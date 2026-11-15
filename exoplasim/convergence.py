"""
Spin-up and energy-budget diagnostics for slowly-equilibrating runs.

``Model.runtobalance()`` judges equilibrium from the drift of the top-of-atmosphere and
surface net energy fluxes alone. On planets with large permanent ice cover (for example
synchronously-rotating planets) those fluxes can stay non-zero for decades because the
surface is still storing or releasing energy, mostly as latent heat of thickening sea ice.
The functions here separate *storage* from *leak* and test the slow reservoirs directly.

Typical use::

    >>> import exoplasim.convergence as conv
    >>> hist = conv.history("path/to/workdir")          # one row per output year
    >>> print(conv.report(hist))                        # table plus a verdict
    >>> ok, checks = conv.isconverged(hist)

or from the command line::

    python -m exoplasim.convergence path/to/workdir

Nothing here changes how the model runs.
"""
import glob
import sys

import numpy as np

import exoplasim.gcmt as gcmt

RHO_SEA = 1030.0    # kg m-3, as in oceanmod.f90
CP_SEA = 4180.0     # J kg-1 K-1
RHO_ICE = 920.0     # kg m-3, as in icemod.f90
L_FUSION = 3.28e5   # J kg-1
RHO_WATER = 1000.0  # kg m-3 (snow depth is in liquid water equivalent)

#: Default convergence criteria. Each is the largest acceptable absolute value over the window.
DEFAULT_CRITERIA = {"toa": 0.5,           # W m-2, |mean net top-of-atmosphere flux|
                    "surface": 0.5,       # W m-2, |mean net surface flux|
                    "atmosphere": 1.0,    # W m-2, |mean (toa - surface)|, the atmosphere's residual
                    "ts_trend": 0.05,     # K yr-1, global-mean surface temperature
                    "seaice_trend": 0.002,  # yr-1, fractional change of sea-ice volume per year
                    "snow_trend": 0.01,   # m yr-1 liquid equivalent, global-mean land snow depth
                    "unexplained": 0.5}   # W m-2, |surface flux - diagnosed storage|


def _fields(filename, keys):
    data = gcmt.load(filename)
    lat = np.array(data.variables["lat"][:], dtype=float)
    out = {}
    for key in keys:
        if key in data.variables:
            var = np.array(data.variables[key][:], dtype=float)
            if var.ndim == 3:
                var = var.mean(axis=0)
            out[key] = var
    try:
        data.close()
    except Exception:
        pass
    return lat, out


def _gmean(field, lat):
    weights = np.cos(np.radians(lat))[:, None] * np.ones(field.shape)
    return float((field * weights).sum() / weights.sum())


def yearly(filename):
    """Global annual means of the quantities the convergence test uses, from one output file.

    Parameters
    ----------
    filename : str
        A postprocessed annual output file (any format :py:func:`exoplasim.gcmt.load` reads).

    Returns
    -------
    dict
        ``ts`` (K), ``toa`` and ``surface`` (W m-2, positive into the planet/surface),
        ``seaice_area`` (fraction of the globe), ``seaice_volume`` (m, global-mean thickness),
        ``snow`` (m liquid equivalent, global mean), ``open_area``, ``open_ts`` and ``open_mld``
        (area fraction, mean temperature and mean depth of the ice-free mixed layer), and ``ohdf_abs`` (W m-2, global mean
        magnitude of ocean diffusion heating, if present).
    """
    keys = ["ts", "ntr", "hfns", "sic", "sit", "snd", "mld", "lsm", "ohdf"]
    lat, f = _fields(filename, keys)
    row = {"file": filename}
    row["ts"] = _gmean(f["ts"], lat)
    row["toa"] = _gmean(f["ntr"], lat)
    row["surface"] = _gmean(f["hfns"], lat)
    lsm = f.get("lsm", np.zeros_like(f["ts"]))
    sea = lsm < 0.5
    sic = np.where(sea, f.get("sic", np.zeros_like(lsm)), 0.0)
    sit = np.where(sea, f.get("sit", np.zeros_like(lsm)), 0.0)
    row["seaice_area"] = _gmean(sic, lat)
    row["seaice_volume"] = _gmean(sit, lat)
    row["snow"] = _gmean(np.where(sea, 0.0, f.get("snd", np.zeros_like(lsm))), lat)
    if "mld" in f:
        openwater = sea & (sic < 0.5)
        row["open_area"] = _gmean(openwater * 1.0, lat)
        if row["open_area"] > 0:
            row["open_ts"] = _gmean(np.where(openwater, f["ts"], 0.0), lat) / row["open_area"]
            row["open_mld"] = _gmean(np.where(openwater, f["mld"], 0.0), lat) / row["open_area"]
        else:
            row["open_ts"] = np.nan
            row["open_mld"] = 0.0
    else:
        row["open_area"] = np.nan
        row["open_ts"] = np.nan
        row["open_mld"] = 0.0
    if "ohdf" in f:
        row["ohdf_abs"] = _gmean(np.abs(f["ohdf"]), lat)
        row["ohdf_net"] = _gmean(f["ohdf"], lat)
    return row


def history(workdir, pattern="MOST.[0-9]*", extension=None):
    """Collect :py:func:`yearly` for every annual output file in a run directory.

    Parameters
    ----------
    workdir : str
        The model's working directory.
    pattern : str, optional
        Glob pattern for annual output files, without the extension.
    extension : str, optional
        File extension such as ``".nc"``. By default every supported extension is tried.

    Returns
    -------
    list of dict
        One entry per year, in order.
    """
    exts = [extension] if extension else [".nc", ".npz", ".hdf5", ".he5", ".h5"]
    files = []
    for ext in exts:
        files += [f for f in glob.glob("%s/%s%s" % (workdir, pattern, ext)) if "_metadata" not in f]
    return [yearly(f) for f in sorted(files)]


def _trend(values):
    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return np.nan
    return float(np.polyfit(np.arange(len(values)), values, 1)[0])


def budget(hist, window=10, yearseconds=360.0 * 86400.0):
    """Energy budget of the last ``window`` years: fluxes, diagnosed storage and what is left.

    Storage is estimated from the trends of sea-ice volume, land snow and the heat content of
    the ice-free mixed layer. Soil heat and the heat content of ice and snow are not included,
    nor is the heat released where open water freezes over, so a remainder is expected while
    the ice edge is still moving.

    Parameters
    ----------
    hist : list of dict
        Output of :py:func:`history`.
    window : int, optional
        Number of years at the end of the run to assess.
    yearseconds : float, optional
        Length of one output year in seconds (ExoPlaSim's default year is 360 days of 86400 s).

    Returns
    -------
    dict
        All fluxes in W m-2 as global means. ``storage_*`` are positive when the reservoir
        gains energy.
    """
    rows = hist[-window:]
    b = {"years": len(rows)}
    b["toa"] = float(np.mean([r["toa"] for r in rows]))
    b["surface"] = float(np.mean([r["surface"] for r in rows]))
    b["atmosphere"] = b["toa"] - b["surface"]
    b["ts_trend"] = _trend([r["ts"] for r in rows])
    icetrend = _trend([r["seaice_volume"] for r in rows])       # m yr-1
    snowtrend = _trend([r["snow"] for r in rows])               # m yr-1
    b["seaice_volume"] = float(np.mean([r["seaice_volume"] for r in rows]))
    b["seaice_trend"] = icetrend / b["seaice_volume"] if b["seaice_volume"] > 0 else 0.0
    b["seaice_growth"] = icetrend
    b["snow_trend"] = snowtrend
    # growing ice or snow releases latent heat, i.e. the surface loses energy
    b["storage_seaice"] = -RHO_ICE * L_FUSION * icetrend / yearseconds
    b["storage_snow"] = -RHO_WATER * L_FUSION * snowtrend / yearseconds
    # mixed layer: warming of the ice-free ocean, weighted by its mean area and depth
    oceantrend = _trend([r["open_ts"] for r in rows])           # K yr-1
    if np.isnan(oceantrend):
        b["storage_ocean"] = 0.0
    else:
        area = float(np.mean([r["open_area"] for r in rows]))
        depth = float(np.mean([r["open_mld"] for r in rows]))
        b["storage_ocean"] = RHO_SEA * CP_SEA * depth * area * oceantrend / yearseconds
    b["storage"] = b["storage_seaice"] + b["storage_snow"] + b["storage_ocean"]
    b["unexplained"] = b["surface"] - b["storage"]
    return b


def isconverged(hist, window=10, criteria=None, yearseconds=360.0 * 86400.0):
    """Test the slow reservoirs as well as the energy fluxes.

    Parameters
    ----------
    hist : list of dict
        Output of :py:func:`history`.
    window : int, optional
        Number of years at the end of the run to assess. Fewer available years is a failure.
    criteria : dict, optional
        Overrides for entries of :py:data:`DEFAULT_CRITERIA`. Set an entry to None to skip it.

    Returns
    -------
    (bool, dict)
        The verdict, and for each criterion a tuple ``(value, limit, passed)``.
    """
    limits = dict(DEFAULT_CRITERIA)
    if criteria:
        limits.update(criteria)
    if len(hist) < window:
        return False, {"years": (len(hist), window, False)}
    b = budget(hist, window=window, yearseconds=yearseconds)
    checks = {}
    for key, limit in limits.items():
        if limit is None:
            continue
        value = b[key]
        checks[key] = (value, limit, bool(abs(value) <= limit))
    return all(c[2] for c in checks.values()), checks


def report(hist, window=10, criteria=None, yearseconds=360.0 * 86400.0):
    """Human-readable table of the yearly history, the budget and the verdict."""
    lines = ["%4s %8s %8s %8s %8s %9s %9s" % ("year", "ts", "toa", "surface", "toa-sfc",
                                               "ice area", "ice vol")]
    for n, r in enumerate(hist):
        lines.append("%4d %8.2f %8.2f %8.2f %8.2f %9.3f %9.3f" % (
            n, r["ts"], r["toa"], r["surface"], r["toa"] - r["surface"],
            r["seaice_area"], r["seaice_volume"]))
    if len(hist) >= 2:
        w = min(window, len(hist))
        b = budget(hist, window=w, yearseconds=yearseconds)
        lines.append("")
        lines.append("Energy budget over the last %d years (W m-2, global mean):" % w)
        lines.append("  top of atmosphere net   %7.2f" % b["toa"])
        lines.append("  surface net             %7.2f" % b["surface"])
        lines.append("  atmosphere residual     %7.2f   (toa - surface)" % b["atmosphere"])
        lines.append("  storage: sea ice        %7.2f   (ice volume %+.3f m/yr)" % (
            b["storage_seaice"], b["seaice_growth"]))
        lines.append("  storage: land snow      %7.2f" % b["storage_snow"])
        lines.append("  storage: open mixed layer %5.2f" % b["storage_ocean"])
        lines.append("  surface net - storage   %7.2f   (unexplained)" % b["unexplained"])
        lines.append("  surface temperature trend %+.3f K/yr" % b["ts_trend"])
    ok, checks = isconverged(hist, window=window, criteria=criteria, yearseconds=yearseconds)
    lines.append("")
    for key, (value, limit, passed) in checks.items():
        lines.append("  %-13s %10.4f  limit %-8g %s" % (key, value, limit, "ok" if passed else "FAIL"))
    lines.append("CONVERGED" if ok else "NOT CONVERGED")
    return "\n".join(lines)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: python -m exoplasim.convergence WORKDIR [WINDOW]")
    _window = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    print(report(history(sys.argv[1]), window=_window))
