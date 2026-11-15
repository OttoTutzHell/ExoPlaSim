"""Build the embedded data for the Ki-Engir planet sim.

Ki frame: lat phi (N+), lon lam (east+), SSP at (0,0). Grid: 2 deg, NX=180 lon x NY=90 lat,
cell centres lon -179..179, lat -89..89 (row 0 = south).
Relief: z_ki (canon rescaled, 6 arcmin, Earth frame) from
'Eversun/Ki-Engir De-terraced Relief (6-min).nc', sampled 6x6 per cell via the
Earth-Ki-Engir Map Transform.
GCM: ExoPlaSim T42 year-70 (rc067) and year-65 (rc070) surface temperature, already Ki frame.
"""
import json, base64, sys
import numpy as np
from scipy.io import netcdf_file

V = sys.argv[1]  # vault root
OUT = sys.argv[2]
RES = 2.0
NX, NY = int(360 / RES), int(180 / RES)

f = netcdf_file(f"{V}/Eversun/Ki-Engir De-terraced Relief (6-min).nc", "r", mmap=False)
Z = f.variables["z_ki"][:].copy().astype(np.float32)  # (1801, 3600) lat -90..90, lon -180..179.9


def ki_to_earth(phi, lam):
    p, l = np.radians(phi), np.radians(lam)
    pe = -np.arcsin(np.cos(p) * np.cos(l))
    le = np.arctan2(np.cos(p) * np.sin(l), np.sin(p))
    return np.degrees(pe), np.degrees(le)


SUB = 6
off = (np.arange(SUB) + 0.5) / SUB * RES - RES / 2
lonc = -180 + RES / 2 + RES * np.arange(NX)
latc = -90 + RES / 2 + RES * np.arange(NY)
zmean = np.zeros((NY, NX), np.float32)
land = np.zeros((NY, NX), np.float32)
zland = np.zeros((NY, NX), np.float32)
zmax = np.zeros((NY, NX), np.float32)
for j, la in enumerate(latc):
    PH = (la + off)[:, None] + 0 * off[None, :]
    for i, lo in enumerate(lonc):
        LA = (lo + off)[None, :] + 0 * off[:, None]
        pe, le = ki_to_earth(PH, LA)
        r = np.clip(np.rint((pe + 90) / 0.1).astype(int), 0, 1800)
        c = np.rint((le + 180) / 0.1).astype(int) % 3600
        zz = Z[r, c]
        zmean[j, i] = zz.mean()
        lm = zz > 0
        land[j, i] = lm.mean()
        zland[j, i] = zz[lm].mean() if lm.any() else 0
        zmax[j, i] = zz.max()

# cold-pool floor: lowest surface (sea counts as 0) within +-3 cells, on the sphere-ish window
surf = np.where(land >= 0.5, np.maximum(zland, 0), 0.0)
floor = np.full_like(surf, 1e9)
for dj in range(-3, 4):
    for di in range(-3, 4):
        jj = np.clip(np.arange(NY) + dj, 0, NY - 1)
        floor = np.minimum(floor, np.roll(surf[jj], -di, axis=1))
zrel = np.maximum(surf - floor, 0)


def b64(a, dtype):
    return base64.b64encode(np.ascontiguousarray(a, dtype=dtype).tobytes()).decode()


# GCM fields: lat descending 87.9..-87.9, lon 0..357.2 (Ki frame, east+). Regrid to our grid bilinear.
def regrid(ts, glat, glon):
    glat = np.asarray(glat)[::-1]
    ts = np.asarray(ts)[::-1]
    glon2 = np.concatenate([glon - 360, glon, glon + 360])
    ts3 = np.concatenate([ts, ts, ts], axis=1)
    out = np.zeros((NY, NX))
    for j, la in enumerate(latc):
        la_c = np.clip(la, glat[0], glat[-1])
        jj = np.searchsorted(glat, la_c) - 1
        jj = np.clip(jj, 0, len(glat) - 2)
        wy = (la_c - glat[jj]) / (glat[jj + 1] - glat[jj])
        row = ts3[jj] * (1 - wy) + ts3[jj + 1] * wy
        out[j] = np.interp(lonc, glon2, row)
    return out


# GCM fields: two ExoPlaSim T21 runs with the stellar near-infrared absorption on and stock clouds
# (exported by anchors.py). Slot "rc067" is the primary comparison run, "rc070" the second.
GA, GB = sys.argv[3], sys.argv[4]          # e.g. gcm_tr_f1100.npz gcm_tg_f1200.npz
ANCH = json.load(open(sys.argv[5]))        # anchors.json
gcm = {}
for key, fn in [("rc067", GA), ("rc070", GB)]:
    d = np.load(fn)
    g = regrid(d["ts"], d["lat"], d["lon"])
    gcm[key] = g
d67 = np.load(GA)
gcm_pr = np.maximum(regrid(d67["pr"], d67["lat"], d67["lon"]), 0)    # mm/yr
gcm_tas = regrid(d67["tas"], d67["lat"], d67["lon"])                 # K, 2 m air

# area weights & summaries
w = np.cos(np.radians(latc))[:, None] * np.ones((1, NX))
LON = lonc[None, :] * np.ones((NY, 1))
LAT = latc[:, None] * np.ones((1, NX))
mu = np.cos(np.radians(LAT)) * np.cos(np.radians(LON))
for k, g in gcm.items():
    print(k, "global %.1f day %.1f night %.1f" % ((g * w).sum() / w.sum(), (g * w)[mu > 0].sum() / w[mu > 0].sum(),
                                                 (g * w)[mu <= 0].sum() / w[mu <= 0].sum()))
print("land frac %.3f" % ((land * w).sum() / w.sum()), "dayside land %.3f" % ((land * w)[mu > 0].sum() / w[mu > 0].sum()))

# ---- routing grid: 0.5 deg in the Ki frame, 4x4 per climate cell, land mean height and land fraction
RR = 0.5
RX, RY = int(360 / RR), int(180 / RR)
SUBR = 4
offr = (np.arange(SUBR) + 0.5) / SUBR * RR - RR / 2
rlon = -180 + RR / 2 + RR * np.arange(RX)
rlat = -90 + RR / 2 + RR * np.arange(RY)
LA = (rlat[:, None, None, None] + offr[None, None, :, None]) * np.ones((RY, RX, SUBR, SUBR))
LO = (rlon[None, :, None, None] + offr[None, None, None, :]) * np.ones((RY, RX, SUBR, SUBR))
pe, le = ki_to_earth(LA, LO)
r = np.clip(np.rint((pe + 90) / 0.1).astype(int), 0, 1800)
c = np.rint((le + 180) / 0.1).astype(int) % 3600
zz = Z[r, c]
lm = zz > 0
rland = lm.mean(axis=(2, 3))
rz = np.where(lm, zz, 0).sum(axis=(2, 3)) / np.maximum(lm.sum(axis=(2, 3)), 1)
rz = np.where(rland >= 0.5, np.maximum(rz, 1), np.minimum(zz.mean(axis=(2, 3)), 0))

data = dict(
    NX=NX, NY=NY, RES=RES,
    zmean=b64(np.rint(zmean), np.int16), land=b64(np.rint(land * 255), np.uint8),
    zland=b64(np.rint(zland), np.int16), zmax=b64(np.rint(zmax), np.int16), zrel=b64(np.rint(zrel), np.int16),
    gcm067=b64(np.rint((gcm["rc067"] - 150) * 100), np.uint16),
    gcm070=b64(np.rint((gcm["rc070"] - 150) * 100), np.uint16),
    gcmPr=b64(np.clip(np.rint(gcm_pr), 0, 65535), np.uint16), gcmTas=b64(np.rint((gcm_tas - 150) * 100), np.uint16),
    RX=RX, RY=RY, RRES=RR, rz=b64(np.rint(rz), np.int16), rland=b64(np.rint(rland * 255), np.uint8),
)
def meta(fn):
    run = fn.split("gcm_")[-1][:-4]; a = [x for x in ANCH if x["run"] == run][0]
    short = ("T-1 " if a["star"].startswith("TRAP") else "Teeg. ") + "%.2f d" % a["period"]
    return dict(short=short, label="%s, %.2f-day orbit, %d W/m²" % (a["star"], a["period"], round(a["flux"])),
                means=[a["G"], a["D"], a["N"], a["maxC"] + 273.15], albedo=a["alb"], iceFree=a["iceFree"])
data["gcmMeta"] = dict(a=meta(GA), b=meta(GB))
data["anchors"] = ANCH
data["stars"] = [
    dict(name="TRAPPIST-1", T=2566, L=5.53e-4, M=0.0898, period=5.0, sweep=[4.5, 6.5, 0.125],
         spec=dict(fvis=0.015, ratio=[1.046, 3.192, 4.331]),
         planets=[["b", 1.511], ["c", 2.422], ["d", 4.049], ["e", 6.101], ["f", 9.208], ["g", 12.352], ["h", 18.773]],
         note="Agol et al. 2021 values; band weights from the TRAPPIST-1 spectrum shipped with ExoPlaSim."),
    dict(name="Teegarden's Star", T=3034, L=7.22e-4, M=0.097, period=5.68, sweep=[5.0, 7.0, 0.125], spec=None,
         planets=[["b", 4.906], ["c", 11.416], ["d", 26.13]],
         note="Dreizler et al. 2024 values; blackbody colour."),
]
with open(OUT, "w") as fh:
    fh.write("const KI_DATA = " + json.dumps(data) + ";\n")
print("bytes", sum(len(v) for v in data.values() if isinstance(v, str)))
