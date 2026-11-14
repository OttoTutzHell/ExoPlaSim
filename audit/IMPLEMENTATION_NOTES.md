# Implementation notes: items 1–3 of the tidally locked audit

Companion to `TIDALLY_LOCKED_AUDIT.md`. Three changes, each behind an option that is off by default.

| item | option | where |
|---|---|---|
| 1. Ocean heat transport | `oceandiffusion`, `oceandiffusionzonal` | `oceanmod.f90`, `seamod.f90`, `outmod.f90`, wrapper |
| 2. Convergence tooling | `Model.convergence()`, `python -m exoplasim.convergence` | `exoplasim/convergence.py` |
| 3. Stellar-spectrum near-infrared absorption | `nirabsorption` | `radmod.f90`, wrapper |

## Default behaviour is unchanged

One model year of the THAI Hab 1 configuration, same restart, same executable flags, unmodified
build against modified build with no new option set: all 100 output variables are bit-identical.
A configuration file written with the new options reloads to the same namelists; a default
configuration writes none of the new namelist entries.

## Benchmark

THAI Hab 1 (Fauchez et al. 2020; Sergeev et al. 2022): TRAPPIST-1e, 900 W/m², 2,600 K, 6.1 d,
0.91 R⊕, 9.12 m/s², 1 bar N₂, 400 ppm CO₂, aquaplanet, slab ocean without heat transport.
Here: T21 L10, 30-minute step, 50 m slab, `gp|exp|sp` filter. Departures from the protocol:
blackbody star instead of BT-Settl, ExoPlaSim's own sea-ice and ocean albedos instead of the
prescribed 0.25 and 0.06.

Reference values, read from the paper's own source (Table 4, Sec. IV.2):

| | ExoCAM | LMD-G | ROCKE-3D | UM |
|---|---|---|---|---|
| global mean T_s (K) | 245.6 | 242.2 | 244.0 | 231.6 |
| planetary albedo | 0.24 | 0.22 | 0.27 | 0.28 |
| ice-free ocean fraction | 0.20–0.23 across the four | | | |

Four branches from the same year-13 state, 12 years each, last 5 averaged
(`results/thai_hab1_t21.md`, `results/thai_hab1_t21.png`):

| quantity | control | `oceandiffusion=5e4` | `nirabsorption` | both |
|---|---|---|---|---|
| T_s global (K) | 230.8 | 233.4 | 232.9 | 234.0 |
| T_s dayside / nightside (K) | 258.7 / 202.8 | 261.2 / 205.5 | 255.0 / 210.8 | 256.1 / 212.0 |
| T_s warmest / coldest cell (K) | 295.6 / 141.5 | 291.1 / 149.9 | 289.2 / 146.5 | 285.6 / 150.8 |
| ice-free ocean fraction | 0.219 | 0.300 | 0.173 | 0.227 |
| planetary albedo | 0.307 | 0.301 | 0.308 | 0.300 |
| absorbed shortwave: surface / atmosphere (W/m²) | 112.2 / 43.8 | 112.7 / 44.6 | 88.0 / 67.7 | 89.2 / 68.4 |
| atmosphere's share of absorbed shortwave | 0.28 | 0.28 | 0.44 | 0.43 |
| precipitation, global (mm/yr) | 492 | 503 | 321 | 329 |
| precipitation, substellar cap (mm/yr) | 5,685 | 5,148 | 3,454 | 3,311 |
| TOA net / surface net (W/m²) | −0.3 / −0.1 | −0.4 / −0.1 | −0.6 / −0.3 | −0.4 / −0.1 |
| ocean transport heating, mean magnitude / global mean (W/m²) | — | 11.7 / −0.001 | — | 6.9 / 0.000 |

How to read it against the benchmark:

- Every branch stays inside or at the edge of the four-model temperature range (231.6–245.6 K). The control is 0.8 K below the coldest model; the other three are inside.
- Planetary albedo is 0.30–0.31 in every branch, above all four models (0.22–0.28). Neither change addresses that. It is a cloud and ice-albedo property of this configuration.
- The benchmark cannot say which branch is "right". The four models differ by 14 K with no ocean transport, and none has the true answer.

## Item 1: ocean heat transport by diffusion

**What changed.**
- `hdiffo` used the Earth's radius for every planet. It now uses the model's planet radius (`oceansetradius`, called from `seaini`).
- New namelist factor `HDIFFZON` scales the diffusivity in the east–west direction only (default 1).
- With diffusion on, the heating of the mixed layer by the transport is written to the standard output as `ohdf` (code 330, W/m²).
- Wrapper options `oceandiffusion` (m²/s) and `oceandiffusionzonal`, in `configure()`, `modify()` and the configuration file.

**What it fixes.**
- The model can now carry heat from the substellar ocean to the ice margin. In the benchmark it opens the ice-free area from 0.22 to 0.30, cools the warmest cell by 4.5 K and warms the coldest by 8 K, with the global mean of `ohdf` zero to 10⁻³ W/m².
- The direction matches published dynamic-ocean results for tidally locked planets (Hu & Yang 2014; Del Genio et al. 2019; Yang et al. 2019: more open ocean, warmer nightside). I have not reproduced any of those papers' configurations, so the comparison is qualitative.

**What it does not fix.**
- It is diffusion, not circulation. No equatorial jet in the ocean, no boundary currents, no upwelling, no east–west asymmetry unless `oceandiffusionzonal` is used to stretch it.
- The diffusivity is a free parameter. Transport scales with `oceandiffusion × mldepth`, so fix the mixed-layer depth first. 5×10⁴ m²/s at 50 m was chosen to give a visible response, not derived.
- A simple test that would tell you the value is too high: `ohdf` larger than the local absorbed starlight, or an ice-free nightside.
- The solver is explicit with 100 sub-steps on the root process; I estimate it is stable to about 10⁶ m²/s at T42 with a 30-minute step but tested only 2–5×10⁴ at T21.
- The one-step lag in `mkiflux` noted in the audit is left alone, because fixing it would change default results.

## Item 2: convergence and energy-budget tooling

**What changed.** New module `exoplasim.convergence`:
- `history(workdir)`: one row per output year of global-mean surface temperature, TOA and surface net flux, sea-ice area and volume, land snow, and the ice-free mixed layer's temperature.
- `budget(hist, window)`: mean fluxes, and the energy stored in sea ice (latent heat), land snow and the open mixed layer, from their trends. What is left is reported as unexplained.
- `isconverged(hist, window, criteria)`: tests TOA and surface net flux, the atmosphere's residual, the trends of surface temperature, sea-ice volume and snow, and the unexplained term. Defaults are in `DEFAULT_CRITERIA` and can be overridden.
- `Model.convergence()` and `python -m exoplasim.convergence WORKDIR [WINDOW]`.

**Before and after.** `Model.runtobalance()` cannot return a verdict before year 50–75 and looks only at the slopes of two fluxes. On the control run, model years 8–13, the new report gives:

```
  top of atmosphere net     -0.69
  surface net               -0.44
  storage: sea ice          -0.42   (ice volume +0.043 m/yr)
  surface net - storage     -0.02   (unexplained)
  seaice_trend      0.0065  limit 0.002    FAIL
NOT CONVERGED
```

and on model years 18–25 of the same configuration (`results/convergence_report.txt`): surface net −0.10, sea-ice storage −0.09, unexplained 0.00, all criteria passed.

**What it fixes.** It answers the audit's balance question with a measurement: the surface imbalance in an icy tidally locked run is the latent heat of sea ice still thickening toward the model's 9 m cap, to within 0.02–0.04 W/m². It also gives a spin-up test that works at year 10–20 instead of year 75.

**What it does not fix.**
- It does not change `runtobalance()`; that would alter existing behaviour.
- Soil heat and the heat content of ice and snow are not in the budget. While the ice edge is still moving the unexplained term is large and should be read as "not settled", not as a leak.
- The open-ocean term uses surface temperature, since sea-surface temperature is not in the default output.
- It does not address restarts across resolutions, the `pyfft` build failure under numpy 2, or two smaller wrapper faults found on the way: the executable cache ignores `precision` (a precision-8 restart crashes on a cached precision-4 executable unless `recompile=True`), and two models created at once race on compiling.

## Item 3: near-infrared water-vapour absorption that sees the star

**What changed.** With `nirabsorption=True` (namelist `NSWRNIR=1`) and a stellar temperature or spectrum set:
- The Lacis & Hansen (1974) curve, which is the fraction of *solar* flux absorbed, is mapped onto the near-infrared water-vapour k-distribution of Chou & Lee (1996) in the three bands of Chou & Suarez (1999): 0.70–1.22, 1.22–2.27 and 2.27–10 µm. The table values were copied from the public NASA/GSFC shortwave code as distributed in WRF (`module_ra_gsfcsw.F`); its stated checksum (0.52926) matches.
- Each band is re-weighted by the ratio of the star's flux fraction in that band to a 5,772 K blackbody's.
- A look-up table is built once at start-up (`mknirtab`). The water amount keeps PlaSim's own pressure and temperature scaling.

**Check of the table.** A Python replica (`scripts/nirtab.py`) gives the same band ratios and absorbed fraction as the model prints at start-up (2,600 K: ratios 0.80, 2.75, 5.84; 27.7% of stellar flux absorbed by 1 cm of water against 9.9% for the solar fit). For 5,772 K the ratio to the original curve is 1.00 at every water amount.

| star | absorbed fraction relative to the solar fit, at 0.01 / 0.1 / 1 / 10 cm of water |
|---|---|
| 5,772 K blackbody | 1.00 / 1.00 / 1.00 / 1.00 |
| 3,800 K blackbody | 2.15 / 2.00 / 1.84 / 1.70 |
| 2,825 K blackbody (Ki-Engir, 1,600 W/m² star) | 3.40 / 2.98 / 2.58 / 2.23 |
| 2,600 K blackbody | 3.81 / 3.28 / 2.79 / 2.36 |
| TRAPPIST-1 spectrum shipped with ExoPlaSim (`trap`) | 3.29 / 2.99 / 2.64 / 2.31 |

**What it fixes.** In the benchmark the atmosphere's share of absorbed starlight rises from 28% to 44%, the surface receives 22% less, global precipitation falls by a third, the dayside cools by 3.7 K and the nightside warms by 8 K. The audit's guess of "roughly double" was low: for a 2,600–2,900 K star the factor is 2.2 to 3.8.

**What it does not fix, and where it can be wrong.**
- Only the weighting between three bands changes. Inside each band the spectrum is assumed solar-shaped. A real M dwarf has its own water absorption in the same bands, so the starlight is already depleted exactly where the planet's water vapour absorbs. This approximation therefore tends to **overestimate** absorption. Using a real spectrum through `starspec` corrects the band totals (the TRAPPIST-1 row above) but not the alignment inside bands.
- The benchmark does not settle it. Sergeev et al. report that the four models differ in gas absorption among themselves and that one over-absorbs. I found no tabulated Hab 1 surface shortwave to compare against; the paper gives it only as maps.
- Still no shortwave absorption by CO₂, and the cloud optics are untouched.
- Planetary albedo does not move (0.31), so this is a redistribution of absorbed energy between surface and atmosphere, not a change in how much is absorbed.
- It lowers the ice-free fraction (0.22 → 0.17) because less starlight reaches the ocean. With ocean diffusion on as well it returns to 0.23.

## What this means for a tuned state

The vault's cloud settings (`tswr2` tripled) were chosen with the old absorption. With `nirabsorption` on, less light reaches the surface for the same clouds, so those settings will over-cool the dayside surface and should be relaxed toward stock before anything else is retuned. The Ki-Engir demonstration below shows the size of the effect.

## Ki-Engir demonstration (T21, in progress)

Not a validation. Settings: the vault's `f1800_rc070` (1,800 W/m², 2,879 K, `rcritmod` 0.70, `tswr1/2/3` 0.3/0.18/0.002), driver `scripts/kiengir_t21.py`.

**From the warm, ice-free restart (`ki_t21_clouds_hard_y20.rest`, 311 K, 195 kg/m² of column water), 10 years, last 3 averaged:**

| quantity | reference | `nirabsorption` |
|---|---|---|
| T_s global (K) | 286.6, cooling ~3 K/yr | 318.0, near balance |
| planetary albedo | 0.507 | 0.257 |
| absorbed shortwave, global (W/m²) | 222 | 335 |
| atmosphere's share of absorbed shortwave | 0.41 | 0.66 |
| column water vapour (kg/m²) | 48 | 398 |
| dayside cloud cover | 0.77 | 0.55 |

From this start the two branches part: the reference slides toward the icy state as the vault found, and the near-infrared branch settles on a hot, very moist state that absorbs more than canon's runaway bound (295 ± 15 W/m²). In a moist column the stronger absorption takes the light before it reaches the cloud deck, the albedo collapses, and the warming feeds itself. The approximation's bias (see item 3) is in the direction of exaggerating this.

**From the reference run's cooler year-30 state (284 K, 10% sea ice, ~35 kg/m² of water), 10 years, last 3 averaged** (`results/kiengir_t21_demo.md`, `.png`):

| quantity | reference | `nirabsorption` |
|---|---|---|
| T_s global (K) | 272.7, cooling 0.3 K/yr | 313.9, still warming 1.4 K/yr |
| T_s dayside / nightside (K) | 286.7 / 258.7 | 321.2 / 306.7 |
| warmest cell (K) | 321.6 | 360.3 |
| planetary albedo | 0.525 | 0.258 |
| absorbed shortwave, global (W/m²) | 214 | 334 |
| ice-free ocean fraction | 0.57 | 1.00 |
| precipitation, global / substellar cap (mm/yr) | 883 / 6,320 | 1,189 / 1,636 |
| column water vapour (kg/m²) | 24 | 337 |

The reference reproduces the vault's `f1800_rc070` track. The near-infrared branch turns around within a year and warms. So this is not start-state hysteresis: **with `nirabsorption` on, the tuned Ki-Engir settings have no cold state at 1,800 W/m², at least down to 284 K.** Whether a fully iced 270 K state would hold was not tested.

Why the benchmark stayed moderate and Ki-Engir did not: twice the stellar flux and a much moister dayside. Vapour above the cloud deck absorbs the light the tuned clouds were reflecting, so the tripled `tswr2` loses most of its leverage.

How far to trust it:
- The sign is physical and published for cool stars (lower albedo, weaker substellar cloud, earlier moist greenhouse). The size is uncertain and biased high by the in-band approximation.
- 1,800 W/m² around a ~2,900 K star is, from memory of the Kopparapu et al. (2016, 2017) inner-edge results, close to where three-dimensional models put the moist-greenhouse transition for tidally locked planets. I have not re-read those papers for this. If that is right, the flux is the parameter under strain, not a model fault.

**A Ki-Engir ocean-diffusion branch (`oceandiffusion=2e4`, no near-infrared change) was started from the same state and had not finished when this was written.** Run `scripts/compare.py` on `runs/ki_ohd2` if it is needed.

## Star temperature

`nirabsorption` takes whatever `startemp` is set, or a `starspec` file; nothing is tied to a particular value. `scripts/nirtab.py TEFF ...` prints the factors for any temperature. 2,566 K and 2,600 K differ by about 1% (3.87/3.32/2.82/2.38 against 3.81/3.28/2.79/2.36 at 0.01/0.1/1/10 cm of water). The table is capped so the band transmissivity stays positive for any star.

## Recommended run configuration for Ki-Engir

| setting | T21 | T42 |
|---|---|---|
| grid, levels, step | T21 L10, 30 min | T42 L10, 30 min (20 if it crashes near the terminator range) |
| `physicsfilter` | `gp|exp|sp` | `gp|exp|sp` |
| `th2oc`, `gamma` | defaults | pin to 0.024 and 0.01, as the vault does |
| `precision` | 8, with `recompile=True` on first use | same |
| `mldepth` | 50 m, held fixed | same |
| `oceandiffusion` | 2×10⁴ m²/s as a starting value; bracket with 1×10⁴ and 5×10⁴ | the value chosen at T21 |
| `oceandiffusionzonal` | 1 unless a lobster-shaped open ocean is wanted | same |
| `nirabsorption` | **decide first; see below** | same as T21 |
| `glaciers`, `maxsnow` | off, default; mask `mrro` where snow is at the 5 m cap | same |
| start | branch from an existing restart | transplant from T21 |
| stop | `Model.convergence(window=10)` passes | same |
| averaging | last 10 years | last 10 years |

**The `nirabsorption` decision.** Two coherent options; mixing them is not one.

| | keep it off | turn it on |
|---|---|---|
| what you keep | every tuned setting and run in the vault | — |
| what you accept | a shortwave scheme that treats a 2,900 K star like the Sun, with clouds tuned to stand in for the missing absorption; rain about a third too strong on the benchmark | the tuned state is gone; clouds, and probably the 1,800 W/m² flux, must be re-chosen from scratch |
| failure mode | canon rests on compensating errors | no setting inside the canon bounds may give a temperate dayside; the absorption is itself biased high |
| first test | none needed | T21, `nirabsorption=True`, stock `tswr2/3`, at 1,600 and 1,400 W/m², 15 years each from the year-30 state: does a cold or temperate branch exist |

Either way the convergence tool and the ocean diffusion can be used as they are.
