# ExoPlaSim audit for tidally locked, slow-rotating planets with extreme topography

Audited tree: `5f8dff4` (ExoPlaSim 3.4.2). No model code is changed by this document.
The reference planet is Ki-Engir (synchronous, 6.097 d, 0.92 R⊕, 7.99 m/s², 1 bar N₂,
1,000 ppm CO₂, 2,825–2,879 K host, 26% land, relief to 8–10 km), with the run history in
the Eversun vault's `tools/exoplasim/PLAN.txt`.

Every claim below is tagged:

- **[code]** read directly in the source, with file and line.
- **[measured]** computed here, from a run or an offline calculation in `audit/scripts/`.
- **[literature]** a published result. The reference checker was unavailable during this
  audit, so citations are from memory and unverified unless marked otherwise.
- **[inferred]** my reasoning from the above. Treat as a hypothesis.

## 1. Which outputs carry the weight

From the vault (`Ki-Engir GCM Run (Scope)` §5 and §7, `Open for Otto` K60):

1. Dayside land surface temperature at the polity sites. This is the score a canon run is judged on.
2. Precipitation pattern and P−E over land. It feeds aridity, climate classes, runoff routing, lakes, soils. It is the worst-fitting field (land pattern r = 0.18 against the hand-built map).
3. Near-surface wind direction, and the equatorial-versus-mid-latitude jet structure (K58).
4. Sea-ice edge and the frozen share of the nightside (a bound, not a target).
5. Nightside temperatures are explicitly low priority.

The ranking in §3 weights problems by how hard they hit items 1–3.

## 2. Findings by subsystem

### 2.1 Ocean heat transport

- **[code]** The slab ocean already contains a horizontal heat diffusion scheme, switched off by default: `nhdiff` and `hdiffk` in `oceanmod.f90:41,51`, applied in `mksst` (`:844`) through `hdiffo` (`:1295`). The vault's note that ExoPlaSim has no ocean heat transport option is wrong for the Fortran core and right for the Python wrapper, which exposes neither parameter. They can be reached today with `otherargs={"NHDIFF@oceanmod_namelist": 1, "HDIFFK@oceanmod_namelist": ...}`.
- **[code]** The scheme is flux-form on the Gaussian grid with zero flux across coasts, so it conserves heat. Heat diffused under sea ice is passed to the ice as basal melt (`:943–949`), which is the physically relevant pathway for a tidally locked ice margin.
- **[code]** Defects:
  1. It uses a hard-coded Earth radius (`PLARAD=6.371E6`, `:24`, used at `:1317`). On a 0.92 R⊕ planet the effective diffusivity is 15% low; on a 1.6 R⊕ planet 2.6 times too high.
  2. The default `hdiffk = 1e3 m²/s` is negligible for a 50 m slab. Scaling: transport ≈ ρ c_p h K ∇T × perimeter; an Earth-like ~1–2 PW needs K ~ 1–2 × 10⁴ m²/s at 50 m **[inferred]**. Only the product h·K matters, so changing `mldepth` silently changes the transport.
  3. The diffusive flux `yqhd` is accumulated but written only to the separate ocean output file, not to the main output, so the transport cannot be diagnosed from a normal run.
  4. It is solved explicitly in 100 sub-steps on the root process. It is stable for K up to roughly 10⁶ m²/s at T42 with a 30-minute step **[inferred from the CFL estimate]**.
  5. `mkiflux` tests the previous step's SST rather than the updated one (`:997`), a one-step lag in ice formation. Harmless, but it is an error.
- **[measured]** The switch works. See §4.1.
- **[literature]** Dynamic oceans on tidally locked planets carry heat mainly zonally along the equator, widen the open-ocean region into a "lobster" shape, warm the nightside and thin its ice (Hu & Yang 2014, PNAS; Del Genio et al. 2019, Astrobiology; Yang et al. 2019, ApJ).
- **What a patch cannot do.** Isotropic diffusion spreads heat down-gradient in every direction. It cannot produce the east–west asymmetry of real equatorial currents, western boundary currents along Ki-Engir's coasts, or upwelling. K is a tuning parameter with no first-principles value. Whether any ocean transport wets the western heartland is unknown; the mechanism there may be atmospheric.

### 2.2 Shortwave radiation and the stellar spectrum

- **[code]** The star enters the shortwave scheme through four things only: the split of flux between two bands at 0.75 µm (`radmod.f90:367–368`), band-weighted surface albedos for snow, ice and ground (`:377–520`), a spectrum-weighted Rayleigh coefficient (`:302`), and the cloud optics, which differ by band.
- **[code]** Clear-sky near-infrared absorption by water vapour is the Lacis & Hansen (1974) broadband function, `2.9y/((1+141.5y)^0.635+5.925y)`, divided by the band-2 fraction (`:1906–1907`, `:2001–2003`). That function is the fraction of **total solar** flux absorbed, fitted for the Sun's spectrum. Dividing by `zsolar2` and multiplying by band-2 flux returns the same fraction of total stellar flux whatever the star. The absorption does not respond to the stellar spectrum at all.
- **[code]** There is no shortwave absorption by CO₂, and none by water vapour in band 1.
- **[inferred]** For a 2,825 K star, about 90% of the flux is beyond 0.75 µm and most of it lies in or near the 0.94, 1.1, 1.4, 1.9 and 2.7 µm water bands. Real clear-sky absorption should be roughly twice the solar-spectrum fraction. The model therefore delivers too much starlight to the dayside surface and too little to the dayside atmosphere. Expected consequences: dayside surface too warm for a given absorbed flux, evaporation and precipitation too strong, static stability too weak.
- **[inferred]** This overlaps with the vault's tuning history. The stock model overheated with 578–715 W/m² reaching the substellar surface, and the fix was to triple the band-2 cloud back-scatter (`tswr2` 0.065 → 0.18). Part of what the clouds are being asked to reflect should have been absorbed aloft instead. Reflection and absorption have opposite effects on planetary albedo and on atmospheric heating, so the tuned state gets the surface budget roughly right for the wrong reason.
- **[measured]** In the THAI Hab 1 analogue (§4), year 1, the atmosphere absorbs 43 W/m², 19% of incident, clouds included.
- **[literature]** Redder stars reduce surface shortwave and weaken the hydrological cycle in models with resolved spectra (Eager-Nash et al. 2020, A&A; Shields et al. 2013, Astrobiology).
- **Longwave.** **[code]** Sasamori (1968) broadband emissivities for H₂O, CO₂ and O₃, with optional pressure broadening by surface pressure (`npbroaden`, `:2176`). At 1 bar and 1,000 ppm this is inside its design envelope. It is less trustworthy in the 0.5 bar branch and in the hot stock state (220 kg/m² of column water). No CO₂ condensation; the vault already checks the frost point by hand.

### 2.3 Sea ice, land ice, and the calving question

- **[code]** Sea ice is thermodynamic only. There is no ice advection, no rheology, no export.
- **[code]** Sea-ice thickness is capped at `xmaxd = 9 m` (`icemod.f90:58`). At the cap, conduction through the ice is set to zero (`:1751`), so capped ice is thermally decoupled from the ocean, which behaves like very thick ice. Any thickness above the cap is removed and its latent heat is taken from every ocean cell with thinner ice, weighted by (cap − thickness), globally (`getiflx`, `:1948–1978`). This is a crude, non-local stand-in for ice export and melt.
- **[code]** Land snow is capped at `dsmax = 5 m` water equivalent (`landmod.f90:56`, `:839–848`). Excess is sent to soil water and from there to runoff, with a compensating cooling. With `glaciers` on, snow deeper than 30 m water equivalent becomes a glacier with no flow (`glaciermod.f90:199`).
- **[code]** The model has no ocean mass budget. Sea level and ocean volume are fixed.
- **Answer to the canon question.** The model cannot represent terminator ice calving under load, and it cannot test whether the dayside ocean is cold-trapped. It assumes it is not: whatever water is deposited on the nightside is returned instantly and invisibly. Canon's ruling is therefore hard-wired, not simulated. That is consistent with canon but is no evidence for it. **[literature]** Yang, Liu, Hu & Abbot (2014, ApJL) found with a model that does have ocean and sea-ice dynamics that trapping needs large nightside continents and a small water inventory.
- **A consequence for the hydrology pipeline.** **[inferred]** On any land cell sitting at the snow cap, `mrro` equals the snowfall rate and is ice export, not liquid river flow. Routing it as rivers on the nightside and terminator land invents discharge. The pipeline should mask or relabel it.

### 2.4 Topography in a spectral model

- **[code]** The "physics filter" is not confined to physics. `skgpsp` and `skspgp` multiply every grid-to-spectral and spectral-to-grid transform (`legmod.f90:197–416`), including the orography's (`surfmod.f90:403–408`) and the dynamics'. With the tidally locked default `gp|exp|sp`, the exponential filter exp(−8(n/N)⁸) is applied in both directions.
- **[measured]** Offline transform of the Ki-Engir boundary files (`audit/scripts/gibbs.py`), which reproduces the model's own output topography to within 2%:

  | | T21 | T42 |
  |---|---|---|
  | highest cell in the boundary file | 5,726 m | 6,568 m |
  | after truncation only | 4,221 m | 5,866 m |
  | what the dynamics see (exp filter once) | 3,388 m | 4,947 m |
  | ocean cells: rms height / min / max | 212 / −349 / +1,946 m | 146 / −444 / +1,786 m |
  | ocean area with \|z\| > 200 m | 11.3% | 6.5% |
  | land rms height error | 604 m | 491 m |

- **[measured]** The "ocean" errors are dominated by mountains smeared across coasts (up to +1.8 km over sea cells), not by oscillating ripples. Lanczos or Cesàro filtering removes the negative lobes but lowers the peaks further.
- Canon's 8–10 km summits become ≤ 5 km at T42 and ≤ 3.4 km at T21. The substellar shield volcano is one or two cells of 1–2 km. Blocking by the terminator range and orographic rain on single ranges are not resolved.
- **[inferred]** The filter also lowers the effective resolution of the flow itself to roughly two-thirds of the nominal truncation.
- **What a patch can do.** Decouple the orography filter from the physics filter, recovering up to ~900 m of peak height at T42 at the price of more ringing in rain and cloud. **What it cannot do.** Represent 8–10 km relief. That needs T85 or higher (§5). An aquaplanet benchmark cannot validate any topography change.

### 2.5 Energy balance

- **[measured]** From the vault's saved diagnostics, TOA net minus surface net:

  | run | grid | TOA − surface (W/m²) | TOA net | sea ice (ocean fraction) |
  |---|---|---|---|---|
  | t42_rc067, y70 | T42 | −0.06 | −0.19 | 0.64 |
  | t42_rc070, y65 | T42 | −0.16 | −0.08 | 0.63 |
  | f1800_rc070, y50 | T21 | −0.40 | −0.05 | 0.64 |
  | f1800_nir3x, y40 | T21 | −0.52 | −0.35 | 0.66 |
  | base, y15 | T21 | −0.84 | +1.06 | 0.00 |
  | nir_3x, y30 | T21 | −0.89 | −5.38 | 0.70 |
  | nir_2x, y30 | T21 | −1.92 | −10.22 | 0.55 |
  | clouds_hard, y20 | T21 | −2.42 | −5.05 | 0.00 |

- The atmosphere closes its budget to about 0.1 W/m² at T42 and 0.4–0.9 W/m² at T21 in settled states. There is no multi-W/m² leak.
- **[inferred]** Runs that "do not balance" at similar insolation are storing or releasing real energy, not leaking it. Two stores are slow:
  1. Sea-ice latent heat. A −5.4 W/m² global imbalance with half the planet ice-covered is about 1.1 m/yr of ice growth; reaching the 9 m cap takes most of a decade after the ice edge stops moving.
  2. The ice-albedo branch. The vault's own bracket shows a cliff between 1× and 2× near-infrared cloud optics. Runs near it drift for decades.
- **[code]** The wrapper's balance test (`_isbalanced`, `__init__.py:975`) looks only at the slope of |10-year mean| of TOA and surface net flux, needs at least 50–75 years, and never looks at temperature, ice thickness or snow. It can pass a run whose sea ice is still thickening, and it cannot say why a run fails.

### 2.6 Land surface, boundary layer, clouds

- **[code]** Land is a single bucket (`wsmax`), one snow layer, uniform ground albedo, roughness 2 m. Adjustable, not improvable without a new scheme. Desert-river results depend on it.
- **[literature]** Nightside surface temperature depends on stable-boundary-layer mixing to ±20 K across models (Boutle et al. 2017, A&A). ExoPlaSim's ten levels cannot resolve the inversion. This is a model limit, and it matters little for the weighted outputs.
- Cloud cover and optics are Earth-tuned diagnostics and are being used as free parameters. That is a legitimate use; it is not a prediction.

### 2.7 Wrapper and run tooling

- **[measured]** `sysconfigure()` fails to build `pyfft`/`pyfft991` under numpy 2 (return code 1); a manual `numpy.f2py -c` works. Already in the vault's notes; still unfixed here.
- **[measured]** Two `Model()` objects created at the same time for a not-yet-built executable race on compiling in the shared package directory and one of them crashes.
- **[code]** Restarts are bound to resolution; there is no supported way to carry a T21 state to T42. The vault built its own.
- **[code]** No per-year convergence log beyond TOA and surface means.

## 3. Ranked list

Score = impact on the weighted outputs × tractability. "Limit" means no patch of reasonable size fixes it.

| rank | problem | impact | tractability | verdict |
|---|---|---|---|---|
| 1 | Ocean heat transport: fix, expose and diagnose the existing slab diffusion | high: ice edge, terminator and dayside SST, nightside warmth; the one thing K60 names as missing | high: about 100 lines, default off | **patchable, partially.** Diffusion is a stand-in, not ocean dynamics |
| 2 | Convergence and energy-budget tooling | medium–high: decides whether any run can be called canon | high: Python only | **patchable** |
| 3 | Spectrum-aware near-infrared gas absorption | high: dayside surface energy, evaporation, rain, and it undoes part of the cloud over-tuning | medium: needs published band data and careful benchmarking | **patchable**, largest physics change of the set |
| 4 | Hydrology artefact: capped-snow runoff is ice export | medium: nightside and terminator discharge | high: a diagnostic flag | **patchable** as a diagnostic; the missing ice flow is a limit |
| 5 | Orography filter decoupled from physics filter | medium: orographic rain, blocking | medium to write, **poor to validate** | marginal gain; resolution is the **limit** |
| 6 | Sea-ice dynamics, calving, ocean mass budget | canon-relevant, low for weighted outputs | none | **limit** |
| 7 | Land bucket, nightside boundary layer, cloud scheme | medium (rain over land), low, tuned | none | **limit** |
| 8 | Longwave scheme at 1 bar / 1,000 ppm | low | — | no action needed |

### Rationale and failure modes

**1. Ocean heat transport.**
- For: it is the gap the vault already identified; most of the code exists; it is cheap; it has a clear published comparison (slab versus dynamic ocean).
- Against: K is a free parameter, so this adds a dial rather than removing one. It may move the ice edge without touching the western-heartland rain, in which case K60's hard choice is unchanged. Isotropic diffusion can give the right ice area with the wrong shape.
- Validation: THAI Hab 1 configuration, slab without transport against the four-model range (global mean 232–246 K, ice-free ocean 0.20–0.23), then with transport against the published dynamic-ocean changes for Proxima b and TRAPPIST-1e-like planets.

**2. Convergence tooling.**
- For: no physics risk; directly answers "why do some runs balance"; needed for the K60 acceptance test whatever else is chosen.
- Against: it improves confidence, not the climate.
- Validation: show on the benchmark run that the budget closes when storage terms are included.

**3. Near-infrared absorption.**
- For: it is a first-order error for any M-dwarf planet, and upstream would want it.
- Against: it changes the model's climate for every non-solar star, so it must be a flag. Coefficients must come from a published band model, which I have not yet sourced or verified. After it, the cloud tuning and probably the 1,800 W/m² choice have to be redone. The four benchmark models disagree by 14 K, so validation can show "inside the published range", not "correct".
- Validation: THAI Hab 1 surface shortwave and cloud radiative effect per model. I have those numbers only second-hand and unverified.

**4–5.** Small. Item 5 cannot meet the validation rule on a published aquaplanet case, so I would not implement it under that rule.

**6–7. Limits.** Stated in §2.3, §2.4, §2.6.

### Canon parameters

No canon parameter looks like the cause of a model failure. Two things to know:
- The cloud settings and the move to 1,800 W/m² were chosen against a model with the §2.2 error. They are tuning, not canon, but they would shift if item 3 is implemented.
- `mldepth` is un-pinned in the vault because it "only sets spin-up time". Once ocean diffusion is on, that stops being true: transport scales with depth × K.

## 4. Measurements made for this audit

Configuration: THAI Hab 1 values (900 W/m², 2,600 K, 6.1 d, 0.91 R⊕, 9.12 m/s², 1 bar N₂, 400 ppm CO₂, aquaplanet, 50 m slab), T21 L10, 30-minute step, `gp|exp|sp` filter. Driver: `audit/scripts/thai_hab1.py`. Differences from the protocol: blackbody star instead of BT-Settl, ExoPlaSim's own ice and ocean albedos instead of 0.25 / 0.06.

### 4.1 Does the existing ocean diffusion do anything?

Two runs from the same year-1 restart, 2 MPI ranks each, differing only in `NHDIFF@oceanmod_namelist=1`, `HDIFFK@oceanmod_namelist=50000.` (nominal; 4.1 × 10⁴ m²/s effective on this radius because of the Earth-radius error). Means of model years 9–11 (the diffusion run stopped at year 11 when the session's background job ended; the control ran to year 13).

| | control | ocean diffusion on | THAI Hab 1, four GCMs |
|---|---|---|---|
| global mean T_s (K) | 231.1 | 233.2 | 232–246 |
| dayside / nightside mean (K) | 257.2 / 203.4 | 259.3 / 205.5 | not tabulated |
| warmest / coldest cell (K) | 295.8 / 141.6 | 291.4 / 148.6 | UM minimum ≈ 125, others 20–30 K higher |
| ice-free ocean fraction | 0.22 | 0.29–0.30 | 0.20–0.23 |
| mean sea-ice thickness (m), still rising | 6.7 | 6.2 | — |
| TOA net / surface net (W/m²) | −0.6 / −0.4 | −0.7 / −0.4 | — |
| global precipitation (mm/yr) | 494 | 502 | — |
| planetary albedo | 0.31 | 0.30 | 0.22–0.28 (unverified table) |

Reading:
- The unmodified slab reproduces the benchmark's ice-free fraction and sits at the cold edge of its temperature range, with a higher albedo than any of the four models. That is a usable baseline.
- The existing diffusion works and conserves energy to the precision visible here (TOA and surface net are unchanged). At this K it opens a third more ocean, cools the substellar maximum by 4 K, and warms the coldest cell by 7 K. The direction matches the dynamic-ocean literature; the size is whatever K makes it.
- The control's sea ice is still thickening by ~0.05 m/yr with TOA at −0.6 W/m². This is the storage term of §2.5 seen directly.
- Ten years is short. These are differences between branches, not equilibria.

## 5. Compute

Measured: T21 takes 2.2 minutes per model year on 4 cores here; the vault measured 5–6 minutes (T21) and 17–21 minutes (T42) on the home PC.

- Items 1–4 need only T21 benchmark pairs of 30–40 years: one to four hours each on the home PC.
- A T42 confirmation of any change is 15–20 hours.
- Out of reach on a home PC: resolving 8–10 km relief (T85 is roughly 8–10 times T42 per year, days per decade) and a dynamic ocean (a different model).

## 6. Provisional run configuration for Ki-Engir

Provisional because the top-ranked change is not implemented. Values not listed stay as in the vault's `run_branch.py`.

| setting | T21 | T42 | note |
|---|---|---|---|
| `timestep` | 30 min | 30 min; drop to 20 if it crashes near the terminator range | |
| `physicsfilter` | `gp|exp|sp` | `gp|exp|sp` | keep; it costs ~900 m of peak height at T42 but removing it is unvalidated |
| `th2oc`, `gamma` | defaults | pin to 0.024 and 0.01 | as the vault does, so resolution is the only change |
| cloud and flux settings | vault's `f1800_rc070` | vault's `t42_rc067` | tuning, to be redone if item 3 lands |
| `mldepth` | 50 m | 50 m | pin it before using ocean diffusion |
| ocean diffusion | off until item 1 is done | off | the Earth-radius error makes current values non-portable |
| `glaciers` | off, ice surface in the topography | same | no ice flow in the model |
| `maxsnow` | default | default | but mask `mrro` where snow is at the cap |
| start | cold, or branch | transplant from T21 | as the vault does |
| stop | TOA and surface net within ±0.3 W/m² of each other and of zero over 10 years, global T_s drift < 0.05 K/yr, mean sea-ice thickness drift < 1%/decade | same | the wrapper's own test is not sufficient (§2.5) |
| averaging | last 10 years | last 10 years | |

## 7. References (from memory; not machine-verified unless stated)

- Boutle et al. 2017, A&A 601, A120. Proxima b with the UM.
- Del Genio et al. 2019, Astrobiology 19, 99. Proxima b with a dynamic ocean.
- Eager-Nash et al. 2020, A&A 639, A99. Stellar spectra and tidally locked climates.
- Fauchez et al. 2020, GMD 13, 707. THAI protocol.
- Hu & Yang 2014, PNAS 111, 629. Ocean heat transport on tidally locked planets.
- Lacis & Hansen 1974, J. Atmos. Sci. 31, 118.
- Paradise et al. 2022, MNRAS 511, 3272. ExoPlaSim.
- Sasamori 1968, J. Appl. Meteorol. 7, 721.
- Sergeev et al. 2022, Planet. Sci. J. 3, 212, doi:10.3847/PSJ/ac6cf2. THAI moist cases. Bibliographic details and the two quoted ranges were read from the arXiv page during this audit.
- Shields et al. 2013, Astrobiology 13, 715.
- Yang, Liu, Hu & Abbot 2014, ApJL 796, L22. Water trapping.
- Yang et al. 2019, ApJ 871, 29. Ocean dynamics and the inner edge.
