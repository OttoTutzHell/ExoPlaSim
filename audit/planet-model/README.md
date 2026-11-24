# Ki-Engir Planet Model: updated source

Source for the republished artifact. The original lives in the vault at
`tools/kiengir-planet-model/`; the vault was read-only for this work, so the changed files are kept here.

- `model.js`, `template.html`: the page. Star and orbit, spectrum-aware starlight, orbit sweep.
- `anchors.py`: summarises ExoPlaSim runs into `anchors.json` and exports two fields (`gcm_*.npz`).
- `prep.py VAULT kidata.js gcm_A.npz gcm_B.npz anchors.json`: builds the embedded data from the vault's relief.
- `build.py DIR`: assembles `ki-engir-planet-model.html`.
- `fit.js KEYS START ITERS`: coordinate search of the fitted parameters against `anchors.json`; result in `fit.json`.
