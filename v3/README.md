# v3: Reef Wave 3D

A 3D breaking wave where the water motion comes from a physics solver, not a hand-drawn shape. Two presets: a Teahupo'o-style reef (thick lip, trough below sea level) and a Cloudbreak-style reef (thinner lip, rounder tube).

## Run it

Open `index.html` in a browser. It loads three.js r128 from cdnjs and the two `data_*.js` files next to it.

Drag to orbit, pinch or scroll to zoom. Camera buttons: Channel, Beach, Aerial, Tube, Lineup. Settings: break, wave size, peel speed, playback speed, follow the peel or stay fixed on the reef, cross-section inset, spray.

## How it works

1. **Solver** (`sim/solver.py`, Python + numpy/scipy, run offline). A 2D cross-section of the water: inviscid, irrotational, fully nonlinear free surface. Boundary-integral method in dipole / vortex-sheet form (after Longuet-Higgins & Cokelet, Baker-Meiron-Orszag). Surface markers move with the water, Bernoulli updates the potential, RK4 in time. The domain is periodic (440 m) and the reef is put in with a conformal map from a flat strip, so the sea-bed condition is exact. Markers are redistributed every 4 steps, down to about 6 cm on the face, lip and tube. The run stops when the lip touches the trough.
2. **Mesh** (`sim/process.py`). A wave meeting a straight reef at an angle is the same 2D event delayed along the reef, so solver time becomes distance along the reef. Mesh columns are water particles, mesh rows are time steps. Normals use the solver's particle velocities. Lighting helpers (how much sky each point sees, water thickness toward the sun) are ray-cast in the cross-section.
3. **Viewer** (`index.html`). Fresnel sky reflection, refraction to a reef bottom with absorption, sunlight through thin water for the glowing lip, surface texture carried by the water particles.

## What is not physics

- Everything after the lip lands. Potential flow ends at touchdown; the whitewater bore, foam and spray are a model.
- Each slice along the reef is treated as independent (fine for small angles).
- No wind, surface tension or air in the tube.
- The solver runs in advance; the page plays its output back.

## Reproduce the data

```
cd sim
python3 runcfg.py teahupoo   '{"Nmax":900,"ds_min":0.025,"ds_face":0.06,"N0":256,"cfl_g":0.18,"cfl_s":0.1,"dto":0.1,"dto_fine":0.02,"a":2.4,"x1":200,"hr":0.9}'
python3 runcfg.py cloudbreak '{"Nmax":900,"ds_min":0.025,"ds_face":0.06,"N0":256,"cfl_g":0.18,"cfl_s":0.1,"dto":0.1,"dto_fine":0.02,"a":2.0,"hr":1.0}'
python3 process.py teahupoo
python3 process.py cloudbreak
```

Each solver run takes about 80 s and ends with "NaN, stop" at touchdown, which is expected. `a` is swell amplitude (m), `x1` the start of the reef ramp (it ends at x = 285 m), `hr` the reef depth (m). `plot_profiles.py NAME.npz 14 8` draws the breaking sequence.

## Checks

`sim/tests/test_exact_fields.py`: velocities match exact harmonic test fields to 1e-9 or better, with and without the reef map.
`sim/tests/test_linear_wave.py`: a small wave over a flat bed has a phase error of 5e-4 rad after two periods; energy is conserved over a reef.
In the production runs total energy stays within about 3e-5 of its starting value through to touchdown.

## Known weak spots

- Whitewater is a smooth mound plus sprites, with foam "pillars" where the lip lands.
- The tube is lower and wider than photos of Teahupo'o. Very abrupt reef steps make the wave collapse at the toe instead of throwing, so the reef shape needs more tuning.
- Not yet checked on a real phone GPU (resolution drops automatically if frames are slow).
