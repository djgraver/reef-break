# v4: Reef Wave 3D — one wave, a ride, a place

Same solver data as v3 (see `../v3/README.md` and `../v3/sim/` for the physics and how to regenerate it). The changes are all in the viewer, `index.html`.

## Run it

Open `index.html` in a browser. It loads three.js r128 from cdnjs and the two `data_*.js` files next to it (identical to the v3 ones).

## What changed from v3

- **Sections.** The lip throws slightly early or late along the reef (smooth noise fixed to the reef, about ±0.05 s at normal size). Only the crest and lip shift; the face stays put.
- **Size.** The slider runs from about 1.6 m to 12 m faces. As size goes up the sections get bigger (up to about ±0.2 s), crest height varies along the wave, the surface gets rougher, and there is more spray and whitewater turbulence.
- **Start and end.** The wave breaks along a finite stretch of reef and the ride takes 10 s at any size. The swell arrives, peels from the start point, fades into a channel at the far end, and the cycle repeats. Replay button in Settings.
- **One wave.** The repeated waves behind are gone. A single low swell ribbon sits about 230 m further out to sea.
- **Land.** An exposed reef platform at the take-off, and a sand beach with vegetation inshore of it, with shallow water and a wash line around them.

## What is not physics

Everything listed in the v3 README, plus: the sections, the size-driven unruliness, and the start and end of the wave are shaping applied in the viewer, not solver output.

## Known weak spots

- The wave switches on beside the exposed reef by growing from flat rather than wrapping around it.
- The exposed reef is a plain round platform.
- The far swell is a short low hump, not a full swell line.
- Whitewater is still a smooth mound plus sprites.
- Not yet checked on a real phone GPU.
