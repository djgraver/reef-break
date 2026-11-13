# Study 1: Breaking Wave

A fresh start on the wave, focused only on how it looks and how it crashes (physics is not the goal). Single self-contained HTML/Canvas2D file, no dependencies. Open `index.html` in a browser.

## How it works

- **Sheet profile.** Each column across the wave is a rolled sheet of water: flat sea, one straight low ramp, then a curve that rolls up the face and over the top as the lip. The radius of curvature changes gradually (about phi^2 x the tightest radius at the start of the roll, tightest in the barrel, looser along the thrown lip), so there are no sharp bends. The ramp is built from a noise-free state, so it is identical along the whole wave; all irregularity lives above it.
- **Life of a section.** swell, pitch, barrel, impact, whitewater, fade. In peel mode one frame shows the whole sequence along the wave.
- **Software 3D.** Per-column skin polylines are meshed, shaded per quad (Fresnel sky reflection, sun glint, translucent glowing lip, foam) and drawn painter's-algorithm style. Pixel mode quantises to a ~50 colour ramp with an 8x8 Bayer dither.
- **Foam and spray** are stateless procedural sprites: foam clumps, the impact burst where the lip meets the rising water, ballistic spray, drips, spindrift.
- **Variation.** Multi-scale noise on height, lip thickness, hollowness, break speed and the crest line, deterministic per seed.
- **Model scale.** The Size slider scales the whole model (geometry, foam, spray), not the camera.
- **Flow.** Soft bands travel along the sheet so the water appears to run up the ramp, up the face and over the lip.
- **Sea.** Darker gradient matching the wave, with a chop layer of small ripples.

## Camera presets

Beach 3/4 (front-on from the shoulder side), beach front-on, down the line, and a camera inside the barrel.
