# Reef Break

A 90s-style (EGA/VGA pixel look) breaking-wave simulator, the first step toward a surfing game inspired by California Games.

Waves are computed, not modelled by hand. Each wave is a crest line of 520 columns, and every column carries its own state (position, height, curl phase, crash timer). Columns move at sqrt(g * depth) over a sampled seabed depth grid, so the crest bends over the reef, shoals, and breaks once height over depth passes a breaker index that depends on bed slope. A broken neighbour makes the break run along the crest (the peel). The wave surface is a mesh of 520 x 31 computed profile points, shaded per pixel into a 30-colour palette at 480x270. Whitewater and spray use pre-generated sprites.

## Run it

Open index.html in a browser. It is a single self-contained file (it loads two web fonts from Google Fonts).

## Controls

Space: pause. S: send a wave. D: data view (crest lines over seabed depth). Buttons: speed, reef preset (peels left to right, right to left, A-frame, closeout sandbar), swell size.

## Status

Branch `v1` holds the first tuned version (deeper barrel, longer-lived whitewater). `main` is v2.

**v2** (current): a more realistic wave, built from real barrel photos.
- The tube interior glows green-turquoise where sunlight comes through the lip, over a deep navy face.
- A thicker, rounder lip, and a curl that scales up with how hard the wave plunges.
- Wind mist blows off the lip as soft, drifting, fading plumes; the face has fine chop and sun glints.
- The lip lands and splashes within a second, and broken water keeps rolling in as clumpy whitewater instead of vanishing.

## v3: Reef Wave 3D

`v3/` is a new start in real 3D. The wave shape comes from a 2D free-surface flow solver run over a reef (Teahupo'o-style and Cloudbreak-style presets), swept along the reef to make a peeling barrel, and shown with three.js. Whitewater after the lip lands is modelled, not solved. See `v3/README.md`.

## v4: one wave, a ride, a place

`v4/` is the v3 viewer with sections along the lip, a wider and wilder size range, a 10-second ride with a start and an end, a single wave with a distant swell line, and an exposed reef and beach at the take-off. Same solver data as v3. See `v4/README.md`.

**Next**: better whitewater and reef tuning, then gameplay (the surfer riding the line).

## Study 1 (restart, wave only)

`study-1/` is a fresh take on the wave alone, with a different approach: a rolled-sheet profile (flat sea, straight ramp, gradually tightening spiral up and over the lip), a bigger impact explosion, rolling whitewater tail, and flowing water. See `study-1/README.md`.
