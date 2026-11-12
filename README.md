# Reef Break

A 90s-style (EGA/VGA pixel look) breaking-wave simulator, the first step toward a surfing game inspired by California Games.

Waves are computed, not modelled by hand. Each wave is a crest line of 520 columns, and every column carries its own state (position, height, curl phase, crash timer). Columns move at sqrt(g * depth) over a sampled seabed depth grid, so the crest bends over the reef, shoals, and breaks once height over depth passes a breaker index that depends on bed slope. A broken neighbour makes the break run along the crest (the peel). The wave surface is a mesh of 520 x 31 computed profile points, shaded per pixel into a 30-colour palette at 480x270. Whitewater and spray use pre-generated sprites.

## Run it

Open index.html in a browser. It is a single self-contained file (it loads two web fonts from Google Fonts).

## Controls

Space: pause. S: send a wave. D: data view (crest lines over seabed depth). Buttons: speed, reef preset (peels left to right, right to left, A-frame, closeout sandbar), swell size.

## Status

**v1** (current): fixed beach camera, four reef presets, a tuned barrel (plunging waves pivot lower and throw a thicker lip, so the tube is deeper and stays open longer), and longer-lasting whitewater.

**v2** (next): the surfer riding the line.
