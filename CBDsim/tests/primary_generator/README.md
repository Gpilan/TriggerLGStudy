# Beam request validation

Build this CMake project after `source envset.sh`. `primary_generator_test` uses
the production geometry/generator and checks /gun position, normalized direction,
saved kinematics, both rotated tile intersections, and nonaccumulating x/y spreads.

Prefer `/gun/position x y z mm` and `/gun/direction dx dy dz`. Legacy y0/z0 and
theta/phi immediately update the gun; the last command wins. `randx`/`randy` mean
world x/y full widths (old y/z mapping was misleading and is no longer used).
Past position scans are not validated by this new implementation.
