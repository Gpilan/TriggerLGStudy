# Nominal composite sensor verification

This compiles production proto geometry, materials, SD and optical diagnostics.
The sensor is an idealized **conditional wafer-QE model**, not a validated R2076
photocathode. QE is applied once on the first wafer step, every accepted/rejected
track is killed, and time is recorded at the pre-step point (wafer entry).

```bash
source envset.sh
cmake -S CBDsim/tests/sensor_response -B /tmp/trigger-sensor-build
cmake --build /tmp/trigger-sensor-build -j2
# Run from repository root to resolve the QE CSV:
/tmp/trigger-sensor-build/sensor_response_test 3000 42 /tmp/sensor.txt 450 45 1 0
```

Arguments: photons, seed, diagnostic output, wavelength nm, angle degrees,
start inside wafer (1) or just inside front of window (0), noLG (1) or LG (0).
Polarization is perpendicular to the incidence plane. The test covers real
nominal dimensions for both geometries and counts wafer steps independently of
terminal fates. Window-start losses include refraction/reflection/absorption
and possible returns through the complete detector.

Conditional detection probability is CSV QE; it is not detection probability
per incident photon at the window. Source outside CSV range currently uses
endpoint clamping, an unvalidated extrapolation inherited from the old model.
Do not call the n=3.8 silicon layer a bialkali photocathode, or infer real PMT
collection efficiency from it. LG aperture is circular (176.715 mm²), noLG is
15×5 mm² rectangular (75 mm²); this is a configuration comparison.
