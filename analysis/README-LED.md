# LED selection with CFD timing analysis

`LED.py` implements fixed-threshold event selection on the existing reconstructed waveforms. It does not replace CFD or change the SPE model.

- Absolute threshold: `k * SPE_peak`, with analytic peak-normalized SPE reference1 a.u. The reference is fixed across all events and geometries. Sampled one-photon peaks may be slightly below1 due to sampling and fractional-bin placement; do not renormalize each event.
- LED time: first upward crossing, linearly interpolated between adjacent samples. Empty signals and signals below threshold fail selection. Already-above-at-start waveforms have no known crossing.
- Optional timed coincidence: `abs(tLED1-tLED2-offset) <= window_ns`, where window is a **half-width**. `none` applies the two thresholds without a timing cut.
- CFD remains fraction0.3,10ps waveform sampling, local-fit half-width2 with the established rise1.3ns/FWHM3ns/transit14ns response. Both channel CFD times are computed before LED cuts and retained unchanged.
- All generated events remain in the incident denominator. Out-of-range photon bins are marked unknown for LED, not silently treated as a measured trigger failure. The reported selected fraction is then the known-selected fraction; unknown counts must be inspected.
- No primary-truth cut is used to form electronic coincidence. Primary entry categories are separate diagnostic fields when available.

Example (exploratory settings, not an approved operating point):

```sh
python3 analysis/LED.py /path/to/canonical-events.json --output /new/path/LED.json \
  --k 1 10 50 100 200 400 --window-ns none 0.5 1 2 --offset-ns 0
```

The output retains pre-window CFD values, per-event selection masks and selected CFD arrays for subsequent ROOT fitting. No LED Δt fit is substituted for the CFD resolution. Zero-signal events and low-statistics samples are not replaced. ROOT fitting/efficiency plots for full beam production will use these arrays after pilot calibration.
