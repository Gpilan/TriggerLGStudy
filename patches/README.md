# Study-specific patches

`2026-09-23-tape-entry-study.patch` records the additional boundary recovery used for the final noLG W10,s=−25mm production. It is **not applied to the shared source** in this checkpoint.

It extends only the named fully absorbing gel-tape entry test: when4 tolerances is still on the surface, accept two consecutive strict-inside probes among8,16,32,64 tolerances; any outside probe rejects recovery. The original guards and standard Geant4 absorption remain. Geometry, physical step length/time and material properties are unchanged.

Validation previously completed: exact original failing event payload reproduced before the repair; repaired event has no abnormal fates; paired64 local photons change10 NoRINDEX to tape absorption, with54 other records unchanged; normal200-photon budgets and RNG state identical. Existing24 world exits in the extreme corner fixture remain unchanged. The final3000-event production passes strict validation.

Inspect with `git apply --check patches/2026-09-23-tape-entry-study.patch`. Applying or promoting it is a separate model change and requires the relevant regression validation. Do not apply it blindly to an already-modified header.
