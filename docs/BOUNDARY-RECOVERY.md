# Optical boundary recovery — 2026-09-20

This change addresses rare production failures at packaging corners in Geant4
11.2.0, with a bounded 11.4.0 migration check recorded below. It does not change the intended tile, gel, window or LG dimensions,
material properties, reflectivity, QE, beam, or global geometry tolerance.

## Reproduced failures

- LG40: the inlet gel and a nested-Boolean Al corner rim repeatedly relocated
  between one another with zero steps, ending in `GeomNav1002`.
- noLG40: a gel/window corner relocation entered the tape with a zero step. The
  boundary process returned `StepTooSmall` before applying the tape surface;
  subsequent travel through PVC ended in `NoRINDEX`.
- noLG15: after reflection from one foil face, a second, perpendicular foil face
  was only 5.54e-10 mm away. Its boundary was skipped as `StepTooSmall`. Applying
  reflection alone was insufficient because the Boolean corner supplied an
  averaged normal that still pointed the reflected photon into the foil.

These are path-dependent edge cases, not evidence that 40 mm uniquely fails.
The other widths have different trajectories and encounter these edges with
different random samples. A separate old noLG30 failure also remains recorded.

## Geometry representation

The LG Al inlet rim is represented by four touching, non-overlapping boxes.
The noLG five-face foil shell is represented by nine non-overlapping boxes,
including the original 1 µm bottom notch around the tape. The occupied volume,
material and existing Al skin surface are preserved. These changes remove the
nested Boolean navigation ambiguity and provide the face normal for the
reproduced noLG corner. The earlier finite-facet intersection guard is retained.

## Narrow boundary-process adaptation

`CBDsimBoundaryProcess` derives from `G4OpBoundaryProcess`. The installer replaces
only the optical boundary process and checks that its effective GPIL/DoIt
ordering remains unchanged. The compile-time guard permits only Geant4 11.2.0 and 11.4.0.

The standard process handles ordinary steps. Recovery requires all of:

1. A geometry boundary with a step no greater than the existing surface tolerance.
2. An incident medium with RINDEX and a post-material without RINDEX. The last
   physical incident medium can be retained across a local zero-step relocation
   chain, including a transient, incorrectly reported pre-volume inside Al.
3. The existing `AluminumSurf` or `protoGelTapeAbsorber` surface, configured as
   unified/polished/dielectric-metal with zero transmission and efficiency.
4. Forward probes at 4 and 8 surface tolerances strictly inside that post-solid,
   and a valid navigator normal. Ordinary reflection-return relocations fail
   this entry predicate and are not processed twice.

For a qualifying entry, the original Geant4 surface algorithm is called once
using persistent process-owned track/step copies. Only the copy's track-length
dispatch marker is set to 16 tolerances, above both Geant4 short-step guards.
The actual step length, coordinates, times, track and global tolerance remain
unchanged. The copied G4Step preserves the true path length used to initialize
the particle change; this is checked at runtime. The marker is not a physical
displacement and does not add propagation time. Proxies remain alive until the
particle change has been consumed, and their touchable handles are released at
EndTracking before worker geometry/allocator teardown.

Existing Geant4 reflection, polarization, absorption and random draws determine
the outcome. If a reflective normal would return a photon into the opaque bulk,
the run fails explicitly instead of silently accepting the path. No fake bulk
RINDEX, forced world kill, additional absorber, or relaxed validation gate is
introduced. A compact log line records each recovered boundary; this is not
full photon step logging.

## Validation scope

- Exact selected-event replays from the original fixed seeds: LG40 event133,
  noLG40 event84, noLG15 event155. Before repair, all persisted target ROOT members
  matched their original production event. After repair all three exit normally
  with no navigation warning, NoRINDEX, unknown/unfinished fate, world exit or
  photon-balance error in these events.
- LG known-corner navigation rays and 64 adversarial optical edge probes pass.
  The latter had 16 NoRINDEX terminations with the original boundary handling.
- One million occupied-volume samples each for the LG rim and noLG foil show
  no mismatch; noLG box intersections show no double occupancy.
- All three production conditions pass 20,000-sample placement-overlap checks,
  interface navigation checks, 200 optical coupling probes and two full electron
  events. The two-event canonical ROOT physics digests match the earlier normal
  controls exactly in each condition (120,995 photons total).
- Normal optical control runs with/without the adapter have identical diagnostic
  summaries and complete final RNG states, for 200 photons in each mode.
- A separate two-event actual Condor worker run matches the new local reference
  digest exactly before production is allowed to start.

This is a bounded workaround for the reproduced geometry/model and version,
not a general replacement for Geant4 boundary tracking. Full production remains
subject to strict geometry, NoRINDEX, photon-budget, ROOT and provenance checks.
Do not mix old chunks with reruns under this implementation.

## Upstream context

The Geant4 11.2 [boundary implementation](https://github.com/Geant4/geant4/blob/v11.2.0/source/processes/optical/src/G4OpBoundaryProcess.cc)
uses a short-step return before optical-surface handling. Maintainer discussions
describe the difficulty of distinguishing genuine close boundaries from
[relocation steps](https://geant4-forum.web.cern.ch/t/issues-with-optical-simulations-at-boundary/11823/6)
and the corresponding [corner NoRINDEX problem](https://geant4-forum.web.cern.ch/t/error-occurs-when-an-optical-photon-hits-the-edge-of-a-cubic-scintillator/8748/2).

## noLG30 rim-corner extension (2026-09-21)

The noLG30 production seed181830200/event118/track18167 exposed a distinct
long-step corner case. The scintillator exit is within tolerance of two faces;
the selected exit normal reflected the photon into the receiving Al rim. It
then traveled0.0341mm in Al and terminated as NoRINDEX. The prior short-step
adapter alone did not correct this event.

The noLG rim is now partitioned at the existing axis-aligned solid boundaries,
classified against the original unplaced Boolean reference, and merged along
y first. The normal configuration has14 boxes with the same material and skin
surface. A candidate that split only at y=-30mm was rejected: zero-length
reflections at that artificial seam eventually absorbed the photon.

The boundary adapter adds a restricted receiving-face correction for a regular
scintillator→noLG rim-box step only when:

- Both forward tests at4/8 surface tolerances are strictly inside the opaque box,
  and the point8 tolerances before the boundary is inside the incident scintillator.
- The supported polished, opaque surface and refractive incident medium pass
  the existing checks; the ordinary specular reflection points inside that box.
- Exactly one receiving box face lies within tolerance, all other faces are more
  than8 tolerances away, and that face's reflection exits the opaque box.

A scoped navigator facade supplies only the normal queried synchronously by
Geant4 11.2's unchanged optical boundary process. The active-navigator entry is
restored before return; the transportation tracking navigator is not replaced.
This path requires one mass-world navigator and no parallel-world hyperstep.
The actual step length, points, time, optical properties and standard random
reflection/absorption calculation remain unchanged. Ordinary paths allocate no
facade and invoke the existing boundary process.

Validation: exact original event ROOT payload equality before repair; repaired
event NoRINDEX/unknown/unfinished/world exit/balance error0, with track18167
reaching gel/window/wafer and being detected.64 nearby optical probes have0
NoRINDEX (standard boundary:24),0 loops and at most2 consecutive zero steps.
Normal200-photon summaries and full RNG states match the standard process.
One million sampled points comparing actual placed boxes against the Boolean
reference show no volume mismatch or double occupancy. noLG30 placement overlap,
interface navigation,200 coupling photons and2 production electron events pass;
normal electron ROOT digest equals the earlier control exactly.

This validation is scoped to these cases. Production still rejects navigation,
NoRINDEX, budget or provenance failures; it does not bypass them to fill the scan.
Submitted-job status is retained in the separately archived production campaign
`nolg30_retry_20260921`.

## Geant4 11.4 migration check (2026-10-06)

An isolated candidate changed only the version guard before promotion. The
current shared W40 geometry, materials and wafer QE were preserved; the separate
September 23 tape candidate was not imported. With Geant4 11.4.0 and seed 42:

- The standard process reproduced 16 NoRINDEX terminations in 64 LG inlet edge
  probes. The existing recovery algorithm passed all 64 with zero NoRINDEX,
  unknown/unfinished fates or loops, and a closed photon budget.
- Normal LG and noLG controls, 200 photons each, had byte-identical diagnostic
  summaries, path histograms and complete final RNG states with/without recovery.
- The original short-step dispatch and normal-recovery algorithm is unchanged.
  Its true-path-length runtime assertion and GPIL/DoIt ordering checks remain.

This checks the stated current-geometry probes. It does not revalidate every
historical width/length/thickness or the selected old production failure events.
An LG-specific edge probe was also attempted against noLG; its entry condition
is inapplicable there and it is not counted as a passing noLG edge test.
Version-to-version timing equivalence and production readiness remain separate
gates. See [build and validation notes](GEANT4-MIGRATION.md).
