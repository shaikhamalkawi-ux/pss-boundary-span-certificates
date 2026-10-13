# Claims that must not cross from this engineering example into the paper

1. The MovingAI city/street maps are discretized benchmark rasters of roads/buildings, not verified as-built engineering floor plans or measured obstacle clearances.
2. Native benchmark cell side length is not established as a physically calibrated metre value. All `r` values are **cells**.
3. The closed-square union connects on a single common vertex. That may not allow actual traversal, even for an idealized point robot under collision-free motion. The 4-neighbour edge-only component count is reported separately.
4. The full source case computes direct cubical homology. The exact triangle-cover PSS proof is built **only on bounded 8x8 windows**, not all 256x256 cells.
5. Generating an offline smaller proof does not imply fewer discovery queries or lower wall-clock cost. The generator evaluates all pair and candidate-triple statuses.
6. Same-source clearance perturbations and nine chosen V23 windows are dependent; no statistical confidence or claim of representative urban engineering scenarios follows.
7. The 2x2 conservative adapter retains only all-passable source blocks; it may remove physical paths and alter topology. The encoded admitted union is exact, not the original full-world union.
8. The official upstream MovingAI archive has not been byte-compared against the mirrored source maps, although the locally preserved V23 SHA-256 locks are verified.
9. No full urban navigation plan, robot clearance, industrial operating improvement, material fatigue, emergency response or safety claim is admitted from these outputs.
10. Real engineering transfer requires: source data license, survey-calibrated spatial resolution, physical constraints, motion model, start-goal tasks, operating condition variation, validation against direct route planning, and human engineering sign-off.
