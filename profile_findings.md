# Profile Analysis: `achilles_geom_driver` (2026-06-15)

Two new samply profiles (1 ms interval):

- `build/achilles_geom_driver 2026-06-15 18.18 profile.json.gz` — **Retry** mode (TotalXSecRetry), ~137.8 s, 137,786 samples.
- `build/achilles_geom_driver 2026-06-15 18.26 profile.json.gz` — **NoRetry** mode, ~924.9 s, 924,860 samples.

The earlier `Material` value-copy churn (previously the #1 hotspot, ~35–40%) is **gone**. The run is
now a compute-bound geometry-traversal problem. NoRetry is the cleanest signal: event generation is
only ~4%, and `World::GetLineSegments` is **88%** of wall time. Retry shares the identical kernel and
just adds generator cost on top (GetLineSegments ~52%, event gen ~26%), so every fix here helps both.

## NoRetry self-time, grouped

| Group | Functions (self %) | ~Self % |
|---|---|---|
| **AABB / BVH kernel** | `BoundingBox::IntersectT` (19.7), `BVH::Traverse` (9.1), `std::isfinite` (2.3), `vector<Node>::operator[]` (1.3) | **~32%** |
| **Per-iteration ray transforms** | `PhysicalVolume::TransformRay` (6.4), `Transform3D::Apply` (3.2), `ApplyRayDirect` (3.2), `Ray::Ray` ctor (3.4), `TranslateRay` (1.3), `Transform3D::operator*` (1.3) | **~19%** |
| **Shape intersection** | `Box::Intersect2Impl` (6.9), `Shape::Intersect2` (2.4), `CombinedShape::Intersect2Impl` (1.7), `PhysicalVolume::Intersect` (1.4) | **~12%** |

### Root structural inefficiency

In `PhysicalVolume::GetLineSegments` / `LogicalVolume::GetLineSegments` the ray is carried in the
**global frame** and re-transformed into the volume's local frame **on every loop iteration**
(`ApplyRayDirect` + `TransformRay`), and a fresh `Ray` is reconstructed at each step — each
reconstruction redoing 3 divisions for `inv_direction` and re-rotating the direction — even though
**only the origin advances along an otherwise-fixed ray**. Transforms are rigid (rotation+translation),
so the ray parameter `t` is preserved across frames; the local direction/inv-direction are invariant
within a volume.

---

## Plan (priority order)

**#1 — Branchless `BoundingBox::IntersectT`** (biggest single hotspot, ~20% + the 2.3% `isfinite`).
Replace the per-axis `std::isfinite(inv_d)` branch + `std::swap` with the robust min/max slab
formulation (PBRT / Tavian Barnes); IEEE handles infinite `inv_d` without branching, ordered to be
NaN-safe for the on-face zero-direction case. `Ray::InvDirection()` is already precomputed.

**#2 — Hoist ray transforms out of the traversal loop** (~19%).
Carry the ray in the volume's local frame and advance only the origin by `t` each iteration instead of
re-deriving the local ray from the global ray every step. Moves the transforms from per-iteration to
per-volume-entry.

**#3 — Stop recomputing `inv_direction` on every `Ray` reconstruction** (~3.4% ctor).
Add a cheap "advance origin / reuse direction+inv_direction" path. Falls out of #2.

**#4 — BVH leaf/node micro-optimizations** (~10%).
Cache `m_nodes.data()` raw pointer in `Traverse`; output a raw `PhysicalVolume*`/index instead of
copying a `shared_ptr` per hit; reuse the already-local ray for the leaf `Intersect`.

**#5 — Shape-intersect cleanup** (~12%, re-profile first).
`Box::Intersect2Impl` should reuse the optimized slab kernel from #1; consolidate redundant
entry+exit queries per step.

### Validation between steps
Re-run `[BVH]` regression tests (`test/test_bvh.cc`: adjacent-daughter / exact-fill / AABB-grazing
cases guard exactly these numerics) and re-profile to confirm the share moved as predicted. Behavior
must stay bit-stable — these are correctness-sensitive paths with a long bug-fix history.

### TODO afterwards — investigate "ray-march"
The traversal does a full root-to-leaf BVH descent **plus** a shape exit-test at every segment
boundary, so an N-segment ray pays N descents. Investigate replacing this with a single ordered
ray-march that collects all boundary crossings in one pass. Larger algorithmic win but a substantial
rewrite; evaluate only after #1–#4 are measured.

---

## Ray-march investigation (boundary-sweep vs Geant4-optimized sequential)

Geant4/ROOT both use *stateful sequential tracking* (ComputeStep/DistanceToOut+DistanceToIn,
voxelization = our BVH, safety cache, navigation-history incremental transforms). They deliberately
avoid an analytic all-boundaries solve because solids are non-convex (CSG). Our neutrinos fly straight
and we want the whole path up front, so we *can* try the analytic sweep — but this geometry has 336
subtractions, and our `CombinedShape::Intersect2` returns only the **first** [enter,exit] interval.

Built **both paths side by side**, runtime-selectable via `NUGEOM_TRAVERSAL=sweep` (default =
sequential): `World::GetLineSegmentsSequential` (the Geant4-optimized stepper) and
`World::GetLineSegmentsSweep` (one BVH `CollectHits` descent collecting all volume `[enter,exit]`
events via `PhysicalVolume::CollectIntervals`, then a depth-tracking sweep where the deepest active
volume wins). Cross-validated by the hidden `[.][benchmark]` test in `test/bench_segments.cc` (500 rays
through the real ND GDML, sequential as oracle) + the whole unit suite re-run with the sweep enabled.

**Result (2026-06-15):**
- Unit suite under `NUGEOM_TRAVERSAL=sweep`: all 272 [BVH] + 857 assertions PASS (sweep matches
  sequential on every synthetic geometry incl. adjacent/exact-fill/AABB-grazing/CSG cases).
- Real ND GDML, 500 rays: **sweep ~1.6-1.8x faster** (~20.6 ms vs ~36.6 ms; 24.3k vs 13.6k rays/s).
- BUT sweep mis-segments **29/500 rays** (27 axis-aligned, 2 angled). Example: seq `Vac Rock Air Rock
  Vac` vs sweep `Vac Rock Steel Rock Air Rock Vac` -- the sweep inserts spurious Steel/Rock where the
  ray is actually in a hole. Root cause = single-interval `Intersect2` on non-convex CSG (and/or
  overlapping same-depth placements): a subtraction's outer interval marks the volume active across its
  hole. **Confirms the sweep needs `Shape::IntersectAll` (true crossing lists, CSG = boolean merge of
  child crossings) to be correct** -- exactly why Geant4/ROOT step sequentially.

**Decision input:** sweep has real ~1.7x upside but is not yet correct on this geometry; the
sequential-cached path (committed 5d20fe2 + 051eeb7) is correct now.

**Mismatch diagnosed (2026-06-15) — it is NOT a CSG Intersect2 bug.** Instrumented CollectIntervals
(env `NUGEOM_SWEEP_DEBUG`) on the first mismatching ray (#5): spurious `Steel` came from
`elevatorBlock_lv` (depth 3) with interval t=[30988.8,31019.7], SDF@mid=-15cm (genuinely INSIDE its
raw Steel solid), yet its mother `volDetEnclosure` (depth 2) spans only t=[31215.6,32241], and
`FindMaterial` (point-in-volume oracle) returns Rock there. Across all 536 collected intervals, **0
were "OUTSIDE"** their raw solid -> `CombinedShape::Intersect2` is correct. Root cause:
`elevatorBlock_lv` **protrudes outside its mother** (or GDML overlap). The geometry's material
semantics are *hierarchical* (deepest volume contained at every ancestor level = what sequential &
FindMaterial do); the sweep takes "deepest raw solid active," so a daughter poking outside its mother
leaks. **Fix = hierarchical interval clipping (child intervals ∩ parent intervals), NOT IntersectAll**
(intervals are already complete here). IntersectAll only needed later if a CSG *hole* (multi-interval
non-convex mother) appears after clipping. Next: implement clipping in CollectIntervals, re-validate to
500/500, re-benchmark.

**Clipping implemented + re-diagnosed (2026-06-15).** CollectIntervals now clips each volume to its
parent's window (child intervals ∩ parent). Result: 29 -> 25 mismatches (fixed the protrusion class
e.g. elevatorBlock); sweep still ~1.65x faster (22 vs 36 ms). **Residual ~25 (5%) are a deeper model
mismatch, NOT CSG and NOT fixable by IntersectAll:** e.g. ray 36 misses two 0.2 cm SSteel304 cryostat
membrane layers. The sweep's deep stack is volUllage(GAr)->volTopInsulation(GAr)->volTubs(GAr)->
volTub(Polyurethane d15) so depth-override picks Polyurethane, but the steel wall belongs to a
DIFFERENT tree branch whose box overlaps the region (volTub raw box [402.4,929.6] overlaps the steel at
[403.39,403.59]). Root cause = **non-strictly-nested / overlapping volumes** (common in DUNE GDMLs):
the sweep's "deepest TREE-DEPTH wins" != the sequential's exact step-by-step containment (resolves
overlaps by traversal order). The missed features are thin cryostat steel = real interaction targets.
**Conclusion: depth-override sweep is an approximate fast path (~1.65x) that is NOT bit-exact on ~5% of
rays for this overlapping geometry; exact parity would need per-region containment resolution, eroding
the speed advantage.** Decision pending (accept approximate optional fast path / invest in exact
containment / shelve sweep & keep committed sequential-cached).

**CORRECTION (2026-06-16, external ROOT check) — the geometry is NOT overlapping.**
`tools/check_overlaps_root.C` (ROOT 6.38 TGeoManager, `/usr/bin/root.exe`, NOT the broken ~/.local
python-alpha ROOT) imports the same GDML (194447 nodes) and reports **ZERO overlaps** at 1e-4 cm (mesh
+ sampling). So the sweep's disagreements are **our sweep's own approximation, not GDML overlaps**:
an incomplete single-interval `Intersect2` on a *non-convex (CSG) ancestor* over-clips a legitimately
nested deep child (e.g. volTPCActive[LAr] gets clipped out under a CSG ancestor, leaving the shallower
volHalfDetector[G10] as deepest-active). Implications: (1) the [overlapscan] is a *sweep-accuracy*
probe, not a geometry-overlap detector; (2) **IntersectAll (complete CSG crossing lists) IS the fix
after all** — it would give ancestors their full multi-interval extent so clipping stops dropping
nested children (earlier dismissal was based on the pre-clipping protrusion case, since fixed by
clipping); (3) still TODO: validate our SEQUENTIAL (production) path against ROOT's TGeoNavigator on
identical rays to confirm our parser + sequential are correct and only the sweep is approximate.
Volume-named scan signatures + the ROOT script committed (5ef10e6).

## ✅ PRODUCTION BUG (found by ROOT comparison) — FIXED 2026-06-16 (commit 2de1e94)

Root cause: `daughter_fg = from_global * m_transform` built each child's world→local transform in the
**wrong operand order** — `(A*B).Apply(p) == A.Apply(B.Apply(p))`, so it applied m_transform (parent→
local) before from_global (world→parent). Correct = `m_transform * from_global`. Translations commute,
so it was latent until the rotated exact-fill chain `volHalfDetector_R → volFieldcage → volTPCActive`,
where the rotated daughter's frame was off by ~850 cm (NOT an eps/boundary issue — the first guess),
the descent skipped it, and ~5 m of LAr was reported as the G10 shell. Fix = swap operands at both
daughter_fg sites in `src/geom/Volume.cc` (sequential + sweep). Ray 7 LAr now 523.43 cm = ROOT (was
45). Regression tests in `test/test_bvh.cc` ([BVH]): "Transform3D composition order" and "GetLineSegments
TPC traversal matches containment" (loads the ND GDML, checks vs FindMaterial oracle, requires metres of
LAr); both fail on the pre-fix order. A synthetic rotated nesting did NOT reproduce it — the real
exact-fill/dt==0 chain was required.

### Second bug, also found via ROOT — FIXED 2026-06-16 (commit 00ef333)
`CombinedShape::Intersect2` walked ALL child boundary crossings, including ones behind the ray origin.
For a non-convex two-wall CSG (window frame = box minus central slab) the ray, once past the first
wall, has the earlier crossings at negative t; the function returned that behind-origin interval and
never reported the wall ahead. This dropped the exit-side composite-window steel (volCompositeWindow),
so GetLineSegments reported ~1 cm SSteel304 as Air. Fix: skip non-positive crossing times (origin
in/out state is already in inside(0)). **After both fixes ALL 12/12 navigator-comparison rays match
ROOT's TGeoNavigator exactly.** Regression test: `test/test_shape.cc` "CombinedShape::Intersect2 returns
the forward wall (two-wall subtraction)".

### IntersectAll + final re-comparison (2026-06-16, commit 8060f37)
Added `Shape::IntersectAll` (all solid spans; CombinedShape = interval-set union/intersect/subtract of
children; convex = single Intersect2 span). `CollectIntervals` now emits an event pair per span, so the
sweep captures multi-interval CSG (the two-wall composite window) instead of just the first wall.
**Result: sweep vs sequential cross-validation 400/500 -> 500/500; navcompare 12/12 vs ROOT; the
[navmismatch] test now writes 0 mismatch rays.** Timing: **sweep ~1.9x faster** (~23.6 ms vs ~44.6 ms /
500 rays). Sweep-vs-containment overlap scan: 12 -> ~7/2000, the residual being the parser hall bug
(elevatorBlock) + a couple thin TPC-PCB depth-override edges that the sequential shares (sweep ≡
sequential now). **Net: both traversals agree and match ROOT; the sweep is a correct ~1.9x-faster
alternative.** Only open item: the parser nested-CSG hall shape (volDetEnclosure) -- a geometry-
construction bug, not traversal, lower priority.

### Sequential vs sweep re-comparison (2026-06-16, after the daughter_fg + Intersect2 fixes)
Both fixes touch `Intersect2`/`daughter_fg`, which the sweep's `CollectIntervals` also uses, so the
sweep improved a lot: the `[overlapscan]` (sweep vs FindMaterial containment, 2000 rays) dropped from
**154 -> 12** disagreements. Timing unchanged: **sweep ~1.9-2x faster** (~22 ms vs ~45 ms / 500 rays).
The `[benchmark]` strict cross-check (sweep vs the now-ROOT-correct sequential) is 400/500; the
remaining mismatches are almost entirely the **two-wall composite window**: `CollectIntervals` does a
*single-interval* `Intersect2` per volume, so on a multi-interval CSG it captures only the first wall
(seq: `...Air SSteel304 Air...`, sweep: `...Air...`). **The sweep's last correctness gap is multi-
interval CSG volumes -> needs `IntersectAll` (all crossings) in CollectIntervals.** Sequential is exact.

### ⚠ THIRD bug found (open) — angled rays into CSG volDetEnclosure (2026-06-16)
Took the 100 rays where sweep != sequential and did a 3-way compare vs ROOT (`[.][navmismatch]` test +
`navigate_rays`):
- **Sweep: all 100 mismatches are a pure subsequence of ROOT** (it only ever *drops* walls, never
  reports a wrong material) -> confirms the sweep's only gap is the single-interval CSG wall miss
  (needs IntersectAll), NOT a new bug.
- **Sequential: 19/100 also disagree with ROOT** — DIAGNOSED 2026-06-16 as a **GEOMETRY-CONSTRUCTION
  (parser) bug, NOT a traversal bug.** With z computed correctly for the angled ray (×dir_z), our
  `FindMaterial` (point-in-volume containment) returns Rock through z≈863 and Air at z≈950, i.e. our
  parsed `volDetEnclosure` (Air) starts ~475 cm LATER than ROOT's (z≈388). The sequential traversal
  MATCHES our FindMaterial (~863), so the traversal is correct — our *parsed shape* differs from ROOT's.
  `volDetEnclosure`'s solid `NDHallAirVolShape` is a **6-level left-nested union of boxes** with large
  positioned operands (e.g. NDHallAirVolSpace7 at x=-2758). `GDMLParser::ParseSolids` builds each
  boolean's second operand as a `TransformedShape` (ParseCSGTransform = `(rot*Translation3D(pos))
  .Inverse()`); the chained-union result lands ~475 cm off ROOT's TGeoManager. **Open parser/nested-CSG
  positioning issue, hall Air volume only** (the detector region parses correctly — the 12 navcompare
  rays match). Lower physics priority (hall/enclosure edges, not interaction targets). Separate effort:
  verify how chained boolean operand frames compose vs ROOT. Repro: `[.][navmismatch]`.

### Original report (for context)
## ⚠ PRODUCTION BUG found by ROOT TGeoNavigator comparison (2026-06-16)

External validation against ROOT's `TGeoNavigator` (`tools/check_overlaps_root.C` `navigate_rays`,
diffed by `tools/compare_segments.py`, rays from the `[navcompare]` test) found a **correctness bug in
the production sequential `World::GetLineSegments` traversal**, independent of the perf/sweep work:

- Total path length matches ROOT to <1e-4 cm on all 12 rays; bulk materials (Vac/Rock/Air/Steel/
  Polyurethane/Scintillator/SteelTMS/Nitrogen) match exactly. So the **parser/units/transforms and
  bulk traversal are correct**.
- BUT in the **TPC** the sequential traversal emits a pathological **`99 cm G10 / 1 cm LAr` pattern
  (~5 m)** where ROOT sees continuous LAr with thin 0.635 cm G10 boards. Affected ray positions only
  (e.g. ray 7 = x@0.65 of ArgonCube). 6/12 sample rays differ, all in this region.
- `[nav7probe]` proves it: inside each bogus 99 cm G10 span, **`FindMaterial`/`FindVolume` (our own
  point-in-volume containment) returns `LAr`/`volTPCActive`** — so the geometry is right and the
  **ray traversal is mis-assigning the LAr active target to G10**.
- Both ray-traversal paths share the bug: the sweep's overlap scan flagged the same `G10 (sweep) vs
  LAr (contained)` 600x. So `FindMaterial` (point containment) is the only one of our methods that is
  correct here; **both `GetLineSegmentsSequential` and `GetLineSegmentsSweep` are wrong in the TPC**.

Impact: the LAr interaction target is mislabeled as G10 over meters → wrong materials/targets for
events whose ray crosses the affected TPC region. This is a **higher-priority correctness issue than
the traversal optimizations**. Likely related to the documented TPC-region descent bugs (exact-fill /
AABB-grazing; volTPCActive/volPixelPlane/volTPCPCB are the G10 pixel planes). Next: debug the
sequential daughter-descent in the TPC pixel-plane stack (why a thin G10 board's exit is computed as
~99 cm and LAr between boards is skipped).

## Status log

- 2026-06-15: Findings recorded, plan agreed. Implementing #1 and #2, stopping after each to run
  `[BVH]` tests and re-profile.
- 2026-06-15: **#1 done** (branchless `IntersectT`). NoRetry profile "19.49": `IntersectT` self
  19.7%->18.1%, `std::isfinite` self (total) 2.3%->0.65%, wall 924.9s->910.1s (~1.6%). Modest:
  `IntersectT` is inherently slab arithmetic and is still #1 because it runs for every BVH node on
  every segment step (call *frequency* is the real lever -> see ray-march TODO). `[BVH]` 272 asserts
  + full suite 857 asserts pass, bit-stable. Next: #2.
- 2026-06-15: **#2 done** (hoist invariant ray transforms out of `PhysicalVolume::GetLineSegments`;
  added no-division `Ray(origin, dir, inv_dir, pot)` ctor = the #3 piece). NoRetry profile "20.34":
  wall 910.1s->842.5s (~7.5%; cumulative vs baseline ~8.9%). `ApplyRayDirect` −36%, normalizing
  `Ray::Ray` divisions moved off the hot path. `[BVH]` + full suite pass, bit-stable (change is
  numerically identical by construction).
  **Residual finding:** `TransformRay` only dropped 12% because its dominant caller is now the **BVH
  leaf intersection** (`BVH::Traverse` -> `PhysicalVolume::Intersect` -> `TransformRay`), which
  re-transforms the world ray into each candidate daughter's frame per leaf test. That is exactly
  what **#4** targets and is now the clearest next transform win. `LogicalVolume::GetLineSegments` was
  left untouched — it is dead on the live path (root is a PV; all recursion stays in the PV version).
- 2026-06-15: **#4 done** (raw `PhysicalVolume*` through `BVH::TraverseIndex` + `GetLineSegments`;
  cached `m_nodes.data()`/`m_daughter_ptrs`; kept `shared_ptr` API for World/tests). Profile
  `profile_after_4.json.gz` (envelope/NoRetry, POT=1e16, 815.7s on-CPU). Targets eliminated:
  `vector<BVH::Node>::operator[]` 1.11%->gone, `__shared_ptr::operator bool` 0.61%->gone, per-Traverse
  shared_ptr copy gone; `BVH` traverse self 7.67%->7.13%. ~1.7-2% self removed, as predicted (small
  fixed overhead). Unchanged (not targeted): `IntersectT` ~18%, `Box::Intersect2Impl` ~7.7%,
  `TransformRay` ~6.5% (BVH-leaf re-transform -> #5 / per-daughter dir cache / ray-march). `[BVH]` 272
  + full 857 asserts pass.
  **Profiling caveats for this run:** (1) this samply build emitted the per-thread, *unsymbolicated*
  schema (earlier 4 were symbolicated Firefox exports) -> use `build/analyze_profile_sym.py` /
  `agg.py` which symbolicate native offsets via addr2line and aggregate by name; (2) ~3.5%
  `_dl_mcount_wrapper` tax from an instrumented (-pg) lib (libgeom itself is clean, 0 mcount) absent
  from earlier profiles, so absolute totals aren't comparable -> trust libgeom-internal composition;
  (3) seed still unfixed (run-to-run noise). Also: an idle-machine *suspend* once inflated a run's
  wall/CPU% -- samply samples on-CPU only, so suspend doesn't corrupt composition, just wall clock.
  **Post-geometry-opt note:** with the kernel faster, flux-file I/O (HepMC ASCII, `strtof64` ~1.2%
  on-CPU + off-CPU read/parse) is a growing share of *wall* time; see plan item #7 (binary flux
  format) in the original hotspots report.

Analysis scripts: `build/analyze_profile.py`, `build/analyze_callers.py`, `build/analyze_callees.py`
(usage: `python3 analyze_profile.py "<profile>.json.gz"`).
