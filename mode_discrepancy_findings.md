# Why the two interaction-determining mechanisms give different event counts

`DetectorSim::SamplingMode` offers two ways to decide whether a ray interacts:

| mode | layer 1 (NuGeometry) | layer 2 (generator) |
|---|---|---|
| `TotalXSecRetry` (`Mode: total`) | accept with `w·A·N_col·σ_tot(E) / M` | retry `GenerateEvent` until it emits |
| `EnvelopeNoRetry` (`Mode: envelope`) | accept with `w·A·N_col·σ_env / M` | one shot; a rejection produces no event |

They are supposed to give the same physics. On the DUNE ND hall they do not:
**`envelope` produces ~16% fewer events than `total` for the same POT.**

## Root cause: Achilles' unweighter is partial, so the modes' *counts* measure different things

`data/default/OptionDefaults.yml` in Achilles:

```yaml
Unweighting:
  Name: Percentile
  percentile: 99
```

`PercentileUnweighter::MaxValue()` is the **99th percentile** of the weight
distribution, not its maximum, and `AcceptEvent` is a *partial* unweighter:

```cpp
double prob = |w| / max_wgt;
if(prob < U(0,1)) return 0;          // P(accept) = min(1, prob)
return prob > 1 ? prob : 1;          // emitted weight = max(1, prob)
```

Write `r = w/max_w` and `⟨W⟩ = E[r] / E[min(1,r)]` — the mean weight of the
emitted events. Then for one trial:

* `P(emit)   = E[min(1,r)] = σ_tot / (max_w · ⟨W⟩)`
* `E[weight] = E[max(1,r)·1{accept}] = E[r] = σ_tot / max_w`

So the weights are unbiased but the **counts are not**, and the two modes put
the physical rate in different places:

```
EnvelopeNoRetry:  events/POT  = (1/⟨W⟩) · rate     sum(w)/POT = rate      ✓ weights
TotalXSecRetry:   events/POT  =           rate     sum(w)/POT = ⟨W⟩·rate  ✓ counts
```

`⟨W⟩` is exactly the factor by which the raw event counts disagree.

### Measured on the existing `build/achilles_geom_total.hepmc` (3644 events)

```
mean E.C.1 weight     = 1.1963
events with weight>1  = 31%          (enriched: overweight trials always accept)
max weight            = 9.77
overweight excess     = 1 − N/Σw = 16.4% of the cross section
```

⟹ predicted `total`/`envelope` event-count ratio **1.196**, i.e. `envelope`
short by 16.4%. That is far outside the ~1.7% statistical error of a
3644-event sample.

Note the 31%: only ~1% of *trials* are overweight, but they are accepted with
probability 1 while typical trials are accepted with probability `r ≪ 1`, so
they are heavily enriched in the emitted sample.

### ⟨W⟩ depends strongly on the beam energy range

That 1.196 was produced with `AutoEnergyRange` widening Achilles'
`Beam/EnergyRange` to the flux file's extremes — **[0.535, 80161] MeV**, set by a
handful of outlier rays. VEGAS, the 200-bin σ spline and `max_w` were all
calibrated over that range instead of the 0.1–10 GeV where essentially all of
the flux lives, which fattens the weight distribution and hence the overweight
tail. Re-running with the run card's [100, 10000] MeV gives

```
⟨W⟩ = 1.024 – 1.045      (2.4% – 4.3% overweight excess)
```

so the mode-to-mode discrepancy shrinks from ~16% to ~3–4% on the same geometry
and flux. The energy range now comes from the run card unless the flux file
records `NuGeom.Flux.EnergyRange_GeV` (or `Driver/FluxEnergyRange` is set).

### End-to-end, POT = 1e16, `nd_hall_with_lar_tms_nosand.gdml`, `rays.hepmc`

```
EnvelopeNoRetry   536 events / 2 624 035 rays, POT = 1.000e+16
  events/POT       = 5.360e-14        <- LOW by <weight>
  sum(weights)/POT = 5.602e-14        <- the physical rate
  <weight>         = 1.045
```

The prediction to check against is therefore that `total` mode's `events/POT`
lands on 5.602e-14, i.e. 4.5% above `envelope`'s raw count. **That run could not
be completed** at the time — see the driver crash below, since fixed.

### The A/B, now that the driver no longer crashes (2026-08-07)

Both modes ran to completion on `nd_hall_with_lar_tms_nosand.gdml` + `rays.hepmc`
with run cards differing **only** in `Driver/Mode`
(`InitRays: 500`, `BurnInRays: 1000`, `SafetyFactor: 1.5`, `AutoEnergyRange: true`).
Two independent pairs — full physics at matched POT, and a cheaper cascade-off
pair at matched event count. The RNG is entropy-seeded, so the pairs are
independent measurements.

| pair | mode | N | rays | POT | events/POT | Σw/POT | ⟨w⟩ |
|---|---|---|---|---|---|---|---|
| full physics | total    | 5303 | 2.25e6 | 1.000e16 | 5.303e-13 | 5.528e-13 | 1.0423 |
| full physics | envelope | 4887 | 2.88e7 | 1.000e16 | 4.887e-13 | 5.070e-13 | 1.0375 |
| cascade off  | total    | 3000 | 8.95e5 | 5.952e15 | 5.041e-13 | 5.258e-13 | 1.0430 |
| cascade off  | envelope | 3000 | 1.72e7 | 6.087e15 | 4.928e-13 | 5.126e-13 | 1.0402 |

**The physical-rate test** — `envelope`'s `Σw/POT` against `total`'s `events/POT`:

```
R = 0.9561 +/- 1.98%   full physics   (2.2 sigma)
R = 1.0170 +/- 2.58%   cascade off    (0.7 sigma)
------------------------------------------------
R = 0.9787 +/- 1.57%   combined       (1.4 sigma from 1)
```

So the two modes agree on the physical rate within statistics, with the right
estimator used for each mode. The raw **counts** differ as predicted
(total/envelope = 1.085 ± 2.0% and 1.023 ± 2.6%, against the predicted ⟨w⟩ ≈
1.04); the two pairs scatter to either side of 1.04, consistent with statistics.

Clipping is confirmed negligible at `SafetyFactor: 1.5`: 2 clipped rays per run,
0.000–0.003% of the rate — so the envelope is not biasing either mode.

**The one caveat**: the full-physics pair sits 2.2σ low on its own. That is not
significant, but it is the larger of the two samples, so a residual few-percent
effect is not excluded by this data. Re-running that pair (different seeds)
would settle it.

### Cost: `envelope` is 14–19× more expensive per event

Not a correctness issue, but the dominant practical difference, and it was not
previously measured:

```
rays per emitted event    total 424       envelope 5884    (full physics, 13.9x)
rays per emitted event    total 298       envelope 5746    (cascade off, 19.3x)
layer-1 envelope value    total 3.479e-10 envelope 3.063e-09  (8.8x)
```

`envelope` mode's layer 1 uses σ_env = max_w instead of σ_tot, an 8.8× larger
envelope here, so it must throw ~9× more rays to accumulate the same POT
(`m_pot += ray.POT()/envelope`) and then discards ~90% of its layer-1 accepts on
the single-shot trial (49 140 accepts → 4887 events, 9.9% efficiency, matching
Achilles' quoted ~15% unweighting efficiency). `total` mode retries instead, so
essentially every layer-1 accept becomes an event (5305 accepts → 5303 events).
For a fixed event budget, prefer `total`.

### Why the existing unit test did not catch it

`AlternatingGenerator` in `test_detsim.cc` is a *strict* unweighter — every
emitted event has weight exactly 1, so `⟨W⟩ = 1` and the modes agree exactly.
The new test **"A partial unweighter makes the two modes' raw event counts
differ by ⟨weight⟩"** adds a generator with an explicit overweight tail and
reproduces the divergence (measured ratio 1.555 against a predicted 1.554),
while confirming that `envelope`'s `Σw/POT` equals `total`'s `events/POT`.

## What to do about it

Nothing in NuGeometry can invent the missing events — the information is in the
weights. Three options, in order of preference:

1. **Use the right estimator per mode**: `Σw/POT` in `envelope`, `events/POT`
   in `total`. Both are now printed by `ReportRunSummary()` and by the driver,
   with a warning whenever `⟨W⟩` deviates from 1 by more than 2%.
2. **Make the envelope a true bound** — Achilles
   `Options/Unweighting/percentile: 100` (or `Name: None`). Then `⟨W⟩ = 1`,
   both modes' raw counts are the physical rate, and they agree.
3. If unweighted event *counts* must be right in `envelope` mode, Achilles'
   `VertexEnvelope` has to include the overweight excess, i.e. return
   `max_w · ⟨W⟩` rather than `max_w`.

## Secondary effect: envelope clipping (small here, but mode-asymmetric)

`ProduceEvent` accepts with probability 1 whenever `accept_w > max_prob`, while
still charging `POT/max_prob` — a downward bias on the rate. It is
mode-asymmetric because the two modes' `accept_w` distributions differ: `σ_env`
is energy-independent, so `envelope`'s spread comes only from the flux weight,
while `total` carries the extra `σ_tot(E)` tail and clips more.

Measured on `build/rays.hepmc` (961 060 rays) with the old `InitRays: 50000`,
`SafetyFactor: 2`:

```
max weight over the first 50 000 rays : 1.245e-05   (0.225 × the global max)
rays above 2 × that                   : 3 of 961 060, carrying 0.11% of Σw
```

An end-to-end cross-check agrees it is small: the same run with a clipping-free
envelope gave `events/POT = 5.530e-14` against `5.550e-14` with clipping — 0.4%
apart, well inside the ~4% statistical error of a 550-event sample. Clipping is
now counted and reported (`RunStats::clipped`, `clipped_excess`) so it can no
longer hide.

## RESOLVED: the combined Achilles driver crash

No longer reproduces. Four consecutive driver runs (both modes, two Achilles
configs, up to 2.9e7 traced rays and 5303 events) all exited 0 with a clean
`EventGen: finalizing run` and no `free(): invalid pointer`. The A/B above is
the run that this crash used to block. Original report kept below for context.

### Original report

`achilles_geom_driver` segfaults inside yaml-cpp's `YAML::Load`, while Achilles
parses its own config, immediately after NuGeometry's GDML parse:

```
#0  free () from /usr/lib/libc.so.6
#9  YAML::detail::node_data::push_back (...) at node_data.cpp:189
#12 YAML::NodeBuilder::Pop (...) at nodebuilder.cpp:112
#24 YAML::Load (input=...) at parse.cpp:25
```

i.e. heap corruption from somewhere earlier. It is **not** mode-dependent —
`envelope` and `total` both crash on the same card, and an identical card ran to
completion (14 min, 536 events) half an hour earlier, so it is intermittent.
It blocks a direct `envelope`-vs-`total` A/B at matched settings.

What has been ruled out:

* **NuGeometry.** The GDML parse plus the probe-ray `Init()` on the same
  geometry, isolated from Achilles, is clean under AddressSanitizer — no
  overflow, no use-after-free, no invalid free. (LeakSanitizer does report
  ~128 MB in 225 k allocations for the volume hierarchy; a leak, not a
  corruption, and pre-existing.) The 82-test suite is green.
* **An HepMC3 ABI skew.** `LD_LIBRARY_PATH` puts `~/.local/lib/libHepMC3.so.4`
  (1.1 MB, Oct 2025) ahead of the bundled 21 MB build the driver compiles
  against. Forcing the bundled one does not change the crash.
* **A yaml-cpp node-lifetime bug in the driver.** `FromRunCard` stored
  `card["Achilles"]` while the root document `card` died at return. That is
  genuinely unsafe (a sub-node shares the root's node arena) and is now fixed by
  keeping the root in `DriverConfig::runcard_root`, but it was not the cause.

Remaining suspect: the Achilles side, or the two builds' duplicate yaml-cpp /
HepMC3 stacks being linked into one process.

## Related: the auto-widened beam energy range

`AutoEnergyRange` used to widen Achilles' `Beam/EnergyRange` to the flux file's
extremes, which for `rays.hepmc` are **[0.535, 80161] MeV** — driven by a
handful of outlier rays. The VEGAS warm-up, the 200-bin σ spline and `max_w`
were all calibrated over that range instead of the 0.1–10 GeV where essentially
all of the flux lives. That inflates `max_w` (hence the unweighter's overweight
tail) and coarsens the spline. The range now comes from the run card unless the
flux file explicitly records `NuGeom.Flux.EnergyRange_GeV`.

## The E.C.4 cross section is NOT consistent between the modes (2026-08-08)

Final `GenCrossSection` from the four A/B files (pb):

| pair | total | envelope | total/env | significance |
|---|---|---|---|---|
| full physics | 0.111563 ± 0.000166 | 0.102921 ± 0.001428 | 1.0840 | 6.0σ |
| cascade off  | 0.111718 ± 0.000227 | 0.104447 ± 0.001875 | 1.0696 | 3.8σ |

This is by construction — `AchillesAdapter::GenerateEvent` uses a *different*
estimator per mode — but the two do not land on the same number. Decomposition:

**~2%: `total` mode's estimator is a σ²-weighted mean, not a flux-averaged one.**
Layer-1 accepts are drawn ∝ w·N_col·σ_tot(E), and each emitted event then
contributes `TotalXSec(E)` *again*, so `m_results.Mean()` converges to
E[σ²]/E[σ] rather than E[σ]. Backing out the per-event σ_tot(E_i) from the
running mean in `ab_tot.hepmc` (w_n = n·M_n − (n−1)·M_{n−1}):

```
arithmetic mean = 0.111563 pb   <- what GenCrossSection reports
harmonic  mean  = 0.109301 pb   <- = E[sigma], the flux-averaged value
arith/harm      = 1.0207
```

(The harmonic mean of a p ∝ w·N·σ sample is exactly E_{w·N}[σ], which is what
`envelope` mode estimates — so the harmonic mean is the apples-to-apples number.)

**~6%: `envelope` mode multiplies by a stale `VertexEnvelope`.**
`AchillesAdapter::CachedVertexEnvelope` caches M once per (nu, target) —
"the value is fixed for the run" — but Achilles' `PercentileUnweighter` is
adaptive: `Process::AddWeight` → `Unweighter::AddEvent` keeps feeding the
percentile *during* generation (`Process.cc:308`), and `AcceptEvent` reads
`m_percentile.Get()` fresh on every call. The estimator
`w_trial = Weight()·M_cached` then converges to `M_cached·σ/M_current`, i.e.
biased LOW whenever the percentile grows during the run — the direction observed
(0.1029 vs the 0.1093 harmonic mean, −5.8%).

**This does not affect the event rate**, which comes from NuGeometry's own POT
accounting and layer-1 acceptance. It affects only the E.C.4 number written into
the HepMC file, i.e. any downstream analysis that normalizes with it.

## Absolute normalization vs DUNE (2026-08-08)

The naive scaling `5.303e-13 events/POT x 1e21 = 5.3e8` is not comparable to a
DUNE detector rate, because **97% of those events are in the rock**:

```
ab_tot.hepmc, 5303 events:  Rock 97.04%, Scintillator 2.32%, EJ280WLS 0.34%,
                            Polyurethane 0.17%, G10 0.13%
```

and every event is on 12C (`TargetPDG 1000060120`) — the run card's only
nucleus — so the LAr detector contributes exactly nothing (argon has no
configured target, hence zero cross section, hence layer 1 never selects a LAr
segment). Detector materials only: 2.96% → 1.6e7 per 1e21 POT.

**Argon-only re-run** (`ar_tot.yml`, build/example_runcard.yml's Ar config:
QE-only, cascade off, EnergyRange [10, 100000], total mode, 3000 events):

```
3000 events / 1404126 rays, POT = 5.077e16
  events/POT = 5.909e-14   ->  5.9e7 events per 1e21 POT
  100% of vertices in LAr, all on Ar40 (1000180400)
```

The normalization is **internally consistent**. Independent check from the flux
file (961060 rays, Σw = 0.11578 cm^-2 over 1e6 POT → 1.158e-7 ν/cm²/POT, i.e.
1.16e14 ν/cm² per 1e21 POT — the right ballpark for the DUNE ND):

```
implied sigma from the measured rate = 0.188 pb / Ar40
Achilles' own GenCrossSection        = 0.225 pb / Ar40
```

agreeing to ~17%, which is within the crudeness of treating the LAr as a uniform
539.6 cm slab uniformly illuminated. So the engine's absolute rate is defensible;
the apparent 100–1000x excess over the quoted DUNE ~1e6 was the rock.

Caveats on the remaining argon-vs-DUNE gap (5.9e7 vs ~1e6):
* **No fiducial cut** — this is all LAr in the hall (bath, cryostat argon, ullage),
  not an active/fiducial volume.
* The vertices' transverse extent is **exactly the flux window** (x 698.7 cm,
  y 342.3 cm vs `FluxWindow.Size_cm 700.0, 342.424`), so any LAr outside the
  window is unilluminated — the rate is a *lower* bound in that respect.
* QE-only, so adding RES/SIS/DIS moves it further *up*, not down.

### Two defects surfaced by the argon run

1. **`kMaxTries` is being hit**: 55 of 3055 layer-1 accepts (1.8%) logged
   "Generator failed to emit after 100000 retries". In `TotalXSecRetry` those
   vertices charge POT but emit no event, biasing the rate ~1.8% LOW — and each
   costs 1e5 wasted `GenerateEvent` calls. Did not occur with carbon.
2. **⟨w⟩ = 1.41 for argon** (29% overweight tail) against 1.04 for carbon. The
   argon card uses `EnergyRange: [10, 100000]` where the carbon card uses
   [100, 10000] — exactly the max_w inflation documented above, now much worse.
