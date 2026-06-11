#pragma once

#include <memory>

namespace HepMC3 {
class GenEvent;
class GenRunInfo;
} // namespace HepMC3

namespace NuGeom {

// Adapter base class for neutrino event generators.
//
// NuGeometry owns flux sampling, ray tracing, vertex/material selection.  The
// generator owns channel selection and full final-state kinematics.  All data
// exchange is via HepMC3/NuHepMC records so that generator-specific metadata
// rides along as attributes without changes to this interface.
class GeneratorInterface {
  public:
    virtual ~GeneratorInterface() = default;

    // Called once before event generation.  Generator populates run-level
    // NuHepMC metadata on `run_info` (process list, citations, etc.).
    virtual void InitRun(std::shared_ptr<HepMC3::GenRunInfo> run_info) = 0;

    // NuGeometry's layer 1 accepts a vertex with probability
    //     p_accept = W(E, nu, tgt) * n * L / max_geom
    // where W is one of two valid choices.  Implementations must supply
    // both so results can be cross-checked against each other.
    //
    // (1) TotalXSec — standard cross-section.
    //     Used with retry-to-emit in layer 2 so the generator produces
    //     exactly one event per accepted vertex.  Event weights follow the
    //     LHC partial-unweight convention
    //         w_emit = max(max_w, raw_w)
    //     with N_trials (including rejections) tracked generator-side and
    //     surfaced via E.C.4 / GR.6.  sigma estimator:
    //         Sum w_emit / N_trials -> sigma
    //     Unbiased for strict or partial unweighting.
    //
    // (2) EnvelopeXSec — unweighter's max_w (or equivalent upper bound).
    //     Used without retry; the generator's single per-vertex call may
    //     return no event, which is physically correct because max_w
    //     absorbs the unweighter inefficiency into the layer-1 envelope.
    //     Event weights are close to 1 (partial-unweighter tails at
    //     max(1, raw_w/max_w)).  Also exact.
    //
    // The two paths must yield the same physical observables; they differ
    // in where the unweighter inefficiency is bookkept (per-event weights
    // vs. layer-1 envelope).  DetectorSim selects between them via
    // SetSamplingMode.
    virtual double TotalXSec(double energy, int nu_pdg, int target_pdg) = 0;
    virtual double EnvelopeXSec(double energy, int nu_pdg, int target_pdg) = 0;

    // DetectorSim pre-fills `evt` with:
    //   * primary vertex at the chosen position
    //   * incoming neutrino (status 4) with the sampled 4-momentum
    //   * target nucleus (status 20) attached to that vertex
    //   * attribute "NuGeom.material" on the vertex
    // Generator adds outgoing particles, channel id (E.R.5), and per-event
    // attributes / weights (E.C.1).
    // Generate one phase-space sample and try to emit it onto `evt`.
    // Returns true if an event was emitted (evt populated with outgoing
    // particles and E.C.1 weight), false if the partial-unweighter rejected
    // this sample (evt untouched beyond what DetectorSim pre-filled).
    // DetectorSim decides whether to call this once or retry, based on the
    // SamplingMode (see DetectorSim.hh).
    virtual bool GenerateEvent(HepMC3::GenEvent &evt) = 0;

    virtual void FinalizeRun() {}
};

} // namespace NuGeom
