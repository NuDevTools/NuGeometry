#pragma once

#include "geom/Ray.hh"

#include <functional>

namespace HepMC3 {
class GenEvent;
}

namespace NuGeom {

// One sample from the incoming neutrino flux.  flux_weight is the per-sample
// importance weight (e.g. dk2nu `wgt`), distinct from Ray::POT() which carries
// the file-level POT normalisation.  The product
//     flux_weight * ray.POT() / max_prob
// is the POT contribution of a single thrown sample under accept/reject.
//
// `decorate` (optional) writes flux-specific NuHepMC attributes onto the event
// that eventually results from this sample (parent hadron info, flux window,
// etc.).  Left null for toy flux sources.
struct FluxSample {
    double energy = 0.0;
    int pdg = 0;
    Ray ray;
    double flux_weight = 1.0;
    std::function<void(HepMC3::GenEvent &)> decorate;
};

using FluxCallback = std::function<FluxSample()>;

} // namespace NuGeom
