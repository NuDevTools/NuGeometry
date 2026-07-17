#pragma once

#include "geom/Ray.hh"

#include <functional>

namespace HepMC3 {
class GenEvent;
}

namespace NuGeom {

// One sample from the incoming neutrino flux.  flux_weight is the per-sample
// importance weight (e.g. dk2nu `wgt`, per cm^2), distinct from Ray::POT() which
// carries the file-level POT normalisation.
//
// window_area [cm^2] is the transverse area of the flux window over which the
// ray origins were spread when the file was produced (NuGeom.FluxWindow.Area_cm2
// in the NuHepMC run info; 1.0 if absent).  Because the per-cm^2 flux_weight is
// only the rate *per unit area*, the detector event rate needs the area
// integral over the window; with the origins already Monte-Carlo spread across
// it, that integral is simply
//     accept_w = flux_weight * window_area * interaction_probability
// summed over rays.  Without the window_area factor (all rays converging to a
// point) the rate comes out low by ~window_area (e.g. ~2.4e5 cm^2 for ND-LAr).
//
// `decorate` (optional) writes flux-specific NuHepMC attributes onto the event
// that eventually results from this sample (parent hadron info, flux window,
// etc.).  Left null for toy flux sources.
struct FluxSample {
    double energy = 0.0;
    int pdg = 0;
    Ray ray;
    double flux_weight = 1.0;
    double window_area = 1.0;
    std::function<void(HepMC3::GenEvent &)> decorate;
};

using FluxCallback = std::function<FluxSample()>;

} // namespace NuGeom
