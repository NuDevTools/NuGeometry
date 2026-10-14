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
// Provenance of one flux ray: what produced the neutrino.  dk2nu records the
// parent hadron and the point at which it decayed, and the NuHepMC flux file
// carries them (parent as the status-4 particle, decay point as its production
// vertex).  Forwarding them onto the generated event is what lets an analysis
// tie an interaction back to the beamline decay that made its neutrino.
//
// `valid` is false for toy flux sources and for CSV files, which carry no
// parent information; nothing is written in that case.
struct FluxProvenance {
    bool valid = false;
    int parent_pdg = 0;     ///< dk2nu ptype: the decaying hadron
    double parent_px = 0.0; ///< parent momentum at decay, event momentum units
    double parent_py = 0.0;
    double parent_pz = 0.0;
    double parent_e = 0.0;
    double decay_x = 0.0; ///< decay point, event length units, DETECTOR frame
    double decay_y = 0.0;
    double decay_z = 0.0;
    double nimpwt = 0.0;    ///< dk2nu importance weight
    double nuray_wgt = 0.0; ///< dk2nu per-ray weight at this location
    int tgtexit_tptype = 0; ///< species exiting the target (0 = not recorded)
};

struct FluxSample {
    double energy = 0.0;
    int pdg = 0;
    Ray ray;
    double flux_weight = 1.0;
    double window_area = 1.0;
    std::function<void(HepMC3::GenEvent &)> decorate;
    FluxProvenance provenance{};
};

using FluxCallback = std::function<FluxSample()>;

} // namespace NuGeom
