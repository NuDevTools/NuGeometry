#pragma once

#include "geom/FluxSource.hh"
#include "geom/Vector3D.hh"

#include <string>
#include <vector>

namespace NuGeom {

// Result of loading a NuHepMC (HepMC3 ASCII) flux file produced by
// dune_dk2nu/dk2nu_to_nuhepmc.py.
struct HepMCFlux {
    std::vector<FluxSample> samples; // rays, beam->detector transform applied
    double pot = 0.0;                // NuHepMC.Exposure.POT from the file
    double window_area = 1.0;        // NuGeom.FluxWindow.Area_cm2 (1.0 if absent)
    Vector3D offset;                 // translation (cm) actually applied
    std::string frame;               // NuGeom.CoordinateFrame ("beam")
    std::string length_unit;         // NuGeom.LengthUnit ("cm")
};

// Read a NuHepMC flux file into FluxSamples.
//
// Each event carries a decay vertex with the parent hadron incoming
// (status IncomingBeam) and the neutrino outgoing (status UndecayedPhysical).
// We take the outgoing neutrino: its production-vertex position is the ray
// origin (in the BEAM frame), its momentum the direction, and its energy E.
//
// The beam->detector translation recorded in the run-level attribute
// NuGeom.BeamToDetector.Translation_cm is added to every ray origin at load,
// since the rays are stored in their native beam frame (see the converter's
// coordinate-frame contract). Pass `offset_override` to ignore the recorded
// value and apply a different translation (cm) instead.
//
// Each FluxSample's Ray carries POT = total_pot / n_loaded so that summing
// over the thrown rays reproduces the file's beam exposure; flux_weight holds
// the per-ray dk2nu flux weight (nimpwt * nuray.wgt).
//
// `max_samples` > 0 stops after that many rays (a convenience for interactive
// previews such as BeamViz); with a partial load the per-ray POT is normalized
// to the loaded subset, so use the full file (max_samples = 0) for any POT-
// dependent result.
HepMCFlux LoadHepMCFlux(const std::string &path, const Vector3D *offset_override = nullptr,
                        std::size_t max_samples = 0);

} // namespace NuGeom
