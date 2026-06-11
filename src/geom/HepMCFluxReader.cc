#include "geom/HepMCFluxReader.hh"
#include "geom/Ray.hh"

// HepMC3 and NuHepMC headers trip our strict warnings; silence them here
// (same pattern as DetectorSim.cc).
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "HepMC3/Attribute.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/ReaderAscii.h"

#include "NuHepMC/Constants.hxx"
#pragma GCC diagnostic pop

#include "geom/Logging.hh"

#include <array>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

namespace {

constexpr int kStatusFinal = NuHepMC::ParticleStatus::UndecayedPhysical; // 1

bool IsNeutrino(int pdg) {
    const int a = std::abs(pdg);
    return a == 12 || a == 14 || a == 16;
}

// Read a run-level string attribute, empty if absent.
std::string RunAttr(const std::shared_ptr<HepMC3::GenRunInfo> &ri, const std::string &name) {
    if(!ri) return {};
    auto a = ri->attribute<HepMC3::StringAttribute>(name);
    return a ? a->value() : std::string{};
}

// Parse "x,y,z" (cm) into a Vector3D; throws on malformed input.
NuGeom::Vector3D ParseTriplet(const std::string &s, const std::string &ctx) {
    std::array<double, 3> v{0.0, 0.0, 0.0};
    std::size_t start = 0;
    for(int i = 0; i < 3; ++i) {
        std::size_t comma = s.find(',', start);
        const std::string tok =
            s.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        try {
            v[static_cast<std::size_t>(i)] = std::stod(tok);
        } catch(const std::exception &) {
            throw std::runtime_error("HepMCFluxReader: malformed " + ctx + ": '" + s + "'");
        }
        if(comma == std::string::npos) {
            if(i != 2)
                throw std::runtime_error("HepMCFluxReader: expected 3 values in " + ctx + ": '" +
                                         s + "'");
            break;
        }
        start = comma + 1;
    }
    return {v[0], v[1], v[2]};
}

} // namespace

namespace NuGeom {

HepMCFlux LoadHepMCFlux(const std::string &path, const Vector3D *offset_override,
                        std::size_t max_samples) {
    HepMC3::ReaderAscii reader(path);
    if(reader.failed()) throw std::runtime_error("HepMCFluxReader: cannot open flux file: " + path);

    HepMCFlux out;
    bool metadata_read = false;

    // Raw per-ray data; rays get their POT after we know the total count.
    struct Raw {
        int pdg;
        double energy;
        double weight;
        Vector3D origin; // beam frame + applied offset (cm)
        Vector3D direction;
    };
    std::vector<Raw> raws;

    std::size_t n_no_neutrino = 0;
    std::size_t n_no_vertex = 0;

    HepMC3::GenEvent evt;
    while(true) {
        reader.read_event(evt);
        if(reader.failed()) break; // clean EOF or read error -> stop

        if(!metadata_read) {
            const auto &ri = evt.run_info();
            out.frame = RunAttr(ri, "NuGeom.CoordinateFrame");
            out.length_unit = RunAttr(ri, "NuGeom.LengthUnit");

            const std::string pot_s = RunAttr(ri, "NuHepMC.Exposure.POT");
            if(!pot_s.empty()) {
                try {
                    out.pot = std::stod(pot_s);
                } catch(const std::exception &) {
                    NuGeom::Log().warn("HepMCFluxReader: could not parse NuHepMC.Exposure.POT='{}'",
                                       pot_s);
                }
            }

            if(offset_override) {
                out.offset = *offset_override;
                NuGeom::Log().info("HepMCFluxReader: overriding recorded transform with "
                                   "({}, {}, {}) cm",
                                   out.offset.X(), out.offset.Y(), out.offset.Z());
            } else {
                const std::string t_s = RunAttr(ri, "NuGeom.BeamToDetector.Translation_cm");
                if(t_s.empty())
                    NuGeom::Log().warn(
                        "HepMCFluxReader: no NuGeom.BeamToDetector.Translation_cm in "
                        "'{}'; applying identity transform",
                        path);
                else
                    out.offset = ParseTriplet(t_s, "NuGeom.BeamToDetector.Translation_cm");
            }

            if(!out.length_unit.empty() && out.length_unit != "cm")
                NuGeom::Log().warn("HepMCFluxReader: file length unit is '{}', not 'cm'; the "
                                   "transform is in cm and positions are read in the event's "
                                   "native unit -- check consistency",
                                   out.length_unit);
            metadata_read = true;
        }

        // Find the outgoing (status 1) neutrino.
        HepMC3::ConstGenParticlePtr nu;
        for(const auto &p : evt.particles()) {
            if(p->status() == kStatusFinal && IsNeutrino(p->pid())) {
                nu = p;
                break;
            }
        }
        if(!nu) {
            ++n_no_neutrino;
            continue;
        }
        const auto vtx = nu->production_vertex();
        if(!vtx) {
            ++n_no_vertex;
            continue;
        }

        const auto &pos = vtx->position();
        const auto &mom = nu->momentum();

        double weight = 1.0;
        if(!evt.weights().empty())
            weight = evt.weights().front();
        else {
            auto wa = evt.attribute<HepMC3::DoubleAttribute>("NuGeom.flux_weight");
            if(wa) weight = wa->value();
        }

        Raw r;
        r.pdg = nu->pid();
        r.energy = mom.e();
        r.weight = weight;
        r.origin = Vector3D{pos.x(), pos.y(), pos.z()} + out.offset;
        r.direction = Vector3D{mom.px(), mom.py(), mom.pz()};
        raws.push_back(std::move(r));

        if(max_samples > 0 && raws.size() >= max_samples) break;
    }

    if(raws.empty())
        throw std::runtime_error("HepMCFluxReader: no neutrino rays parsed from " + path);

    // POT carried per ray so that summing over thrown rays reproduces the beam
    // exposure recorded in the file.  If the file lacks POT, fall back to 1.0.
    const double pot_per_ray = out.pot > 0.0 ? out.pot / static_cast<double>(raws.size()) : 1.0;

    out.samples.reserve(raws.size());
    for(const auto &r : raws) {
        FluxSample fs;
        fs.energy = r.energy;
        fs.pdg = r.pdg;
        fs.ray = Ray(r.origin, r.direction, pot_per_ray);
        fs.flux_weight = r.weight;
        out.samples.push_back(std::move(fs));
    }

    NuGeom::Log().info("HepMCFluxReader: loaded {} rays from '{}' (frame='{}', POT={:.6e}, "
                       "offset=({}, {}, {}) cm)",
                       out.samples.size(), path, out.frame.empty() ? "?" : out.frame, out.pot,
                       out.offset.X(), out.offset.Y(), out.offset.Z());
    if(n_no_neutrino || n_no_vertex)
        NuGeom::Log().warn("HepMCFluxReader: skipped {} events without a status-1 neutrino, {} "
                           "without a production vertex",
                           n_no_neutrino, n_no_vertex);

    return out;
}

} // namespace NuGeom
