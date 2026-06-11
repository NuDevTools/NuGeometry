#include "geom/FluxStreamer.hh"
#include "geom/Logging.hh"
#include "geom/Ray.hh"

// HepMC3 headers trip our strict warnings; silence them here (same pattern as
// HepMCFluxReader.cc).
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

#include <array>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

namespace {

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
            throw std::runtime_error("FluxStreamer: malformed " + ctx + ": '" + s + "'");
        }
        if(comma == std::string::npos) {
            if(i != 2)
                throw std::runtime_error("FluxStreamer: expected 3 values in " + ctx + ": '" + s +
                                         "'");
            break;
        }
        start = comma + 1;
    }
    return {v[0], v[1], v[2]};
}

// Parse one CSV row (pid,wgt,E,px,py,pz,x,y,z); false on blank/malformed rows.
bool ParseCSVRow(const std::string &line, NuGeom::FluxSample &fs) {
    std::istringstream ss(line);
    std::string field;
    double cols[9];
    std::size_t n = 0;
    for(; n < 9 && std::getline(ss, field, ','); ++n) {
        try {
            cols[n] = std::stod(field);
        } catch(const std::exception &) { return false; }
    }
    if(n < 9) return false;

    fs.energy = cols[2];
    fs.pdg = static_cast<int>(cols[0]);
    fs.ray = NuGeom::Ray(NuGeom::Vector3D{cols[6], cols[7], cols[8]},
                         NuGeom::Vector3D{cols[3], cols[4], cols[5]}, 1.0);
    fs.flux_weight = cols[1];
    fs.decorate = nullptr;
    return true;
}

// Pull the ray out of one NuHepMC flux event: the outgoing (status 1)
// neutrino's production vertex is the origin (beam frame), its momentum the
// direction.  Returns false for events without a usable neutrino.
// Mirrors the extraction in HepMCFluxReader.cc (LoadHepMCFlux).
bool ExtractHepMCRay(const HepMC3::GenEvent &evt, const NuGeom::Vector3D &offset,
                     double pot_per_ray, NuGeom::FluxSample &fs) {
    constexpr int kStatusFinal = NuHepMC::ParticleStatus::UndecayedPhysical; // 1

    HepMC3::ConstGenParticlePtr nu;
    for(const auto &p : evt.particles()) {
        if(p->status() == kStatusFinal && IsNeutrino(p->pid())) {
            nu = p;
            break;
        }
    }
    if(!nu) return false;
    const auto vtx = nu->production_vertex();
    if(!vtx) return false;

    const auto &pos = vtx->position();
    const auto &mom = nu->momentum();

    double weight = 1.0;
    if(!evt.weights().empty())
        weight = evt.weights().front();
    else {
        auto wa = evt.attribute<HepMC3::DoubleAttribute>("NuGeom.flux_weight");
        if(wa) weight = wa->value();
    }

    fs.energy = mom.e();
    fs.pdg = nu->pid();
    fs.ray = NuGeom::Ray(NuGeom::Vector3D{pos.x(), pos.y(), pos.z()} + offset,
                         NuGeom::Vector3D{mom.px(), mom.py(), mom.pz()}, pot_per_ray);
    fs.flux_weight = weight;
    fs.decorate = nullptr;
    return true;
}

} // namespace

namespace NuGeom {

FluxSample FluxStreamer::Next() {
    FluxSample fs;
    if(TryNext(fs)) return fs;
    if(!m_loop)
        throw std::runtime_error("FluxStreamer: end of flux file '" + m_path + "' reached after " +
                                 std::to_string(m_count) +
                                 " rays; enable ray looping to rewind instead");
    Rewind();
    ++m_loops;
    NuGeom::Log().debug("FluxStreamer: rewound '{}' (pass {})", m_path, m_loops + 1);
    if(!TryNext(fs))
        throw std::runtime_error("FluxStreamer: no rays found in '" + m_path + "' after rewind");
    return fs;
}

// --- CSV ---------------------------------------------------------------

CSVFluxStreamer::CSVFluxStreamer(const std::string &path, bool loop) : FluxStreamer(path, loop) {
    std::ifstream pass(path);
    if(!pass) throw std::runtime_error("FluxStreamer: cannot open rays file: " + path);

    // Summary pass: count rays and record the energy range.
    std::string line;
    bool header_skipped = false;
    FluxSample fs;
    while(std::getline(pass, line)) {
        if(!line.empty() && line.back() == '\r') line.pop_back();
        if(line.empty()) continue;
        if(!header_skipped) {
            header_skipped = true;
            continue;
        }
        if(!ParseCSVRow(line, fs)) continue;
        ++m_count;
        m_emin = std::min(m_emin, fs.energy);
        m_emax = std::max(m_emax, fs.energy);
    }
    if(m_count == 0) throw std::runtime_error("FluxStreamer: no rays parsed from " + path);

    m_in.open(path);
    Rewind();
}

void CSVFluxStreamer::Rewind() {
    m_in.clear();
    m_in.seekg(0);
    // Skip the header row (first non-empty line).
    std::string line;
    while(std::getline(m_in, line)) {
        if(!line.empty() && line.back() == '\r') line.pop_back();
        if(!line.empty()) break;
    }
}

bool CSVFluxStreamer::TryNext(FluxSample &fs) {
    std::string line;
    while(std::getline(m_in, line)) {
        if(!line.empty() && line.back() == '\r') line.pop_back();
        if(line.empty()) continue;
        if(ParseCSVRow(line, fs)) return true;
    }
    return false;
}

// --- NuHepMC -----------------------------------------------------------

HepMCFluxStreamer::HepMCFluxStreamer(const std::string &path, bool loop,
                                     const Vector3D *offset_override)
    : FluxStreamer(path, loop) {
    Summarize(offset_override);
    Rewind();
}

HepMCFluxStreamer::~HepMCFluxStreamer() = default;

void HepMCFluxStreamer::Summarize(const Vector3D *offset_override) {
    HepMC3::ReaderAscii reader(m_path);
    if(reader.failed()) throw std::runtime_error("FluxStreamer: cannot open flux file: " + m_path);

    bool metadata_read = false;
    std::size_t n_skipped = 0;
    HepMC3::GenEvent evt;
    FluxSample fs;
    while(true) {
        reader.read_event(evt);
        if(reader.failed()) break; // clean EOF or read error -> stop

        if(!metadata_read) {
            const auto &ri = evt.run_info();
            const std::string pot_s = RunAttr(ri, "NuHepMC.Exposure.POT");
            if(!pot_s.empty()) {
                try {
                    m_pot = std::stod(pot_s);
                } catch(const std::exception &) {
                    NuGeom::Log().warn("FluxStreamer: could not parse NuHepMC.Exposure.POT='{}'",
                                       pot_s);
                }
            }

            if(offset_override) {
                m_offset = *offset_override;
                NuGeom::Log().info(
                    "FluxStreamer: overriding recorded transform with ({}, {}, {}) cm",
                    m_offset.X(), m_offset.Y(), m_offset.Z());
            } else {
                const std::string t_s = RunAttr(ri, "NuGeom.BeamToDetector.Translation_cm");
                if(t_s.empty())
                    NuGeom::Log().warn("FluxStreamer: no NuGeom.BeamToDetector.Translation_cm in "
                                       "'{}'; applying identity transform",
                                       m_path);
                else
                    m_offset = ParseTriplet(t_s, "NuGeom.BeamToDetector.Translation_cm");
            }

            const std::string unit = RunAttr(ri, "NuGeom.LengthUnit");
            if(!unit.empty() && unit != "cm")
                NuGeom::Log().warn("FluxStreamer: file length unit is '{}', not 'cm'; the "
                                   "transform is in cm and positions are read in the event's "
                                   "native unit -- check consistency",
                                   unit);
            metadata_read = true;
        }

        if(!ExtractHepMCRay(evt, m_offset, 1.0, fs)) {
            ++n_skipped;
            continue;
        }
        ++m_count;
        m_emin = std::min(m_emin, fs.energy);
        m_emax = std::max(m_emax, fs.energy);
    }

    if(m_count == 0)
        throw std::runtime_error("FluxStreamer: no neutrino rays parsed from " + m_path);

    // POT carried per ray so a full pass over the file reproduces the recorded
    // beam exposure (same normalization as LoadHepMCFlux on a full load).
    m_pot_per_ray = m_pot > 0.0 ? m_pot / static_cast<double>(m_count) : 1.0;

    NuGeom::Log().info("FluxStreamer: '{}' has {} rays (E in [{}, {}] GeV, POT={:.6e}, "
                       "offset=({}, {}, {}) cm)",
                       m_path, m_count, m_emin, m_emax, m_pot, m_offset.X(), m_offset.Y(),
                       m_offset.Z());
    if(n_skipped)
        NuGeom::Log().warn("FluxStreamer: skipped {} events without a usable neutrino ray",
                           n_skipped);
}

void HepMCFluxStreamer::Rewind() {
    // HepMC3 readers are forward-only; reopen the file for each pass.
    m_reader = std::make_unique<HepMC3::ReaderAscii>(m_path);
    if(m_reader->failed())
        throw std::runtime_error("FluxStreamer: cannot reopen flux file: " + m_path);
}

bool HepMCFluxStreamer::TryNext(FluxSample &fs) {
    HepMC3::GenEvent evt;
    while(true) {
        m_reader->read_event(evt);
        if(m_reader->failed()) return false;
        if(ExtractHepMCRay(evt, m_offset, m_pot_per_ray, fs)) return true;
    }
}

// --- Factory -----------------------------------------------------------

namespace {
bool EndsWith(const std::string &s, const std::string &suffix) {
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}
} // namespace

std::unique_ptr<FluxStreamer> OpenFluxStreamer(const std::string &path, bool loop,
                                               const Vector3D *offset_override) {
    if(EndsWith(path, ".hepmc") || EndsWith(path, ".hepmc3"))
        return std::make_unique<HepMCFluxStreamer>(path, loop, offset_override);
    return std::make_unique<CSVFluxStreamer>(path, loop);
}

} // namespace NuGeom
