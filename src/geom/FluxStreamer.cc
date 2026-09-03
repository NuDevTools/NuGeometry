#include "geom/FluxStreamer.hh"
#include "geom/Logging.hh"
#include "geom/Random.hh"
#include "geom/Ray.hh"
#include <algorithm>

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
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

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
                     double pot_per_ray, double window_area, NuGeom::FluxSample &fs) {
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
    fs.window_area = window_area;
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
    ReadMetadata(offset_override);
    Rewind();
}

HepMCFluxStreamer::~HepMCFluxStreamer() = default;

namespace {

// Count HepMC3 ASCII event records without parsing them.  Every event begins a
// line with "E "; no other record type does, and the file's own header/footer
// lines start with "HepMC::".  Reading raw blocks and scanning for the marker
// costs a single sequential pass at I/O speed, against the ~20x more expensive
// alternative of materializing a GenEvent (with its vertices and particles)
// per ray only to throw it away.
std::size_t CountHepMCEvents(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    if(!in) throw std::runtime_error("FluxStreamer: cannot open flux file: " + path);

    constexpr std::size_t kBlock = 1 << 20;
    std::vector<char> buf(kBlock);
    std::size_t count = 0;
    bool at_line_start = true;
    // "E " can straddle a block boundary: remember a trailing 'E' seen at the
    // start of a line so the following block's first byte can complete it.
    bool pending_e = false;

    while(in) {
        in.read(buf.data(), static_cast<std::streamsize>(kBlock));
        const std::size_t n = static_cast<std::size_t>(in.gcount());
        for(std::size_t i = 0; i < n; ++i) {
            const char c = buf[i];
            if(pending_e) {
                if(c == ' ') ++count;
                pending_e = false;
                at_line_start = (c == '\n');
                continue;
            }
            if(at_line_start && c == 'E') {
                pending_e = true;
                continue;
            }
            at_line_start = (c == '\n');
        }
    }
    return count;
}

} // namespace

void HepMCFluxStreamer::ReadMetadata(const Vector3D *offset_override) {
    // One event is enough: HepMC3 attaches the run info to every event it
    // reads, so the run-level attributes are available after the first record.
    HepMC3::ReaderAscii reader(m_path);
    if(reader.failed()) throw std::runtime_error("FluxStreamer: cannot open flux file: " + m_path);

    HepMC3::GenEvent evt;
    reader.read_event(evt);
    if(reader.failed())
        throw std::runtime_error("FluxStreamer: no events found in flux file: " + m_path);

    const auto &ri = evt.run_info();
    const std::string pot_s = RunAttr(ri, "NuHepMC.Exposure.POT");
    if(!pot_s.empty()) {
        try {
            m_pot = std::stod(pot_s);
        } catch(const std::exception &) {
            NuGeom::Log().warn("FluxStreamer: could not parse NuHepMC.Exposure.POT='{}'", pot_s);
        }
    }

    if(offset_override) {
        m_offset = *offset_override;
        NuGeom::Log().info("FluxStreamer: overriding recorded transform with ({}, {}, {}) cm",
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

    const std::string area_s = RunAttr(ri, "NuGeom.FluxWindow.Area_cm2");
    if(!area_s.empty()) {
        try {
            m_window_area = std::stod(area_s);
        } catch(const std::exception &) {
            NuGeom::Log().warn("FluxStreamer: could not parse NuGeom.FluxWindow.Area_cm2='{}'",
                               area_s);
        }
    } else {
        NuGeom::Log().warn("FluxStreamer: no NuGeom.FluxWindow.Area_cm2 in '{}'; using 1.0 cm^2 "
                           "(event rate will be the per-cm^2 flux rate, NOT integrated over the "
                           "detector face -- regenerate the flux with a flux window to fix the "
                           "normalization)",
                           m_path);
    }

    const std::string unit = RunAttr(ri, "NuGeom.LengthUnit");
    if(!unit.empty() && unit != "cm")
        NuGeom::Log().warn("FluxStreamer: file length unit is '{}', not 'cm'; the "
                           "transform is in cm and positions are read in the event's "
                           "native unit -- check consistency",
                           unit);

    // Ray count: only needed to turn the file POT into a per-ray exposure.
    const std::string n_s = RunAttr(ri, "NuGeom.Flux.NRays");
    if(!n_s.empty()) {
        try {
            m_count = static_cast<std::size_t>(std::stoull(n_s));
        } catch(const std::exception &) {
            NuGeom::Log().warn("FluxStreamer: could not parse NuGeom.Flux.NRays='{}'", n_s);
        }
    }
    if(m_count == 0) {
        NuGeom::Log().info("FluxStreamer: '{}' has no NuGeom.Flux.NRays attribute; counting "
                           "event records (regenerate the flux to record the count and skip "
                           "this pass)",
                           m_path);
        m_count = CountHepMCEvents(m_path);
    }
    if(m_count == 0) throw std::runtime_error("FluxStreamer: no rays found in " + m_path);

    // Energy range: reported only if the file records it.  Deriving it means
    // reading every ray, which is what this streamer exists to avoid.
    const std::string e_s = RunAttr(ri, "NuGeom.Flux.EnergyRange_GeV");
    if(!e_s.empty()) {
        const auto comma = e_s.find(',');
        try {
            m_emin = std::stod(e_s.substr(0, comma));
            m_emax = std::stod(e_s.substr(comma + 1));
        } catch(const std::exception &) {
            NuGeom::Log().warn("FluxStreamer: could not parse NuGeom.Flux.EnergyRange_GeV='{}'",
                               e_s);
        }
    }

    // Flux-weight summary, if the converter recorded it.  Knowing mean(w)
    // without reading the file lets the driver report the importance-sampling
    // payoff (max/mean) up front, and cross-checks the cache once it is built.
    const std::string mw_s = RunAttr(ri, "NuGeom.Flux.MeanWeight");
    if(!mw_s.empty()) {
        try {
            m_mean_weight = std::stod(mw_s);
        } catch(const std::exception &) {
            NuGeom::Log().warn("FluxStreamer: could not parse NuGeom.Flux.MeanWeight='{}'", mw_s);
        }
    }
    const std::string tw_s = RunAttr(ri, "NuGeom.Flux.TotalWeight");
    if(!tw_s.empty()) {
        try {
            m_total_weight = std::stod(tw_s);
        } catch(const std::exception &) {
            NuGeom::Log().warn("FluxStreamer: could not parse NuGeom.Flux.TotalWeight='{}'", tw_s);
        }
    }
    // Either one implies the other given the ray count.
    if(m_mean_weight <= 0.0 && m_total_weight > 0.0 && m_count > 0)
        m_mean_weight = m_total_weight / static_cast<double>(m_count);
    if(m_total_weight <= 0.0 && m_mean_weight > 0.0 && m_count > 0)
        m_total_weight = m_mean_weight * static_cast<double>(m_count);

    // POT carried per ray so a full pass over the file reproduces the recorded
    // beam exposure (same normalization as LoadHepMCFlux on a full load).
    m_pot_per_ray = m_pot > 0.0 ? m_pot / static_cast<double>(m_count) : 1.0;

    if(HasEnergyRange())
        NuGeom::Log().info("FluxStreamer: '{}' has {} rays (E in [{}, {}] GeV, POT={:.6e}, "
                           "window_area={:.6e} cm^2, offset=({}, {}, {}) cm)",
                           m_path, m_count, m_emin, m_emax, m_pot, m_window_area, m_offset.X(),
                           m_offset.Y(), m_offset.Z());
    else
        NuGeom::Log().info("FluxStreamer: '{}' has {} rays (energy range not recorded, POT={:.6e}, "
                           "window_area={:.6e} cm^2, offset=({}, {}, {}) cm)",
                           m_path, m_count, m_pot, m_window_area, m_offset.X(), m_offset.Y(),
                           m_offset.Z());
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
        if(ExtractHepMCRay(evt, m_offset, m_pot_per_ray, m_window_area, fs)) return true;
    }
}

// --- Cached / importance-sampling decorator ----------------------------

CachedFluxStreamer::CachedFluxStreamer(std::unique_ptr<FluxStreamer> inner, Options opts)
    : FluxStreamer("", inner->Looping()), m_inner{std::move(inner)}, m_caching{opts.cache},
      m_importance{opts.importance}, m_max_bytes{opts.max_bytes} {
    m_count = m_inner->Count();
    m_emin = m_inner->EMin();
    m_emax = m_inner->EMax();
    if(m_importance && !m_caching) {
        NuGeom::Log().warn("FluxStreamer: importance sampling needs random access to the rays; "
                           "enabling the in-memory cache for it.");
        m_caching = true;
    }
    // The envelope only ever grows, so a run that starts on real weights locks
    // in a max(w)-sized envelope; arming importance sampling later then buys
    // nothing.  Preload so the very first thrown ray already carries the mean.
    bool preload = opts.preload || m_importance;
    // Report the payoff up front when the file recorded its weight summary.
    if(m_inner->HasWeightSummary()) {
        m_mean_weight = m_inner->MeanFluxWeight();
        m_total_weight = m_inner->TotalFluxWeight();
    }
    if(m_caching && m_count > 0 && m_max_bytes > 0) {
        const std::size_t need = m_count * sizeof(Cached);
        if(need > m_max_bytes)
            NuGeom::Log().warn("FluxStreamer: caching {} rays would need {:.1f} MB but the cap is "
                               "{:.1f} MB; the cache will fill until the cap and then fall back to "
                               "streaming.",
                               m_count, static_cast<double>(need) / 1048576.0,
                               static_cast<double>(m_max_bytes) / 1048576.0);
        else
            m_rays.reserve(m_count);
    }
    if(preload && m_caching) Preload();
}

void CachedFluxStreamer::Preload() {
    NuGeom::Log().info("FluxStreamer: preloading rays into memory...");
    FluxSample fs;
    while(TryNext(fs)) {
        if(!m_caching) { // hit the cap: nothing more to gain from preloading
            NuGeom::Log().warn("FluxStreamer: preload aborted at the cache cap; importance "
                               "sampling is unavailable and the run will stream.");
            m_importance = false;
            return;
        }
    }
    // TryNext() returning false on the first pass already flipped m_filled and
    // built the sampler; the cursor just needs resetting for the run proper.
    m_pos = 0;
}

void CachedFluxStreamer::Emit(const Cached &c, FluxSample &fs, double weight_override) const {
    fs.energy = c.energy;
    fs.pdg = c.pdg;
    fs.ray = NuGeom::Ray(NuGeom::Vector3D{c.ox, c.oy, c.oz}, NuGeom::Vector3D{c.dx, c.dy, c.dz},
                         m_pot_per_ray);
    fs.flux_weight = weight_override > 0.0 ? weight_override : c.weight;
    fs.window_area = m_window_area;
    fs.decorate = nullptr;
}

void CachedFluxStreamer::BuildSampler() {
    m_cumulative.clear();
    m_cumulative.reserve(m_rays.size());
    double running = 0.0;
    for(const auto &r : m_rays) {
        running += r.weight > 0.0 ? r.weight : 0.0;
        m_cumulative.push_back(running);
    }
    m_total_w = running;
    m_mean_w = m_rays.empty() ? 0.0 : m_total_w / static_cast<double>(m_rays.size());
    if(m_mean_weight > 0.0 && m_mean_w > 0.0) {
        const double ratio = m_mean_w / m_mean_weight;
        if(ratio < 0.99 || ratio > 1.01)
            NuGeom::Log().warn("FluxStreamer: mean flux weight from the rays ({:.6e}) disagrees "
                               "with the file's recorded NuGeom.Flux.MeanWeight ({:.6e}) by "
                               "{:.2f}%. The metadata may be stale.",
                               m_mean_w, m_mean_weight, 100.0 * (ratio - 1.0));
    }
    m_mean_weight = m_mean_w;
    m_total_weight = m_total_w;
    if(!(m_total_w > 0.0)) {
        NuGeom::Log().warn("FluxStreamer: total flux weight is {}; importance sampling disabled.",
                           m_total_w);
        m_importance = false;
        return;
    }
    NuGeom::Log().info("FluxStreamer: importance sampling armed over {} rays -- mean weight {:.6e} "
                       "vs max {:.6e} (factor {:.0f} fewer rays for the same POT).",
                       m_rays.size(), m_mean_w,
                       m_rays.empty() ? 0.0
                                      : std::max_element(m_rays.begin(), m_rays.end(),
                                                         [](const Cached &a, const Cached &b) {
                                                             return a.weight < b.weight;
                                                         })
                                            ->weight,
                       m_mean_w > 0.0 ? (std::max_element(m_rays.begin(), m_rays.end(),
                                                          [](const Cached &a, const Cached &b) {
                                                              return a.weight < b.weight;
                                                          })
                                             ->weight /
                                         m_mean_w)
                                      : 0.0);
}

bool CachedFluxStreamer::TryNext(FluxSample &fs) {
    // Serving from memory.
    if(m_filled) {
        if(m_rays.empty()) return false;
        if(m_importance) {
            // Draw i with probability w_i / sum(w); hand it the mean weight.
            const double u = NuGeom::Random::Instance().Uniform(0.0, m_total_w);
            auto it = std::upper_bound(m_cumulative.begin(), m_cumulative.end(), u);
            std::size_t idx = static_cast<std::size_t>(it - m_cumulative.begin());
            if(idx >= m_rays.size()) idx = m_rays.size() - 1;
            Emit(m_rays[idx], fs, m_mean_w);
            // No end-of-file in sampling mode; count pseudo-passes so the run
            // report still means "exposures worth of rays drawn".
            if(++m_draws % m_rays.size() == 0) ++m_loops;
            return true;
        }
        if(m_pos >= m_rays.size()) return false;
        Emit(m_rays[m_pos++], fs, 0.0);
        return true;
    }

    // First pass: pull through the wrapped streamer, caching as we go.
    if(!m_inner->TryNext(fs)) {
        if(m_caching && !m_rays.empty()) {
            m_filled = true;
            m_pos = 0;
            m_rays.shrink_to_fit();
            NuGeom::Log().info("FluxStreamer: cached {} rays ({:.1f} MB); later passes are served "
                               "from memory.",
                               m_rays.size(), static_cast<double>(CacheBytes()) / 1048576.0);
            if(m_importance) BuildSampler();
        }
        return false;
    }
    if(m_caching) {
        // File-level constants, captured from the first ray.
        m_window_area = fs.window_area;
        m_pot_per_ray = fs.ray.POT();
        if(m_max_bytes > 0 && (m_rays.size() + 1) * sizeof(Cached) > m_max_bytes) {
            NuGeom::Log().warn("FluxStreamer: cache cap of {:.1f} MB reached after {} rays; "
                               "falling back to streaming for the rest of the run.",
                               static_cast<double>(m_max_bytes) / 1048576.0, m_rays.size());
            m_caching = false;
            m_importance = false;
            m_rays.clear();
            m_rays.shrink_to_fit();
        } else {
            const auto &o = fs.ray.Origin();
            const auto &d = fs.ray.Direction();
            m_rays.push_back(Cached{fs.energy, o.X(), o.Y(), o.Z(), d.X(), d.Y(), d.Z(),
                                    fs.flux_weight, fs.pdg});
        }
    }
    return true;
}

void CachedFluxStreamer::Rewind() {
    if(m_filled) {
        m_pos = 0;
        return;
    }
    // Still on the first pass (or not caching): the inner streamer owns the
    // position.  A rewind here means the first pass ended, so whatever was
    // cached is the whole file.
    if(m_caching && !m_rays.empty() && !m_filled) {
        m_filled = true;
        m_pos = 0;
        if(m_importance) BuildSampler();
    }
    m_inner->Rewind();
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

std::unique_ptr<FluxStreamer> OpenFluxStreamer(const std::string &path, bool loop,
                                               const Vector3D *offset_override,
                                               const CachedFluxStreamer::Options &opts) {
    auto inner = OpenFluxStreamer(path, loop, offset_override);
    if(!opts.cache && !opts.importance) return inner;
    return std::make_unique<CachedFluxStreamer>(std::move(inner), opts);
}

} // namespace NuGeom
