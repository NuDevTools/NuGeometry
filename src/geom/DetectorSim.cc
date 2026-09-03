#include "geom/DetectorSim.hh"
#include "geom/GeneratorInterface.hh"
#include "geom/LineSegment.hh"
#include "geom/Material.hh"
#include "geom/Parser.hh"
#include "geom/Random.hh"
#include <cstdlib>

// HepMC3 and NuHepMC headers trip our strict warnings; silence them here.
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
#include "HepMC3/WriterAscii.h"

#include "NuHepMC/Constants.hxx"
#include "NuHepMC/EventUtils.hxx"
#include "NuHepMC/WriterUtils.hxx"
#pragma GCC diagnostic pop

#include "geom/BoundingBox.hh"
#include "geom/Logging.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>

using NuGeom::DetectorSim;

namespace {

constexpr int kStatusFinal = NuHepMC::ParticleStatus::UndecayedPhysical;
constexpr int kStatusIncoming = NuHepMC::ParticleStatus::IncomingBeam;
constexpr int kStatusTarget = NuHepMC::ParticleStatus::Target;

// Minimal PDG ion code for a nucleus Z,A (1000000000 + Z*10000 + A*10).
int NucleusPDG(const NuGeom::Element &elm) {
    // Element::PDG already returns the ion code for composite nuclei;
    // for protons it returns 2212.  Pass through.
    return static_cast<int>(elm.PDG());
}

} // namespace

DetectorSim::DetectorSim(double safety_factor) : m_safety_factor{safety_factor} {
    // Two-stage layer-1 rejection is exact and on by default; NUGEOM_TWO_STAGE=0
    // forces the single-stage path for cross-checking (same pattern as
    // NUGEOM_TRAVERSAL=sequential).  Read once.
    static const bool disabled = [] {
        const char *e = std::getenv("NUGEOM_TWO_STAGE");
        return e && std::string(e) == "0";
    }();
    if(disabled) m_two_stage = false;
}
DetectorSim::~DetectorSim() = default;

void DetectorSim::SetEventFile(const std::string &outfile) {
    if(!m_run_info) m_run_info = std::make_shared<HepMC3::GenRunInfo>();
    m_writer = std::make_unique<HepMC3::WriterAscii>(outfile, m_run_info);
}

void DetectorSim::Setup(const std::string &geometry) {
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(geometry.c_str());
    if(!result) throw std::runtime_error("GDMLParser: Invalid file");
    NuGeom::GDMLParser parser(doc);
    world = parser.GetWorld();
}

void DetectorSim::SetGenerator(std::shared_ptr<GeneratorInterface> gen) {
    m_generator = std::move(gen);
}

void DetectorSim::SetEnergyRange(double emin, double emax) {
    if(!(emin < emax) || !std::isfinite(emin) || !std::isfinite(emax))
        throw std::invalid_argument(
            "DetectorSim::SetEnergyRange: need a finite emin < emax, got [" + std::to_string(emin) +
            ", " + std::to_string(emax) + "] GeV");
    m_emin = emin;
    m_emax = emax;
}

std::vector<NuGeom::Ray> DetectorSim::GenerateProbeRays(size_t nrays) const {
    const BoundingBox world_box = world.GetWorldBox();
    if(!world_box.IsValid())
        throw std::runtime_error("DetectorSim::Init: world bounding box is invalid");

    std::vector<NuGeom::Ray> rays;
    rays.reserve(nrays);

    // A probe ray is a chord: pick two points, then back the origin out beyond
    // the world so the ray enters from outside exactly as a flux ray does.
    const double reach = (world_box.max - world_box.min).Norm();
    auto add = [&rays, reach](const Vector3D &from, const Vector3D &to) {
        const Vector3D dir = to - from;
        if(dir.Norm() <= 0) return;
        rays.emplace_back(from - dir.Unit() * reach, dir, 1.0);
    };

    // The chords that maximize a box's traversal: the three face-to-face lines
    // through its centre and the four body diagonals.
    auto add_extremes = [&add](const BoundingBox &b) {
        const Vector3D lo = b.min, hi = b.max, mid = (b.min + b.max) * 0.5;
        add({lo.X(), mid.Y(), mid.Z()}, {hi.X(), mid.Y(), mid.Z()});
        add({mid.X(), lo.Y(), mid.Z()}, {mid.X(), hi.Y(), mid.Z()});
        add({mid.X(), mid.Y(), lo.Z()}, {mid.X(), mid.Y(), hi.Z()});
        add(lo, hi);
        add({lo.X(), lo.Y(), hi.Z()}, {hi.X(), hi.Y(), lo.Z()});
        add({lo.X(), hi.Y(), lo.Z()}, {hi.X(), lo.Y(), hi.Z()});
        add({hi.X(), lo.Y(), lo.Z()}, {lo.X(), hi.Y(), hi.Z()});
    };
    auto add_random_chord = [&add](const BoundingBox &b) {
        auto rng = NuGeom::Random::Instance();
        auto point = [&rng, &b]() {
            return Vector3D{rng.Uniform(b.min.X(), b.max.X()), rng.Uniform(b.min.Y(), b.max.Y()),
                            rng.Uniform(b.min.Z(), b.max.Z())};
        };
        add(point(), point());
    };

    // Chords of the world box alone almost never traverse the detectors: on the
    // DUNE ND hall the dense volumes are a small fraction of the hall's volume,
    // and hitting one lengthwise at random is rarer still, so a world-box scan
    // lands well below the column density an ordinary beam ray sees.  Aim the
    // probes at the placed volumes instead, densest first — the maximum column
    // density is by construction a long chord through dense material.
    auto bounds = world.GetVolumeBounds(/*max_depth=*/4);
    auto potential = [](const World::VolumeBounds &vb) {
        double n_total = 0;
        for(const auto &elm : vb.material.Elements()) n_total += vb.material.NumberDensity(elm);
        return n_total * (vb.bb.max - vb.bb.min).Norm();
    };
    std::sort(bounds.begin(), bounds.end(),
              [&potential](const auto &a, const auto &b) { return potential(a) > potential(b); });

    add_extremes(world_box);
    for(const auto &vb : bounds) {
        if(rays.size() + 7 > nrays) break;
        add_extremes(vb.bb);
    }
    // Spend whatever budget is left on random chords through the ranked
    // volumes, so off-centre paths that stack more material still get sampled.
    for(size_t i = 0; rays.size() < nrays; ++i) {
        if(bounds.empty()) {
            add_random_chord(world_box);
            continue;
        }
        add_random_chord(bounds[i % bounds.size()].bb);
    }
    return rays;
}

void DetectorSim::Init(size_t nrays) {
    if(!m_generator) throw std::runtime_error("DetectorSim::Init: generator not set");
    if(!(m_emin < m_emax))
        throw std::runtime_error("DetectorSim::Init: no energy range set; call SetEnergyRange() "
                                 "with the generator's supported beam range (GeV) first, or use "
                                 "InitFromFlux() to calibrate from flux samples instead");

    if(!m_run_info) m_run_info = std::make_shared<HepMC3::GenRunInfo>();
    m_generator->InitRun(m_run_info);

    // Energy grid over the supported range.  The cross section need not be
    // monotone in E (thresholds, resonances), so scan rather than assuming the
    // top of the range wins.  Log-spaced: sigma varies fastest near threshold.
    constexpr size_t kEnergyPoints = 32;
    std::vector<double> energies;
    energies.reserve(kEnergyPoints);
    const bool log_grid = m_emin > 0;
    for(size_t i = 0; i < kEnergyPoints; ++i) {
        const double f = static_cast<double>(i) / (kEnergyPoints - 1);
        energies.push_back(log_grid ? m_emin * std::pow(m_emax / m_emin, f)
                                    : m_emin + f * (m_emax - m_emin));
    }

    const auto rays = GenerateProbeRays(nrays);
    double best = 0;
    double best_e = 0;
    int best_pdg = 0;
    m_max_col_elem.clear();
    for(size_t i = 0; i < rays.size(); ++i) {
        if(i % 100 == 0)
            NuGeom::Log().info("Probe ray {} / {}. Max interaction probability = {}", i,
                               rays.size(), best);
        // One trace per chord, and one walk of its segments: the grid is then
        // scanned against the per-element column densities, so the cost is
        // (grid points x elements) cross-section lookups rather than
        // (grid points x segments x elements).
        const auto cols = ColumnDensity(world.GetLineSegments(rays[i]));
        // Per-element running maximum: the ingredients of GeomBound().
        for(const auto &[elm, n_col] : cols) {
            auto it = m_max_col_elem.find(elm);
            if(it == m_max_col_elem.end())
                m_max_col_elem.emplace(elm, n_col);
            else if(n_col > it->second)
                it->second = n_col;
        }
        for(const double energy : energies) {
            for(const int pdg : m_species) {
                double prob = 0;
                for(const auto &[elm, n_col] : cols) prob += n_col * ElementXSec(energy, pdg, elm);
                if(prob > best) {
                    best = prob;
                    best_e = energy;
                    best_pdg = pdg;
                }
            }
        }
    }

    // A zero or non-finite probability means the generator's cross section
    // vanishes over the whole scan: ProduceEvent would divide by a zero
    // envelope (0/0 -> NaN) and spin forever, so fail loudly with the likely
    // causes instead.
    if(!std::isfinite(best) || best <= 0.0)
        throw std::runtime_error(
            "DetectorSim::Init: max interaction probability = " + std::to_string(best) + " over " +
            std::to_string(rays.size()) + " probe rays and " + std::to_string(energies.size()) +
            " energies in [" + std::to_string(m_emin) + ", " + std::to_string(m_emax) +
            "] GeV — the generator's cross section is zero everywhere in the geometry. "
            "Check that the cross section is non-zero over this energy range (e.g. matching "
            "energy units, GeV vs MeV), that the scanned neutrino species match the flux, and "
            "that the geometry materials map onto the generator's targets.");

    // The flux factor (flux_weight * window_area) is *not* known here and is
    // deliberately not sampled from the flux file: ProduceEvent learns it, and
    // the tighter running maximum of the product, from the rays it throws.
    m_geom_envelope = best * m_safety_factor;
    m_flux_scale = 0.0;
    max_prob = 0.0;
    m_adaptive_envelope = true;

    // ---- GeomBound table -------------------------------------------------
    // sum_e max(n_col_e) * sigma_e(E, pdg): the largest column density of each
    // element combined into one synthetic ray.  No real ray can beat it, and
    // unlike `best` it resolves energy and species instead of collapsing to the
    // worst case over both.  Tabulated per coarse energy bin, each bin holding
    // the max over a finer sub-grid so a lookup stays an upper bound.
    constexpr size_t kSub = 4;
    m_bound_energies = energies;
    m_bound_table.assign(m_species.size(), std::vector<double>(energies.size(), 0.0));
    m_bound_global = 0.0;
    for(size_t si = 0; si < m_species.size(); ++si) {
        const int pdg = m_species[si];
        for(size_t i = 0; i + 1 < energies.size(); ++i) {
            double bin_max = 0.0;
            for(size_t k = 0; k <= kSub; ++k) {
                const double f = static_cast<double>(k) / kSub;
                const double e = energies[i] * std::pow(energies[i + 1] / energies[i], f);
                double p = 0.0;
                for(const auto &[elm, n_col] : m_max_col_elem)
                    p += n_col * ElementXSec(e, pdg, elm);
                bin_max = std::max(bin_max, p);
            }
            m_bound_table[si][i] = bin_max;
            m_bound_global = std::max(m_bound_global, bin_max);
        }
        m_bound_table[si].back() = m_bound_table[si][energies.size() - 2];
    }
    NuGeom::Log().info("GeomBound: summed per-element column maxima over {} elements give a "
                       "worst-case bound of {} (vs {} for the largest single probe ray) -- "
                       "energy-resolved, so most rays are held to a much tighter bound.",
                       m_max_col_elem.size(), m_bound_global, best);

    NuGeom::Log().info("Envelope calibrated on {} probe rays x {} energies x {} species: max "
                       "interaction probability {} (at E = {} GeV, pdg = {}); geometry envelope "
                       "with safety factor {} = {}. The flux-weight scale is learned per ray.",
                       rays.size(), energies.size(), m_species.size(), best, best_e, best_pdg,
                       m_safety_factor, m_geom_envelope);
}

void DetectorSim::InitFromFlux(size_t nrays) {
    if(!m_generator) throw std::runtime_error("DetectorSim::Init: generator not set");
    if(!flux_callback) throw std::runtime_error("DetectorSim::Init: flux callback not set");

    if(!m_run_info) m_run_info = std::make_shared<HepMC3::GenRunInfo>();
    m_generator->InitRun(m_run_info);

    size_t last_update = 0;
    for(size_t i = 0; i < nrays; ++i) {
        if(i % 1000 == 0)
            NuGeom::Log().info("Shot {} / {} rays. Max prob = {}", i, nrays, max_prob);
        FluxSample fs = flux_callback();
        auto [probs, segments] = HandleRay(fs);
        // Layer-1 acceptance weight = flux importance weight * flux-window area *
        // interaction probability.  The flux weight and window area MUST ride in
        // the envelope (and the acceptance below), exactly like GENIE's
        // fMaxWeight over nimpwt*wgt; see ProduceEvent for the correspondence and
        // why window_area is needed.
        double accept_w =
            fs.flux_weight * fs.window_area * std::accumulate(probs.begin(), probs.end(), 0.0);
        if(accept_w > max_prob) {
            NuGeom::Log().debug("Updating max prob: {} -> {}", max_prob, accept_w);
            max_prob = accept_w;
            last_update = i;
        }
        if(i - last_update > 10000) break;
    }

    max_prob *= m_safety_factor;
    // No geometry factor to separate out here, but the envelope still adapts:
    // a warm-up over `nrays` flux samples only ever sees the part of the weight
    // tail it happened to draw, so let ProduceEvent raise it on overflow.
    m_geom_envelope = 0.0;
    m_adaptive_envelope = true;

    // A zero or non-finite max_prob means every sampled ray had a zero
    // acceptance weight (flux_weight * interaction probability): either the
    // generator's cross section vanished over the whole flux, or every ray
    // carried a zero flux weight. ProduceEvent would then divide by it (0/0 ->
    // NaN) and spin forever, so fail loudly with the likely causes instead.
    if(!std::isfinite(max_prob) || max_prob <= 0.0)
        throw std::runtime_error(
            "DetectorSim::InitFromFlux: max_prob = " + std::to_string(max_prob) +
            " after sampling " + std::to_string(nrays) +
            " rays — every ray had zero acceptance weight (flux_weight * interaction prob). "
            "Check that the flux weights are non-zero, that the generator's cross section is "
            "non-zero over the flux energies (e.g. matching energy units, GeV vs MeV), and that "
            "the geometry materials map onto the generator's targets.");

    NuGeom::Log().info("Maximum probability found with safety factor ({}): {}", m_safety_factor,
                       max_prob);
}

bool DetectorSim::ProduceEvent() {
    FluxSample fs = flux_callback();

    // Layer-1 acceptance weight = flux importance weight * flux-window area *
    // interaction probability.  This mirrors GENIE GNuMIFlux's unweighted/retry
    // path:
    //   * the flux weight (dk2nu nimpwt*wgt) rides in the acceptance, accepting
    //     a ray with probability accept_w / max_prob — GENIE's w_i / fMaxWeight;
    //   * window_area [cm^2] turns the per-cm^2 flux weight into a rate over the
    //     detector face.  The converter spreads the ray origins uniformly across
    //     a transverse flux window of this area (GENIE samples positions over its
    //     flux window for the same reason), so summing flux_weight*window_area*
    //     tot_probs Monte-Carlo integrates the rate over the face.  Omitting it
    //     (all rays converging to one point) leaves the rate low by ~window_area
    //     (~2.4e5 cm^2 for ND-LAr);
    //   * max_prob is the envelope max_i(flux_weight_i * window_area * tot_probs_i),
    //     GENIE's fMaxWeight (here it also folds the interaction prob, since
    //     DetectorSim couples flux and interaction in one rejection step).
    const double flux_scale = fs.flux_weight * fs.window_area;
    ++m_stats.thrown;

    // The envelope in force for THIS throw.  It is deliberately a function of
    // the rays already seen and NOT of the current one: making max_prob track
    // the current ray's flux weight would cancel that weight out of the
    // acceptance (accept_w/max_prob loses its w) and move it into the POT
    // charge (POT/max_prob picks up a spurious 1/w), which is exactly the
    // normalization error the flux-weight tests guard against.  With a
    // past-only envelope, events and POT both scale as 1/max_prob for the same
    // throw, so the per-ray rate accept_w/POT_ray is envelope-independent and
    // the envelope may move freely between throws.
    if(m_geom_envelope > 0.0) {
        // Ceiling from Init()'s probe rays x the flux weights seen so far, and
        // the much tighter running maximum of the product.  The ceiling alone
        // overshoots by the factor between max(w)*max(P) and max(w*P) -- ~45x
        // on the DUNE ND hall -- because the heaviest ray and the longest path
        // are different rays; the running maximum recovers that, and the
        // ceiling keeps it honest before enough rays have been seen.
        const double ceiling = m_geom_envelope * m_flux_scale;
        const double running = m_joint_max * m_safety_factor;
        if(running > 0.0)
            max_prob = ceiling > 0.0 ? std::min(ceiling, running) : running;
        else if(ceiling > 0.0)
            max_prob = ceiling;
    }
    const double envelope = max_prob;

    // The flux half of the envelope needs no traversal, so it is updated from
    // EVERY ray -- including the ones stage 1 rejects below.  That keeps the
    // `ceiling` bound valid for the whole flux, which is what makes the stage-1
    // pre-rejection safe.
    if(m_geom_envelope > 0.0 && flux_scale > m_flux_scale) m_flux_scale = flux_scale;

    // Burn-in rays seed the envelope (both halves) and are never accepted, so
    // they always traverse -- there are only m_envelope_burn_in of them.
    const bool in_burn_in = (m_geom_envelope > 0.0 && m_stats.thrown <= m_envelope_burn_in);

    // --- Stage 1: reject on the flux weight alone, before any traversal ------
    //
    // accept_w = flux_scale * total_prob and total_prob <= m_geom_envelope (the
    // probe-ray calibration), so the acceptance probability
    //     p = accept_w / envelope
    // is bounded by
    //     q = flux_scale * m_geom_envelope / envelope,
    // which depends only on this ray's flux weight.  Throwing against q here
    // and against p/q = total_prob/m_geom_envelope after the traversal gives
    // exactly p, so this is an exact restructuring, not an approximation.
    //
    // A ray that fails stage 1 has accept_w <= flux_scale*m_geom_envelope
    // < envelope, so it could never have been accepted AND could never have
    // clipped -- only its (traversal-free) POT charge and flux scale matter.
    const bool two_stage = m_two_stage && !in_burn_in && m_geom_envelope > 0.0 && envelope > 0.0 &&
                           std::isfinite(envelope);
    double q = 1.0;
    double bound = 0.0;
    if(two_stage) {
        // Energy- and species-resolved bound on total_prob (see GeomBound):
        // tighter than the flat m_geom_envelope, so stage 1 rejects more rays
        // before paying for a traversal.
        bound = GeomBound(fs.energy, fs.pdg);
        if(!(bound > 0.0) || !std::isfinite(bound)) bound = m_geom_envelope;
        q = flux_scale * bound / envelope;
        if(q < 1.0) {
            if(NuGeom::Random::Instance().Uniform(0.0, 1.0) > q) {
                m_pot += fs.ray.POT() / envelope;
                m_stats.pot = m_pot;
                m_stats.max_prob = envelope;
                ++m_stats.stage1_rejected;
                return false;
            }
        } else {
            q = 1.0; // stage 1 cannot help this ray; stage 2 carries the full p
        }
    }

    // --- Geometry traversal (the expensive part) -----------------------------
    // Layer 1 needs only the total interaction probability, so take the
    // column-density-only traversal here.  The full LineSegment list is built
    // below, once, and only for rays that survive the acceptance throw.
    const double total_prob = ColumnProb(fs.energy, fs.pdg, fs.ray);
    const double accept_w = flux_scale * total_prob;
    ++m_stats.traversed;

    // Adaptive half of the bound: if a real ray beats the tabulated maximum,
    // raise the inflation so subsequent rays are held to a bound that would
    // have covered it.  Past-only, monotone, and it converges -- the summed
    // per-element maxima make this rare by construction.
    if(two_stage && bound > 0.0 && total_prob > bound) {
        m_bound_inflation *= total_prob / bound;
        ++m_bound_raises;
        NuGeom::Log().debug("GeomBound exceeded ({} > {}); inflation now {}", total_prob, bound,
                            m_bound_inflation);
    }

    // Fold this ray into the joint envelope for SUBSEQUENT throws.
    if(m_geom_envelope > 0.0) {
        if(accept_w > m_joint_max) {
            NuGeom::Log().debug("ProduceEvent: raising the running envelope {} -> {}", m_joint_max,
                                accept_w);
            m_joint_max = accept_w;
            ++m_stats.envelope_growths;
        }
    } else if(m_adaptive_envelope && accept_w > max_prob) {
        // Fixed envelope (SetMaxProb / InitFromFlux): raise it for the next
        // throw so a warm-up that missed the tail self-corrects.
        max_prob = accept_w * m_safety_factor;
        ++m_stats.envelope_growths;
    }

    // Burn-in: the first rays after a custom-ray Init() only establish the
    // envelope above.  Dropping them keeps the accounting clean -- a ray thrown
    // against a barely-seeded envelope charges POT/envelope, an enormous
    // exposure, while producing at most one event, which would dilute
    // events/POT by several percent for the whole run.
    if(in_burn_in) return false;

    // Bootstrap: with no envelope yet there is nothing to throw against.
    if(!(envelope > 0.0) || !std::isfinite(envelope)) {
        NuGeom::Log().debug("ProduceEvent: envelope not yet defined (accept_w = {}); "
                            "seeding it from this ray and skipping the throw",
                            accept_w);
        return false;
    }

    // Clipping is the only way the envelope can bias the rate: a ray with
    // accept_w > envelope is accepted with probability 1 instead of
    // accept_w/envelope, so it under-counts.  With the geometry factor pinned
    // by Init() this can only happen on a new record flux weight — rare, and
    // bounded — but it is accounted for here so the bias is measurable rather
    // than silent.
    // The probability actually thrown against.  Without two-stage this is
    // accept_w/envelope as before; with it, stage 1 has already consumed a
    // factor q, so stage 2 must carry p/q = total_prob/m_geom_envelope.
    const double p_full = accept_w / envelope;
    const double p_stage2 = q < 1.0 ? p_full / q : p_full;

    if(p_stage2 > 1.0) {
        ++m_stats.clipped;
        m_stats.clipped_excess += p_stage2 - 1.0;
        NuGeom::Log().debug("ProduceEvent: acceptance probability {} exceeds 1 (clip #{}); "
                            "accept_w = {}, envelope = {}, q = {}",
                            p_stage2, m_stats.clipped, accept_w, envelope, q);
    }

    // POT charged per thrown ray, flat and independent of accept/reject — the
    // analogue of GENIE's "fAccumPOTs += fEffPOTsPerNu / fMaxWeight".  ray.POT()
    // is filePOT/Nrays (fEffPOTsPerNu); dividing by the envelope makes m_pot the
    // *effective/used* POT (GENIE UsedPOTs), so that the unweighted events
    // satisfy events/POT == Σ w_i*tot_probs_i / POT_total (the true rate).  The
    // flux weight must NOT divide POT per ray (that gives a spurious Σ 1/w_i) —
    // it belongs in accept_w above.
    m_pot += fs.ray.POT() / envelope;
    m_stats.pot = m_pot;
    m_stats.max_prob = envelope;

    NuGeom::Log().trace("ProduceEvent: E = {} GeV, pdg = {}, P = {} w = {} "
                        "/ envelope = {} (accept prob {})",
                        fs.energy, fs.pdg, total_prob, fs.flux_weight, envelope,
                        accept_w / envelope);

    if(p_stage2 < NuGeom::Random::Instance().Uniform(0.0, 1.0)) return false;
    ++m_stats.accepted;

    // Accepted: now pay for the full segment list, which PickVertex needs to
    // place the vertex.  Draw against this traversal's own total so a last-bit
    // difference in summation order cannot walk off the end of `probs`.
    auto [probs, segments] = HandleRay(fs);
    const double seg_total = std::accumulate(probs.begin(), probs.end(), 0.0);
    double r = NuGeom::Random::Instance().Uniform(0.0, seg_total);
    VertexPick vp = PickVertex(segments, probs, r);

    Element target = PickElement(vp.material, fs.energy, fs.pdg);
    NuGeom::Log().debug("ProduceEvent: vertex at ({}, {}, {}) in segment {} ({}), target = {}",
                        vp.position.X(), vp.position.Y(), vp.position.Z(), vp.segment_idx,
                        vp.material.Name(), target.PDG());

    // Build the HepMC event skeleton.
    HepMC3::GenEvent evt(m_run_info, HepMC3::Units::GEV, HepMC3::Units::CM);
    evt.set_event_number(static_cast<int>(m_event_number++));

    // Incoming neutrino four-momentum: treat as massless, direction from ray.
    auto dir = fs.ray.Direction();
    HepMC3::FourVector nu_p4(fs.energy * dir.X(), fs.energy * dir.Y(), fs.energy * dir.Z(),
                             fs.energy);

    auto nu = std::make_shared<HepMC3::GenParticle>(nu_p4, fs.pdg, kStatusIncoming);

    auto vtx = std::make_shared<HepMC3::GenVertex>(
        HepMC3::FourVector(vp.position.X(), vp.position.Y(), vp.position.Z(), 0.0));
    vtx->add_particle_in(nu);
    evt.add_vertex(vtx);

    // Communicate the target nucleus to the generator by PDG only. We do NOT
    // attach a target particle to the vertex: where the struck nucleon sits
    // inside the nucleus (and the nucleus -> nucleon separation) is the
    // generator's internal detail, which it fills onto this event.
    evt.add_attribute("NuGeom.TargetPDG",
                      std::make_shared<HepMC3::IntAttribute>(NucleusPDG(target)));

    // NuGeom-side attributes for audit / POT reconstruction.
    vtx->add_attribute("NuGeom.material",
                       std::make_shared<HepMC3::StringAttribute>(vp.material.Name()));
    evt.add_attribute("NuGeom.Flux.Weight",
                      std::make_shared<HepMC3::DoubleAttribute>(fs.flux_weight));
    evt.add_attribute("NuGeom.Flux.WindowArea_cm2",
                      std::make_shared<HepMC3::DoubleAttribute>(fs.window_area));
    evt.add_attribute("NuGeom.Flux.POT", std::make_shared<HepMC3::DoubleAttribute>(fs.ray.POT()));
    evt.add_attribute("NuGeom.MaxProbNorm",
                      std::make_shared<HepMC3::DoubleAttribute>(1.0 / envelope));

    if(fs.decorate) fs.decorate(evt);

    bool emitted = false;
    if(m_mode == SamplingMode::TotalXSecRetry) {
        // sigma-based layer-1 envelope: generator must produce exactly one
        // event per accepted vertex to keep the per-process rate unbiased.
        // Safety cap so a badly-tuned unweighter can't hang the run.
        constexpr size_t kMaxTries = 100000;
        size_t tries = 0;
        for(; tries < kMaxTries && !emitted; ++tries) emitted = m_generator->GenerateEvent(evt);
        if(!emitted)
            NuGeom::Log().warn("Generator failed to emit after {} retries", kMaxTries);
        else if(tries > 1)
            NuGeom::Log().trace("ProduceEvent: generator emitted after {} trials", tries);
    } else {
        // max_w-based envelope absorbs the unweighter inefficiency into
        // layer 1; a single-shot rejection is physically correct here.
        emitted = m_generator->GenerateEvent(evt);
        if(!emitted) NuGeom::Log().trace("ProduceEvent: generator rejected the single-shot trial");
    }

    if(emitted) {
        // Track the E.C.1 weight the generator attached.  A partial unweighter
        // emits max(1, raw_w/max_w), so <weight> > 1 means part of the rate is
        // carried by the weights rather than by the event count -- the factor
        // by which raw event counts disagree between the two sampling modes.
        // See ReportRunSummary().
        ++m_stats.emitted;
        m_stats.sum_weight += evt.weights().empty() ? 1.0 : evt.weights().front();
        if(m_writer) m_writer->write_event(evt);
    }
    return emitted;
}

void DetectorSim::ReportRunSummary() const {
    auto &log = NuGeom::Log();
    if(m_stats.emitted == 0 || m_pot <= 0.0) {
        log.warn("Run summary: {} events over {} thrown rays, POT = {}", m_stats.emitted,
                 m_stats.thrown, m_pot);
        return;
    }

    const double n = static_cast<double>(m_stats.emitted);
    const double count_rate = n / m_pot;
    const double weighted_rate = m_stats.sum_weight / m_pot;
    const double mean_w = m_stats.sum_weight / n;

    log.info("Run summary ({} mode): {} events / {} thrown rays, POT = {:.6e}",
             m_mode == SamplingMode::TotalXSecRetry ? "TotalXSecRetry" : "EnvelopeNoRetry",
             m_stats.emitted, m_stats.thrown, m_pot);
    log.info("  events/POT       = {:.6e}", count_rate);
    log.info("  sum(weights)/POT = {:.6e}", weighted_rate);
    log.info("  <weight>         = {:.6f}", mean_w);
    log.info("  envelope         = {:.6e} ({} growths)", max_prob, m_stats.envelope_growths);
    if(m_stats.stage1_rejected > 0) {
        const double frac =
            static_cast<double>(m_stats.traversed) / static_cast<double>(m_stats.thrown);
        log.info("  traversals       = {} of {} rays ({:.4f}%, 1 in {:.0f}) -- two-stage "
                 "rejection skipped the geometry for {} rays",
                 m_stats.traversed, m_stats.thrown, 100.0 * frac, frac > 0 ? 1.0 / frac : 0.0,
                 m_stats.stage1_rejected);
        if(m_bound_raises > 0)
            log.info("  GeomBound        = raised {} times, inflation {:.4f} (a real ray beat the "
                     "summed per-element column maxima)",
                     m_bound_raises, m_bound_inflation);
    }

    // Clipped rays are the one envelope-induced bias, and it is always a
    // shortfall: `clipped_excess` is the number of extra layer-1 accepts they
    // should have produced beyond the one each was capped at.  It is a
    // *layer-1* deficit, so it must be compared against the layer-1 accept
    // count -- not against the emitted events, which are fewer by the
    // generator's single-shot efficiency and would inflate the quoted bias.
    if(m_stats.clipped > 0 && m_stats.accepted > 0)
        log.warn("Run summary: {} of {} rays ({:.3f}%) exceeded the envelope and were clipped, "
                 "costing an estimated {:.1f} of {} layer-1 accepts ({:.3f}% of the rate). "
                 "Raise the safety factor if this is not negligible.",
                 m_stats.clipped, m_stats.thrown,
                 100.0 * static_cast<double>(m_stats.clipped) / static_cast<double>(m_stats.thrown),
                 m_stats.clipped_excess, m_stats.accepted,
                 100.0 * m_stats.clipped_excess / static_cast<double>(m_stats.accepted));

    // The two sampling modes only agree on the *physical* rate, and which
    // estimator carries it depends on the mode.  With a strict unweighter
    // (every emitted weight exactly 1) the two estimators coincide and the
    // modes' raw event counts agree.  With a partial unweighter -- e.g.
    // Achilles' Percentile unweighter, whose max_w is a percentile of the
    // weight distribution and therefore leaves an overweight tail -- they do
    // not, and the two modes' raw counts differ by exactly this <weight>:
    //
    //   EnvelopeNoRetry: layer 1 uses sigma_env = max_w and the single-shot
    //     trial emits with probability E[min(1, w/max_w)] = sigma/(max_w
    //     <weight>).  sum(weights)/POT is the physical rate; the raw count is
    //     LOW by a factor <weight>.
    //   TotalXSecRetry: layer 1 uses sigma_tot and retries until emit, so
    //     there is exactly one event per accepted vertex.  events/POT is the
    //     physical rate; sum(weights)/POT is HIGH by a factor <weight>.
    //
    // Comparing raw event counts between the modes therefore shows a
    // systematic <weight> discrepancy that is not a statistical fluctuation.
    constexpr double kWeightTolerance = 0.02;
    if(std::abs(mean_w - 1.0) > kWeightTolerance)
        log.warn(
            "Run summary: <weight> = {:.4f} != 1, so the generator's unweighter left an "
            "overweight tail carrying {:.1f}% of the cross section. events/POT and "
            "sum(weights)/POT therefore differ by that factor, and the two sampling modes will "
            "NOT produce the same number of events for the same POT (they differ by exactly "
            "{:.4f}). Use sum(weights)/POT in EnvelopeNoRetry and events/POT in TotalXSecRetry, "
            "or make the generator's envelope a true bound (e.g. Achilles "
            "Options/Unweighting/percentile: 100) to make the raw counts agree.",
            mean_w, 100.0 * (1.0 - 1.0 / mean_w), mean_w);
}

std::set<NuGeom::Material> DetectorSim::GetMaterials(const LineSegments &segments) const {
    std::set<NuGeom::Material> mats;
    for(const auto &segment : segments) mats.insert(segment.GetMaterial());
    return mats;
}

std::vector<double> DetectorSim::EvaluateProbs(const LineSegments &segments,
                                               const std::map<NuGeom::Element, double> &xsecsmaps) {
    std::vector<double> probs;
    probs.reserve(segments.size());
    for(const auto &segment : segments) {
        const auto &material = segment.GetMaterial();
        double inv_mfp = 0;
        for(const auto &elm : material.Elements()) {
            auto cross_section = xsecsmaps.at(elm);
            inv_mfp += cross_section * material.NumberDensity(elm);
        }
        probs.push_back(segment.Length() * inv_mfp);
    }
    double total_prob = std::accumulate(probs.begin(), probs.end(), 0.0);
    if(total_prob > max_prob) max_prob = total_prob;
    return probs;
}

void DetectorSim::GenerateEvents(size_t nevents) {
    if(!m_generator) throw std::runtime_error("DetectorSim::GenerateEvents: generator not set");
    size_t accepted = 0;
    while(accepted < nevents) {
        // Report only when an event is actually accepted; checking outside the
        // branch would re-print the same count for every rejected ray.
        if(ProduceEvent()) {
            ++accepted;
            if(accepted % 100 == 0)
                NuGeom::Log().info("Generated {} / {} events", accepted, nevents);
        }
    }
    ReportRunSummary();
    m_generator->FinalizeRun();
    if(m_writer) m_writer->close();
}

void DetectorSim::GenerateEvents(double pot) {
    if(!m_generator) throw std::runtime_error("DetectorSim::GenerateEvents: generator not set");
    size_t nhits = 0;
    while(m_pot < pot) {
        if(ProduceEvent()) {
            ++nhits;
            if(nhits % 100 == 0)
                NuGeom::Log().info("Generated {} events; POT = {:.3e} / {:.3e}", nhits, m_pot, pot);
        }
    }
    NuGeom::Log().info("Accumulated {} events with {} POT", nhits, m_pot);
    ReportRunSummary();
    m_generator->FinalizeRun();
    if(m_writer) m_writer->close();
}

std::pair<NuGeom::Vector3D, NuGeom::Material>
DetectorSim::Interaction(const LineSegments &segments,
                         const std::map<NuGeom::Element, double> &xsecsmaps) {
    auto prob = EvaluateProbs(segments, xsecsmaps);
    double r1 = NuGeom::Random::Instance().Uniform(0.0, 1.0);
    double total_prob = std::accumulate(prob.begin(), prob.end(), 0.0);
    if(total_prob / max_prob < r1) return {NuGeom::Vector3D(9e9, 9e9, 9e9), {}};
    double r2 = NuGeom::Random::Instance().Uniform(0.0, total_prob);
    VertexPick vp = PickVertex(segments, prob, r2);
    return {vp.position, vp.material};
}

NuGeom::DetectorSim::VertexPick DetectorSim::PickVertex(const LineSegments &segments,
                                                        const std::vector<double> &probs,
                                                        double r) const {
    double cumulative = 0.0;
    size_t idx = probs.size() - 1;
    for(size_t i = 0; i < probs.size(); ++i) {
        cumulative += probs[i];
        if(cumulative >= r) {
            idx = i;
            break;
        }
    }
    // Fraction of the way into segment idx: (r - cumulative_before_idx) / probs[idx].
    double into = (r - (cumulative - probs[idx])) / probs[idx];
    Vector3D pos = segments[idx].Start() + (segments[idx].End() - segments[idx].Start()) * into;
    return {pos, segments[idx].GetMaterial(), idx};
}

NuGeom::Element DetectorSim::PickElement(const Material &mat, double energy, int nu_pdg) const {
    const auto &elms = mat.Elements();
    if(elms.size() == 1) return elms.front();

    // Weight by n_e * sigma_e (the same mode-appropriate cross section used for
    // the mean free path), so the chosen element is distributed by interaction
    // rate rather than by abundance alone.
    double total = 0;
    std::vector<double> weights;
    weights.reserve(elms.size());
    for(const auto &elm : elms) {
        double w = mat.NumberDensity(elm) * ElementXSec(energy, nu_pdg, elm);
        weights.push_back(w);
        total += w;
    }
    double r = NuGeom::Random::Instance().Uniform(0.0, total);
    double cum = 0;
    for(size_t i = 0; i < elms.size(); ++i) {
        cum += weights[i];
        if(cum >= r) return elms[i];
    }
    return elms.back();
}

std::map<NuGeom::Element, double> DetectorSim::ColumnDensity(const LineSegments &segments) const {
    std::map<NuGeom::Element, double> cols;
    for(const auto &segment : segments) {
        const auto &material = segment.GetMaterial();
        const double length = segment.Length();
        for(const auto &elm : material.Elements())
            cols[elm] += length * material.NumberDensity(elm);
    }
    return cols;
}

NuGeom::HandledRay DetectorSim::HandleRay(const FluxSample &fs) const {
    return HandleRay(fs.energy, fs.pdg, fs.ray);
}

NuGeom::HandledRay DetectorSim::HandleRay(double energy, int nu_pdg, const NuGeom::Ray &ray) const {
    LineSegments segments = world.GetLineSegments(ray);
    std::vector<double> probs;
    probs.reserve(segments.size());
    for(const auto &segment : segments) {
        const auto &material = segment.GetMaterial();
        auto meanfreepath = CalculateMeanFreePath(energy, nu_pdg, material);
        probs.push_back(segment.Length() / meanfreepath);
    }
    return {std::move(probs), std::move(segments)};
}

double DetectorSim::GeomBound(double energy, int nu_pdg) const {
    if(m_bound_table.empty() || m_bound_energies.size() < 2)
        return m_geom_envelope * m_bound_inflation;
    // Outside the calibrated range the table says nothing: fall back to the
    // global maximum rather than extrapolating into an under-estimate.
    if(energy < m_bound_energies.front() || energy > m_bound_energies.back())
        return m_bound_global * m_safety_factor * m_bound_inflation;
    size_t si = 0;
    for(size_t i = 0; i < m_species.size(); ++i)
        if(m_species[i] == nu_pdg) {
            si = i;
            break;
        }
    // Log-spaced grid: the bin index is a closed form.
    const double lo = m_bound_energies.front(), hi = m_bound_energies.back();
    const size_t nbin = m_bound_energies.size() - 1;
    const double f = std::log(energy / lo) / std::log(hi / lo);
    size_t i = static_cast<size_t>(f * static_cast<double>(nbin));
    if(i >= nbin) i = nbin - 1;
    return m_bound_table[si][i] * m_safety_factor * m_bound_inflation;
}

double DetectorSim::ColumnProb(double energy, int nu_pdg, const NuGeom::Ray &ray) const {
    double total = 0.0;
    for(const auto &[mat, len] : world.GetColumnLengths(ray))
        total += len / CalculateMeanFreePath(energy, nu_pdg, *mat);
    return total;
}

double DetectorSim::ElementXSec(double energy, int nu_pdg, const NuGeom::Element &elm) const {
    if(m_generator)
        return (m_mode == SamplingMode::TotalXSecRetry)
                   ? m_generator->TotalXSec(energy, nu_pdg, static_cast<int>(elm.PDG()))
                   : m_generator->EnvelopeXSec(energy, nu_pdg, static_cast<int>(elm.PDG()));
    if(xsec_callback) return xsec_callback(energy, elm.PDG());
    return 0;
}

double DetectorSim::CalculateMeanFreePath(double energy, int nu_pdg,
                                          const NuGeom::Material &mat) const {
    // This loop runs per element, per segment, per ray. Hoist the logger out
    // and gate the trace on the level so disabled traces cost nothing (no
    // per-element format-argument capture). See profile-2026-06 notes.
    auto &log = NuGeom::Log();
    const bool trace_on = log.should_log(spdlog::level::trace);
    double inv_mfp = 0;
    for(const auto &elm : mat.Elements()) {
        if(trace_on) log.trace("Element: {}", elm.PDG());
        inv_mfp += mat.NumberDensity(elm) * ElementXSec(energy, nu_pdg, elm);
    }
    return 1.0 / inv_mfp;
}
