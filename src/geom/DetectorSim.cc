#include "geom/DetectorSim.hh"
#include "geom/GeneratorInterface.hh"
#include "geom/LineSegment.hh"
#include "geom/Material.hh"
#include "geom/Parser.hh"
#include "geom/Random.hh"

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

#include "geom/Logging.hh"
#include <numeric>

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

DetectorSim::DetectorSim(double safety_factor) : m_safety_factor{safety_factor} {}
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

void DetectorSim::Init(size_t nrays) {
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
        double prob_tot = std::accumulate(probs.begin(), probs.end(), 0.0) * fs.flux_weight;
        if(prob_tot > max_prob) {
            NuGeom::Log().debug("Updating max prob: {} -> {}", max_prob, prob_tot);
            max_prob = prob_tot;
            last_update = i;
        }
        if(i - last_update > 10000) break;
    }

    max_prob *= m_safety_factor;
    NuGeom::Log().info("Maximum probability found with safety factor ({}): {}", m_safety_factor,
                       max_prob);
}

bool DetectorSim::ProduceEvent() {
    FluxSample fs = flux_callback();

    // Charge POT per thrown ray, independent of accept/reject outcome.
    const double pot_contribution = fs.flux_weight * fs.ray.POT() / max_prob;
    m_pot += pot_contribution;

    auto [probs, segments] = HandleRay(fs);
    double tot_probs = std::accumulate(probs.begin(), probs.end(), 0.0) * fs.flux_weight;

    NuGeom::Log().trace("ProduceEvent: E = {} GeV, pdg = {}, {} segments, P = {} / max_prob = {} "
                        "(accept prob {})",
                        fs.energy, fs.pdg, segments.size(), tot_probs, max_prob,
                        tot_probs / max_prob);
    if(tot_probs > max_prob)
        NuGeom::Log().warn(
            "ProduceEvent: ray probability {} exceeds max_prob {}; increase the safety "
            "factor or Init() sample size",
            tot_probs, max_prob);

    if(tot_probs / max_prob < NuGeom::Random::Instance().Uniform(0.0, 1.0)) return false;

    double total_prob = std::accumulate(probs.begin(), probs.end(), 0.0);
    double r = NuGeom::Random::Instance().Uniform(0.0, total_prob);
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
    evt.add_attribute("NuGeom.Flux.POT", std::make_shared<HepMC3::DoubleAttribute>(fs.ray.POT()));
    evt.add_attribute("NuGeom.MaxProbNorm",
                      std::make_shared<HepMC3::DoubleAttribute>(1.0 / max_prob));

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

    if(emitted && m_writer) m_writer->write_event(evt);
    return emitted;
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
        auto material = segment.GetMaterial();
        double inv_mfp = 0;
        for(auto elm : material.Elements()) {
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
                NuGeom::Log().info("Generated {} events; POT = {} / {}", nhits, m_pot, pot);
        }
    }
    NuGeom::Log().info("Accumulated {} events with {} POT", nhits, m_pot);
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

NuGeom::HandledRay DetectorSim::HandleRay(const FluxSample &fs) const {
    return HandleRay(fs.energy, fs.pdg, fs.ray);
}

NuGeom::HandledRay DetectorSim::HandleRay(double energy, int nu_pdg, const NuGeom::Ray &ray) const {
    LineSegments segments = world.GetLineSegments(ray);
    std::vector<double> probs;
    probs.reserve(segments.size());
    for(const auto &segment : segments) {
        auto material = segment.GetMaterial();
        auto meanfreepath = CalculateMeanFreePath(energy, nu_pdg, material);
        probs.push_back(segment.Length() / meanfreepath);
    }
    return {probs, segments};
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
    double inv_mfp = 0;
    for(const auto &elm : mat.Elements()) {
        NuGeom::Log().trace("Element: {}", elm.PDG());
        inv_mfp += mat.NumberDensity(elm) * ElementXSec(energy, nu_pdg, elm);
    }
    return 1.0 / inv_mfp;
}
