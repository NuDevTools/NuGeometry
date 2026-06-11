#include "adapters/achilles/AchillesAdapter.hh"

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

// Achilles.
#include "Achilles/Constants.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventGen.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/ParticleInfo.hh" // PID
// FillHepMC3Event merges the generated Achilles event onto NuGeom's GenEvent.
#include "Achilles/Writers/NuHepMCWriter.hh"
#pragma GCC diagnostic pop

#include "geom/Logging.hh"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace NuGeom {

AchillesAdapter::AchillesAdapter(const std::string &yaml_config) : m_config_path{yaml_config} {}

AchillesAdapter::AchillesAdapter(const YAML::Node &config, std::string run_card_name)
    : m_config_path{std::move(run_card_name)}, m_config_node{config} {}

AchillesAdapter::~AchillesAdapter() = default;

YAML::Node &AchillesAdapter::Config() {
    // Constructed from a path: load it so the node can be edited.  InitRun
    // then takes the node branch; achilles::Settings resolves any !include
    // tags either way.
    if(!m_config_node) m_config_node = YAML::LoadFile(m_config_path);
    return *m_config_node;
}

std::optional<std::pair<double, double>> AchillesAdapter::BeamEnergyRange() {
    const YAML::Node range = Config()["Beam"]["EnergyRange"];
    if(!range) return std::nullopt;
    return std::make_pair(range[0].as<double>(), range[1].as<double>());
}

void AchillesAdapter::SetBeamEnergyRange(double emin, double emax) {
    if(m_evgen)
        throw std::logic_error("AchillesAdapter: SetBeamEnergyRange must be called before InitRun");
    if(!(emin < emax))
        throw std::invalid_argument("AchillesAdapter: SetBeamEnergyRange needs emin < emax, got [" +
                                    std::to_string(emin) + ", " + std::to_string(emax) + "] MeV");
    Config()["Beam"]["EnergyRange"] = std::vector<double>{emin, emax};
}

bool AchillesAdapter::EnsureBeamEnergyCoverage(double emin, double emax) {
    const auto current = BeamEnergyRange();
    if(current && current->first <= emin && current->second >= emax) {
        NuGeom::Log().debug(
            "[adapter] Beam/EnergyRange [{}, {}] MeV already covers the flux [{}, {}] MeV",
            current->first, current->second, emin, emax);
        return false;
    }

    double lo = emin, hi = emax;
    if(current) {
        lo = std::min(current->first, lo);
        hi = std::max(current->second, hi);
        NuGeom::Log().warn("[adapter] widening Beam/EnergyRange [{}, {}] -> [{}, {}] MeV to "
                           "cover the flux energies [{}, {}] MeV",
                           current->first, current->second, lo, hi, emin, emax);
    } else {
        NuGeom::Log().info("[adapter] setting Beam/EnergyRange from the flux: [{}, {}] MeV", lo,
                           hi);
    }

    // GeometryBeam divides by (max - min); give a monoenergetic flux a
    // finite support.
    if(hi - lo < 1.0) {
        lo = std::max(0.0, lo - 0.5);
        hi = lo + 1.0;
    }

    SetBeamEnergyRange(lo, hi);
    return true;
}

void AchillesAdapter::InitRun(std::shared_ptr<HepMC3::GenRunInfo> run_info) {
    // Achilles does everything it normally does for a run.
    if(m_config_node)
        m_evgen = std::make_unique<achilles::EventGen>(
            *m_config_node, /*plugins=*/std::vector<std::string>{}, m_config_path);
    else
        m_evgen = std::make_unique<achilles::EventGen>(m_config_path,
                                                       /*plugins=*/std::vector<std::string>{});
    m_evgen->Initialize(); // MultiChannel warm-up, unweighter calibration, splines

    // Match Achilles' per-event weight convention to the DetectorSim mode.
    m_evgen->SetSamplingMode(m_total_xsec_retry
                                 ? achilles::EventGen::SamplingMode::TotalXSecRetry
                                 : achilles::EventGen::SamplingMode::EnvelopeNoRetry);

    // Extend NuGeometry's run_info with the Achilles GR.* run-level metadata
    // (version, tools, process list, status codes, conventions, units).
    m_evgen->WriteRunInfo(std::move(run_info));
}

double AchillesAdapter::TotalXSec(double energy, int nu_pdg, int target_pdg) {
    // Cross-section spline value at this (E, nu, target).  Used when
    // DetectorSim is in TotalXSecRetry mode; DetectorSim then retries
    // GenerateEvent until emit per accepted vertex.
    return m_evgen->TotalXSec(energy, achilles::PID{nu_pdg}, achilles::PID{target_pdg});
}

double AchillesAdapter::EnvelopeXSec(double energy, int nu_pdg, int target_pdg) {
    // Unweighter max_w for the (nu, target) process group.  Achilles'
    // integrand weight already absorbs V_ps as the Jacobian, so max_w has
    // the units of sigma_tot and is the correct envelope as-is.
    (void)energy; // max_w is a single post-warmup scalar per ProcessGroup.
    double v = m_evgen->VertexEnvelope(achilles::PID{nu_pdg}, achilles::PID{target_pdg});
    NuGeom::Log().debug("[adapter] EnvelopeXSec(E={}, nu={}, tgt={}) = {}", energy, nu_pdg,
                        target_pdg, v);
    return v;
}

bool AchillesAdapter::GenerateEvent(HepMC3::GenEvent &evt) {
    // NuGeom pre-filled the primary vertex with the incoming neutrino, and the
    // target nucleus PDG as an event attribute.
    auto vtxs = evt.vertices();
    if(vtxs.empty()) throw std::runtime_error("AchillesAdapter: no vertex");
    auto ins = vtxs.front()->particles_in();
    if(ins.empty()) throw std::runtime_error("AchillesAdapter: vertex missing incoming neutrino");
    auto nu_in = ins[0];

    auto tgt_attr = evt.attribute<HepMC3::IntAttribute>("NuGeom.TargetPDG");
    if(!tgt_attr) throw std::runtime_error("AchillesAdapter: event missing NuGeom.TargetPDG");

    // NuGeom builds momenta in the event's units (GeV); Achilles works in MeV.
    const double to_mev = (evt.momentum_unit() == HepMC3::Units::GEV) ? 1000.0 : 1.0;
    achilles::FourVector nu_p4{nu_in->momentum().e() * to_mev, nu_in->momentum().px() * to_mev,
                               nu_in->momentum().py() * to_mev, nu_in->momentum().pz() * to_mev};
    // Note: NuGeom-side flux_weight is intentionally NOT forwarded here. It is
    // already folded into layer-1 acceptance and POT accounting in DetectorSim;
    // Achilles' internal weight must stay on the warm-up normalization so that
    // max_w remains a valid envelope. InjectRay also switches the GeometryBeam
    // to Generate mode and restricts selection to this (nu, target).
    m_evgen->InjectRay(nu_p4, achilles::PID{nu_in->pid()}, achilles::PID{tgt_attr->value()});

    // One Achilles call.  DetectorSim (not the adapter) decides whether to loop
    // on false according to SamplingMode:
    //   * TotalXSecRetry  : layer 1 used sigma_tot, so DetectorSim retries until
    //                       emit and the event carries max(max_w, raw_w).
    //   * EnvelopeNoRetry : layer 1 used max_w, single-shot rejection is correct.
    //
    // Either way, every call is one trial of the sigma estimator
    //     sigma = sum(w_trial) / N_trials,
    // so rejected trials must enter m_results with weight zero; otherwise the
    // E.C.4 mean is inflated by 1/efficiency.
    if(!m_evgen->GenerateSingleEvent()) {
        m_results += 0.0;
        NuGeom::Log().trace(
            "[adapter] trial rejected (E = {}, nu = {}, tgt = {}); running sigma = {} "
            "+/- {} pb over {} trials",
            nu_in->momentum().e(), nu_in->pid(), tgt_attr->value(), m_results.Mean(),
            m_results.Error(), m_results.Calls());
        return false;
    }

    // Per-trial weight in sigma units (pb).  The event itself carries the
    // partial-unweight weight max(1, raw_w/max_w_p) ~ 1 (the E.C.1 "CV"
    // value) in both modes.  The sigma scale is the per-(nu, target) vertex
    // envelope M = sum_p max_w_p in both modes too: a trial selects process p
    // with probability max_w_p / M, so the importance-sampling factor that
    // makes sum(w_trial)/N_trials estimate sum_p sigma_p is M, independent of
    // which process was drawn.
    constexpr double cm2_to_pb = 1.0 / achilles::Constant::TO_CM2;
    const double w_trial =
        m_evgen->LastEvent().Weight() *
        m_evgen->VertexEnvelope(achilles::PID{nu_in->pid()}, achilles::PID{tgt_attr->value()}) *
        cm2_to_pb;
    m_results += w_trial;
    NuGeom::Log().debug(
        "[adapter] trial emitted with w = {} pb; running sigma = {} +/- {} pb over {} "
        "trials",
        w_trial, m_results.Mean(), m_results.Error(), m_results.Calls());

    // Merge the generated (lab-frame) Achilles event onto NuGeom's GenEvent:
    // appends the nucleon-separation, primary outgoing, and cascade onto the
    // existing primary vertex (reusing the incoming nu), and stamps E.R.3
    // ProcessID, E.R.5 lab position, E.C.4 cross section, and E.C.1 weight.
    achilles::FillHepMC3Event(m_evgen->LastEvent(), evt, m_results);
    return true;
}

void AchillesAdapter::FinalizeRun() {
    if(m_evgen) m_evgen->Finalize();
}

} // namespace NuGeom
