#pragma once

#include "geom/GeneratorInterface.hh"

#include "Achilles/Statistics.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include <memory>
#include <optional>
#include <string>
#include <utility>

namespace achilles {
class EventGen;
} // namespace achilles

namespace NuGeom {

// NuGeom::GeneratorInterface implementation driven by Achilles.
//
// Achilles' own EventGen (constructor + Initialize) owns the full setup:
// Settings, NuclearModel, ProcessGroups, MultiChannel, VEGAS, Unweighter,
// splines.  This adapter only:
//
//   1. Requires the YAML to select a new Achilles beam type
//      ``Type: Geometry`` (achilles::GeometryBeam) so warm-up is FlatFlux-
//      like and event-time sampling is delta-function on the caller's nu
//      4-momentum.  VEGAS/MultiChannel are warm-up-transparent — the
//      integrator sees a FlatFlux during Optimize().
//
//   2. Switches that beam to Generate mode after Initialize() and injects
//      the per-vertex neutrino 4-momentum from the HepMC event NuGeometry
//      has pre-filled.
//
//   3. Asks Achilles for *one* phase-space sample per NuGeom-accepted
//      vertex (no retry).  The resulting achilles::Event carries the
//      partial-unweighter's max(1, W/max_w) in Event::m_wgt; the adapter
//      stamps it on the HepMC GenEvent as E.C.1.  Achilles' existing
//      E.C.4 TotalCrossSection stamp rides along via the translation.
//
//   4. Layer-1 vertex envelope returned from EnvelopeXSec is the Achilles
//      unweighter's V_ps * max_w for the relevant process group, not
//      sigma_tot.  This makes the two-layer composition exact under
//      partial unweighting with no retry and no adapter-side correction.
//
// See AchillesAdapter.cc for the (small) list of Achilles-side additions
// required to wire this up.
class AchillesAdapter : public GeneratorInterface {
  public:
    explicit AchillesAdapter(const std::string &yaml_config);
    // Construct from an already-loaded (and possibly programmatically
    // updated) YAML node -- e.g. the "Achilles" section of a geometry run
    // card.  `run_card_name` is recorded for provenance (Achilles.RunCard).
    explicit AchillesAdapter(const YAML::Node &config,
                             std::string run_card_name = "<geometry run card>");
    ~AchillesAdapter() override;

    void InitRun(std::shared_ptr<HepMC3::GenRunInfo> run_info) override;
    double TotalXSec(double energy, int nu_pdg, int target_pdg) override;
    double EnvelopeXSec(double energy, int nu_pdg, int target_pdg) override;
    bool GenerateEvent(HepMC3::GenEvent &evt) override;
    void FinalizeRun() override;

    // Keep Achilles' per-event weight convention in sync with the DetectorSim
    // SamplingMode. Call before InitRun(); default is EnvelopeNoRetry.
    void SetTotalXSecRetry(bool on) { m_total_xsec_retry = on; }

    // --- Beam energy support (MeV) -------------------------------------
    // The Achilles GeometryBeam treats injected rays outside its configured
    // Beam/EnergyRange as zero cross section, so the range must cover the
    // flux energies for every ray to be generatable.  These helpers let the
    // caller align the configured range with the flux support; they edit the
    // YAML before EventGen is constructed, so they must be called before
    // InitRun().

    // The configured Beam/EnergyRange, or nullopt if the run card omits it.
    std::optional<std::pair<double, double>> BeamEnergyRange();

    // Overwrite Beam/EnergyRange with [emin, emax].
    void SetBeamEnergyRange(double emin, double emax);

    // Widen (or set, if absent) Beam/EnergyRange to cover [emin, emax].
    // Returns true if the configuration was modified.
    bool EnsureBeamEnergyCoverage(double emin, double emax);

  private:
    // The Achilles config as a mutable YAML node; loads m_config_path on
    // first use when constructed from a file path.
    YAML::Node &Config();
    std::string m_config_path; // also used as the recorded run-card name
    std::optional<YAML::Node> m_config_node;

    // Owns the whole Achilles run.  EventGen's ctor + Initialize() do all
    // setup; this adapter never peeks inside at ProcessGroup, NuclearModel,
    // Sherpa, etc.
    std::unique_ptr<achilles::EventGen> m_evgen;

    // Running cross-section accumulator stamped onto each event as E.C.4.
    achilles::StatsData m_results;

    // Whether DetectorSim is in TotalXSecRetry (true) or EnvelopeNoRetry (false).
    bool m_total_xsec_retry = false;
};

} // namespace NuGeom
