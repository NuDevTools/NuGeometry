#pragma once

#include "geom/FluxSource.hh"
#include "geom/LineSegment.hh"
#include "geom/Ray.hh"
#include "geom/Volume.hh"
#include "geom/World.hh"

#include <memory>

namespace HepMC3 {
class GenRunInfo;
class WriterAscii;
} // namespace HepMC3

namespace NuGeom {

class GeneratorInterface;

using LineSegments = std::vector<NuGeom::LineSegment>;
using EnergyRay = std::pair<double, NuGeom::Ray>;
using HandledRay = std::pair<std::vector<double>, LineSegments>;

// Legacy callback kept for callers driving the simulation themselves
// (generator-led mode via EvaluateProbs/Interaction).
using GeneratorCallback = std::function<double(double, size_t)>;
using RayGenCallback = std::function<EnergyRay()>;

class DetectorSim {
  public:
    // Which envelope layer-1 uses, and whether layer-2 retries.  The two
    // modes must produce the same physical observables on the same input
    // — use both to cross-check a generator adapter.  See
    // GeneratorInterface::TotalXSec / EnvelopeXSec for the math.
    enum class SamplingMode {
        TotalXSecRetry,  // layer 1 = sigma_tot; retry Achilles until emit
        EnvelopeNoRetry, // layer 1 = max_w-style envelope; single shot
    };

    explicit DetectorSim(double safety_factor = 1.5);
    ~DetectorSim();

    void SetSamplingMode(SamplingMode m) { m_mode = m; }
    SamplingMode GetSamplingMode() const { return m_mode; }

    void SetEventFile(const std::string &outfile);

    void Setup(const std::string &geometry);
    void Setup(NuGeom::World world_) { world = world_; }

    // --- Layer-1 envelope calibration ----------------------------------------
    // Layer 1 accepts a ray with probability accept_w / max_prob, where
    //     accept_w = flux_weight * window_area * P_int(ray)
    // so max_prob must bound the product.  It factorizes cleanly:
    //     max_prob = [safety * max_ray P_int] * [max_ray flux_weight*window_area]
    //                 \___ geometry x xsec ___/  \______ flux scale ________/
    // Init() pins the first factor with rays *it generates itself* -- chords
    // aimed at the placed volumes, densest first, scanned over an energy grid --
    // so it needs neither the flux file nor a representative sample of it.  The
    // second factor is a running maximum over the rays already thrown.
    //
    // That product is a genuine bound, but a loose one: the highest-weight ray
    // and the longest-path ray are different rays, so max(w)*max(P) overshoots
    // max(w*P) by a large factor (~45x on the DUNE ND hall, i.e. ~45x more
    // throws per event).  So it is used as a *ceiling* only, and the envelope
    // actually thrown against is the tighter
    //     min( ceiling, safety * max over previous rays of accept_w )
    //
    // Moving the envelope mid-run does not bias events/POT: both the accept
    // probability and the POT charge divide by the same max_prob, so their
    // ratio is envelope-independent for every throw.  What that argument needs
    // is that max_prob for a throw depends only on *previous* rays -- letting it
    // track the current ray's flux weight would cancel the weight out of the
    // acceptance and move it into the POT charge.  So a ray that beats the
    // running maximum is thrown against the old envelope and, if it overflows
    // it, counted as clipped (RunStats::clipped): clipping caps the acceptance
    // at 1 and is the one residual bias, always a shortfall, and it is measured
    // rather than assumed negligible.
    //
    // Energy range for the scan; required by Init().  In GeV, matching
    // FluxSample::energy.  Use the generator's supported beam range.
    void SetEnergyRange(double emin, double emax);
    // Neutrino species to scan (default: nu_mu, nu_mu_bar, nu_e, nu_e_bar).
    // Species the generator has no process for contribute zero.
    void SetFluxSpecies(std::vector<int> pdgs) { m_species = std::move(pdgs); }
    // Rays consumed after Init() purely to settle the envelope's running
    // maximum, before any exposure is charged or any event produced.  They are
    // traced (that is how accept_w is known) but no generator call is made, so
    // they cost a fraction of a normal throw, and keeping them out of the
    // accounting removes the transient where the envelope is still far from its
    // steady-state value.  This is what the old flux warm-up did, except inline
    // and without a separate pass over the file.
    void SetEnvelopeBurnIn(size_t nrays) { m_envelope_burn_in = nrays; }

    // Calibrate the envelope from `nrays` self-generated chords.  Requires the
    // generator and an energy range; does NOT touch the flux callback.
    void Init(size_t nrays);
    // Legacy calibration: stream `nrays` real flux samples and take the max of
    // their acceptance weights.  Needs the flux callback, consumes those rays
    // from it, and only ever sees the part of the tail it happened to draw.
    void InitFromFlux(size_t nrays);

    std::vector<NuGeom::Material> GetMaterials() const { return world.GetMaterials(); }
    void GenerateEvents(size_t nevents);
    void GenerateEvents(double POT);

    // --- Geometry-driven mode -------------------------------------------------
    // Preferred path: NuGeometry owns the loop; the generator supplies xsec
    // and fills in final states via the GeneratorInterface.
    void SetGenerator(std::shared_ptr<GeneratorInterface> gen);
    void SetFluxCallback(FluxCallback flux) { flux_callback = std::move(flux); }

    // --- Generator-driven mode ------------------------------------------------
    // Legacy helpers kept for callers (e.g. external generators) that own the
    // event loop themselves and only need geometry/xsec weighting.
    void SetGeneratorCallback(GeneratorCallback xsec) { xsec_callback = std::move(xsec); }
    void SetRayGenCallback(RayGenCallback ray_gen) { ray_gen_callback = std::move(ray_gen); }

    // Pin the envelope by hand.  Unlike the calibrated paths this one is never
    // adapted at runtime: a caller who states the envelope owns it, and a ray
    // that overflows it is reported as clipped rather than silently rescaled.
    void SetMaxProb(double prob) {
        max_prob = prob * m_safety_factor;
        m_geom_envelope = 0.0;
        m_adaptive_envelope = false;
    }
    // The envelope in force.  After Init() this is 0 until the first ray is
    // thrown: Init() fixes only the geometry x cross-section factor, and the
    // flux-weight factor is not known until a flux sample has been seen.
    double GetMaxProb() const { return max_prob; }
    double GetAccumulatedPOT() const { return m_pot; }

    // End-of-run bookkeeping.  `sum_weight` accumulates the emitted events'
    // E.C.1 weights; comparing it with `emitted` is what exposes a generator
    // whose partial unweighter leaves an overweight tail -- see
    // ReportRunSummary() and the SamplingMode notes above.
    struct RunStats {
        double max_prob = 0;         // envelope in force at the end of the run
        size_t thrown = 0;           // rays drawn from the flux
        size_t accepted = 0;         // rays that passed layer 1
        size_t emitted = 0;          // events written
        double sum_weight = 0;       // sum of the emitted events' E.C.1 weights
        size_t envelope_growths = 0; // times the envelope had to be raised
        size_t clipped = 0;          // rays whose accept_w exceeded the envelope
        double clipped_excess = 0;   // sum of (accept_w/envelope - 1) over those
        double pot = 0;
    };
    const RunStats &GetRunStats() const { return m_stats; }
    // Log events/POT, sum(w)/POT and their ratio, warning when they disagree.
    void ReportRunSummary() const;

    std::set<NuGeom::Material> GetMaterials(const LineSegments &segments) const;
    LineSegments GetLineSegments(const NuGeom::Ray &ray) const {
        return world.GetLineSegments(ray);
    }

    std::vector<double> EvaluateProbs(const LineSegments &segments,
                                      const std::map<NuGeom::Element, double> &xsecsmaps);

    std::pair<NuGeom::Vector3D, Material>
    Interaction(const LineSegments &segments, const std::map<NuGeom::Element, double> &xsecsmaps);

  private:
    // Shared core of vertex picking.  Given per-segment probabilities and a
    // uniform draw r in [0, total_prob), returns the segment index and the
    // vertex position within it.  Used by both Interaction() and the
    // NuHepMC event-generation path.
    struct VertexPick {
        Vector3D position;
        Material material;
        size_t segment_idx;
    };
    VertexPick PickVertex(const LineSegments &segments, const std::vector<double> &probs,
                          double r) const;

    HandledRay HandleRay(const FluxSample &fs) const;
    HandledRay HandleRay(double energy, int nu_pdg, const NuGeom::Ray &ray) const;
    // Per-element column density (atoms/cm^2) along already-traced segments.
    // Init()'s scan collapses onto this: the interaction probability at any
    // (energy, species) is then sum_e N_e * sigma_e, so a chord is traced once
    // and its segments walked once for the whole grid instead of once per grid
    // point.
    std::map<NuGeom::Element, double> ColumnDensity(const LineSegments &segments) const;
    // The chords Init() scans: the world box's face-to-face and corner-to-corner
    // extremes, then random surface-point pairs.  Long chords through dense
    // material are what maximizes the column density, and a real flux ray can
    // never beat the longest chord of the world box.
    std::vector<NuGeom::Ray> GenerateProbeRays(size_t nrays) const;
    double CalculateMeanFreePath(double energy, int nu_pdg, const NuGeom::Material &material) const;
    // Cross section for one element, using whichever quantity the sampling mode is
    // built on: sigma_tot (TotalXSecRetry) or the max-weight envelope
    // (EnvelopeNoRetry). The same quantity drives the mean free path and the
    // per-element selection, so the two stay consistent.
    double ElementXSec(double energy, int nu_pdg, const NuGeom::Element &elm) const;

    // Draws one flux sample, charges POT, returns either a populated
    // HepMC3::GenEvent (interaction occurred) or std::nullopt.  Generator is
    // consulted to fill the final state.
    bool ProduceEvent();

    // Element sampling within a chosen material, weighted by n_e * sigma_e (the
    // mode-appropriate cross section) so the per-element rate is unbiased.
    Element PickElement(const Material &mat, double energy, int nu_pdg) const;

    NuGeom::World world;
    std::vector<std::shared_ptr<NuGeom::Shape>> shapes;
    std::vector<NuGeom::Material> m_mats;
    std::map<NuGeom::Material, double> m_mfp;

    std::shared_ptr<GeneratorInterface> m_generator;
    FluxCallback flux_callback;

    // Legacy generator-led mode callbacks.
    GeneratorCallback xsec_callback;
    RayGenCallback ray_gen_callback;

    double max_prob = 0;
    double m_safety_factor;
    SamplingMode m_mode = SamplingMode::EnvelopeNoRetry;

    // Factorized envelope (see SetEnergyRange / Init).  m_geom_envelope == 0
    // means "fixed envelope" -- set by SetMaxProb or InitFromFlux -- in which
    // case max_prob is used as-is and only the overflow guard can raise it.
    double m_geom_envelope = 0;
    double m_flux_scale = 0;
    double m_joint_max = 0; // max accept_w over the rays already thrown
    bool m_adaptive_envelope = false;
    size_t m_envelope_burn_in = 1000;
    double m_emin = 0, m_emax = 0;
    std::vector<int> m_species{14, -14, 12, -12};

    RunStats m_stats;

    double m_pot = 0;
    size_t m_event_number = 0;

    std::shared_ptr<HepMC3::GenRunInfo> m_run_info;
    std::unique_ptr<HepMC3::WriterAscii> m_writer;
};

} // namespace NuGeom
