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

    // Sample rays to fix max_prob with the safety factor applied.  Uses the
    // flux callback and the generator's TotalXSec, so both must be set first.
    void Init(size_t nrays);

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

    void SetMaxProb(double prob) { max_prob = prob * m_safety_factor; }
    double GetMaxProb() const { return max_prob; }
    double GetAccumulatedPOT() const { return m_pot; }

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

    double m_pot = 0;
    size_t m_event_number = 0;

    std::shared_ptr<HepMC3::GenRunInfo> m_run_info;
    std::unique_ptr<HepMC3::WriterAscii> m_writer;
};

} // namespace NuGeom
