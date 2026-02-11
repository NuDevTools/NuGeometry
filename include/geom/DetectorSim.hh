#pragma once

#include "geom/LineSegment.hh"
#include "geom/Ray.hh"
#include "geom/Volume.hh"
#include "geom/World.hh"

#include <fstream>

namespace NuGeom {

using LineSegments = std::vector<NuGeom::LineSegment>;
using EnergyRay = std::pair<double, NuGeom::Ray>;
using HandledRay = std::pair<std::vector<double>, LineSegments>;
using GeneratorCallback = std::function<double(double, size_t)>;
using RayGenCallback = std::function<EnergyRay()>;

class DetectorSim {
  public:
    DetectorSim(double safety_factor = 1.5) : m_safety_factor{safety_factor} {}

    void SetEventFile(const std::string &outfile) { m_outfile = std::ofstream(outfile); }

    void Setup(const std::string &geometry);
    void Setup(NuGeom::World world_) { world = world_; }
    void Init(size_t nrays);
    std::vector<NuGeom::Material> GetMaterials() const { return world.GetMaterials(); }
    void GenerateEvents(size_t nevents) const;
    void GenerateEvents(double POT) const;

    // Expects a function that returns the total cross section per nucleus given the energy and the
    // PDG code of the target
    void SetGeneratorCallback(GeneratorCallback xsec) { xsec_callback = xsec; }
    // Expects a function that returns the next ray to propagate
    void SetRayGenCallback(RayGenCallback ray_gen) { ray_gen_callback = ray_gen; }
    void SetMaxProb(double prob) { max_prob = prob * m_safety_factor; }
    double GetMaxProb() const { return max_prob; }

    std::set<NuGeom::Material> GetMaterials(const LineSegments &segments) const;
    // Gets linesegments given a ray
    LineSegments GetLineSegments(const NuGeom::Ray &ray) const {
        return world.GetLineSegments(ray);
    }

    std::vector<double> EvaluateProbs(const LineSegments &segments,
                                      const std::map<NuGeom::Element, double> &xsecsmaps);

    std::pair<NuGeom::Vector3D, Material>
    Interaction(const LineSegments &segments, const std::map<NuGeom::Element, double> &xsecsmaps);

  private:
    HandledRay HandleRay(double energy, const NuGeom::Ray &ray) const;
    double CalculateMeanFreePath(double energy, const NuGeom::Material &material) const;
    std::tuple<Vector3D, Material, EnergyRay> GetInteraction() const;

    NuGeom::World world;
    std::vector<std::shared_ptr<NuGeom::Shape>> shapes;
    std::vector<NuGeom::Material> m_mats;
    std::map<NuGeom::Material, double> m_mfp;
    GeneratorCallback xsec_callback;
    RayGenCallback ray_gen_callback;
    double max_prob = 0;
    double m_safety_factor;

    mutable double m_pot = 0; // This is a code smell. Figure out how to handle better
    mutable std::ofstream m_outfile;
};

} // namespace NuGeom
