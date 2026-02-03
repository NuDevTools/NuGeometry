#include <cstdlib>
#include <cmath>
#include <ctime>
#include "geom/Material.hh"
#include "geom/Vector3D.hh"
#include <map>
#include <string>
#include "geom/Volume.hh"
#include "geom/Ray.hh"
#include "geom/Random.hh"
#include "geom/World.hh"
#include "geom/DetectorSim.hh"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"
#include "CLI/CLI.hpp"


class TestRayGen {
public:
    TestRayGen(double emin, double emax,
               NuGeom::Vector3D corner1, NuGeom::Vector3D corner2) :
        m_emin{emin}, m_emax{emax} {
        auto [xmin, xmax] = std::minmax(corner1.X(), corner2.X());
        auto [ymin, ymax] = std::minmax(corner1.Y(), corner2.Y());
        auto [zmin, zmax] = std::minmax(corner1.Z(), corner2.Z());
        m_xmin = xmin;
        m_xmax = xmax;
        m_ymin = ymin;
        m_ymax = ymax;
        m_zmin = zmin;
        m_zmax = zmax;
    }

    NuGeom::EnergyRay GetRay() const {
        auto ray = ShootRay();
        double energy = NuGeom::Random::Instance().Uniform(m_emin, m_emax);
        return {energy, ray};
    }

private:
    NuGeom::Ray ShootRay() const {
        auto rand = NuGeom::Random::Instance();

        NuGeom::Vector3D position{rand.Uniform(m_xmin, m_xmax),
                                  rand.Uniform(m_ymin, m_ymax),
                                  rand.Uniform(m_zmin, m_zmax)};
        double costheta = rand.Uniform(0.0, 1.0);
        double sintheta = std::sqrt(1-costheta*costheta);
        double phi = rand.Uniform(0.0, 2*M_PI);
        NuGeom::Vector3D direction{sintheta*cos(phi), sintheta*sin(phi), costheta};

        return NuGeom::Ray(position, direction, 1e2);
    }

    double m_emin, m_emax;
    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;
};

class TestEventGen {
public:
    TestEventGen(std::map<size_t, double> xsec) : m_xsec{xsec} {}

    double CrossSection(double, size_t pdg) const {
        return m_xsec.at(pdg);
    }

    std::map<NuGeom::Material, double> EvaluateCrossSections(double energy, const std::set<NuGeom::Material> &mats) const {
        std::map<NuGeom::Material,double> result;
        for (const auto &mat : mats) {
            for (size_t i = 0; i < mat.NElements(); ++i) {
                auto elem = mat.Elements()[i];
                size_t pdg = 1000000000 + elem.Z() * 10000 + elem.A() * 10;
                double xsec = CrossSection(energy, pdg);
                result[mat] += xsec;
            }
        }
        return result;
    }

    std::map<size_t, double> GetXSecs(double energy, const std::set<NuGeom::Element> &mats) {
        std::map<size_t, double> result;
        for (const auto &elm : mats) {
            if (result.find(elm.PDG()) == result.end()) {
                result[elm.PDG()] = CrossSection(energy, elm.PDG());
            }
    }
    return result;
    }

    

private:
    std::map<size_t, double> m_xsec;
};

int main(int argc, char **argv){
    auto console = spdlog::stdout_color_mt("NuGeom");
    spdlog::set_default_logger(console);
    spdlog::set_pattern("[%n] [%^%l%$] %v");
    CLI::App app("Neutrino Geometry Driver");
    argv = app.ensure_utf8(argv);
    int verbosity = 0;
    std::string outfile = "hit_locations1.txt";
    size_t nevents = 1000;
    double pot = 1e23;
    app.add_flag("-v,--verbose", verbosity, "Increase the verbosity level");
    app.add_option("-f,--file", outfile, "File to write out the locations of the hits.");
    app.add_option("-e,--events", nevents, "Number of events to generate. (Default 1000)");
    app.add_option("-p,--pot", pot, "Number of protons on target (POT). (Default 10^23)");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    if(verbosity == 1) spdlog::set_level(spdlog::level::debug);
    else if(verbosity == 2) spdlog::set_level(spdlog::level::trace);


    // Define materials in the detector
    NuGeom::Material mat("Water", 1.0, 2);
    NuGeom::Material mat1("Argon", 9.0, 1);
    mat.AddElement(NuGeom::Element("Hydrogen", 1, 1), 2);
    mat.AddElement(NuGeom::Element("Oxygen", 8, 16), 1);
    mat1.AddElement(NuGeom::Element("Argon", 18, 40), 1);

    // Define the interaction geometry
    // Define the inner detector
    auto inner_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{1, 1, 1}); // Define a 1x1x1 box

    auto inner_vol = std::make_shared<NuGeom::LogicalVolume>(mat, inner_box); 
    NuGeom::RotationX3D rot(45*M_PI/180.0);
    auto inner_pvol = std::make_shared<NuGeom::PhysicalVolume>(inner_vol, NuGeom::Transform3D{}, rot);

    // Define the outer detector
    auto outer_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2});
    auto outer_vol = std::make_shared<NuGeom::LogicalVolume>(mat1, outer_box); 
    outer_vol->AddDaughter(inner_pvol);
    inner_vol->SetMother(outer_vol);
    NuGeom::RotationX3D rot2(30*M_PI/180.0);
    auto outer_pvol = std::make_shared<NuGeom::PhysicalVolume>(outer_vol, NuGeom::Transform3D{}, rot2);
    inner_pvol->SetMother(outer_pvol);

    // Define the "World"
    auto world_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 4});
    auto world_vol = std::make_shared<NuGeom::LogicalVolume>(mat, world_box);
    outer_vol->SetMother(world_vol); 
    world_vol -> AddDaughter(outer_pvol);
    NuGeom::World world(world_vol);


    NuGeom::DetectorSim sim;
    sim.Setup(world);
    sim.SetEventFile(outfile);


    std::map<size_t, double> xsec_map = {{1000010010, 1e-38},
                                         {1000080160, 1e-38},
                                         {1000180400, 1e-38}};

    double max_prob = 0;
    auto raygen = std::make_shared<TestRayGen>(0, 10, NuGeom::Vector3D{-2, -2, -2},
                                               NuGeom::Vector3D{2, 2, -2});
    TestEventGen gen1(xsec_map);

    
    std::vector<NuGeom::EnergyRay>  rays;
    size_t ntest_rays = 1 << 20;
    for (size_t i = 0; i < ntest_rays; ++i) {
        rays.push_back(raygen->GetRay());
    }
   for (const auto &ray : rays) {
    
    std::vector<NuGeom::LineSegment> segments = sim.GetLineSegments(ray.second);
     std::set<NuGeom::Material> materials = sim.GetMaterials(segments);
     std::map<NuGeom::Material, double> xsecs= gen1.EvaluateCrossSections(ray.first, materials);
     std::vector<double> probs = sim.EvaluateProbs(segments, xsecs);
     double prob_sum = std::accumulate(probs.begin(), probs.end(), 0.0);
     if (prob_sum > max_prob) {
         max_prob = prob_sum;
     }

   }
   sim.SetMaxProb(max_prob);
   spdlog::info("Set max probability to {}", max_prob);
    std::ofstream hist;
    hist.open("hit_locations2.txt", std::ios::app);
    if (!hist.is_open()) {
    throw std::runtime_error("Failed to open hit_locations2.txt");
}
//for rays const 
size_t ntest_hits = 0;


   while (ntest_hits < nevents) {
        auto ray= raygen->GetRay();
    std::vector<NuGeom::LineSegment> segments0 = sim.GetLineSegments(ray.second);
    std::set<NuGeom::Material> materials = sim.GetMaterials(segments0);
    std::map<NuGeom::Material, double> material_xsec_map = gen1.EvaluateCrossSections(ray.first, materials);
    auto [interaction_point, interaction_material] = sim.Interaction(segments0, material_xsec_map);
     if (interaction_point.X() < 9e9) {
            ntest_hits++;
         // Valid interaction point
         hist << interaction_point.X() << " "
              << interaction_point.Y() << " "
              << interaction_point.Z() << " "
              << interaction_material.Name() << "\n";
    }
        else {
            // No interaction occurred
            spdlog::info("No interaction occurred for this ray.");
     }
     
   }
   hist.close();
        

   return 0;
}