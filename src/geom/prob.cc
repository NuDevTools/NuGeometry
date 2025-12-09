#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "geom/Material.hh"
#include "geom/Vector3D.hh"
#include <numeric>
#include <map>
#include <string>
#include "geom/Volume.hh"
#include "geom/Ray.hh"
#include "geom/LineSegment.hh"
#include "geom/Parser.hh"
#include "geom/Random.hh"
#include "spdlog/spdlog.h"

using NuGeom::LogicalVolume;
using NuGeom::PhysicalVolume;
using NuGeom::Vector3D;

using LineSegments = std::vector<NuGeom::LineSegment>;
using HandledRay = std::pair<std::vector<double>, LineSegments>;
using GeneratorCallback = std::function<double(double, size_t)>;
using RayGenCallback = std::function<std::pair<double, NuGeom::Ray>()>;

class DetectorSim {
public:
    void Setup(const std::string &geometry) {
        // Create world from geometry
        pugi::xml_document doc;
        pugi::xml_parse_result result = doc.load_file(geometry.c_str());
        if(!result)
            throw std::runtime_error("GDMLParser: Invalid file");
        NuGeom::GDMLParser parser(doc);
        world = parser.GetWorld();
    }

    void Setup(NuGeom::World world_) { world = world_; }

    void Init(const std::vector<std::pair<double, NuGeom::Ray>> &energy_rays, double safety_factor) {
        // Calculate the max_prob and store it with the safety factor
        for(const auto &[energy, ray] : energy_rays) {
            auto [probs, segments] = HandleRay(energy, ray);
            auto prob_tot = std::accumulate(probs.begin(), probs.end(), 0.0);
            if(prob_tot > max_prob) max_prob = prob_tot;
        }

        max_prob *= safety_factor;
    }

    std::vector<NuGeom::Material> GetMaterials() const { return world.GetMaterials(); }

    void GetMeanFreePath(const std::vector<double> &cross_section) {
        // Fill result mfp
        if(cross_section.size() != m_mats.size())
            throw "ERROR";

        for(size_t i = 0; i < m_mats.size(); ++i) {
            m_mfp[m_mats[i]] = cross_section[i];
        }
    }

    std::pair<Vector3D, NuGeom::Material> GetInteraction(double energy, const NuGeom::Ray &ray) const {
        auto [probs, segments] = HandleRay(energy, ray);
        auto tot_probs = std::accumulate(probs.begin(), probs.end(), 0.0);
        if(tot_probs / max_prob < NuGeom::Random::Instance().Uniform(0.0, 1.0)) {
            // Return a dummy material if no interaction occurs
            return {{}, {}};
        }

        // Determine hit location and hit material
        double rand = NuGeom::Random::Instance().Uniform(0.0, tot_probs);
        Vector3D position;
        NuGeom::Material mat;
        double sum = 0;
        for(size_t i = 0; i < probs.size(); ++i) {
            sum += probs[i];
            if(sum > rand) {
                position = segments[i].Start() + (segments[i].End() - segments[i].Start()) * (sum - rand)/probs[i];
                mat = segments[i].GetMaterial();
                break;
            }
        }
        return {position, mat};
    }

    // Expects a function that returns the total cross section per nucleus given the energy and the PDG code of the target
    void SetGeneratorCallback(GeneratorCallback &xsec) { xsec_callback = xsec; }
    // Expects a function that returns the next ray to propagate 
    void SetRayGenCallback(RayGenCallback &ray_gen) { ray_gen_callback = ray_gen; }

    //Event Generator is in control 
    //Materials GetMaterials(Ray)

    std::set<NuGeom::Material> GetMaterials(const NuGeom::Ray &rays) {
        std::set<NuGeom::Material> mats;
        auto segments = world.GetLineSegments(rays);
        for(const auto &segment : segments) {
            mats.insert(segment.GetMaterial());
        }
        return mats;
    }

    //gets linesegments given a ray
    LineSegments GetLineSegments(const NuGeom::Ray &ray) {
        return world.GetLineSegments(ray);
    }

    std::vector<double> EvaluateProbs(const LineSegments &segments, const std::map<NuGeom::Material, double> &xsecsmaps) {
        std::vector<double> probs;

        for(const auto &segment : segments) {
            auto material = segment.GetMaterial();
            auto cross_section = xsecsmaps.at(material);
            auto meanfreepath = 1.0 / (cross_section * material.Density());
            probs.push_back(segment.Length() / meanfreepath);

        }
        auto total_prob = std::accumulate(probs.begin(), probs.end(), 0.0);
        if (total_prob > max_prob) {
            max_prob = total_prob;
        }

        return probs;
    }

    NuGeom::Vector3D Interaction(const LineSegments &segments, const std::map<NuGeom::Material, double> &xsecsmaps) {

        // Evaluate interaction probability
        auto prob = EvaluateProbs(segments, xsecsmaps);

        // Draw random uniform [0,1)
        double r1 = NuGeom::Random::Instance().Uniform(0.0, 1.0);

        // If interaction occurs
        auto total_prob = std::accumulate(prob.begin(), prob.end(), 0.0);
        if (total_prob/(max_prob * 1.5) < r1) {
            // No interaction occurred
            return NuGeom::Vector3D(9e9, 9e9, 9e9);
        }
        double r2 = NuGeom::Random::Instance().Uniform(0.0, 1.0);
        // Select interaction segment
        double cumulative_prob = 0.0;
        size_t idx = 0;
        for (size_t i = 0; i < prob.size(); ++i) {
            cumulative_prob += prob[i];
            if (cumulative_prob >= r2) {
                idx = i;
                break;
            }
        }
        // Determine interaction location within the segment
        double segment_length = segments[idx].Length();
        double interaction_distance = (r2 - (cumulative_prob - prob[idx])) / prob[idx] * segment_length;
        NuGeom::Vector3D interaction_point = segments[idx].Start() + (segments[idx].End() - segments[idx].Start()).Unit() * interaction_distance;
        return interaction_point;
    }

private:
    HandledRay HandleRay(double energy, const NuGeom::Ray &ray) const;
    double CalculateMeanFreePath(double energy, const NuGeom::Material &material) const;
    double CalculateCrossSection(double energy, const NuGeom::Material mat) const {
        double cross_section = 0;
        for(const auto &elm : mat.Elements()) {
            cross_section += xsec_callback(energy, elm.PDG()); 
        }
        return cross_section;
    }

    NuGeom::World world;
    std::vector<std::shared_ptr<NuGeom::Shape>> shapes;
    std::vector<NuGeom::Material> m_mats;
    std::map<NuGeom::Material, double> m_mfp;
    GeneratorCallback xsec_callback;
    RayGenCallback ray_gen_callback;
    double max_prob = 0;
};

HandledRay DetectorSim::HandleRay(double energy, const NuGeom::Ray &ray) const {
    LineSegments segments = world.GetLineSegments(ray);

    // Calculate probability to interact for each line segment
    std::vector<double> probs;
    for(const auto &segment : segments) {
        auto material = segment.GetMaterial();
        auto meanfreepath = CalculateMeanFreePath(energy, material);
        probs.push_back(segment.Length()/meanfreepath);
    }

    return {probs, segments};
}

double DetectorSim::CalculateMeanFreePath(double energy, const NuGeom::Material &mat) const {
    auto xsec = CalculateCrossSection(energy, mat);
    return 1.0 / (xsec * mat.Density());
}

NuGeom::Ray ShootRay(const NuGeom::Vector3D &corner1, const NuGeom::Vector3D &corner2) {
    auto [xmin, xmax] = std::minmax(corner1.X(), corner2.X());
    auto [ymin, ymax] = std::minmax(corner1.Y(), corner2.Y());
    auto [zmin, zmax] = std::minmax(corner1.Z(), corner2.Z());

    auto rand = NuGeom::Random::Instance();

    NuGeom::Vector3D position{rand.Uniform(xmin, xmax), rand.Uniform(ymin, ymax), rand.Uniform(zmin, zmax)};
    double costheta = rand.Uniform(0.0, 1.0);
    double sintheta = std::sqrt(1-costheta*costheta);
    double phi = rand.Uniform(0.0, 2*M_PI);
    NuGeom::Vector3D direction{sintheta*cos(phi), sintheta*sin(phi), costheta};

    return NuGeom::Ray(position, direction);
}

int main(){
    // Define materials in the detector
    NuGeom::Material mat("Water", 1.0, 2);
    NuGeom::Material mat1("Argon", 9.0, 1);
    mat.AddElement(NuGeom::Element("Hydrogen", 1, 1), 2);
    mat.AddElement(NuGeom::Element("Oxygen", 8, 16), 1);
    mat1.AddElement(NuGeom::Element("Argon", 18, 40), 1);

    // Define the interaction geometry
    // Define the inner detector
    auto inner_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{1, 1, 1}); // Define a 1x1x1 box

    auto inner_vol = std::make_shared<LogicalVolume>(mat, inner_box); 
    NuGeom::RotationX3D rot(45*M_PI/180.0);
    auto inner_pvol = std::make_shared<PhysicalVolume>(inner_vol, NuGeom::Transform3D{}, rot);


    // Define the outer detector
    auto outer_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2});
    auto outer_vol = std::make_shared<LogicalVolume>(mat1, outer_box); 
    outer_vol->AddDaughter(inner_pvol);
    inner_vol->SetMother(outer_vol);
    NuGeom::RotationX3D rot2(30*M_PI/180.0);
    auto outer_pvol = std::make_shared<PhysicalVolume>(outer_vol, NuGeom::Transform3D{}, rot2);
    inner_pvol->SetMother(outer_pvol);

    // Define the "World"
    auto world_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 4});
    auto world = std::make_shared<LogicalVolume>(mat, world_box);
    outer_vol->SetMother(world); 
    world -> AddDaughter(outer_pvol);

    // Calculate interaction location
    std::map<std::string,double> meanfreepaths{
        {"Water",1e+7},
        {"Air",1e+7},
        {"Argon",1e+7}
    };

    // Shoot Rays and get LineSegments
    double prob_max = 0;
    size_t ntrials = 1 << 22;
    for(size_t i = 0; i < ntrials; ++i) {
        NuGeom::Ray ray = ShootRay({-2, -2, -2}, {2, 2, -2});
        std::vector<NuGeom::LineSegment> segments0;
        world->GetLineSegments(ray, segments0);

        std::vector<double> probs0(segments0.size());
        std::vector<double> probs1(segments0.size());
        std::vector<double> seglength0(segments0.size());

        std::vector<std::string> material0(segments0.size());

        // Calculate probability to interact for each line segment
        for (size_t j=0; j<segments0.size(); ++j){
            material0[j]=segments0[j].GetMaterial().Name();

            // NOTE: This only works for l/meanfreepath tiny
            probs0[j]=segments0[j].Length()/meanfreepaths[material0[j]];
            probs1[j]=1-exp(-segments0[j].Length()/meanfreepaths[material0[j]]);

            // For testing only:
            seglength0[j] = segments0[j].Length();
        }
        double normconst0=std::accumulate(probs0.begin(), probs0.end(), 0.0);
        if (normconst0>prob_max) {
            std::cout << "Trial: " << i << " Ray: O(" << ray.Origin() << ") D(" << ray.Direction() << ")\n";
            std::cout << "Hit " << segments0.size() << " segments\n";
            double normconst1=std::accumulate(probs1.begin(), probs1.end(), 0.0);
            std::cout << "Exact = " << normconst1 << " Approx = " << normconst0 << "\n";
            std::cout << "Old prob max = " << prob_max << "\n";
            std::cout << "New prob max = " << normconst0 << "\n";
            prob_max=normconst0;
        }
    }

    constexpr double safety_factor = 1.5;
    prob_max *= safety_factor;

    auto rand = NuGeom::Random::Instance();
    double nrays = 0;
    size_t nhits = 0;
    std::ofstream hist;
    hist.open("hit_locations.txt");

    // Send Ray for check to hit
    for(size_t i = 0; i < ntrials; ++i) {
        nrays += 1.0/prob_max;
        if(i % 10000 == 0) {
            std::cout << "Shot " << nrays << " rays\n";
        }
        NuGeom::Ray ray = ShootRay({-2, -2, -2}, {2, 2, -2});
        std::vector<NuGeom::LineSegment> segments0;

        world->GetLineSegments(ray, segments0);

        std::vector<double> probs0(segments0.size());
        std::vector<double> probs1(segments0.size());
        std::vector<double> seglength0(segments0.size());

        std::vector<std::string> material0(segments0.size());

        // Calculate probability to interact for each line segment
        for (size_t j=0; j<segments0.size(); ++j){
            material0[j]=segments0[j].GetMaterial().Name();
            // NOTE: This only works for l/meanfreepath tiny
            probs0[j]=segments0[j].Length()/meanfreepaths[material0[j]];
            // std::cout << -segments0[i].Length()/meanfreepaths[material0[i]] << std::endl;
            probs1[j]=1-exp(-segments0[j].Length()/meanfreepaths[material0[j]]);

            // For testing only:
            seglength0[j] = segments0[j].Length();
        }
        double normconst0=std::accumulate(probs0.begin(), probs0.end(), 0.0);
        if(normconst0 / prob_max < rand.Uniform(0.0, 1.0)) {
            continue;
        }
        nhits++;
        std::cout << nhits << " / " << nrays << "\n";

        for(size_t j=0; j<segments0.size(); ++j) probs0[j]/=normconst0;
        double sum = 0;
        size_t idx = 0;
        double r2 = rand.Uniform(0.0, 1.0);
        NuGeom::Vector3D position;
        for(size_t j = 0; j < probs0.size(); ++j) {
            sum += probs0[j]; 
            if(sum > r2) {
                idx = j;
                position = segments0[j].Start() + (segments0[j].End() - segments0[j].Start()) * (sum - r2)/probs0[j]; 
                break;
            }
        }

        std::cout << "Start: " << segments0[idx].Start() << "\n";
        std::cout << "End: " << segments0[idx].End() << "\n";
        std::cout << "Sum: " << sum << " R2: " << r2 << " Prob: " << probs0[idx] << " Dist: " << (sum - r2)/probs0[idx] << "\n";
        std::cout << "Hit Location: " << position << "\n";

        hist << position.X() << "," << position.Y() << "," << position.Z() << "\n";
        if(nhits == 10000) break;
    } 

    hist.close();

    return 0;
}
