#include "geom/DetectorSim.hh"
#include "geom/LineSegment.hh"
#include "geom/Parser.hh"
#include "geom/Random.hh"

#include "spdlog/spdlog.h"
#include <numeric>

using NuGeom::DetectorSim;

void DetectorSim::Setup(const std::string &geometry) {
    // Create world from geometry
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(geometry.c_str());
    if(!result) throw std::runtime_error("GDMLParser: Invalid file");
    NuGeom::GDMLParser parser(doc);
    world = parser.GetWorld();
}

void DetectorSim::Init(size_t nrays) {
    // Calculate the max_prob and store it with the safety factor
    for(size_t i = 0; i < nrays; ++i) {
        auto [energy, ray] = ray_gen_callback();
        auto [probs, segments] = HandleRay(energy, ray);
        auto prob_tot = std::accumulate(probs.begin(), probs.end(), 0.0);
        if(prob_tot > max_prob) {
            spdlog::debug("Updating max prob: {} -> {}", max_prob, prob_tot);
            max_prob = prob_tot;
        }
    }

    max_prob *= m_safety_factor;
    spdlog::info("Maximum probability found with safety factor ({}): {}", m_safety_factor,
                 max_prob);
}

std::pair<NuGeom::Vector3D, NuGeom::Material> DetectorSim::GetInteraction() const {
    auto [energy, ray] = ray_gen_callback();
    auto [probs, segments] = HandleRay(energy, ray);
    auto tot_probs = std::accumulate(probs.begin(), probs.end(), 0.0);
    m_pot += ray.POT() / max_prob;
    spdlog::debug("Probability / Max_Prop = {}", tot_probs / max_prob);
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
            position = segments[i].Start() +
                       (segments[i].End() - segments[i].Start()) * (sum - rand) / probs[i];
            mat = segments[i].GetMaterial();
            break;
        }
    }
    spdlog::debug("Interaction occured! Placing vertex at {} on {}", position, mat.Name());
    return {position, mat};
}

std::set<NuGeom::Material> DetectorSim::GetMaterials(const LineSegments &segments) const {
    std::set<NuGeom::Material> mats;
    for(const auto &segment : segments) { mats.insert(segment.GetMaterial()); }
    return mats;
}

std::vector<double> DetectorSim::EvaluateProbs(const LineSegments &segments,
                                               const std::map<NuGeom::Element, double> &xsecsmaps) {
    std::vector<double> probs;

    for(const auto &segment : segments) {
        auto material = segment.GetMaterial();
        double mean_free_path = 0;
        for(auto elm : material.Elements()) {
            auto cross_section = xsecsmaps.at(elm);
            mean_free_path += cross_section * material.NumberDensity(elm);
        }
        mean_free_path = 1 / mean_free_path;
        probs.push_back(segment.Length() / mean_free_path);
    }
    auto total_prob = std::accumulate(probs.begin(), probs.end(), 0.0);
    if(total_prob > max_prob) { max_prob = total_prob; }

    return probs;
}

void DetectorSim::GenerateEvents(size_t nevents) const {
    for(size_t i = 0; i < nevents; ++i) {
        NuGeom::Vector3D hit_loc;
        NuGeom::Material hit_mat;
        while(hit_mat == NuGeom::Material()) {
            auto hit = GetInteraction();
            hit_loc = hit.first;
            hit_mat = hit.second;
        }
        if(m_outfile) { m_outfile << hit_loc << ", " << hit_mat.Name() << "\n"; }
    }
}

void DetectorSim::GenerateEvents(double pot) const {
    size_t nhits = 0;
    while(m_pot < pot) {
        NuGeom::Vector3D hit_loc;
        NuGeom::Material hit_mat;
        auto hit = GetInteraction();
        hit_loc = hit.first;
        hit_mat = hit.second;
        if(hit_mat != NuGeom::Material()) {
            nhits++;
            if(m_outfile) { m_outfile << hit_loc << ", " << hit_mat.Name() << "\n"; }
        }
    }
    spdlog::info("Accumulated {} events with {} POT", nhits, m_pot);
}

// std::vector<double> Evaluate(const std::set<NuGeom::Material> &mats, double energy) {
//     std::vector<double> probs;
//     for(const auto &mat : mats) {
//         double inverse_mfp = 0;
//         for(const auto &element : mat.Elements()) {
//             inverse_mfp += xsec_callback(energy, element.PDG()) * mat.NumberDensity(element);
//         }
//         probs.push_back(1.0 / inverse_mfp);
//     }
//     return probs;
// }

std::pair<NuGeom::Vector3D, NuGeom::Material>
DetectorSim::Interaction(const LineSegments &segments,
                         const std::map<NuGeom::Element, double> &xsecsmaps) {
    // Evaluate interaction probability
    auto prob = EvaluateProbs(segments, xsecsmaps);

    // Draw random uniform [0,1)
    double r1 = NuGeom::Random::Instance().Uniform(0.0, 1.0);

    // If interaction occurs
    auto total_prob = std::accumulate(prob.begin(), prob.end(), 0.0);
    if(total_prob / max_prob < r1) {
        // No interaction occurred
        return {NuGeom::Vector3D(9e9, 9e9, 9e9), {}};
    }
    double r2 = NuGeom::Random::Instance().Uniform(0.0, total_prob);
    // Select interaction segment
    double cumulative_prob = 0.0;
    size_t idx = 0;
    for(size_t i = 0; i < prob.size(); ++i) {
        cumulative_prob += prob[i];
        if(cumulative_prob >= r2) {
            idx = i;
            break;
        }
    }
    // Determine interaction location within the segment
    double segment_length = segments[idx].Length();
    double interaction_distance = (r2 - (cumulative_prob - prob[idx])) / prob[idx] * segment_length;
    NuGeom::Vector3D interaction_point =
        segments[idx].Start() +
        (segments[idx].End() - segments[idx].Start()).Unit() * interaction_distance;
    return {interaction_point, segments[idx].GetMaterial()};
}

NuGeom::HandledRay DetectorSim::HandleRay(double energy, const NuGeom::Ray &ray) const {
    LineSegments segments = world.GetLineSegments(ray);

    // Calculate probability to interact for each line segment
    std::vector<double> probs;
    for(const auto &segment : segments) {
        auto material = segment.GetMaterial();
        auto meanfreepath = CalculateMeanFreePath(energy, material);
        probs.push_back(segment.Length() / meanfreepath);
    }

    return {probs, segments};
}

double DetectorSim::CalculateMeanFreePath(double energy, const NuGeom::Material &mat) const {
    double mean_free_path = 0;
    for(const auto &elm : mat.Elements()) {
        spdlog::trace("Element: {}", elm.PDG());
        mean_free_path += mat.NumberDensity(elm) * xsec_callback(energy, elm.PDG());
    }
    return 1.0 / mean_free_path;
}
