#include "geom/NuRay.hh"
#include "geom/Ray.hh"
#include "geom/Tokenizer.hh"
#include "geom/Vector3D.hh"
#include "spdlog/spdlog.h"
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

std::vector<NuGeom::NuRay> ParseCSVFlux(const std::string &filename) {
    std::ifstream file(filename);
    if(!file.is_open()) {
        spdlog::error("Could not open file: {}", filename);
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::vector<NuGeom::NuRay> nu_rays;
    nu_rays.reserve(10000);

    std::string line;
    while(std::getline(file, line)) {
        std::vector<std::string> tokens;
        tokenize(line, tokens, ",");

        if(tokens.size() < 9) {
            spdlog::warn("Skipping malformed line: {}", line);
            continue;
        }

        try {
            int flav = std::stoi(tokens[0]);
            double WGT = std::stod(tokens[1]);

            NuGeom::Vector3D pos(std::stod(tokens[2]), std::stod(tokens[3]), std::stod(tokens[4]));

            NuGeom::Vector3D dir(std::stod(tokens[5]), std::stod(tokens[6]), std::stod(tokens[7]));

            double Enu = std::stod(tokens[8]);

            if(Enu < 0) {
                spdlog::warn("Negative energy: {}", line);
                continue;
            }

            NuGeom::Ray ray(pos, dir, WGT);
            nu_rays.emplace_back(std::move(ray), flav, Enu, WGT);

        } catch(const std::exception &e) {
            spdlog::warn("Skipping invalid line: {} ({})", line, e.what());
        }
    }

    return nu_rays;
}

double ComputeTotalWeight(const std::vector<NuGeom::NuRay> &nu_rays) {
    return std::accumulate(
        nu_rays.begin(), nu_rays.end(), 0.0,
        [](double sum, const NuGeom::NuRay &ray) { return sum + ray.GetWeight(); });
}

double ComputeWeightedEnergy(const std::vector<NuGeom::NuRay> &nu_rays) {
    return std::accumulate(nu_rays.begin(), nu_rays.end(), 0.0,
                           [](double sum, const NuGeom::NuRay &ray) {
                               return (sum + ray.GetEnergy() * ray.GetWeight());
                           });
}

int main(int argc, char *argv[]) {
    if(argc < 2) {
        spdlog::error("Usage: {} <csv_file>", argv[0]);
        return -1;
    }

    try {
        auto nu_rays = ParseCSVFlux(argv[1]);
        std::ofstream outfile("parsed_nu_rays.txt");
        for(const auto &nu_ray : nu_rays) {
            const auto &ray = nu_ray.GetRay();
            outfile << nu_ray.GetFlavor() << " " << nu_ray.GetEnergy() << " " << nu_ray.GetWeight()
                    << " " << ray.Origin().X() << " " << ray.Origin().Y() << " " << ray.Origin().Z()
                    << ray.Direction().X() << " " << ray.Direction().Y() << " "
                    << ray.Direction().Z() << "\n";
        }
        outfile.close();
        spdlog::info("Parsed {} neutrino rays from CSV", nu_rays.size());
        double total_weight = ComputeTotalWeight(nu_rays);
        spdlog::info("Total weight of parsed neutrino rays: {}", total_weight);
        double weighted_energy = ComputeWeightedEnergy(nu_rays);
        weighted_energy /= total_weight;
        spdlog::info("Weighted energy of parsed neutrino rays: {}", weighted_energy);

    } catch(const std::exception &e) {
        spdlog::error("Error parsing CSV: {}", e.what());
        return -1;
    }

    return 0;
}