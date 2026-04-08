#include <iostream>
#include <vector>
#include <fstream>
#include <stdexcept>
#include "geom/Ray.hh"
#include "spdlog/spdlog.h"
#include "geom/Tokenizer.hh"
#include "geom/Vector3D.hh"

namespace NuGeom {

    class NuRay {
        Ray ray_data;
        int flav;
        double Enu;
        double WGT;

    public:
        NuRay(const Ray &ray, int flav, double Enu, double WGT);
    };

    inline NuRay::NuRay(const Ray &ray, int flav, double Enu, double WGT)
        : ray_data(ray), flav(flav), Enu(Enu), WGT(WGT) {}
};

std::vector<NuGeom::NuRay> ParseCSVFlux(const std::string &filename) {
    std::ifstream file(filename);
    if(!file.is_open()) {
        spdlog::error("Could not open file: {}", filename);
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::vector<NuGeom::NuRay> nu_rays;
    std::string line;
    while(std::getline(file, line)) {
        std::vector<std::string> tokens;
        tokenize(line, tokens, ",");

        if(tokens.size() < 9) {
        spdlog::warn("Skipping malformed line: {}", line);
        continue;
    }

        int flav = std::stoi(tokens[0]);
        double Enu = std::stod(tokens[1]);
        double WGT = std::stod(tokens[2]);
        NuGeom::Ray ray = NuGeom::Ray(NuGeom::Vector3D(std::stod(tokens[3]), 
        std::stod(tokens[4]), std::stod(tokens[5])), 
        NuGeom::Vector3D(std::stod(tokens[6]), 
        std::stod(tokens[7]), std::stod(tokens[8])), WGT); 
        NuGeom::NuRay nu_ray(ray, flav, Enu, WGT);
        nu_rays.push_back(nu_ray);
    }
    return nu_rays;

}

int main() {
    auto rays = ParseCSVFlux("file.csv");
    std::cout << "Parsed " << rays.size() << " rays\n";
    return 0;
}