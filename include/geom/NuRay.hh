#include <iostream>
#include <vector>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include "geom/Ray.hh"
#include "spdlog/spdlog.h"
#include "geom/Tokenizer.hh"
#include "geom/Vector3D.hh"

namespace NuGeom {

class NuRay {
public:
    NuRay() = delete;

    explicit NuRay(Ray ray, int flavor, double energy, double weight)
        : ray_(std::move(ray)), flavor_(flavor), energy_(energy), weight_(weight)
    {
        if (energy_ < 0) {
            throw std::invalid_argument("Energy must be non-negative");
        }
    }

    const Ray& GetRay() const { return ray_; }
    int GetFlavor() const { return flavor_; }
    double GetEnergy() const { return energy_; }
    double GetWeight() const { return weight_; }

private:
    Ray ray_;
    int flavor_;
    double energy_;
    double weight_;
};
}