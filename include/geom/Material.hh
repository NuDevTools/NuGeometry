#pragma once

#include "geom/Element.hh"

#include <unordered_map>
#include <vector>

namespace NuGeom {

class Material {
public:
    Material() = default; 
    Material(const std::string &name, double density, size_t ncomponents)
    : m_name{name}, m_density{density}, m_ncomponents{ncomponents} {}

    size_t NComponents() const { return m_ncomponents; }
    std::vector<Element> Elements() const { return m_elements; }
    std::vector<double> MassFractions() const { return m_fractions; }
    size_t NElements() const { return m_elements.size(); }
    void AddElement(const Element&, int);
    void AddElement(const Element&, double);
    void AddMaterial(const Material&, double);
    Element SelectElement(double) const;
    double Density() const { return m_density; }
    double NumberDensity(const Element&) const;
    std::string Name() const { return m_name; }

    template<typename OStream>
    friend OStream& operator<<(OStream &os, const Material &material) {
        os << "Material:\n";
        os << "  Name: " << material.Name() << "\n";
        os << "  Density: " << material.m_density << "\n"; 
        os << "  NComponents: " << material.m_ncomponents << "\n";
        os << "  Elements:\n";
        bool component = material.m_natoms.size() != 0;
        size_t idx = 0;
        for(const auto &elm : material.m_elements) {
            ++idx;
            if(component) {
                os << "    - " << idx << ": " << elm << " " << material.m_natoms[idx-1] << "\n";
            } else {
                os << "    - " << idx << ": " << elm << " " << material.m_fractions[idx-1] << "\n";
            }
        }
        return os;
    }

    // Comparisons
    bool operator<(const Material &other) const {
        return m_name < other.m_name;
    }
    bool operator==(const Material &other) const {
        // Probably should add some additional validation
        return m_name == other.m_name; 
    }
    bool operator!=(const Material &other) const {
        return !(*this == other);
    }

private:
    void ComputeNumberDensities();

    std::string m_name = "Dummy";
    std::vector<Element> m_elements;
    std::vector<double> m_fractions;
    std::vector<int> m_natoms;
    double m_density;
    size_t m_ncomponents;
    std::unordered_map<size_t, double> m_number_densities;

    static std::map<std::string, Material> s_materials;
};

}
