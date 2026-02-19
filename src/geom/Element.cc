#include "geom/Element.hh"
#include "geom/Units.hh"
#include "spdlog/spdlog.h"
// #include "yaml-cpp/yaml.h"

#include <cmath>
#include <stdexcept>

using NuGeom::Element;
using NuGeom::Isotope;

Isotope::Isotope(const std::string &name) {
    if(CommonIsotopes().find(name) == CommonIsotopes().end()) {
        throw std::runtime_error("Invalid isotope " + name);
    }
    *this = CommonIsotopes().at(name);
}

Isotope::Isotope(const std::string &name, size_t z, size_t a, double mass)
    : m_name{name}, m_Z{z}, m_A{a}, m_mass{mass} {
    if(CommonIsotopes().find(name) == CommonIsotopes().end()) { CommonIsotopes()[name] = *this; }
}

Element::Element(const std::string &name) {
    if(CommonElements().find(name) == CommonElements().end()) {
        throw std::runtime_error("Invalid element " + name);
    }
    *this = CommonElements().at(name);
}

Element::Element(const std::string &name, size_t Z, double mass, size_t A)
    : m_name{name}, m_Z{Z}, m_mass{mass} {
    if(A == 0) {
        m_A = static_cast<size_t>(mass);
    } else {
        m_A = A;
    }

    if(CommonElements().find(name) == CommonElements().end()) { CommonElements()[name] = *this; }
}

Element::Element(const std::string &name, const std::string &symbol, size_t Z, double mass,
                 size_t A)
    : m_name{name}, m_symbol{symbol}, m_Z{Z}, m_mass{mass} {
    if(A == 0) {
        m_A = static_cast<size_t>(mass);
    } else {
        m_A = A;
    }

    if(CommonElements().find(name) == CommonElements().end()) { CommonElements()[name] = *this; }

    if(CommonElements().find(symbol) == CommonElements().end()) {
        CommonElements()[symbol] = *this;
    }
}

void Element::AddIsotope(const std::string &name, double fraction) {
    Isotope iso(name);
    if(m_Z == 0)
        m_Z = iso.m_Z;
    else if(m_Z != iso.m_Z) { throw std::runtime_error("Isotopes must all have the same Z!"); }
    m_isotopes.push_back({Isotope(name), fraction});
    m_mass = 0;
    double frac_sum = 0;
    spdlog::trace("Element: {}", m_name);
    for(const auto &iso_frac : m_isotopes) {
        spdlog::trace("Adding mass = {} with frac = {}", iso_frac.first.m_mass, iso_frac.second);
        m_mass += iso_frac.first.m_mass * iso_frac.second;
        frac_sum += iso_frac.second;
    }

    if(frac_sum - 1e-4 > 1) { throw std::runtime_error("Isotope fractions sum larger than 1!"); }

    // Update element database
    CommonElements()[m_name] = *this;
}

size_t Element::PDG() const {
    static constexpr size_t base = 1000000000;
    static constexpr size_t zbase = 10000;
    static constexpr size_t abase = 10;
    return base + m_Z * zbase + m_A * abase;
}

// void NuGeom::LoadElements(const YAML::Node &node) {
//     for(const auto &subnode : node) {
//         auto name = subnode[0].as<std::string>();
//         auto symbol = subnode[1].as<std::string>();
//         auto z = subnode[2].as<size_t>();
//         auto a = subnode[3].as<size_t>();
//         auto m = subnode[4].as<double>();
//
//         Element elm(name, symbol, z, a, m);
//         Element::CommonElements()[name] = elm;
//         Element::CommonElements()[symbol] = elm;
//     }
// }
