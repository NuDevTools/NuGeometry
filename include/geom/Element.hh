#pragma once

#include <map>
#include <string>
#include <vector>

namespace YAML {
class Node;
}

namespace NuGeom {

struct Isotope {
    Isotope() = default;
    Isotope(const std::string &name);
    Isotope(const std::string &name, size_t z, size_t a, double mass);
    std::string m_name{};
    size_t m_Z{};
    size_t m_A{};
    double m_mass{};
    static std::map<std::string, Isotope> &CommonIsotopes() {
        static std::map<std::string, Isotope> common_isotopes;
        return common_isotopes;
    }
};

class Element {
  private:
    std::string m_name{};
    std::string m_symbol{};
    size_t m_Z{};
    size_t m_A{};
    double m_mass{};
    std::vector<std::pair<Isotope, double>> m_isotopes;

  public:
    Element() = default;
    Element(const std::string &);
    Element(const std::string &, size_t, double, size_t = 0);
    Element(const std::string &, const std::string &, size_t, double, size_t = 0);
    Element(const Element &) = default;
    Element &operator=(const Element &) = default;
    Element(Element &&) = default;
    Element &operator=(Element &&) = default;

    bool operator<(const Element &other) const {
        if(m_Z != other.m_Z)
            return m_Z < other.m_Z;
        else
            return m_A < other.m_A;
    }
    bool operator==(const Element &other) const { return m_Z == other.m_Z; }
    bool operator!=(const Element &other) const { return !(*this == other); }

    std::string Name() const { return m_name; }
    std::string Symbol() const { return m_symbol; }
    size_t PDG() const;
    size_t Z() const { return m_Z; }
    size_t A() const { return m_A; }
    double Mass() const { return m_mass; }
    void AddIsotope(const std::string &name, double fraction);

    static std::map<std::string, Element> &CommonElements() {
        static std::map<std::string, Element> common_elements;
        return common_elements;
    }

    template <typename OStream> friend OStream &operator<<(OStream &os, const Element &elm) {
        os << "Element(" << elm.m_name << ", " << elm.m_symbol << ", " << elm.m_Z << ", "
           << elm.m_mass << ")";
        return os;
    }
};

void LoadElements(const YAML::Node &);

} // namespace NuGeom
