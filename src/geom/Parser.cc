#include "geom/Parser.hh"
#include "spdlog/spdlog.h"
#include <cmath>
#include <cstring>

using NuGeom::GDMLParser;

GDMLParser::GDMLParser(const std::string &filename) {
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filename.c_str());
    if(!result) throw std::runtime_error("GDMLParser: Invalid geometry file");

    *this = GDMLParser(doc);
}

GDMLParser::GDMLParser(const pugi::xml_document &doc) {
    auto root = doc.child("gdml");
    ParseDefines(root.child("define"));
    ParseMaterials(root.child("materials"));
    ParseSolids(root.child("solids"));
    ParseStructure(root.child("structure"));

    auto setup = root.child("setup");
    m_world = World(m_volumes[setup.child("world").attribute("ref").value()]);

    spdlog::info("Number of constants defined: {}", m_def_constants.size());
    spdlog::info("Number of positions defined: {}", m_def_positions.size());
    spdlog::info("Number of rotations defined: {}", m_def_rotations.size());
    spdlog::info("Number of materials: {}", m_materials.size());
    spdlog::info("Number of solids: {}", m_shapes.size());
    spdlog::info("Number of volumes: {}", m_volumes.size());
    spdlog::info("Number of physical volumes: {}", m_phys_vols.size());
}

double GDMLParser::GetConstant(const std::string &name) const {
    if(m_def_constants.find(name) == m_def_constants.end())
        throw std::runtime_error(fmt::format("GDMLParser: Undefined constant {}", name));

    return m_def_constants.at(name);
}

NuGeom::Vector3D GDMLParser::GetPosition(const std::string &name) const {
    if(m_def_positions.find(name) == m_def_positions.end())
        throw std::runtime_error(fmt::format("GDMLParser: Undefined position {}", name));

    return m_def_positions.at(name);
}

std::shared_ptr<NuGeom::Shape> GDMLParser::GetShape(const std::string &name) const {
    if(m_shapes.find(name) == m_shapes.end())
        throw std::runtime_error(fmt::format("GDMLParser: Undefined shape {}", name));

    return m_shapes.at(name);
}

NuGeom::Transform3D GDMLParser::GetTransform(const std::string &name) const {
    if(m_def_rotations.find(name) == m_def_rotations.end())
        throw std::runtime_error(fmt::format("GDMLParser: Undefined transform {}", name));

    return m_def_rotations.at(name);
}

NuGeom::Material GDMLParser::GetMaterial(const std::string &name) const {
    if(m_materials.find(name) == m_materials.end())
        throw std::runtime_error(fmt::format("GDMLParser: Undefined material {}", name));

    return m_materials.at(name);
}

std::vector<NuGeom::Material> GDMLParser::GetMaterials() const {
    std::vector<NuGeom::Material> mats;
    for(const auto &mat : m_materials) mats.push_back(mat.second);
    return mats;
}

void GDMLParser::ParseDefines(const pugi::xml_node &define) {
    for(const auto &node : define.children("constant")) {
        std::string name = node.attribute("name").value();
        double value = node.attribute("value").as_double();
        m_def_constants[name] = value;
    }

    for(const auto &node : define.children("position")) {
        // Load the position information
        std::string name = node.attribute("name").value();
        double x = node.attribute("x").as_double();
        double y = node.attribute("y").as_double();
        double z = node.attribute("z").as_double();
        Vector3D position(x, y, z);

        // Convert the units
        std::string unit = node.attribute("unit").value();
        if(unit == "m") {
            position *= 100;
        } else if(unit == "mm") {
            position /= 10;
        } else if(unit != "cm") {
            std::cout << name << std::endl;
            throw std::runtime_error(
                fmt::format("GDMLParser: Invalid position unit found ({})", unit));
        }
        m_def_positions[name] = position;
    }

    for(const auto &node : define.children("rotation")) {
        // Load the rotation information
        std::string name = node.attribute("name").value();
        std::string unit = node.attribute("unit").value();
        double xRot = node.attribute("x").as_double();
        double yRot = node.attribute("y").as_double();
        double zRot = node.attribute("z").as_double();

        // Convert if needed
        double convert{};
        if(unit == "deg" || unit == "degree")
            convert = M_PI / 180;
        else if(unit == "rad")
            convert = 1;
        else {
            std::cout << name << std::endl;
            throw std::runtime_error(
                fmt::format("GDMLParser: Invalid angle unit found ({})", unit));
        }
        auto rotX = RotationX3D(xRot * convert);
        auto rotY = RotationY3D(yRot * convert);
        auto rotZ = RotationZ3D(zRot * convert);
        Transform3D rot = rotZ * rotY * rotX;
        m_def_rotations[name] = rot;
    }
}

void GDMLParser::ParseMaterials(const pugi::xml_node &materials) {
    // Parse all isotopes
    for(const auto &node : materials.children("isotope")) {
        std::string name = node.attribute("name").value();
        size_t z = node.attribute("Z").as_ullong();
        size_t a = node.attribute("N").as_ullong();
        double mass = node.child("atom").attribute("value").as_double();
        Isotope iso(name, z, a, mass);
    }

    // Parse all elements
    for(const auto &node : materials.children("element")) {
        std::string name = node.attribute("name").value();
        std::string symbol = node.attribute("formula").value();
        size_t z = node.attribute("Z").as_ullong();
        if(node.child("atom")) {
            double mass = node.child("atom").attribute("value").as_double();
            Element elm(name, symbol, z, mass);
        } else if(node.child("fraction")) {
            Element elm(name, 0, 0);
            for(const auto &isotope : node.children("fraction")) {
                elm.AddIsotope(isotope.attribute("ref").value(),
                               isotope.attribute("n").as_double());
            }
        } else {
            throw std::runtime_error("Element: Invalid format");
        }
    }

    // Parse all materials
    for(const auto &node : materials.children("material")) {
        std::string name = node.attribute("name").value();
        spdlog::debug("Parsing material {}", name);
        double density = node.child("D").attribute("value").as_double();
        std::string unit = "g/cm3";
        if(node.child("D").attribute("unit")) unit = node.child("D").attribute("unit").value();

        // TODO: Refactor this to make it cleaner
        if(node.child("fraction")) {
            auto nelements = static_cast<size_t>(
                std::distance(node.children("fraction").begin(), node.children("fraction").end()));

            Material material(name, density, nelements);
            for(const auto &element : node.children("fraction")) {
                double fraction = element.attribute("n").as_double();
                // TODO: Allow other materials to be added to a new material
                if(m_materials.find(element.attribute("ref").as_string()) != m_materials.end()) {
                    material.AddMaterial(m_materials[element.attribute("ref").as_string()],
                                         fraction);
                } else {
                    Element elm(element.attribute("ref").as_string());
                    material.AddElement(elm, fraction);
                }
            }
            m_materials[name] = material;
        } else if(node.child("composite")) {
            auto nelements = static_cast<size_t>(std::distance(node.children("composite").begin(),
                                                               node.children("composite").end()));
            Material material(name, density, nelements);
            for(const auto &element : node.children("composite")) {
                auto natoms = element.attribute("n").as_double();
                // TODO: Allow other materials to be added to a new material
                if(m_materials.find(element.attribute("ref").as_string()) != m_materials.end()) {
                    throw std::runtime_error(
                        "GDMLParser: Using composite materials requires mass fractions");
                } else {
                    Element elm(element.attribute("ref").as_string());
                    if(natoms < 1) {
                        material.AddElement(elm, natoms);
                    } else {
                        material.AddElement(elm, static_cast<int>(std::floor(natoms)));
                    }
                }
            }
            m_materials[name] = material;
        } else {
            throw std::runtime_error("GDMLParser: Invalid material");
        }
    }
}

/// Parse the position/rotation of the second operand in a CSG solid.
/// Returns a parent-to-local transform for the second shape.
static NuGeom::Transform3D
ParseCSGTransform(const pugi::xml_node &solid,
                  const std::map<std::string, NuGeom::Vector3D> &positions,
                  const std::map<std::string, NuGeom::Transform3D> &rotations) {
    using namespace NuGeom;
    Vector3D pos;
    Transform3D rot;

    if(solid.child("positionref")) {
        std::string ref = solid.child("positionref").attribute("ref").value();
        auto it = positions.find(ref);
        if(it != positions.end()) pos = it->second;
    } else if(solid.child("position")) {
        auto pos_node = solid.child("position");
        double x = pos_node.attribute("x").as_double();
        double y = pos_node.attribute("y").as_double();
        double z = pos_node.attribute("z").as_double();
        pos = Vector3D(x, y, z);
        std::string unit = pos_node.attribute("unit").value();
        if(unit.empty()) unit = pos_node.attribute("lunit").value();
        if(unit == "m")
            pos *= 100;
        else if(unit == "mm")
            pos /= 10;
    }

    if(solid.child("rotationref")) {
        std::string ref = solid.child("rotationref").attribute("ref").value();
        auto it = rotations.find(ref);
        if(it != rotations.end()) rot = it->second;
    } else if(solid.child("rotation")) {
        auto rot_node = solid.child("rotation");
        double xRot = rot_node.attribute("x").as_double();
        double yRot = rot_node.attribute("y").as_double();
        double zRot = rot_node.attribute("z").as_double();
        double convert = 1;
        std::string unit = rot_node.attribute("unit").value();
        if(unit.empty()) unit = rot_node.attribute("aunit").value();
        if(unit == "deg" || unit == "degree") convert = M_PI / 180;
        auto rotX = RotationX3D(xRot * convert);
        auto rotY = RotationY3D(yRot * convert);
        auto rotZ = RotationZ3D(zRot * convert);
        rot = rotZ * rotY * rotX;
    }

    return (rot * Translation3D(pos)).Inverse();
}

void GDMLParser::ParseSolids(const pugi::xml_node &solids) {
    for(const auto &solid : solids) {
        std::string name = solid.attribute("name").value();
        spdlog::debug("Parsing solid {} of type {}", name, solid.name());
        if(std::strcmp(solid.name(), "subtraction") == 0 ||
           std::strcmp(solid.name(), "union") == 0 ||
           std::strcmp(solid.name(), "intersection") == 0) {
            std::string first_name = solid.child("first").attribute("ref").value();
            std::string second_name = solid.child("second").attribute("ref").value();
            auto first_shape = m_shapes[first_name];
            auto second_shape = m_shapes[second_name];

            // Parse position/rotation of second operand relative to first
            auto second_transform = ParseCSGTransform(solid, m_def_positions, m_def_rotations);
            bool has_transform = !second_transform.IsIdentity();
            std::shared_ptr<Shape> positioned_second =
                has_transform ? std::make_shared<TransformedShape>(second_shape, second_transform)
                              : second_shape;

            ShapeBinaryOp op;
            if(std::strcmp(solid.name(), "subtraction") == 0) {
                op = ShapeBinaryOp::kSubtraction;
                // CombinedShape convention: right - left. GDML: first - second.
                // So pass (second, first) to get first - second.
                m_shapes[name] =
                    std::make_shared<CombinedShape>(positioned_second, first_shape, op);
            } else if(std::strcmp(solid.name(), "union") == 0) {
                op = ShapeBinaryOp::kUnion;
                m_shapes[name] =
                    std::make_shared<CombinedShape>(first_shape, positioned_second, op);
            } else {
                op = ShapeBinaryOp::kIntersect;
                m_shapes[name] =
                    std::make_shared<CombinedShape>(first_shape, positioned_second, op);
            }
        } else {
            std::shared_ptr<Shape> shape = ShapeFactory::Initialize(solid.name(), solid);
            m_shapes[name] = shape;
        }
    }
}

void GDMLParser::ParseStructure(const pugi::xml_node &structure) {
    for(const auto &node : structure.children("volume")) {
        std::string name = node.attribute("name").value();
        std::string material_ref = node.child("materialref").attribute("ref").value();
        std::string solid_ref = node.child("solidref").attribute("ref").value();
        Material material = m_materials[material_ref];
        auto shape = m_shapes[solid_ref];
        auto volume = std::make_shared<LogicalVolume>(name, material, shape);

        // Check for sub-volumes
        for(const auto &subnode : node.children("physvol")) {
            std::string volume_ref = subnode.child("volumeref").attribute("ref").value();

            Vector3D translation;
            if(subnode.child("positionref")) {
                std::string position_ref = subnode.child("positionref").attribute("ref").value();
                translation = m_def_positions[position_ref];
            } else if(subnode.child("position")) {
                auto pos_node = subnode.child("position");
                double x = pos_node.attribute("x").as_double();
                double y = pos_node.attribute("y").as_double();
                double z = pos_node.attribute("z").as_double();
                translation = Vector3D(x, y, z);

                // Convert the units
                std::string unit = pos_node.attribute("unit").value();
                if(unit.empty()) unit = pos_node.attribute("lunit").value();
                if(unit == "m") {
                    translation *= 100;
                } else if(unit == "mm") {
                    translation /= 10;
                }
            }

            Transform3D rotation;
            if(subnode.child("rotationref")) {
                rotation = m_def_rotations[subnode.child("rotationref").attribute("ref").value()];
            } else if(subnode.child("rotation")) {
                auto rot_node = subnode.child("rotation");
                double xRot = rot_node.attribute("x").as_double();
                double yRot = rot_node.attribute("y").as_double();
                double zRot = rot_node.attribute("z").as_double();

                // Convert if needed
                double convert = 1;
                std::string unit = rot_node.attribute("unit").value();
                if(unit.empty()) unit = rot_node.attribute("aunit").value();
                if(unit == "deg" || unit == "degree")
                    convert = M_PI / 180;
                else if(unit != "rad" && !unit.empty())
                    throw std::runtime_error("GDMLParser: Invalid angle unit: " + unit);
                auto rotX = RotationX3D(xRot * convert);
                auto rotY = RotationY3D(yRot * convert);
                auto rotZ = RotationZ3D(zRot * convert);
                rotation = rotZ * rotY * rotX;
            }

            auto subvolume = m_volumes[volume_ref];
            auto phys_vol = std::make_shared<PhysicalVolume>(volume_ref, subvolume,
                                                             Translation3D(translation), rotation);
            m_phys_vols.push_back(phys_vol);
            volume->AddDaughter(m_phys_vols.back());
        }

        spdlog::info("Volume: {}", name);
        // spdlog::info("  Mass = {}", volume->Mass());
        //  Store volume information
        m_volumes[name] = volume;
    }
}
