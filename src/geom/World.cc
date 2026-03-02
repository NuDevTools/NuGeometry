#include "geom/World.hh"
#include "geom/BVH.hh"
#include "geom/BoundingBox.hh"
#include "geom/LineSegment.hh"
#include "geom/Ray.hh"
#include "spdlog/spdlog.h"
#include <limits>
#include <map>
#include <set>

using NuGeom::World;

// Recursively build a unique PV tree: for each daughter of pv's LogicalVolume,
// create a new PhysicalVolume that wraps the same LogicalVolume + transform
// but is a unique node with its own parent pointer.  This avoids the
// shared-LV bug where multiple placements overwrite each other's m_mother.
static void BuildUniqueTree(const std::shared_ptr<NuGeom::PhysicalVolume> &pv,
                            const std::shared_ptr<NuGeom::PhysicalVolume> &parent) {
    pv->SetMother(parent);
    for(const auto &lv_daughter : pv->GetLogicalVolume()->Daughters()) {
        auto own_child = std::make_shared<NuGeom::PhysicalVolume>(
            lv_daughter->Name(), lv_daughter->GetLogicalVolume(), lv_daughter->GetTransform(),
            0 /* tag: pre-computed transform */);
        pv->AddOwnDaughter(own_child);
        BuildUniqueTree(own_child, pv);
    }
}

World::World(std::shared_ptr<NuGeom::LogicalVolume> volume) : m_volume{std::move(volume)} {
    // Create a synthetic root PV wrapping the world LV with identity transform.
    m_root_pv = std::make_shared<NuGeom::PhysicalVolume>("WorldPV", m_volume, NuGeom::Transform3D{},
                                                         0 /* tag: pre-computed transform */);
    // Build unique PV tree: top-level daughters get m_root_pv as their mother.
    BuildUniqueTree(m_root_pv, nullptr);
}

NuGeom::Shape *World::GetShape(size_t idx) const {
    if(idx == 0)
        return m_volume->GetShape();
    else {
        const auto &daughters = m_root_pv->Daughters();
        return daughters[idx - 1]->GetLogicalVolume()->GetShape();
    }
}

NuGeom::BoundingBox World::GetWorldBox() const {
    return m_volume->GetShape()->GetBoundingBox();
}

NuGeom::Material World::GetMaterial(size_t idx) const {
    if(idx == 0)
        return m_volume->GetMaterial();
    else {
        const auto &daughters = m_root_pv->Daughters();
        return daughters[idx - 1]->GetLogicalVolume()->GetMaterial();
    }
}

std::vector<NuGeom::Material> World::GetMaterials() const {
    std::set<NuGeom::Material> material_set;
    material_set.insert(m_volume->GetMaterial());
    for(const auto &daughter : m_root_pv->Daughters()) {
        daughter->GetLogicalVolume()->GetMaterials(material_set);
    }
    std::vector<NuGeom::Material> materials(material_set.begin(), material_set.end());
    return materials;
}

bool World::InWorld(const Vector3D &pos) const {
    return m_volume->GetShape()->SignedDistance(pos) <= 0;
}

bool World::SphereTrace(const Ray &ray, double &distance, size_t &step, size_t &idx) const {
    step = 0;
    Vector3D pos = ray.Origin();
    distance = 0;
    auto res = GetSDF(pos);
    while(step < m_max_steps && std::abs(res.first) > m_epsilon && InWorld(pos)) {
        pos = ray.Propagate(distance);
        res = GetSDF(pos);
        distance += res.first;
        step++;
    }
    if(step == m_max_steps || !InWorld(pos) || res.second == 0) return false;
    idx = res.second;
    return true;
}

bool World::RayTrace(const Ray &ray, double &distance, size_t &idx) const {
    std::shared_ptr<PhysicalVolume> pvol;
    if(!m_root_pv->RayTrace(ray, distance, pvol)) return false;
    const auto &daughters = m_root_pv->Daughters();
    for(size_t i = 0; i < daughters.size(); ++i) {
        if(daughters[i] == pvol) {
            idx = i + 1;
            return true;
        }
    }
    return false;
}

// Remove sub-threshold segments and merge adjacent same-material segments.
// Traversal can produce eps-sized (≈1e-8 cm) artifacts at volume boundaries;
// any segment shorter than kPruneThreshold is an artifact, not a real feature.
static std::vector<NuGeom::LineSegment> PruneSegments(std::vector<NuGeom::LineSegment> raw) {
    static constexpr double kPruneThreshold = 1e-4; // 1 µm — below any real feature
    std::vector<NuGeom::LineSegment> out;
    for(const auto &s : raw) {
        if(s.Length() < kPruneThreshold) continue;
        if(!out.empty() && out.back().GetMaterial() == s.GetMaterial())
            out.back() = NuGeom::LineSegment(out.back().Start(), s.End(), s.GetMaterial());
        else
            out.push_back(s);
    }
    return out;
}

std::vector<NuGeom::LineSegment> World::GetLineSegments(const Ray &ray) const {
    std::vector<NuGeom::LineSegment> segments;
    if(!m_volume) throw std::runtime_error("NuGeom: World does not have a volume!");
    // Use the root PV's GetLineSegments which handles daughter traversal
    // and recursive unwinding via physical mothers.
    m_root_pv->GetLineSegments(ray, segments, {});
    spdlog::trace("World: Line Segments found -> {}", fmt::join(segments, ", "));
    return PruneSegments(std::move(segments));
}

std::pair<double, size_t> World::GetSDF(const Vector3D &pos) const {
    double distance = std::numeric_limits<double>::max();
    size_t idx = 0;
    const auto &daughters = m_root_pv->Daughters();
    for(size_t i = 0; i < daughters.size(); ++i) {
        double tmp = daughters[i]->SignedDistance(pos);
        if(tmp < distance) {
            distance = tmp;
            idx = i + 1;
        }
    }

    return {distance, idx};
}

// ---------------------------------------------------------------------------

static void collect_bounds(const std::vector<std::shared_ptr<NuGeom::PhysicalVolume>> &daughters,
                           const NuGeom::Transform3D &world_to_parent, int depth, int max_depth,
                           std::vector<NuGeom::World::VolumeBounds> &out) {
    if(depth > max_depth) return;
    auto parent_to_world = world_to_parent.Inverse();
    for(const auto &pv : daughters) {
        // GetParentBoundingBox returns AABB in the parent's local frame.
        auto bb = pv->GetParentBoundingBox();
        const NuGeom::Vector3D cx[2] = {bb.min, bb.max};
        auto c0 = parent_to_world.Apply(cx[0]);
        NuGeom::Vector3D mn = c0, mx = c0;
        for(size_t ix = 0; ix < 2; ++ix)
            for(size_t iy = 0; iy < 2; ++iy)
                for(size_t iz = 0; iz < 2; ++iz) {
                    if(!ix && !iy && !iz) continue;
                    NuGeom::Vector3D corner{cx[ix].X(), cx[iy].Y(), cx[iz].Z()};
                    auto c = parent_to_world.Apply(corner);
                    mn = {std::min(mn.X(), c.X()), std::min(mn.Y(), c.Y()),
                          std::min(mn.Z(), c.Z())};
                    mx = {std::max(mx.X(), c.X()), std::max(mx.Y(), c.Y()),
                          std::max(mx.Z(), c.Z())};
                }
        out.push_back({{mn, mx},
                       pv->GetLogicalVolume()->GetMaterial(),
                       depth,
                       pv->GetLogicalVolume()->Name()});
        // world_to_pv_local = pv->GetTransform() * world_to_parent
        collect_bounds(pv->Daughters(), pv->GetTransform() * world_to_parent, depth + 1, max_depth,
                       out);
    }
}

std::vector<NuGeom::World::VolumeBounds> World::GetVolumeBounds(int max_depth) const {
    if(!m_volume) return {};
    std::vector<VolumeBounds> result;
    collect_bounds(m_root_pv->Daughters(), NuGeom::Transform3D{}, 0, max_depth, result);
    return result;
}

// ---------------------------------------------------------------------------
// FindMaterial: walk the volume hierarchy to find the deepest volume
// containing the query point.
// ---------------------------------------------------------------------------

static NuGeom::Material
find_material_recursive(const NuGeom::Vector3D &local_point,
                        const std::vector<std::shared_ptr<NuGeom::PhysicalVolume>> &daughters,
                        const NuGeom::Material &parent_material) {
    for(const auto &pv : daughters) {
        if(pv->SignedDistance(local_point) <= 0) {
            // Point is inside this daughter — transform to daughter-local frame
            auto child_point = pv->GetTransform().Apply(local_point);
            return find_material_recursive(child_point, pv->Daughters(),
                                           pv->GetLogicalVolume()->GetMaterial());
        }
    }
    return parent_material;
}

NuGeom::Material World::FindMaterial(const Vector3D &point) const {
    if(!m_volume) return Material{};
    if(m_volume->GetShape()->SignedDistance(point) > 0) return Material{}; // outside world
    return find_material_recursive(point, m_root_pv->Daughters(), m_volume->GetMaterial());
}

// ---------------------------------------------------------------------------
// Voxelize
// ---------------------------------------------------------------------------

NuGeom::World::VoxelGrid World::Voxelize(int resolution) const {
    VoxelGrid grid;
    if(!m_volume) return grid;

    grid.bounds = GetWorldBox();
    double dx = grid.bounds.max.X() - grid.bounds.min.X();
    double dy = grid.bounds.max.Y() - grid.bounds.min.Y();
    double dz = grid.bounds.max.Z() - grid.bounds.min.Z();
    double longest = std::max({dx, dy, dz});

    grid.nx = std::max(1, static_cast<int>(resolution * dx / longest));
    grid.ny = std::max(1, static_cast<int>(resolution * dy / longest));
    grid.nz = std::max(1, static_cast<int>(resolution * dz / longest));
    grid.data.resize(static_cast<size_t>(grid.nx) * static_cast<size_t>(grid.ny) *
                     static_cast<size_t>(grid.nz));

    std::map<std::string, int16_t> mat_index;

    double sx = dx / grid.nx, sy = dy / grid.ny, sz = dz / grid.nz;
    for(int iz = 0; iz < grid.nz; ++iz)
        for(int iy = 0; iy < grid.ny; ++iy)
            for(int ix = 0; ix < grid.nx; ++ix) {
                Vector3D pt(grid.bounds.min.X() + (ix + 0.5) * sx,
                            grid.bounds.min.Y() + (iy + 0.5) * sy,
                            grid.bounds.min.Z() + (iz + 0.5) * sz);
                Material mat = FindMaterial(pt);
                auto it = mat_index.find(mat.Name());
                int16_t idx;
                if(it != mat_index.end()) {
                    idx = it->second;
                } else {
                    idx = static_cast<int16_t>(grid.materials.size());
                    mat_index[mat.Name()] = idx;
                    grid.materials.push_back(mat);
                }
                grid.data[static_cast<size_t>(ix + grid.nx * (iy + grid.ny * iz))] = idx;
            }

    spdlog::info("Voxelize: {}x{}x{} grid, {} materials", grid.nx, grid.ny, grid.nz,
                 grid.materials.size());
    return grid;
}

// ---------------------------------------------------------------------------

std::pair<double, size_t> World::GetSDFNonNeg(const Vector3D &pos) const {
    double distance = std::numeric_limits<double>::max();
    size_t idx = 0;
    const auto &daughters = m_root_pv->Daughters();
    for(size_t i = 0; i < daughters.size(); ++i) {
        double tmp = std::abs(daughters[i]->SignedDistance(pos));
        if(tmp < distance) {
            distance = tmp;
            idx = i + 1;
        }
    }

    return {distance, idx};
}
