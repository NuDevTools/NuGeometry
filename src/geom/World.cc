#include "geom/World.hh"
#include "geom/BVH.hh"
#include "geom/BoundingBox.hh"
#include "geom/LineSegment.hh"
#include "geom/Logging.hh"
#include "geom/Ray.hh"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <string>

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

namespace {

/// Number of placed volumes in a subtree, counting the root.
size_t SubtreeNodes(const std::shared_ptr<NuGeom::PhysicalVolume> &pv) {
    size_t n = 1;
    for(const auto &d : pv->Daughters()) n += SubtreeNodes(d);
    return n;
}

/// Rebuild `pv`'s daughter list, dropping the ones that match, and recurse into
/// the survivors.  Decisions are taken against the *unpruned* total mass so the
/// threshold means the same thing at every depth.
void PruneRecursive(const std::shared_ptr<NuGeom::PhysicalVolume> &pv,
                    const NuGeom::World::PruneOptions &opts, double total_mass,
                    NuGeom::World::PruneReport &report) {
    std::vector<std::shared_ptr<NuGeom::PhysicalVolume>> kept;
    kept.reserve(pv->Daughters().size());
    const auto &mother = pv->GetLogicalVolume()->GetMaterial();

    for(const auto &d : pv->Daughters()) {
        const auto &lv = d->GetLogicalVolume();
        const std::string name = d->Name().empty() ? lv->Name() : d->Name();
        const std::string mat = lv->GetMaterial().Name();

        std::string reason;
        if(opts.keep_volumes.count(name) == 0) {
            if(opts.drop_materials.count(mat) > 0) {
                reason = "material";
            } else if(opts.min_mass_fraction > 0.0 &&
                      lv->Mass() < opts.min_mass_fraction * total_mass) {
                reason = "mass fraction";
            }
        }

        if(reason.empty()) {
            kept.push_back(d);
            PruneRecursive(d, opts, total_mass, report);
            continue;
        }

        NuGeom::World::PruneReport::Entry e;
        e.volume = name;
        e.material = mat;
        e.replaced_by = mother.Name();
        e.mass = lv->Mass();
        e.fraction = total_mass > 0 ? e.mass / total_mass : 0.0;
        // The hole left behind is filled by the mother, so the geometry's mass
        // changes by (displaced volume x mother density) - (mass removed).
        e.mass_delta = lv->GetShape()->Volume() * mother.Density() - e.mass;
        e.nodes = SubtreeNodes(d);
        e.reason = reason;

        report.removed_subtrees += 1;
        report.removed_nodes += e.nodes;
        report.removed_mass += e.mass;
        report.mass_delta += e.mass_delta;
        report.entries.push_back(std::move(e));
    }

    pv->SetOwnDaughters(std::move(kept));
}

} // namespace

NuGeom::World::PruneReport NuGeom::World::Prune(const PruneOptions &opts) {
    PruneReport report;
    if(!m_volume || !m_root_pv) return report;

    report.total_mass = m_volume->Mass();
    if(!(report.total_mass > 0.0)) {
        NuGeom::Log().warn("World::Prune: world mass is {} -- nothing pruned. A mass-fraction "
                           "threshold is meaningless without densities.",
                           report.total_mass);
        return report;
    }
    if(opts.min_mass_fraction <= 0.0 && opts.drop_materials.empty()) return report;

    PruneRecursive(m_root_pv, opts, report.total_mass, report);

    NuGeom::Log().info(
        "World::Prune: removed {} subtrees ({} placed volumes) carrying {:.4e} g = {:.4f}% of the "
        "{:.4e} g world mass; net geometry mass change {:+.4e} g ({:+.4f}%) after the mothers "
        "fill the gaps.",
        report.removed_subtrees, report.removed_nodes, report.removed_mass,
        100.0 * report.RemovedMassFraction(), report.total_mass, report.mass_delta,
        100.0 * report.MassDeltaFraction());
    for(const auto &e : report.entries) {
        NuGeom::Log().debug("World::Prune:   dropped {} ({}, {} nodes, {:.3e} g = {:.4f}%, {}) -> "
                            "space now {}",
                            e.volume, e.material, e.nodes, e.mass, 100.0 * e.fraction, e.reason,
                            e.replaced_by);
    }
    if(std::abs(report.MassDeltaFraction()) > 0.01)
        NuGeom::Log().warn("World::Prune: pruning changed the geometry mass by {:+.2f}%. That is "
                           "not a small perturbation -- the dropped volumes' space is now mother "
                           "material, which changes interaction rates.",
                           100.0 * report.MassDeltaFraction());
    return report;
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
    // The analytic boundary-sweep is the default: it matches the sequential
    // traversal ray-for-ray (cross-validated) and ROOT's TGeoNavigator, and is
    // ~1.9x faster.  Set NUGEOM_TRAVERSAL=sequential to fall back to the
    // step-by-step traversal (kept for cross-checking / debugging).  Read once.
    static const bool use_sequential = [] {
        const char *e = std::getenv("NUGEOM_TRAVERSAL");
        return e && std::string(e) == "sequential";
    }();
    return use_sequential ? GetLineSegmentsSequential(ray) : GetLineSegmentsSweep(ray);
}

std::vector<NuGeom::LineSegment> World::GetLineSegmentsSequential(const Ray &ray) const {
    std::vector<NuGeom::LineSegment> segments;
    if(!m_volume) throw std::runtime_error("NuGeom: World does not have a volume!");
    // Use the root PV's GetLineSegments which handles daughter traversal
    // and recursive unwinding via physical mothers.
    m_root_pv->GetLineSegments(ray, segments, {});
    NuGeom::Log().trace("World: Line Segments found -> {}", fmt::join(segments, ", "));
    return PruneSegments(std::move(segments));
}

namespace {

/// Shared core of the boundary sweep.  Collects every volume crossing along the
/// ray in one tree descent, orders them, and calls `emit(material, t0, t1)` for
/// each stretch of constant material.  Both the LineSegment traversal and the
/// column-density-only traversal are thin wrappers around this, so they cannot
/// drift apart.
template <typename Emit>
void SweepCore(const std::shared_ptr<NuGeom::PhysicalVolume> &root, const NuGeom::Ray &ray,
               Emit &&emit) {
    // 1. Collect every volume-crossing event along the ray in one tree descent.
    //    The root spans the whole forward ray; descendants are clipped to their
    //    parents so a protruding daughter cannot leak outside its mother.
    std::vector<NuGeom::IntervalEvent> events;
    constexpr double inf = std::numeric_limits<double>::infinity();
    root->CollectIntervals(ray, NuGeom::Transform3D{}, 0, 0.0, inf, events);
    if(events.empty()) return;

    // 2. Order boundaries along the ray.
    std::sort(
        events.begin(), events.end(),
        [](const NuGeom::IntervalEvent &a, const NuGeom::IntervalEvent &b) { return a.t < b.t; });

    // 3. Sweep.  `active` maps nesting depth -> material of the volume active at
    // that depth; properly nested geometry has at most one volume per depth
    // active at any t, so the material is the deepest active entry.  A stretch
    // spans from one boundary to the next, carrying the material in force there.
    std::map<int, const NuGeom::Material *> active;
    const NuGeom::Material *cur_mat = nullptr;
    double seg_start = 0.0;
    const size_t n = events.size();
    for(size_t i = 0; i < n;) {
        const double t = events[i].t;
        // Apply all events at this boundary before deciding the next material.
        for(; i < n && events[i].t == t; ++i) {
            const auto &e = events[i];
            if(e.delta > 0)
                active[e.depth] = e.material;
            else {
                auto it = active.find(e.depth);
                if(it != active.end() && it->second == e.material) active.erase(it);
            }
        }
        if(cur_mat && t > seg_start) emit(*cur_mat, seg_start, t);
        seg_start = t;
        cur_mat = active.empty() ? nullptr : active.rbegin()->second;
    }
}

} // namespace

std::vector<NuGeom::LineSegment> World::GetLineSegmentsSweep(const Ray &ray) const {
    if(!m_volume) throw std::runtime_error("NuGeom: World does not have a volume!");
    std::vector<NuGeom::LineSegment> segments;
    SweepCore(m_root_pv, ray, [&](const Material &mat, double t0, double t1) {
        segments.emplace_back(ray.Propagate(t0), ray.Propagate(t1), mat);
    });
    if(segments.empty()) return {};
    NuGeom::Log().trace("World(sweep): Line Segments found -> {}", fmt::join(segments, ", "));
    return PruneSegments(std::move(segments));
}

std::vector<std::pair<const NuGeom::Material *, double>>
World::GetColumnLengths(const Ray &ray) const {
    if(!m_volume) throw std::runtime_error("NuGeom: World does not have a volume!");
    // Layer-1 acceptance needs only sum(n_i sigma_i L_i), so skip building the
    // LineSegment vector and its per-boundary Propagate/Norm entirely, and
    // accumulate path length per material instead.  A ray crosses only a
    // handful of distinct materials, so the linear search beats a map.
    std::vector<std::pair<const Material *, double>> out;
    const double dir_norm = ray.Direction().Norm();
    SweepCore(m_root_pv, ray, [&](const Material &mat, double t0, double t1) {
        const double len = (t1 - t0) * dir_norm;
        for(auto &e : out) {
            if(e.first == &mat) {
                e.second += len;
                return;
            }
        }
        out.emplace_back(&mat, len);
    });
    return out;
}

std::vector<World::SweepOverlap> World::CheckSweepConsistency(const Ray &ray, bool warn) const {
    std::vector<SweepOverlap> out;
    if(!m_volume) return out;

    // Collect + sweep with volume tracking (mirrors GetLineSegmentsSweep, but
    // keeps the deepest-active *volume* so disagreements can name the culprits).
    std::vector<NuGeom::IntervalEvent> events;
    constexpr double inf = std::numeric_limits<double>::infinity();
    m_root_pv->CollectIntervals(ray, NuGeom::Transform3D{}, 0, 0.0, inf, events);
    if(events.empty()) return out;
    std::sort(
        events.begin(), events.end(),
        [](const NuGeom::IntervalEvent &a, const NuGeom::IntervalEvent &b) { return a.t < b.t; });

    std::map<int, const NuGeom::IntervalEvent *> active; // depth -> entering event
    const NuGeom::IntervalEvent *cur = nullptr;
    double seg_start = 0.0;
    const size_t n = events.size();
    for(size_t i = 0; i < n;) {
        const double t = events[i].t;
        for(; i < n && events[i].t == t; ++i) {
            const auto &e = events[i];
            if(e.delta > 0)
                active[e.depth] = &e;
            else {
                auto it = active.find(e.depth);
                if(it != active.end() && it->second->material == e.material) active.erase(it);
            }
        }
        if(cur && t > seg_start) {
            const Vector3D start = ray.Propagate(seg_start), end = ray.Propagate(t);
            if((end - start).Norm() >= 1e-4) { // ignore prune-threshold slivers
                const Vector3D mid = start + 0.5 * (end - start);
                const std::string truth_mat = FindMaterial(mid).Name();
                if(truth_mat != cur->material->Name()) {
                    const std::string sweep_vol = cur->volume ? cur->volume->Name() : "";
                    const std::string truth_vol = FindVolume(mid);
                    out.push_back(
                        {start, end, cur->material->Name(), truth_mat, sweep_vol, truth_vol});
                    if(warn)
                        NuGeom::Log().warn("Sweep/containment mismatch (possible GDML overlap): "
                                           "z=[{:.3f}, {:.3f}] "
                                           "sweep='{}'({}) but containment='{}'({})",
                                           start.Z(), end.Z(), cur->material->Name(), sweep_vol,
                                           truth_mat, truth_vol);
                }
            }
        }
        seg_start = t;
        cur = active.empty() ? nullptr : active.rbegin()->second;
    }
    return out;
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

static std::string
find_volume_recursive(const NuGeom::Vector3D &local_point,
                      const std::vector<std::shared_ptr<NuGeom::PhysicalVolume>> &daughters,
                      const std::string &parent_name) {
    for(const auto &pv : daughters) {
        if(pv->SignedDistance(local_point) <= 0) {
            auto child_point = pv->GetTransform().Apply(local_point);
            return find_volume_recursive(child_point, pv->Daughters(), pv->Name());
        }
    }
    return parent_name;
}

std::string World::FindVolume(const Vector3D &point) const {
    if(!m_volume || m_volume->GetShape()->SignedDistance(point) > 0) return "";
    return find_volume_recursive(point, m_root_pv->Daughters(), m_volume->Name());
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

    NuGeom::Log().info("Voxelize: {}x{}x{} grid, {} materials", grid.nx, grid.ny, grid.nz,
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
