/// GeoDisplay — Interactive GDML Geometry Viewer
///
/// Loads a GDML file and displays it in 3D with outline rendering,
/// overlap checking, configurable visibility depth, and volume/material info.
///
/// Controls:
///   Right-mouse drag   — orbit camera
///   Scroll wheel       — zoom
///   W/A/S/D / Q/E      — pan target point
///   ImGui right panel   — vis level, overlaps, volume tree, material legend

#include "ViewerCommon.hh"
#include "geom/Parser.hh"
#include "geom/World.hh"

#include <SFML/System/Clock.hpp>
#include <imgui-SFML.h>
#include <imgui.h>
#include <misc/cpp/imgui_stdlib.h>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <map>
#include <numeric>
#include <vector>

// ---------------------------------------------------------------------------
// Overlap detection
// ---------------------------------------------------------------------------
struct OverlapPair {
    size_t idx_a, idx_b;
    NuGeom::BoundingBox overlap_bb;
};

/// Check if two AABBs overlap, returning the intersection region.
static bool aabb_overlap(const NuGeom::BoundingBox &a, const NuGeom::BoundingBox &b,
                         NuGeom::BoundingBox &out) {
    double xmin = std::max(a.min.X(), b.min.X());
    double ymin = std::max(a.min.Y(), b.min.Y());
    double zmin = std::max(a.min.Z(), b.min.Z());
    double xmax = std::min(a.max.X(), b.max.X());
    double ymax = std::min(a.max.Y(), b.max.Y());
    double zmax = std::min(a.max.Z(), b.max.Z());
    if(xmin >= xmax || ymin >= ymax || zmin >= zmax) return false;
    out = {{xmin, ymin, zmin}, {xmax, ymax, zmax}};
    return true;
}

/// Detect overlapping volumes at the same depth level using AABB checks.
static std::vector<OverlapPair>
detect_overlaps(const std::vector<NuGeom::World::VolumeBounds> &vbs) {
    std::vector<OverlapPair> result;
    for(size_t i = 0; i < vbs.size(); ++i) {
        for(size_t j = i + 1; j < vbs.size(); ++j) {
            if(vbs[i].depth != vbs[j].depth) continue;
            NuGeom::BoundingBox overlap;
            if(!aabb_overlap(vbs[i].bb, vbs[j].bb, overlap)) continue;
            result.push_back({i, j, overlap});
        }
    }
    return result;
}

// ---------------------------------------------------------------------------
// Volume tree
// ---------------------------------------------------------------------------
struct VolumeNode {
    std::string display_name;
    std::string logical_name;
    std::string material_name;
    NuGeom::BoundingBox bb;
    std::array<double, 12> transform;
    int depth;
    std::vector<VolumeNode> children;
};

/// Compute world-frame AABB from a PhysicalVolume given parent→world transform.
static NuGeom::BoundingBox compute_world_aabb(const std::shared_ptr<NuGeom::PhysicalVolume> &pv,
                                              const NuGeom::Transform3D &parent_to_world) {
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
                mn = {std::min(mn.X(), c.X()), std::min(mn.Y(), c.Y()), std::min(mn.Z(), c.Z())};
                mx = {std::max(mx.X(), c.X()), std::max(mx.Y(), c.Y()), std::max(mx.Z(), c.Z())};
            }
    return {mn, mx};
}

static VolumeNode build_tree(const std::shared_ptr<NuGeom::PhysicalVolume> &pv,
                             const NuGeom::Transform3D &parent_to_world, int depth,
                             const std::string &display_name) {
    VolumeNode node;
    node.display_name = display_name;
    node.logical_name = pv->GetLogicalVolume()->Name();
    node.material_name = pv->GetLogicalVolume()->GetMaterial().Name();
    node.depth = depth;

    node.bb = compute_world_aabb(pv, parent_to_world);

    // Global transform (local→world) = parent_to_world * placement.
    // Placement (local→parent) = m_transform.Inverse().
    auto local_to_world = parent_to_world * pv->GetTransform().Inverse();
    node.transform = local_to_world.GetTransform();

    // Count daughter name occurrences for disambiguation
    const auto &daughters = pv->Daughters();
    std::map<std::string, int> name_counts;
    for(const auto &d : daughters) name_counts[d->GetLogicalVolume()->Name()]++;

    std::map<std::string, int> name_idx;
    for(const auto &d : daughters) {
        const auto &lname = d->GetLogicalVolume()->Name();
        std::string dname = lname;
        if(name_counts[lname] > 1) { dname = lname + "_" + std::to_string(name_idx[lname]++); }
        node.children.push_back(build_tree(d, local_to_world, depth + 1, dname));
    }
    return node;
}

static VolumeNode build_root_tree(const NuGeom::World &world) {
    VolumeNode root;
    const auto &lv = world.GetLogicalVolume();
    root.display_name = lv->Name().empty() ? "World" : lv->Name();
    root.logical_name = root.display_name;
    root.material_name = lv->GetMaterial().Name();
    root.bb = world.GetWorldBox();
    root.transform = NuGeom::Transform3D{}.GetTransform();
    root.depth = -1;

    // Use the unique PV tree (root PV's own daughters) rather than the
    // original LV daughters, so the tree reflects BuildUniqueTree's structure.
    const auto &root_pv = world.GetRootPV();
    if(!root_pv) return root;
    const auto &daughters = root_pv->Daughters();
    std::map<std::string, int> name_counts;
    for(const auto &d : daughters) name_counts[d->GetLogicalVolume()->Name()]++;

    std::map<std::string, int> name_idx;
    NuGeom::Transform3D identity;
    for(const auto &d : daughters) {
        const auto &lname = d->GetLogicalVolume()->Name();
        std::string dname = lname;
        if(name_counts[lname] > 1) { dname = lname + "_" + std::to_string(name_idx[lname]++); }
        root.children.push_back(build_tree(d, identity, 0, dname));
    }
    return root;
}

// ---------------------------------------------------------------------------
// Case-insensitive substring search
// ---------------------------------------------------------------------------
static bool icontains(const std::string &haystack, const std::string &needle) {
    if(needle.empty()) return true;
    auto it = std::search(haystack.begin(), haystack.end(), needle.begin(), needle.end(),
                          [](char a, char b) { return std::tolower(a) == std::tolower(b); });
    return it != haystack.end();
}

/// Returns true if this node or any descendant matches the query.
static bool tree_matches(const VolumeNode &node, const std::string &query) {
    if(icontains(node.display_name, query) || icontains(node.material_name, query)) return true;
    for(const auto &child : node.children)
        if(tree_matches(child, query)) return true;
    return false;
}

// ---------------------------------------------------------------------------
// ImGui tree rendering (recursive)
// ---------------------------------------------------------------------------
static void draw_volume_tree(const VolumeNode &node, const std::string &search_query,
                             const VolumeNode *&selected_node,
                             const std::map<std::string, sf::Color> &palette) {
    bool is_searching = !search_query.empty();
    bool matches = is_searching && tree_matches(node, search_query);
    if(is_searching && !matches) return;

    bool self_matches = is_searching && (icontains(node.display_name, search_query) ||
                                         icontains(node.material_name, search_query));
    bool is_leaf = node.children.empty();
    bool is_selected = (selected_node == &node);

    ImGuiTreeNodeFlags flags = ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_SpanAvailWidth;
    if(is_leaf) flags |= ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;
    if(is_selected) flags |= ImGuiTreeNodeFlags_Selected;
    if(is_searching && matches) flags |= ImGuiTreeNodeFlags_DefaultOpen;

    // Material colour chip
    auto pal_it = palette.find(node.material_name);
    if(pal_it != palette.end()) {
        const auto &col = pal_it->second;
        ImGui::ColorButton(("##tc" + node.display_name).c_str(),
                           {col.r / 255.f, col.g / 255.f, col.b / 255.f, 1.f},
                           ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoBorder, {8, 8});
        ImGui::SameLine();
    }

    // Highlight matching text
    if(self_matches) ImGui::PushStyleColor(ImGuiCol_Text, {1.f, 1.f, 0.4f, 1.f});

    std::string label = node.display_name + " (" + node.material_name + ")";
    bool open = ImGui::TreeNodeEx(("##" + node.display_name + std::to_string(node.depth)).c_str(),
                                  flags, "%s", label.c_str());

    if(self_matches) ImGui::PopStyleColor();

    if(ImGui::IsItemClicked() && !ImGui::IsItemToggledOpen())
        selected_node = is_selected ? nullptr : &node;

    if(open && !is_leaf) {
        for(const auto &child : node.children)
            draw_volume_tree(child, search_query, selected_node, palette);
        ImGui::TreePop();
    }
}

/// Draw the detail panel for a selected node.
static void draw_node_info(const VolumeNode &node,
                           const std::map<std::string, NuGeom::Material> &mat_registry) {
    ImGui::TextColored({1.f, 1.f, 0.4f, 1.f}, "%s", node.display_name.c_str());
    if(node.display_name != node.logical_name)
        ImGui::Text("Logical: %s", node.logical_name.c_str());

    ImGui::Separator();
    ImGui::Text("Material: %s", node.material_name.c_str());
    auto it = mat_registry.find(node.material_name);
    if(it != mat_registry.end()) {
        const auto &mat = it->second;
        ImGui::Text("Density: %.4g g/cm\xc2\xb3", mat.Density());
        const auto &elems = mat.Elements();
        const auto &fracs = mat.MassFractions();
        for(size_t i = 0; i < elems.size(); ++i) {
            if(i < fracs.size() && fracs[i] > 0.0)
                ImGui::Text("  %s  Z=%-3zu A=%-3zu  %.2f%%", elems[i].Name().c_str(), elems[i].Z(),
                            elems[i].A(), fracs[i] * 100.0);
            else
                ImGui::Text("  %s  Z=%-3zu A=%-3zu", elems[i].Name().c_str(), elems[i].Z(),
                            elems[i].A());
        }
    }

    ImGui::Separator();
    double dx = node.bb.max.X() - node.bb.min.X();
    double dy = node.bb.max.Y() - node.bb.min.Y();
    double dz = node.bb.max.Z() - node.bb.min.Z();
    ImGui::Text("Bounding Box:");
    ImGui::Text("  min: (%.1f, %.1f, %.1f)", node.bb.min.X(), node.bb.min.Y(), node.bb.min.Z());
    ImGui::Text("  max: (%.1f, %.1f, %.1f)", node.bb.max.X(), node.bb.max.Y(), node.bb.max.Z());
    ImGui::Text("  size: %.1f x %.1f x %.1f cm", dx, dy, dz);

    ImGui::Separator();
    ImGui::Text("Global Transform (local -> MARS):");
    const auto &m = node.transform;
    ImGui::Text("  | %8.4f %8.4f %8.4f | %8.2f |", m[0], m[1], m[2], m[3]);
    ImGui::Text("  | %8.4f %8.4f %8.4f | %8.2f |", m[4], m[5], m[6], m[7]);
    ImGui::Text("  | %8.4f %8.4f %8.4f | %8.2f |", m[8], m[9], m[10], m[11]);
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    constexpr int WIN_W = 1400;
    constexpr int WIN_H = 900;
    constexpr float PANEL_W = 360.f;

    sf::ContextSettings ctx;
    ctx.depthBits = 24;
    ctx.stencilBits = 8;
    sf::RenderWindow window(sf::VideoMode({(unsigned)WIN_W, (unsigned)WIN_H}), "GeoDisplay",
                            sf::State::Windowed, ctx);
    window.setFramerateLimit(60);
    (void)window.setActive(true);
    (void)ImGui::SFML::Init(window);

    // --- persistent state ---
    std::string geomfile;
    bool has_world = false;
    NuGeom::World world;
    NuGeom::BoundingBox bbox;
    std::vector<NuGeom::World::VolumeBounds> vol_bounds;

    int vis_level = 3;
    bool outline_mode = true;
    float fill_alpha = 0.06f;
    float line_width = 1.0f;

    std::map<std::string, sf::Color> palette;
    std::map<std::string, NuGeom::Material> mat_registry;
    std::string selected_material;
    char err_buf[256] = {};

    std::vector<OverlapPair> overlaps;
    bool overlaps_checked = false;
    bool highlight_overlaps = true;

    int volume_filter_depth = -1;

    // Volume tree browser state
    VolumeNode volume_tree;
    bool has_tree = false;
    std::string tree_search;
    const VolumeNode *selected_tree_node = nullptr;

    OrbitalCamera cam;
    bool right_down = false;
    sf::Vector2i prev_mouse;
    sf::Clock delta_clock;

    // Load geometry from command line if provided
    if(argc > 1) geomfile = argv[1];

    auto rebuild_palette = [&]() {
        palette.clear();
        mat_registry.clear();
        selected_material.clear();
        for(const auto &vb : vol_bounds) {
            assign_color(vb.material.Name(), palette);
            mat_registry.emplace(vb.material.Name(), vb.material);
        }
    };

    auto load_geometry = [&]() {
        err_buf[0] = '\0';
        try {
            NuGeom::GDMLParser parser(geomfile);
            world = parser.GetWorld();
            bbox = world.GetWorldBox();
            has_world = true;
            vol_bounds = world.GetVolumeBounds(vis_level);
            overlaps.clear();
            overlaps_checked = false;
            rebuild_palette();
            centre_camera_on_bbox(cam, bbox);
            volume_tree = build_root_tree(world);
            has_tree = true;
            selected_tree_node = nullptr;
        } catch(const std::exception &e) {
            std::snprintf(err_buf, sizeof(err_buf), "%s", e.what());
            has_world = false;
            has_tree = false;
        }
    };

    if(!geomfile.empty()) load_geometry();

    while(window.isOpen()) {
        // ----------------------------------------------------------------
        // Events
        // ----------------------------------------------------------------
        while(const auto event = window.pollEvent()) {
            ImGui::SFML::ProcessEvent(window, *event);
            if(event->is<sf::Event::Closed>()) window.close();

            bool consumed = handle_camera_event(*event, cam, window, right_down, prev_mouse);
            // Show/hide cursor on right-button state change
            if(consumed) {
                if(const auto *mb = event->getIf<sf::Event::MouseButtonPressed>())
                    if(mb->button == sf::Mouse::Button::Right) window.setMouseCursorVisible(false);
                if(const auto *mb = event->getIf<sf::Event::MouseButtonReleased>())
                    if(mb->button == sf::Mouse::Button::Right) window.setMouseCursorVisible(true);
            }
        }

        ImGui::SFML::Update(window, delta_clock.restart());

        // ----------------------------------------------------------------
        // ImGui control panel
        // ----------------------------------------------------------------
        ImGui::SetNextWindowPos({WIN_W - PANEL_W, 0}, ImGuiCond_Always);
        ImGui::SetNextWindowSize({PANEL_W, static_cast<float>(WIN_H)}, ImGuiCond_Always);
        ImGui::Begin("GeoDisplay", nullptr,
                     ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize |
                         ImGuiWindowFlags_NoCollapse);

        // --- GDML file loading ---
        ImGui::Text("GDML File:");
        ImGui::SetNextItemWidth(PANEL_W - 20.f);
        ImGui::InputText("##geom", &geomfile);
        if(ImGui::Button("Load", {PANEL_W - 20.f, 0})) load_geometry();
        if(err_buf[0]) ImGui::TextColored({1.f, 0.4f, 0.4f, 1.f}, "%s", err_buf);

        if(has_world) {
            ImGui::Separator();
            ImGui::Text("World: %.0f x %.0f x %.0f cm", bbox.max.X() - bbox.min.X(),
                        bbox.max.Y() - bbox.min.Y(), bbox.max.Z() - bbox.min.Z());
            ImGui::Text("%zu volumes, %zu materials", vol_bounds.size(), palette.size());

            // --- Vis Level ---
            ImGui::Separator();
            ImGui::Text("Visibility:");
            if(ImGui::SliderInt("Vis Level", &vis_level, 0, 10)) {
                vol_bounds = world.GetVolumeBounds(vis_level);
                overlaps.clear();
                overlaps_checked = false;
                rebuild_palette();
            }
            ImGui::Checkbox("Outline mode", &outline_mode);
            if(outline_mode) ImGui::SliderFloat("Fill alpha", &fill_alpha, 0.f, 0.3f, "%.3f");
            ImGui::SliderFloat("Line width", &line_width, 0.5f, 3.f, "%.1f");

            // --- Overlap Detection ---
            ImGui::Separator();
            if(ImGui::CollapsingHeader("Overlaps", ImGuiTreeNodeFlags_DefaultOpen)) {
                if(ImGui::Button("Check Overlaps", {PANEL_W - 20.f, 0})) {
                    overlaps = detect_overlaps(vol_bounds);
                    overlaps_checked = true;
                }
                if(overlaps_checked) {
                    if(overlaps.empty()) {
                        ImGui::TextColored({0.4f, 1.f, 0.4f, 1.f}, "No AABB overlaps found.");
                    } else {
                        ImGui::TextColored({1.f, 0.7f, 0.3f, 1.f},
                                           "%zu AABB overlap(s):", overlaps.size());
                        ImGui::Checkbox("Highlight overlaps", &highlight_overlaps);
                        for(size_t i = 0; i < overlaps.size(); ++i) {
                            const auto &op = overlaps[i];
                            const auto &ma = vol_bounds[op.idx_a].material.Name();
                            const auto &mb = vol_bounds[op.idx_b].material.Name();
                            double dx = op.overlap_bb.max.X() - op.overlap_bb.min.X();
                            double dy = op.overlap_bb.max.Y() - op.overlap_bb.min.Y();
                            double dz = op.overlap_bb.max.Z() - op.overlap_bb.min.Z();
                            ImGui::Text("  [%zu] %s x %s", i, ma.c_str(), mb.c_str());
                            ImGui::Text("       %.1f x %.1f x %.1f cm", dx, dy, dz);
                        }
                    }
                }
            }

            // --- Volume Tree Browser ---
            ImGui::Separator();
            if(ImGui::CollapsingHeader("Volume Tree", ImGuiTreeNodeFlags_DefaultOpen)) {
                ImGui::SetNextItemWidth(PANEL_W - 20.f);
                ImGui::InputTextWithHint("##treesearch", "Search volumes...", &tree_search);

                ImGui::BeginChild("##voltree", {0, 250}, ImGuiChildFlags_Borders);
                if(has_tree)
                    draw_volume_tree(volume_tree, tree_search, selected_tree_node, palette);
                ImGui::EndChild();

                // Node info panel
                if(selected_tree_node) {
                    ImGui::Separator();
                    draw_node_info(*selected_tree_node, mat_registry);
                }
            }

            // --- Material Legend ---
            ImGui::Separator();
            if(ImGui::CollapsingHeader("Materials", ImGuiTreeNodeFlags_DefaultOpen)) {
                draw_material_legend(palette, mat_registry, selected_material);
            }
        }
        ImGui::End();

        // ----------------------------------------------------------------
        // 3-D OpenGL rendering
        // ----------------------------------------------------------------
        (void)window.setActive(true);

        float vw = WIN_W - PANEL_W;
        float vh = static_cast<float>(WIN_H);
        glViewport(0, 0, static_cast<GLsizei>(vw), static_cast<GLsizei>(vh));
        glClearColor(0.07f, 0.07f, 0.10f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

        cam.apply_gl(vw, vh);

        if(has_world) {
            float ex, ey, ez;
            cam.eye_pos(ex, ey, ez);
            auto order = sort_back_to_front(vol_bounds, ex, ey, ez);

            // Draw volumes
            glDepthMask(GL_FALSE);
            glLineWidth(line_width);
            for(size_t i : order) {
                const auto &vb = vol_bounds[i];
                sf::Color col = assign_color(vb.material.Name(), palette);
                float r = col.r / 255.f, g = col.g / 255.f, b = col.b / 255.f;
                if(outline_mode) {
                    float alpha = fill_alpha + 0.02f * static_cast<float>(vb.depth);
                    gl_box(vb.bb, r, g, b, alpha);
                } else {
                    gl_box(vb.bb, r, g, b, 0.f, /*wire_only=*/true);
                }
            }
            glDepthMask(GL_TRUE);

            // World box outline
            glLineWidth(1.5f);
            gl_box(bbox, 0.5f, 0.5f, 0.6f, 0.f, /*wire_only=*/true);

            // Highlight overlap regions
            if(highlight_overlaps && overlaps_checked) {
                glDepthMask(GL_FALSE);
                glLineWidth(2.f);
                for(const auto &op : overlaps) gl_box(op.overlap_bb, 1.f, 0.2f, 0.2f, 0.15f);
                glDepthMask(GL_TRUE);
            }

            // Highlight selected tree node
            if(selected_tree_node) {
                glDepthMask(GL_FALSE);
                glLineWidth(3.f);
                gl_box(selected_tree_node->bb, 1.f, 1.f, 0.2f, 0.12f);
                glDepthMask(GL_TRUE);
            }
        }

        ImGui::SFML::Render(window);
        window.display();
    }

    ImGui::SFML::Shutdown();
    return 0;
}
