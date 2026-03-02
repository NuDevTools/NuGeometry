/// Beam Ray Inspector — 3-D SFML + OpenGL + ImGui visualizer
///
/// Loads a GDML geometry, draws each daughter volume as a transparent AABB
/// box, then lets the user shoot beam rays one at a time and see each segment
/// coloured by material.
///
/// Controls:
///   Right-mouse drag   — orbit camera
///   Scroll wheel       — zoom
///   W/A/S/D / Q/E      — pan target point
///   ImGui right panel  — beam parameters, shoot/clear rays, depth filter

#include "ViewerCommon.hh"
#include "geom/DetectorSim.hh" // defines NuGeom::EnergyRay
#include "geom/LineSegment.hh"
#include "geom/Parser.hh"
#include "geom/TestGen.hh"

#include <SFML/System/Clock.hpp>
#include <imgui-SFML.h>
#include <imgui.h>
#include <misc/cpp/imgui_stdlib.h>

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <numeric>

// ---------------------------------------------------------------------------
// Ray log I/O  (format: ox oy oz dx dy dz pot, one ray per line)
// ---------------------------------------------------------------------------
static void save_ray(const std::string &path, const NuGeom::Ray &ray) {
    std::ofstream f(path, std::ios::app);
    if(!f) return;
    const auto &o = ray.Origin(), &d = ray.Direction();
    f << std::setprecision(17) << o.X() << " " << o.Y() << " " << o.Z() << " " << d.X() << " "
      << d.Y() << " " << d.Z() << " " << ray.POT() << "\n";
}

static std::vector<NuGeom::Ray> load_rays(const std::string &path, char *err_buf, size_t err_sz) {
    std::vector<NuGeom::Ray> out;
    std::ifstream f(path);
    if(!f) {
        std::snprintf(err_buf, err_sz, "Cannot open: %s", path.c_str());
        return out;
    }
    double ox, oy, oz, dx, dy, dz, pot;
    while(f >> ox >> oy >> oz >> dx >> dy >> dz >> pot)
        out.emplace_back(NuGeom::Vector3D{ox, oy, oz}, NuGeom::Vector3D{dx, dy, dz}, pot,
                         /*normalize=*/false);
    return out;
}

// ---------------------------------------------------------------------------
// Per-ray record
// ---------------------------------------------------------------------------
struct RayRecord {
    NuGeom::Ray ray{{0, 0, 0}, {0, 0, 1}, 1.0};
    std::vector<NuGeom::LineSegment> segs;
};

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main() {
    constexpr int WIN_W = 1400;
    constexpr int WIN_H = 900;
    constexpr float PANEL_W = 340.f;

    sf::ContextSettings ctx;
    ctx.depthBits = 24;
    ctx.stencilBits = 8;
    sf::RenderWindow window(sf::VideoMode({(unsigned)WIN_W, (unsigned)WIN_H}),
                            "Beam Ray Inspector 3D", sf::State::Windowed, ctx);
    window.setFramerateLimit(60);
    (void)window.setActive(true);
    (void)ImGui::SFML::Init(window);

    // --- persistent state ---
    std::string geomfile;
    bool has_world = false;
    NuGeom::World world;
    NuGeom::BoundingBox bbox;
    std::vector<NuGeom::World::VolumeBounds> vol_bounds;

    float beam_cx = 0.f, beam_cy = 0.f;
    float beam_radius = 200.f;
    float sigma_theta = 0.003f;
    int max_depth = 3;
    int nrays = 1;

    std::vector<RayRecord> ray_history;
    bool show_all_rays = true;
    bool fill_volumes = true;
    std::string raylog_path{"rays.log"};
    bool auto_log_rays = false;

    std::map<std::string, sf::Color> palette;
    std::map<std::string, NuGeom::Material> mat_registry;
    std::string selected_material;
    char err_buf[256] = {};

    // Voxelization state
    bool has_voxels = false;
    NuGeom::World::VoxelGrid voxel_grid;
    int voxel_resolution = 64;
    int slice_axis = 2; // 0=X, 1=Y, 2=Z
    int slice_index = 0;
    float voxel_opacity = 0.3f;

    OrbitalCamera cam;
    bool right_down = false;
    sf::Vector2i prev_mouse;

    sf::Clock delta_clock;

    while(window.isOpen()) {
        // ----------------------------------------------------------------
        // Events
        // ----------------------------------------------------------------
        while(const auto event = window.pollEvent()) {
            ImGui::SFML::ProcessEvent(window, *event);

            if(event->is<sf::Event::Closed>()) window.close();

            bool consumed = handle_camera_event(*event, cam, window, right_down, prev_mouse);
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
        ImGui::Begin("Controls", nullptr,
                     ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize |
                         ImGuiWindowFlags_NoCollapse);

        ImGui::Text("GDML file:");
        ImGui::SetNextItemWidth(PANEL_W - 20.f);
        ImGui::InputText("##geom", &geomfile);
        if(ImGui::Button("Load geometry", {PANEL_W - 20.f, 0})) {
            err_buf[0] = '\0';
            try {
                NuGeom::GDMLParser parser(geomfile);
                world = parser.GetWorld();
                bbox = world.GetWorldBox();
                has_world = true;
                vol_bounds = world.GetVolumeBounds(max_depth);
                ray_history.clear();
                palette.clear();
                mat_registry.clear();
                selected_material.clear();
                has_voxels = false;
                voxel_grid = {};
                // Pre-fill palette and material registry from volumes
                for(const auto &vb : vol_bounds) {
                    assign_color(vb.material.Name(), palette);
                    mat_registry.emplace(vb.material.Name(), vb.material);
                }

                centre_camera_on_bbox(cam, bbox);
                float diag =
                    static_cast<float>(std::sqrt(std::pow(bbox.max.X() - bbox.min.X(), 2) +
                                                 std::pow(bbox.max.Y() - bbox.min.Y(), 2) +
                                                 std::pow(bbox.max.Z() - bbox.min.Z(), 2)));
                beam_cx = 0.5f * static_cast<float>(bbox.min.X() + bbox.max.X());
                beam_cy = 0.5f * static_cast<float>(bbox.min.Y() + bbox.max.Y());
                beam_radius = 0.05f * diag;
            } catch(const std::exception &e) {
                std::snprintf(err_buf, sizeof(err_buf), "%s", e.what());
                has_world = false;
            }
        }
        if(err_buf[0]) ImGui::TextColored({1.f, 0.4f, 0.4f, 1.f}, "%s", err_buf);

        if(has_world) {
            ImGui::Separator();
            ImGui::Text("World: %.0f × %.0f × %.0f cm", bbox.max.X() - bbox.min.X(),
                        bbox.max.Y() - bbox.min.Y(), bbox.max.Z() - bbox.min.Z());

            if(ImGui::SliderInt("Max depth##vd", &max_depth, 0, 15)) {
                vol_bounds = world.GetVolumeBounds(max_depth);
                palette.clear();
                mat_registry.clear();
                selected_material.clear();
                for(const auto &vb : vol_bounds) {
                    assign_color(vb.material.Name(), palette);
                    mat_registry.emplace(vb.material.Name(), vb.material);
                }
                for(const auto &rv : ray_history)
                    for(const auto &s : rv.segs) {
                        assign_color(s.GetMaterial().Name(), palette);
                        mat_registry.emplace(s.GetMaterial().Name(), s.GetMaterial());
                    }
            }
            ImGui::Checkbox("Fill volumes", &fill_volumes);

            // --- Voxelization ---
            ImGui::Separator();
            ImGui::Text("Voxelization:");
            ImGui::SliderInt("Resolution##vox", &voxel_resolution, 32, 128);
            if(ImGui::Button("Voxelize", {(PANEL_W - 20.f) * 0.5f, 0})) {
                err_buf[0] = '\0';
                try {
                    voxel_grid = world.Voxelize(voxel_resolution);
                    has_voxels = true;
                    slice_index = 0;
                    for(const auto &m : voxel_grid.materials) {
                        assign_color(m.Name(), palette);
                        mat_registry.emplace(m.Name(), m);
                    }
                } catch(const std::exception &e) {
                    std::snprintf(err_buf, sizeof(err_buf), "Voxelize: %s", e.what());
                    has_voxels = false;
                }
            }
            ImGui::SameLine();
            if(ImGui::Button("Clear voxels", {-1.f, 0})) {
                has_voxels = false;
                voxel_grid = {};
            }
            if(has_voxels) {
                ImGui::RadioButton("X slice", &slice_axis, 0);
                ImGui::SameLine();
                ImGui::RadioButton("Y slice", &slice_axis, 1);
                ImGui::SameLine();
                ImGui::RadioButton("Z slice", &slice_axis, 2);
                int max_slice = (slice_axis == 0   ? voxel_grid.nx
                                 : slice_axis == 1 ? voxel_grid.ny
                                                   : voxel_grid.nz) -
                                1;
                if(slice_index > max_slice) slice_index = max_slice;
                ImGui::SliderInt("Slice##sl", &slice_index, 0, max_slice);
                ImGui::SliderFloat("Voxel opacity", &voxel_opacity, 0.f, 1.f);
                ImGui::Text("Grid: %dx%dx%d, %zu materials", voxel_grid.nx, voxel_grid.ny,
                            voxel_grid.nz, voxel_grid.materials.size());
            }

            ImGui::Separator();
            ImGui::Text("Beam Parameters:");
            ImGui::SliderFloat("Centre X", &beam_cx, static_cast<float>(bbox.min.X()),
                               static_cast<float>(bbox.max.X()));
            ImGui::SliderFloat("Centre Y", &beam_cy, static_cast<float>(bbox.min.Y()),
                               static_cast<float>(bbox.max.Y()));
            ImGui::SliderFloat("Radius (cm)", &beam_radius, 0.f,
                               0.3f * static_cast<float>(bbox.max.X() - bbox.min.X()));
            ImGui::SliderFloat("σ_θ (rad)", &sigma_theta, 0.f, 0.05f, "%.4f");
            ImGui::SliderInt("nrays", &nrays, 0, 1000);

            ImGui::Separator();
            if(ImGui::Button("Shoot one ray", {-1.f, 0})) {
                NuGeom::Vector3D center{beam_cx, beam_cy, static_cast<double>(bbox.min.Z())};
                BeamRayGen gen(0.0, 10.0, center, beam_radius, sigma_theta,
                               static_cast<double>(bbox.min.Z()));
                auto [energy, ray] = gen.GetRay();
                if(auto_log_rays) save_ray(raylog_path, ray);
                RayRecord rec;
                rec.ray = ray;
                try {
                    rec.segs = world.GetLineSegments(ray);
                    for(const auto &s : rec.segs) {
                        assign_color(s.GetMaterial().Name(), palette);
                        mat_registry.emplace(s.GetMaterial().Name(), s.GetMaterial());
                    }
                } catch(const std::exception &e) {
                    std::snprintf(err_buf, sizeof(err_buf), "GetLineSegments: %s", e.what());
                }
                ray_history.push_back(std::move(rec));
            }
            if(ImGui::Button("Shoot rays", {-1.f, 0})) {
                NuGeom::Vector3D center{beam_cx, beam_cy, static_cast<double>(bbox.min.Z())};
                BeamRayGen gen(0.0, 10.0, center, beam_radius, sigma_theta,
                               static_cast<double>(bbox.min.Z()));
                for(int i = 0; i < nrays; ++i) {
                    auto [energy, ray] = gen.GetRay();
                    if(auto_log_rays) save_ray(raylog_path, ray);
                    RayRecord rec;
                    rec.ray = ray;
                    try {
                        rec.segs = world.GetLineSegments(ray);
                        for(const auto &s : rec.segs) {
                            assign_color(s.GetMaterial().Name(), palette);
                            mat_registry.emplace(s.GetMaterial().Name(), s.GetMaterial());
                        }
                    } catch(const std::exception &e) {
                        std::snprintf(err_buf, sizeof(err_buf), "GetLineSegments: %s", e.what());
                    }
                    ray_history.push_back(std::move(rec));
                }
            }
            if(ImGui::Button("Clear all rays", {-1.f, 0})) ray_history.clear();
            ImGui::Checkbox("Show all rays", &show_all_rays);

            // --- Ray log / replay ---
            ImGui::Separator();
            ImGui::Text("Ray log file:");
            ImGui::SetNextItemWidth(PANEL_W - 70.f);
            ImGui::InputText("##rlogpath", &raylog_path);
            ImGui::SameLine();
            ImGui::Checkbox("Auto", &auto_log_rays);
            if(ImGui::IsItemHovered())
                ImGui::SetTooltip("Append every shot ray to the log file automatically");
            if(!ray_history.empty())
                if(ImGui::Button("Save latest ray", {-1.f, 0}))
                    save_ray(raylog_path, ray_history.back().ray);
            if(ImGui::Button("Replay rays from file", {-1.f, 0}) && has_world) {
                auto loaded = load_rays(raylog_path, err_buf, sizeof(err_buf));
                for(const auto &r : loaded) {
                    RayRecord rec;
                    rec.ray = r;
                    try {
                        rec.segs = world.GetLineSegments(r);
                        for(const auto &s : rec.segs) {
                            assign_color(s.GetMaterial().Name(), palette);
                            mat_registry.emplace(s.GetMaterial().Name(), s.GetMaterial());
                        }
                    } catch(const std::exception &e) {
                        std::snprintf(err_buf, sizeof(err_buf), "Replay: %s", e.what());
                    }
                    ray_history.push_back(std::move(rec));
                }
            }

            if(!ray_history.empty()) {
                const RayRecord &latest = ray_history.back();
                ImGui::Separator();
                ImGui::Text("Latest ray (%zu segs):", latest.segs.size());
                for(size_t i = 0; i < latest.segs.size(); ++i) {
                    std::string nm = latest.segs[i].GetMaterial().Name();
                    double len = latest.segs[i].Length();
                    bool bad = !std::isfinite(len);
                    sf::Color c = assign_color(nm, palette);
                    if(bad)
                        ImGui::TextColored({1.f, 0.3f, 0.3f, 1.f}, "  [%zu] %-14s *** %s ***", i,
                                           nm.c_str(), std::isinf(len) ? "inf" : "nan");
                    else
                        ImGui::TextColored({c.r / 255.f, c.g / 255.f, c.b / 255.f, 1.f},
                                           "  [%zu] %-14s %.2f cm", i, nm.c_str(), len);
                }
            }

            // Material legend
            ImGui::Separator();
            ImGui::Text("Materials (click for composition):");
            draw_material_legend(palette, mat_registry, selected_material);
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

        cam.apply_gl(vw, vh);

        if(has_world) {
            float ex, ey, ez;
            cam.eye_pos(ex, ey, ez);
            auto order = sort_back_to_front(vol_bounds, ex, ey, ez);

            // Draw volumes
            glDepthMask(GL_FALSE); // don't write depth for transparent volumes
            for(size_t i : order) {
                const auto &vb = vol_bounds[i];
                sf::Color col = assign_color(vb.material.Name(), palette);
                float r = col.r / 255.f, g = col.g / 255.f, b = col.b / 255.f;
                float alpha = fill_volumes ? (0.04f + 0.03f * static_cast<float>(vb.depth)) : 0.f;
                gl_box(vb.bb, r, g, b, alpha, !fill_volumes);
            }
            glDepthMask(GL_TRUE);

            // World box outline
            glLineWidth(1.5f);
            gl_box(bbox, 0.5f, 0.5f, 0.6f, 0.f, /*wire_only=*/true);

            // Draw voxel slice
            if(has_voxels && voxel_opacity > 0.f) {
                const auto &vb = voxel_grid.bounds;
                double vdx = (vb.max.X() - vb.min.X()) / voxel_grid.nx;
                double vdy = (vb.max.Y() - vb.min.Y()) / voxel_grid.ny;
                double vdz = (vb.max.Z() - vb.min.Z()) / voxel_grid.nz;
                int world_mat_idx = 0;

                auto draw_voxel = [&](int ix, int iy, int iz) {
                    int16_t mi = voxel_grid.at(ix, iy, iz);
                    if(mi == world_mat_idx) return;
                    const auto &mat = voxel_grid.materials[static_cast<size_t>(mi)];
                    sf::Color col = assign_color(mat.Name(), palette);
                    NuGeom::BoundingBox cell{
                        {vb.min.X() + ix * vdx, vb.min.Y() + iy * vdy, vb.min.Z() + iz * vdz},
                        {vb.min.X() + (ix + 1) * vdx, vb.min.Y() + (iy + 1) * vdy,
                         vb.min.Z() + (iz + 1) * vdz}};
                    gl_box(cell, col.r / 255.f, col.g / 255.f, col.b / 255.f, voxel_opacity);
                };

                if(slice_axis == 0) {
                    for(int iy = 0; iy < voxel_grid.ny; ++iy)
                        for(int iz = 0; iz < voxel_grid.nz; ++iz) draw_voxel(slice_index, iy, iz);
                } else if(slice_axis == 1) {
                    for(int ix = 0; ix < voxel_grid.nx; ++ix)
                        for(int iz = 0; iz < voxel_grid.nz; ++iz) draw_voxel(ix, slice_index, iz);
                } else {
                    for(int ix = 0; ix < voxel_grid.nx; ++ix)
                        for(int iy = 0; iy < voxel_grid.ny; ++iy) draw_voxel(ix, iy, slice_index);
                }
            }

            // Draw ray segments
            glLineWidth(2.5f);
            size_t first = show_all_rays ? 0 : (ray_history.empty() ? 0 : ray_history.size() - 1);
            for(size_t ri = first; ri < ray_history.size(); ++ri) {
                bool is_latest = (ri + 1 == ray_history.size());
                glLineWidth(is_latest ? 3.f : 1.5f);
                for(const auto &seg : ray_history[ri].segs) {
                    sf::Color col = assign_color(seg.GetMaterial().Name(), palette);
                    float a = is_latest ? 1.f : 0.4f;
                    glColor4f(col.r / 255.f, col.g / 255.f, col.b / 255.f, a);
                    glBegin(GL_LINES);
                    glVertex3f(static_cast<float>(seg.Start().X()),
                               static_cast<float>(seg.Start().Y()),
                               static_cast<float>(seg.Start().Z()));
                    glVertex3f(static_cast<float>(seg.End().X()), static_cast<float>(seg.End().Y()),
                               static_cast<float>(seg.End().Z()));
                    glEnd();
                }
                // Ray origin marker (small point)
                if(is_latest) {
                    const auto &o = ray_history[ri].ray.Origin();
                    glPointSize(8.f);
                    glColor4f(1.f, 1.f, 0.f, 1.f);
                    glBegin(GL_POINTS);
                    glVertex3f(static_cast<float>(o.X()), static_cast<float>(o.Y()),
                               static_cast<float>(o.Z()));
                    glEnd();
                }
            }
        }

        // ----------------------------------------------------------------
        // ImGui overlay — imgui-sfml manages its own GL state in SFML 3
        // ----------------------------------------------------------------
        ImGui::SFML::Render(window);

        window.display();
    }

    ImGui::SFML::Shutdown();
    return 0;
}
