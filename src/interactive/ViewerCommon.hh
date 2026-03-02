#pragma once
/// ViewerCommon.hh — Shared utilities for SFML+OpenGL+ImGui geometry viewers.
///
/// Provides: colour palette, OrbitalCamera, gl_box, and SFML event helpers.

#include "geom/BoundingBox.hh"
#include "geom/World.hh"

#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>
#include <SFML/Window/Event.hpp>
#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Colour palette
// ---------------------------------------------------------------------------
inline const sf::Color kPalette[] = {
    {0x4F, 0x9D, 0xFF}, {0xFF, 0x6B, 0x6B}, {0x6B, 0xFF, 0x9D}, {0xFF, 0xDB, 0x58},
    {0xDA, 0x70, 0xD6}, {0xFF, 0xA5, 0x00}, {0x00, 0xE5, 0xFF}, {0xFF, 0xC0, 0xCB},
};
inline constexpr int kNColors = 8;

inline sf::Color assign_color(const std::string &name, std::map<std::string, sf::Color> &palette) {
    auto it = palette.find(name);
    if(it != palette.end()) return it->second;
    sf::Color col = kPalette[static_cast<int>(palette.size()) % kNColors];
    palette[name] = col;
    return col;
}

// ---------------------------------------------------------------------------
// Camera — orbital
// ---------------------------------------------------------------------------
struct OrbitalCamera {
    float az{30.f}, el{20.f}, dist{1000.f};
    float tx{0.f}, ty{0.f}, tz{0.f};

    void orbit(float daz, float del) {
        az += daz;
        el = std::max(-89.f, std::min(89.f, el + del));
    }
    void zoom(float delta) {
        dist *= std::pow(0.85f, delta);
        dist = std::max(1.f, dist);
    }
    void pan(float right, float up) {
        float azr = az * static_cast<float>(M_PI) / 180.f;
        tx -= right * std::cos(azr);
        tz -= right * std::sin(azr);
        ty += up;
    }

    void eye_pos(float &ex, float &ey, float &ez) const {
        float azr = az * static_cast<float>(M_PI) / 180.f;
        float elr = el * static_cast<float>(M_PI) / 180.f;
        ex = tx + dist * std::cos(elr) * std::sin(azr);
        ey = ty + dist * std::sin(elr);
        ez = tz + dist * std::cos(elr) * std::cos(azr);
    }

    /// Upload a perspective projection and look-at modelview to the current GL context.
    void apply_gl(float w, float h) const {
        float ex, ey, ez;
        eye_pos(ex, ey, ez);

        // Projection
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        float fov = 60.f * static_cast<float>(M_PI) / 180.f;
        float f = 1.f / std::tan(fov * 0.5f);
        float asp = w / h;
        float zn = dist * 0.001f, zf = dist * 10.f;
        float pm[16] = {f / asp,
                        0,
                        0,
                        0,
                        0,
                        f,
                        0,
                        0,
                        0,
                        0,
                        (zf + zn) / (zn - zf),
                        -1,
                        0,
                        0,
                        2 * zf * zn / (zn - zf),
                        0};
        glLoadMatrixf(pm);

        // Modelview (look-at)
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        float fx = tx - ex, fy = ty - ey, fz = tz - ez;
        float fn = std::sqrt(fx * fx + fy * fy + fz * fz);
        fx /= fn;
        fy /= fn;
        fz /= fn;
        float sx = fy * 0.f - fz * 1.f, sy = fz * 0.f - fx * 0.f, sz = fx * 1.f - fy * 0.f;
        float sn = std::sqrt(sx * sx + sy * sy + sz * sz);
        sx /= sn;
        sy /= sn;
        sz /= sn;
        float ux = sy * fz - sz * fy, uy = sz * fx - sx * fz, uz = sx * fy - sy * fx;
        float mv[16] = {sx,
                        ux,
                        -fx,
                        0,
                        sy,
                        uy,
                        -fy,
                        0,
                        sz,
                        uz,
                        -fz,
                        0,
                        -(sx * ex + sy * ey + sz * ez),
                        -(ux * ex + uy * ey + uz * ez),
                        fx * ex + fy * ey + fz * ez,
                        1};
        glLoadMatrixf(mv);
    }
};

// ---------------------------------------------------------------------------
// OpenGL drawing helpers
// ---------------------------------------------------------------------------

/// Draw a box as wireframe edges + optional transparent fill.
inline void gl_box(const NuGeom::BoundingBox &bb, float r, float g, float b, float alpha,
                   bool wire_only = false) {
    float x0 = static_cast<float>(bb.min.X()), x1 = static_cast<float>(bb.max.X());
    float y0 = static_cast<float>(bb.min.Y()), y1 = static_cast<float>(bb.max.Y());
    float z0 = static_cast<float>(bb.min.Z()), z1 = static_cast<float>(bb.max.Z());

    if(!wire_only) {
        glColor4f(r, g, b, alpha);
        glBegin(GL_QUADS);
        // -z face
        glVertex3f(x0, y0, z0);
        glVertex3f(x1, y0, z0);
        glVertex3f(x1, y1, z0);
        glVertex3f(x0, y1, z0);
        // +z face
        glVertex3f(x0, y0, z1);
        glVertex3f(x0, y1, z1);
        glVertex3f(x1, y1, z1);
        glVertex3f(x1, y0, z1);
        // -x face
        glVertex3f(x0, y0, z0);
        glVertex3f(x0, y1, z0);
        glVertex3f(x0, y1, z1);
        glVertex3f(x0, y0, z1);
        // +x face
        glVertex3f(x1, y0, z0);
        glVertex3f(x1, y0, z1);
        glVertex3f(x1, y1, z1);
        glVertex3f(x1, y1, z0);
        // -y face
        glVertex3f(x0, y0, z0);
        glVertex3f(x0, y0, z1);
        glVertex3f(x1, y0, z1);
        glVertex3f(x1, y0, z0);
        // +y face
        glVertex3f(x0, y1, z0);
        glVertex3f(x1, y1, z0);
        glVertex3f(x1, y1, z1);
        glVertex3f(x0, y1, z1);
        glEnd();
    }

    // Wireframe
    glColor4f(r, g, b, std::min(1.f, alpha + 0.3f));
    glBegin(GL_LINES);
    // bottom ring
    glVertex3f(x0, y0, z0);
    glVertex3f(x1, y0, z0);
    glVertex3f(x1, y0, z0);
    glVertex3f(x1, y0, z1);
    glVertex3f(x1, y0, z1);
    glVertex3f(x0, y0, z1);
    glVertex3f(x0, y0, z1);
    glVertex3f(x0, y0, z0);
    // top ring
    glVertex3f(x0, y1, z0);
    glVertex3f(x1, y1, z0);
    glVertex3f(x1, y1, z0);
    glVertex3f(x1, y1, z1);
    glVertex3f(x1, y1, z1);
    glVertex3f(x0, y1, z1);
    glVertex3f(x0, y1, z1);
    glVertex3f(x0, y1, z0);
    // verticals
    glVertex3f(x0, y0, z0);
    glVertex3f(x0, y1, z0);
    glVertex3f(x1, y0, z0);
    glVertex3f(x1, y1, z0);
    glVertex3f(x1, y0, z1);
    glVertex3f(x1, y1, z1);
    glVertex3f(x0, y0, z1);
    glVertex3f(x0, y1, z1);
    glEnd();
}

// ---------------------------------------------------------------------------
// SFML event handling helpers
// ---------------------------------------------------------------------------

/// Process orbital camera events (right-drag orbit, scroll zoom, WASDQE pan).
/// Call this inside your event polling loop. Returns true if the event was consumed.
inline bool handle_camera_event(const sf::Event &event, OrbitalCamera &cam,
                                const sf::RenderWindow &window, bool &right_down,
                                sf::Vector2i &prev_mouse) {
    if(!ImGui::GetIO().WantCaptureMouse) {
        if(const auto *mb = event.getIf<sf::Event::MouseButtonPressed>()) {
            if(mb->button == sf::Mouse::Button::Right) {
                right_down = true;
                prev_mouse = sf::Mouse::getPosition(window);
                return true;
            }
        }
        if(const auto *mb = event.getIf<sf::Event::MouseButtonReleased>()) {
            if(mb->button == sf::Mouse::Button::Right) {
                right_down = false;
                return true;
            }
        }
        if(event.is<sf::Event::MouseMoved>() && right_down) {
            auto cur = sf::Mouse::getPosition(window);
            cam.orbit(static_cast<float>(cur.x - prev_mouse.x) * 0.4f,
                      static_cast<float>(cur.y - prev_mouse.y) * 0.4f);
            prev_mouse = cur;
            return true;
        }
        if(const auto *mw = event.getIf<sf::Event::MouseWheelScrolled>()) {
            cam.zoom(mw->delta);
            return true;
        }
    }

    if(!ImGui::GetIO().WantCaptureKeyboard) {
        if(const auto *kp = event.getIf<sf::Event::KeyPressed>()) {
            float step = cam.dist * 0.02f;
            if(kp->code == sf::Keyboard::Key::W) {
                cam.pan(0, step);
                return true;
            }
            if(kp->code == sf::Keyboard::Key::S) {
                cam.pan(0, -step);
                return true;
            }
            if(kp->code == sf::Keyboard::Key::A) {
                cam.pan(step, 0);
                return true;
            }
            if(kp->code == sf::Keyboard::Key::D) {
                cam.pan(-step, 0);
                return true;
            }
            if(kp->code == sf::Keyboard::Key::Q) {
                cam.ty -= step;
                return true;
            }
            if(kp->code == sf::Keyboard::Key::E) {
                cam.ty += step;
                return true;
            }
            if(kp->code == sf::Keyboard::Key::Z) {
                cam.zoom(5);
                return true;
            }
            if(kp->code == sf::Keyboard::Key::X) {
                cam.zoom(-5);
                return true;
            }
        }
    }
    return false;
}

// ---------------------------------------------------------------------------
// Back-to-front volume sorting
// ---------------------------------------------------------------------------

/// Return indices into vol_bounds sorted back-to-front from the camera eye position.
inline std::vector<size_t>
sort_back_to_front(const std::vector<NuGeom::World::VolumeBounds> &vol_bounds, float ex, float ey,
                   float ez) {
    auto dist_sq = [&](const NuGeom::BoundingBox &bb) {
        float cx = 0.5f * static_cast<float>(bb.min.X() + bb.max.X()) - ex;
        float cy = 0.5f * static_cast<float>(bb.min.Y() + bb.max.Y()) - ey;
        float cz = 0.5f * static_cast<float>(bb.min.Z() + bb.max.Z()) - ez;
        return cx * cx + cy * cy + cz * cz;
    };
    std::vector<size_t> order(vol_bounds.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return dist_sq(vol_bounds[a].bb) > dist_sq(vol_bounds[b].bb);
    });
    return order;
}

// ---------------------------------------------------------------------------
// Material legend ImGui widget
// ---------------------------------------------------------------------------

/// Draw a clickable colour-coded material legend. Updates selected_material on click.
inline void draw_material_legend(std::map<std::string, sf::Color> &palette,
                                 const std::map<std::string, NuGeom::Material> &mat_registry,
                                 std::string &selected_material) {
    for(const auto &[nm, col] : palette) {
        ImGui::ColorButton(("##c" + nm).c_str(), {col.r / 255.f, col.g / 255.f, col.b / 255.f, 1.f},
                           ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoBorder, {10, 10});
        ImGui::SameLine();
        bool is_sel = (selected_material == nm);
        if(ImGui::Selectable(nm.c_str(), is_sel, 0, {0, 0})) selected_material = is_sel ? "" : nm;
    }

    // Element composition for selected material
    if(!selected_material.empty()) {
        auto it = mat_registry.find(selected_material);
        if(it != mat_registry.end()) {
            const NuGeom::Material &mat = it->second;
            ImGui::Separator();
            ImGui::TextColored({1.f, 1.f, 0.4f, 1.f}, "%s", selected_material.c_str());
            ImGui::Text("Density: %.4g g/cm\xc2\xb3", mat.Density());
            const auto &elems = mat.Elements();
            const auto &fracs = mat.MassFractions();
            for(size_t i = 0; i < elems.size(); ++i) {
                const auto &el = elems[i];
                if(i < fracs.size() && fracs[i] > 0.0)
                    ImGui::Text("  %s  Z=%-3zu A=%-3zu  %.2f%%", el.Name().c_str(), el.Z(), el.A(),
                                fracs[i] * 100.0);
                else
                    ImGui::Text("  %s  Z=%-3zu A=%-3zu", el.Name().c_str(), el.Z(), el.A());
            }
        }
    }
}

/// Centre camera on a bounding box and return a reasonable beam radius.
inline void centre_camera_on_bbox(OrbitalCamera &cam, const NuGeom::BoundingBox &bbox) {
    float cx = 0.5f * static_cast<float>(bbox.min.X() + bbox.max.X());
    float cy = 0.5f * static_cast<float>(bbox.min.Y() + bbox.max.Y());
    float cz = 0.5f * static_cast<float>(bbox.min.Z() + bbox.max.Z());
    float diag = static_cast<float>(std::sqrt(std::pow(bbox.max.X() - bbox.min.X(), 2) +
                                              std::pow(bbox.max.Y() - bbox.min.Y(), 2) +
                                              std::pow(bbox.max.Z() - bbox.min.Z(), 2)));
    cam = OrbitalCamera{30.f, 20.f, diag * 0.8f, cx, cy, cz};
}
