#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>

#include "Particle.h"
#include "World.h"

enum class ParticleFate { Unhandled, Removed, Reflected };

struct ParticleFateResult {
    ParticleFate fate = ParticleFate::Unhandled;
    Face face = Face::UNAVAILABLE;
};

// Грани, условие на которых удаляет частицы: депозит Esirkepov-тока должен
// заканчивать траекторию ровно на такой грани (заряд уходит через стенку).
struct PathPolicy {
    bool can_reflect = false;
    // Pure classifier shared with the boundary handler. Empty for geometric-only tests.
    std::function<ParticleFateResult(const Particle&)> classify;
    std::array<bool, static_cast<size_t>(Face::Count)> clip{};
    bool active_any = false;
    static PathPolicy none() noexcept {
        return {};
    }
    bool clips(Face f) const noexcept {
        const auto i = static_cast<size_t>(f);
        return i < clip.size() && clip[i];
    }
    bool active() const noexcept {
        return active_any;
    }
    void enable(Face f) noexcept {
        const auto i = static_cast<size_t>(f);
        if (i < clip.size()) {
            clip[i] = true;
            active_any = true;
        }
    }
};

struct PathResult {
    bool clipped = false;
    Face face = Face::UNAVAILABLE;
};

// Привязывает пересечённую плоскую грань точно к координате стенки.
inline void snap_path_end_to_plane(Face face, const Domain& domain, Vector3R& p) noexcept {
    switch (face) {
        case Face::XMIN:
            p.x() = domain.geom.box_min.x();
            break;
        case Face::XMAX:
            p.x() = domain.geom.box_max.x();
            break;
        case Face::YMIN:
            p.y() = domain.geom.box_min.y();
            break;
        case Face::YMAX:
            p.y() = domain.geom.box_max.y();
            break;
        case Face::ZMIN:
            p.z() = domain.geom.box_min.z();
            break;
        case Face::ZMAX:
            p.z() = domain.geom.box_max.z();
            break;
        default:
            break;
    }
}

// Обрезает отрезок [start, trial_end] в ПЕРВОЙ удаляющей грани, через которую
// он выходит. end_for_current получает trial_end, если обрезки не было.
inline PathResult shape_particle_path(const Vector3R& start, const Vector3R& trial_end, const Domain& domain,
                                      const PathPolicy& policy, Vector3R& end_for_current) {
    end_for_current = trial_end;
    if (!policy.active())
        return {};

    const auto& box_min = domain.geom.box_min;
    const auto& box_max = domain.geom.box_max;
    const Vector3R d = trial_end - start;

    double best_t = 1.0;
    Face best_face = Face::UNAVAILABLE;

    const auto consider_plane = [&](Face min_face, Face max_face, double s, double e, double lo, double hi) {
        if (policy.clips(min_face) && s >= lo && e < lo) {
            const double t = (lo - s) / (e - s);
            if (t < best_t) {
                best_t = t;
                best_face = min_face;
            }
        }
        if (policy.clips(max_face) && s < hi && e >= hi) {
            const double t = (hi - s) / (e - s);
            if (t < best_t) {
                best_t = t;
                best_face = max_face;
            }
        }
    };

    consider_plane(Face::XMIN, Face::XMAX, start.x(), trial_end.x(), box_min.x(), box_max.x());
    consider_plane(Face::YMIN, Face::YMAX, start.y(), trial_end.y(), box_min.y(), box_max.y());
    consider_plane(Face::ZMIN, Face::ZMAX, start.z(), trial_end.z(), box_min.z(), box_max.z());

    // Цилиндр: |(p.xy(t) - c.xy)|^2 = R^2.
    if (policy.clips(Face::CYLINDER) && domain.geom.use_cylinder) {
        const Vector3R& c = domain.geom.cyl_center;
        const double R = domain.geom.cyl_radius;
        const double sx = start.x();
        const double sy = start.y();
        const double dx = d.x();
        const double dy = d.y();
        const double rs2 = (sx - c.x()) * (sx - c.x()) + (sy - c.y()) * (sy - c.y());
        // Старт уже снаружи цилиндра — обрезать нечего (защитный случай).
        if (rs2 <= R * R) {
            const double a = dx * dx + dy * dy;
            if (a >= 1e-30) {
                const double b = 2.0 * ((sx - c.x()) * dx + (sy - c.y()) * dy);
                const double cc = (sx - c.x()) * (sx - c.x()) + (sy - c.y()) * (sy - c.y()) - R * R;
                const double discr = b * b - 4.0 * a * cc;
                if (discr >= 0.0) {
                    const double sq = std::sqrt(discr);
                    const double t1 = (-b - sq) / (2.0 * a);
                    const double t2 = (-b + sq) / (2.0 * a);
                    // Берём наименьший корень в [0,1] с ПОЛОЖИТЕЛЬНОЙ
                    // производной dr2/dt: там отрезок выходит из цилиндра
                    // (касание стенки изнутри при t=0 не обрезается).
                    const auto try_root = [&](double t) {
                        if (t < 0.0 || t > 1.0)
                            return;
                        const double dr2dt = b + 2.0 * a * t;
                        if (dr2dt > 0.0 && t < best_t) {
                            best_t = t;
                            best_face = Face::CYLINDER;
                        }
                    };
                    try_root(t1);
                    try_root(t2);
                }
            }
        }
    }

    if (best_face != Face::UNAVAILABLE) {
        end_for_current = start + best_t * (trial_end - start);
        snap_path_end_to_plane(best_face, domain, end_for_current);
        return {true, best_face};
    }
    return {};
}

struct ParticleCurrentPath {
    Vector3R turn;
    Vector3R end;
    bool reflected = false;
};

inline ParticleCurrentPath current_particle_path(const Vector3R& start, const Particle& trial,
                                                 const Domain& domain, const PathPolicy& policy) {
    ParticleCurrentPath result{trial.coord, trial.coord, false};
    if (!policy.active())
        return result;
    const auto crossing = shape_particle_path(start, trial.coord, domain, policy, result.end);
    result.turn = result.end;
    if (!crossing.clipped)
        return result;
    const auto fate = policy.classify ? policy.classify(trial) : ParticleFateResult{};
    if (fate.fate == ParticleFate::Reflected) {
        PathPolicy wall;
        wall.enable(fate.face);
        const auto crossing = shape_particle_path(start, trial.coord, domain, wall, result.turn);
        if (crossing.clipped) {
            result.end = domain.geom.reflect_from_face(fate.face, trial.coord);
            result.reflected = true;
        }
    }
    return result;
}
