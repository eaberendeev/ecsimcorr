#pragma once

#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "Particle.h"
#include "World.h"
#include "timer.h"

enum class ParticleFate { Unhandled, Removed, Reflected };

struct ParticleFateResult {
    ParticleFate fate = ParticleFate::Unhandled;
    Face face = Face::UNAVAILABLE;
};

// One geometric hit and the interaction selected at its position and incoming velocity.
struct BoundaryCrossing {
    bool clipped = false;
    Face face = Face::UNAVAILABLE;
    double t_hit = 1.0;
    Vector3R point{0, 0, 0};
    ParticleFate fate = ParticleFate::Unhandled;
};

// Sparse snapshots: only particles touching a wall consume entries. Each deposit
// thread owns its vector; lookup happens after the deposit's OpenMP barrier.
// sort_for_lookup() must be called once before the first find() of a step; it
// orders each thread vector by particle pointer so find() can use binary search.
// Do not call record() between sort_for_lookup() and the last find(): find()
// returns a pointer into a per-thread vector, which push_back may reallocate.
// The deposit classifies the exact wall hit, while the boundary handler only sees
// the trial endpoint, so the fate must be looked up by particle pointer. The
// lookup is sparse and only consulted for out-of-domain particles.
// Particle addresses stay valid until apply_to_particles, before update_cells.
class ParticleBoundaryEvents {
    struct Entry {
        const Particle* particle;
        ParticleFateResult fate;
    };
    std::vector<std::vector<Entry>> threads_;

   public:
    ParticleBoundaryEvents() : threads_(omp_get_max_threads()) {
    }
    void record(const Particle& p, ParticleFateResult fate) {
        threads_[omp_get_thread_num()].push_back({&p, fate});
    }
    void sort_for_lookup() {
        RECORD_TIMER;
        for (auto& entries : threads_)
            std::sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) {
                return std::less<const Particle*>{}(a.particle, b.particle);
            });
    }
    const ParticleFateResult* find(const Particle& p) const {
        const Particle* key = &p;
        for (const auto& entries : threads_) {
            const auto it = std::lower_bound(
                entries.begin(), entries.end(), key,
                [](const Entry& entry, const Particle* q) { return std::less<const Particle*>{}(entry.particle, q); });
            if (it != entries.end() && it->particle == key)
                return &it->fate;
        }
        return nullptr;
    }
    size_t size() const {
        size_t count = 0;
        for (const auto& entries : threads_) count += entries.size();
        return count;
    }
};

// Грани, условие на которых удаляет частицы: депозит Esirkepov-тока должен
// заканчивать траекторию ровно на такой грани (заряд уходит через стенку).
struct PathPolicy {
    // Pure classifier shared with the boundary handler. Empty for geometric-only tests.
    std::function<ParticleFateResult(const Particle&, Face)> classify;
    ParticleBoundaryEvents* events = nullptr;
    std::array<bool, static_cast<size_t>(Face::Count)> clip{};
    bool active_any = false;
    static PathPolicy none() noexcept {
        return {};
    }
    bool clips(Face f) const noexcept {
        const size_t i = static_cast<size_t>(f);
        return i < clip.size() && clip[i];
    }
    bool active() const noexcept {
        return active_any;
    }
    void enable(Face f) noexcept {
        const size_t i = static_cast<size_t>(f);
        if (i < clip.size()) {
            clip[i] = true;
            active_any = true;
        }
    }
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
inline BoundaryCrossing shape_particle_path(const Vector3R& start, const Vector3R& trial_end, const Domain& domain,
                                            const PathPolicy& policy, Vector3R& end_for_current) {
    end_for_current = trial_end;
    if (!policy.active())
        return {};

    const Vector3R& box_min = domain.geom.box_min;
    const Vector3R& box_max = domain.geom.box_max;
    const Vector3R d = trial_end - start;

    double best_t = 1.0;
    Face best_face = Face::UNAVAILABLE;

    const auto consider_plane = [&](Face min_face, Face max_face, double s, double e, double lo, double hi) {
        if (policy.clips(min_face) && s >= lo && e < lo) {
            const double t = (lo - s) / (e - s);
            if (t < best_t || best_face == Face::UNAVAILABLE) {
                best_t = t;
                best_face = min_face;
            }
        }
        if (policy.clips(max_face) && s <= hi && e >= hi && e > s) {
            const double t = (hi - s) / (e - s);
            if (t < best_t || best_face == Face::UNAVAILABLE) {
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
        return {true, best_face, best_t, end_for_current};
    }
    return {};
}

struct ParticleCurrentPath {
    Vector3R turn;
    Vector3R end;
    bool reflected = false;
    BoundaryCrossing crossing;
};

inline ParticleCurrentPath current_particle_path(const Vector3R& start, const Particle& trial, const Domain& domain,
                                                 const PathPolicy& policy) {
    ParticleCurrentPath result{trial.coord, trial.coord, false, {}};
    if (!policy.active())
        return result;
    result.crossing = shape_particle_path(start, trial.coord, domain, policy, result.end);
    result.turn = result.end;
    if (!result.crossing.clipped)
        return result;

    // Select the interaction on the first wall, using the hit position for
    // spatially selective reflection (the endpoint may be beyond another wall).
    Particle at_hit = trial;
    at_hit.coord = result.crossing.point;
    const ParticleFateResult fate = policy.classify ? policy.classify(at_hit, result.crossing.face)
                                                    : ParticleFateResult{ParticleFate::Removed, result.crossing.face};
    result.crossing.fate = fate.fate;
    if (fate.fate == ParticleFate::Reflected) {
        result.end = domain.geom.reflect_from_face(result.crossing.face, trial.coord);
        result.reflected = true;
    }
    return result;
}
