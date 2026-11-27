#pragma once

#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

#include "Particle.h"
#include "grid/domain.h"
#include "runtime/profiler.h"

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
    std::array<bool, 3> periodic{};

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

// Конец тока поглощённой частицы: на полъячейки за плоской стенкой (по нормали,
// в точке удара). На самой стенке квадратичное ядро оставляет 1/8 заряда на первом
// внутреннем узле; этот заряд остаётся в поле фантомом (div E != rho у стенки).
// С полъячейки за стенкой вес первого внутреннего узла нулевой (для линейного
// ядра тоже), весь заряд уходит в узлы стенки и за ней, где его снимает граничное
// условие. Цилиндр: cylinder_absorbed_path_end.
inline void extend_absorbed_path_end(Face face, const Domain& domain, Vector3R& end) noexcept {
    const Vector3R half_cell = 0.5 * domain.cell_size();
    switch (face) {
        case Face::XMIN:
            end.x() -= half_cell.x();
            break;
        case Face::XMAX:
            end.x() += half_cell.x();
            break;
        case Face::YMIN:
            end.y() -= half_cell.y();
            break;
        case Face::YMAX:
            end.y() += half_cell.y();
            break;
        case Face::ZMIN:
            end.z() -= half_cell.z();
            break;
        case Face::ZMAX:
            end.z() += half_cell.z();
            break;
        default:
            break;
    }
}

// Узел решается полем изнутри цилиндра: все шесть рёбер E при нём лежат внутри
// (как маска открытого цилиндра с gap = 0). На узле с ребром снаружи закон Гаусса
// не держится, его заряд снимает стенка.
inline bool node_inside_cylinder_field(const Domain& domain, int i, int j, int k) noexcept {
    constexpr double mask_eps = 1e-10;   // eps открытого условия
    const std::array<Vector3I, 5> edges{Vector3I(i, j, k), Vector3I(i - 1, j, k), Vector3I(i, j, k),
                                        Vector3I(i, j - 1, k), Vector3I(i, j, k)};
    const std::array<int, 5> components{0, 0, 1, 1, 2};
    for (int e = 0; e < 5; ++e) {
        const Vector3R position = domain.grid.node_position(edges[e], FieldType::ELECTRIC, components[e]);
        if (domain.geom.is_outside_face(Face::CYLINDER, position, mask_eps))
            return false;
    }
    return true;
}

// Конец тока поглощённой на цилиндре частицы: ближайшая к точке удара точка на
// радиальном луче наружу, где форма радиуса support_radius (в ячейках) не задевает
// ни одного узла, решаемого полем изнутри. Иначе часть заряда остаётся на таком
// узле фантомом (div E != rho), как на плоской стенке без сдвига. Узлы стенки
// лежат на сетке ступенькой, поэтому сдвиг не постоянный: для каждого
// внутреннего узла рядом берётся отрезок луча, где форма его задевает, и конец
// ставится в первую точку луча вне всех отрезков (до ~(r + 0.5)·sqrt(2) ячеек).
inline Vector3R cylinder_absorbed_path_end(const Domain& domain, const Vector3R& hit, double support_radius) {
    const Vector3R& center = domain.geom.cyl_center;
    const double rx = hit.x() - center.x();
    const double ry = hit.y() - center.y();
    const double r = std::hypot(rx, ry);
    if (r == 0.0)
        return hit;
    const Vector3R c = domain.to_cell_coordinates(hit);
    // Луч в ячейках: c + s * u, s — длина в мировых единицах.
    const Vector3R u(rx / r / domain.cell_size().x(), ry / r / domain.cell_size().y(), 0.0);
    const int k = static_cast<int>(std::floor(c.z()));
    constexpr int window = 6;   // больше максимального сдвига плюс радиус формы
    const int i0 = static_cast<int>(std::floor(c.x()));
    const int j0 = static_cast<int>(std::floor(c.y()));

    std::array<double, (2 * window + 1) * (2 * window + 1)> lo{};
    std::array<double, (2 * window + 1) * (2 * window + 1)> hi{};
    int count = 0;
    for (int i = i0 - window; i <= i0 + window; ++i)
        for (int j = j0 - window; j <= j0 + window; ++j) {
            if (!node_inside_cylinder_field(domain, i, j, k))
                continue;
            double a = -std::numeric_limits<double>::infinity();
            double b = std::numeric_limits<double>::infinity();
            const std::array<double, 2> offset{i - c.x(), j - c.y()};
            const std::array<double, 2> direction{u.x(), u.y()};
            for (int d = 0; d < 2; ++d) {
                if (direction[d] == 0.0) {
                    if (std::abs(offset[d]) >= support_radius)
                        b = a;
                    continue;
                }
                const double t1 = (offset[d] - support_radius) / direction[d];
                const double t2 = (offset[d] + support_radius) / direction[d];
                a = std::max(a, std::min(t1, t2));
                b = std::min(b, std::max(t1, t2));
            }
            if (b > 0.0 && b > a) {
                lo[count] = a;
                hi[count] = b;
                ++count;
            }
        }

    // Отрезки открытые: на их конце вес узла уже нулевой.
    double s = 0.0;
    bool moved = true;
    while (moved) {
        moved = false;
        for (int n = 0; n < count; ++n)
            if (lo[n] < s && s < hi[n]) {
                s = hi[n];
                moved = true;
            }
    }
    return Vector3R(hit.x() + s * rx / r, hit.y() + s * ry / r, hit.z());
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

// Current path of one particle: start -> turn, then turn -> end for a reflected
// leg or for the absorption tail beyond the cylinder (tail: turn is the hit).
// reflected and tail are mutually exclusive: a particle is either reflected
// from the wall or absorbed by the cylinder.
struct ParticleCurrentPath {
    Vector3R turn;
    Vector3R end;
    BoundaryCrossing crossing;
    bool reflected = false;
    bool tail = false;
};

// support_radius: half-width of the deposit shape in cells (1 linear, 1.5 quadratic).
inline ParticleCurrentPath current_particle_path(const Vector3R& start, const Particle& trial, const Domain& domain,
                                                 const PathPolicy& policy, double support_radius) {
    ParticleCurrentPath result{trial.coord, trial.coord, {}, false, false};
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
    } else if (fate.fate == ParticleFate::Removed && result.crossing.face == Face::CYLINDER) {
        result.end = cylinder_absorbed_path_end(domain, result.turn, support_radius);
        result.tail = true;
    } else if (fate.fate == ParticleFate::Removed) {
        extend_absorbed_path_end(result.crossing.face, domain, result.end);
        result.turn = result.end;
    }
    return result;
}
