#pragma once

#include <array>
#include <cmath>

#include "core/Vector3.h"

enum class Face { XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, CYLINDER, Count, UNAVAILABLE };
inline constexpr std::array<Face, 7> ALL_FACES{Face::XMIN, Face::XMAX, Face::YMIN,    Face::YMAX,
                                               Face::ZMIN, Face::ZMAX, Face::CYLINDER};

// -----------------------------------------------------------------------------
// Геометрия домена: прямоугольный параллелепипед + цилиндр (внутренняя область)
// -----------------------------------------------------------------------------
struct Geometry {
    Vector3R box_min, box_max;   // границы прямоугольной области
    bool use_cylinder = false;   // активен ли цилиндр
    Vector3R cyl_center;
    double cyl_radius;

    Geometry() : box_min(0, 0, 0), box_max(0, 0, 0), use_cylinder(false), cyl_center(0, 0, 0), cyl_radius(0.0) {
    }

    Geometry(const Vector3R& bmin, const Vector3R& bmax, const bool cyl = false) {
        init(bmin, bmax, cyl);
    }
    void init(const Vector3R& bmin, const Vector3R& bmax, const bool cyl = false) {
        box_min = bmin;
        box_max = bmax;
        use_cylinder = cyl;
    }
    bool in_cylinder(const Vector3R& p, double eps = 0.0) const {
        if (!use_cylinder)
            return false;   // цилиндр отсутствует ⇒ точка не внутри
        double dx = p.x() - cyl_center.x();
        double dy = p.y() - cyl_center.y();
        return (dx * dx + dy * dy <= cyl_radius * cyl_radius - eps);
    }
    // Проверка, находится ли точка внутри области (box ∩ cylinder)
    bool contains(const Vector3R& p, double eps = 0.0) const {
        // Проверка прямоугольного параллелепипеда
        if (p.x() < box_min.x() + eps || p.x() >= box_max.x() - eps)
            return false;
        if (p.y() < box_min.y() + eps || p.y() >= box_max.y() - eps)
            return false;
        if (p.z() < box_min.z() + eps || p.z() >= box_max.z() - eps)
            return false;

        // Если цилиндр не используется, точка внутри
        if (!use_cylinder)
            return true;

        // Проверка цилиндра
        if (!in_cylinder(p, eps))
            return false;

        return true;
    }

    bool contains(const Vector3R& p, bool periodic[3], double eps = 0.0) const {
        // Проверка прямоугольного параллелепипеда
        if (!periodic[0] && (p.x() < box_min.x() + eps || p.x() >= box_max.x() - eps))
            return false;
        if (!periodic[1] && (p.y() < box_min.y() + eps || p.y() >= box_max.y() - eps))
            return false;
        if (!periodic[2] && (p.z() < box_min.z() + eps || p.z() >= box_max.z() - eps))
            return false;

        // Если цилиндр не используется, точка внутри
        if (!use_cylinder)
            return true;

        // Проверка цилиндра
        if (!in_cylinder(p, eps))
            return false;

        return true;
    }

    // Функция, возвращающая true, если точка p находится вне указанной грани
    bool is_outside_face(const Face face, const Vector3R& p, double eps = 0.0) const {
        switch (face) {
            case Face::XMIN:
                return p.x() < box_min.x() + eps;
            case Face::XMAX:
                return p.x() >= box_max.x() - eps;
            case Face::YMIN:
                return p.y() < box_min.y() + eps;
            case Face::YMAX:
                return p.y() >= box_max.y() - eps;
            case Face::ZMIN:
                return p.z() < box_min.z() + eps;
            case Face::ZMAX:
                return p.z() >= box_max.z() - eps;
            case Face::CYLINDER:
                if (!use_cylinder) {
                    return false;
                }
                return !in_cylinder(p, eps);
            default:
                return false;   // на случай добавления новых граней
        }
    }
    // Проверяет, находится ли точка внутри области, но игнорируя указанную
    // грань
    // (как будто эта грань не ограничивает область)
    bool contains_ignoring_face(Face face, const Vector3R& p, double eps = 0.0) const {
        // Прямоугольные грани – проверяем только те, которые не совпадают с
        // face
        if (face != Face::XMIN && p.x() < box_min.x() + eps)
            return false;
        if (face != Face::XMAX && p.x() >= box_max.x() - eps)
            return false;
        if (face != Face::YMIN && p.y() < box_min.y() + eps)
            return false;
        if (face != Face::YMAX && p.y() >= box_max.y() - eps)
            return false;
        if (face != Face::ZMIN && p.z() < box_min.z() + eps)
            return false;
        if (face != Face::ZMAX && p.z() >= box_max.z() - eps)
            return false;

        // Цилиндр – проверяем, только если цилиндр активен и не является
        // игнорируемой гранью
        if (face != Face::CYLINDER && use_cylinder) {
            if (!in_cylinder(p, eps))
                return false;
        }
        return true;
    }
    // Точка строго вне указанной грани, но внутри всех остальных граней
    // области.
    // Возвращает true, если одновременно:
    //   - is_outside_face(face, p) == true
    //   - contains_ignoring_face(face, p) == true
    bool is_outside_only_face(Face face, const Vector3R& p, double eps = 0.0) const {
        return is_outside_face(face, p, eps) && contains_ignoring_face(face, p, eps);
    }
    // Отражает точку p симметрично относительно заданной грани.
    // Предполагается, что точка находится вне этой грани.
    Vector3R reflect_from_face(Face face, const Vector3R& p) const {
        Vector3R reflected = p;
        switch (face) {
            case Face::XMIN:
                reflected.x() = 2.0 * box_min.x() - p.x();
                break;
            case Face::XMAX:
                reflected.x() = 2.0 * box_max.x() - p.x();
                break;
            case Face::YMIN:
                reflected.y() = 2.0 * box_min.y() - p.y();
                break;
            case Face::YMAX:
                reflected.y() = 2.0 * box_max.y() - p.y();
                break;
            case Face::ZMIN:
                reflected.z() = 2.0 * box_min.z() - p.z();
                break;
            case Face::ZMAX:
                reflected.z() = 2.0 * box_max.z() - p.z();
                break;
            case Face::CYLINDER: {
                double dx = p.x() - cyl_center.x();
                double dy = p.y() - cyl_center.y();
                double dist = std::sqrt(dx * dx + dy * dy);
                if (dist > 0.0) {
                    double new_dist = 2.0 * cyl_radius - dist;
                    double scale = new_dist / dist;
                    reflected.x() = cyl_center.x() + dx * scale;
                    reflected.y() = cyl_center.y() + dy * scale;
                }
                break;
            }
            default:
                break;
        }
        return reflected;
    }

    // Простое отражение точки от границы области.
    // Последовательно проверяет все грани и отражает точку от тех,
    // за которыми она оказалась. Результат гарантированно лежит внутри
    // области (численная погрешность ~ 1e-15).
    Vector3R reflect_from_boundary(const Vector3R& p) const {
        Vector3R result = p;
        bool changed = true;
        int max_iter = 10;
        while (changed && max_iter-- > 0) {
            changed = false;
            for (auto face : ALL_FACES) {
                if (is_outside_face(face, result, 0.0)) {
                    result = reflect_from_face(face, result);
                    changed = true;
                }
            }
        }
        return result;
    }
};
