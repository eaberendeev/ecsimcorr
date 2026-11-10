#pragma once


#include <stdexcept>

#include <nlohmann/json.hpp>
#include "grid/domain.h"
#include "grid/grid.h"
#include "particles/Particle.h"
#include "utils/containers.h"
#include "simulation/boundary_conditions.h"

// ECCAPFIM formulas from vakurshakov/xpic, implemented with Beren3D shapes/buffers.
// See XPIC_LICENSE and docs/eccapfim.md for provenance and time levels.
namespace algorithmsECCAPFIM {
int periodic_node(int physical_index, int axis, const Domain& domain);
int node_owner(int physical_index, int axis, const Domain& domain, const BoundaryConditionHandler* bc = nullptr);
Vector3R wrap_position(const Vector3R& position, const Domain& domain, const BoundaryConditionHandler* bc = nullptr);
void copy_periodic_nodes(Field3d& field, const Domain& domain, const BoundaryConditionHandler* bc = nullptr);

struct OrbitFields {
    Vector3R electric{0, 0, 0};
    Vector3R magnetic{0, 0, 0};
};

// Integrate along a straight orbit, split at half-integer cell coordinates.
// The electric gather and current scatter use exactly the same weights.
OrbitFields gather_orbit(const Field3d& electric, const Field3d& magnetic, const Vector3R& start, const Vector3R& end,
                         const Domain& domain, const BoundaryConditionHandler* bc = nullptr);
void deposit_orbit(Field3d& current, const Vector3R& start, const Vector3R& end, const Vector3R& velocity,
                   double charge_weight, const Domain& domain, const BoundaryConditionHandler* bc = nullptr);
// Point quadratic interpolation of B at the Yee B positions.
Vector3R gather_magnetic_point(const Field3d& magnetic, const Vector3R& position, const Domain& domain,
                               const BoundaryConditionHandler* bc = nullptr);
void deposit_charge(Field3d& density, const Vector3R& position, double charge_weight, const Domain& domain,
                    const BoundaryConditionHandler* bc = nullptr);

// Numerical failure data are assembled only on failure; kernels perform no I/O.
struct ParticleSolveFailure : std::runtime_error {
    nlohmann::json details;
    explicit ParticleSolveFailure(nlohmann::json report)
        : std::runtime_error("ECCAPFIM: particle solve did not converge: " + report.at("reason").get<std::string>()),
          details(std::move(report)) {
    }
};

struct OrbitBoundaryEvents {
    Face removed = Face::UNAVAILABLE;
    std::array<int, static_cast<size_t>(Face::Count)> reflected{};
};

// x and v are colocated at integer times. Each unwrapped substep is deposited
// before wrapping; the returned endpoint is periodic. Current uses dtau/dt.
Particle advance_particle(const Particle& initial, const Field3d& electric, const Field3d& magnetic, Field3d& current,
                          double charge_mass, double charge_weight, double dt, const Domain& domain, double tolerance,
                          int max_iterations, const BoundaryConditionHandler* bc = nullptr,
                          const PathPolicy& policy = PathPolicy::none(), OrbitBoundaryEvents* events = nullptr);
}   // namespace algorithmsECCAPFIM
