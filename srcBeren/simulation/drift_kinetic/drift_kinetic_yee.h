#pragma once

#include <stdexcept>
#include <vector>

#include "drift_kinetic_push.h"
#include "grid/domain.h"
#include "simulation/boundary_conditions.h"
#include "utils/containers.h"

// Guiding-centre coupling to the Yee grid, ported from xpic
// `drift_kinetic::DriftKineticEsirkepov` (src/algorithms/implicit_drift_kinetic.cpp)
// and `drift_kinetic::cell_segments` (src/impls/drift_kinetic/segments.cpp).
//
// Coordinates are grid units r / h. E_c lives at the node shifted by +1/2
// along c, B_c at the node shifted by +1/2 across c (grid/yee.h). Along each
// straight piece of the orbit the electric field and the current use the same
// weights: linear S1 at the piece midpoint along the component and the
// Esirkepov average of quadratic S2 at the piece ends across it. Therefore
// sum_g E_g . J_g = q V . E_p exactly, which is the discrete work of the
// field on the guiding centre. The magnetic moment is deposited as the
// magnetization M = -mu b at the B positions with the same point S2 used to
// interpolate |B| at the orbit ends.
//
// Time levels: E is the midpoint field E^{n+1/2}. B is given at three levels
// of the step, B^n, B^{n+1/2}, B^{n+1}; a sub-interval [theta0, theta1] of the
// step uses B(theta) = (1 - theta) B^n + theta B^{n+1} at its ends and at its
// midpoint for the derivatives of B.
namespace driftKinetic {

// Straight piece of an orbit in grid units, inside the cell centred on the
// node `cell` (pieces are split at half-integer planes).
struct TrackSegment {
    Vector3R r0, r1, rh, dr;
    double length = 0.0;
    Vector3I cell;
};

// Splits the orbit r0 -> r1 (grid units) at half-integer planes. Endpoints
// keep their unwrapped coordinates; node indices are wrapped by the caller.
std::vector<TrackSegment> split_at_half_planes(const Vector3R& r0, const Vector3R& r1);

// Magnetic field levels of one time step. All fields are full fields
// (self-consistent plus external).
struct MagneticLevels {
    const Field3d* begin = nullptr;   // B^n
    const Field3d* half = nullptr;    // B^{n+1/2}
    const Field3d* end = nullptr;     // B^{n+1}
};

// Time sub-interval of the step, as fractions of dt.
struct StepInterval {
    double begin = 0.0;
    double end = 1.0;
};

class YeeCoupling {
   public:
    YeeCoupling(const Domain& domain, const BoundaryConditionHandler* bc);

    // Precomputes the forward differences of B at the three levels; call once
    // per field evaluation, before `interpolate`, whenever the levels change.
    void prepare(const MagneticLevels& magnetic);

    // Fields of the guiding-centre equations for the orbit start -> end
    // (world coordinates): segment-averaged E, grad|B| and curl B, the mean
    // |B| and the normalized mean direction b of the orbit ends.
    // `start_magnetic`: B at the orbit start, fixed during the orbit iteration.
    // Uses the prepared differences when `prepare` was called for `magnetic`.
    Fields interpolate(const Field3d& electric, const MagneticLevels& magnetic, const Vector3R& start,
                       const Vector3R& end, StepInterval interval, const Vector3R* start_magnetic = nullptr) const;
    // Reference evaluation with differences taken on the fly.
    Fields interpolate_reference(const Field3d& electric, const MagneticLevels& magnetic, const Vector3R& start,
                                 const Vector3R& end, StepInterval interval) const;
    Vector3R start_magnetic(const MagneticLevels& magnetic, const Vector3R& start, StepInterval interval) const;

    // Deposits the free current q w V along the orbit (V = (end - start) /
    // tau of the sub-interval) and the magnetization -mu w b, half at each
    // orbit end. `charge_weight` = q w tau / dt, `moment_weight` = mu w tau / dt.
    // Not thread-safe: every thread deposits into its own arrays.
    void deposit(Field3d& current, Field3d& magnetization, const MagneticLevels& magnetic, const Vector3R& start,
                 const Vector3R& end, const Vector3R& velocity, double charge_weight, double moment_weight,
                 StepInterval interval) const;

    // Point S2 interpolation of a full magnetic field at the B positions.
    Vector3R magnetic_at(const Field3d& magnetic, const Vector3R& position) const;

   private:
    Vector3R magnetic_at(const MagneticLevels& magnetic, const Vector3R& position, double theta) const;
    void deposit_magnetization(Field3d& magnetization, const Vector3R& position, const Vector3R& direction,
                               double moment) const;
    Fields interpolate_impl(const Field3d& electric, const MagneticLevels& magnetic, const Vector3R& start,
                            const Vector3R& end, StepInterval interval, const Vector3R* start_magnetic,
                            bool prepared) const;

    // Storage index of the node with physical index `p` (periodic wrap).
    int owner(int p, int axis) const {
        if (periodic_[axis]) {
            const int n = cells_[axis];
            // Stencils reach at most a few nodes past a face: avoid the division.
            const int m = p < 0 ? p + n : (p >= n ? p - n : p);
            if (m >= 0 && m < n)
                return m + GHOST_CELLS;
            const int r = p % n;
            return (r < 0 ? r + n : r) + GHOST_CELLS;
        }
        const int stored = p + GHOST_CELLS;
        if (stored < 0 || stored >= storage_[axis])
            throw std::runtime_error("drift_kinetic: orbit stencil outside nonperiodic storage");
        return stored;
    }

    const Domain& domain_;
    const BoundaryConditionHandler* bc_;
    bool periodic_[3];
    int cells_[3], storage_[3];
    Field3d differences_[3];   // d_a B_c of B^n, B^{n+1/2}, B^{n+1}
    MagneticLevels prepared_;
};

}   // namespace driftKinetic
