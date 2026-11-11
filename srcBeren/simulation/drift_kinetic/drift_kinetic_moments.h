#pragma once

#include <map>
#include <memory>
#include <string>

#include "diagnostics/SpeciesMomentSource.h"
#include "drift_kinetic_particles.h"
#include "drift_kinetic_yee.h"

namespace driftKinetic {

// Grid moments of one guiding-centre species for the diagnostic outputs, at
// the accepted step (positions r^{n+1}, v_par^{n+1}, B^{n+1}).
// - density: guiding-centre charge density q w S(r); it differs from the
//   particle density by the polarization and magnetization corrections of
//   order (rho / L)^2;
// - current: free guiding-centre current plus the magnetization current,
//   J_gc + curl_B M, on the E positions; the sum over species is the current
//   that enters Maxwell's equations;
// - second moment: the gyrotropic tensor m w [v_par^2 b b + (mu |B| / m) (I - b b)],
//   i.e. p_par b b + p_perp (I - b b) plus the parallel flow; the drift flow
//   is omitted, like the drift energy in the guiding-centre energy.
class GuidingCenterMoments : public SpeciesMomentSource {
   public:
    using Species = std::map<std::string, std::unique_ptr<DriftKineticParticles>>;

    // The references must outlive the adapter; `species` is looked up by name
    // on every call, since the simulation replaces the containers each step.
    GuidingCenterMoments(std::string name, const Species& species, const YeeCoupling& coupling, const Field3d& magnetic,
                         const Field3d& external_magnetic, const Field3d& current)
        : name_(std::move(name)),
          species_(species),
          coupling_(coupling),
          magnetic_(magnetic),
          external_magnetic_(external_magnetic),
          current_(current) {
    }

    bool is_neutral() const override {
        return particles().is_neutral();
    }

    void deposit_density(Field3d& density, ShapeType shape) const override;

    const Field3d& current() const override {
        return current_;
    }

    void deposit_second_moment(Field3d& moment, int a, int b) const override;

   private:
    const DriftKineticParticles& particles() const {
        return *species_.at(name_);
    }

    std::string name_;
    const Species& species_;
    const YeeCoupling& coupling_;
    const Field3d& magnetic_;            // self-consistent B^{n+1}
    const Field3d& external_magnetic_;   // prescribed B
    const Field3d& current_;
};

}   // namespace driftKinetic
