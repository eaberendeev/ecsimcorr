#include <functional>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include "ParticlesArray.h"
#include "World.h"
#include "boundary_conditions.h"
#include "containers.h"
#include "timer.h"

void Er0Condition::init_electric_field(Field3d& fieldE, const Domain& domain) const {
    RECORD_TIMER;

    if (face_ != Face::ZMIN && face_ != Face::ZMAX)
        return;
    if (width_ <= 0.0 || inner_radius_ < 0.0)
        return;

    const auto size_x = fieldE.sizes().x();
    const auto size_y = fieldE.sizes().y();
    const auto size_z = fieldE.sizes().z();
    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();

    const double center_x = 0.5 * domain.num_cells().x() * dx;
    const double center_y = 0.5 * domain.num_cells().y() * dy;

    const double outer_radius = inner_radius_ + width_;

    int k;
    if (face_ == Face::ZMIN) {
        k = domain.grid.ghost_cells();   // k=1 for ghost=1 — z=0 for Ex/Ey
    } else {
        k = size_z - domain.grid.ghost_cells() - 1;   // k = Nz - 2 for ghost=1 — z = box_max for Ex/Ey
    }

    const int ghost = domain.grid.ghost_cells();

    for (int i = 0; i < size_x; ++i) {
        for (int j = 0; j < size_y; ++j) {
            // Ex node (x-face: i+0.5, j, k)
            double xx = (i + 0.5) * dx - center_x - dx * ghost;
            double yy = j * dy - center_y - dy * ghost;
            double rr = std::hypot(xx, yy);
            if (rr >= inner_radius_ && rr <= outer_radius) {
                fieldE(i, j, k, 0) = -potential_drop_ / width_ * xx / rr;
            }

            // Ey node (y-face: i, j+0.5, k)
            xx = i * dx - center_x - dx * ghost;
            yy = (j + 0.5) * dy - center_y - dy * ghost;
            rr = std::hypot(xx, yy);
            if (rr >= inner_radius_ && rr <= outer_radius) {
                fieldE(i, j, k, 1) = -potential_drop_ / width_ * yy / rr;
            }
        }
    }
}
