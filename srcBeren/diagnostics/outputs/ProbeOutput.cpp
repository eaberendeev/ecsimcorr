#include "ProbeOutput.h"

#include <iomanip>
#include <stdexcept>

#include "interpolation.h"
#include "logger.h"
#include "util.h"

// ============================================================
//  ProbeOutput  — interpolate fields at probe points, write probes.txt
// ============================================================
namespace {

double interpolate_scalar_cell(const Field3d& field, int ix, int iy, int iz, double wx, double wy, double wz) {
    const double w000 = (1 - wx) * (1 - wy) * (1 - wz);
    const double w001 = (1 - wx) * (1 - wy) * wz;
    const double w010 = (1 - wx) * wy * (1 - wz);
    const double w011 = (1 - wx) * wy * wz;
    const double w100 = wx * (1 - wy) * (1 - wz);
    const double w101 = wx * (1 - wy) * wz;
    const double w110 = wx * wy * (1 - wz);
    const double w111 = wx * wy * wz;
    const Field3d& f = field;
    return w000 * f(ix, iy, iz, 0) + w001 * f(ix, iy, iz + 1, 0) + w010 * f(ix, iy + 1, iz, 0) +
           w011 * f(ix, iy + 1, iz + 1, 0) + w100 * f(ix + 1, iy, iz, 0) + w101 * f(ix + 1, iy, iz + 1, 0) +
           w110 * f(ix + 1, iy + 1, iz, 0) + w111 * f(ix + 1, iy + 1, iz + 1, 0);
}

Vector3R interpolate_vector_at(const Field3d& field, const Vector3R& world, const Domain& domain, bool is_b) {
    const Vector3R cell = domain.cell_size();
    if (is_b) {
        return interpolateB_linear(field, Vector3R(world.x() / cell.x(), world.y() / cell.y(), world.z() / cell.z()));
    } else {
        return interpolateE_linear(field, Vector3R(world.x() / cell.x(), world.y() / cell.y(), world.z() / cell.z()));
    }
}

double interpolate_scalar_at(const Field3d& field, const Vector3R& world, const Domain& domain) {
    const Vector3R cell = domain.cell_size();
    const double x = world.x() / cell.x() + GHOST_CELLS;
    const double y = world.y() / cell.y() + GHOST_CELLS;
    const double z = world.z() / cell.z() + GHOST_CELLS;
    const int ix = int(x), iy = int(y), iz = int(z);
    return interpolate_scalar_cell(field, ix, iy, iz, x - ix, y - iy, z - iz);
}

}   // namespace

ProbeOutput::ProbeOutput(const Domain& domain, const Field3d& fieldEn, const Field3d& fieldBn,
                         const Field3d& fieldE_external, const Field3d& fieldBInit, const Species& species,
                         const std::vector<std::string>& species_names, const std::vector<Vector3R>& points,
                         int interval)
    : domain_(domain),
      fieldEn_(fieldEn),
      fieldBn_(fieldBn),
      fieldE_external_(fieldE_external),
      fieldBInit_(fieldBInit),
      species_(species),
      species_names_(species_names),
      points_(points),
      interval_(interval),
      header_written_(false) {
    ensure_directory(".//Probes");
    fProbes_.open(".//Probes//probes.txt");
    if (!fProbes_.is_open()) {
        const std::string message = "Cannot open .//Probes//probes.txt for writing";
        logger::error(message);
        throw std::runtime_error(message);
    }
}

void ProbeOutput::output(int timestep, const Diagnostics&) {
    if (timestep % interval_ != 0)
        return;
    if (points_.empty())
        return;

    if (!header_written_) {
        fProbes_ << "# Time x y z";
        fProbes_ << " Ex Ey Ez Bx By Bz";
        fProbes_ << " Enx Eny Enz Bnx Bny Bnz";
        fProbes_ << " E_ex_x E_ex_y E_ex_z B_init_x B_init_y B_init_z";
        for (const auto& name : species_names_) fProbes_ << " dens_" << name;
        fProbes_ << "\n";
        header_written_ = true;
    }

    auto field_e_full = [&]() {
        Field3d f(fieldEn_.sizes(), fieldEn_.nd());
        f.data() = fieldEn_.data() + fieldE_external_.data();
        return f;
    };
    auto field_b_full = [&]() {
        Field3d f(fieldBn_.sizes(), fieldBn_.nd());
        f.data() = fieldBn_.data() + fieldBInit_.data();
        return f;
    };

    for (const auto& pt : points_) {
        fProbes_ << timestep << " " << pt.x() << " " << pt.y() << " " << pt.z();

        const Field3d eTmp = field_e_full();
        const Field3d bTmp = field_b_full();
        fProbes_ << std::scientific << std::setprecision(6);
        Vector3R w = interpolate_vector_at(eTmp, pt, domain_, false);
        fProbes_ << " " << w.x() << " " << w.y() << " " << w.z();
        w = interpolate_vector_at(bTmp, pt, domain_, true);
        fProbes_ << " " << w.x() << " " << w.y() << " " << w.z();
        w = interpolate_vector_at(fieldEn_, pt, domain_, false);
        fProbes_ << " " << w.x() << " " << w.y() << " " << w.z();
        {
            Field3d bOnly(fieldBn_.sizes(), fieldBn_.nd());
            bOnly.data() = fieldBn_.data();
            w = interpolate_vector_at(bOnly, pt, domain_, true);
        }
        fProbes_ << " " << w.x() << " " << w.y() << " " << w.z();
        w = interpolate_vector_at(fieldE_external_, pt, domain_, false);
        fProbes_ << " " << w.x() << " " << w.y() << " " << w.z();
        w = interpolate_vector_at(fieldBInit_, pt, domain_, true);
        fProbes_ << " " << w.x() << " " << w.y() << " " << w.z();

        for (const auto& name : species_names_) {
            const auto it = species_.find(name);
            double d = 0;
            if (it != species_.end()) {
                d = interpolate_scalar_at(it->second->densityOnGrid, pt, domain_);
            }
            fProbes_ << " " << d;
        }
        fProbes_ << "\n";
    }
    fProbes_.flush();
}

void ProbeOutput::finalize() {
    if (fProbes_.is_open()) {
        fProbes_.flush();
        fProbes_.close();
    }
}
