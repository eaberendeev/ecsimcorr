// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "particle_coupling.h"

#include <omp.h>

#include <iostream>
#include <source_location>

#include "algorithms_ecsim.h"
#include "algorithms_ecsim_corr.h"
#include "env_options.h"
#include "sparse.h"
#include "timer.h"

LinearEcsimCoupling::LinearEcsimCoupling() : response_row_blocks_(omp_get_max_threads()) {
    constexpr int reserve_bytes = 32 * 1024 * 1024;
    for (SimpleArrayBuffer<RowBlock<36>>& row_blocks : response_row_blocks_)
        row_blocks.reserve(reserve_bytes / sizeof(RowBlock<36>));
}

void LinearEcsimCoupling::assemble_response(Species& species, const Field3d& B_full, const Domain& domain, double dt,
                                            Operator& L) {
    Operator tmpMat;

    if (assembly_check_counter_ % envOptions::validationPeriodicity() == 0) {
        timer::commonTimer timerCopy("copy matrix");
        tmpMat = L;
    }
    timer::commonTimer timerTestAssemble("new optimized assemble");
    ecsim_response::clear_row_blocks(response_row_blocks_);

    for (auto& kv : species) {
        ParticlesArray& sp = *kv.second;
        ecsim_response::accumulate_species<SHAPE>(sp, B_full, domain, dt, response_row_blocks_);
    }
    ecsim_response::assemble_matrix(L, response_row_blocks_, response_workspace_);
    timerTestAssemble.finish();

    if (assembly_check_counter_ % envOptions::validationPeriodicity() == 0) {
        timer::commonTimer timerRefAssemble("old assemble");
        prepare_reference_blocks(species, domain);

        for (auto& kv : species) {
            ParticlesArray& sp = *kv.second;
            ecsim_response::accumulate_species_reference<SHAPE>(sp, B_full, domain, dt, reference_response_blocks_);
        }
        ecsim_response::assemble_reference(tmpMat, domain, reference_response_blocks_);
        timerRefAssemble.finish();

        const bool isSameShape = checkMatrixPortraitCoincidence(L, tmpMat);
        if (!isSameShape) {
            const std::source_location location = std::source_location::current();
            std::cerr << location.file_name() << ":" << location.line()
                      << " Warning: optimized and reference assembly produced different matrix portraits" << std::endl;
        }
        timer::commonTimer timerNorms("compute norms");
        const double diffNorm = (L - tmpMat).norm();
        const double refNorm = tmpMat.norm();
        timerNorms.finish();
        const double normalizedErr = refNorm == 0.0 && diffNorm == 0.0 ? 0.0 : diffNorm / refNorm;
        // Optimized/reference sums differ by round-off. This warning measures
        // Frobenius-relative deviation without changing the solver tolerance.
        constexpr double expectedTol = 1e-14;
        if (normalizedErr >= expectedTol) {
            const std::source_location location = std::source_location::current();
            std::cerr << location.file_name() << ":" << location.line() << " Warning: assembly deviation "
                      << normalizedErr << " exceeds expectedTol=" << expectedTol << std::endl;
        }
    }
    assembly_check_counter_ += 1;
}

void LinearEcsimCoupling::prepare_reference_blocks(const Species& species, const Domain& domain) {
    RECORD_TIMER;
    Array3D<int> countInCell(domain.size());
    countInCell.setZero();
    for (const auto& kv : species) {
        const ParticlesArray& sp = *kv.second;
#pragma omp parallel for schedule(static, 32)
        for (int i = 0; i < sp.cells().capacity(); i++) {
            countInCell(i) += sp(i).size();
        }
    }

    static_assert(SHAPE == ShapeType::Linear, "LinearEcsimCoupling supports Linear response only");
    reference_response_blocks_.resize(domain.total_size());
    reference_response_blocks_.prepare(countInCell);
}

void LinearEcsimCoupling::predict_current(const ParticlesArray& species, const Field3d& B_full, Field3d& J,
                                          double dt) const {
    algorithmsECSIM::predict_current<SHAPE>(species, B_full, J, dt);
}

void LinearEcsimCoupling::push_velocities(ParticlesArray& species, const Field3d& E_full, const Field3d& B_full,
                                          double dt) const {
    algorithmsECSIM::predict_velocity<SHAPE>(species, E_full, B_full, dt);
}

void LinearEcsimCoupling::push_and_deposit(ParticlesArray& species, const Field3d& E_full, const Field3d& B_full,
                                           double dt, double half_dt, Field3d& J, double& predictor_power,
                                           const PathPolicy& policy) const {
    ::fused_push_and_deposit<SHAPE_CH>(species, E_full, B_full, dt, half_dt, J, predictor_power, policy);
}
