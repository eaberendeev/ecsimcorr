// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "simulation_ecsim.h"

#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "Coil.h"
#include "Damping.h"
#include "Diagnostic.h"
#include "DiagnosticOutput.h"
#include "Mesh.h"
#include "ParticlesArray.h"
#include "ParticlesDiagnostic.h"
#include "Read.h"
#include "World.h"
#include "collision.h"
#include "containers.h"
#include "external_fieldsB.h"
#include "external_fieldsE.h"
#include "log_macros.h"
#include "mesh/aux.h"
#include "recovery.h"
#include "solverSLE.h"
#include "timer.h"

void SimulationEcsim::first_push() {
    RECORD_TIMER;

    const double dt = get_checked<double>(system_config, "Dt");

    globalTimer.start("particles1");
    // parallel allocation-free fieldBFull = fieldB + fieldBInit;
    blas::sum(1.0, fieldB.data(), 1.0, fieldBInit.data(), fieldBFull.data());

    for (auto &kv : species) {
        auto &sp = *kv.second;
        //  +++ x_{n-1/2} -> x_{n+1/2}
        sp.move(dt);

        sp.update_cells(domain);
        bc_handler.apply_to_particles(sp, species, domain);

        // +++ get J(x_{n+1/2},v_n)_predict
        algorithmsECSIM::predict_current(sp, fieldBFull, fieldJp, dt, SHAPE);
    }
    globalTimer.finish("particles1");

    globalTimer.start("particlesLmat2");

    // for (auto &sp : species) {
    //     sp->fill_matrixL(mesh, fieldBFull, domain, dt, SHAPE);
    // }

    prepare_block_matrix(SHAPE);

    std::vector<std::vector<RowBlock<12>>> rowBlocksGlobal;
    rowBlocksGlobal.resize(omp_get_max_threads());

    for (auto &kv : species) {
        ParticlesArray &sp = *kv.second;
        sp.fill_matrixL2(mesh, fieldBFull, domain, dt, SHAPE);
    }

    for (auto &kv : species) {
        ParticlesArray &sp = *kv.second;
        sp.fill_matrixL2_Optimized(mesh, fieldBFull, domain, dt, SHAPE, rowBlocksGlobal);
    }
    // todo: zeros Lmat + current
    globalTimer.finish("particlesLmat2");

    bc_handler.apply_to_fields(fieldJp, FieldType::CURRENT, domain);

    globalTimer.start("bound1");
    // mesh.apply_boundaries(mesh.LmatX, domain);
    globalTimer.finish("bound1");

    globalTimer.start("stencilLmat2");

    Operator tmpMat = mesh.Lmat2;
    mesh.stencil_Lmat2(mesh.Lmat2, domain, mesh.workspacePtr);
    mesh.stencil_Lmat2_Optimized_V2(tmpMat, domain, rowBlocksGlobal, mesh.workspacePtr);

    checkMatrixCoincidence(mesh.Lmat2, tmpMat, 1e-10);

    // convert_block_matrix(SHAPE);

    globalTimer.finish("stencilLmat2");

    globalTimer.start("bound2");
    bc_handler.apply_to_operator(mesh.Lmat2, domain);
    globalTimer.finish("bound2");
}

void SimulationEcsim::second_push() {
    RECORD_TIMER;

    const double dt = get_checked<double>(system_config, "Dt");

    globalTimer.start("particles2");

    blas::sum(1.0, fieldB.data(), 1.0, fieldBInit.data(), fieldBFull.data());
    Field3d fieldE_full = fieldEp + fieldE_external;
    for (auto &kv : species) {
        auto &sp = *kv.second;
        // +++ get v'_{n+1} from v_{n} and E'_{n+1/2}
        algorithmsECSIM::predict_velocity(sp, fieldE_full, fieldBFull, dt, SHAPE);
    }
    globalTimer.finish("particles2");
}

// Particles have ccordinates and velocities. Mesh have 3D fields in nodes (each
// field stored in 1D array with 4d index x,y,z,d)
void SimulationEcsim::make_step([[maybe_unused]] const int timestep) {
    RECORD_TIMER;
    globalTimer.start("Total");

    first_push();

    globalTimer.start("FieldsPredict");
    // --- solve A*E'_{n+1/2}=f(E_n, B_n, J(x_{n+1/2})).
    predict_electric_field(fieldEp, fieldE, fieldE_external, fieldB, fieldJp);
    bc_handler.apply_to_fields(fieldEp, FieldType::ELECTRIC, domain);

    globalTimer.finish("FieldsPredict");

    second_push();

    for (auto &kv : species) {
        auto &sp = *kv.second;
        sp.density_on_grid_update(SHAPE);
        bc_handler.apply_to_fields(sp.densityOnGrid, FieldType::DENSITY, domain);
    }

    globalTimer.start("computeB");
    // calculate fieldB
    blas::sum(2.0, fieldEp.data(), -1.0, fieldE.data(), fieldEn.data());
    bc_handler.apply_to_fields(fieldEn, FieldType::ELECTRIC, domain);

    mesh.compute_fieldB(fieldBn, fieldB, fieldE, fieldEn, get_checked<double>(system_config, "Dt"));
    bc_handler.apply_to_fields(fieldBn, FieldType::MAGNETIC, domain);
    globalTimer.finish("computeB");

    globalTimer.finish("Total");
}

void SimulationEcsim::prepare_block_matrix(ShapeType type) {
    RECORD_TIMER;
    Array3D<int> countInCell(domain.size());
    countInCell.setZero();
    for (auto &kv : species) {
        auto &sp = *kv.second;
#pragma omp parallel for schedule(static, 32)
        for (int i = 0; i < sp.particlesData.capacity(); i++) {
            countInCell(i) += sp.particlesData(i).size();
        }
    }

    switch (type) {
        case ShapeType::NGP:
            mesh.LmatX_NGP.prepare(countInCell);
            break;
        case ShapeType::Linear:
            mesh.LmatX2.prepare(countInCell);
            break;
        case ShapeType::Quadratic:
            std::cout << "Fill Lmatrix for quadratic shape function is not "
                         "implemented"
                      << std::endl;
            exit(-1);
    }
}

void SimulationEcsim::convert_block_matrix(ShapeType type) {
    switch (type) {
        case ShapeType::NGP:
            mesh.convert_block_to_crs_format<XIndexerNGP, YIndexerNGP, ZIndexerNGP>(mesh.LmatX_NGP, mesh.Lmat2, domain);
            break;
        case ShapeType::Linear:
            mesh.convert_block_to_crs_format<XIndexer, YIndexer, ZIndexer>(mesh.LmatX2, mesh.Lmat2, domain);
            break;
        case ShapeType::Quadratic:
            std::cout << "Fill Lmatrix for quadratic shape function is not "
                         "implemented"
                      << std::endl;
            exit(-1);
    }
}

void SimulationEcsim::predict_electric_field(Field3d &Ep, const Field3d &E, const Field3d &B, Field3d &J) {
    const double dt = get_checked<double>(system_config, "Dt");

    Operator A = mesh.IMmat + mesh.Lmat2;
    mesh.Lmat2.makeCompressed();

    Field3d rhs = E - 0.5 * dt * J + 0.5 * dt * mesh.curlB * B;

    // E(n+1/2) = (M-L) * E(n+1/2) + E - 0.5*dt*(J + rotB)
    solve_linear_system<BicgstabSolver<Field3d>>(A, rhs, Ep, E);
    LOG_STEP("  solver error=" << (A * Ep - rhs).norm() << "\n");
}

void SimulationEcsim::predict_electric_field(Field3d &Ep, const Field3d &E, const Field3d &E_ex, const Field3d &B,
                                             Field3d &J) {
    RECORD_TIMER;

    const double dt = get_checked<double>(system_config, "Dt");

    timer::flatTimer timerDestructors(timer::NoStart{});
    {
        timer::commonTimer timerA("construct A");
        Operator A = parallelSparseSum(mesh.IMmat, mesh.Lmat2);
        timerA.finish();

        timer::commonTimer compressTimer("compress Lmat2");
        mesh.Lmat2.makeCompressed();
        compressTimer.finish();

        timer::commonTimer timerRhs("make rhs");
        Field3d rhs = E + 0.5 * dt * (mesh.curlB * B - J) - mesh.Lmat2 * E_ex;
        timerRhs.finish();

        // E(n+1/2) = (M-L) * E(n+1/2)  - L*E_ex + E - 0.5*dt*(J + rotB)
        // (M*Ex = 0)
        const double err = solve_linear_system<BicgstabSolver<Field3d>>(A, rhs, Ep, E);
        LOG_STEP("  solver error=" << err << "\n");

        // A и rhs уничтожаются при выходе из этого scope — замеряем их деструкторы
        timerDestructors.start("destructor operator A and rhs");
    }
}

void SimulationEcsim::init_operators() {
    Simulation::init_operators();

    // const double dt = get_checked<double>(system_config, "Dt");
    // Mmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    // IMmat.resize(domain.total_size() * 3, domain.total_size() * 3);

    // Mmat = -0.25 * dt * dt * mesh.curlB * mesh.curlE;
    // IMmat = mesh.Imat - mesh.Mmat;
    // IMmat.makeCompressed();
}
void SimulationEcsim::init_fields() {
    fieldJp.resize(domain.size(), 3);
    fieldJp_full.resize(domain.size(), 3);
    fieldJe.resize(domain.size(), 3);

    fieldE.resize(domain.size(), 3);
    fieldEn.resize(domain.size(), 3);
    fieldEp.resize(domain.size(), 3);
    fieldB.resize(domain.size(), 3);
    fieldBn.resize(domain.size(), 3);
    fieldBInit.resize(domain.size(), 3);
    fieldBFull.resize(domain.size(), 3);
    fieldE_external.resize(domain.size(), 3);

    fieldJp.setZero();
    fieldJe.setZero();

    fieldE.setZero();
    fieldB.setZero();

    if (get_checked<int>(system_config, "StartFromTime") > 0) {
        read_fields_from_recovery(fieldE, fieldB);
    }

    fieldEn = fieldE;
    fieldBn = fieldB;

    bc_handler.init_electric_field(fieldE, domain);

    fieldE_external.setZero();
    if (auto e_cfg = create_electric_field_config(system_config, "ExternalFieldE")) {
        e_cfg->apply(fieldE_external, domain);
        std::cout << "Electric field config: " << fieldE_external.norm() << std::endl;
    }

    fieldBInit.setZero();
    if (auto b_cfg = create_magnetic_field_config(system_config, "ExternalFieldB")) {
        b_cfg->apply(fieldBInit, domain);
    }
    LmatX.resize(domain.total_size() * 3);
}

void SimulationEcsim::prepare_step(const int timestep) {
    RECORD_TIMER;
    const double dt = get_checked<double>(system_config, "Dt");
    for (auto &kv : species) {
        auto &sp = *kv.second;
        sp.diag.injection_energy = sp.inject_particles_step(sp.get_injection_distributions(), timestep, domain, dt,
                                                            system_config.value("co_locate_species", true));
    }

    dampingEnergy_ = damping_fields(fieldEn, fieldBn, domain, system_config);
    fieldE = fieldEn;
    fieldB = fieldBn;
    fieldJp.setZero();
    fieldJe.setZero();

    // #pragma omp parallel for
    //   for ( size_t i = 0; i < LmatX.size(); i++){
    //       for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it){
    //         it->second = 0.;
    //     }
    //   }

    for (auto &kv : species) {
        auto &sp = *kv.second;
        sp.prepare();   // save start coord for esirkepov current
    }
}

SimulationEcsim::~SimulationEcsim() = default;

void SimulationEcsim::make_diagnostic(const int timestep) {
    RECORD_TIMER;

    if (!diagnostic_ptr_) {
        nlohmann::json diagnostic_config =
            system_config.contains("diagnostics") ? system_config["diagnostics"] : nlohmann::json::object();
        diagnostic_ptr_ = std::make_unique<Diagnostics>(diagnostic_config, domain, species);
        outputs_ = OutputFactory::create(diagnostic_config, domain, species, system_config, fieldEn, fieldBn,
                                         fieldE_external, fieldBInit);
    }

    diagnostic_energy(*diagnostic_ptr_);

    for (auto &out : outputs_) {
        out->output(timestep, *diagnostic_ptr_);
    }
}

void SimulationEcsim::collect_per_species_diagnostics(Diagnostics &diagnostic, double &kineticEnergy,
                                                      double &kineticEnergyNew, double &totalLostEnergy) {
    static const Face all_faces[] = {Face::XMIN, Face::XMAX, Face::YMIN,    Face::YMAX,
                                     Face::ZMIN, Face::ZMAX, Face::CYLINDER};
    static const char *face_names[] = {"XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX", "CYLINDER"};
    static const PerFaceStats zero_stats;
    double totalInjectEnergy = 0;
    double totalEmitEnergy = 0;

    for (auto &kv : species) {
        auto &sp = *kv.second;
        diagnostic.addEnergy(sp.name(), sp.get_kinetic_energy());
        diagnostic.addEnergy(sp.name() + "Particles", sp.get_total_num_of_particles());
        diagnostic.addEnergy(sp.name() + "Inject", sp.diag.injection_energy);
        diagnostic.addEnergy(sp.name() + "Z", sp.get_kinetic_energy(Z));
        diagnostic.addEnergy(sp.name() + "XY", sp.get_kinetic_energy(X, Y));
        kineticEnergy += sp.kineticEnergy;
        kineticEnergyNew += diagnostic.energy[sp.name()];
        totalInjectEnergy += sp.diag.injection_energy;

        diagnostic.addBoundary(sp.name() + "_InjectE", sp.diag.injection_energy);
        diagnostic.addBoundary(sp.name() + "_InjectN", static_cast<double>(sp.diag.injection_count));

        double sp_lost_energy = 0;
        for (int fi = 0; fi < 7; ++fi) {
            Face f = all_faces[fi];
            auto it = sp.diag.boundary.find(f);
            const auto &s = (it != sp.diag.boundary.end()) ? it->second : zero_stats;
            std::string prefix = sp.name() + "_" + face_names[fi];
            diagnostic.addBoundary(prefix + "LostE", s.lost_energy);
            diagnostic.addBoundary(prefix + "LostN", static_cast<double>(s.lost_count));
            diagnostic.addBoundary(prefix + "ReflN", static_cast<double>(s.reflected_count));
            diagnostic.addBoundary(prefix + "EmitE", s.emitted_energy);
            diagnostic.addBoundary(prefix + "EmitN", static_cast<double>(s.emitted_count));
            sp_lost_energy += s.lost_energy;
            totalEmitEnergy += s.emitted_energy;
        }
        totalLostEnergy += sp_lost_energy;

        sp.kineticEnergy = sp.get_kinetic_energy();
        sp.diag.clear();
    }
    diagnostic.addEnergy("totalInjectEnergy", totalInjectEnergy);
    diagnostic.addEnergy("totalEmitEnergy", totalEmitEnergy);
}

void SimulationEcsim::compute_field_energy_and_conservation(Diagnostics &diagnostic, const IndexRange &irange,
                                                            double dt, double kineticEnergy, double kineticEnergyNew,
                                                            double totalLostEnergy, double totalInjectEnergy,
                                                            double energyJe_ex, double dampingEnergy) {
    diagnostic.addEnergy("totalLostEnergy", totalLostEnergy);
    diagnostic.addEnergy("energyDamping", dampingEnergy);
    diagnostic.addEnergy("energyFieldE", calc_energy_field(fieldEn, irange));
    diagnostic.addEnergy("energyFieldB", calc_energy_field(fieldBn, irange));
    fieldBFull.data() = fieldBn.data() + fieldBInit.data();
    diagnostic.addEnergy("energyFieldBFull", calc_energy_field(fieldBFull, irange));
    double energyFieldEold = calc_energy_field(fieldE, irange);
    double energyFieldBold = calc_energy_field(fieldB, irange);

    double energyFieldDifference =
        diagnostic.energy["energyFieldB"] + diagnostic.energy["energyFieldE"] - energyFieldBold - energyFieldEold;

    diagnostic.addEnergy("energyConserve", std::abs(kineticEnergyNew - kineticEnergy - totalInjectEnergy +
                                                    energyFieldDifference - dt * energyJe_ex + totalLostEnergy));
}

void SimulationEcsim::diagnostic_energy(Diagnostics &diagnostic) {
    double kineticEnergy = 0;
    double kineticEnergyNew = 0;
    double energyJe_ex = 0;
    double energyJe = 0;
    double totalLostEnergy = 0;
    IndexRange irange = bc_handler.active_range(domain.grid);

    collect_per_species_diagnostics(diagnostic, kineticEnergy, kineticEnergyNew, totalLostEnergy);

    for (auto &kv : species) {
        auto &sp = *kv.second;
        algorithmsECSIM::calculate_current(sp, sp.currentOnGrid);
        bc_handler.apply_to_fields(sp.currentOnGrid, FieldType::CURRENT, domain);
        energyJe_ex += dot_product_sum(fieldE_external, sp.currentOnGrid, irange);
        energyJe += dot_product_sum(fieldEp, sp.currentOnGrid, irange);
    }

    const double dt = get_checked<double>(system_config, "Dt");
    compute_field_energy_and_conservation(diagnostic, irange, dt, kineticEnergy, kineticEnergyNew, totalLostEnergy,
                                          diagnostic.energy["totalInjectEnergy"], energyJe_ex, dampingEnergy_);

    fieldJp_full.data() = fieldJp.data() + mesh.Lmat2 * (fieldE.data() + fieldEn.data()) / dt;
    double energyJe2 = dot_product_sum(fieldEp, fieldJp_full, irange);
    double energyFieldDifference = diagnostic.energy["energyFieldB"] + diagnostic.energy["energyFieldE"] -
                                   calc_energy_field(fieldE, irange) - calc_energy_field(fieldB, irange);

    diagnostic.addEnergy("deltaKinetic", kineticEnergyNew - kineticEnergy);
    diagnostic.addEnergy("deltaField", energyFieldDifference);
    diagnostic.addEnergy("JePredWork", dt * energyJe2);
    diagnostic.addEnergy("JeRealWork", dt * energyJe);
    diagnostic.addEnergy("JeExWork", dt * energyJe_ex);
}

void update_Lmat(std::vector<IndexMap> &LmatX, const Vector3R &coord, const Domain &domain, double charge, double mass,
                 double mpw, const Field3d &fieldB, const double dt) {
    const int SMAX = SHAPE_SIZE;
    double wx, wy, wz;
    int cellLocX, cellLocY, cellLocZ, cellLocX05, cellLocY05, cellLocZ05;
    double coordLocX, coordLocY, coordLocZ;
    double coordLocX05, coordLocY05, coordLocZ05;
    int i, j, k;
    int indx, indy, indz;
    int indx05, indy05, indz05;
    alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
    alignas(64) double sx05[SMAX], sy05[SMAX], sz05[SMAX];

    coordLocX = coord.x() / domain.cell_size().x() + GHOST_CELLS;
    coordLocY = coord.y() / domain.cell_size().y() + GHOST_CELLS;
    coordLocZ = coord.z() / domain.cell_size().z() + GHOST_CELLS;
    coordLocX05 = coordLocX - 0.5;
    coordLocY05 = coordLocY - 0.5;
    coordLocZ05 = coordLocZ - 0.5;

    cellLocX = int(coordLocX);
    cellLocY = int(coordLocY);
    cellLocZ = int(coordLocZ);
    cellLocX05 = int(coordLocX05);
    cellLocY05 = int(coordLocY05);
    cellLocZ05 = int(coordLocZ05);

    sx[1] = (coordLocX - cellLocX);
    sx[0] = 1 - sx[1];
    sy[1] = (coordLocY - cellLocY);
    sy[0] = 1 - sy[1];
    sz[1] = (coordLocZ - cellLocZ);
    sz[0] = 1 - sz[1];

    sx05[1] = (coordLocX05 - cellLocX05);
    sx05[0] = 1 - sx05[1];
    sy05[1] = (coordLocY05 - cellLocY05);
    sy05[0] = 1 - sy05[1];
    sz05[1] = (coordLocZ05 - cellLocZ05);
    sz05[0] = 1 - sz05[1];

    Vector3R B = Vector3R(0, 0, 0);

    for (i = 0; i < SMAX; ++i) {
        indx = cellLocX + i;
        indx05 = cellLocX05 + i;
        for (j = 0; j < SMAX; ++j) {
            indy = cellLocY + j;
            indy05 = cellLocY05 + j;
            for (k = 0; k < SMAX; ++k) {
                indz = cellLocZ + k;
                indz05 = cellLocZ05 + k;
                wx = sx[i] * sy05[j] * sz05[k];
                wy = sx05[i] * sy[j] * sz05[k];
                wz = sx05[i] * sy05[j] * sz[k];
                B.x() += (wx * fieldB(indx, indy05, indz05, 0));
                B.y() += (wy * fieldB(indx05, indy, indz05, 1));
                B.z() += (wz * fieldB(indx05, indy05, indz, 2));
            }
        }
    }

    const double q_m = charge / mass;
    const Vector3R b = 0.5 * dt * q_m * B;

    const double betaI = mpw * charge / (1.0 + b.squared());
    const double betaL = 0.5 * dt * q_m * betaI;

    const double matB[3][3] = {{1.0 + b.x() * b.x(), +b.z() + b.x() * b.y(), -b.y() + b.x() * b.z()},
                               {-b.z() + b.y() * b.x(), 1.0 + b.y() * b.y(), +b.x() + b.y() * b.z()},
                               {+b.y() + b.z() * b.x(), -b.x() + b.z() * b.y(), 1.0 + b.z() * b.z()}};

    constexpr double eps = 1.e-16;

    for (int i = 0; i < SMAX; ++i) {
        for (int j = 0; j < SMAX; ++j) {
            for (int k = 0; k < SMAX; ++k) {
                // веса и индексы для (i,j,k)
                const double s1[3] = {
                    sx05[i] * sy[j] * sz[k],   // X
                    sx[i] * sy05[j] * sz[k],   // Y
                    sx[i] * sy[j] * sz05[k]    // Z
                };

                const int id1[3] = {domain.vind(cellLocX05 + i, cellLocY + j, cellLocZ + k, 0),
                                    domain.vind(cellLocX + i, cellLocY05 + j, cellLocZ + k, 1),
                                    domain.vind(cellLocX + i, cellLocY + j, cellLocZ05 + k, 2)};

                for (int i1 = 0; i1 < SMAX; ++i1) {
                    for (int j1 = 0; j1 < SMAX; ++j1) {
                        for (int k1 = 0; k1 < SMAX; ++k1) {
                            const double s2[3] = {sx05[i1] * sy[j1] * sz[k1], sx[i1] * sy05[j1] * sz[k1],
                                                  sx[i1] * sy[j1] * sz05[k1]};

                            const int id2[3] = {domain.vind(cellLocX05 + i1, cellLocY + j1, cellLocZ + k1, 0),
                                                domain.vind(cellLocX + i1, cellLocY05 + j1, cellLocZ + k1, 1),
                                                domain.vind(cellLocX + i1, cellLocY + j1, cellLocZ05 + k1, 2)};

                            const double common = betaL;

                            // 3×3 вместо 9 копипаст
                            for (int c1 = 0; c1 < 3; ++c1) {
                                const int row = id1[c1];
                                const double w1 = s1[c1];
                                if (w1 == 0.0)
                                    continue;

                                for (int c2 = 0; c2 < 3; ++c2) {
                                    const double value = common * w1 * s2[c2] * matB[c1][c2];

                                    if (fabs(value) > eps) {
                                        LmatX[row][id2[c2]] += value;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void update_LmatNGP(std::vector<IndexMap> &LmatX, const Vector3R &coord, const Domain &domain, double charge,
                    double mass, double mpw, const Field3d &fieldB, const double dt) {
    const double coordLocX = coord.x() / domain.cell_size().x() + GHOST_CELLS;
    const double coordLocY = coord.y() / domain.cell_size().y() + GHOST_CELLS;
    const double coordLocZ = coord.z() / domain.cell_size().z() + GHOST_CELLS;
    const double coordLocX05 = coordLocX - 0.5;
    const double coordLocY05 = coordLocY - 0.5;
    const double coordLocZ05 = coordLocZ - 0.5;

    const int cellLocX = ngp(coordLocX);
    const int cellLocY = ngp(coordLocY);
    const int cellLocZ = ngp(coordLocZ);
    const int cellLocX05 = ngp(coordLocX05);
    const int cellLocY05 = ngp(coordLocY05);
    const int cellLocZ05 = ngp(coordLocZ05);

    const int indx = domain.vind(cellLocX05, cellLocY, cellLocZ, 0);
    const int indy = domain.vind(cellLocX, cellLocY05, cellLocZ, 1);
    const int indz = domain.vind(cellLocX, cellLocY, cellLocZ05, 2);
    Vector3R B = Vector3R(0.);

    B.x() = fieldB(cellLocX, cellLocY05, cellLocZ05, 0);
    B.y() = fieldB(cellLocX05, cellLocY, cellLocZ05, 1);
    B.z() = fieldB(cellLocX05, cellLocY05, cellLocZ, 2);

    const double q_m = charge / mass;
    const Vector3R b = 0.5 * dt * q_m * B;

    const double betaI = mpw * charge / (1.0 + b.squared());
    const double betaL = 0.5 * dt * q_m * betaI;

    const double matB[3][3] = {{1.0 + b.x() * b.x(), +b.z() + b.x() * b.y(), -b.y() + b.x() * b.z()},
                               {-b.z() + b.y() * b.x(), 1.0 + b.y() * b.y(), +b.x() + b.y() * b.z()},
                               {+b.y() + b.z() * b.x(), -b.x() + b.z() * b.y(), 1.0 + b.z() * b.z()}};

    constexpr double eps = 1.e-16;
    const double common = betaL;

    const int id[3] = {indx, indy, indz};

    for (int c1 = 0; c1 < 3; ++c1) {
        const int row = id[c1];
        for (int c2 = 0; c2 < 3; ++c2) {
            const double value = common * matB[c1][c2];
            if (fabs(value) > eps) {
                LmatX[row][id[c2]] += value;
            }
        }
    }
}
