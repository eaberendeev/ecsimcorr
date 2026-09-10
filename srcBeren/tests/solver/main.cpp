// Тесты солвера (этап 1: темы 1+4) — быстрые, без полного PIC.
// Малая сетка 6x6x6 (9^3 узлов, 3*729 = 2187 неизвестных), периодика по всем осям.
// Проверяют: сопряжённость роторов, Jacobi vs identity, дефектную переформулировку
// corrector-шага, экстраполяцию x0, mixed precision vs эталон, и энергетическое
// тождество (MATH.md §7): dU_E + dU_B + dt*E_avg*J_e = E_avg*r.
#include <Eigen/SparseLU>
#include <cmath>
#include <iostream>
#include <random>
#include <string>

#include "Mesh.h"
#include "World.h"
#include "boundary_conditions.h"
#include "containers.h"
#include "solverSLE.h"

namespace test {
static int passed = 0;
static int failed = 0;

void assert_true(bool cond, const std::string& name) {
    if (cond) {
        passed++;
    } else {
        std::cout << "[FAIL] " << name << "\n";
        failed++;
    }
}

void assert_le(double actual, double limit, const std::string& name) {
    if (actual <= limit) {
        passed++;
    } else {
        std::cout << "[FAIL] " << name << " (got " << actual << ", limit " << limit << ")\n";
        failed++;
    }
}

void print_summary() {
    std::cout << "\n========================================\n";
    std::cout << "Test Summary: " << passed << " passed, " << failed << " failed\n";
    std::cout << "========================================\n";
}
}   // namespace test

namespace {

const int NC = 6;              // ячеек по осям
const double DX = 0.25;        // шаг сетки
const double DT = 1.0;         // временной шаг

struct Setup {
    Domain domain;
    BoundaryConditionHandler bc;
    Mesh mesh;
    int n;   // 3 * число узлов

    Setup() {
        domain.init(Vector3I(NC, NC, NC), Vector3R(DX, DX, DX));
        nlohmann::json cfg = R"({
            "Boundary_conditions": [
                {"periodic": [{"face": "XMIN"}, {"face": "YMIN"}, {"face": "ZMIN"}]}
            ]
        })"_json;
        bc.load_from_json(cfg, domain);
        mesh.init(domain, DT, bc);
        n = 3 * domain.total_size();
    }
};

void fill_random(Field3d& f, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> ud(-1.0, 1.0);
    for (int i = 0; i < f.size(); ++i) f(i) = ud(rng);
}

void fill_random_small(Field3d& f, std::mt19937_64& rng, double amp) {
    std::uniform_real_distribution<double> ud(-1.0, 1.0);
    for (int i = 0; i < f.size(); ++i) f(i) = amp * ud(rng);
}

// Решение системы через SparseLU (эталон для малой задачи; SparseLU требует column-major)
Field3d solve_dense(const Operator& A, const Field3d& rhs) {
    Eigen::SparseMatrix<double, Eigen::ColMajor> Ac = A;
    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> lu;
    lu.compute(Ac);
    Field3d x(rhs.size());
    x.data() = lu.solve(rhs.data());
    return x;
}

Field3d ones_diag(int n) {
    Field3d d(n);
    for (int i = 0; i < n; ++i) d(i) = 1.0;
    return d;
}

// Истинная диагональ IMmat: 1 + 0.25*dt^2*sum_j curlB(i,j)^2
Field3d imat_diag(const Operator& curlB, double dt, int n) {
    Field3d d(n);
    for (int i = 0; i < n; ++i) d(i) = 1.0;
    for (int k = 0; k < curlB.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(curlB, k); it; ++it) {
            d(it.row()) += 0.25 * dt * dt * it.value() * it.value();
        }
    }
    return d;
}

double field_energy(const Field3d& f) {
    return 0.5 * f.dot(f);
}

// --- Тест 1: сопряжённость роторов ------------------------------------------
// В этом коде curlB = curlE^T (entry-wise) для строк, чьи столбцы имеют строки
// в curlE. Проверяем поэлементно на активном диапазоне периодической сетки.
// (Причина: curlB использует разности назад, curlE — минус разности вперёд;
// на ghost-строках периодической сетки wrap_index ломает парность — они исключены.)
void test_adjointness(const Setup& s, const std::vector<int>& active) {
    double max_diff = 0.0;
    long checked = 0;
    for (int r : active) {
        for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(s.mesh.curlB, r); it; ++it) {
            const double expect = it.value();
            const double got = s.mesh.curlE.coeff(it.col(), r);
            max_diff = std::max(max_diff, std::abs(got - expect));
            ++checked;
        }
    }
    std::cout << "  pairs checked = " << checked << ", max |(curlB)_rb - (curlE)_br| = " << max_diff << "\n";
    test::assert_le(max_diff, 1e-13, "adjointness: curlB == curlE^T (entry-wise, active rows)");
}

// --- Тест 2: Jacobi vs identity (решения совпадают, итерации печатаем) -------
void test_jacobi(const Setup& s) {
    const auto& A = s.mesh.IMmat;
    Field3d rhs(s.n), x0(s.n), x_id(s.n), x_jac(s.n);
    std::mt19937_64 rng(123);
    fill_random(rhs, rng);
    fill_random_small(x0, rng, 0.1);

    size_t iters_id = 1000, iters_jac = 1000;
    double err_id = 1e-10, err_jac = 1e-10;
    x_id = x0;
    x_jac = x0;
    bicgstab_iteration<Field3d>(A, rhs, x_id, ones_diag(s.n), iters_id, err_id);
    Field3d diag = imat_diag(s.mesh.curlB, DT, s.n);
    bicgstab_iteration<Field3d>(A, rhs, x_jac, diag, iters_jac, err_jac);

    std::cout << "  iters: identity=" << iters_id << " jacobi=" << iters_jac << "\n";
    test::assert_le((x_id - x_jac).norm() / x_jac.norm(), 1e-8, "jacobi: solution equals identity solution");
    Field3d resid(s.n);
    resid.data() = rhs.data() - A * x_jac.data();
    test::assert_le(resid.norm() / rhs.norm(), 1e-8, "jacobi: true residual below tolerance");
}

// --- Тест 3: дефектная переформулировка corrector ---------------------------
void test_defect(const Setup& s) {
    const auto& M = s.mesh.Mmat;
    const auto& IM = s.mesh.IMmat;
    const auto& curlB = s.mesh.curlB;

    std::mt19937_64 rng(456);
    Field3d E(s.n), B(s.n), Je(s.n), Jp(s.n), Eex(s.n);
    fill_random(E, rng);
    fill_random(B, rng);
    fill_random(Je, rng);
    fill_random_small(Eex, rng, 0.1);
    Jp = Je;   // J' близко к J_e (как в реальной симуляции)
    Field3d pert(s.n);
    fill_random_small(pert, rng, 0.1);
    Jp += pert;

    // L = положительная диагональная матрица
    Operator L(s.n, s.n);
    std::uniform_real_distribution<double> ud(0.0, 1.0);
    for (int i = 0; i < s.n; ++i) L.insert(i, i) = 0.1 + 0.2 * ud(rng);
    L.makeCompressed();

    // predictor: (I - M + L) Ep = E + 0.5dt(curlB*B - J') - L*E_ex
    Operator A = IM + L;
    Field3d rhs_pred(s.n), Ep(s.n);
    rhs_pred.data() = E.data() + 0.5 * DT * (curlB * B.data() - Jp.data()) - L * Eex.data();
    size_t iters = 1000;
    double err = 1e-11;
    Ep = E;
    bicgstab_iteration<Field3d>(A, rhs_pred, Ep, ones_diag(s.n), iters, err);
    std::cout << "  predictor iters = " << iters << ", err = " << err << "\n";

    // direct corrector: (I-M) E1 = E + dt(curlB*B - J_e) + M*E
    Field3d rhs_dir(s.n), E_dir(s.n);
    rhs_dir.data() = E.data() + DT * (curlB * B.data() - Je.data()) + M * E.data();
    iters = 1000;
    err = 1e-11;
    E_dir = E;
    bicgstab_iteration<Field3d>(IM, rhs_dir, E_dir, ones_diag(s.n), iters, err);
    std::cout << "  direct corrector iters = " << iters << ", err = " << err << "\n";

    // defect corrector: (I-M) delta = dt(J' - J_e) + 2 L (Ep + E_ex)
    Field3d rhs_def(s.n), delta(s.n), E_def(s.n);
    rhs_def.data() = DT * (Jp.data() - Je.data()) + 2.0 * L * (Ep.data() + Eex.data());
    iters = 1000;
    err = 1e-11;
    delta.setZero();
    bicgstab_iteration<Field3d>(IM, rhs_def, delta, ones_diag(s.n), iters, err);
    E_def = 2.0 * Ep - E + delta;
    std::cout << "  defect corrector iters = " << iters << ", err = " << err << "\n";
    std::cout << "  |defect rhs| / |direct rhs| = " << rhs_def.norm() / rhs_dir.norm() << "\n";

    test::assert_le((E_def - E_dir).norm() / E_dir.norm(), 1e-9, "defect: E_dir == 2Ep - E + delta");
    test::assert_le(rhs_def.norm() / rhs_dir.norm(), 0.5, "defect: rhs smaller than direct rhs");
}

// --- Тест 4: экстраполяция x0 (решения совпадают) ---------------------------
void test_extrapolated_x0(const Setup& s) {
    const auto& IM = s.mesh.IMmat;
    std::mt19937_64 rng(789);
    Field3d E(s.n), E_prev(s.n), rhs(s.n);
    fill_random(E, rng);
    fill_random_small(E_prev, rng, 0.1);
    fill_random(rhs, rng);

    size_t iters_plain = 1000, iters_extr = 1000;
    double err_plain = 1e-11, err_extr = 1e-11;
    Field3d x_plain(s.n), x_extr(s.n);
    x_plain = E;
    bicgstab_iteration<Field3d>(IM, rhs, x_plain, ones_diag(s.n), iters_plain, err_plain);
    // линейная экстраполяция на полушаг: 1.5 E - 0.5 E_prev
    x_extr = 1.5 * E - 0.5 * E_prev;
    bicgstab_iteration<Field3d>(IM, rhs, x_extr, ones_diag(s.n), iters_extr, err_extr);

    std::cout << "  iters: plain=" << iters_plain << " extrapolated=" << iters_extr << "\n";
    test::assert_le((x_extr - x_plain).norm() / x_plain.norm(), 1e-8, "extrapolated x0: same solution");
}

// --- Тест 5: mixed precision == эталон (SparseLU) ----------------------------
void test_mixed_precision(const Setup& s) {
    const auto& IM = s.mesh.IMmat;
    std::mt19937_64 rng(101112);
    Field3d rhs(s.n), x(s.n);
    fill_random(rhs, rng);
    fill_random_small(x, rng, 0.1);

    size_t iters = 1000;
    double err = 1e-10;
    bicgstab_iteration<Field3d>(IM, rhs, x, ones_diag(s.n), iters, err);

    Field3d x_ref = solve_dense(IM, rhs);
    double rel = (x - x_ref).norm() / x_ref.norm();
    std::cout << "  iters = " << iters << ", err = " << err << ", |x - x_ref|/|x_ref| = " << rel << "\n";
    test::assert_le(rel, 2e-8, "mixed precision: matches dense reference");
}

// --- Тест 6: энергетическое тождество (MATH.md §7) --------------------------
// Для активных E-строк R_E (оси 1..Ncells+1) и b-set (B-столбцы, на которые
// они ссылаются, включая ghost-B) выполнено ТОЧНО при ЛЮБОЙ невязке r:
//   dU_E + dU_B + dt*E_avg*J_e = SUM_b G(b)*[0.5*dt^2*(curlE E_avg)(b) - dt*B(b)] - E_avg*r,
// где G(b) = SUM_{r не в R_E} (curlE)_{b,r}*E_avg(r) — вклад ghost-E строк
// (на периодической сетке их строки curlB ломаются wrap'ом, поэтому G берём
// через столбцы curlE). В частности, дрейф полной энергии при неточном
// солве определяется слагаемым -E_avg*r.
void test_energy_identity(const Setup& s, const std::vector<int>& active) {
    const auto& M = s.mesh.Mmat;
    const auto& IM = s.mesh.IMmat;
    const auto& curlB = s.mesh.curlB;
    const auto& curlE = s.mesh.curlE;

    // b-set: B-узлы, на которые ссылаются активные E-строки
    std::vector<char> is_b(s.n, 0);
    for (int r : active) {
        for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(curlB, r); it; ++it) is_b[it.col()] = 1;
    }
    std::vector<int> bset;
    for (int i = 0; i < s.n; ++i)
        if (is_b[i]) bset.push_back(i);

    std::vector<char> is_active(s.n, 0);
    for (int r : active) is_active[r] = 1;

    std::mt19937_64 rng(131415);
    Field3d E(s.n), B(s.n), Je(s.n);
    fill_random(E, rng);
    fill_random(B, rng);
    fill_random(Je, rng);

    Field3d rhs(s.n);
    rhs.data() = E.data() + DT * (curlB * B.data() - Je.data()) + M * E.data();

    auto dot_restricted = [](const Field3d& a, const Field3d& b, const std::vector<int>& idx) {
        double sum = 0.0;
        for (int i : idx) sum += a(i) * b(i);
        return sum;
    };

    for (double tol : {1e-13, 1e-3}) {
        size_t iters = 1000;
        double err = tol;
        Field3d E1 = E;
        bicgstab_iteration<Field3d>(IM, rhs, E1, ones_diag(s.n), iters, err);

        // истинная невязка r
        Field3d r(s.n);
        r.data() = rhs.data() - IM * E1.data();

        // B_{n+1} = B_n - 0.5 dt curlE (E_n + E_{n+1})  (как в коде)
        Field3d Esum = E + E1;
        Field3d B1(s.n);
        B1.data() = B.data() - 0.5 * DT * (curlE * Esum.data());

        Field3d Eavg = 0.5 * Esum;
        Field3d curlE_Eavg(s.n);
        curlE_Eavg.data() = curlE * Eavg.data();

        // ghost-коррекция: G(b) = сумма (curlE)_{b,r}*E_avg(r) по r вне R_E.
        // Точное тождество (проверено: свёртка Ампера с E_avg и преобразования
        // парности совпадают с машинной точностью):
        // dU_E + dU_B + dt*E_avg*J = SUM_b G(b)*[0.5*dt^2*(curlE E_avg)(b) - dt*B(b)] - E_avg*r
        Field3d G(s.n);
        G.setZero();
        for (int b = 0; b < s.n; ++b) {
            for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(curlE, b); it; ++it) {
                const int r = it.col();
                if (!is_active[r]) G(b) += it.value() * Eavg(r);
            }
        }
        double corr = 0.0;
        for (int b : bset) corr += G(b) * (0.5 * DT * DT * curlE_Eavg(b) - DT * B(b));

        double dUE = 0.0, dUB = 0.0;
        for (int i : active) dUE += 0.5 * (E1(i) * E1(i) - E(i) * E(i));
        for (int b : bset) dUB += 0.5 * (B1(b) * B1(b) - B(b) * B(b));
        double work = DT * dot_restricted(Eavg, Je, active);
        double lhs = dUE + dUB + work;
        double rhs_term = dot_restricted(Eavg, r, active);
        double scale = std::max(1.0, std::abs(dUE) + std::abs(dUB) + std::abs(work));

        std::cout << "  tol=" << tol << ": |r|/|rhs|=" << r.norm() / rhs.norm()
                  << ", |lhs + E_avg*r - corr|/scale = " << std::abs(lhs + rhs_term - corr) / scale << "\n";
        test::assert_le(std::abs(lhs + rhs_term - corr) / scale, 1e-12,
                        "energy identity: dU_E + dU_B + dt*E_avg*J_e == corr - E_avg*r");
    }
}

// --- Тест 7: предобусловленный BiCGSTAB (M3) --------------------------------
void test_precond(const Setup& s) {
    const auto& IM = s.mesh.IMmat;
    std::mt19937_64 rng(161718);
    Field3d rhs(s.n), x0(s.n);
    fill_random(rhs, rng);
    fill_random_small(x0, rng, 0.1);

    Field3d x_ref = solve_dense(IM, rhs);

    // (a) точный LU-предобуславливатель: должен сходиться за 1-2 итерации
    SparseLUPreconditioner pc(IM);
    {
        // Прямая проверка факторизации: lu.solve должен давать решение сразу
        Field3d x_direct(s.n);
        pc.apply(rhs, x_direct);
        double rel_direct = (x_direct - x_ref).norm() / x_ref.norm();
        std::cout << "  lu direct solve rel = " << rel_direct << "\n";

        size_t iters = 1000;
        double err = 1e-10;
        Field3d x = x0;
        bicgstab_iteration_precond(IM, rhs, x, pc, iters, err);
        double rel = (x - x_ref).norm() / x_ref.norm();
        std::cout << "  lu precond: iters=" << iters << " err=" << err << " rel=" << rel << "\n";
        test::assert_le(rel, 1e-8, "precond lu: matches dense reference");
        test::assert_le((double)iters, 5.0, "precond lu: converges in <=5 iterations");
    }
    // (b) Jacobi через тот же путь
    {
        Field3d diag = imat_diag(s.mesh.curlB, DT, s.n);
        JacobiPreconditioner pj(diag);
        size_t iters = 1000;
        double err = 1e-10;
        Field3d x = x0;
        bicgstab_iteration_precond(IM, rhs, x, pj, iters, err);
        double rel = (x - x_ref).norm() / x_ref.norm();
        std::cout << "  jacobi precond: iters=" << iters << " err=" << err << " rel=" << rel << "\n";
        test::assert_le(rel, 1e-8, "precond jacobi: matches dense reference");
    }
}

// --- Тест 8: энерго-осведомлённый (адаптивный) допуск солвера ---------------
// Следствие теоремы (T): |дрейф| <= |E_avg|*|r| <= |E_avg|*tol*|rhs|.
// Для заданного бюджета сохранения eps достаточно tol = eps/(|E_avg|*|rhs|).
// Решаем с таким tol и проверяем: фактический дрейф <= eps.
void test_adaptive_tolerance(const Setup& s, const std::vector<int>& active) {
    const auto& M = s.mesh.Mmat;
    const auto& IM = s.mesh.IMmat;
    const auto& curlB = s.mesh.curlB;
    const auto& curlE = s.mesh.curlE;

    std::vector<char> is_b(s.n, 0);
    for (int r : active)
        for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(curlB, r); it; ++it) is_b[it.col()] = 1;
    std::vector<int> bset;
    for (int i = 0; i < s.n; ++i)
        if (is_b[i]) bset.push_back(i);
    std::vector<char> is_active(s.n, 0);
    for (int r : active) is_active[r] = 1;

    std::mt19937_64 rng(192021);
    Field3d E(s.n), B(s.n), Je(s.n);
    fill_random(E, rng);
    fill_random(B, rng);
    fill_random(Je, rng);
    Field3d rhs(s.n);
    rhs.data() = E.data() + DT * (curlB * B.data() - Je.data()) + M * E.data();

    const double U0 = 0.5 * (E.dot(E) + B.dot(B));

    for (double eps_rel : {1e-3, 1e-6}) {
        const double eps = eps_rel * U0;
        // консервативная априорная оценка |E_avg| <= |E_n| + dt*|rhs|/2
        const double Eavg_bound = E.norm() + 0.5 * DT * rhs.norm();
        const double tol = eps / (Eavg_bound * rhs.norm());

        size_t iters = 1000;
        double err = tol;
        Field3d E1 = E;
        bicgstab_iteration<Field3d>(IM, rhs, E1, ones_diag(s.n), iters, err);

        Field3d r(s.n);
        r.data() = rhs.data() - IM * E1.data();

        Field3d Esum = E + E1;
        Field3d B1(s.n);
        B1.data() = B.data() - 0.5 * DT * (curlE * Esum.data());
        Field3d Eavg = 0.5 * Esum;
        Field3d curlE_Eavg(s.n);
        curlE_Eavg.data() = curlE * Eavg.data();

        Field3d G(s.n);
        G.setZero();
        for (int b = 0; b < s.n; ++b)
            for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(curlE, b); it; ++it)
                if (!is_active[it.col()]) G(b) += it.value() * Eavg(it.col());
        double corr = 0.0;
        for (int b : bset) corr += G(b) * (0.5 * DT * DT * curlE_Eavg(b) - DT * B(b));

        double dUE = 0.0, dUB = 0.0;
        for (int i : active) dUE += 0.5 * (E1(i) * E1(i) - E(i) * E(i));
        for (int b : bset) dUB += 0.5 * (B1(b) * B1(b) - B(b) * B(b));
        double work = 0.0;
        for (int i : active) work += DT * Eavg(i) * Je(i);

        const double drift = dUE + dUB + work - corr;   // = -E_avg*r по теореме
        std::cout << "  eps_rel=" << eps_rel << ": tol=" << tol << " iters=" << iters
                  << " |drift|/U0=" << std::abs(drift) / U0 << " (limit " << eps_rel << ")\n";
        test::assert_le(std::abs(drift), eps * 1.05 + 1e-13,
                        "adaptive tol: |energy drift| <= eps (theorem bound)");
    }
}

// --- Тест 9: пропуск corrector (P1) -----------------------------------------
// Если E_{n+1} := 2E' − E_n (корректор пропущен, δ = 0), то истинная невязка
// corrector-системы РАВНА дефекту rhs_def = dt(J'−J_e) + 2L(E'+E_ex). Отсюда
// дрейф энергии = −⟨E_avg, rhs_def⟩ и |дрейф| <= |E_avg|·|rhs_def| — бюджет,
// которым управляется пропуск (SKIP_CORR_TOL в коде).
void test_corrector_skip(const Setup& s, const std::vector<int>& active) {
    const auto& M = s.mesh.Mmat;
    const auto& IM = s.mesh.IMmat;
    const auto& curlB = s.mesh.curlB;
    const auto& curlE = s.mesh.curlE;

    std::mt19937_64 rng(222324);
    Field3d E(s.n), B(s.n), Je(s.n), Jp(s.n), Eex(s.n);
    fill_random(E, rng);
    fill_random(B, rng);
    fill_random(Je, rng);
    fill_random_small(Eex, rng, 0.1);
    Jp = Je;
    Field3d pert(s.n);
    fill_random_small(pert, rng, 0.1);
    Jp += pert;

    Operator L(s.n, s.n);
    std::uniform_real_distribution<double> ud(0.0, 1.0);
    for (int i = 0; i < s.n; ++i) L.insert(i, i) = 0.1 + 0.2 * ud(rng);
    L.makeCompressed();

    // predictor: ((I−M)+L)·E' = E + ½dt(curlB·B − J') − L·E_ex
    Operator A = IM + L;
    Field3d rhs_pred(s.n), Ep(s.n);
    rhs_pred.data() = E.data() + 0.5 * DT * (curlB * B.data() - Jp.data()) - L * Eex.data();
    size_t iters = 1000;
    double err = 1e-11;
    Ep = E;
    bicgstab_iteration<Field3d>(A, rhs_pred, Ep, ones_diag(s.n), iters, err);

    // corrector rhs и дефект
    Field3d rhs_corr(s.n);
    rhs_corr.data() = E.data() + DT * (curlB * B.data() - Je.data()) + M * E.data();
    Field3d rhs_def(s.n);
    rhs_def.data() = DT * (Jp.data() - Je.data()) + 2.0 * L * (Ep.data() + Eex.data());

    // ПРОПУСК: E1 = 2E' − E
    Field3d E1 = 2.0 * Ep - E;
    Field3d r(s.n);
    r.data() = rhs_corr.data() - IM * E1.data();

    const double rel_identity = (r - rhs_def).norm() / rhs_corr.norm();
    std::cout << "  skip: |r - rhs_def|/|rhs| = " << rel_identity
              << ", |defect|/|rhs| = " << rhs_def.norm() / rhs_corr.norm() << "\n";
    // тождество верно с точностью до невязки predictor-солва (~1e-11 => здесь ~2.8e-12)
    test::assert_le(rel_identity, 1e-10, "corrector skip: residual == defect (identity)");

    // бюджет: |дрейф| <= |E_avg|·|rhs_def| (с ghost-коррекцией)
    std::vector<char> is_b(s.n, 0);
    for (int rr : active)
        for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(curlB, rr); it; ++it) is_b[it.col()] = 1;
    std::vector<int> bset;
    for (int i = 0; i < s.n; ++i)
        if (is_b[i]) bset.push_back(i);
    std::vector<char> is_active(s.n, 0);
    for (int rr : active) is_active[rr] = 1;

    Field3d Esum = E + E1;
    Field3d B1(s.n);
    B1.data() = B.data() - 0.5 * DT * (curlE * Esum.data());
    Field3d Eavg = 0.5 * Esum;
    Field3d curlE_Eavg(s.n);
    curlE_Eavg.data() = curlE * Eavg.data();
    Field3d G(s.n);
    G.setZero();
    for (int b = 0; b < s.n; ++b)
        for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(curlE, b); it; ++it)
            if (!is_active[it.col()]) G(b) += it.value() * Eavg(it.col());
    double corr = 0.0;
    for (int b : bset) corr += G(b) * (0.5 * DT * DT * curlE_Eavg(b) - DT * B(b));
    double dUE = 0.0, dUB = 0.0, work = 0.0;
    for (int i : active) {
        dUE += 0.5 * (E1(i) * E1(i) - E(i) * E(i));
        work += DT * Eavg(i) * Je(i);
    }
    for (int b : bset) dUB += 0.5 * (B1(b) * B1(b) - B(b) * B(b));
    const double drift = dUE + dUB + work - corr;

    double Eavg_norm = 0.0, def_norm = 0.0;
    for (int i : active) Eavg_norm += Eavg(i) * Eavg(i);
    for (int i = 0; i < s.n; ++i) def_norm += rhs_def(i) * rhs_def(i);
    const double bound = sqrt(Eavg_norm * def_norm);
    std::cout << "  skip: |drift|=" << std::abs(drift) << ", bound |E_avg|*|defect|=" << bound << "\n";
    test::assert_le(std::abs(drift), bound + 1e-12, "corrector skip: |drift| <= |E_avg|*|defect|");
}

// --- Тест 10: замкнутая форма (λ, u) E+p коррекции --------------------------
// v' = λ·v + u должна давать заданные K_t и P_t (4 уравнения, 4 неизвестных).
// λ² = (K_t − |P_t|²/(2mN)) / (½m·W2),  u = (P_t/m − λ·V1)/N,
// N = Σw, V1 = Σw·v, V2 = Σw|v|², W2 = V2 − |V1|²/N.
void test_lambda_u_math() {
    std::mt19937_64 rng(252627);
    std::uniform_real_distribution<double> ud(-1.0, 1.0);

    for (int trial = 0; trial < 20; ++trial) {
        const int np = 1000 + int(1000 * ud(rng) * 0 + 500);
        const double m = 0.5 + ud(rng);
        double N = 0, V2 = 0, V1x = 0, V1y = 0, V1z = 0, K = 0, Px = 0, Py = 0, Pz = 0;
        std::vector<Vector3R> vs(np);
        std::vector<double> ws(np);
        for (int i = 0; i < np; ++i) {
            Vector3R v(ud(rng), ud(rng), ud(rng));
            double w = 0.5 + ud(rng);
            vs[i] = v;
            ws[i] = w;
            N += w;
            V1x += w * v.x();
            V1y += w * v.y();
            V1z += w * v.z();
            V2 += w * v.squared();
            K += 0.5 * m * w * v.squared();
            Px += m * w * v.x();
            Py += m * w * v.y();
            Pz += m * w * v.z();
        }
        // цели: малые возмущения
        const double K_t = K + 0.05 * ud(rng) * K + 1e-9;
        const double Ptx = Px + 0.05 * ud(rng) * Px + 1e-9;
        const double Pty = Py + 0.05 * ud(rng) * Py + 1e-9;
        const double Ptz = Pz + 0.05 * ud(rng) * Pz + 1e-9;

        const double W2 = V2 - (V1x * V1x + V1y * V1y + V1z * V1z) / N;
        const double num = K_t - 0.5 * (Ptx * Ptx + Pty * Pty + Ptz * Ptz) / (m * N);
        const double lambda = sqrt(num / (0.5 * m * W2));
        const double ux = (Ptx / m - lambda * V1x) / N;
        const double uy = (Pty / m - lambda * V1y) / N;
        const double uz = (Ptz / m - lambda * V1z) / N;

        double K_new = 0, Px_new = 0, Py_new = 0, Pz_new = 0;
        for (int i = 0; i < np; ++i) {
            Vector3R vn = lambda * vs[i] + Vector3R(ux, uy, uz);
            K_new += 0.5 * m * ws[i] * vn.squared();
            Px_new += m * ws[i] * vn.x();
            Py_new += m * ws[i] * vn.y();
            Pz_new += m * ws[i] * vn.z();
        }
        const double scale = std::abs(K_t) + std::abs(Ptx) + std::abs(Pty) + std::abs(Ptz) + 1e-30;
        const double rel = (std::abs(K_new - K_t) + std::abs(Px_new - Ptx) + std::abs(Py_new - Pty) +
                            std::abs(Pz_new - Ptz)) / scale;
        test::assert_le(rel, 1e-12, "lambda-u closed form: K and P match targets");
    }
    std::cout << "  lambda-u closed form: 20 random cases OK\n";
}

}   // namespace

int main() {
    std::cout << "solver stage-1 tests\n";
    Setup s;
    std::cout << "system size N = " << s.n << "\n";

    // Активные E-строки периодической сетки: оси в [1, Ncells+1]
    // (индексы сетки 0..Ncells+2; ghost-строки исключены из сумм)
    const int Ntot = NC + 3;
    std::vector<int> active;
    for (int i = 1; i <= NC + 1; ++i)
        for (int j = 1; j <= NC + 1; ++j)
            for (int k = 1; k <= NC + 1; ++k)
                for (int d = 0; d < 3; ++d) active.push_back(d + 3 * (i * Ntot * Ntot + j * Ntot + k));

    test_adjointness(s, active);
    test_jacobi(s);
    test_defect(s);
    test_extrapolated_x0(s);
    test_mixed_precision(s);
    test_energy_identity(s, active);
    test_precond(s);
    test_adaptive_tolerance(s, active);
    test_corrector_skip(s, active);
    test_lambda_u_math();

    test::print_summary();
    return test::failed == 0 ? 0 : 1;
}
