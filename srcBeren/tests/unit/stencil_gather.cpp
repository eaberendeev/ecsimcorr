#include "coupling/current.h"
#include "coupling/gather.h"
#include "raw_gather_reference.h"
#include "test_case.h"

namespace {
// Absolute roundoff bound for the small tensor sums (field magnitudes below 100).
constexpr double tolerance = 1e-13;

void vector_near(const Vector3R& actual, const Vector3R& expected, const std::string& name) {
    for (int d = 0; d < 3; ++d)
        check::near(actual[d], expected[d], tolerance, name + " component " + std::to_string(d));
}

template <ShapeType ST>
void shape_moments() {
    // Coordinates are in cell units. Reproduce a constant and a linear function
    // on both integer and half-integer nodes, including the origin crossings.
    for (double x : {-1.500001, -1.5, -1.499999, -0.500001, -0.5, -0.499999, -0.000001, 0.0, 0.000001, 0.499999, 0.5,
                     0.500001, 2.37}) {
        for (double shift : {0.0, 0.5}) {
            const Vector3R position(x, x + 0.17, x - 0.23);
            const auto stencil = make_point_stencil<ST>(position, Vector3R(shift));
            constexpr int N = ST == ShapeType::Linear ? 2 : 3;
            const auto fixed = make_stencil_at<ShapeKernel<ST>::weight, N>(position, stencil.origin, Vector3R(shift));
            for (int axis = 0; axis < 3; ++axis) {
                double sum = 0, moment = 0;
                for (int n = 0; n < N; ++n) {
                    const double w = stencil(n, static_cast<Axis>(axis));
                    check::require(w >= 0 && std::isfinite(w), "finite nonnegative weight");
                    sum += w;
                    moment += (stencil.node(n, static_cast<Axis>(axis)) + shift) * w;
                    check::near(fixed(n, static_cast<Axis>(axis)), w, 0, "given-origin factory applies shift");
                }
                check::near(sum, 1, tolerance, "partition of unity");
                check::near(moment, position[axis], tolerance, "first spatial moment");
            }
        }
    }
    // An explicit wider origin must be used as supplied, with the shift applied
    // once. Expected weights below follow the polynomial definition directly.
    const auto wide = make_stencil_at<ShapeKernel<ST>::weight, 4>(Vector3R(0.25, 0.5, 0.75), Vector3I(-1));
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
        double sum = 0;
        for (int local = 0; local < 4; ++local) sum += wide(local, axis);
        check::near(sum, 1, tolerance, "four-node stencil preserves unity");
    }
    const auto fixed = make_stencil_at<ShapeKernel<ST>::weight, 4>(Vector3R(0.75), Vector3I(-1), Vector3R(0.5));
    check::require(fixed.origin == Vector3I(-1), "explicit logical origin is preserved");
    const std::array<double, 4> expected = ST == ShapeType::Linear ? std::array<double, 4>{0, 0.75, 0.25, 0}
                                                                   : std::array<double, 4>{0.03125, 0.6875, 0.28125, 0};
    for (int axis = 0; axis < 3; ++axis)
        for (int n = 0; n < 4; ++n)
            check::near(fixed(n, static_cast<Axis>(axis)), expected[n], tolerance, "given-origin polynomial weights");
}

double affine(const Vector3R& p, int component) {
    return 0.7 + component * 1.1 + (component + 1) * p.x() - 0.3 * p.y() + 0.2 * p.z();
}

// Independent Yee locations: E has half an interval in its component direction;
// B has half an interval in the other two directions.
Vector3R yee_position(const Grid& grid, Vector3I node, int component, bool magnetic) {
    Vector3R p;
    for (int a = 0; a < 3; ++a)
        p[a] = (node[a] + ((magnetic ? a != component : a == component) ? 0.5 : 0.0)) * grid.cell_size(a);
    return p;
}

template <ShapeType ST>
void gather_fields() {
    constexpr int ghost = GHOST_CELLS;
    {
        const Grid grid(Vector3R(0.5, 1.25, 2), Vector3I(5, 6, 7), ghost);
        Field3d electric(grid.size(), 3), magnetic(grid.size(), 3);
        for (bool constant : {true, false}) {
            for (int i = -ghost; i <= 5 + ghost; ++i)
                for (int j = -ghost; j <= 6 + ghost; ++j)
                    for (int k = -ghost; k <= 7 + ghost; ++k) {
                        const auto raw = grid.to_storage_index({i, j, k});
                        for (int d = 0; d < 3; ++d) {
                            electric(raw.x(), raw.y(), raw.z(), d) =
                                constant ? d + 2.0 : affine(yee_position(grid, {i, j, k}, d, false), d);
                            magnetic(raw.x(), raw.y(), raw.z(), d) =
                                constant ? d + 2.0 : affine(yee_position(grid, {i, j, k}, d, true), d);
                        }
                    }
            for (const Vector3R coord : {Vector3R(0.0), Vector3R(0.5), Vector3R(1.13, 2.72, 3.41),
                                         Vector3R(4.9, 5.8, 6.7), Vector3R(-0.2, 0.1, -0.1)}) {
                if (ghost == 1 && coord.x() < 0)
                    continue;
                const auto no = make_point_stencil<ST>(coord);
                const auto sh = make_point_stencil<ST>(coord, Vector3R(0.5));
                Vector3R expected;
                const Vector3R physical(coord.x() * 0.5, coord.y() * 1.25, coord.z() * 2);
                for (int d = 0; d < 3; ++d) expected[d] = constant ? d + 2.0 : affine(physical, d);
                vector_near(interpolateE(electric, no, sh), expected, "electric affine reproduction");
                vector_near(interpolateB(magnetic, no, sh), expected, "magnetic affine reproduction");
            }
        }
    }
}

template <ShapeType ST>
void point_current_work() {
    constexpr int ghost = GHOST_CELLS;
    constexpr int N = ST == ShapeType::Linear ? 2 : 3;
    {
        const Grid grid(Vector3R(0.5, 1.25, 2), Vector3I(5, 6, 7), ghost);
        Field3d field(grid.size(), 3), current(grid.size(), 3);
        for (int i = 0; i < grid.size().x(); ++i)
            for (int j = 0; j < grid.size().y(); ++j)
                for (int k = 0; k < grid.size().z(); ++k)
                    for (int d = 0; d < 3; ++d) field(i, j, k, d) = std::sin(0.3 * i - 0.4 * j + 0.17 * k + d);
        for (const Vector3R coord :
             {Vector3R(0.0), Vector3R(0.5), Vector3R(1.13, 2.72, 3.41), Vector3R(5, 6, 7), Vector3R(-0.2, 0.1, -0.1)}) {
            if (ghost == 1 && coord.x() < 0)
                continue;
            current.setZero();
            const auto no = make_point_stencil<ST>(coord);
            const auto sh = make_point_stencil<ST>(coord, Vector3R(0.5));
            Vector3I origin;
            for (int a = 0; a < 3; ++a) origin[a] = std::min(no.origin[a], sh.origin[a]);
            CurrentBuffer<N + 1> buffer;
            buffer.zero();
            const Vector3R value(0.7, -1.2, 2.3);
            // Accumulating twice also checks that deposition adds to the buffer.
            decompose_point_current(no, sh, value * 0.25, buffer, origin);
            decompose_point_current(no, sh, value * 0.75, buffer, origin);
            flush_point_current(current, buffer, origin, grid);
            Vector3R total(0.0);
            double grid_work = 0;
            for (int i = 0; i < grid.size().x(); ++i)
                for (int j = 0; j < grid.size().y(); ++j)
                    for (int k = 0; k < grid.size().z(); ++k)
                        for (int d = 0; d < 3; ++d) {
                            total[d] += current(i, j, k, d);
                            grid_work += current(i, j, k, d) * field(i, j, k, d);
                        }
            vector_near(total, value, "total deposited current");
            const auto gathered = interpolateE(field, no, sh);
            double particle_work = 0;
            for (int d = 0; d < 3; ++d) particle_work += value[d] * gathered[d];
            check::near(grid_work, particle_work, tolerance, "gather/scatter work equality");
        }
    }
}

void invalid_stencil_coordinates_are_rejected() {
    check::equal(checked_floor_index(double(std::numeric_limits<int>::max()) + 1, -1), std::numeric_limits<int>::max(),
                 "logical offset is applied before narrowing");
    bool underflow_rejected = false;
    try {
        (void) checked_floor_index(double(std::numeric_limits<int>::min()), -1);
    } catch (const std::out_of_range&) {
        underflow_rejected = true;
    }
    check::require(underflow_rejected, "logical stencil offset cannot underflow int");

    for (double invalid : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), 1e300}) {
        bool rejected = false;
        try {
            (void) make_point_stencil<ShapeType::Linear>(Vector3R(invalid, 0, 0));
        } catch (const std::out_of_range&) {
            rejected = true;
        }
        check::require(rejected, "invalid point coordinate rejected before index conversion");
    }
    bool rejected = false;
    try {
        (void) make_stencil_at<Shape, 4>(Vector3R(0.0), Vector3I(std::numeric_limits<int>::max(), 0, 0));
    } catch (const std::out_of_range&) {
        rejected = true;
    }
    check::require(rejected, "explicit stencil support cannot overflow its last node");
}

void legacy_linear_agreement() {
    const Grid grid(Vector3R(1), Vector3I(5, 6, 7), GHOST_CELLS);
    Field3d field(grid.size(), 3);
    for (int i = 0; i < grid.size().x(); ++i)
        for (int j = 0; j < grid.size().y(); ++j)
            for (int k = 0; k < grid.size().z(); ++k)
                for (int d = 0; d < 3; ++d) field(i, j, k, d) = std::cos(0.2 * i + 0.3 * j - 0.1 * k + d);
    for (const Vector3R coord : {Vector3R(0.0), Vector3R(0.5), Vector3R(1.13, 2.72, 3.41), Vector3R(4.9, 5.8, 6.7)}) {
        vector_near(interpolateE(field, make_point_stencil<ShapeType::Linear>(coord),
                                 make_point_stencil<ShapeType::Linear>(coord, Vector3R(0.5))),
                    raw_gather_reference::interpolateE_linear(field, coord), "legacy linear electric gather");
        vector_near(interpolateB(field, make_point_stencil<ShapeType::Linear>(coord),
                                 make_point_stencil<ShapeType::Linear>(coord, Vector3R(0.5))),
                    raw_gather_reference::interpolateB_linear(field, coord), "legacy linear magnetic gather");
    }
}
}   // namespace

int main() {
    return check::run([] {
        check::scenario("Invalid stencil coordinates and integer overflow", invalid_stencil_coordinates_are_rejected);
        check::scenario("Linear shape preserves unity and first moment across shifted origin",
                        shape_moments<ShapeType::Linear>);
        check::scenario("Quadratic shape preserves unity and first moment across shifted origin",
                        shape_moments<ShapeType::Quadratic>);
        check::scenario("Linear gathers reproduce constant and affine Yee fields",
                        [] { gather_fields<ShapeType::Linear>(); });
        check::scenario("Quadratic gathers reproduce constant and affine Yee fields",
                        [] { gather_fields<ShapeType::Quadratic>(); });
        check::scenario("Linear current conserves sum and grid-particle work",
                        [] { point_current_work<ShapeType::Linear>(); });
        check::scenario("Quadratic current conserves sum and grid-particle work",
                        [] { point_current_work<ShapeType::Quadratic>(); });
        check::scenario("Logical Linear gather agrees with legacy gather", legacy_linear_agreement);
    });
}
