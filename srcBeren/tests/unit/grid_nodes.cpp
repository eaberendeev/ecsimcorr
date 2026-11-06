#include <type_traits>
#include <utility>
#ifndef NDEBUG
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

#include <csignal>
#endif

#include "grid/grid_node_view.h"
#include "test_case.h"
#include "utils/containers.h"

using Scalar = Array3D<double>;
using Vector = Field3dBase<double>;
static_assert(!std::is_constructible_v<GridNodeView<Scalar>, Scalar&&, const Grid&>);
static_assert(!std::is_constructible_v<GridNodeView<const Scalar>, const Scalar&&, const Grid&>);
static_assert(std::is_same_v<decltype(std::declval<GridNodeView<const Scalar>>()(0, 0, 0)), const double&>);
static_assert(std::is_same_v<decltype(std::declval<GridNodeView<const Vector>>()(0, 0, 0, 0)), const double&>);

void logical_nodes_and_storage() {
    for (int ghost : {1, 2}) {
        const Grid grid(Vector3R(0.5, 1, 2), Vector3I(3, 4, 5), ghost);
        Scalar scalar(grid.size());
        Vector vector(grid.size(), 3);
        const GridNodeView scalar_nodes(scalar, grid);
        const GridNodeView vector_nodes(vector, grid);
        const GridNodeView scalar_const(std::as_const(scalar), grid);
        const GridNodeView vector_const(std::as_const(vector), grid);
        for (int i = -ghost; i <= 3 + ghost; ++i)
            for (int j = -ghost; j <= 4 + ghost; ++j)
                for (int k = -ghost; k <= 5 + ghost; ++k) {
                    const Vector3I node(i, j, k);
                    const Vector3I raw(i + ghost, j + ghost, k + ghost);
                    check::require(grid.to_storage_index(node) == raw, "logical to storage offset");
                    check::require(grid.to_node_index(raw) == node, "storage to logical offset");
                    scalar_nodes(i, j, k) = grid.sind(raw.x(), raw.y(), raw.z());
                    check::require(&scalar_nodes(i, j, k) == &scalar(raw.x(), raw.y(), raw.z()), "scalar alias");
                    check::equal(scalar_const(i, j, k), double(grid.sind(raw.x(), raw.y(), raw.z())), 0,
                                 "const scalar value");
                    for (int d = 0; d < 3; ++d) {
                        vector_nodes(i, j, k, d) = grid.vind(raw.x(), raw.y(), raw.z(), d);
                        check::require(&vector_const(i, j, k, d) == &vector(raw.x(), raw.y(), raw.z(), d),
                                       "component alias");
                        for (FieldType field :
                             {FieldType::ELECTRIC, FieldType::MAGNETIC, FieldType::CURRENT, FieldType::DENSITY}) {
                            // Independent Yee coordinates: E/J have one half-step;
                            // B has the other two, density has none.
                            Vector3R expected(i * 0.5, j, k * 2.0);
                            for (int axis = 0; axis < 3; ++axis) {
                                const bool shifted =
                                    field == FieldType::MAGNETIC ? axis != d : field != FieldType::DENSITY && axis == d;
                                if (shifted)
                                    expected[axis] += 0.5 * grid.cell_size(axis);
                            }
                            const auto position = grid.node_position(node, field, d);
                            const auto legacy = grid.get_node_position(raw.x(), raw.y(), raw.z(), field, d);
                            for (int axis = 0; axis < 3; ++axis) {
                                check::equal(position[axis], expected[axis], 0, "logical Yee position");
                                check::equal(legacy[axis], expected[axis], 0, "legacy stored Yee position");
                            }
                            check::require(grid.field_node_index(expected, field, d) == node, "logical Yee inverse");
                            check::require(grid.get_field_node_index(expected, field, d) == raw, "legacy Yee inverse");
                            const auto before = expected - Vector3R(0.125, 0.25, 0.5);
                            check::require(grid.field_node_index(before, field, d) == node - Vector3I(1),
                                           "floor below node, including negative coordinates");
                        }
                    }
                }
        // Legacy cell lookup truncates after the ghost offset, unlike floor.
        const auto cell = grid.get_cell_index(Vector3R((-ghost - 0.25) * 0.5, 0, 0));
        check::equal(cell.x(), 0, "legacy cell truncation outside ghost region");
    }
}

#ifndef NDEBUG
template <typename Action>
void requires_assertion(Action action) {
    const pid_t child = fork();
    check::require(child >= 0, "fork for debug assertion test");
    if (child == 0) {
        const rlimit no_core{0, 0};
        setrlimit(RLIMIT_CORE, &no_core);
        action();
        _exit(0);
    }
    int status = 0;
    check::require(waitpid(child, &status, 0) == child, "wait for debug assertion");
    check::require(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT, "invalid access must assert");
}
void invalid_view_access() {
    const Grid grid(Vector3R(1), Vector3I(3, 4, 5), 2);
    Scalar scalar(grid.size());
    Vector vector(grid.size(), 3);
    GridNodeView scalar_nodes(scalar, grid);
    GridNodeView vector_nodes(vector, grid);
    for (int axis = 0; axis < 3; ++axis) {
        for (bool upper : {false, true}) {
            Vector3I node(0);
            node[axis] = upper ? grid.num_cells(axis) + 3 : -3;
            requires_assertion([&] { scalar_nodes(node.x(), node.y(), node.z()); });
        }
    }
    requires_assertion([&] { vector_nodes(0, 0, 0, -1); });
    requires_assertion([&] { vector_nodes(0, 0, 0, 3); });
    requires_assertion([&] {
        Scalar wrong(Vector3I(2));
        GridNodeView invalid(wrong, grid);
    });
    requires_assertion([&] {
        Vector wrong(Vector3I(2), 3);
        GridNodeView invalid(wrong, grid);
    });
    requires_assertion([&] {
        scalar.resize(Vector3I(2));
        scalar_nodes(0, 0, 0);
    });
}
#endif

int main() {
    return check::run([] {
        check::scenario("Logical nodes alias raw scalar/vector storage for ghost widths 1 and 2",
                        logical_nodes_and_storage);
#ifndef NDEBUG
        check::scenario("Invalid dimensions, components and logical nodes fail in Debug", invalid_view_access);
#endif
    });
}
