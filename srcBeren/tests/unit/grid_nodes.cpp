#include <type_traits>
#include <utility>

#include "assert_contract.h"
#include "grid/grid.h"
#include "grid/grid_node_view.h"
#include "test_case.h"
#include "utils/containers.h"

using Scalar = Array3D<double>;
using Vector = Field3dBase<double>;
static_assert(!std::is_constructible_v<GridNodeView<Scalar>, Scalar&&>);
static_assert(!std::is_constructible_v<GridNodeView<const Scalar>, const Scalar&&>);
static_assert(std::is_convertible_v<Scalar&, GridNodeView<Scalar>>);
static_assert(std::is_convertible_v<const Scalar&, GridNodeView<const Scalar>>);
static_assert(std::is_same_v<decltype(std::declval<GridNodeView<const Scalar>>()(0, 0, 0)), const double&>);
static_assert(std::is_same_v<decltype(std::declval<GridNodeView<const Vector>>()(0, 0, 0, 0)), const double&>);

template <int ghost>
void logical_nodes_and_storage() {
    {
        const Grid grid(Vector3R(0.5, 1, 2), Vector3I(3, 4, 5), ghost);
        Scalar scalar(grid.size());
        Vector vector(grid.size(), 3);
        const GridNodeView<Scalar, ghost> scalar_nodes(scalar);
        const GridNodeView<Vector, ghost> vector_nodes(vector);
        const GridNodeView<const Scalar, ghost> scalar_const(std::as_const(scalar));
        const GridNodeView<const Vector, ghost> vector_const(std::as_const(vector));
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
                            const auto legacy = grid.storage_node_position(raw.x(), raw.y(), raw.z(), field, d);
                            for (int axis = 0; axis < 3; ++axis) {
                                check::equal(position[axis], expected[axis], 0, "logical Yee position");
                                check::equal(legacy[axis], expected[axis], 0, "legacy stored Yee position");
                            }
                            check::require(grid.field_node_index(expected, field, d) == node, "logical Yee inverse");
                            check::require(grid.storage_field_node_index(expected, field, d) == raw,
                                           "legacy Yee inverse");
                            const auto before = expected - Vector3R(0.125, 0.25, 0.5);
                            check::require(grid.field_node_index(before, field, d) == node - Vector3I(1),
                                           "floor below node, including negative coordinates");
                        }
                    }
                }
        // Particle bins use the same logical floor before the storage offset.
        const auto cell = grid.storage_cell_index(Vector3R((-ghost - 0.25) * 0.5, 0, 0));
        check::equal(cell.x(), -1, "cell below lower halo must not alias storage zero");
    }
}

#ifndef NDEBUG
void invalid_view_access() {
    const Grid grid(Vector3R(1), Vector3I(3, 4, 5), 2);
    Scalar scalar(grid.size());
    Vector vector(grid.size(), 3);
    GridNodeView<Scalar, 2> scalar_nodes(scalar);
    GridNodeView<Vector, 2> vector_nodes(vector);
    for (int axis = 0; axis < 3; ++axis) {
        for (bool upper : {false, true}) {
            Vector3I node(0);
            node[axis] = upper ? grid.num_cells(axis) + 3 : -3;
            check::assertion_fails([&] { scalar_nodes(node.x(), node.y(), node.z()); });
        }
    }
    check::assertion_fails([&] { vector_nodes(0, 0, 0, -1); });
    check::assertion_fails([&] { vector_nodes(0, 0, 0, 3); });
    check::assertion_fails([&] { grid.cell_size(Axis::Comp); });
    check::assertion_fails([&] { grid.num_cells(Axis::Comp); });
    check::assertion_fails([&] { grid.size(Axis::Comp); });
    check::assertion_fails([&] { grid.pos_sind(0, Axis::Comp); });
    check::assertion_fails([&] {
        scalar.resize(Vector3I(2));
        scalar_nodes(0, 0, 0);
    });
}
#endif

int main() {
    return check::run([] {
        check::scenario("Logical nodes alias raw scalar/vector storage for ghost widths 1 and 2", [] {
            logical_nodes_and_storage<1>();
            logical_nodes_and_storage<2>();
        });
#ifndef NDEBUG
        check::scenario("Invalid components and logical nodes fail in Debug", invalid_view_access);
#endif
    });
}
