#pragma once

#include <cassert>

#include "grid/grid.h"

// A non-owning logical-node interface to raw storage. It only adds the ghost
// offset: no Yee shift, periodic wrapping, clipping or atomic update is implied.
// Storage must outlive the view. Recreate it after resize/reallocation or a grid
// change; it is a transient access helper and is not part of restart state.
template <typename Storage>
class GridNodeView {
   public:
    GridNodeView(Storage& storage, const Grid& grid) : storage_(&storage), grid_(grid) {
        check_dimensions();
    }
    GridNodeView(Storage&&, const Grid&) = delete;

    decltype(auto) operator()(int i, int j, int k) const
        requires requires(Storage& s) { s(i, j, k); }
    {
        const Vector3I raw = storage_index(i, j, k);
        return (*storage_)(raw.x(), raw.y(), raw.z());
    }
    decltype(auto) operator()(int i, int j, int k, int component) const
        requires requires(Storage& s) { s(i, j, k, component); }
    {
        const Vector3I raw = storage_index(i, j, k);
        assert(component >= 0 && component < storage_->nd());
        return (*storage_)(raw.x(), raw.y(), raw.z(), component);
    }

   private:
    void check_dimensions() const {
#ifndef NDEBUG
        // Scalar Array3D exposes size(); multi-component Field3dBase uses sizes().
        if constexpr (requires { storage_->sizes(); }) {
            assert(storage_->sizes() == grid_.size());
        } else {
            assert(storage_->size() == grid_.size());
        }
#endif
    }
    Vector3I storage_index(int i, int j, int k) const {
        check_dimensions();
        const Vector3I node(i, j, k);
        for (int axis = 0; axis < 3; ++axis) {
            assert(node[axis] >= -grid_.ghost_cells());
            assert(node[axis] <= grid_.num_cells(axis) + grid_.ghost_cells());
        }
        return grid_.to_storage_index(node);
    }

    Storage* storage_;
    Grid grid_;
};
