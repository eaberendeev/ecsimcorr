#pragma once

#include <cassert>

#include "core/types.h"

// A non-owning logical-node interface to raw storage. The compile-time offsets
// add ghost layers only: no Yee shift, wrapping, clipping or atomic update.
// Storage must outlive the view; this transient helper is not restart state.
template <typename Storage, int OffsetX = GHOST_CELLS, int OffsetY = OffsetX, int OffsetZ = OffsetX>
class GridNodeView {
    static constexpr bool IsScalar = requires(Storage& s) {
        s(0, 0, 0);
        s.size();
    };
    static constexpr bool IsVector = requires(Storage& s) {
        s(0, 0, 0, 0);
        s.sizes();
        s.nd();
    };
    static_assert(IsScalar != IsVector, "GridNodeView requires scalar or component storage");

   public:
    // Intentionally implicit: the target type specifies logical indexing.
    GridNodeView(Storage& storage) noexcept : storage_(storage) {
    }
    GridNodeView(Storage&&) = delete;

    decltype(auto) operator()(int i, int j, int k) const
        requires(IsScalar)
    {
        return storage_(i + OffsetX, j + OffsetY, k + OffsetZ);
    }
    decltype(auto) operator()(int i, int j, int k, int component) const
        requires(IsVector)
    {
        assert(component >= 0 && component < storage_.nd());
        return storage_(i + OffsetX, j + OffsetY, k + OffsetZ, component);
    }
    decltype(auto) operator()(int i, int j, int k, Axis component) const
        requires(IsVector)
    {
        return (*this)(i, j, k, static_cast<int>(component));
    }

   private:
    Storage& storage_;
};
