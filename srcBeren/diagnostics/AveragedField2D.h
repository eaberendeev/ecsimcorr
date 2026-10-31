#pragma once

#include <array>
#include <map>
#include <string>
#include <vector>

#include "FieldSlices.h"
#include "containers.h"
#include "grid/domain.h"

// An averaging window stores only requested planes, including ghost nodes just
// like instantaneous files. All named fields share the output-call count.
class AveragedField2D {
   public:
    void set_planes(const std::vector<double>& px, const std::vector<double>& py, const std::vector<double>& pz,
                    const Domain& domain);
    void add(const std::string& key, const Field3d& field, const Field3d* extra = nullptr);
    void add(const std::string& key, const FieldSlices& samples);
    void set_planes(const SliceIndices& positions) {
        positions_ = positions;
    }
    const SliceIndices& planes() const {
        return positions_;
    }
    void write(const std::string& key, const std::string& basename, int sequence, bool averaged = true) const;
    std::vector<std::string> keys() const;
    long long sequence() const;
    size_t stored_values() const;
    void reset();

    // Samples in the current window; updated by the output owner and cleared by reset().
    int count = 0;

   private:
    SliceIndices positions_;
    std::map<std::string, FieldSlices> accum_;
    long long sequence_ = 0;
};
