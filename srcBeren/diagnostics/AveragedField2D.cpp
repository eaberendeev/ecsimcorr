#include "AveragedField2D.h"

#include <algorithm>

#include "outputs/SliceFiles.h"
#include "timer.h"
#include "util.h"

// ============================================================
//  AveragedField2D
// ============================================================
void AveragedField2D::set_planes(const std::vector<double>& px, const std::vector<double>& py,
                                 const std::vector<double>& pz, const Domain& domain) {
    positions_ = slice_indices(px, py, pz, domain);
}

void AveragedField2D::add(const std::string& key, const Field3d& field, const Field3d* extra) {
    RECORD_TIMER;
    if (field.capacity() == 0)
        return;
    auto [it, inserted] = accum_.try_emplace(key);
    auto& planes = it->second;
    if (inserted)
        planes = make_field_slices(field.sizes(), field.nd(), positions_);
    add_field_to_slices(planes, field, extra);
}

void AveragedField2D::add(const std::string& key, const FieldSlices& samples) {
    auto [it, inserted] = accum_.try_emplace(key);
    auto& planes = it->second;
    const auto find_sample = [&](int axis, int position) {
        return std::find_if(samples.begin(), samples.end(),
                            [&](const auto& sample) { return sample.axis == axis && sample.position == position; });
    };
    if (inserted)
        for (int axis = 0; axis < 3; ++axis)
            for (int position : positions_[axis]) {
                const auto sample = find_sample(axis, position);
                if (sample == samples.end())
                    continue;
                planes.push_back({axis, position, Field3d(sample->values.sizes(), sample->values.nd())});
                planes.back().values.setZero();
            }
    for (auto& plane : planes) {
        const auto sample = find_sample(plane.axis, plane.position);
        if (sample != samples.end())
            plane.values.data() += sample->values.data();
    }
}

void AveragedField2D::write(const std::string& key, const std::string& basename, int sequence, bool averaged) const {
    RECORD_TIMER;
    const auto it = accum_.find(key);
    if (it == accum_.end() || count == 0)
        return;
    const std::string suffix = std::string(averaged ? "avg" : "") + to_string(sequence, 4);
    for (const auto& plane : it->second) {
        // Preserve the binary layout: three component slots, even for scalar
        // density/pressure files. Unused slots remain zero.
        const auto sizes = plane.values.sizes();
        const int first = plane.axis == 0 ? 1 : 0;
        const int second = plane.axis == 2 ? 1 : 2;
        const int n1 = sizes[first], n2 = sizes[second];
        std::vector<float> values(3 * n1 * n2, 0.0f);
        for (int i = 0; i < n1; ++i)
            for (int j = 0; j < n2; ++j) {
                Vector3I node(0, 0, 0);
                node[first] = i;
                node[second] = j;
                for (int d = 0; d < plane.values.nd(); ++d)
                    values[output_field_index2d(i, j, d, n1, n2)] =
                        static_cast<float>(plane.values(node.x(), node.y(), node.z(), d) / count);
            }
        const std::string axis(1, "XYZ"[plane.axis]);
        output_array2d(values, n1, n2, basename + "_Plane" + axis + "_" + to_string(plane.position, 3) + "_" + suffix);
    }
}

size_t AveragedField2D::stored_values() const {
    size_t result = 0;
    for (const auto& [key, planes] : accum_)
        for (const auto& plane : planes) result += plane.values.capacity();
    return result;
}

std::vector<std::string> AveragedField2D::keys() const {
    std::vector<std::string> result;
    for (const auto& entry : accum_) result.push_back(entry.first);
    return result;
}

long long AveragedField2D::sequence() const {
    return sequence_;
}

void AveragedField2D::reset() {
    RECORD_TIMER;
    for (auto& [key, planes] : accum_)
        for (auto& plane : planes) plane.values.setZero();
    count = 0;
    ++sequence_;
}
