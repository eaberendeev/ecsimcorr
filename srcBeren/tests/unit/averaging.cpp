#include <filesystem>
#include <fstream>

#include "DiagnosticOutput.h"
#include "output_util.h"
#include "test_support.h"

namespace {
void verify_plane(const std::string& filename, const Field3d& field, int axis, int position, double offset) {
    std::ifstream file(filename, std::ios::binary);
    check::require(file.is_open(), "missing averaged plane: " + filename);
    // The binary slice format stores both dimensions and samples as float.
    float widthRaw = 0, heightRaw = 0;
    file.read(reinterpret_cast<char*>(&widthRaw), sizeof(float));
    file.read(reinterpret_cast<char*>(&heightRaw), sizeof(float));
    check::require(bool(file) && std::isfinite(widthRaw) && std::isfinite(heightRaw), "invalid plane header");
    const int n1 = static_cast<int>(widthRaw);
    const int n2 = static_cast<int>(heightRaw);
    const int first = axis == 0 ? 1 : 0;
    const int second = axis == 2 ? 1 : 2;
    check::near(n1, field.sizes()[first], 0, "plane width");
    check::near(n2, field.sizes()[second], 0, "plane height");
    for (int d = 0; d < 3; ++d)
        for (int i = 0; i < n1; ++i)
            for (int j = 0; j < n2; ++j) {
                float stored = 0;
                file.read(reinterpret_cast<char*>(&stored), sizeof(float));
                check::require(bool(file), "truncated averaged plane");
                const double actual = stored;
                Vector3I node(0, 0, 0);
                node[axis] = position;
                node[first] = i;
                node[second] = j;
                check::near(actual, field(node.x(), node.y(), node.z(), d) + offset, 1e-6,
                            "arithmetic mean at plane node");
            }
    check::require(file.peek() == std::char_traits<char>::eof(), "extra data in averaged plane");
    file.close();
    std::filesystem::remove(filename);
}

void selected_planes_have_correct_mean_and_storage() {
    Domain domain;
    domain.init(Vector3I(4, 6, 8), Vector3R(0.5, 1, 2));
    Field3d field(domain.size(), 3), extra(domain.size(), 3);
    for (int i = 0; i < field.sizes().x(); ++i)
        for (int j = 0; j < field.sizes().y(); ++j)
            for (int k = 0; k < field.sizes().z(); ++k)
                for (int d = 0; d < 3; ++d) field(i, j, k, d) = 1000 * d + 100 * i + 10 * j + k;
    extra.data().setConstant(2);
    AveragedField2D avg;
    avg.set_planes({0.5}, {2}, {6}, domain);   // Storage indices 1, 2, 3; same convention as write_slices.
    check::near(avg.stored_values(), 0, 0, "empty window has no field storage");
    check::require(avg.count == 0, "empty window is not ready");
    avg.add("E", field);
    avg.add("B", field, &extra);
    ++avg.count;
    extra.data().setConstant(4);
    avg.add("E", field, &extra);
    avg.add("B", field, &extra);
    ++avg.count;
    check::require(avg.count == 2, "window readiness counts samples");
    const auto n = field.sizes();
    const size_t expected = 2 * 3 * (n.y() * n.z() + n.x() * n.z() + n.x() * n.y());
    check::near(avg.stored_values(), expected, 0, "storage scales with requested plane areas");
    for (const auto& key : {"E", "B"}) {
        avg.write(key, key, 0);
        for (int axis = 0; axis < 3; ++axis) {
            const int pos = axis + 1;
            verify_plane(std::string(key) + "_Plane" + "XYZ"[axis] + "_00" + std::to_string(pos) + "_avg0000", field,
                         axis, pos, key == std::string("E") ? 2 : 3);
        }
    }
    avg.reset();
    check::equal(avg.count, 0, "reset clears sample count");
    check::near(avg.sequence(), 1, 0, "reset advances file sequence");
    avg.add("E", field);
    ++avg.count;
    avg.write("E", "partial", 1);
    for (int axis = 0; axis < 3; ++axis)
        verify_plane(std::string("partial_Plane") + "XYZ"[axis] + "_00" + std::to_string(axis + 1) + "_avg0001", field,
                     axis, axis + 1, 0);
    avg.reset();
    check::near(avg.stored_values(), expected, 0, "storage is reused between windows");
}
}   // namespace

int main() {
    return check::run([] {
        check::scenario("Selected X/Y/Z planes: arithmetic means, residual window and storage size",
                        selected_planes_have_correct_mean_and_storage);
    });
}
