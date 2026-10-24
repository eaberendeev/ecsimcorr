#include <filesystem>
#include <fstream>
#include <limits>

#include "output_util.h"
#include "test_support.h"

void recovery_fields() {
    const std::string path = "recovery-field-test.bin";
    Field3d expected(Vector3I(2, 3, 4), 3), actual(expected.sizes(), expected.nd());
    for (int i = 0; i < expected.capacity(); ++i) expected[i] = 0.125 * i;
    write_field_to_file(path, expected);
    read_field_from_file(path, actual);
    check::near((actual - expected).norm(), 0, 0, "recovery field round trip");

    const auto reject = [&] {
        actual.data().setConstant(7);
        bool rejected = false;
        try {
            read_field_from_file(path, actual);
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        check::require(rejected, "damaged recovery field must fail");
        for (int i = 0; i < actual.capacity(); ++i)
            check::near(actual[i], 7, 0, "failed recovery must not partially replace field");
    };
    for (auto bytes : {size_t{0}, sizeof(int), 3 * sizeof(int),
                       3 * sizeof(int) + size_t(expected.capacity()) * sizeof(double) - 1}) {
        write_field_to_file(path, expected);
        std::filesystem::resize_file(path, bytes);
        reject();
    }
    Field3d wrong_shape(Vector3I(3, 2, 4), 3);
    wrong_shape.setZero();
    write_field_to_file(path, wrong_shape);
    reject();
    for (double value : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        expected[0] = value;
        write_field_to_file(path, expected);
        reject();
    }
    std::filesystem::remove(path);
    reject();
}

int main() {
    return check::run(recovery_fields);
}
