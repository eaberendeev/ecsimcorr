// Unit tests for the shared diagnostic 2D averaging helper (AveragedField2D).
// The helper is used by Fields2DOutput and Species2DOutput; these tests cover
// window accumulation, the readiness threshold, reset/sequence bookkeeping and
// the partial-window mean used by finalize().
#include <cmath>
#include <iostream>
#include <string>

#include "DiagnosticOutput.h"
#include "containers.h"

namespace test {

static int passed = 0;
static int failed = 0;

void assert_true(bool condition, const std::string& test_name) {
    if (condition) {
        passed++;
    } else {
        std::cout << "[FAIL] " << test_name << std::endl;
        failed++;
    }
}

void assert_eq(double actual, double expected, double eps, const std::string& test_name) {
    if (std::abs(actual - expected) < eps) {
        passed++;
    } else {
        std::cout << "[FAIL] " << test_name << " (expected: " << expected << ", got: " << actual << ")" << std::endl;
        failed++;
    }
}

void assert_eq(int actual, int expected, const std::string& test_name) {
    if (actual == expected) {
        passed++;
    } else {
        std::cout << "[FAIL] " << test_name << " (expected: " << expected << ", got: " << actual << ")" << std::endl;
        failed++;
    }
}

void print_summary() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test Summary: " << passed << " passed, " << failed << " failed" << std::endl;
    std::cout << "========================================\n" << std::endl;
}

int get_failed() {
    return failed;
}

}   // namespace test

namespace {

Field3d constant_field(double value) {
    Field3d field(Vector3I(2, 2, 2), 1);
    field.data().setConstant(value);
    return field;
}

void assert_field_value(const Field3d& field, double expected, const std::string& test_name) {
    bool all_match = field.capacity() > 0;
    for (int i = 0; i < field.capacity(); ++i) {
        if (!std::isfinite(field(i)) || std::abs(field(i) - expected) >= 1e-15) {
            all_match = false;
            break;
        }
    }
    test::assert_true(all_match, test_name);
}

void test_empty_helper() {
    AveragedField2D avg;
    test::assert_eq(avg.count(), 0, "empty: count is zero");
    test::assert_true(!avg.ready(1), "empty: not ready for a 1-step window");
    test::assert_true(!avg.ready(0), "empty: average_steps=0 is never ready");
    test::assert_true(avg.mean("E").capacity() == 0, "empty: no mean available");
    test::assert_eq(static_cast<int>(avg.sequence()), 0, "empty: sequence starts at zero");
}

void test_two_field_mean() {
    AveragedField2D avg;
    avg.add("E", constant_field(1.0));
    avg.increment();
    avg.add("E", constant_field(3.0));
    avg.increment();

    test::assert_eq(avg.count(), 2, "mean: count tracks output calls");
    test::assert_true(avg.ready(2), "mean: ready when window is full");
    test::assert_true(!avg.ready(3), "mean: not ready before the window is full");
    assert_field_value(avg.mean("E"), 2.0, "mean: arithmetic mean of two fields");

    // Unknown key must not produce a field.
    test::assert_true(avg.mean("B").capacity() == 0, "mean: unknown key is empty");
}

void test_independent_keys() {
    AveragedField2D avg;
    avg.add("E", constant_field(2.0));
    avg.add("B", constant_field(10.0));
    avg.increment();

    assert_field_value(avg.mean("E"), 2.0, "keys: E mean is independent");
    assert_field_value(avg.mean("B"), 10.0, "keys: B mean is independent");
}

void test_reset_and_sequence() {
    AveragedField2D avg;
    avg.add("E", constant_field(4.0));
    avg.increment();
    avg.increment();

    test::assert_eq(static_cast<int>(avg.sequence()), 0, "reset: initial sequence");
    test::assert_true(avg.ready(2), "reset: window is full before reset");

    avg.reset();
    test::assert_eq(avg.count(), 0, "reset: count cleared");
    test::assert_eq(static_cast<int>(avg.sequence()), 1, "reset: sequence advanced");
    test::assert_true(avg.mean("E").capacity() == 0, "reset: mean unavailable after reset");

    // A new window continues from the next sequence number.
    avg.add("E", constant_field(6.0));
    avg.increment();
    assert_field_value(avg.mean("E"), 6.0, "reset: new window accumulates from zero");
    test::assert_eq(static_cast<int>(avg.sequence()), 1, "reset: sequence kept for next window");
}

void test_partial_window_mean() {
    // Mirrors finalize(): a non-empty window with count < average_steps still
    // exposes a mean.
    AveragedField2D avg;
    avg.add("E", constant_field(5.0));
    avg.increment();

    test::assert_eq(avg.count(), 1, "partial: one accumulated call");
    test::assert_true(!avg.ready(4), "partial: not a full window");
    assert_field_value(avg.mean("E"), 5.0, "partial: mean of the residual window");

    avg.reset();
    test::assert_eq(avg.count(), 0, "partial: reset clears residual window");
}

}   // namespace

int main() {
    test_empty_helper();
    test_two_field_mean();
    test_independent_keys();
    test_reset_and_sequence();
    test_partial_window_mean();

    test::print_summary();

    return test::get_failed();
}
