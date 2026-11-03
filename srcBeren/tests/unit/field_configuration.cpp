#include <algorithm>
#include <cmath>
#include <numbers>

#include "Coil.h"
#include "config/field_config.h"
#include "test_case.h"

namespace {
void external_fields_add_the_configured_contributions() {
    Domain domain;
    domain.init(Vector3I(6, 6, 6), Vector3R(0.5, 0.75, 0.25));
    const nlohmann::json uniform = {{"uniform_field", {{"value", {1, 2, 3}}}}};
    const nlohmann::json cylinder = {{"uniformly_charged_cylinder", {{"radius", 0.7}, {"value", 4}}}};
    for (bool as_array : {false, true}) {
        nlohmann::json config;
        if (as_array)
            config["ExternalFieldE"] = nlohmann::json::array({uniform, cylinder});
        else
            config["ExternalFieldE"] = {{"uniform_field", uniform.at("uniform_field")},
                                        {"uniformly_charged_cylinder", cylinder.at("uniformly_charged_cylinder")}};
        Field3d field(domain.size(), 3);
        field.data().setConstant(0.5);
        create_electric_field_config(config)->apply(field, domain);
        // Ex nodes: r=0.25 inside the cylinder and r=0.75 outside it.
        // Er = value*r/2 inside; Er = value*R^2/(2*r) outside.
        check::near(field(3, 4, 2, 0), 1.5 - 2 * 0.25, 1e-13, "cylinder inner Ex plus uniform and initial E");
        check::near(field(5, 4, 2, 0), 1.5 + 2 * 0.7 * 0.7 / 0.75, 1e-13,
                    "cylinder outer Ex plus uniform and initial E");
        check::near(field(3, 4, 2, 2), 3.5, 0, "cylinder preserves the uniform axial E");
    }
    Field3d magnetic(domain.size(), 3);
    magnetic.data().setConstant(0.5);
    create_magnetic_field_config({{"ExternalFieldB", nlohmann::json::array({uniform, uniform})}})
        ->apply(magnetic, domain);
    for (int d = 0; d < 3; ++d)
        check::near(magnetic(2, 2, 2, d), 0.5 + 2 * (d + 1), 0, "uniform B contributions add to the initial field");
    check::require(!create_electric_field_config(nlohmann::json::object()), "absent external E is disabled");
    check::require(!create_magnetic_field_config({{"ExternalFieldB", nlohmann::json::array({3, "unknown"})}}),
                   "unrecognized magnetic array entries remain ignored");
    Field3d zero_contribution(domain.size(), 3);
    zero_contribution.data().setConstant(0.5);
    create_electric_field_config({{"ExternalFieldE", {{"uniform_field", {{"value", {1, 2}}}}}}})
        ->apply(zero_contribution, domain);
    check::near((zero_contribution.data().array() - 0.5).matrix().norm(), 0, 0,
                "legacy incomplete uniform vector contributes zero");
}

void coil_parameters_reproduce_the_axis_field_and_superposition() {
    const std::vector<Coil> parameters{{0.5, 2, 0.3}, {-1, 1.5, -0.1}};
    CoilsArray coils(parameters);
    for (double z : {-2.0, 0.5, 3.0}) {
        double expected = 0;
        for (const Coil& coil : parameters) {
            const double distance = z - coil.z0;
            expected +=
                2 * std::numbers::pi * coil.I * coil.R * coil.R / std::pow(coil.R * coil.R + distance * distance, 1.5);
        }
        check::near(coils.get_Bz(z, 0), expected, 1e-12, "circular-loop axial B from Biot-Savart");
        check::near(coils.get_Br(z, 0), 0, 1e-12, "axial symmetry cancels radial B");
    }
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(0.5, 0.5, 0.5));
    const nlohmann::json coil_json =
        nlohmann::json::array({{{"z", 0.5}, {"R", 2}, {"I", 0.3}}, {{"z", -1}, {"R", 1.5}, {"I", -0.1}}});
    Field3d expected(domain.size(), 3), actual(domain.size(), 3);
    expected.data().setConstant(0.5);
    actual.data().setConstant(0.5);
    set_coils(expected, domain, parameters);
    const nlohmann::json config = {
        {"ExternalFieldB",
         nlohmann::json::array({{{"coils", coil_json}}, {{"uniform_field", {{"value", {1, 1, 1}}}}}})}};
    create_magnetic_field_config(config)->apply(actual, domain);
    check::require(actual.data().allFinite() && expected.data().allFinite(),
                   "coil fields must be finite at every node");
    check::near((actual.data() - expected.data()).minCoeff(), 1, 1e-12, "coils add to initial and uniform B");
    check::near((actual.data() - expected.data()).maxCoeff(), 1, 1e-12, "all magnetic nodes retain superposition");
    bool rejected = false;
    try {
        create_magnetic_field_config({{"ExternalFieldB", {{"coils", nlohmann::json::array({{{"z", 0}, {"R", 2}}})}}}});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check::require(rejected, "missing coil current must be rejected before applying a field");
}

void damping_matches_the_spatial_profile_and_removed_energy() {
    Domain domain;
    domain.init(Vector3I(6, 6, 6), Vector3R(0.5, 0.75, 0.25));
    for (const std::string type : {"None", "Rectangle", "CircleXY"}) {
        nlohmann::json config = {{"DampingType", type},
                                 {"DampingCoeff", 0.8},
                                 {"DampCellsX_glob", {2, 1}},
                                 {"DampCellsY_glob", {1, 2}},
                                 {"DampCellsZ_glob", {2, 2}}};
        const DampingParameters parameters = damping_parameters_from_config(config);
        // A parsed record owns its values; subsequent JSON edits cannot change the layer.
        config["DampingCoeff"] = 0;
        Field3d e(domain.size(), 3), b(domain.size(), 3);
        e.data().setConstant(1);
        b.data().setConstant(-2);
        const double removed = damping_fields(e, b, domain, parameters);
        double expected_removed = 0;
        for (int i = 0; i < domain.size().x(); ++i)
            for (int j = 0; j < domain.size().y(); ++j)
                for (int k = 0; k < domain.size().z(); ++k) {
                    double factor = 1;
                    if (type == "Rectangle") {
                        const Vector3I index(i, j, k);
                        const std::array<std::array<int, 2>, 3> widths{{{2, 1}, {1, 2}, {2, 2}}};
                        for (int axis = 0; axis < 3; ++axis) {
                            for (int side = 0; side < 2; ++side) {
                                const int distance = side == 0 ? index[axis] : domain.size()[axis] - 1 - index[axis];
                                const int width = widths[axis][side];
                                if (distance < width)
                                    factor *= 0.8 + 0.2 * std::pow(double(distance) / width, 2);
                            }
                        }
                    } else if (type == "CircleXY") {
                        // Preserve the circular layer's existing radial convention: X spacing on both axes.
                        const int center = (domain.size().x() - 1) / 2;
                        const double r = 0.5 * std::hypot(i - center, j - center);
                        const double outer = 0.5 * 0.5 * (domain.size().x() - 1);
                        const double width = 2 * 0.5;
                        if (r > outer - width)
                            factor = r >= outer ? 0.8 : 0.8 + 0.2 * std::pow((outer - r) / width, 2);
                    }
                    for (int d = 0; d < 3; ++d) {
                        check::near(e(i, j, k, d), factor, 1e-14, "E follows the quadratic damping profile");
                        check::near(b(i, j, k, d), -2 * factor, 1e-14, "B follows the same damping profile");
                    }
                    expected_removed += 0.5 * 3 * (1 + 4) * (1 - factor * factor);
                }
        check::near(removed, expected_removed, 1e-10, "damping returns the independently integrated energy loss");
    }
    check::require(damping_parameters_from_config({{"DampingCoeff", "ignored"}}).type == DampingType::None,
                   "disabled damping ignores unused parameters");
    const DampingParameters circle = damping_parameters_from_config({{"DampingType", "CircleXY"}});
    check::equal(circle.x_cells[0], 0, "missing circular width defaults to zero");
    check::equal(circle.x_cells[1], 0, "missing circular width defaults to zero on both sides");
    bool rejected = false;
    try {
        damping_parameters_from_config({{"DampingType", "unknown"}});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check::require(rejected, "unknown damping type is rejected during initialization");
}
}   // namespace

int main() {
    return check::run([] {
        check::scenario("External E/B: additive JSON composition and analytic cylinder field",
                        external_fields_add_the_configured_contributions);
        check::scenario("Coils: Biot-Savart axis field, parsed superposition and incomplete input",
                        coil_parameters_reproduce_the_axis_field_and_superposition);
        check::scenario("Damping: spatial profiles and energy loss from owned parameters",
                        damping_matches_the_spatial_profile_and_removed_energy);
    });
}
