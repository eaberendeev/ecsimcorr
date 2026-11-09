import math
import json
from pathlib import Path
import struct
import tempfile
import unittest

from validate_physics import (component_errors, configuration, read_field, read_particles,
                              snapshot_binary, conserved_counts, check_continuity, write_payload, read_payload,
                              reference_input, save_reference, compare_reference)


class PhysicsValidationTest(unittest.TestCase):
    def test_detects_local_field_error_even_if_integral_is_unchanged(self):
        errors = component_errors([1.01, 0.99], [1, 1], 1)
        self.assertGreater(errors[0], 1e-3)

    def test_nonfinite_and_missing_values_cannot_pass(self):
        for actual, expected in (([math.nan], [1]), ([1], [math.inf]), ([], []), ([1], [1, 2])):
            with self.subTest(actual=actual), self.assertRaises(ValueError):
                component_errors(actual, expected, 1)

    def test_bitwise_repeatable_reference_has_a_roundoff_floor(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reference"
            arrays = {"FieldE": {"shape": [2], "components": 1, "values": [0., 2.]}}
            save_reference(path, arrays, {"FieldE": [0.]}, "scenario", {}, 5, 3)
            arrays["FieldE"]["values"] = [1e-15, 2 + 1e-15]
            compare_reference(path, arrays, "scenario")

    def test_truncated_recovery_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "data"
            path.write_bytes(struct.pack("=3i", 2, 2, 2) + bytes(8))
            with self.assertRaises(ValueError):
                read_field(path)
            path.write_bytes(struct.pack("=i", 2) + bytes(48))
            with self.assertRaises(ValueError):
                read_particles(path)

    def test_particle_redistribution_order_does_not_change_comparison(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "particles"
            path.write_bytes(struct.pack("=i12d", 2, *([2.] * 6 + [1.] * 6)))
            self.assertEqual(read_particles(path), (2, [1.] * 6 + [2.] * 6))

    def test_executable_snapshot_survives_a_working_tree_rebuild(self):
        with tempfile.TemporaryDirectory() as directory:
            source, target = Path(directory) / "binary", Path(directory) / "snapshot"
            source.write_bytes(b"original executable")
            source.chmod(0o755)
            info = snapshot_binary(source, target)
            source.write_bytes(b"rebuilt executable")
            self.assertEqual(target.read_bytes(), b"original executable")
            self.assertTrue(target.stat().st_mode & 0o111)
            self.assertEqual(info["snapshot"], str(target))

    def test_common_particle_loss_is_detected_without_using_the_other_binary(self):
        with self.assertRaises(ValueError):
            conserved_counts([{"ElectronsParticles": 20.0}, {"ElectronsParticles": 19.0}], ["Electrons"])

    def test_continuity_checks_the_physical_residual(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "continuity.txt"
            path.write_text("Time L2 Linf\n0.1 2e-14 1e-14\n")
            self.assertLess(check_continuity(Path(directory), 1, .1)["max_Linf"], 1e-10)
            path.write_text("Time L2 Linf\n0.1 2e-4 1e-4\n")
            with self.assertRaises(ValueError):
                check_continuity(Path(directory), 1, .1)

    def test_portable_state_payload_checks_integrity(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "field.f64"
            descriptor = write_payload(path, [1., 2., 3., 4.], [2, 2])
            self.assertEqual(read_payload(Path(directory), descriptor), [1., 2., 3., 4.])
            path.write_bytes(path.read_bytes()[:-1] + b"x")
            with self.assertRaises(ValueError):
                read_payload(Path(directory), descriptor)

    def test_reference_uses_noise_with_margin_and_rejects_other_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reference"
            arrays = {"FieldE": {"shape": [2], "components": 1, "values": [1., 2.]}}
            save_reference(path, arrays, {"FieldE": [1e-12]}, "scenario", {}, 5, 3)
            arrays["FieldE"]["values"][0] += 2e-12
            compare_reference(path, arrays, "scenario")
            arrays["FieldE"]["values"][0] += 1e-8
            with self.assertRaises(ValueError):
                compare_reference(path, arrays, "scenario")
            with self.assertRaises(ValueError):
                compare_reference(path, arrays, "another scenario")

    def test_invalid_tolerances_cannot_turn_a_reference_check_into_a_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reference"
            arrays = {"FieldE": {"shape": [1], "components": 1, "values": [1.]}}
            save_reference(path, arrays, {"FieldE": [0.]}, "scenario", {}, 5, 3)
            metadata = json.loads((path / "manifest.json").read_text())
            metadata["arrays"]["FieldE"]["abs_limits"] = [float("nan")]
            (path / "manifest.json").write_text(json.dumps(metadata))
            with self.assertRaises(ValueError):
                compare_reference(path, arrays, "scenario")

    def test_profiles_keep_distribution_inside_the_requested_box(self):
        for profile in ("smoke", "local", "volume50", "dense50"):
            system, particles = configuration("ecsim_corr", profile)
            self.assertEqual(system["Scheme"], "ecsim_corr")
            space = particles["particles"][0]["distribution"][0]["dist_space"]
            for index, axis in enumerate("XYZ"):
                extent = system["NumCells" + axis] * system["D" + axis.lower()]
                self.assertEqual(2 * space["half_length"][index], extent)
                self.assertEqual(space["center"][index], extent / 2)


if __name__ == "__main__":
    unittest.main()
