import math
from pathlib import Path
import struct
import tempfile
import unittest

from compare_simulations import compare_values, configuration, read_field, read_particles, read_plane, snapshot_binary, conserved_counts, check_continuity


class SimulationComparisonTest(unittest.TestCase):
    def test_detects_local_field_error_even_if_integral_is_unchanged(self):
        with self.assertRaises(ValueError):
            compare_values([1.01, 0.99], [1, 1], 1e-8, 1e-10, "field")

    def test_nonfinite_and_missing_values_cannot_pass(self):
        for actual, expected in (([math.nan], [1]), ([1], [math.inf]), ([], []), ([1], [1, 2])):
            with self.subTest(actual=actual), self.assertRaises(ValueError):
                compare_values(actual, expected, 1e-8, 1e-10, "sample")

    def test_tolerance_allows_roundoff_at_zero(self):
        result = compare_values([1e-12, 2 + 1e-12], [0, 2], 1e-8, 1e-10, "sample")
        self.assertGreater(result["limit"], result["max_abs_error"])

    def test_truncated_recovery_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "data"
            path.write_bytes(struct.pack("=3i", 2, 2, 2) + bytes(8))
            with self.assertRaises(ValueError):
                read_field(path)
            path.write_bytes(struct.pack("=i", 2) + bytes(48))
            with self.assertRaises(ValueError):
                read_particles(path)

    def test_density_plane_checks_header_and_payload(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "density"
            for payload in (struct.pack("=2f", math.nan, 2), struct.pack("=2f", 2, 2) + bytes(4)):
                path.write_bytes(payload)
                with self.assertRaises(ValueError):
                    read_plane(path)
            path.write_bytes(struct.pack("=2f12f", 2, 2, *([0.5] * 12)))
            self.assertEqual(read_plane(path), ((2., 2.), [0.5] * 12))

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
