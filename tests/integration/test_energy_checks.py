import json
import tempfile
import unittest
from pathlib import Path

from run_integration import check_energy_conservation, check_nan_inf


class EnergyValidationTest(unittest.TestCase):
    def result(self, contents, species="Beam", finite_only=False):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / "energy.txt").write_text(contents)
            (path / "particles_config.json").write_text(json.dumps({"particles": [{"Name": species}]}))
            return (check_nan_inf if finite_only else check_energy_conservation)(directory)[0]

    def test_custom_species_is_included_in_balance(self):
        self.assertTrue(self.result("Time Beam energyConserve\n0 10 0\n1 10 0.00001\n"))
        self.assertFalse(self.result("Time Beam energyConserve\n0 10 0\n1 10 2\n"))

    def test_nonfinite_overflow_and_malformed_tables_fail(self):
        for row in ("0 nan 0", "0 inf 0", "0 1e999 0", "0 invalid 0", "0 1", "0 1 0 extra"):
            with self.subTest(row=row):
                self.assertFalse(self.result("Time Beam energyConserve\n" + row + "\n"))
                self.assertFalse(self.result("Time Beam energyConserve\n" + row + "\n", finite_only=True))

    def test_zero_energy_does_not_skip_conservation(self):
        self.assertFalse(self.result("Time Beam energyConserve\n0 0 1\n"))
        self.assertTrue(self.result("Time Beam energyConserve\n0 0 0\n"))

    def test_field_energy_provides_scale_without_particles(self):
        self.assertTrue(self.result("Time Beam energyFieldE energyConserve\n0 0 10 0\n1 0 10 0.0001\n"))

    def test_headers_are_not_numeric_data(self):
        self.assertTrue(self.result("Time InfiniteBeam energyConserve\n0 1 0\n", "InfiniteBeam", True))
        self.assertFalse(self.result("Time Beam Beam energyConserve\n0 1 1 0\n"))
        self.assertFalse(self.result("Time Other energyConserve\n0 1 0\n"))

    def test_initial_snapshot_is_not_a_completed_calculation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / "energy.txt").write_text("Time Beam energyConserve\n0 1 0\n")
            (path / "particles_config.json").write_text(json.dumps({"particles": [{"Name": "Beam"}]}))
            (path / "system_config.json").write_text(json.dumps({"LastTimestep": 4, "Dt": 0.5}))
            self.assertFalse(check_energy_conservation(directory)[0])
            (path / "energy.txt").write_text("Time Beam energyConserve\n0 1 0\n2 1 0\n")
            self.assertTrue(check_energy_conservation(directory)[0])


if __name__ == "__main__":
    unittest.main()
