"""Validate the PASS/FAIL reader; these synthetic tables are not PIC validation."""

import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

CHECKER = Path(__file__).resolve().parents[2] / "srcBeren/tests/coulomb/check_relaxation.py"
SPEC = importlib.util.spec_from_file_location("coulomb_relaxation_check", CHECKER)
CHECK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECK)


class RelaxationReportTest(unittest.TestCase):
    def fixture(self, directory, frozen=False):
        parameters = dict(
            test_type="isotropization",
            mass=1,
            charge=-1,
            n_density=4000 / 6,
            N=4000,
            npc=6,
            n0=1e13,
            coulomb_log=15,
            T_par0=4e-4,
            T_perp0=2e-4,
            dt=1,
            num_steps=1200,
            diag_interval=20,
        )
        # main.cpp normalization: Etotal=0.5*(N/npc)*(T_par+2*T_perp).
        initial_energy = 0.5 * parameters["n_density"] * (parameters["T_par0"] + 2 * parameters["T_perp0"])
        times = list(range(0, 1201, 20))
        temperatures = CHECK.reference_temperatures(
            CHECK.isotropization_rhs(parameters), [4e-4, 2e-4], times
        )
        if frozen:
            temperatures = [[4e-4, 2e-4] for _ in times]
        rows = [
            dict(step=t, time=t, T_par=a, T_perp=b, Etotal=initial_energy, Etotal_rel_err=0)
            for t, (a, b) in zip(times, temperatures)
        ]
        (directory / "sample_params.json").write_text(json.dumps(parameters))
        self.write_rows(directory, rows)
        return rows

    def write_rows(self, directory, rows):
        with (directory / "sample.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=["step", "time", "T_par", "T_perp", "Etotal", "Etotal_rel_err"]
            )
            writer.writeheader()
            writer.writerows(rows)

    def test_theory_table_passes_and_frozen_collision_fails(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            self.fixture(path)
            self.assertTrue(CHECK.check_case(path, "sample")["pass"])
            self.fixture(path, frozen=True)
            self.assertFalse(CHECK.check_case(path, "sample")["pass"])

    def test_actual_energy_drift_fails_even_if_reported_residual_is_zero(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            rows = self.fixture(path)
            rows[-1]["Etotal"] *= 1.01
            self.write_rows(path, rows)
            self.assertFalse(CHECK.check_case(path, "sample")["pass"])

    def test_inconsistent_metadata_is_rejected_independently(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            self.fixture(path)
            parameters = json.loads((path / "sample_params.json").read_text())
            parameters["n_density"] = 1.0
            (path / "sample_params.json").write_text(json.dumps(parameters))
            with self.assertRaisesRegex(ValueError, "disagrees"):
                CHECK.check_case(path, "sample")
            self.fixture(path)
            with (path / "sample.csv").open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            rows[0]["Etotal"] = str(float(rows[0]["Etotal"]) * 1.01)
            self.write_rows(path, rows)
            with self.assertRaisesRegex(ValueError, "normalization"):
                CHECK.check_case(path, "sample")

    def equilibration_fixture(self, directory):
        parameters = dict(
            test_type="equilibration",
            n0=1e13,
            coulomb_log=15,
            dt=0.5,
            num_steps=100,
            diag_interval=10,
            charge1=-1,
            mass1=1,
            npc1=1,
            N1=4000,
            n1_density=4000.0,
            T1_0=4e-4,
            charge2=1,
            mass2=10,
            npc2=1,
            N2=4000,
            n2_density=4000.0,
            T2_0=4e-5,
        )
        times = list(range(0, 101, 10))
        temperatures = CHECK.reference_temperatures(
            CHECK.equilibration_rhs(parameters), [parameters["T1_0"], parameters["T2_0"]], times
        )
        # main.cpp normalization: Etotal=1.5*sum((N/npc)*T_total).
        initial_energy = 1.5 * (
            parameters["N1"] / parameters["npc1"] * parameters["T1_0"]
            + parameters["N2"] / parameters["npc2"] * parameters["T2_0"]
        )
        rows = [
            dict(
                step=step,
                time=step * parameters["dt"],
                T1=a,
                T2=b,
                Etotal=initial_energy,
                Etotal_rel_err=0,
            )
            for step, (a, b) in zip(times, temperatures)
        ]
        (directory / "equilibration_params.json").write_text(json.dumps(parameters))
        with (directory / "equilibration.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=["step", "time", "T1", "T2", "Etotal", "Etotal_rel_err"]
            )
            writer.writeheader()
            writer.writerows(rows)
        return rows

    def test_equilibration_density_metadata_is_rejected_for_both_species(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            for key in ("n1_density", "n2_density"):
                with self.subTest(key=key):
                    self.equilibration_fixture(path)
                    parameters = json.loads((path / "equilibration_params.json").read_text())
                    parameters[key] = 1.0
                    (path / "equilibration_params.json").write_text(json.dumps(parameters))
                    with self.assertRaisesRegex(ValueError, "disagrees"):
                        CHECK.check_case(path, "equilibration")

    def test_nan_missing_rows_and_wrong_time_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            for problem in ("nan", "missing", "truncated", "time"):
                with self.subTest(problem=problem):
                    rows = self.fixture(path)
                    if problem == "nan":
                        rows[4]["T_par"] = float("nan")
                    if problem == "missing":
                        rows.pop(4)
                    if problem == "truncated":
                        rows.pop()
                    if problem == "time":
                        rows[4]["time"] += 1
                    self.write_rows(path, rows)
                    with self.assertRaises(ValueError):
                        CHECK.check_case(path, "sample")

    def test_success_without_new_output_cannot_reuse_previous_results(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            for name in CHECK.QUICK_CASES:
                for suffix in (".csv", "_params.json"):
                    (directory / (name + suffix)).write_text("previous result")
            # A process may exit successfully without producing a calculation.
            # Existing files in the destination must not hide that failure.
            with mock.patch.object(CHECK.subprocess, "run", return_value=None):
                with self.assertRaisesRegex(ValueError, "did not create fresh results"):
                    CHECK.generate_quick_data(directory / "no-output-binary", directory)
            for name in CHECK.QUICK_CASES:
                self.assertEqual((directory / (name + ".csv")).read_text(), "previous result")

    def test_temperature_ode_conserves_energy_and_converges(self):
        parameters = dict(mass=1, charge=-1, n_density=4000 / 6, n0=1e13, coulomb_log=15)
        rhs = CHECK.isotropization_rhs(parameters)
        initial = [4e-4, 2e-4]
        coarse = CHECK.reference_temperatures(rhs, initial, [0, 1200])[-1]
        # Halve the reference ODE step at the same observation time.
        fine = CHECK.reference_temperatures(rhs, initial, [0, 1200], substeps=16384)[-1]
        self.assertLess(abs(coarse[0] + 2 * coarse[1] - 8e-4), 1e-16)
        self.assertLess(coarse[0] - coarse[1], 0.7 * (initial[0] - initial[1]))
        for a, b in zip(coarse, fine):
            self.assertLess(abs(a - b), 1e-14)


if __name__ == "__main__":
    unittest.main()
