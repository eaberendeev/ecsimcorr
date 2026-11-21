#!/usr/bin/env python3
"""Numerical PASS/FAIL for Coulomb temperature relaxation, without plotting.

Units: v/c, m/m_e, T/(m_e*c^2), t*omega_pe(n0). Temperature curves come
from the NRL isotropization equation and finite-mass thermal equilibration,
not from the production pair-scattering routine.
"""
import argparse
import csv
import json
import math
import shutil
from pathlib import Path
import subprocess
import tempfile

ME = 9.10938356e-28
QE = 4.80320427e-10
C = 2.99792458e10
QUICK_CASES = (
    "quick_isotropization",
    "quick_equilibration",
    "quick_equilibration_weighted",
    "quick_hot_electrons_w10",
)
ENERGY_LIMIT = 1e-10
# Fixed before running the validation: finite N=4000 and finite-step TA bias.
# Errors are normalized by the INITIAL temperature difference, not a late
# anisotropy that can vanish in statistical noise.
CURVE_LIMIT = 0.10
FINAL_DIFFERENCE_LIMIT = 0.20
# The stored density must be exactly N/npc (main.cpp uses Density=1).
DENSITY_TOLERANCE = 1e-12
# Etotal and T_total come from the same initial velocities but are accumulated
# in different orders; after the exact bulk-velocity subtraction they agree to
# floating-point roundoff.
ENERGY_NORMALIZATION_TOLERANCE = 1e-10


def check_density_metadata(parameters, isotropization):
    """Recompute the kernel density N/npc instead of trusting the metadata.

    main.cpp sets Density=1 and mpw=Density/NumPartPerCell, so the density
    entering the collision rate is N/npc (then multiplied by n0).
    """
    if parameters["test_type"] == "multispecies":
        records = [(species, "N", "npc", "n_density") for species in parameters["species"]]
    else:
        records = [
            (parameters, *keys)
            for keys in (
                (("N", "npc", "n_density"),)
                if isotropization
                else (("N1", "npc1", "n1_density"), ("N2", "npc2", "n2_density"))
            )
        ]
    for source, count_key, npc_key, density_key in records:
        count = source[count_key]
        npc = source[npc_key]
        if count <= 0 or npc <= 0:
            raise ValueError(f"{count_key}/{npc_key} must be positive")
        density = source[density_key]
        expected = count / npc
        if not math.isfinite(density) or abs(density - expected) > DENSITY_TOLERANCE * max(1.0, abs(density)):
            raise ValueError(f"{density_key} disagrees with {count_key}/{npc_key}")


def check_initial_energy_normalization(parameters, isotropization, energy0):
    """Check the macroparticle energy normalization independently of Etotal_rel_err.

    main.cpp uses Density=1 and mpw=Density/NumPartPerCell.  For an isotropic
    species Etotal=0.5*mpw*m*sum(v^2)=1.5*(N/npc)*T_total because
    T_total=m*sum(v^2)/(3N).  The initial bulk velocity is subtracted exactly,
    and the two species add independently.
    """
    if parameters["test_type"] == "multispecies":
        expected = 1.5 * sum(s["N"] / s["npc"] * s["T0"] for s in parameters["species"])
    elif isotropization:
        density = parameters["N"] / parameters["npc"]
        expected = 0.5 * density * (parameters["T_par0"] + 2 * parameters["T_perp0"])
    else:
        expected = 1.5 * (
            parameters["N1"] / parameters["npc1"] * parameters["T1_0"]
            + parameters["N2"] / parameters["npc2"] * parameters["T2_0"]
        )
    if not math.isfinite(expected) or expected <= 0:
        raise ValueError("invalid initial energy normalization")
    if abs(energy0 - expected) > ENERGY_NORMALIZATION_TOLERANCE * expected:
        raise ValueError("initial Etotal disagrees with the N/npc * T normalization")


def plasma_frequency(n0):
    return math.sqrt(4 * math.pi * n0 * QE**2 / ME)


def isotropization_rhs(parameters):
    """dT_parallel/dt=-2*nu*(T_parallel-T_perp); dT_perp/dt=+nu*D."""
    mass = parameters["mass"] * ME
    density = parameters["n_density"] * parameters["n0"]
    prefactor = 2 * math.sqrt(math.pi) * (QE * abs(parameters["charge"])) ** 4 * density
    prefactor *= parameters["coulomb_log"] / (math.sqrt(mass) * plasma_frequency(parameters["n0"]))

    def rhs(parallel, perpendicular):
        anisotropy = perpendicular / parallel - 1
        if abs(anisotropy) < 1e-3:
            angular_factor = (
                4 / 15 - 8 / 35 * anisotropy + 4 / 21 * anisotropy**2 - 16 / 99 * anisotropy**3
            )
        else:
            root = math.sqrt(abs(anisotropy))
            integral = math.atan(root) / root if anisotropy > 0 else math.atanh(root) / root
            angular_factor = (-3 + (anisotropy + 3) * integral) / anisotropy**2
        rate = prefactor * angular_factor / (parallel * ME * C * C) ** 1.5
        exchange = rate * (parallel - perpendicular)
        return -2 * exchange, exchange

    return rhs


def equilibration_rhs(parameters):
    """Finite-mass Glinskij Eq.17 differentiated into dT_e/dt.

    Valid here for singly charged equal-density species, m_light=m_e.
    T_inf=(Te0+Ti0)/2; y=Te/T_inf;
    dy/dt=(1-y)/(tau*[2/M+y*(1-1/M)]^(3/2)).
    """
    if not (
        parameters["mass1"] == 1
        and parameters["mass2"] > 1
        and abs(parameters["charge1"]) == abs(parameters["charge2"]) == 1
        and parameters["n1_density"] == parameters["n2_density"]
    ):
        raise ValueError("thermal theory requires equal densities, |q|=1 and mass1=1<mass2")
    mass_ratio = parameters["mass2"]
    equilibrium = (parameters["T1_0"] + parameters["T2_0"]) / 2
    density = parameters["n1_density"] * parameters["n0"]
    tau_seconds = 3 * (2 * math.pi) ** 1.5 / (4 * parameters["coulomb_log"])
    tau_seconds *= C**3 * ME**2 / (16 * math.pi**2 * density * QE**4)
    tau = tau_seconds * mass_ratio * equilibrium**1.5 * plasma_frequency(parameters["n0"])

    def rhs(electron, ion):
        factor = (2 / mass_ratio + electron / equilibrium * (1 - 1 / mass_ratio)) ** 1.5
        exchange = (equilibrium - electron) / (tau * factor)
        return exchange, -exchange

    return rhs


def multispecies_rhs(parameters):
    """Maxwellian NRL temperature exchange for every pair of species a != b:
    dT_a/dt = sum_b nu_ab (T_b - T_a),
    nu_ab = (8 sqrt(2 pi)/3) n_b q_a^2 q_b^2 lnL sqrt(m_a m_b) / (m_a T_b + m_b T_a)^(3/2)
    (Gaussian units, T in erg). For two equal-density species it reduces to the
    finite-mass equilibration used above."""
    species = parameters["species"]
    unit_time = 1 / plasma_frequency(parameters["n0"])
    coefficients = []
    for a in species:
        row = []
        for b in species:
            density = b["n_density"] * parameters["n0"]
            row.append(
                8 * math.sqrt(2 * math.pi) / 3 * density * (QE * a["charge"]) ** 2 * (QE * b["charge"]) ** 2
                * parameters["coulomb_log"] * math.sqrt(a["mass"] * b["mass"]) * ME
            )
        coefficients.append(row)

    def rhs(*temperatures):
        energy = [t * ME * C * C for t in temperatures]   # erg
        result = []
        for i, a in enumerate(species):
            total = 0.0
            for j, b in enumerate(species):
                if i == j:
                    continue
                rate = coefficients[i][j] / (ME * (a["mass"] * energy[j] + b["mass"] * energy[i])) ** 1.5
                total += rate * (temperatures[j] - temperatures[i])
            result.append(total * unit_time)
        return result

    return rhs


def reference_temperatures(rhs, initial, times, substeps=8192):
    """Integrate a two-temperature physical ODE, independently of PIC kernels."""
    state = list(initial)
    time = 0.0
    result = []
    max_step = times[-1] / substeps
    for target in times:
        while time < target:
            h = min(max_step, target - time)
            k1 = rhs(*state)
            k2 = rhs(*(y + h * k / 2 for y, k in zip(state, k1)))
            k3 = rhs(*(y + h * k / 2 for y, k in zip(state, k2)))
            k4 = rhs(*(y + h * k for y, k in zip(state, k3)))
            state = [
                y + h * (a + 2 * b + 2 * c + d) / 6 for y, a, b, c, d in zip(state, k1, k2, k3, k4)
            ]
            time += h
        result.append(tuple(state))
    return result


# Weighted and equal-weight runs of one mixture start from different random samples.
REFERENCE_INITIAL_TOLERANCE = 0.05


def reference_run(directory, parameters, columns, times):
    """Temperatures of the equal-weight run of the same physical mixture, sample by sample.

    Used where the Maxwellian exchange theory does not hold: the weighted
    scattering must then reproduce the unweighted algorithm (already checked
    against theory) within sampling noise."""
    name = parameters["reference"]
    with (directory / f"{name}_params.json").open() as stream:
        reference = json.load(stream)
    if len(reference["species"]) != len(parameters["species"]) or any(
        not math.isclose(a["n_density"], b["n_density"], rel_tol=1e-12)
        or a["mass"] != b["mass"]
        or a["charge"] != b["charge"]
        or not math.isclose(a["T0"], b["T0"], rel_tol=REFERENCE_INITIAL_TOLERANCE)
        for a, b in zip(parameters["species"], reference["species"])
    ):
        raise ValueError(f"{name} is not the same physical mixture")
    if len({s["npc"] for s in reference["species"]}) != 1:
        raise ValueError(f"{name} must use equal weights")
    with (directory / f"{name}.csv").open(newline="") as stream:
        rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
    if [row["time"] for row in rows] != times:
        raise ValueError(f"{name} is sampled at different times")
    return [tuple(row[key] for key in columns) for row in rows]


def check_case(directory, name):
    with (directory / f"{name}_params.json").open() as stream:
        parameters = json.load(stream)
    if any(
        not math.isfinite(value) for value in parameters.values() if isinstance(value, (int, float))
    ):
        raise ValueError("nonfinite physical parameter")
    with (directory / f"{name}.csv").open(newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("missing/duplicate columns")
        records = [{key: float(value) for key, value in row.items()} for row in reader]
    if len(records) < 8 or any(not math.isfinite(x) for row in records for x in row.values()):
        raise ValueError("need >=8 finite samples")
    times = [row["time"] for row in records]
    if times[0] != 0 or any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError("time must start at zero and increase")
    cadence = parameters["diag_interval"]
    if cadence <= 0 or parameters["num_steps"] <= 0 or parameters["dt"] <= 0:
        raise ValueError("positive steps, cadence and dt are required")
    expected_steps = list(range(0, parameters["num_steps"] + 1, cadence))
    if len(records) != len(expected_steps) or any(
        row["step"] != step or not math.isclose(row["time"], step * parameters["dt"], abs_tol=1e-12)
        for row, step in zip(records, expected_steps)
    ):
        raise ValueError("missing samples or inconsistent step/time columns")
    last_step = parameters["num_steps"] // cadence * cadence
    if not math.isclose(times[-1], last_step * parameters["dt"], rel_tol=1e-12):
        raise ValueError("incomplete relaxation calculation")
    isotropization = parameters["test_type"] == "isotropization"
    multispecies = parameters["test_type"] == "multispecies"
    if parameters["test_type"] not in ("isotropization", "equilibration", "multispecies"):
        raise ValueError("unknown physical problem")
    check_density_metadata(parameters, isotropization)
    if multispecies:
        columns = tuple(f"T{k}" for k in range(len(parameters["species"])))
    else:
        columns = ("T_par", "T_perp") if isotropization else ("T1", "T2")
    if any(row[key] <= 0 for row in records for key in columns):
        raise ValueError("temperatures must stay positive")
    initial = [records[0][key] for key in columns]
    if multispecies:
        initial_metadata = [s["T0"] for s in parameters["species"]]
    else:
        initial_metadata = [parameters[key] for key in (("T_par0", "T_perp0") if isotropization else ("T1_0", "T2_0"))]
    if any(not math.isclose(actual, expected, rel_tol=1e-12) for actual, expected in zip(initial, initial_metadata)):
        raise ValueError("initial temperatures disagree with experiment metadata")
    scale = max(initial) - min(initial)
    if scale <= 0 or min(initial) <= 0:
        raise ValueError("need positive unequal initial temperatures")
    if multispecies:
        rhs = multispecies_rhs(parameters)
    elif isotropization:
        rhs = isotropization_rhs(parameters)
    else:
        rhs = equilibration_rhs(parameters)
    theory = reference_temperatures(rhs, initial, times)
    if parameters.get("reference"):
        theory = reference_run(directory, parameters, columns, times)
    # Recompute from the measured energy column; do not trust Etotal_rel_err.
    energy0 = records[0]["Etotal"]
    if energy0 <= 0:
        raise ValueError("initial energy must be positive")
    check_initial_energy_normalization(parameters, isotropization, energy0)
    energy_error = max(abs(row["Etotal"] / energy0 - 1) for row in records)
    curve_error = math.sqrt(
        sum(
            ((row[key] - expected[d]) / scale) ** 2
            for row, expected in zip(records, theory)
            for d, key in enumerate(columns)
        )
        / (len(columns) * len(records))
    )
    final_spread = max(records[-1][key] for key in columns) - min(records[-1][key] for key in columns)
    theory_spread = max(theory[-1]) - min(theory[-1])
    final_difference_error = abs(final_spread - theory_spread) / scale
    theoretical_relaxation = 1 - theory_spread / scale
    if not all(
        math.isfinite(value)
        for value in (energy_error, curve_error, final_difference_error, theoretical_relaxation)
    ):
        raise ValueError("nonfinite error metric")
    # An almost unchanged curve would make a no-collision implementation pass.
    informative = theoretical_relaxation >= 0.30
    passed = (
        informative
        and energy_error <= ENERGY_LIMIT
        and curve_error <= CURVE_LIMIT
        and final_difference_error <= FINAL_DIFFERENCE_LIMIT
    )
    result = {
        "case": name,
        "pass": passed,
        "max_relative_energy_error": energy_error,
        "rms_temperature_error_over_initial_difference": curve_error,
        "final_temperature_difference_error": final_difference_error,
        "predicted_relaxed_fraction": theoretical_relaxation,
        "limits": {
            "energy": ENERGY_LIMIT,
            "curve": CURVE_LIMIT,
            "final_difference": FINAL_DIFFERENCE_LIMIT,
            "minimum_predicted_relaxation": 0.30,
        },
    }
    return result


def generate_quick_data(binary, directory):
    """A successful process must produce new data, never reuse old CSV files."""
    with tempfile.TemporaryDirectory(prefix="relaxation-", dir=directory) as temporary:
        fresh = Path(temporary)
        subprocess.run([str(binary.resolve()), "--quick"], cwd=fresh, check=True, timeout=240)
        filenames = [name + suffix for name in QUICK_CASES for suffix in (".csv", "_params.json")]
        missing = [name for name in filenames if not (fresh / name).is_file()]
        if missing:
            raise ValueError("calculation did not create fresh results: " + ", ".join(missing))
        # Copy every result, including equal-weight reference runs.
        for path in fresh.iterdir():
            if path.suffix in (".csv", ".json"):
                shutil.copyfile(path, directory / path.name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="*")
    parser.add_argument("--directory", type=Path, default=Path.cwd())
    parser.add_argument("--binary", type=Path, help="run the bounded --quick experiment first")
    args = parser.parse_args()
    if args.binary:
        if any(name not in QUICK_CASES for name in args.cases):
            parser.error("--binary supports only the quick relaxation cases")
        try:
            generate_quick_data(args.binary, args.directory)
        except (OSError, ValueError, subprocess.SubprocessError) as error:
            result = [{"case": "quick_calculation", "pass": False, "error": str(error)}]
            print(f"FAIL quick calculation: {error}")
            (args.directory / "collision_report.json").write_text(
                json.dumps(result, indent=2) + "\n"
            )
            return 1
    cases = args.cases or (
        list(QUICK_CASES)
        if args.binary
        else [
            "iso_electrons",
            "iso_electrons_r5",
            "eq_two_species",
            "eq_two_species_dt1",
            "eq_two_species_real",
            "eq_two_species_weighted",
            "hot_electrons_w10",
            "hot_electrons_w100",
        ]
    )
    results = []
    for name in cases:
        try:
            result = check_case(args.directory, name)
            print(
                ("PASS" if result["pass"] else "FAIL")
                + f' {name}: energy={result["max_relative_energy_error"]:.3g}'
                + f' (limit {ENERGY_LIMIT:g}), curve RMS/DeltaT0={result["rms_temperature_error_over_initial_difference"]:.3g}'
                + f' (limit {CURVE_LIMIT:g}), final difference error={result["final_temperature_difference_error"]:.3g}'
                + f' (limit {FINAL_DIFFERENCE_LIMIT:g}), predicted relaxation={result["predicted_relaxed_fraction"]:.1%}'
            )
        except (
            OSError,
            ValueError,
            KeyError,
            TypeError,
            OverflowError,
            ZeroDivisionError,
        ) as error:
            result = {"case": name, "pass": False, "error": str(error)}
            print(f"FAIL {name}: {error}")
        results.append(result)
    (args.directory / "collision_report.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n"
    )
    return 0 if all(item["pass"] for item in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
