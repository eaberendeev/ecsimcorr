"""Short desktop workloads; see README.md beside this file."""

import argparse
import copy
import json
import math
from pathlib import Path


SCHEMES = ("ecsim", "ecsim_corr", "symplectic")
CASES = (
    "dense", "sparse", "mixed", "cylinder", "cylinder-sparse",
    "cylinder-coulomb", "cylinder-neutral", "cylinder-ionization",
)


def make_config(case="dense", scheme="ecsim_corr", cells=(50, 50, 30),
                ppc=None, steps=5):
    if case not in CASES or scheme not in SCHEMES:
        raise ValueError("unknown case or scheme")
    if len(cells) != 3 or any(n < 4 for n in cells):
        raise ValueError("use at least four cells on each axis")
    if steps < 2 or (ppc is not None and ppc < 1):
        raise ValueError("steps must be >= 2 and ppc must be positive")
    cylinder = case.startswith("cylinder")
    if cylinder and scheme == "symplectic":
        raise ValueError("symplectic supports only collisionless periodic boxes")
    if ppc is None:
        ppc = 2 if case in ("sparse", "cylinder-sparse") else 200
    if case == "mixed" and ppc < 3:
        raise ValueError("mixed requires ppc >= 3 (two background particles plus core)")

    spacing = 0.5
    lengths = [n * spacing for n in cells]
    center = [length / 2 for length in lengths]
    box = {"type": "rectangle", "center": center, "half_length": center}
    space = box
    if cylinder:
        space = {"type": "cylinder_z", "center": center,
                 "radius": 0.4 * min(lengths[:2]), "half_length": 0.5 * lengths[2]}

    def population(region, density, temperature):
        return {"type": "initial", "dist_space": copy.deepcopy(region),
                "dist_pulse": {"type": "gaussian", "mean": [0, 0, 0],
                               "sigma": [temperature] * 3}, "density": density}

    def species(name, charge, mass, temperature):
        distributions = [population(space, 1.0, temperature)]
        if case == "mixed":
            core = {"type": "rectangle", "center": center,
                    "half_length": [length / 4 for length in lengths]}
            # Equal macro weights across the background and core, and e/i.
            distributions = [population(box, 2 / ppc, temperature),
                             population(core, (ppc - 2) / ppc, temperature)]
        return {"Name": name, "Charge": charge, "Mass": mass, "Density": 1.0,
                "NumPartPerCell": ppc, "distribution": distributions}

    particles = {"particles": [species("Electrons", -1.0, 1.0, 0.1),
                               species("Ions", 1.0, 1836.0, 0.001)]}
    boundaries = [{"periodic": {"face": face}} for face in ("XMIN", "YMIN", "ZMIN")]
    if cylinder:
        boundaries = [{"open": [{"face": face} for face in
                                ("CYLINDER", "XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX")]}]
    system = {
        "Scheme": scheme, "Dx": spacing, "Dy": spacing, "Dz": spacing,
        "Dt": 0.25 if scheme == "symplectic" else 1.2,
        "NumCellsX": cells[0], "NumCellsY": cells[1], "NumCellsZ": cells[2],
        "StartTimeStep": 0, "LastTimestep": steps, "StartFromTime": 0,
        "n0": 1.e18, "NumPartPerCell": ppc, "k_particles_reservation": -1.0,
        "Boundary_conditions": boundaries, "DampingType": "None",
        "DampCellsX_glob": [0, 0], "DampCellsY_glob": [0, 0], "DampCellsZ_glob": [0, 0],
        "Collisions": [], "VerboseStep": True,
        "diagnostics": {"outputs": []},
    }
    if cylinder:
        system["CylinderDomain"] = {"radius": 0.48 * min(lengths[:2]), "center": center[:2]}
    if case == "cylinder-coulomb":
        system["Collisions"] = [
            {"type": "coulomb", "species1": a, "species2": b, "coulomb_log": 15.0}
            for a, b in (("Electrons", "Electrons"), ("Ions", "Ions"), ("Electrons", "Ions"))
        ]
    if case in ("cylinder-neutral", "cylinder-ionization"):
        particles["particles"].append(species("Neutrals", 0.0, 1837.0, 0.001))
    if case == "cylinder-ionization":
        system["Collisions"] = [{
            "type": "neutral_ionization", "charged_species": "Electrons",
            "neutral_species": "Neutrals", "electron_product": "Electrons", "ion_product": "Ions",
            "scheme": "physical_only", "electron_ionization": True,
            "proton_ionization": False, "proton_charge_exchange": False,
        }]
    return system, particles


def write_config(directory, system, particles):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for name, config in (("system_config.json", system), ("particles_config.json", particles)):
        (directory / name).write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    wp = math.sqrt(4 * math.pi * system["n0"] * 4.80320427e-10 ** 2 / 9.10938356e-28)
    (directory / "phys.par").write_text(f"w_p = {wp}\n1/w_p = {1 / wp}\n", encoding="utf-8")


# Compatible with build.py --config configs/desktop/gen_config.py.
NumProcs = 4
system_config, particles_config = make_config()


def generate_config():
    write_config(Path.cwd(), system_config, particles_config)
    return "Res_desktop_dense_ecsim_corr"


def add_config_arguments(parser):
    parser.add_argument("--case", choices=CASES, default="dense")
    parser.add_argument("--scheme", choices=SCHEMES, default="ecsim_corr")
    parser.add_argument("--cells", type=int, nargs=3, default=(50, 50, 30), metavar=("NX", "NY", "NZ"))
    parser.add_argument("--ppc", type=int, help="per species, in occupied cells; mixed: core value")
    parser.add_argument("--steps", type=int, default=5)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_config_arguments(parser)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        system, particles = make_config(args.case, args.scheme, args.cells, args.ppc, args.steps)
    except ValueError as error:
        parser.error(str(error))
    args.output.mkdir(parents=True, exist_ok=False)
    write_config(args.output, system, particles)
    print(args.output)
