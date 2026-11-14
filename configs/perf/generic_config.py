"""Short performance workloads; see README.md beside this file.

Used through build.py, with arguments after the file name:
    ./run.sh --config "configs/perf/generic_config.py -type dense -scheme ecsim_corr -collisions 1"
"""

import argparse
import copy
import json
import math
import os


SCHEMES = ("ecsim", "ecsim_corr", "symplectic")
TYPES = ("dense", "sparse", "mixed", "cylinder", "cylinder-sparse")
# Machine class: (cells, threads, minimal physical cores, minimal RAM in GiB); threads=0
# means all physical cores. socket is one 64-core CPU of a stark node, node is the
# whole 2 x 64-core, 1 TB node; their grids are those of working runs (Glinskij
# 120x120x200, mirror trap 140x140x400). Dense 200 PPC e/i needs about 95 KB per
# cell: 24, 270 and 750 GB; the cylinders occupy about half of it.
MACHINES = {
    "desktop": ((50, 50, 100), 0, 0, 0),
    "socket": ((120, 120, 200), 64, 64, 384),
    "node": ((140, 140, 400), 128, 128, 768),
}

NumProcs = 4  # set by generate_config; read by build.py


def physical_cores():
    cores = set()
    physical_id = "0"
    try:
        with open("/proc/cpuinfo") as cpuinfo:
            for line in cpuinfo:
                key, _, value = line.partition(":")
                if key.strip() == "physical id":
                    physical_id = value.strip()
                elif key.strip() == "core id":
                    cores.add((physical_id, value.strip()))
    except OSError:
        pass
    return len(cores) or os.cpu_count() or 4


def memory_gib():
    return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 2**30


def detect_machine():
    cores, memory = physical_cores(), memory_gib()
    fitting = [name for name, (_, _, min_cores, min_memory) in MACHINES.items()
               if cores >= min_cores and memory >= min_memory]
    return fitting[-1]


def make_config(kind, scheme, cells, ppc, steps, collisions, neutrals, ionization):
    if any(n < 4 for n in cells) or steps < 2 or (ppc is not None and ppc < 1):
        raise ValueError("need at least 4 cells per axis, steps >= 2 and positive ppc")
    cylinder = kind.startswith("cylinder")
    neutrals = neutrals or ionization
    if scheme == "symplectic" and (cylinder or collisions or neutrals):
        raise ValueError("symplectic supports only collisionless periodic boxes")
    if ppc is None:
        ppc = 2 if kind.endswith("sparse") else 200
    if kind == "mixed" and ppc < 3:
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
        if kind == "mixed":
            core = {"type": "rectangle", "center": center,
                    "half_length": [length / 4 for length in lengths]}
            # Equal macro weights across the background and core, and e/i.
            distributions = [population(box, 2 / ppc, temperature),
                             population(core, (ppc - 2) / ppc, temperature)]
        return {"Name": name, "Charge": charge, "Mass": mass, "Density": 1.0,
                "NumPartPerCell": ppc, "distribution": distributions}

    particles = {"particles": [species("Electrons", -1.0, 1.0, 0.1),
                               species("Ions", 1.0, 1836.0, 0.001)]}
    if neutrals:
        particles["particles"].append(species("Neutrals", 0.0, 1837.0, 0.001))
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
    if collisions:
        system["Collisions"] += [
            {"type": "coulomb", "species1": a, "species2": b, "coulomb_log": 15.0}
            for a, b in (("Electrons", "Electrons"), ("Ions", "Ions"), ("Electrons", "Ions"))
        ]
    if ionization:
        system["Collisions"].append({
            "type": "neutral_ionization", "charged_species": "Electrons",
            "neutral_species": "Neutrals", "electron_product": "Electrons", "ion_product": "Ions",
            "scheme": "physical_only", "electron_ionization": True,
            "proton_ionization": False, "proton_charge_exchange": False,
        })
    return system, particles


def parse_arguments(argv):
    flag = {"type": int, "choices": (0, 1), "default": 0}
    parser = argparse.ArgumentParser(prog="generic_config.py", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-type", choices=TYPES, default="dense")
    parser.add_argument("-scheme", choices=SCHEMES, default="ecsim_corr")
    parser.add_argument("-collisions", **flag, help="Coulomb e-e, i-i, e-i")
    parser.add_argument("-neutrals", **flag, help="add neutral hydrogen, 'ppc' per cell")
    parser.add_argument("-ionization", **flag, help="electron impact ionization (adds neutrals)")
    parser.add_argument("-machine", choices=("auto", *MACHINES), default="auto",
                        help="desktop 50x50x100 on all cores; socket 120x120x200, 64 threads "
                             "(64+ cores, 384+ GiB); node 140x140x400, 128 threads (128+ cores, 768+ GiB)")
    parser.add_argument("-cells", type=int, nargs=3, metavar=("NX", "NY", "NZ"))
    parser.add_argument("-ppc", type=int, help="per species in occupied cells; mixed: core value")
    parser.add_argument("-steps", type=int, default=15)
    parser.add_argument("-threads", type=int, help="OpenMP threads (default: of the machine class)")
    args = parser.parse_args(argv)
    if args.machine == "auto":
        args.machine = detect_machine()
    cells, threads, _, _ = MACHINES[args.machine]
    args.cells = args.cells or cells
    args.threads = args.threads or threads or physical_cores()
    return parser, args


def generate_config(argv=()):
    global NumProcs
    parser, args = parse_arguments(list(argv))
    try:
        system, particles = make_config(args.type, args.scheme, args.cells, args.ppc, args.steps,
                                        args.collisions, args.neutrals, args.ionization)
    except ValueError as error:
        parser.error(str(error))
    NumProcs = args.threads
    for name, config in (("system_config.json", system), ("particles_config.json", particles)):
        with open(name, "w", encoding="utf-8") as file:
            json.dump(config, file, indent=2)
            file.write("\n")
    wp = math.sqrt(4 * math.pi * system["n0"] * 4.80320427e-10 ** 2 / 9.10938356e-28)
    with open("phys.par", "w", encoding="utf-8") as file:
        file.write(f"w_p = {wp}\n1/w_p = {1 / wp}\n")
    suffix = "".join(f"_{name}" for name in ("collisions", "neutrals", "ionization")
                     if getattr(args, name))
    workdir = f"Res_{args.type}_{args.scheme}{suffix}"
    print(f"{workdir}: {args.machine}, cells {args.cells}, ppc {system['NumPartPerCell']}, "
          f"{args.steps} steps, {NumProcs} threads")
    return workdir


if __name__ == "__main__":
    import sys
    generate_config(sys.argv[1:])
