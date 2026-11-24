"""Loss cone of a magnetic mirror, test electrons, for any kinetic scheme.

Two coils (R = 12, on axis at z = 8 and z = 56) make a mirror in a box of
16 x 16 x 64 cells, h = 1 c/omega_pe, with open walls: B = 1 at the centre,
5.67 in the throats. A 5 keV isotropic electron cloud starts in a 4-cell cube
at the centre (density 1e-6, so its own field is negligible). Electrons in the
loss cone leave through the throats within one transit; the others bounce, as
long as the scheme keeps their magnetic moment. The time step gives
Omega_e dt = 2 at the centre and 11 in the throats with a thermal path of 0.2
cell per step.

loss_cone_check.py compares the escaped fraction at the end with
1 - sqrt(1 - B0 / Bmax), averaged over the starting points, and requires that
almost no electron leaves in the second half of the run.
"""

import json
import math

Cells, CellSize = (16, 16, 64), 1.0
CoilRadius, CoilZ = 12.0, (8.0, 56.0)
CentreField = 1.0
TemperatureKev = 5.0
CloudHalfLength = 2.0
ParticlesPerCell = 256
TimeStep = 2.0
EndTime = 1500.0


def coil_field(z):
    """On-axis field of the coils per unit current, as fields/Coil.cpp at r = 0."""
    return sum(2 * math.pi * CoilRadius**2 / (CoilRadius**2 + (z - z0) ** 2) ** 1.5 for z0 in CoilZ)


def coil_current():
    return CentreField / coil_field(0.5 * Cells[2] * CellSize)


def configuration(scheme, dt=TimeStep, **extra):
    box = [Cells[a] * CellSize for a in range(3)]
    centre = [0.5 * b for b in box]
    steps = round(EndTime / dt)
    current = coil_current()
    system = {
        "Scheme": scheme,
        "Dt": dt,
        "Dx": CellSize,
        "Dy": CellSize,
        "Dz": CellSize,
        "NumCellsX": Cells[0],
        "NumCellsY": Cells[1],
        "NumCellsZ": Cells[2],
        "StartFromTime": 0,
        "StartTimeStep": 0,
        "LastTimestep": steps,
        "n0": 1e10,
        "k_particles_reservation": 0,
        "Collisions": [],
        "DampingType": "None",
        "Boundary_conditions": [{"open": [{"face": f} for f in ("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX")]}],
        "ExternalFieldB": [{"coils": [{"z": z0, "R": CoilRadius, "I": current} for z0 in CoilZ]}],
        "diagnostics": {"outputs": [{"type": "energy_balance", "interval": max(1, steps // 300)},
                                    {"type": "console_summary"}]},
    }
    system.update(extra)
    electrons = {
        "Name": "Electrons",
        "Charge": -1.0,
        "Mass": 1.0,
        "Density": 1e-6,
        "NumPartPerCell": ParticlesPerCell,
        "distribution": [{
            "type": "initial",
            "sampling": "halton_quiet",
            "dist_space": {"type": "rectangle", "center": centre, "half_length": [CloudHalfLength] * 3},
            "dist_pulse": {"type": "gaussian", "mean": [0, 0, 0], "sigma": [TemperatureKev] * 3},
            "density": 1.0,   # markers per cell relative to NumPartPerCell; the weight is from Density
        }],
    }
    return system, {"particles": [electrons]}


def write_configuration(scheme, **options):
    system, particles = configuration(scheme, **options)
    for filename, config in (("system_config.json", system), ("particles_config.json", particles)):
        with open(filename, "w", encoding="utf-8") as stream:
            json.dump(config, stream, indent=2)
    wp = math.sqrt(4 * math.pi * system["n0"] * 4.80320427e-10 ** 2 / 9.10938356e-28)
    with open("phys.par", "w", encoding="utf-8") as stream:
        stream.write(f"w_p = {wp}\n1/w_p = {1 / wp}\n")
