#!/usr/bin/env python3
"""Ion-acoustic wave ringdown with drift-kinetic electrons and ions.

The reduced analogue of xpic drift_kinetic_ex5_iaw: a periodic column
1 x 1 x 32 cells along B = z, Te = 5 keV, Ti = 0.1 keV, mi / me = 25,
T_perp = 0 (one-dimensional dynamics along the field), quiet-start markers
with a mode-1 density modulation dn_i / n = 0.03. The electrons start in
Boltzmann equilibrium with the initial field E_z = -Te dn_e/dz / n (as xpic
SetCosineField), with dn_e / dn_i = 1 / (1 + k^2 lambda_De^2) from the Gauss
law; without that field the initial charge separation would launch Langmuir
oscillations (omega_pe dt = 20, aliased). Both species get the mean flow
u = (omega / k) dn / n of the forward wave (continuity), which suppresses the
backward wave and most of the ballistic ion response of a density-only load.
check_iaw.py fits the mode-1 ion density to exp(-gamma t) (P cos omega t +
Q sin omega t) and compares omega and gamma with the root of the kinetic
dispersion relation (iaw_theory.py): omega = 3.9304e-3, gamma = 4.7291e-4.
"""

import json
import math

IntegrationHeavy = True   # 33k markers x 160 steps, about six minutes on 8 cores
NumProcs = 8
LastTimestep = 160        # two wave periods (2 pi / omega = 1599), omega dt = 0.079
DirName = "Res_IAW_Drift_Kinetic"
Cells = (1, 1, 32)
Amplitude = 0.03                       # dn_i / n
ElectronTemperatureKev, IonTemperatureKev, IonMass = 5.0, 0.1, 25.0
_k = 2 * math.pi / Cells[2]            # dz = 1
_te = ElectronTemperatureKev / 511.0   # Te / (m_e c^2) = lambda_De^2 for n = 1
ElectronAmplitude = Amplitude / (1 + _k**2 * _te)
# Field seen by the markers -Te a_e k cos(kz); the S1 gather of the Yee E_z and
# the Esirkepov charge shape reduce the harmonic by sinc^2(k dz / 2).
_sinc = math.sin(0.5 * _k) / (0.5 * _k)
InitialEz = -_te * ElectronAmplitude * _k / _sinc**2
# Phase velocity of the kinetic root, iaw_theory.ion_acoustic_root: omega = 3.9304083806759e-3.
PhaseVelocity = 3.9304083806759e-3 / _k

system_config = {
    "Scheme": "drift_kinetic",
    "Dt": 20.0,
    "Dx": 1.0,
    "Dy": 1.0,
    "Dz": 1.0,
    "NumCellsX": Cells[0],
    "NumCellsY": Cells[1],
    "NumCellsZ": Cells[2],
    "StartFromTime": 0,
    "StartTimeStep": 0,
    "LastTimestep": LastTimestep,
    "n0": 1e10,
    "k_particles_reservation": 0,
    "Collisions": [],
    "DampingType": "None",
    "Boundary_conditions": [{"periodic": {"face": face}} for face in ("XMIN", "YMIN", "ZMIN")],
    "DriftKinetic": {"magnetic": [0, 0, 1], "electric": [0, 0, 0]},
    "InitialFieldE": [{"cosine_field": {"amplitude": [0, 0, InitialEz], "mode": [0, 0, 1]}}],
    "DriftKineticPush": {"coord_is_gc": True},
    "Eccapfim": {"coupled_fields": True, "preconditioner": "mass_matrix", "coupled_max_outer": 40},
    "diagnostics": {
        "outputs": [
            {"type": "energy_balance", "interval": 1},
            {"type": "console_summary"},
            {"type": "probes", "interval": 1, "species": ["Electrons", "Ions"],
             "points": [[0.5, 0.5, k + 0.5] for k in range(Cells[2])]},
        ]
    },
}


def species(name, charge, mass, temperature_kev, amplitude):
    return {
        "Name": name,
        "Charge": charge,
        "Mass": mass,
        "Density": 1.0,
        "NumPartPerCell": 512,
        "model": "drift_kinetic",
        "distribution": [{
            "type": "initial",
            "sampling": "halton_quiet",
            "stratified_axis": 2,
            "perturbation": {"amplitude": amplitude, "mode": 1, "phase": 0.0, "flow": PhaseVelocity * amplitude},
            "dist_space": {"type": "rectangle", "center": [0.5, 0.5, 16], "half_length": [0.5, 0.5, 16]},
            "dist_pulse": {"type": "gaussian", "mean": [0, 0, 0], "sigma": [0, 0, temperature_kev]},
            "density": 1.0,
        }],
    }


particles_config = {"particles": [
    species("Electrons", -1.0, 1.0, ElectronTemperatureKev, ElectronAmplitude),
    species("Ions", 1.0, IonMass, IonTemperatureKev, Amplitude),
]}


def generate_config():
    for filename, config in (("system_config.json", system_config),
                             ("particles_config.json", particles_config)):
        with open(filename, "w", encoding="utf-8") as stream:
            json.dump(config, stream, indent=2)
    wp = math.sqrt(4 * math.pi * system_config["n0"] * 4.80320427e-10 ** 2 / 9.10938356e-28)
    with open("phys.par", "w", encoding="utf-8") as stream:
        stream.write(f"w_p = {wp}\n1/w_p = {1 / wp}\n")
    return DirName


if __name__ == "__main__":
    generate_config()
