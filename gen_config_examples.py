# =============================================================================
#   gen_config_examples.py — reference for all available options
#   Copy-paste from this file into gen_config.py as needed.
# =============================================================================

# ---- 1. SCHEME --------------------------------------------------------------
Scheme_name = "ecsim"          # implicit moment method, single push dt
Scheme_name = "ecsim_corr"     # ecsim + energy-conserving correction

# ---- 2. GRID & TIME ---------------------------------------------------------
Dx = Dy = Dz = 0.3            # cell size
Dt = 1.0                      # timestep (1/omega_p). Keep v*dt < cell size
NumCellsX = 30                # number of cells in X
NumCellsY = 30                # number of cells in Y
NumCellsZ = 5                 # number of cells in Z (thin = 2D-like)

# ---- 3. DOMAIN GEOMETRY -----------------------------------------------------
# Box (default — just set NumCellsX/Y/Z):
#   Domain is [0, NumCellsX*Dx] x [0, NumCellsY*Dy] x [0, NumCellsZ*Dz]

# Cylinder in XY cross-section (truncates corners):
# Cyl = {"radius": 50.0, "center": [50.0, 50.0]}
# system_config["CylinderDomain"] = Cyl

# ---- 4. BOUNDARY CONDITIONS -------------------------------------------------
# Faces: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, CYLINDER
# Multiple types can be combined (e.g. periodic in X/Y, open in Z).

# Open — particles lost, fields zeroed:
# {"open": {"face": "ZMIN"}}

# Periodic — must pair MIN/MAX on same axis:
# {"periodic": {"face": "ZMIN"}}

# Secondary electron emission — ions hitting Z-face produce electrons:
# {"second_emisson": {"face": "ZMAX", "mean": [0,0,0], "sigma": [Tx,Ty,Tz]}}

# Er0 — radial electric field ring at Z boundary:
# {"er0": {"face": "ZMIN", "inner_radius": 5.0, "width": 2.0,
#          "potential_drop": 0.1 / MC2}}

# Electron reflection — reflect electrons at Z within radius R
# if z-kinetic energy below threshold (replaces open for that face):
# {"electron_reflection": {"face": "ZMIN", "radius": 5.0,
#                          "energy_threshold": 0.1 / MC2}}

# ---- 5. DAMPING -------------------------------------------------------------
DampingType = "None"           # "None" | "CircleXY" | "Rectangle"
DampCellsX_glob = [0, 0]      # [left, right] damping cells in X
DampCellsY_glob = [0, 0]      # [left, right] damping cells in Y
DampCellsZ_glob = [0, 0]      # [left, right] damping cells in Z
DampingCoeff = 0.8            # damping factor at the outer edge of the layer

# ---- 6. EXTERNAL FIELDS -----------------------------------------------------
# Electric:
# system_config["ExternalFieldE"] = [
#     {"uniform_field": {"value": [ex, ey, ez]}},
#     {"uniformly_charged_cylinder": {"radius": 25.0, "value": -0.025}},
# ]

# Magnetic (summed):
# system_config["ExternalFieldB"] = [
#     {"uniform_field": {"value": [bx, by, bz]}},
#     {"coils": [{"z": z0, "R": R0, "I": I0}, ...]},  # Biot-Savart sum
# ]

# ---- 7. PARTICLE SPECIES ----------------------------------------------------
# Each species is a dict appended to particles_config["particles"].

# Electrons example:
# {
#     "Name": "Electrons",
#     "Charge": -1.0,
#     "Density": 1.0,
#     "Mass": 1.0,
#     "NumPartPerCell": 100,
#     "distribution": [
#         {
#             "type": "initial",         # populate once at t=0
#             "density": 1.0,
#             "dist_space": {...},
#             "dist_pulse": {...},
#         },
#     ],
# }

# dist_space — spatial shape of the particle cloud:
#   rectangle:       {"type": "rectangle", "center": [x,y,z],
#                     "half_length": [dx,dy,dz]}
#   cylinder_z:      {"type": "cylinder_z", "center": [x,y,z],
#                     "radius": R, "half_length": L}
#   cylinder_x:      {"type": "cylinder_x", "center": [x,y,z],
#                     "radius": R, "half_length": L}
#   cylinder_ring_z: {"type": "cylinder_ring_z", "center": [x,y,z],
#                     "r1": R1, "r2": R2, "half_length": L}

# dist_pulse — velocity/momentum distribution:
#   gaussian:        {"type": "gaussian", "mean": [vx,vy,vz],
#                     "sigma": [Tx,Ty,Tz]}     # sigma in keV
#   rigid_rotation:  {"type": "rigid_rotation",
#                     "rotation_center": [x,y,z], "omega": w_rad_s}
#   tangential:      {"type": "tangential",
#                     "rotation_center": [x,y,z],
#                     "mean_speed": V, "sigma_speed": dV,
#                     "thermal_sigma": [Tx,Ty,Tz]}  # keV

# Distribution types:
#   "initial"         — populate once at t=0
#   "injection"       — continuous injection each timestep
#   "injection_bound" — injection from boundary with velocity*dt offset

# Neutral species (omit distribution or set Charge=0):
# { "Name": "Neutrals", "Charge": 0, "Density": 1.0,
#   "Mass": 1837.0, "NumPartPerCell": 100 }

# ---- 8. COLLISIONS ----------------------------------------------------------
# Collider = "None"  # legacy; ignored if Collisions[] is non-empty

# Coulomb collisions:
# {"type": "coulomb", "species1": "Electrons", "species2": "Electrons",
#  "coulomb_log": 15}
# {"type": "coulomb", "species1": "Electrons", "species2": "Ions",
#  "coulomb_log": 15}

# Neutral ionization / charge exchange:
# {"type": "neutral_ionization",
#  "charged_species": "Electrons", "neutral_species": "Neutrals",
#  "electron_product": "Electrons", "ion_product": "Ions",
#  "scheme": "physical_only",       # or "null_collision"
#  "electron_ionization": True,
#  "proton_ionization": True,
#  "proton_charge_exchange": True}

# ---- 9. DIAGNOSTICS ---------------------------------------------------------
# Diagnostics are configured via the "outputs" array in DiagDict.
# Each entry is a dict with:
#
#   type        — one of the types below
#   interval    — output every N timesteps (default 1)
#   species     — "all" or ["Electrons", "Ions", ...]  (per-species outputs)
#   fields      — ["E", "B", "En", "Bn", "E_ex", "B_init"]  (field outputs)
#   planes      — {"x": [...], "y": [...], "z": [...]}  (2D slice positions)
#   timesteps   — [0, 100, 1000]  (for fields_3d, full3d_all: specific timesteps)
#   points      — [[x1,y1,z1],...]  (for probes: probe coordinates)
#
# Available types:
#   energy_balance      — energy.txt (per-step energy balance)
#   boundary_stats      — boundary.txt (per-face particle statistics)
#   console_summary     — per-step terminal output
#   recovery            — checkpoint/restart
#   fields_2d           — 2D slices of E/B fields (uses planes + fields)
#   density_2d          — 2D slices of particle density per species
#   current_2d          — 2D slices of particle current per species
#   pressure_2d         — 2D slices of pressure tensor (6 components) per species
#   energy_spectrum     — energy spectrum text file per species
#   fields_3d           — full 3D field dump at specific timesteps
#   full3d_all          — full 3D dump of all fields + per-species density/current
#   charge_density_2d   — 2D slice of total charge density (all species sum)
#   total_current_2d    — 2D slice of total current (all species sum)
#   probes              — point probes: all fields at given (x,y,z) coordinates

# --- Example 9a: Default (matches original hardcoded behavior) ---
# diag_outputs = [
#     {"type": "energy_balance",    "interval": 1},
#     {"type": "boundary_stats",    "interval": 1},
#     {"type": "console_summary"},
#     {"type": "recovery",          "interval": RecoveryInterval},
#     {"type": "fields_2d",         "fields": ["E", "B"],
#      "planes": {"x": [cx], "y": [cy], "z": [cz]},
#      "interval": TimeStepDelayDiag2D},
#     {"type": "density_2d",        "species": "all",
#      "planes": {"x": [cx], "y": [cy], "z": [cz]},
#      "interval": TimeStepDelayDiag2D},
#     {"type": "current_2d",        "species": "all",
#      "planes": {"x": [cx], "y": [cy], "z": [cz]},
#      "interval": TimeStepDelayDiag2D},
#     {"type": "pressure_2d",       "species": "all",
#      "planes": {"x": [cx], "y": [cy], "z": [cz]},
#      "interval": TimeStepDelayDiag2D},
#     {"type": "energy_spectrum",   "species": "all",
#      "interval": TimeStepDelayDiag2D},
# ]

# --- Example 9b: Minimal (skip expensive pressure & spectrum) ---
# diag_outputs = [
#     {"type": "energy_balance",    "interval": 1},
#     {"type": "boundary_stats",    "interval": 1},
#     {"type": "console_summary"},
#     {"type": "recovery",          "interval": 500},
#     {"type": "fields_2d",         "fields": ["E", "B"],
#      "planes": {"x": [cx], "z": [cz]}, "interval": 10},
#     {"type": "density_2d",        "species": "all",
#      "planes": {"z": [cz]},       "interval": 10},
#     {"type": "current_2d",        "species": "all",
#      "planes": {"z": [cz]},       "interval": 10},
# ]

# --- Example 9c: Custom intervals per output ---
# diag_outputs = [
#     {"type": "energy_balance",    "interval": 1},
#     {"type": "boundary_stats",    "interval": 1},
#     {"type": "console_summary"},
#     {"type": "recovery",          "interval": 200},
#     {"type": "fields_2d",         "fields": ["E", "B"],
#      "planes": {"z": [0.0]},      "interval": 5},
#     {"type": "density_2d",        "species": "all",
#      "planes": {"z": [0.0]},      "interval": 5},
#     {"type": "current_2d",        "species": ["Electrons"],
#      "planes": {"z": [0.0]},      "interval": 5},
#     {"type": "pressure_2d",       "species": ["Electrons"],
#      "planes": {"z": [0.0]},      "interval": 50},   # every 50 steps
#     {"type": "energy_spectrum",   "species": ["Electrons"],
#      "interval": 100},                                # every 100 steps
# ]

# --- Example 9d: With 3D field dumps and charge diagnostics ---
# diag_outputs = [
#     {"type": "energy_balance",    "interval": 1},
#     {"type": "boundary_stats",    "interval": 1},
#     {"type": "console_summary"},
#     {"type": "recovery",          "interval": 500},
#     {"type": "fields_2d",         "fields": ["E", "B", "E_ex"],
#      "planes": {"z": [0.0]},      "interval": 10},
#     {"type": "density_2d",        "species": "all",
#      "planes": {"z": [0.0]},      "interval": 10},
#     {"type": "current_2d",        "species": "all",
#      "planes": {"z": [0.0]},      "interval": 10},
#     {"type": "fields_3d",         "fields": ["E", "B"],
#      "timesteps": [0, 500, 1000]},                   # full 3D dump
#     {"type": "charge_density_2d",
#      "planes": {"z": [0.0]},      "interval": 10},
#     {"type": "total_current_2d",
#      "planes": {"z": [0.0]},      "interval": 10},
# ]

# --- Example 9e: Different planes for fields vs particles ---
# diag_outputs = [
#     {"type": "energy_balance",    "interval": 1},
#     {"type": "boundary_stats",    "interval": 1},
#     {"type": "console_summary"},
#     {"type": "fields_2d",         "fields": ["E", "B"],
#      "planes": {"x": [cx], "y": [cy], "z": [cz]},
#      "interval": 10},
#     {"type": "density_2d",        "species": "all",
#      "planes": {"z": [cz]},       "interval": 10},   # only XY slices
#     {"type": "current_2d",        "species": "all",
#      "planes": {"z": [cz]},       "interval": 10},
#     # No pressure, no spectrum
# ]
# 
# --- Example 9f: Probes at specific points (Probes/probes.txt) ---
# diag_outputs = [
#     {"type": "energy_balance",    "interval": 1},
#     {"type": "console_summary"},
#     {"type": "probes",
#      "points": [[30.0, 40.0, 5.0], [50.0, 60.0, 10.0]],
#      "interval": 10},
# ]
#
# --- Example 9g: Full 3D dump at specific timesteps (Full3D/ directory) ---
# diag_outputs = [
#     {"type": "energy_balance",    "interval": 1},
#     {"type": "console_summary"},
#     {"type": "full3d_all",
#      "species": "all",
#      "timesteps": [0, 500, 1000, 5000]},
# ]

# ---- 10. OTHER PARAMETERS ---------------------------------------------------
# n0 = 1e13                    # reference density (cm^-3)
# Tau = 4998                   # injection timescale (1/wp)
# StartFromTime = 0            # restart time (0 = fresh start)
# RecoveryInterval = 200       # checkpoint every N timesteps
# VerboseStep = True           # extra solver diagnostics to beren3d.log
# k_particles_reservation = -1 # reserve factor (<=0 disables)
# LastTime = 750               # total simulation time (1/wp)
# LastTimestep = int(round(LastTime / Dt + 1))

# ---------------------------------------------------------------------------
# In gen_config.py, wire up DiagDict like this:
#
# DiagDict = {}
# DiagDict["outputs"] = diag_outputs
# system_config["diagnostics"] = DiagDict
# ---------------------------------------------------------------------------
