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
# DiagDict = {
#     "outTime3D": [5, 150, 200],
#     "zondCoordsLineX": [(x0,y0,z0), ...],   # 1D probes along X
#     "zondCoordsLineY": [(x0,y0,z0), ...],   # 1D probes along Y
#     "zondCoordsLineZ": [(x0,y0,z0), ...],   # 1D probes along Z
#     "sliceFieldsPlaneX": [x0, ...],          # YZ slice planes (2D fields)
#     "sliceFieldsPlaneY": [y0, ...],          # XZ slice planes
#     "sliceFieldsPlaneZ": [z0, ...],          # XY slice planes
#     "TimeStepDelayDiag2D": 1,                # write 2D every N steps
# }
# Automatically written:
#   energy.txt, boundary.txt,
#   Fields/Diag2D/ (E, B slices),
#   Particles/<Name>/Diag2D/ (Current, Density, Pressure components)
#   Recovery/ (checkpoint / restart)

# ---- 10. OTHER PARAMETERS ---------------------------------------------------
# n0 = 1e13                    # reference density (cm^-3)
# Tau = 4998                   # injection timescale (1/wp)
# StartFromTime = 0            # restart time (0 = fresh start)
# RecoveryInterval = 200       # checkpoint every N timesteps
# VerboseStep = True           # extra solver diagnostics to beren3d.log
# k_particles_reservation = -1 # reserve factor (<=0 disables)
# LastTime = 750               # total simulation time (1/wp)
# LastTimestep = int(round(LastTime / Dt + 1))
