# Beren3D Simulation Framework

This repository contains a 3D particle-in-cell (PIC) simulation code **beren3d** along with Python scripts for building, configuring, and running simulations. The framework is designed for plasma physics simulations with support for complex geometries, boundary conditions, and diagnostics.

## Features

- 3D electromagnetic PIC code with customizable schemes (e.g., `ecsim_corr`).
- Flexible particle injection and boundary conditions.
- OpenMP parallelization + MPI support (future).
- Configurable diagnostics: config-driven output system with 14 output types (field slices, particle density/current/pressure, energy balance, spectra, 3D dumps, full 3D all, probes, checkpoint).
- Python scripts for automated build and run directory preparation.
- Support for SGE (Sun Grid Engine) clusters.

## Dependencies

- C++ compiler with C++17 support (gcc ≥ 8, clang ≥ 7).
- CMake ≥ 3.12.
- Eigen3 (header-only linear algebra library).
- AMGCL (header-only C++ library for solving systems of linear equations with algebraic multigrid).
- Python ≥ 3.6 with standard libraries (`argparse`, `json`, `shutil`, `subprocess`).
- MPI (optional, for distributed runs).
- OpenMP (required for shared-memory parallelism).

## Installation

1) Clone the repository:
   ```bash
   git clone https://github.com/eaberendeev/beren3d.git
   cd beren3d
   ```
2) Install Eigen and AMGCL (if not already in system paths). Set environment variables `EIGEN_PATH` and `AMGCL_PATH` pointing to the root directories of these libraries, or pass them via command-line arguments.

### Build using run.sh wrapper

```bash
./run.sh [options]
```

This script sets the library paths and calls `build.py` with any passed arguments.

Common build options:
- `--type Debug|Release` — build type (default: Release)
- `--rebuild` — clean and rebuild from scratch
- `--rerun` — remove existing work directory (if any) before preparing a new run
- `--jobs N` — number of parallel build jobs (default: number of CPU cores)
- `--eigen PATH` — override Eigen path
- `--amgcl PATH` — override AMGCL path

Example:
```bash
./run.sh --type Debug --rebuild
```

## Configuration

Configure simulation parameters by editing **`gen_config.py`** directly. This file defines:

| Parameter | Description |
|---|---|
| `NumCellsX/Y/Z`, `Dx`, `Dy`, `Dz` | Grid dimensions and cell size |
| `Dt`, `MaxTime`, `Tau` | Time step, duration, injection timescale |
| `Scheme_name` | `"ecsim"` or `"ecsim_corr"` |
| `Boundary_conditions` | Open, periodic, reflecting, secondary emission |
| `ExternalFieldB` | Uniform B-field + coil configuration |
| `CylinderDomain` | Optional cylindrical boundary in XY |
| `diag_outputs` | Array of diagnostic output descriptors |
| `DirName` | Work directory prefix |

After running `build.py`, the script generates three configuration files:
- `system_config.json` – main simulation parameters
- `particles_config.json` – particle species definitions
- `phys.par` – physical constants (`w_p`, `1/w_p`)

**Note:** These JSON files are generated artifacts — edit `gen_config.py`, not the JSON files directly.

### Secondary electron emission (`second_emission`)

Secondary emission is configured per face and is **not** a consuming condition:
it fires for particles *removed* on that face by another condition (`open`,
`bphi`, `electron_reflection`), so a consuming condition must be configured on
the same face (validated at load). The `product` species must exist in
`particles_config.json` (validated after particle initialization).

```json
"Boundary_conditions": [
  {"open": {"face": "ZMAX"}},
  {"second_emission": {"face": "ZMAX", "product": "Electrons", "sources": [
      {"species": "Ions",      "yield": 0.3, "energy": {"type": "fixed",       "kev": 2.0}},
      {"species": "Electrons", "yield": 0.1, "energy": {"type": "temperature", "temperature_kev": [1.0, 1.0, 1.0], "mean": [0.0, 0.0, 0.0]}},
      {"species": "Ions2",     "yield": 0.5, "energy": {"type": "fraction",    "fraction": 0.5}}
  ]}}
]
```

- `yield` — mean number of secondaries per **physical** incident particle
  (fractional yields are sampled stochastically; differing macro-particle
  weights of source and product species are handled).
- Energy types: `fixed` — monoenergetic (`kev`); `temperature` — per-component
  temperature `temperature_kev` in keV, converted as σ_v = sqrt(kT/(mc²)),
  with optional drift `mean` in code velocity units (same convention as the
  `gaussian` injection distribution); `fraction` — fraction (0..1) of the
  incident particle's kinetic energy transferred to each secondary.
- Secondaries are emitted with a Lambertian (cosine-law) angular distribution
  into the domain. Emitted energy is accumulated in `totalEmitEnergy` and
  included in the `energyConserve` balance.
- Boundary particle handling is two-phase and lock-free: conditions classify
  particles in parallel (per-thread buffers), secondary emission runs
  sequentially in canonical cell order — results are deterministic and
  independent of the OpenMP thread count.

### Diagnostics (configurable output system)

All diagnostic output is configured via the `diag_outputs` list in `gen_config.py`. Each entry is a dict with:

```python
{
    "type":     str,      # one of the types below
    "interval": int,      # output every N timesteps
    "species":  ...,      # "all" or ["Electrons", "Ions", ...]
    "fields":   [...],    # ["E", "B", "En", "Bn", "E_ex", "B_init"]
    "planes":   {...},    # {"x": [...], "y": [...], "z": [...]}
    "timesteps": [...],   # [0, 100, 1000] for fields_3d / full3d_all
    "points":   [...],    # [[x1,y1,z1],...] for probes
}
```

**Available output types:**

| type | Description |
|---|---|
| `energy_balance` | `energy.txt` — per-step energy balance |
| `boundary_stats` | `boundary.txt` — per-face particle statistics |
| `console_summary` | Per-step terminal output |
| `recovery` | Checkpoint / restart |
| `fields_2d` | 2D slices of E/B fields at configured planes |
| `density_2d` | 2D slices of particle density per species |
| `current_2d` | 2D slices of particle current per species |
| `pressure_2d` | 2D slices of pressure tensor (6 components) — **expensive** |
| `energy_spectrum` | Energy spectrum text file per species |
| `fields_3d` | Full 3D field dump at specific timesteps |
| `charge_density_2d` | 2D slice of total charge density (all species) |
| `total_current_2d` | 2D slice of total current (all species) |
| `full3d_all` | Full 3D dump of all fields + per-species density/current |
| `probes` | Point probes: all fields at given (x,y,z) coordinates → `probes.txt` |

**Default config** (matches original hardcoded behavior):

```python
_default_planes = {
    "x": sliceFieldsPlaneX,
    "y": sliceFieldsPlaneY,
    "z": sliceFieldsPlaneZ,
}
diag_outputs = [
    {"type": "energy_balance",    "interval": 1},
    {"type": "boundary_stats",    "interval": 1},
    {"type": "console_summary"},
    {"type": "recovery",          "interval": RecoveryInterval},
    {"type": "fields_2d",         "fields": ["E", "B"],
     "planes": _default_planes,   "interval": TimeStepDelayDiag2D},
    {"type": "density_2d",        "species": "all",
     "planes": _default_planes,   "interval": TimeStepDelayDiag2D},
    {"type": "current_2d",        "species": "all",
     "planes": _default_planes,   "interval": TimeStepDelayDiag2D},
    {"type": "pressure_2d",       "species": "all",
     "planes": _default_planes,   "interval": TimeStepDelayDiag2D},
    {"type": "energy_spectrum",   "species": "all",
     "interval": TimeStepDelayDiag2D},
]
```

**Optimization:** pressure and energy spectrum are expensive (iterate all particles). Set higher `interval` or remove the entry to skip them entirely. See `gen_config_examples.py` for complete examples with all options.

**Adding a new output type in C++:**
1. Create a class inheriting `IDiagnosticOutput` (in `srcBeren/diagnostics/DiagnosticOutput.h`)
2. Implement `output(timestep, Diagnostics&)`
3. Register in `OutputFactory::create_from_list()`

## Running a Simulation

### Local Run

After building and configuring, simply execute:
```bash
./run.sh
```
The script will:
- Build the code if necessary.
- Generate configuration files (using `gen_config.py`).
- Create a work directory named according to `DirName` (appended with grid and particle settings).
- Copy the binary, configuration files, source tree, and utility scripts into the work directory.
- `cd` into the work directory and launch `beren3d` with `numactl` for optimal memory placement and OpenMP threading.

Note: For Debug builds, the binary is left in `_build/bin` and run directly from there (no work directory is created).

### Cluster Run

The `run.sh` script also contains SGE headers (lines starting with `#$`). To submit a job to a cluster:
```bash
qsub run.sh [options]
```
The script automatically detects the number of slots requested (`-pe smp N`) and sets `OMP_NUM_THREADS` accordingly. Make sure to adjust the queue name (`-q plasma@en067.binp.gpf`) to match your cluster.

## Work Directory Structure

Upon a successful Release run, a directory like `Res_Jz_m0.01_Dx_0.5_np_1000_Dt_1.5` is created containing:
- `beren3d` – the executable
- `system_config.json`, `particles_config.json`, `phys.par` – configuration files
- `srcBeren/` – copy of the source code (for reproducibility)
- `PlotScripts/` – plotting utilities (if any)
- `run.sh`, `build.py`, `gen_config.py` – scripts used for the run

## Testing

Build and run all tests:
```bash
./run.sh --tests
```

Run a specific test:
```bash
./run.sh --tests --test domain      # unit tests (Vector3, Grid, Geometry, Domain)
./run.sh --tests --test collision   # neutral collision tests
./run.sh --tests --test coulomb     # Coulomb (Takizuka-Abe) collision validation
```

When `--tests` is passed, only tests are built and run — simulation config/workdir generation is skipped.

Available tests:
- `domain` — unit tests for `Vector3`, `Grid`, `Geometry`, `Domain` (Yee grid, node positions, cylinder geometry, reflections, boundary conditions)
- `collision` — Vahedi-Surendra neutral collision model tests (ionization, charge exchange)
- `coulomb` — Takizuka-Abe Coulomb collision validation: single-species temperature
  isotropization (vs NRL/Trubnikov rate) and two-species temperature equilibration
  (vs Glinskiy Eq.17). After running the test binary, the analysis scripts
  (`compare_theory.py`, `plot_coulomb.py`) are invoked automatically.

### Coulomb test output

The `coulomb` test writes all of its artifacts (CSV data, `*_params.json`, and the
comparison PNGs) to a dedicated, gitignored folder `test_results/coulomb/` — the
source tree is never polluted. The analysis/plot scripts live in
`srcBeren/tests/coulomb/` and are run with `test_results/coulomb/` as their working
directory by `build.py`, so the figures land alongside the data.

To regenerate the plots manually after a run, from the results folder:
```bash
cd test_results/coulomb
python3 ../../srcBeren/tests/coulomb/compare_theory.py   # theory vs simulation
python3 ../../srcBeren/tests/coulomb/plot_coulomb.py     # raw diagnostics
```

## Troubleshooting

- Build fails:
  - Ensure Eigen and AMGCL paths are correct. Use `--eigen` and `--amgcl` to specify them.
  - Check CMake output for missing dependencies.
- "Work directory exists" error:
  - Use `--rerun` to remove it, or manually delete the directory.
- Segmentation faults / numerical issues:
  - Verify grid spacing and time step satisfy CFL condition.
  - Check particle densities and injection rates.

## Contributing

Feel free to open issues or pull requests. When contributing, please maintain the existing code style and update documentation accordingly.

## License

enjoy using it
