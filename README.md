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

- GCC with C++20 support (the current CMake configuration selects `g++`).
- CMake ≥ 3.16.
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
  weights of source and product species are handled). Three forms are
  supported:
  - a number (e.g. `0.3`): constant yield per incident particle
    (validated to be `>= 0`);
  - `{"type": "threshold", "yield": 2.0, "threshold_kev": 0.05}`:
    constant yield, applied only when the incident particle's kinetic
    energy is at least `threshold_kev` (keV);
  - `{"type": "vaughan", "delta_max": 2.0, "energy_max_kev": 0.3,
    "threshold_kev": 0.05}`: Vaughan secondary-electron-yield curve
    δ(E) = δ_max·(w·exp(1−w))^k with k = 0.62 for w ≤ 1 and k = 0.25 for
    w > 1, where w = E/E_max (E_max = `energy_max_kev`); emission is zero
    below `threshold_kev` (optional, defaults to 0).
  Guard rails (checked at config load): a constant yield `>= 1` where the
  source species equals the product species is **rejected** — it would be
  self-amplifying and run away; a Vaughan rule with `delta_max > 1` on the
  product species itself only **warns**, since net amplification depends on
  the incident energy spectrum near `energy_max_kev`.
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
    "interval": int,      # output every N timesteps (must be > 0)
    "average_steps": int, # optional 2D averaging window (0/absent = off)
    "species":  ...,      # "all" or ["Electrons", "Ions", ...]
    "fields":   [...],    # ["E", "B", "En", "Bn", "E_ex", "B_init"]
    "planes":   {...},    # {"x": [...], "y": [...], "z": [...]}
    "timesteps": [...],   # [0, 100, 1000] for fields_3d / full3d_all
    "points":   [...],    # [[x1,y1,z1],...] for probes
}
```

`interval` must be a positive integer; zero or negative values are rejected at
startup. For `fields_2d`, `density_2d`, `current_2d` and `pressure_2d`, the
two output modes are enabled independently and only by explicit options:

| Options on a 2D output entry | Files written |
|---|---|
| `"interval": 10` | Ordinary snapshots every 10 timesteps |
| `"average_steps": 10` | Averaged snapshots only |
| `"interval": 10, "average_steps": 20` | Both ordinary and averaged snapshots |
| Neither option, or only `"average_steps": 0` | No snapshots |

Omitting an output type from `diagnostics.outputs` disables that output.
An absent or empty `outputs` list enables no diagnostic outputs.
`average_steps` controls per-output averaging: `0` or absent disables it, while
`N > 0` accumulates the field on every `output()` call (i.e. every timestep,
independently of `interval`) and writes the arithmetic mean every `N` calls as
`..._avg<NNNN>`. A residual partial window is written at the end of the run, so
the last averaged file always reflects the remaining steps. The initial
diagnostic at timestep 0 is included as a sample. Pressure averaging includes
all six tensor components. Average indices count windows, not timesteps.

`energy.txt` and `boundary.txt` are written with 12 significant digits
(`std::setprecision(12)`); the energy-conservation values are often ~1e-13 and
would lose all meaningful digits at the default 6.

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

1. Implement `IDiagnosticOutput` from `srcBeren/diagnostics/IDiagnosticOutput.h` with
   `output(int timestep, const Diagnostics&)` and an idempotent `finalize()` when flushing is needed.
2. Add the implementation to the matching group in `srcBeren/diagnostics/outputs/`;
   prepare its directories in the output constructor with `ensure_directory()` and report file-open failures.
3. Register it in `OutputFactory::create_from_list()` in `OutputFactory.cpp`, parse and validate JSON there,
   then pass ordinary settings to the output. `AveragedField2D.*` stores only the selected averaging planes.

Outputs read prepared diagnostics, fields and particles without changing them.
The legacy `pressure_2d` files contain the raw second moment `sum(m*mpw*S*v_a*v_b)`,
including bulk flow; mean velocity is not subtracted.

Output implementations are grouped as `TableOutputs` (tables/console), `SliceOutputs` (2D fields and moments),
`StateOutputs` (3D snapshots/recovery), `ProbeOutput` and `EnergySpectrumOutput`.
`SliceFiles.*` owns the 2D binary format; generic 3D field serialization remains in `utils/output_util.*`.
Configuration rejects unknown species/field names, malformed probe points and invalid numeric settings
when creating the corresponding output. `fields_3d` supports E, B, En, Bn and E_ex;
`fields_2d` also supports B_init. Probes require a valid linear interpolation stencil, including ghosts;
points on physical walls remain supported. Representable out-of-range slice planes retain the skip behavior.

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

Note: For Debug builds, the binary is left in `_build_debug/bin` and run directly from there (no work directory is created).

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

Build and run the fast unit/regression suite:
```bash
./run.sh --tests
```

The fast suite includes geometry and Yee indexing, boundary conditions,
wall-hit scenarios, diagnostic files and averaging, particle motion in constant
fields, charge continuity, Lmat2 assembly, sparse operators, field solving and
collision invariants. A 5×5×5 periodic plasma experiment checks energy, each
species' charge and the stored density at every step for ECSIM and ECSIM-Corr.

CTest runs current build targets in isolated working directories with two
OpenMP threads and time limits. After building, run
`ctest --test-dir <build-dir> -L unit --output-on-failure` without an install step.

`run_collisions` groups Coulomb conservation, neutral reaction statistics and
the bounded Coulomb temperature-relaxation test. The latter compares measured
temperature curves with theory and checks energy drift and complete output.
Each bounded calculation runs in a fresh directory, so old CSV files cannot
mask missing output. The checker writes `collision_report.json` and returns
failure when a criterion is violated.
The full `run_coulomb` experiment applies the same numerical checks to its CSV data.

`run_validation` includes bounded Coulomb relaxation and the separate
`particles_array` stress regression. The particle fixtures contain millions of
particles and can require substantial time and memory. `run_sparse` is an
explicit large benchmark, outside the automatic correctness suite.

Independent regression tests check Boris invariants, Yee interpolation,
charge continuity (including stencil crossings and reflections), current/mass
matrix consistency, field work, sparse algebra, Maxwell operators, damping
energy, solver breakdown/convergence, neutral transport, and collision
energy/momentum. Plot-producing collision experiments do not replace these
assertions.

To measure executable-line coverage with GCC, without modifying a normal build:

```bash
cmake -S srcBeren -B /tmp/beren-coverage -DBUILD_TESTS=ON -DENABLE_COVERAGE=ON -DUSE_TIMERS=OFF
cmake --build /tmp/beren-coverage --target run_tests -j 3
python3 tests/coverage_summary.py --build-dir /tmp/beren-coverage
```

The resulting JSON excludes third-party libraries and test sources. Coverage
percentages describe this build's executed lines; they do not establish
physical correctness of every model or cover every branch.

Run a specific test:
```bash
./run.sh --tests --test domain             # geometry and Yee indexing
./run.sh --tests --test wall_scenarios     # wall-hit fate and charge/energy bookkeeping
./run.sh --tests --test field_solve        # field assembly/solve checks
./run.sh --tests --test particle_motion    # analytic constant-field trajectories
./run.sh --tests --test collisions         # Coulomb + neutral quantitative checks
./run.sh --tests --test neutral_collision_invariants  # neutral reactions only
./run.sh --tests --test coulomb_relaxation # bounded temperature relaxation
./run.sh --tests --test coulomb            # full relaxation, numerical checks and plots
```

When `--tests` is passed, only tests are built and run — simulation config/workdir generation is skipped.
`--test all` (the default) means the fast unit suite only; the heavier
`run_validation` tests (`particles_array`, `coulomb_relaxation`) are separate.

Available tests:

- `domain` — unit tests for `Vector3`, `Grid`, `Geometry`, `Domain` (Yee grid, node positions, cylinder geometry, reflections, boundary conditions)
- `wall_scenarios` — particle fate (absorption/reflection) at a wall and matching macroparticle charge/energy bookkeeping.
- `field_solve` — assembly and solution of the field system against independent residuals/known solutions.
- `neutral_collision_invariants` — reaction frequency, binomial reaction counts,
  charge exchange and electron-impact ionization energy.
- `coulomb` — Takizuka-Abe relaxation experiments: single-species temperature
  isotropization (vs NRL/Trubnikov rate) and two-species temperature equilibration
  (vs Glinskiy Eq.17). `check_relaxation.py` checks energy drift, temperature
  curves and calculation completion, writes `collision_report.json` and fails
  on an invalid result. Plotting scripts then generate supplementary figures.
  Fast conservation/scattering checks run in `coulomb_invariants`.

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
  - Eigen and AMGCL are vendored in `srcBeren/3rd_party`; no external paths are needed.
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
