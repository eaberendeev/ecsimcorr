# Beren3D Simulation Framework

This repository contains a 3D particle-in-cell (PIC) simulation code **beren3d** along with Python scripts for building, configuring, and running simulations. The framework is designed for plasma physics simulations with support for complex geometries, boundary conditions, and diagnostics.

## Features

- 3D electromagnetic PIC code with customizable schemes (e.g., `ecsim_corr`).
- Flexible particle injection and boundary conditions.
- OpenMP parallelization.
- Configurable diagnostics: config-driven output system with 15 output types (field slices, particle density/current/pressure, energy balance, spectra, 3D dumps, full 3D all, probes, checkpoint).
- Python scripts for automated build and run directory preparation.
- Support for SGE (Sun Grid Engine) clusters.

## Grid and geometry interfaces

`srcBeren/grid/domain.h` combines the fixed `Grid` layout and physical `Geometry`.
`grid.h` owns cell coordinates, ghost offsets and stored node positions; `yee.h` defines
component shifts. `core/vector3.h` and `core/indexing.h` provide the shared vectors and
row-major scalar/vector indexing used by both Grid and Mesh. These headers require
neither JSON nor Eigen, particle storage or field solvers.

`config/domain_config.*` reads the grid dimensions and physical geometry from JSON before
particles and operators are constructed. `config/boundary_config.*` separately creates field
constraints, particle wall laws and secondary-emission models for that domain.
`config/json_settings.h` holds the shared checked settings readers. Simulation owns the domain;
each species borrows it through a const reference.
The domain must outlive its species and its grid layout must remain fixed. Constructing a
species from a temporary domain is rejected at compile time. Recovery reconstructs the domain
from configuration; it does not serialize references or introduce another grid copy.

CMake exports `BerenCore` and `BerenGrid` interface targets and the `Configuration` library.
Clients include the thematic grid/core/config headers directly. Current components share the
electric-field Yee edges: each component is shifted by half a cell along its own axis.
Coordinate conversion and boundary membership use this same placement. Current gather/deposition
stencils are unchanged; the geometry helpers now agree with them.

## Prescribed fields and absorbing layers

`config/field_config.*` groups prescribed fields and absorbing-layer settings. It reads
`ExternalFieldE/B`; `fields/external_fieldsE/B.*` apply the parsed contributions.
Coils carry numerical `z`, `R`, `I` parameters rather than JSON into the field quadrature.
Array composition order, uniform-field defaults and existing Yee placement are preserved.
An incomplete coil is rejected before applying its field.

The same configuration module reads `DampingType`, `DampingCoeff` and `DampCellsX/Y/Z_glob`
once in `SimulationEcsimBase::init_fields()`. The resulting parameters belong to that simulation,
remain fixed while stepping and are reconstructed from configuration on restart.
`fields/Damping.*` applies the existing quadratic profiles and returns the energy removed
from the dynamic E/B endpoints before they become the next step's initial fields.
Unknown damping types are rejected during initialization. The circular layer retains its
existing X-spacing convention and uses the larger of the two X widths.

## ECSIM family: shared predictor and separate temporal steps

`Simulation` owns the domain, particles, sources and collision processes and runs
`prepare_step -> make_step -> collisions -> diagnostics`. `SimulationEcsimBase` is
an abstract shared implementation of field initialization, damping, electric prediction
and output. It owns the E/B levels and field operators; output references remain bound
to those same field objects for the lifetime of the calculation.

`SimulationEcsim` and `SimulationEcsimCorr` are sibling temporal schemes:

- ECSIM drifts `x^(n-1/2) -> x^(n+1/2)`, solves the electric predictor, pushes velocities,
  and advances the endpoint fields. Its endpoint diagnostics deposit a point current.
- Corr deposits the first half-path `x^n -> x^(n+1/2)`, solves the same predictor,
  then pushes velocities and deposits the second half-path. It solves the corrected E,
  rescales velocities, applies the saved wall decisions, redistributes particles and
  advances B. Its diagnostics retain the path current. `fieldJe` and continuity data
  belong only to Corr; wall snapshots remain local to the half-step that consumes them.

Thus the correction reuses the ECSIM predictor, rather than invoking a complete ordinary
ECSIM step. The injection, damping, boundary and collision order is unchanged.

The common base owns response assembly scratch and its validation counter; both schemes
call the existing numerical kernels directly. Whitney migration requires a consistent change
to gather, Jp, L and push, including the fused Corr push. Frozen-orbit weights additionally
need a defined lifetime between assembly and push. This MR keeps the existing Linear kernels
and does not introduce a kernel selector or a new interface for the experimental variants.

The checkpoint format is unchanged. `simulation_restart` compares four continuous steps
with two steps followed by checkpoint/restart for both schemes, including E/B and every
particle at the restored endpoint and each continued step. It excludes collisions, sources,
emission and averaging, whose complete runtime state is not currently checkpointed.

## ECSIM response assembly

`simulation/ecsim/response/` assembles the Linear ECSIM response from species,
the domain, the full magnetic field and the time step. Particle cell grouping
must be current before assembly. `response_kernels.cpp` retains the scalar and
batched particle kernels; `response_assembly.cpp` converts their blocks to CSR.
The dense reference path remains available for validation.

`SimulationEcsimBase` owns separate row contributions and dense reference blocks, and
passes them explicitly from species accumulation to sparse assembly. Row storage is
reserved in the constructor; each assembly clears the previous contributions before
accumulating all species. Dense reference blocks are prepared only
when validation is requested. The workspace contains only temporary CSR merge buffers,
whose contents are overwritten within each assembly call; their allocations are reused.
All response data is reconstructed after restart. The final L matrix remains in `EcsimFieldSystem`;
the scheme applies boundary constraints after assembly. The response uses the same Linear
Yee gather/scatter weights and magnetic rotation as Jp and the particle velocity update.

## Dependencies

- GCC with C++20 support (the current CMake configuration selects `g++`).
- CMake ≥ 3.16.
- Eigen3 (header-only linear algebra library).
- Python ≥ 3.6 with standard libraries (`argparse`, `json`, `shutil`, `subprocess`).
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
| `charge_continuity` | ECSIM-Corr: `continuity.txt`, full-grid L2/Linf norms of the charge continuity residual |
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

Density is deposited only when requested by density slices, charged-species totals, species probes or
`full3d_all`; recovery stores particles and fields without density grids. Density averages request every
sample in their window, including step 0. Charged-species density is charge density; neutral density is
particle number density. Outputs sample the particle state at diagnostic time, including particles
created by emission or reactions after transport.

Charge continuity is optional: add `{"type": "charge_continuity", "interval": 10}` for ECSIM-Corr.
Without that entry no continuity densities or residual are computed. A selected completed step samples
both of its own endpoints, after injection and before secondary emission is flushed, and uses its
trajectory current: `(rho_after - rho_before)/dt + div J`. It does not use the previous diagnostic snapshot.
`continuity.txt` contains `Time L2 Linf`, including ghosts as in the former logged norm. Step 0 has no
completed transport interval and produces no row. ECSIM point current cannot supply this discrete
trajectory check, so that scheme rejects the entry. After restart, endpoints are measured afresh.

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
   `output(int timestep, const Diagnostics&, const DiagnosticMeasurements&)` and an idempotent `finalize()`
   when flushing is needed. Override `request_measurements()` for derived quantities.
2. Add the implementation to the matching group in `srcBeren/diagnostics/outputs/`;
   prepare its directories in the output constructor with `ensure_directory()` and report file-open failures.
3. Register it in `OutputFactory::create_from_list()` in `OutputFactory.cpp`, parse and validate JSON there,
   then pass ordinary settings to the output. `AveragedField2D.*` stores only the selected averaging planes.

Outputs read prepared diagnostics, fields and particles without changing them.
`Simulation` prepares scheme-specific currents and calls the read-only energy measurements in
`EnergyMeasurements.*`. After collisions, `write_diagnostics()` collects the active outputs' requests
through `request_measurements()`, calculates one `DiagnosticMeasurements` frame, then calls the writers.
Densities, pressure moments, spectra, charged-species totals and probe interpolation belong to this measurement phase;
writers only serialize samples and accumulate output windows. Duplicate requests share results.
A frame lives for one diagnostic call and is never saved in recovery. Pressure uses one reusable scalar
3D scratch buffer and retains only the union of requested planes; totals also retain only planes.
ECSIM prepares a point current, while ECSIM-Corr preserves the trajectory current deposited by its step.
The legacy `pressure_2d` files contain the raw second moment `sum(m*mpw*S*v_a*v_b)`,
including bulk flow; mean velocity is not subtracted.

Output implementations are grouped as `TableOutputs` (tables/console), `SliceOutputs` (2D fields and moments),
`StateOutputs` (3D snapshots/recovery), `ProbeOutput` and `EnergySpectrumOutput`.
`SliceFiles.*` owns the 2D binary format; generic 3D field serialization remains in `io/field_io.*`.
Configuration rejects unknown species/field names, malformed probe points and invalid numeric settings
when creating the corresponding output. `fields_3d` supports E, B, En, Bn and E_ex;
`fields_2d` also supports B_init. Probes require a valid linear interpolation stencil, including ghosts;
points on physical walls remain supported. Representable out-of-range slice planes retain the skip behavior.

## Collision processes

`collisions/` owns Coulomb and neutral processes and the ordered `CollisionManager`.
`config/collision_config.*` parses their JSON once during `Simulation::init()`;
`collision_step()` then applies one ordered list after the scheme step and before diagnostics.
The numerical processes receive species and `dt`, without a domain argument, JSON or field operators.
Neutral cross sections and reactions live in `collisions/neutral/`.

Accepted `Collisions` entries run in array order. A missing or empty list disables collisions.
Non-array lists, unknown types and entries without required names retain their warning/skip behavior;
JSON field type errors throw. Settings are fixed after initialization.

The legacy `Collider` selector and its numerical implementations have been removed. Any `Collider`
key, including `"Collider": "None"`, is rejected with a migration message. Remove it from old inputs.
Replace `"Collider": "Coulomb"` with one explicit channel to retain the old selection of species:

```json
{"Collisions": [{"type": "coulomb", "species1": "Electrons", "species2": "Ions", "coulomb_log": 15}]}
```

Same-species channels must be requested separately. Migration does not reproduce legacy trajectories:
the retained Coulomb implementation uses `N*mpw` density, rejects unequal interspecies macro weights,
uses channel/cell/apply RNG streams, and handles both signs of axial relative velocity correctly.

`CoulombCollisionOperator` performs a Takizuka–Abe rotation of relative velocity. In code units,
with reduced mass `mu = m1*m2/(m1+m2)`, its small-angle variance is

```text
variance = f * (wp(n0)^3 / (c^3*n0)) * (lnL*q1^2*q2^2*n*dt / (8*pi*mu^2*|v1-v2|^3))
```

Here `n = N*mpw` for one species and `min(N1*mpw1, N2*mpw2)` for two species;
`mpw = density/particles_per_cell`, without an additional cell-volume factor.
Even populations form random pairs; an odd population first uses three pairs of three particles
with `f = 1/2`, then ordinary pairs with `f = 1`. Different species use the existing shuffled
large/small-population mapping. Each rotation conserves pair momentum and kinetic energy;
interspecies collisions reject unequal macro weights before changing velocities.
The zero-relative-speed guard and the large-variance sampling branch remain in `coulomb.cpp`.

Neutral channels select partners in the same cell, apply ionization or charge exchange,
insert the configured products and remove ionized neutrals. They do not move particles or redistribute
cells. ECSIM transports neutrals in its ordinary drift; ECSIM-Corr calls `advance_neutral_particles(dt)`.
Each scheme transports each neutral once per step, with its normal boundaries and cell update, before
collisions. Adding channels therefore adds no drift. This fixes the previous extra neutral displacement
after each channel; runs with neutral collisions intentionally change their particle trajectories.
Channel order still matters when earlier reactions change the particles available to later channels.
The neutral density passed to the reaction kernel remains `N_neutral/particles_per_cell`.
Cross sections convert normalized speed to physical collision rates using `n0` and `c/wp`.
The reaction kernel supports electron and proton projectiles in hydrogen only. Mass is in units of
`m_e`; the legacy `mass < 2` switch selects electron cross sections within that model. It cannot
identify positrons or other species from mass alone; those projectiles are outside its scope.

Collision state lives for the simulation, without a per-step reset. Coulomb channels retain
an ordered-species seed and an `apply_index` incremented by each positive-time `collide_*` call;
cell streams depend on channel, apply index and cell, independently of OpenMP thread count.
Neutral channels retain their base Lehmer engine (seed 77), with per-cell jump stride `2^16`,
and advance the engine after application. Reaction acceptance and channel selection share the collider
object's `mt19937_64`, seeded from the clock; no engine is constructed per reaction. Scattering still
uses a shared `random_device`-seeded engine; neutral thread-count reproducibility is not guaranteed.
A reproducible restart would need these engine states and Coulomb apply indices, in addition
to the input settings and particle state. Current recovery does not save them; this refactoring
does not change that contract or repair the neutral RNG behavior.

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
- `field_configuration` — analytic cylinder/coil fields, composition, damping profiles and removed field energy.
- `domain_config` — domain JSON and box/cylinder geometry, linked with configuration alone, without the simulation runtime.
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
  - Eigen is vendored in `srcBeren/3rd_party`; no external paths are needed.
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

### Runtime and file I/O

`app/main.cpp` reads the two JSON files and runs the simulation lifecycle.
`config/` decodes physical settings; `io/` owns text input, 3D field files and
checkpoint serialization. The checkpoint format and restart paths are unchanged.
`runtime/` holds logging, environment options, RNG primitives and profiling.
Includes use explicit `io/` and `runtime/` paths; forwarding headers are absent.

`USE_TIMERS=ON` enables the detailed flat/tree profiler in `runtime/profiler.*`.
The former `diagnostics/Timer` and its `Performance/globalFunctions.time` output
have been removed; existing `RECORD_TIMER` measurements cover the step and kernels.
`USE_TIMERS=OFF` compiles those detailed measurements out.

With `VerboseStep=true`, each production step reports `Step wall time: <seconds> s`
through `LOG_STEP`, on the console and in `beren3d.log`, regardless of `USE_TIMERS`.
With `VerboseStep=false`, this message is suppressed with the other step details.
The monotonic wall time includes preparation, particle/field updates, collisions
and diagnostics, including their I/O. Initialization, step-zero diagnostics,
profile serialization and the timing log message itself are excluded. The log
uses actual resumed step numbers and is overwritten when logging is initialized.
Only two clock reads are needed per step; no detailed events are stored when
profiling is off.

The build uses target-scoped include paths, OpenMP and third-party requirements.
Common compiler flags, profiling definitions and coverage settings are set once
in the top-level CMake file before adding any project subdirectories.
Cross-directory project includes spell their paths
from `srcBeren`; local sibling includes remain valid. Executables list their
compiled providers explicitly rather than inheriting an aggregate `LIBRARIES`
list. The existing OBJECT/STATIC library types and numerical compiler options
are preserved. All project targets inherit the common build settings;
`beren_core` carries geometry headers and their C++20 requirement.
Numerical unit tests share the full runtime link set. Only `domain` and
`domain_config` keep minimal link sets to check their standalone API.

A dependency is `PUBLIC` when a public header needs its includes, definitions
or compiler options, and `PRIVATE` when only implementation files need them.
For example, configuration and neutral-collision headers contain inline
OpenMP regions, so their libraries export OpenMP. ECSIM and ECSIM-Corr use
`BicgstabSolver`;
AMGCL and MPI are absent from the build configuration and project includes.
The unused vendored AMGCL directory is retained for a separate removal MR.

Configuration construction still refers to particles, collisions, fields and
the shared boundary API. These existing source dependencies are not a claim of
an acyclic architecture: executable link lists supply the required object files
when resolving the configuration archive. Geometry configuration remains usable
without the numerical runtime, as checked by `domain_config`.

With `BUILD_TESTS=ON`, `architecture_dependencies` checks project include paths,
the established lower-module boundaries, and forbids directory-scoped CMake
commands such as global `include_directories` and `link_libraries`. Common
compile and link settings are allowed in the top-level CMake file only.
It permits only the documented exact
configuration/boundary include exceptions and rejects lower-layer dependencies
on concrete simulation steps. Third-party sources are excluded. This is a
static source-level check: it does not prove that the complete dependency graph
is acyclic and does not inspect effective compiler flags for each target.
`architecture_checker` creates small examples with valid dependencies and deliberate
errors, then checks that the source checker accepts or rejects them as expected.
Run the source check directly with:

```bash
python3 tests/check_architecture.py --source-dir srcBeren
```

Detailed profiling can be restricted without rebuilding:

```bash
TIMER_START_STEP=101 TIMER_END_STEP=120 TIMERED_SIM_STEPS=10 ./beren3d
```

The inclusive window enables collection before each step; outside it, instrumented
scopes skip clock reads and event/tree collection. The tree root is still restarted
by `clearTree()` independently of collection. The flag is changed only in serial
orchestration and read by OpenMP workers. Initialization is collected and retained
with step-zero diagnostics when the window includes zero; otherwise the initial
events are cleared before step zero. Standalone kernels and benchmarks keep
profiling enabled by default.
`TIMERED_SIM_STEPS` controls batch writes, not step sampling: the example writes
`profile_it_110.json` and `profile_it_120.json`. The end of the window flushes an
incomplete batch. Defaults are every 20 steps and an unrestricted window;
`TIMERED_SIM_STEPS=0` rewrites an accumulating `profile.json` each profiled step.
Step wall-time logging is independent of this window.

`Grid` distinguishes logical node indices from array indices: logical `(0,0,0)`
is the unshifted node at the origin, while storage index `(g,g,g)` denotes that
node for `g` ghost layers. `to_storage_index` and `to_node_index` translate between
these indices. `node_position` includes the requested component's Yee shift;
`field_node_index` subtracts that shift and selects the lower logical node using
`floor`. Selecting the support of a particle shape remains a separate operation.
`storage_node_position` and `storage_field_node_index` explicitly use storage indices.
Raw array access, flattened `sind/vind` indices and memory layout are unchanged.

`GridNodeView<Storage, OffsetX, OffsetY, OffsetZ>` provides logical-node access
to scalar arrays (three indices) and component fields (four indices). It holds
only a reference to existing storage. Offsets are compile-time constants,
defaulting to `GHOST_CELLS` in every direction; no `Grid` is stored. Its implicit
constructor accepts lvalues, so gather calls pass fields directly. Storage must
outlive the view. Const storage provides const access. A class-level
`static_assert` rejects unsupported storage types. Component checks and the
underlying container's raw bounds checks apply in Debug; the view does not
validate grid dimensions or logical bounds. It performs no periodic wrapping,
clipping or atomic updates and is not saved in restart files. View tests cover
one and two ghost layers; production gather fixes the default compile-time
offset and does not accept a runtime ghost width.

Point interpolation now uses `coupling/stencil.h`: `Stencil3D<N>` stores logical
node origins and tensor-product weights, without a ghost offset. The point
factory selects Linear support from `floor(coord-shift)` and Quadratic support
from `floor(coord-shift+0.5)-1`. The factory with an explicit origin evaluates
weights on those nodes and also subtracts the Yee shift. Unshifted and half-cell
weights are shared between field components and between E/B gathers.

ECSIM velocity/current prediction, point current, the ECSIM-Corr fused push
and diagnostic E/B probes use these logical stencils. Gather uses the default
compile-time view offset;
point-current scatter uses the runtime grid offset. Linear gather retains
a factored eight-corner kernel; the generic tensor sum remains available; old raw kernels are retained only in
`tests/unit/raw_gather_reference.h` for independent regression comparisons. Unlike the old raw interpolation's integer
truncation, the new point factory selects the correct lower node for negative
coordinates wherever its support fits the allocated halo. This extension is
checked by independent constant/affine-field and gather/scatter-work tests.

ECSIM response remains Linear. Quadratic point-factory tests do not enable
Quadratic response. Response kernels use logical stencil origins while retaining
their storage-coordinate arithmetic for weights;
`Domain::to_storage_cell_coordinates` supplies normalized coordinates with the
runtime ghost offset. Complete simulation support for a different ghost width
requires checking response assembly, gather and boundary operators together.
Point-current cell buffers retain their existing size and atomic scatter;
their logical origin is selected once for the particle group and passed explicitly.

Trajectory current uses the same `Stencil3D` container. Both endpoints of each
segment share a four-node logical origin, selected with `make_stencil_at`;
`make_point_stencil` is not used to select trajectory support. Short paths use
the existing cell buffer, while paths leaving its support use a rebased buffer.
Reflected legs retain the full timestep normalization. `ParticleShape` and its
old gather overloads have been removed; `Shape`, `Shape2` and `ShapeKernel` remain.
Current scatter converts the logical buffer origin through the runtime grid.
The full buffer footprint selects the fast atomic scatter or checked scatter;
this replaces assumptions based on compile-time ghost width and bin position.
The optimized particle traversal includes populated lower ghost bins. Checked
trajectory scatter retains its existing clipping at allocated field bounds;
this migration does not change wall physics or make an insufficient halo valid.

The raw second velocity moment also uses a logical Linear stencil and translates
its origin through the runtime grid. It remains a second moment including bulk
flow, with normalization `m*mpw*sum(v1*v2*S)`, rather than thermal pressure.

Density deposition also uses logical stencils. Both the reference particle
loop and the blocked kernel evaluate a four-node envelope beginning at
`floor(coord/cell_size)-1`, independent of ghost width. The blocked kernel
retains its 2-by-2-by-2 cell grouping and local buffer; bin indices are converted
to logical origins, then the buffer origin is translated to storage once for
atomic scatter. It visits every bin, including populated ghost bins that the
previous blocked loop silently skipped. Species
bins must correspond to the runtime domain grid and contain particles in the
indicated cells.

The deposit overwrites the output and preserves code-density normalization:
`q*mpw` per charged particle, `mpw` per neutral particle, without a cell-volume
factor. Boundary folding remains the caller's operation. Zero padding outside
the allocated halo is skipped; a nonzero contribution outside it stops the run
instead of losing charge. Density tests cover ghost widths 1 and 2 for Linear
and Quadratic shapes. Trajectory tests cover complete local buffer placement, ghost widths 0/1/2
where the nonzero support fits, and charge continuity at walls; these checks do
not establish support for a different ghost width in the complete simulation.

External cylinder and coil fields evaluate physical Yee-node positions through
`Domain` and center their transverse profiles on the physical box, independent
of halo width. The Bphi axial interval is `[startz, endz)` in physical units;
Bx/By nodes have a half-cell z shift. This corrects the previous conversion of
physical bounds directly into storage indices without the halo offset. Coil
radial components vanish on the axis, avoiding division by zero there. The
charged-cylinder radial cutoff is unchanged.

Particle bin selection uses the lower logical cell (`floor(coord/cell_size)`)
and then the runtime storage offset. Coordinates outside allocated particle
bins are rejected rather than being truncated into bin zero. Redistribution
still requires motion between neighboring bins and does not replace boundary
handling. Probe validation checks actual integer and half-integer Linear
stencils; it requires the runtime halo to match the compile-time gather offset.

Response gather and Lmat2 now carry Linear weights in `Stencil3D<2>` with
logical origins. Their private builder deliberately retains the previous
storage-coordinate fractional-weight arithmetic and integer truncation; replacing
it with a shape-function evaluation would change roundoff. Raw indices are
formed at the field/block access boundary. This preserves the existing response
support contract, rather than extending it to arbitrary negative storage
coordinates. Reference, scalar and batched assembly remain independently tested.
