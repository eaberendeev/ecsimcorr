# Beren3D Simulation Framework

This repository contains a 3D particle-in-cell (PIC) simulation code **beren3d** along with Python scripts for building, configuring, and running simulations. The framework is designed for plasma physics simulations with support for complex geometries, boundary conditions, and diagnostics.

## Features

- 3D electromagnetic PIC code with customizable schemes (e.g., `ecsim_corr`).
- Flexible particle injection and boundary conditions.
- OpenMP parallelization + MPI support (future).
- Configurable diagnostics: 1D/2D/3D field outputs, particle probes, radiation diagnostics.
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

Configure simulation parameters by editing `set_params.py`. This file defines:
- Grid dimensions (`NumCellsX_glob`, `Dy`, `Dz`, etc.)
- Time step (`Dt`) and simulation duration (`MaxTime`, `RecTime`)
- Particle species (electrons, ions, neutrals) and their distributions
- Boundary conditions (`BoundTypeX/Y/Z`)
- External fields (`BUniform`, coils)
- Diagnostics (output frequencies, probe positions, radiation planes)
- Work directory name (`DirName`) – used to create a unique folder for each simulation

After editing, the script generates three configuration files:
- `system_config.json` – main simulation parameters
- `particles_config.json` – particle species definitions
- `phys.par` – physical constants (`w_p`, `1/w_p`)

## Running a Simulation

### Local Run

After building and configuring, simply execute:
```bash
./run.sh
```
The script will:
- Build the code if necessary.
- Generate configuration files (using `set_params.py`).
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
- `run.sh`, `build.py`, `set_params.py` – scripts used for the run

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
