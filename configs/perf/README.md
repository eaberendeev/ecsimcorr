# Performance workloads

Short workloads for locating particle, grid, boundary and collision
bottlenecks; not physical validation runs. One generator,
`generic_config.py`, is run through the usual `run.sh`/`build.py`: the first
word of `--config` is the script, the rest are its arguments.

```bash
VALIDATION_PERIODICITY=5 TIMERED_SIM_STEPS=5 TIMER_START_STEP=2 time ./run.sh \
  --jobs 12 --type Release --rerun --timers 1 --timers_memory 0 \
  --config "configs/perf/generic_config.py -type dense -scheme ecsim_corr -collisions 1"
```

The run directory is `Res_<type>_<scheme>[_collisions][_neutrals][_ionization]`;
it holds the inputs, `beren3d.log` (with `Step wall time` per step) and
`profile_it_5.json` (Chrome trace). `TIMER_START_STEP=2` leaves the first
step, with allocations and the reference assembly check, out of the profile;
`VALIDATION_PERIODICITY` equal to the step count keeps reference checks on the
first step only. `python3 configs/perf/generic_config.py -h` lists all options.

## Options

| Option | Values | Default |
|---|---|---|
| `-type` | `dense`, `sparse`, `mixed`, `cylinder`, `cylinder-sparse` | `dense` |
| `-scheme` | `ecsim`, `ecsim_corr`, `symplectic` | `ecsim_corr` |
| `-collisions` | 1: Coulomb e–e, i–i, e–i, ln Λ = 15 | 0 |
| `-neutrals` | 1: add neutral hydrogen with the same PPC | 0 |
| `-ionization` | 1: electron impact ionization of H, `physical_only` (adds neutrals) | 0 |
| `-machine` | `auto`, `desktop`, `server`, `node` (sets the grid) | `auto` |
| `-cells NX NY NZ` | overrides the machine grid | |
| `-ppc` | per species per occupied cell; `mixed`: core value | 200, sparse 2 |
| `-steps` | ≥ 2 | 5 |
| `-procs` | OpenMP threads | physical cores |

`symplectic` supports only the collisionless periodic types (`dense`,
`sparse`, `mixed`). Compare `-ionization 1` against `-neutrals 1`, and
`-collisions 1` against the same type without it.

## Machine classes

`auto` picks the largest class whose minimum the machine meets.
Dense 200 PPC e/i takes about 95 KB per cell (measured peak RSS).

| `-machine` | Minimum | Cells | Dense particles | Memory |
|---|---|---|---:|---:|
| `desktop` | — | 50 × 50 × 30 | 30 million | ~5 GB |
| `server` | 32 cores, 128 GiB | 100 × 100 × 60 | 240 million | ~60 GB |
| `node` | 96 cores, 512 GiB | 160 × 160 × 96 | 983 million | ~240 GB |

On a 128-core node, use `-machine node` for full-node throughput and scaling
runs over `-procs 16/32/64/128`; use `-machine desktop` with the same
`-procs` to see how a small grid scales on many threads.

## Types

dx = dy = dz = 0.5. ECSIM and Corr use dt = 1.2; symplectic uses dt = 0.25.
Electrons and hydrogen ions have masses 1 and 1836, temperatures 0.1 and
0.001 keV, equal macro weights; n0 = 10¹⁸ cm⁻³. No imposed fields, damping,
checkpoints or diagnostic outputs; `VerboseStep` is on.

| `-type` | Initial population | Main target |
|---|---|---|
| `dense` | Periodic, uniform e/i, 200 PPC | Gather, push, deposition, response assembly |
| `sparse` | Same box, 2 PPC | Grid and matrix work against particle work |
| `mixed` | Periodic, 2 PPC background + 198 PPC in the central 1/8 volume | Uneven occupancy, OpenMP load imbalance |
| `cylinder` | Open cylinder, e/i, 200 PPC | Open-end losses, particle removal |
| `cylinder-sparse` | Same cylinder, 2 PPC | Fixed grid costs with open boundaries |

The cylinder wall has radius 0.48 of the smaller transverse box size; plasma
fills radius 0.4 of it over the full length, so particles leave through both
open ends from the first step (about 50% of the box volume is occupied).
`mixed` has density 0.01 n0 in the background and n0 in the core: a neutral
density contrast, not an equilibrium.

Start with `sparse`, then `mixed`, then `dense`. Repeat each run at least
three times with the same binary and thread count and compare medians of
steps 2–5; use `-steps 20` to reduce noise. Timers in the profile are
inclusive and must not be added together. For a quick wiring check use
`-cells 8 8 6 -ppc 8`.
