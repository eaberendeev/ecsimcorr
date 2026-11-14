# Performance workloads

Short workloads for locating particle, grid, boundary and collision
bottlenecks; not physical validation runs. One generator,
`generic_config.py`, is run through the usual `run.sh`/`build.py`: the first
word of `--config` is the script, the rest are its arguments.

```bash
TIMERED_SIM_STEPS=15 TIMER_START_STEP=2 time ./run.sh \
  --jobs 12 --type Release --rerun --timers 1 --timers_memory 0 \
  --config "configs/perf/generic_config.py -type dense -scheme ecsim_corr -collisions 1"
```

The run directory is `Res_<type>_<scheme>[_collisions][_neutrals][_ionization]`;
it holds the inputs, `beren3d.log` (with `Step wall time` per step) and
`profile_it_15.json` (Chrome trace). `TIMER_START_STEP=2` leaves the first
step, with its allocations, out of the profile.
`python3 configs/perf/generic_config.py -h` lists all options.

## Options

| Option | Values | Default |
|---|---|---|
| `-type` | `dense`, `sparse`, `mixed`, `cylinder`, `cylinder-sparse` | `dense` |
| `-scheme` | `ecsim`, `ecsim_corr`, `symplectic` | `ecsim_corr` |
| `-collisions` | 1: Coulomb e–e, i–i, e–i, ln Λ = 15 | 0 |
| `-neutrals` | 1: add neutral hydrogen with the same PPC | 0 |
| `-ionization` | 1: electron impact ionization of H, `physical_only` (adds neutrals) | 0 |
| `-machine` | `auto`, `desktop`, `socket`, `node` (sets grid and threads) | `auto` |
| `-cells NX NY NZ` | overrides the machine grid | |
| `-ppc` | per species per occupied cell; `mixed`: core value | 200, sparse 2 |
| `-steps` | ≥ 2 | 15 |
| `-threads` | OpenMP threads | of the machine class |

`symplectic` supports only the collisionless periodic types (`dense`,
`sparse`, `mixed`). Compare `-ionization 1` against `-neutrals 1`, and
`-collisions 1` against the same type without it.

## Machine classes

`auto` picks the largest class whose minimum the machine meets. `socket` is one
64-core CPU of a stark node (2 × AMD EPYC 7773X, 1 TB), `node` the whole node.
The socket grid is that of a working run, 120 × 120 × 200 (Glinskij article
tests). The node grid is the mirror trap 140 × 140 × 400 cut to 240 cells in z:
a dense response has 123 entries per field row and the matrices use 32-bit
indices, so a dense box is limited to about 5.8 million cells (the full
grid fails with a segmentation fault in the first predictor assembly). Dense
200 PPC e/i takes 80–95 KB per cell (measured peak RSS); the cylinders occupy
about half of it.

| `-machine` | Minimum | Cells | Threads | Dense particles | Dense memory |
|---|---|---|---:|---:|---:|
| `desktop` | — | 50 × 50 × 30 | physical cores | 30 million | ~8 GB |
| `socket` | 64 cores, 384 GiB | 120 × 120 × 200 | 64 | 1.15 billion | ~270 GB |
| `node` | 128 cores, 768 GiB | 140 × 140 × 240 | 128 | 1.9 billion | ~390 GB |

On a node, use `-machine node` for full-node throughput and scaling runs over
`-threads 16/32/64/128`; `-machine socket` keeps one CPU and its memory (pin
with `numactl`); `-machine desktop` with the same `-threads` shows how a small
grid scales on many threads.

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
steps 2–15. Timers in the profile are
inclusive and must not be added together. For a quick wiring check use
`-cells 8 8 6 -ppc 8`.
