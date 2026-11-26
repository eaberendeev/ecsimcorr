"""Read a beren3d run directory: configs, 2D slice catalog and text tables."""
import glob
import json
import os
import re
from dataclasses import dataclass

import numpy as np

GHOST = 1  # BEREN_GHOST_CELLS: storage index i is node i - GHOST
AXES = "xyz"
IN_PLANE = {"x": (1, 2), "y": (0, 2), "z": (0, 1)}

# 2D file basename -> output type that writes it (intervals live there).
OUTPUT_TYPE = {
    "FieldE": "fields_2d", "FieldB": "fields_2d", "FieldEn": "fields_2d", "FieldBn": "fields_2d",
    "FieldE_external": "fields_2d", "FieldBInit": "fields_2d",
    "ChargeDensity": "charge_density_2d", "TotalCurrent": "total_current_2d",
    "Density": "density_2d", "Current": "current_2d",
    **{p: "pressure_2d" for p in ("Prr", "Ppp", "Pzz", "Prp", "Prz", "Pzp")},
}
SLICE_NAME = re.compile(r"^(?P<base>\w+?)_Plane(?P<axis>[XYZ])_(?P<pos>\d{3,})_(?P<avg>avg)?(?P<seq>\d{4,})$")


def decimals_of(dt):
    """Decimal places of dt, so step * dt prints without rounding noise."""
    for decimals in range(13):
        scaled = dt * 10 ** decimals
        if abs(scaled - round(scaled)) <= 1e-9 * max(1.0, scaled):
            return decimals
    return 12


@dataclass(frozen=True)
class Plane:
    axis: str   # normal axis: x, y or z
    index: int  # storage index in the file name

    @property
    def tag(self):
        return f"{self.axis}{self.index}"


@dataclass(frozen=True)
class Frame:
    seq: int
    first: int  # first step covered
    last: int   # last step covered (== first for an instantaneous frame)

    @property
    def averaged(self):
        return self.last != self.first


class Run:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        with open(self._p("system_config.json"), encoding="utf-8") as f:
            self.config = json.load(f)
        species_path = self._p("particles_config.json")
        species = []
        if os.path.isfile(species_path):
            with open(species_path, encoding="utf-8") as f:
                species = [s["Name"] for s in json.load(f).get("particles", [])]
        self.species = species
        c = self.config
        self.dt = float(c["Dt"])
        self.cells = [int(c[f"NumCells{a}"]) for a in "XYZ"]
        self.step = [float(c[f"D{a}"]) for a in "xyz"]
        cyl = c.get("CylinderDomain")
        self.center = cyl["center"] if cyl else [0.5 * n * d for n, d in zip(self.cells[:2], self.step[:2])]
        self.time_decimals = decimals_of(self.dt)
        self.outputs = c.get("diagnostics", {}).get("outputs", [])
        self.slices = self._catalog()

    def _p(self, *parts):
        return os.path.join(self.path, *parts)

    # ---- time -----------------------------------------------------------
    def time(self, step):
        return step * self.dt

    def fmt_time(self, step):
        return f"{self.time(step):.{self.time_decimals}f}"

    def frame_label(self, frame):
        if frame.averaged:
            return (f"⟨t = {self.fmt_time(frame.first)} … {self.fmt_time(frame.last)}⟩"
                    f"   steps {frame.first}…{frame.last}")
        return f"t = {self.fmt_time(frame.first)}   step {frame.first}"

    def _output(self, base):
        kind = OUTPUT_TYPE.get(base)
        return next((o for o in self.outputs if o.get("type") == kind), {})

    def _frame(self, base, seq, averaged):
        out = self._output(base)
        if averaged:
            n = int(out.get("average_steps", 1))
            first = int(self.config.get("StartTimeStep", 0)) + seq * n
            return Frame(seq, first, min(first + n - 1, int(self.config.get("LastTimestep", first + n - 1))))
        interval = int(out.get("interval", self.config.get("diagnostics", {}).get("TimeStepDelayDiag2D", 1)))
        return Frame(seq, seq * interval, seq * interval)

    # ---- 2D slices ------------------------------------------------------
    def _catalog(self):
        """{(owner, base, plane, averaged): {seq: path}}; owner is a species or None."""
        dirs = [(None, self._p("Fields", "Diag2D"))]
        dirs += [(s, self._p("Particles", s, "Diag2D")) for s in self.species]
        catalog = {}
        for owner, directory in dirs:
            if not os.path.isdir(directory):
                continue
            for name in os.listdir(directory):
                m = SLICE_NAME.match(name)
                if not m:
                    continue
                key = (owner, m["base"], Plane(m["axis"].lower(), int(m["pos"])), bool(m["avg"]))
                catalog.setdefault(key, {})[int(m["seq"])] = os.path.join(directory, name)
        return catalog

    def planes(self):
        return sorted({key[2] for key in self.slices}, key=lambda p: ("zyx".index(p.axis), p.index))

    def frames(self, owner, base, plane, averaged):
        seqs = self.slices.get((owner, base, plane, averaged), {})
        return [self._frame(base, s, averaged) for s in sorted(seqs)]

    def slice_path(self, owner, base, plane, averaged, seq):
        return self.slices[(owner, base, plane, averaged)][seq]

    def plane_coordinate(self, plane):
        a = AXES.index(plane.axis)
        return (plane.index - GHOST) * self.step[a]

    def plane_axes(self, plane):
        """Node coordinates of the two in-plane axes for the physical nodes."""
        return [np.arange(self.cells[a] + 1) * self.step[a] for a in IN_PLANE[plane.axis]]

    # ---- text tables ----------------------------------------------------
    def tables(self):
        """{column: (steps, values)} merged from energy, boundary, continuity and probes."""
        columns = {}
        for name in ("energy.txt", "boundary.txt", "continuity.txt"):
            path = self._p(name)
            if os.path.isfile(path):
                header, data = _read_table(path)
                if data.size == 0:
                    continue
                steps = (data[:, header.index("Step")] if "Step" in header
                         else np.rint(data[:, header.index("Time")] / self.dt))
                for i, column in enumerate(header):
                    if column not in ("Step", "Time"):
                        columns[column] = (steps, data[:, i])
        columns.update(self.probes())
        return columns

    def probes(self):
        path = self._p("Probes", "probes.txt")
        if not os.path.isfile(path):
            return {}
        header, data = _read_table(path)
        if data.size == 0:
            return {}
        # probes.txt names its step column "Time".
        points = {tuple(row): i for i, row in reversed(list(enumerate(map(tuple, data[:, 1:4]))))}
        columns = {}
        for p, point in enumerate(sorted(points, key=points.get)):
            rows = np.all(data[:, 1:4] == point, axis=1)
            for i, column in enumerate(header[4:], start=4):
                columns[f"{column}@p{p}"] = (data[rows, 0], data[rows, i])
        return columns

    def probe_points(self):
        path = self._p("Probes", "probes.txt")
        if not os.path.isfile(path):
            return []
        _, data = _read_table(path)
        seen = []
        for row in map(tuple, data[:, 1:4]) if data.size else []:
            if row not in seen:
                seen.append(row)
        return seen

    def spectra(self, species):
        """[(step, energies, counts)] sorted by step."""
        result = []
        for path in glob.glob(self._p("Particles", species, "energy_spectrum_*.txt")):
            step = int(re.search(r"(\d+)\.txt$", path)[1])
            data = np.loadtxt(path, ndmin=2)
            if data.size:
                result.append((step, data[:, 0], data[:, 1]))
        return sorted(result, key=lambda r: r[0])


def _read_table(path):
    with open(path, encoding="utf-8") as f:
        header = f.readline().lstrip("#").split()
    data = np.loadtxt(path, skiprows=1, ndmin=2)
    return header, data[:, :len(header)] if data.size else data


def read_slice(path):
    """Return the 3 component slots of a 2D slice as an array (3, n1, n2)."""
    raw = np.fromfile(path, dtype=np.float32)
    n1, n2 = int(raw[0]), int(raw[1])
    return raw[2:2 + 3 * n1 * n2].reshape(3, n1, n2).astype(np.float64)


def select(frames, time_range, every, run):
    """Frames whose time lies in [a, b] ('a:b', 'a:', ':b', 'a' = nearest, 'last'), then every k-th."""
    if not frames:
        return frames
    if time_range == "last":
        frames = frames[-1:]
    elif time_range:
        if ":" in time_range:
            a, b = (float(v) if v else None for v in time_range.split(":", 1))
            frames = [f for f in frames if (a is None or run.time(f.last) >= a - 1e-12)
                      and (b is None or run.time(f.first) <= b + 1e-12)]
        else:
            t = float(time_range)
            frames = [min(frames, key=lambda f: abs(0.5 * run.time(f.first + f.last) - t))]
    return frames[::max(1, every)]
