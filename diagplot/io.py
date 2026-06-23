import os
import numpy as np
from .transforms import (
    yee_grid_to_common,
    cartesian_to_cylindrical,
)


def read_binary_2d(filepath, nfields):
    with open(filepath, "rb") as f:
        header = np.fromfile(f, dtype=np.float32, count=2)
        nx, ny = int(header[0]), int(header[1])
        data = np.fromfile(f, dtype=np.float32)
    arrays = []
    for i in range(nfields):
        offset = i * nx * ny
        arrays.append(data[offset:offset + nx * ny].reshape(nx, ny))
    return arrays, nx, ny


class DataCache:
    def __init__(self, workdir, plane, slice_num, species_list):
        self.workdir = workdir
        self.plane = plane
        self.slice_num = slice_num
        self.species_list = species_list
        self._cache = {}
        self._rmap = None
        self._nx = None
        self._ny = None

    def clear(self):
        self._cache.clear()

    @property
    def nx(self):
        return self._nx

    @property
    def ny(self):
        return self._ny

    def _field_path(self, base, timestep):
        return os.path.join(
            self.workdir, "Fields", "Diag2D",
            f"{base}_Plane{self.plane}_{self.slice_num}_{timestep}",
        )

    def _particle_path(self, species, base, timestep):
        return os.path.join(
            self.workdir, "Particles", species, "Diag2D",
            f"{base}_Plane{self.plane}_{self.slice_num}_{timestep}",
        )

    def _load_field_e(self, ts):
        key = ("FieldE", ts)
        if key not in self._cache:
            path = self._field_path("FieldE", ts)
            arrays, nx, ny = read_binary_2d(path, 3)
            self._nx, self._ny = nx, ny
            self._cache[key] = {"Ex": arrays[0], "Ey": arrays[1], "Ez": arrays[2]}
        return self._cache[key]

    def _load_field_b(self, ts):
        key = ("FieldB", ts)
        if key not in self._cache:
            path = self._field_path("FieldB", ts)
            arrays, nx, ny = read_binary_2d(path, 3)
            self._nx, self._ny = nx, ny
            self._cache[key] = {"Bx": arrays[0], "By": arrays[1], "Bz": arrays[2]}
        return self._cache[key]

    def _load_cylindrical_e(self, ts):
        key = ("CylE", ts)
        if key not in self._cache:
            fe = self._load_field_e(ts)
            Ex_c, Ey_c = yee_grid_to_common(fe["Ex"], fe["Ey"])
            Er, Ephi = cartesian_to_cylindrical(Ex_c, Ey_c, self._nx, self._ny)
            self._cache[key] = {"Er": Er, "Ephi": Ephi}
        return self._cache[key]

    def _load_cylindrical_b(self, ts):
        key = ("CylB", ts)
        if key not in self._cache:
            fb = self._load_field_b(ts)
            Bx_c, By_c = yee_grid_to_common(fb["Bx"], fb["By"])
            Br, Bphi = cartesian_to_cylindrical(Bx_c, By_c, self._nx, self._ny)
            self._cache[key] = {"Br": Br, "Bphi": Bphi}
        return self._cache[key]

    def _load_current(self, species, ts):
        key = ("Current", species, ts)
        if key not in self._cache:
            path = self._particle_path(species, "Current", ts)
            arrays, nx, ny = read_binary_2d(path, 3)
            self._nx, self._ny = nx, ny
            self._cache[key] = {"Jx": arrays[0], "Jy": arrays[1], "Jz": arrays[2]}
        return self._cache[key]

    def _load_cylindrical_j(self, species, ts):
        key = ("CylJ", species, ts)
        if key not in self._cache:
            j = self._load_current(species, ts)
            Jr, Jphi = cartesian_to_cylindrical(
                j["Jx"], j["Jy"], self._nx, self._ny,
            )
            self._cache[key] = {"Jr": Jr, "Jphi": Jphi}
        return self._cache[key]

    def _load_scalar(self, species, base, ts):
        key = (base, species, ts)
        if key not in self._cache:
            path = self._particle_path(species, base, ts)
            arrays, nx, ny = read_binary_2d(path, 1)
            self._nx, self._ny = nx, ny
            self._cache[key] = arrays[0]
        return self._cache[key]

    def get(self, quantity, species, ts):
        if quantity in ("Ex", "Ey", "Ez"):
            return self._load_field_e(ts)[quantity]
        if quantity in ("Bx", "By", "Bz"):
            return self._load_field_b(ts)[quantity]
        if quantity in ("Er", "Ephi"):
            return self._load_cylindrical_e(ts)[quantity]
        if quantity in ("Br", "Bphi"):
            return self._load_cylindrical_b(ts)[quantity]
        if quantity in ("Jx", "Jy", "Jz"):
            return self._load_current(species, ts)[quantity]
        if quantity in ("Jr", "Jphi"):
            return self._load_cylindrical_j(species, ts)[quantity]
        if quantity == "Density":
            return self._load_scalar(species, "Density", ts)
        if quantity in ("Prr", "Ppp", "Pzz", "Prp", "Prz", "Pzp"):
            return self._load_scalar(species, quantity, ts)
        raise ValueError(f"Unknown quantity: {quantity}")

    def get_rmap(self):
        from .transforms import build_rmap
        if self._rmap is None and self._nx is not None:
            self._rmap = build_rmap(self._nx, self._ny)
        return self._rmap
