"""Guard the physical-node extraction used for the ECCAPFIM growth-rate check."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

HAS_ANALYSIS_PACKAGES = all(importlib.util.find_spec(name) is not None
                            for name in ("numpy", "scipy", "matplotlib"))


@unittest.skipUnless(HAS_ANALYSIS_PACKAGES, "Weibel analysis requires numpy, scipy and matplotlib")
class EccapfimWeibelAnalysisTest(unittest.TestCase):
    def test_physical_nodes_exclude_ghosts_and_periodic_duplicate(self):
        import numpy as np
        from Weibel_eccapfim.check_weibel import read_plane
        for ghost in (1, 2):
            with self.subTest(ghost=ghost), tempfile.TemporaryDirectory() as directory:
                nx, ny = 4, 5
                sx, sy = nx + 2 * ghost + 1, ny + 2 * ghost + 1
                field = np.full((3, sx, sy), -999, dtype=np.float32)
                expected = np.arange(3 * nx * ny, dtype=np.float32).reshape(3, nx, ny)
                field[:, ghost:ghost + nx, ghost:ghost + ny] = expected
                path = Path(directory) / "plane"
                np.concatenate((np.array([sx, sy], dtype=np.float32), field.ravel())).tofile(path)
                np.testing.assert_array_equal(read_plane(path, nx, ny), expected)
                path.write_bytes(path.read_bytes()[:-4])
                with self.assertRaises(ValueError):
                    read_plane(path, nx, ny)

    def test_growth_fit_measures_amplitude_rate_in_prescribed_window(self):
        import numpy as np
        from Weibel_eccapfim.check_weibel import fit_growth
        times = np.arange(0, 181, 6)
        amplitude = 0.0001 * np.exp(0.023 * times)
        fit = fit_growth(times, amplitude)
        self.assertAlmostEqual(fit["gamma"], 0.023, places=12)
        self.assertAlmostEqual(fit["r_squared"], 1, places=12)
        with self.assertRaises(ValueError):
            fit_growth(times, np.zeros_like(times))


if __name__ == "__main__":
    unittest.main()
