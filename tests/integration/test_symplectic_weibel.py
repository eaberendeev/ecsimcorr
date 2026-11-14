import math
import unittest

try:
    from Weibel_symplectic.run import gamma_theory
    from weibel_check import configure_snapshot_times
except ImportError:
    gamma_theory = None


@unittest.skipIf(gamma_theory is None, "Weibel theory requires optional numpy/scipy")
class WeibelTheoryTest(unittest.TestCase):
    def test_isotropic_and_short_wavelength_cutoff(self):
        self.assertEqual(gamma_theory(1.16, 0.03, 1), 0)
        self.assertEqual(gamma_theory(3.0, 0.03, 10), 0)
        self.assertGreater(gamma_theory(1.16, 0.03, 10), 0)

    def test_full_benchmark_preserves_particles_and_physical_snapshot_times(self):
        from Weibel_symplectic import gen_config as config
        source = config.reference
        self.assertEqual(config.particles_config, source.particles_config)
        for key in ("NumCellsX", "NumCellsY", "NumCellsZ", "Dx", "Dy", "Dz", "n0", "Boundary_conditions"):
            self.assertEqual(config.system_config[key], source.system_config[key], key)
        self.assertEqual(config.LastTimestep * config.Dt, source.LastTime)
        self.assertAlmostEqual(configure_snapshot_times(config.system_config),
                               configure_snapshot_times(source.system_config))

    def test_cold_parallel_limit_against_fluid_polynomial(self):
        # For sigma_x -> 0 at fixed sigma_y the transverse susceptibility
        # gives gamma^4+(1+k^2)*gamma^2-k^2*sigma_y^2=0.
        # This independent limiting equation checks signs and sqrt(2) units.
        k, sigma_hot, anisotropy = 1.16, 0.1, 1e6
        a = 1 + k*k
        expected = math.sqrt(2*k*k*sigma_hot**2 / (math.sqrt(a*a+4*k*k*sigma_hot**2)+a))
        actual = gamma_theory(k, sigma_hot / math.sqrt(anisotropy), anisotropy)
        self.assertLess(abs(actual / expected - 1), 1e-3)


if __name__ == "__main__":
    unittest.main()
