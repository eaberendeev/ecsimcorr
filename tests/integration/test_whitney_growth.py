import unittest
from unittest.mock import patch

import numpy as np

from Weibel_ecsim.analyze_whitney import fit_growth
import weibel_check


class WhitneyGrowthFitTests(unittest.TestCase):
    def test_amplitude_and_energy_recover_the_same_physical_increment(self):
        # A growing linear mode has A~exp(gamma*t), W~A^2. The factor of two
        # is physical, independent of the particular least-squares routine.
        time = np.arange(0.0, 300.1, 7.2)
        gamma = 0.0232
        amplitude = 3e-6 * np.exp(gamma * time)
        energy = 0.5 * amplitude ** 2
        amplitude_fit = fit_growth(time, amplitude, (72, 270))
        energy_fit = fit_growth(time, energy, (72, 270), factor=2)
        self.assertAlmostEqual(amplitude_fit["gamma"], gamma, places=13)
        self.assertAlmostEqual(energy_fit["gamma"], gamma, places=13)
        self.assertAlmostEqual(amplitude_fit["r_squared"], 1, places=13)

    def test_short_or_zero_signal_is_not_a_growth_measurement(self):
        time = np.arange(0.0, 40.0, 7.2)
        for values in (np.zeros_like(time), np.ones_like(time)):
            with self.assertRaises(ValueError):
                fit_growth(time, values, (72, 270))

    def test_snapshot_times_follow_the_run_timestep(self):
        plane = [np.zeros((30, 30)) for _ in range(3)]
        with patch.object(weibel_check, "instantaneous_snapshots", return_value=[0, 1, 2]), \
                patch.object(weibel_check, "read_field_plane", return_value=plane):
            for dt in (0.75, 1.2, 1.5):
                time, _, _ = weibel_check.mode_amplitudes("unused", [(1, 1)], interval=12, dt=dt)
                np.testing.assert_allclose(time, np.array([0, 12, 24]) * dt)


if __name__ == "__main__":
    unittest.main()
