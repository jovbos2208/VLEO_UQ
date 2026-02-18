import math
import unittest

import numpy as np

from vleo_uq import (
    AeroAdapter,
    DeterministicPropagator,
    EnvInputs,
    PropagatorConfig,
    VehicleParams,
    default_geometry,
)


class TestDeterministicPropagator(unittest.TestCase):
    def test_two_body_conservation(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)

        vehicle = VehicleParams()
        config = PropagatorConfig()
        prop = DeterministicPropagator(aero, vehicle, config)

        mu = config.mu_earth_m3_s2
        r0 = np.array([7000e3, 0.0, 0.0])
        v0 = np.array([0.0, math.sqrt(mu / np.linalg.norm(r0)), 0.0])

        x0 = np.zeros(prop.state_size)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0  # unit quaternion w

        period = 2.0 * math.pi * math.sqrt(np.linalg.norm(r0) ** 3 / mu)
        t_grid = np.linspace(0.0, period, 201)

        env = []
        for _ in t_grid:
            e = EnvInputs()
            e.density = 0.0
            e.temperature_K = 1000.0
            e.particles_mass_kg = 28.0 * 1.6605390689252e-27
            e.wind_I = np.zeros(3)
            env.append(e)

        X = prop.propagate(x0, t_grid, env)
        rf = X[-1, 0:3]
        vf = X[-1, 3:6]

        energy0 = 0.5 * np.dot(v0, v0) - mu / np.linalg.norm(r0)
        energyf = 0.5 * np.dot(vf, vf) - mu / np.linalg.norm(rf)
        rel_energy_err = abs((energyf - energy0) / energy0)

        h0 = np.linalg.norm(np.cross(r0, v0))
        hf = np.linalg.norm(np.cross(rf, vf))
        rel_h_err = abs((hf - h0) / h0)

        self.assertLess(rel_energy_err, 1e-4)
        self.assertLess(rel_h_err, 1e-4)


if __name__ == "__main__":
    unittest.main()
