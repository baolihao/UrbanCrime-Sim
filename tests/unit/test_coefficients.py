import numpy as np

from urbancrime.coefficients import analytic_profile
from urbancrime.config import CoefficientSpec


def test_diagonal_highway_respects_quadrant_cutoff() -> None:
    profile = analytic_profile(
        CoefficientSpec.from_mapping(
            {
                "kind": "analytic",
                "profile": "diagonal_highway",
                "parameters": {
                    "background": 1.0,
                    "amplitude": 2.0,
                    "sharpness": 20.0,
                    "intercept": 16.0,
                    "x_min": 8.0,
                    "y_max": 8.0,
                },
            }
        )
    )
    values = profile(np.asarray([[9.0, 7.0], [7.0, 9.0]]))
    np.testing.assert_allclose(values, [3.0, 1.0])


def test_piecewise_i90_profile_uses_both_segments() -> None:
    profile = analytic_profile(
        CoefficientSpec.from_mapping(
            {
                "kind": "analytic",
                "profile": "piecewise_linear_ridge",
                "parameters": {
                    "background": 1.0,
                    "amplitude": 1.0,
                    "sharpness": 20.0,
                    "switch_x": 9.5,
                    "left_line": [0.0, -1.0, 20.8],
                    "right_line": [-16.0 / 13.0, -1.0, 32.4923],
                },
            }
        )
    )
    right_y = 32.4923 - 16.0 * 12.0 / 13.0
    values = profile(np.asarray([[5.0, 12.0], [20.8, right_y]]))
    np.testing.assert_allclose(values, [2.0, 2.0])
