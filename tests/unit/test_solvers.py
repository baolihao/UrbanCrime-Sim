import numpy as np

from urbancrime.solvers.partitioned import relative_l2_change


def test_relative_l2_change_handles_zero_reference() -> None:
    value = relative_l2_change(np.array([1.0, 0.0]), np.zeros(2), absolute_floor=0.5)
    assert value == 2.0
