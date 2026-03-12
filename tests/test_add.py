# content of test_logic.py
def add(a, b):
    return a + b

def test_add():
    assert add(1, 2) == 3  # Simple assertion

import numpy as np

def test_magnitude_mismatch():
    # These look close, but the difference is 0.01
    actual = np.array([1.0, 2.0, 3.01])
    expected = np.array([1.0, 2.0, 3.00])
    
    # We set a very strict absolute tolerance (atol) of 0.001
    # Since 0.01 > 0.001, this WILL fail.
    np.testing.assert_allclose(actual, expected, atol=1e-3, err_msg="The arrays are too different!")