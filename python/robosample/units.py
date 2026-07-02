"""units.py

Centralised unit-conversion constants.

All physical quantities in the package are converted to SI-adjacent units on
construction (lengths: nm, energies: kJ/mol, angles: rad).
"""

from __future__ import annotations

import math

# Exact physical/definitional constants (fast-loader Step 4b: this module is
# imported by the parmed-free load path -- amber_loader.py,
# molecule_prototype.py, prmtop_reader.py -- and must not import parmed
# itself). Verified bit-for-bit identical to parmed.unit's conversion
# factors (docs/specs/fast-amber-loader.md Step 4b):
#   1 thermochemical calorie = 4.184 J, exactly (SI definition).
#   1 Angstrom = 0.1 nm, exactly (SI definition).
#   1 degree = pi/180 rad, exactly.
KCAL_TO_KJ: float = 4.184
ANG_TO_NM: float = 0.1
DEG_TO_RAD: float = math.pi / 180.0

# Converts the AMBER r_min (equilibrium pair distance) to the LJ sigma:
#   sigma = r_min / 2^(1/6)   =>   SIGMA_SCALE = 2^(-1/6) ~= 0.8909
SIGMA_SCALE: float = 2.0 ** (-1.0 / 6.0)
