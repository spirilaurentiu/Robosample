"""units.py

Centralised unit-conversion constants.

All physical quantities in the package are converted to SI-adjacent units on
construction (lengths: nm, energies: kJ/mol, angles: rad).
"""

from __future__ import annotations

import parmed as pmd

KCAL_TO_KJ: float = pmd.unit.kilocalories_per_mole.conversion_factor_to(
    pmd.unit.kilojoules_per_mole
)
ANG_TO_NM: float = pmd.unit.angstrom.conversion_factor_to(pmd.unit.nanometer)
DEG_TO_RAD: float = pmd.unit.degree.conversion_factor_to(pmd.unit.radian)

# Converts the AMBER r_min (equilibrium pair distance) to the LJ sigma:
#   sigma = r_min / 2^(1/6)   =>   SIGMA_SCALE = 2^(-1/6) ~= 0.8909
SIGMA_SCALE: float = 2.0 ** (-1.0 / 6.0)
