"""
Utilities for reading AMBER prmtop topology files and detecting NBFIX
non-standard Lennard-Jones pair interactions.

References
----------
Case, D.A. et al. AMBER 2023, University of California, San Francisco.
AMBER file format specification:
  https://ambermd.org/FileFormats.php#topology
"""

import logging
import os
import re
from typing import Any

import numpy as np
import parmed as pmd

logger = logging.getLogger(__name__)

# Compiled once at import time; matches Fortran format strings such as
# "20a4", "5E16.8", "12I6".
_FORMAT_RE = re.compile(r"(\d+)\(?([a-zA-Z]+)(\d+)\.?(\d*)\)?")


def parse_prmtop(prmtop_path: str | os.PathLike[str]) -> dict[str, Any]:
    """Parse an AMBER prmtop (parameter/topology) file into numpy arrays.

    Reads the binary-free, Fortran-formatted ASCII topology file produced by
    AMBER's *tleap* / *parmed* utilities.  All numeric sections are returned as
    ``numpy`` arrays; string sections remain as lists of ``str``.

    Atomic partial charges are stored in the file multiplied by the AMBER
    internal factor 18.2223 (√(332.0636), arising from the conversion of
    kcal mol⁻¹ Å to e²).  This function divides by that factor so that
    ``raw_data["CHARGE"]`` contains charges in units of the proton charge *e*.

    Parameters
    ----------
    prmtop_path : str or os.PathLike[str]
        Path to the ``*.prmtop`` or ``*.top`` file.

    Returns
    -------
    dict with keys:

    ``"version"`` : str or None
        Content of the ``%VERSION`` line.
    ``"flags"`` : list of str
        Ordered list of ``%FLAG`` names found in the file.
    ``"raw_data"`` : dict[str, np.ndarray | list]
        Parsed data arrays keyed by flag name.
    ``"raw_format"`` : dict[str, tuple]
        Fortran format descriptors keyed by flag name, each a 5-tuple
        ``(fmt_str, num_items, item_type, field_width, precision)``.
    ``"chamber"`` : bool
        ``True`` when the file contains CHARMM-style (CHAMBER) sections,
        detected by the presence of the ``CTITLE`` flag.

    Raises
    ------
    FileNotFoundError
        If *prmtop_path* does not exist.
    ValueError
        If a data line is encountered before any ``%FLAG`` directive.
    """
    flags: list[str] = []
    raw_data: dict[str, Any] = {}
    raw_format: dict[str, tuple] = {}
    prmtop_version: str | None = None

    with open(prmtop_path, "r") as fh:
        lines = [line.rstrip("\n") for line in fh]

    for lineno, line in enumerate(lines, start=1):
        if not line:
            continue

        if line.startswith("%"):
            if line.startswith("%VERSION"):
                _, prmtop_version = line.split(None, 1)
            elif line.startswith("%FLAG"):
                _, flag = line.split(None, 1)
                flag = flag.strip()
                flags.append(flag)
                raw_data[flag] = []
                logger.debug("Encountered flag %s", flag)
            elif line.startswith("%FORMAT"):
                fmt_line = line[line.index("(") + 1 : line.index(")")]
                m = _FORMAT_RE.search(fmt_line)
                if m:
                    raw_format[flags[-1]] = (
                        fmt_line,
                        int(m.group(1)),
                        m.group(2),
                        int(m.group(3)),
                        m.group(4),
                    )
                elif "," in fmt_line:
                    # Compound Fortran format (e.g. "i2,a78" in CHAMBER's
                    # FORCE_FIELD_TYPE).  The sub-formats describe fields on a
                    # single line; storing the whole line as a string is correct
                    # for all compound-format flags encountered in practice.
                    logger.debug(
                        "Flag %s: compound FORMAT '%s'; storing lines as strings.",
                        flags[-1],
                        fmt_line,
                    )
                    raw_format[flags[-1]] = (fmt_line, 1, "a", 80, "")
                else:
                    logger.warning(
                        "Line %d: unrecognised FORMAT '%s' for flag %s; "
                        "defaulting to 80-char string.",
                        lineno,
                        fmt_line,
                        flags[-1] if flags else "<none>",
                    )
                    raw_format[flags[-1]] = (fmt_line, 1, "a", 80, "")
            continue

        # ------------------------------------------------------------------ #
        # Data line
        # ------------------------------------------------------------------ #
        if not flags:
            raise ValueError(
                f"Line {lineno}: data encountered before any %FLAG directive."
            )

        flag = flags[-1]
        _fmt, _num_items, item_type, i_length, _item_prec = raw_format[flag]

        # TITLE is a free-format single line; capture it verbatim.
        if flag == "TITLE" and not raw_data[flag]:
            raw_data[flag] = [line]
            continue

        # Split the fixed-width Fortran record into individual fields.
        # Plain string slicing is sufficient and faster than the numpy
        # frombuffer approach for short ASCII lines.
        chunks = [line[i : i + i_length] for i in range(0, len(line), i_length)]
        items = [c.strip() for c in chunks if c.strip()]

        if not items:
            continue

        # Type conversion
        itype_upper = item_type.upper()
        if itype_upper == "A":
            raw_data[flag].extend(items)
        elif itype_upper == "I":
            raw_data[flag].extend(np.array(items, dtype=np.int64))
        elif itype_upper in ("E", "F", "D"):
            # Fortran may use 'D' exponent notation (double precision literals)
            raw_data[flag].extend(
                np.array(
                    [float(x.replace("D", "E")) for x in items],
                    dtype=np.float64,
                )
            )
        else:
            logger.warning(
                "Flag %s: unknown item type '%s'; storing as strings.",
                flag,
                item_type,
            )
            raw_data[flag].extend(items)

    # Consolidate per-flag lists into arrays
    for flag, values in raw_data.items():
        if not isinstance(values, list) or not values:
            continue
        first = values[0]
        if isinstance(first, (int, np.integer)):
            raw_data[flag] = np.array(values, dtype=np.int64)
        elif isinstance(first, (float, np.floating)):
            raw_data[flag] = np.array(values, dtype=np.float64)

    # Rescale charges from AMBER internal units to proton-charge units (e).
    if "CHARGE" in raw_data:
        raw_data["CHARGE"] = raw_data["CHARGE"] / 18.2223
    else:
        logger.warning(
            "No CHARGE section found in '%s'; partial charges not available.",
            prmtop_path,
        )

    chamber_style = "CTITLE" in flags
    if chamber_style:
        logger.debug("CHAMBER-style prmtop detected.")

    return {
        "version": prmtop_version,
        "flags": flags,
        "raw_data": raw_data,
        "raw_format": raw_format,
        "chamber": chamber_style,
    }


def has_nbfix_fast(
    nb_indices: np.ndarray,
    num_types: int,
    acoef: np.ndarray,
    bcoef: np.ndarray,
) -> bool:
    """Detect non-standard (NBFIX) Lennard-Jones pair interactions.

    In AMBER, the Lennard-Jones potential between atom types *i* and *j* is

    .. math::

        U_{ij}(r) = \\frac{A_{ij}}{r^{12}} - \\frac{B_{ij}}{r^6}

    Under the standard Lorentz-Berthelot combining rules the cross-pair
    coefficients are uniquely determined by the diagonal (self-interaction)
    coefficients:

    .. math::

        R_{\\min,i} = \\left(\\frac{2A_{ii}}{B_{ii}}\\right)^{1/6}, \\quad
        \\varepsilon_i = \\frac{B_{ii}^2}{4 A_{ii}}

    .. math::

        A_{ij}^{\\text{comb}} = \\varepsilon_{ij}\\,R_{\\min,ij}^{12}, \\quad
        B_{ij}^{\\text{comb}} = 2\\,\\varepsilon_{ij}\\,R_{\\min,ij}^{6}

    with :math:`R_{\\min,ij} = R_{\\min,i}/2 + R_{\\min,j}/2` and
    :math:`\\varepsilon_{ij} = \\sqrt{\\varepsilon_i \\varepsilon_j}`.

    A pair is flagged as NBFIX when its stored *A*/*B* coefficients deviate
    from the combining-rule prediction by more than a relative tolerance of
    10⁻⁶.

    Parameters
    ----------
    nb_indices : array_like, shape (num_types * num_types,)
        ``NONBONDED_PARM_INDEX`` array from the prmtop file (1-based).
        Internally reshaped to ``(num_types, num_types)``.
    num_types : int
        Number of distinct LJ atom types (``NTYPES`` in the prmtop).
    acoef : np.ndarray, shape (n_pairs,)
        ``LENNARD_JONES_ACOEF`` array from the prmtop (r¹² coefficients).
    bcoef : np.ndarray, shape (n_pairs,)
        ``LENNARD_JONES_BCOEF`` array from the prmtop (r⁶ coefficients).

    Returns
    -------
    bool
        ``True`` if any atom-type pair carries NBFIX coefficients that
        deviate from the Lorentz-Berthelot combining rules;
        ``False`` otherwise.

    Notes
    -----
    Pairs where either *A* or *B* is exactly zero (dummy atoms, virtual
    sites, etc.) are handled separately: a zero coefficient is only
    considered anomalous when the combining-rule prediction is non-zero,
    or when the sibling coefficient (*B* for a zero-*A* pair, and vice
    versa) is itself non-zero.

    The function is silent with respect to floating-point exceptions:
    division by zero arising from zero self-interaction coefficients
    produces ``nan`` / ``inf`` internally, which are mapped to zero via
    ``np.where`` before any comparison.
    """
    nb_indices = (
        np.asarray(nb_indices, dtype=np.int64).reshape(num_types, num_types) - 1
    )

    # ------------------------------------------------------------------ #
    # Diagonal (self) LJ parameters → per-type Rmin and epsilon
    # ------------------------------------------------------------------ #
    diag_idx = nb_indices.diagonal()
    A_ii = acoef[diag_idx]
    B_ii = bcoef[diag_idx]

    with np.errstate(divide="ignore", invalid="ignore"):
        # rmin_i = (2 A_ii / B_ii)^(1/6)  [full Rmin, not Rmin/2]
        rmin = (2.0 * A_ii / B_ii) ** (1.0 / 6.0)
        eps_i = 0.25 * B_ii**2 / A_ii

    # Zero out non-finite values (dummy/virtual-site types with A=B=0)
    ri = np.where(np.isfinite(rmin), rmin / 2.0, 0.0)
    eps_i = np.where(np.isfinite(eps_i), eps_i, 0.0)

    # ------------------------------------------------------------------ #
    # Combining-rule predictions for all type pairs
    # ------------------------------------------------------------------ #
    expected_R = ri[:, None] + ri[None, :]  # Rmin,ij
    expected_E = np.sqrt(eps_i[:, None] * eps_i[None, :])  # eps_ij

    calc_A = expected_E * expected_R**12
    calc_B = 2.0 * expected_E * expected_R**6

    # ------------------------------------------------------------------ #
    # Retrieve actual stored coefficients (0-based index via nb_indices)
    # ------------------------------------------------------------------ #
    mask = nb_indices >= 0  # valid pairs only

    actual_A = np.zeros((num_types, num_types))
    actual_B = np.zeros((num_types, num_types))
    actual_A[mask] = acoef[nb_indices[mask]]
    actual_B[mask] = bcoef[nb_indices[mask]]

    # ------------------------------------------------------------------ #
    # Check 1: inconsistent zeros
    # A zero coefficient is suspicious when the other member of the pair
    # is non-zero, or when the combining-rule prediction is non-zero.
    # ------------------------------------------------------------------ #
    zero_mask = (actual_A == 0) | (actual_B == 0)
    bad_zero = zero_mask & (
        (actual_A != 0) | (actual_B != 0) | ((expected_E != 0) & (expected_R != 0))
    )
    if np.any(bad_zero & mask):
        logger.debug("NBFIX detected via zero-coefficient inconsistency.")
        return True

    # ------------------------------------------------------------------ #
    # Check 2: relative deviation from combining rules
    # Pairs where actual_A == actual_B == 0 (and calc values are also 0)
    # yield 0/0 = nan; suppressing the warning is intentional — nan
    # comparisons evaluate to False, correctly skipping those pairs.
    # ------------------------------------------------------------------ #
    with np.errstate(divide="ignore", invalid="ignore"):
        bad_A = np.abs((actual_A - calc_A) / actual_A) > 1e-6
        bad_B = np.abs((actual_B - calc_B) / actual_B) > 1e-6

    if np.any((bad_A | bad_B) & mask):
        logger.debug("NBFIX detected via combining-rule deviation.")
        return True

    return False


def load_lj_coefs(
    parm_data: dict,
    num_types: int,
    ene_conv: float = pmd.unit.kilocalories_per_mole.conversion_factor_to(
        pmd.unit.kilojoules_per_mole
    ),
    length_conv: float = pmd.unit.angstroms.conversion_factor_to(pmd.unit.nanometers),
) -> tuple[np.ndarray, np.ndarray]:
    """Convert and return Lennard-Jones coefficients in SI-adjacent units.

    ...

    Parameters
    ----------
    parm_data : dict
        ``raw_data`` dict from :func:`parse_prmtop_numpy`.  Must contain
        ``NONBONDED_PARM_INDEX``, ``LENNARD_JONES_ACOEF``, and
        ``LENNARD_JONES_BCOEF``.
    num_types : int
        Number of distinct LJ atom types (``NTYPES`` pointer).
    ene_conv : float, optional
        Multiplicative factor converting kcal mol⁻¹ to the target energy
        unit.  Defaults to the kcal mol⁻¹ → kJ mol⁻¹ factor as reported
        by :mod:`parmed`.  Override only when targeting a non-SI energy
        unit.
    length_conv : float, optional
        Multiplicative factor converting Å to the target length unit.
        Defaults to the Å → nm factor as reported by :mod:`parmed`.
        Override only when targeting a non-SI length unit.

    Returns
    -------
    ...
    """
    """Convert and return Lennard-Jones coefficients in SI-adjacent units.

    Reads the raw LJ *A* and *B* coefficients from an AMBER prmtop data
    dictionary (as produced by :func:`parse_prmtop_numpy`) and converts them
    from AMBER internal units (kcal mol⁻¹ Å¹²/⁶) to the caller-specified
    unit system via *ene_conv* and *length_conv*.

    The coefficients are stored in a flat array of length
    ``num_types * num_types`` indexed as ``i * num_types + j``, mirroring
    the layout of ``NONBONDED_PARM_INDEX``.

    .. note::
        *A* coefficients are returned as their **square root**
        (``sqrt(A) · sqrt(ene_conv) · length_conv^6``), which is the
        natural form for geometric-mean combining rules.
        *B* coefficients are returned linearly
        (``B · ene_conv · length_conv^6``).

    Parameters
    ----------
    parm_data : dict
        ``raw_data`` dict from :func:`parse_prmtop_numpy`.  Must contain
        ``NONBONDED_PARM_INDEX``, ``LENNARD_JONES_ACOEF``, and
        ``LENNARD_JONES_BCOEF``.
    num_types : int
        Number of distinct LJ atom types (``NTYPES`` pointer).
    ene_conv : float
        Multiplicative factor converting kcal mol⁻¹ to the target energy
        unit (e.g. ``pmd.unit.kilocalories_per_mole
        .conversion_factor_to(pmd.unit.kilojoules_per_mole)``).
    length_conv : float
        Multiplicative factor converting Å to the target length unit
        (e.g. ``pmd.unit.angstroms
        .conversion_factor_to(pmd.unit.nanometers)``).

    Returns
    -------
    a_coef : np.ndarray, shape (num_types * num_types,)
        Converted square-root *A* coefficients.
    b_coef : np.ndarray, shape (num_types * num_types,)
        Converted *B* coefficients.

    Raises
    ------
    ValueError
        If any entry in ``NONBONDED_PARM_INDEX`` resolves to a negative
        0-based index, indicating a malformed topology file.
    """
    nb_index = np.asarray(parm_data["NONBONDED_PARM_INDEX"])
    acoef = np.asarray(parm_data["LENNARD_JONES_ACOEF"])
    bcoef = np.asarray(parm_data["LENNARD_JONES_BCOEF"])

    # Full num_types*num_types lookup table (row-major, 0-based)
    idx = nb_index[: num_types * num_types] - 1
    if np.any(idx < 0):
        bad = np.argwhere(idx < 0)
        raise ValueError(
            f"Invalid (negative) nonbonded indices at flat positions: {bad.ravel()}"
        )

    afac = np.sqrt(ene_conv) * length_conv**6
    bfac = ene_conv * length_conv**6

    a_coef = np.sqrt(acoef[idx]) * afac
    b_coef = bcoef[idx] * bfac

    return a_coef, b_coef


def load_cmap(
    parm_data: dict,
    parm: pmd.amber.AmberParm,
    ene_conv: float = pmd.unit.kilocalories_per_mole.conversion_factor_to(
        pmd.unit.kilojoules_per_mole
    ),
) -> dict[str, int | list[float] | list[int]]:
    """Parse CMAP correction grids and torsion assignments from an AMBER prmtop.

    CMAP (Correction MAP) defines a 2D potential-energy surface E(phi, psi)
    applied to pairs of consecutive backbone dihedrals.  Each grid is stored
    in the prmtop with phi as the slow (outer) index and psi as the fast
    (inner) index, with the origin at -180 degrees.  OpenMM expects phi as
    the fast (inner) index with the origin at 0 degrees; this function
    performs both the index transposition and the 180-degree cyclic shift.

    All grids must share the same resolution (number of points per axis),
    which is the case for all standard CHARMM force fields (typically 24).

    Parameters
    ----------
    parm_data : dict
        ``raw_data`` dict produced by :func:`parse_prmtop`.  Must contain
        ``CMAP_RESOLUTION`` and ``CMAP_PARAMETER_XX`` entries for CHAMBER
        topologies; returns empty accumulators silently for plain AMBER files
        that carry no CMAP data.
    parm : pmd.amber.AmberParm
        Loaded parmed structure.  Atom global indices are taken from
        ``pmd.Atom.idx`` to avoid maintaining a separate index mapping.
    ene_conv : float, optional
        Multiplicative factor converting kcal mol⁻¹ to the target energy
        unit.  Defaults to the kcal mol⁻¹ -> kJ mol⁻¹ factor from
        :mod:`parmed`.

    Returns
    -------
    dict with keys matching ``SystemTopology`` attribute names:

    ``"cmap_grid_size"`` : int
        Number of grid points along each axis (same for all grids).
    ``"cmap_grid_energy"`` : list[float]
        Concatenated, reordered energy grids in kJ mol⁻¹.  Grid *k* occupies
        positions ``k * size**2`` to ``(k+1) * size**2``, with phi varying
        fastest (OpenMM convention).
    ``"cmap_torsion_a1"`` .. ``"cmap_torsion_a4"`` : list[int]
        Global atom indices for the first (phi) dihedral of each torsion pair.
    ``"cmap_torsion_b1"`` .. ``"cmap_torsion_b4"`` : list[int]
        Global atom indices for the second (psi) dihedral.  Atoms b1-b3
        overlap with a2-a4 (the two dihedrals share a central triplet).
    ``"cmap_torsion_map_index"`` : list[int]
        0-based index into the grid table for each torsion pair.

    Raises
    ------
    ValueError
        If CMAP grids have inconsistent resolutions, or if any torsion
        references a grid index outside the valid range.
    """
    _EMPTY: dict[str, int | list] = {
        "cmap_grid_size": 0,
        "cmap_grid_energy": [],
        "cmap_torsion_a1": [],
        "cmap_torsion_a2": [],
        "cmap_torsion_a3": [],
        "cmap_torsion_a4": [],
        "cmap_torsion_b1": [],
        "cmap_torsion_b2": [],
        "cmap_torsion_b3": [],
        "cmap_torsion_b4": [],
        "cmap_torsion_map_index": [],
    }

    cmap_resolution = parm_data.get("CMAP_RESOLUTION", [])
    num_grids = len(cmap_resolution)
    if num_grids == 0:
        return _EMPTY

    # ------------------------------------------------------------------ #
    # Validate uniform grid resolution
    # ------------------------------------------------------------------ #
    sizes = np.asarray(cmap_resolution, dtype=np.int64)
    if np.any(sizes != sizes[0]):
        raise ValueError(
            f"All CMAP grids must share the same resolution; "
            f"found {sorted(set(sizes.tolist()))}."
        )
    res = int(sizes[0])
    half = res // 2

    # ------------------------------------------------------------------ #
    # Build reorder indices once, reuse for every grid.
    #
    # AMBER layout : old_index = phi_amber * res + psi_amber
    #                origin at -180 deg, psi fastest
    # OpenMM layout: new_index = phi_omm + res * psi_omm
    #                origin at   0 deg, phi fastest
    #
    # Cyclic shift: phi_amber = (phi_omm + half) % res
    # ------------------------------------------------------------------ #
    phi_amber = (np.arange(res) + half) % res
    psi_amber = (np.arange(res) + half) % res
    old_indices = phi_amber[:, None] * res + psi_amber[None, :]  # (res, res)

    # ------------------------------------------------------------------ #
    # Reorder and convert all grids; concatenate into one flat list
    # ------------------------------------------------------------------ #
    grid_energy: list[float] = []
    for i in range(num_grids):
        cmap = np.asarray(parm_data[f"CMAP_PARAMETER_{i + 1:02d}"], dtype=np.float64)
        # Transpose to [psi_omm, phi_omm] then C-order flatten
        # -> new_index = phi_omm + res * psi_omm  (phi fastest)
        grid_energy.extend((cmap[old_indices] * ene_conv).T.ravel().tolist())

    # ------------------------------------------------------------------ #
    # Torsion assignments.
    # CMAP_INDEX layout (6 integers per entry, 1-based prmtop atom indices):
    #   [a1, a2, a3, a4, b4, map_index]
    # Torsion A = (a1, a2, a3, a4)  -- phi dihedral
    # Torsion B = (a2, a3, a4, b4)  -- psi dihedral (shares a2-a4 with A)
    # ------------------------------------------------------------------ #
    cmap_index = parm_data.get("CMAP_INDEX", [])
    if len(cmap_index) == 0:
        return {**_EMPTY, "cmap_grid_size": res, "cmap_grid_energy": grid_energy}

    entries = np.asarray(cmap_index, dtype=np.int64).reshape(-1, 6)
    map_indices = entries[:, 5] - 1  # 0-based grid index

    invalid = (map_indices < 0) | (map_indices >= num_grids)
    if np.any(invalid):
        raise ValueError(
            f"CMAP map indices out of range [1, {num_grids}]: "
            f"{entries[invalid, 5].tolist()}"
        )

    # Resolve 1-based prmtop indices to global atom indices via pmd.Atom.idx
    global_idx = np.array([atom.idx for atom in parm.atoms])
    atoms = global_idx[entries[:, :5] - 1]  # (n_torsions, 5)

    return {
        "cmap_grid_size": res,
        "cmap_grid_energy": grid_energy,
        "cmap_torsion_a1": atoms[:, 0].tolist(),
        "cmap_torsion_a2": atoms[:, 1].tolist(),
        "cmap_torsion_a3": atoms[:, 2].tolist(),
        "cmap_torsion_a4": atoms[:, 3].tolist(),
        "cmap_torsion_b1": atoms[:, 1].tolist(),  # torsion B shares a2-a4
        "cmap_torsion_b2": atoms[:, 2].tolist(),
        "cmap_torsion_b3": atoms[:, 3].tolist(),
        "cmap_torsion_b4": atoms[:, 4].tolist(),
        "cmap_torsion_map_index": map_indices.tolist(),
    }
