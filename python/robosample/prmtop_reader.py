"""
prmtop_reader.py

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
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .units import ANG_TO_NM, KCAL_TO_KJ

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
    ene_conv: float = KCAL_TO_KJ,
    length_conv: float = ANG_TO_NM,
) -> tuple[np.ndarray, np.ndarray]:
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
    ene_conv : float, optional
        Multiplicative factor converting kcal mol⁻¹ to the target energy
        unit.  Defaults to ``units.KCAL_TO_KJ`` (kcal mol⁻¹ -> kJ mol⁻¹).
        Override only when targeting a non-SI energy unit.
    length_conv : float, optional
        Multiplicative factor converting Å to the target length unit.
        Defaults to ``units.ANG_TO_NM`` (Å -> nm). Override only when
        targeting a non-SI length unit.

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


# NOTE: an earlier `load_cmap(parm_data, parm, ...)` helper used to live
# here. It required a loaded ParmEd structure (to resolve each atom's global
# index) and was never called anywhere in the codebase --
# ``Context.load_amber`` builds the CMAP grid/torsion arrays itself, inline,
# directly from ``raw_data`` + ``prmtop_to_global_index`` (see the "CMAP
# correction maps and torsions" section of ``context.py``). Removed as dead
# code rather than ported, per the fast-loader Step 4b parmed-free
# requirement (docs/specs/fast-amber-loader.md) -- keeping an unused
# ParmEd-typed function around would be the only remaining parmed dependency
# in this module.


# ====================================================================== #
# 1-4 scaling pairs and explicit non-bonded exclusions
#
# These sections are not exposed on the high-level ParmEd Structure API and
# must be reconstructed from the raw prmtop dihedral pointer list and the
# excluded-atoms list.  Records are returned in *local* (prototype-local,
# 0-based ascending) atom indices; orientation into closest->farthest order
# and remapping to compound indices is the caller's responsibility (see
# MoleculePrototype).
# ====================================================================== #

# Mathematical constant relating the AMBER r_min (equilibrium pair distance) to
# the Lennard-Jones sigma:  sigma = r_min / 2^(1/6).
_SIGMA_SCALE: float = 2.0 ** (-1.0 / 6.0)


@dataclass(slots=True)
class Scaling14Record:
    """One unique 1-4 non-bonded pair, in local atom indices."""

    i_local: int
    l_local: int
    charge_product: float
    """q_i * q_l / scee, proton-charge^2 (SCEE already absorbed)."""
    epsilon: float
    """Combined LJ well depth / scnb, kJ/mol (SCNB already absorbed)."""
    sigma: float
    """Combined LJ sigma, nm."""


@dataclass(slots=True)
class ExclusionRecord:
    """One explicit non-bonded exclusion pair (not already a 1-4 pair), local indices."""

    i_local: int
    j_local: int


@dataclass(slots=True)
class NonbondedTables:
    """1-4 scaling pairs and explicit exclusions parsed from a prmtop."""

    scaling14: list[Scaling14Record] = field(default_factory=list)
    exclusions: list[ExclusionRecord] = field(default_factory=list)


def load_nonbonded_exceptions(
    raw_data: dict,
    atom_indices: np.ndarray,
    ene_conv: float = KCAL_TO_KJ,
    length_conv: float = ANG_TO_NM,
) -> NonbondedTables:
    """Parse 1-4 scaling pairs and explicit exclusions for ONE prototype.

    Companion to :func:`load_lj_coefs` / CMAP handling in ``context.py``: where
    those handle the full LJ table and the CMAP grids, this reconstructs the
    per-pair 1-4 scaled interactions and the explicit exclusion list.

    Part of the fast-loader rewrite (``docs/specs/fast-amber-loader.md`` §3(C)):
    reads directly from the WHOLE-SYSTEM ``raw_data`` (prmtop order, GLOBAL
    atom indices) instead of a ParmEd-sliced-per-prototype ``parm_data`` /
    ``Structure``. The 1-4/exclusion atom-index-space work (which pairs exist)
    is scoped to *atom_indices* (one prototype's/instance's GLOBAL 0-based
    prmtop atom indices, ascending) and returned in *local* (0-based,
    position-within-``atom_indices``) indices. The physical VALUES (charges,
    LJ A/B coefficients, SCEE/SCNB) are looked up directly from the
    GLOBAL/whole-system type tables (``NONBONDED_PARM_INDEX``, LJ (14)
    coefficient arrays, ``SCEE_SCALE_FACTOR``/``SCNB_SCALE_FACTOR``, per-atom
    ``CHARGE``): these are keyed by ATOM TYPE / DIHEDRAL TYPE, not by atom
    index, so their values are identical whether read from the whole-system
    table or from ParmEd's own per-prototype-sliced-and-pruned copy of the same
    table -- only the row/column numbering would differ, which this function
    never needs (it always indexes by the atom's own GLOBAL
    ``ATOM_TYPE_INDEX``, not a prototype-local one).

    The 1-4 pairs are read from the dihedral pointer quintuples
    (``DIHEDRALS_INC_HYDROGEN`` + ``DIHEDRALS_WITHOUT_HYDROGEN``, in THAT
    order -- matches the pre-rewrite parmed-slice-based behaviour byte for
    byte): the 3rd pointer < 0 marks a suppressed 1-4 interaction (skipped)
    and the 4th < 0 marks an improper (skipped). Only entries whose first atom
    belongs to *atom_indices* are kept (AMBER bonded-term sections only ever
    reference atoms within a single connected/bonded molecule, so checking the
    first atom is sufficient -- see ``docs/specs/fast-amber-loader.md`` §3(B)).
    ``SCEE``/``SCNB`` factors are absorbed into the returned charge product and
    epsilon respectively. Pairs already accounted for as 1-4 interactions are
    removed from the exclusion list so each pair appears in exactly one table.

    Parameters
    ----------
    raw_data : dict
        The ``raw_data`` dict from :func:`parse_prmtop` (WHOLE system, prmtop
        order). Must contain the LJ 1-4 (or plain LJ) coefficient tables,
        ``NONBONDED_PARM_INDEX``, ``SCEE_SCALE_FACTOR``, ``SCNB_SCALE_FACTOR``,
        the dihedral pointer lists, ``ATOM_TYPE_INDEX``, ``CHARGE``, and the
        excluded-atom lists.
    atom_indices : np.ndarray
        Ascending 0-based GLOBAL prmtop atom indices of the one molecule
        instance/prototype to build tables for (see
        ``amber_loader.MoleculePartition``).
    ene_conv : float, optional
        kcal mol⁻¹ -> target energy unit factor (default: kcal -> kJ).
    length_conv : float, optional
        Å -> target length unit factor (default: Å -> nm).

    Returns
    -------
    NonbondedTables
        ``scaling14`` and ``exclusions`` record lists, in LOCAL (0-based,
        position within *atom_indices*) atom indices.
    """
    natom_total = len(raw_data["ATOM_NAME"])
    local_of_global = np.full(natom_total, -1, dtype=np.int64)
    local_of_global[atom_indices] = np.arange(len(atom_indices), dtype=np.int64)

    num_types = int(raw_data["POINTERS"][1])
    charge = np.asarray(raw_data["CHARGE"], dtype=np.float64)
    nb_type_index = np.asarray(raw_data["ATOM_TYPE_INDEX"], dtype=np.int64)  # 1-based

    lj14_a = raw_data.get("LENNARD_JONES_14_ACOEF")
    if lj14_a is None:
        lj14_a = raw_data["LENNARD_JONES_ACOEF"]
    lj14_b = raw_data.get("LENNARD_JONES_14_BCOEF")
    if lj14_b is None:
        lj14_b = raw_data["LENNARD_JONES_BCOEF"]
    lj14_a = np.asarray(lj14_a, dtype=np.float64)
    lj14_b = np.asarray(lj14_b, dtype=np.float64)
    nb_index = np.asarray(raw_data["NONBONDED_PARM_INDEX"], dtype=np.int64)
    scee_factors = np.asarray(raw_data["SCEE_SCALE_FACTOR"], dtype=np.float64)
    scnb_factors = np.asarray(raw_data["SCNB_SCALE_FACTOR"], dtype=np.float64)

    tables = NonbondedTables()
    seen_14: set[tuple[int, int]] = set()

    # INC_HYDROGEN before WITHOUT_HYDROGEN -- matches the pre-rewrite
    # (parmed-slice-based) function's iteration order exactly. (Note this is
    # the REVERSE of the WITHOUT-then-INC order ParmEd/amber_loader use to
    # build the `.dihedrals` collection for periodic-torsion ENERGY terms --
    # a pre-existing inconsistency in this codebase, preserved here rather
    # than "fixed", since fixing it could silently change which duplicate 1-4
    # pair record (rare: two dihedral paths sharing the same (i, l) atoms with
    # different SCEE/SCNB) wins the `seen_14` dedup.
    dihedral_ptrs = list(raw_data.get("DIHEDRALS_INC_HYDROGEN", [])) + list(
        raw_data.get("DIHEDRALS_WITHOUT_HYDROGEN", [])
    )

    for ii in range(0, len(dihedral_ptrs), 5):
        i_raw, _j_raw, k_raw, l_raw, dtype_idx = dihedral_ptrs[ii : ii + 5]

        if k_raw < 0:  # 1-4 interaction suppressed
            continue
        if l_raw < 0:  # improper dihedral
            continue

        atom_i_g = int(i_raw) // 3
        atom_l_g = int(l_raw) // 3

        local_i = int(local_of_global[atom_i_g])
        if local_i < 0:
            continue  # not this instance
        local_l = int(local_of_global[atom_l_g])
        if local_l < 0:
            continue

        key = (min(local_i, local_l), max(local_i, local_l))
        if key in seen_14:
            continue

        nb_i = int(nb_type_index[atom_i_g]) - 1  # 1-based -> 0-based
        nb_l = int(nb_type_index[atom_l_g]) - 1
        pair_idx = int(nb_index[nb_i * num_types + nb_l]) - 1
        if pair_idx < 0:
            continue

        acoef = lj14_a[pair_idx]
        bcoef = lj14_b[pair_idx]

        if acoef != 0.0 and bcoef != 0.0:
            epsilon_kcal = (bcoef**2) / (4.0 * acoef)
            r_min_ang = (2.0 * acoef / bcoef) ** (1.0 / 6.0)
            epsilon_kj = epsilon_kcal * ene_conv
            sigma_nm = r_min_ang * length_conv * _SIGMA_SCALE
        else:
            epsilon_kj = 0.0
            sigma_nm = 1.0 * length_conv  # placeholder; epsilon is 0

        scee = float(scee_factors[int(dtype_idx) - 1])
        scnb = float(scnb_factors[int(dtype_idx) - 1])

        seen_14.add(key)
        tables.scaling14.append(
            Scaling14Record(
                i_local=key[0],
                l_local=key[1],
                charge_product=float(charge[atom_i_g] * charge[atom_l_g]) / scee,
                epsilon=epsilon_kj / scnb,
                sigma=sigma_nm,
            )
        )

    logger.debug("Parsed %d 1-4 scaling pairs from prmtop.", len(tables.scaling14))

    # ---- Explicit exclusions -------------------------------------------- #
    n_excluded_list = np.asarray(raw_data["NUMBER_EXCLUDED_ATOMS"], dtype=np.int64)
    excluded_atoms = np.asarray(raw_data["EXCLUDED_ATOMS_LIST"], dtype=np.int64)
    offsets = np.concatenate(([0], np.cumsum(n_excluded_list)))

    seen_excl: set[tuple[int, int]] = set(seen_14)
    for local_i, g in enumerate(atom_indices.tolist()):
        n = int(n_excluded_list[g])
        off = int(offsets[g])
        for j_1based in excluded_atoms[off : off + n].tolist():
            if j_1based <= 0:
                continue  # j=0 placeholder for atoms with no exclusions
            j_g = j_1based - 1  # prmtop is 1-based
            local_j = int(local_of_global[j_g])
            if local_j < 0:
                continue
            key = (min(local_i, local_j), max(local_i, local_j))
            if key in seen_excl:
                continue
            seen_excl.add(key)
            tables.exclusions.append(ExclusionRecord(i_local=key[0], j_local=key[1]))

    logger.debug("Parsed %d exclusions from prmtop.", len(tables.exclusions))
    return tables
