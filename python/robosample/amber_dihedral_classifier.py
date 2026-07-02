"""
DihedralClassifier: route any AMBER/CHAMBER dihedral (a real ParmEd
``Dihedral`` or the fast-loader's dihedral shim -- see ``amber_loader.py``)
to its DihedralType.

Architecture
------------
Dispatch is based on the residue names of the **two central-bond atoms**
(atom2 and atom3).  This avoids per-call iteration over all four atoms and
keeps the hot path to a handful of ``frozenset.__contains__`` tests plus one
``dict.get``.

Dispatch order (first applicable block wins):

1.  **Glycan-protein linkage** (OLS/OLT/NLN/HYP central atoms): try protein
    backbone lookup first (these modified AAs share backbone atom names with
    their parent AAs), then fall through to the linkage-specific lookup.
    This block is intentionally *before* the plain protein block to prevent
    the protein branch from silently absorbing cross-class psi dihedrals
    (one central atom in the linkage residue, the other in a GLYCAM sugar).

2.  **Protein**: both central atoms in PROTEIN_RESIDUES; neither in
    LIPID_RESIDUES or GLYCAM_RESIDUES.  Uses ``or`` across the central bond
    to handle cross-residue backbone omega (C-N spans two residues).

3.  **Lipid**: at least one central atom in LIPID_RESIDUES; neither in
    PROTEIN_RESIDUES.

4.  **Nucleic acid**: at least one central atom in NUCLEIC_RESIDUES.

5.  **Glycan** (pure sugar-to-sugar or intra-sugar): at least one central
    atom in GLYCAM_RESIDUES.

6.  **User-defined**: ``extra_lookup`` / ``extra_residue_map`` provided at
    construction.

7.  ``DihedralType.UNKNOWN``.

Why dispatch on central atoms only?
    The central bond (a2-a3) uniquely determines the molecule class for all
    standard AMBER topologies.  Outer atoms (a1, a4) can legitimately span
    residues (backbone cross-residue dihedrals), so checking all four atoms
    would introduce false exclusions.

Topology validation
    Backbone dihedrals (protein phi/psi/omega, nucleic alpha/epsilon/zeta)
    share atom-name tuples with plausible intra-residue patterns.  For these,
    residue-index topology is validated after the name lookup:

        phi / alpha   (i-1, i,   i,   i)
        psi / epsilon (i,   i,   i,   i+1)
        omega / zeta  (i,   i,   i+1, i+1)

Ions and coordinated metals
    Standard AMBER ions (Na+, K+, Ca2+, Mg2+, etc.) are point particles with
    no covalent bonds; the loader generates no dihedrals for them.  MCPB-style
    bonded metal complexes use custom residue names absent from all built-in
    residue sets and therefore return ``UNKNOWN``.  Use ``extra_residue_map``
    to route them if classification is needed.

Known limitations
-----------------
* Furanose tau-1 (C1-C2-C3-C4) returns ``GLYCAN_PYRANOSE_TAU_2`` for both
  pyranoses and furanoses; disambiguation requires ring-topology traversal.
* Sn-2 acyl chi-2 through chi-18 are indistinguishable by atom name from
  their sn-1 counterparts and are not classified.
* GLYCAM_RESIDUES is representative, not exhaustive; add names to
  ``GLYCAM_RESIDUES`` to cover exotic carbohydrates.
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple, Union

from .amber_dihedral_tables import (
    LIPID_DIHEDRAL_LOOKUP,
    LIPID_RESIDUES,
    NUCLEIC_DIHEDRAL_LOOKUP,
    NUCLEIC_RESIDUES,
    PROTEIN_DIHEDRAL_LOOKUP,
    PROTEIN_RESIDUES,
)
from .amber_dihedral_types import DihedralType

# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------

_Atom = Any  # a real ParmEd Atom or the fast-loader's atom shim (see amber_loader.py)
_Quad = Tuple[str, str, str, str]
_UserTable = Dict[_Quad, int]

# ---------------------------------------------------------------------------
# Atom-name helpers
# ---------------------------------------------------------------------------


def _names(a1: _Atom, a2: _Atom, a3: _Atom, a4: _Atom) -> _Quad:
    """Return the 4-tuple of atom names in bond order."""
    return (a1.name, a2.name, a3.name, a4.name)


def _ridx(a1: _Atom, a2: _Atom, a3: _Atom, a4: _Atom) -> Tuple[int, int, int, int]:
    """Return the residue indices of the four dihedral atoms."""
    return (
        a1.residue.idx,
        a2.residue.idx,
        a3.residue.idx,
        a4.residue.idx,
    )


def _carbon_idx(name: str) -> Optional[int]:
    """
    Parse the integer suffix of a carbon atom name.

    Returns ``None`` for names without a pure-integer suffix (e.g. ``CA``,
    ``CB``, ``CG``).  Examples: ``C5`` -> 5, ``C12`` -> 12, ``CA`` -> None.
    """
    if not name.startswith("C"):
        return None
    try:
        return int(name[1:])
    except ValueError:
        return None


def _oxygen_idx(name: str) -> Optional[int]:
    """
    Parse the integer suffix of an oxygen atom name.

    Returns ``None`` for names without a pure-integer suffix (e.g. ``OG``,
    ``OD1``).  Examples: ``O4`` -> 4, ``O`` -> None.
    """
    if not name.startswith("O"):
        return None
    try:
        return int(name[1:])
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# Topology predicates for cross-residue backbone dihedrals
# ---------------------------------------------------------------------------


def _is_phi_topology(r1: int, r2: int, r3: int, r4: int) -> bool:
    """
    Test for phi-like cross-residue pattern ``(i-1, i, i, i)``.

    Atom 1 is in the residue preceding the central-bond residue.
    """
    return (r2 == r3 == r4 and r1 == r2 - 1) or (r1 == r2 == r3 and r4 == r3 - 1)


def _is_psi_topology(r1: int, r2: int, r3: int, r4: int) -> bool:
    """
    Test for psi-like cross-residue pattern ``(i, i, i, i+1)``.

    Atom 4 is in the residue following the central-bond residue.
    """
    return (r1 == r2 == r3 and r4 == r3 + 1) or (r2 == r3 == r4 and r1 == r2 + 1)


def _is_omega_topology(r1: int, r2: int, r3: int, r4: int) -> bool:
    """
    Test for omega-like cross-residue pattern ``(i, i, i+1, i+1)``.

    Atoms 3 and 4 are in the residue following atoms 1 and 2.
    """
    return r1 == r2 and r3 == r4 and abs(r3 - r2) == 1


# ---------------------------------------------------------------------------
# Per-class sub-classifiers
# ---------------------------------------------------------------------------


def _classify_protein(a1: _Atom, a2: _Atom, a3: _Atom, a4: _Atom) -> DihedralType:
    """
    Classify a dihedral whose central atoms belong to protein residues.

    A single forward-direction lookup suffices because the tables are
    bidirectional.  For phi, psi, and omega the residue-index topology is
    validated to reject spurious intra-residue matches of the same
    atom-name pattern.

    This function is also called from the glycan-linkage branch to classify
    backbone dihedrals of modified amino acids (OLS/OLT/NLN/HYP), which share
    the same backbone atom names as their parent residues.
    """
    dtype = PROTEIN_DIHEDRAL_LOOKUP.get(_names(a1, a2, a3, a4), DihedralType.UNKNOWN)
    if dtype is DihedralType.UNKNOWN:
        return DihedralType.UNKNOWN

    if dtype in (
        DihedralType.PROTEIN_PHI,
        DihedralType.PROTEIN_PSI,
        DihedralType.PROTEIN_OMEGA,
    ):
        r1, r2, r3, r4 = _ridx(a1, a2, a3, a4)
        if dtype is DihedralType.PROTEIN_PHI and not _is_phi_topology(r1, r2, r3, r4):
            return DihedralType.UNKNOWN
        if dtype is DihedralType.PROTEIN_PSI and not _is_psi_topology(r1, r2, r3, r4):
            return DihedralType.UNKNOWN
        if dtype is DihedralType.PROTEIN_OMEGA and not _is_omega_topology(
            r1, r2, r3, r4
        ):
            return DihedralType.UNKNOWN

    return dtype


def _classify_lipid(a1: _Atom, a2: _Atom, a3: _Atom, a4: _Atom) -> DihedralType:
    """
    Classify a dihedral whose central atoms belong to lipid residues.

    Pure lookup; LIPID21 headgroup quads are globally unique.
    """
    return LIPID_DIHEDRAL_LOOKUP.get(_names(a1, a2, a3, a4), DihedralType.UNKNOWN)


def _classify_nucleic(a1: _Atom, a2: _Atom, a3: _Atom, a4: _Atom) -> DihedralType:
    """
    Classify a dihedral whose central atoms belong to nucleic-acid residues.

    Sugar ring torsions (nu0-nu4) and glycosidic chi are pure lookups.
    Backbone alpha (phi topology), epsilon (psi topology), and zeta (omega
    topology) additionally require residue-index validation.
    """
    dtype = NUCLEIC_DIHEDRAL_LOOKUP.get(_names(a1, a2, a3, a4), DihedralType.UNKNOWN)
    if dtype is DihedralType.UNKNOWN:
        return DihedralType.UNKNOWN

    if dtype in (
        DihedralType.NUCLEIC_ALPHA,
        DihedralType.NUCLEIC_EPSILON,
        DihedralType.NUCLEIC_ZETA,
    ):
        r1, r2, r3, r4 = _ridx(a1, a2, a3, a4)
        if dtype is DihedralType.NUCLEIC_ALPHA and not _is_phi_topology(r1, r2, r3, r4):
            return DihedralType.UNKNOWN
        if dtype is DihedralType.NUCLEIC_EPSILON and not _is_psi_topology(
            r1, r2, r3, r4
        ):
            return DihedralType.UNKNOWN
        if dtype is DihedralType.NUCLEIC_ZETA and not _is_omega_topology(
            r1, r2, r3, r4
        ):
            return DihedralType.UNKNOWN

    return dtype


# ---------------------------------------------------------------------------
# Main classifier
# ---------------------------------------------------------------------------


class AmberDihedralClassifier:
    """
    Stateless (per instance) classifier mapping an AMBER/CHAMBER dihedral
    (a real ParmEd ``Dihedral`` or the fast-loader's dihedral shim) to a
    ``DihedralType`` (or a user-defined integer).

    Parameters
    ----------
    extra_lookup:
        Additional atom-name quads mapped to integer type values.  Quads are
        bidirectionalised at construction (both forward and reverse stored).
        Values should be >= ``DihedralType.USER_DEFINED`` (100).
    extra_residue_map:
        Maps residue names to molecule-class strings so the corresponding
        built-in sub-classifier handles their dihedrals.  Valid class strings:
        ``"protein"``, ``"lipid"``, ``"nucleic"``, ``"glycan"``,
        ``"linkage"``.

    Examples
    --------
    Standard usage::

        clf = DihedralClassifier()
        dtype = clf.classify(dihedral)

    Custom PTM (methylated backbone)::

        MY_METHYL_CHI = 100
        clf = DihedralClassifier(
            extra_lookup={("N", "CA", "CB", "CM"): MY_METHYL_CHI},
            extra_residue_map={"ACEM": "protein"},
        )

    Range checks::

        is_protein = 1 <= dtype <= 19
        is_custom  = dtype >= DihedralType.USER_DEFINED
    """

    __slots__ = ("_extra_lookup", "_extra_residue_map")

    def __init__(
        self,
        extra_lookup: Optional[Dict[_Quad, int]] = None,
        extra_residue_map: Optional[Dict[str, str]] = None,
    ) -> None:
        # Bidirectionalise the extra lookup at construction time.
        bidi: _UserTable = {}
        for quad, val in (extra_lookup or {}).items():
            rev: _Quad = (quad[3], quad[2], quad[1], quad[0])
            bidi[quad] = val
            bidi[rev] = val
        self._extra_lookup: _UserTable = bidi
        self._extra_residue_map: Dict[str, str] = dict(extra_residue_map or {})

    # ------------------------------------------------------------------

    def classify(
        self, dihedral: Any
    ) -> Union[DihedralType, int]:
        """
        Classify an AMBER/CHAMBER dihedral and return its biochemical role.

        Parameters
        ----------
        dihedral:
            A real ParmEd ``Dihedral`` or the fast-loader's dihedral shim
            (``amber_loader._DihedralShim`` / the throwaway
            ``acyclic_graph._DihedralCandidate``), with non-``None``
            ``.atom1`` through ``.atom4`` attributes and valid
            ``.residue.idx`` values.

        Returns
        -------
        DihedralType or int
            A ``DihedralType`` member for all built-in types.
            A plain ``int >= DihedralType.USER_DEFINED`` for user-defined
            types.  ``DihedralType.UNKNOWN`` (0) when no match is found.

        Notes
        -----
        * Protein backbone phi/psi/omega and nucleic alpha/epsilon/zeta
          require both a name match AND residue-index topology validation.
        * Glycan tau-1 ``(C1, C2, C3, C4)`` returns
          ``GLYCAN_PYRANOSE_TAU_2`` for both pyranoses and furanoses.
        * Ions and MCPB-style metal complexes: residue names absent from all
          built-in sets pass through to UNKNOWN (or to the user lookup if
          registered via ``extra_residue_map``).
        """
        a1: _Atom = dihedral.atom1
        a2: _Atom = dihedral.atom2
        a3: _Atom = dihedral.atom3
        a4: _Atom = dihedral.atom4

        rn2: str = a2.residue.name
        rn3: str = a3.residue.name

        # --- 1. Glycan-protein linkage ------------------------------------------------

        # --- 2. Protein ------------------------------------------------------
        # Use 'or' to handle cross-residue omega (C-N bond spans two residues).
        # Exclude GLYCAM residues as outer atoms to prevent false routing when
        # a custom glycoconjugate creates an unexpected bond pattern.
        if rn2 in PROTEIN_RESIDUES or rn3 in PROTEIN_RESIDUES:
            if rn2 not in LIPID_RESIDUES and rn3 not in LIPID_RESIDUES:
                return _classify_protein(a1, a2, a3, a4)

        # --- 3. Lipid --------------------------------------------------------
        if rn2 in LIPID_RESIDUES or rn3 in LIPID_RESIDUES:
            if rn2 not in PROTEIN_RESIDUES and rn3 not in PROTEIN_RESIDUES:
                return _classify_lipid(a1, a2, a3, a4)

        # --- 4. Nucleic acid -------------------------------------------------
        if rn2 in NUCLEIC_RESIDUES or rn3 in NUCLEIC_RESIDUES:
            return _classify_nucleic(a1, a2, a3, a4)

        # --- 6. User-defined -------------------------------------------------
        if self._extra_lookup or self._extra_residue_map:
            return self._classify_custom(a1, a2, a3, a4, rn2, rn3)

        return DihedralType.UNKNOWN

    # ------------------------------------------------------------------

    def _classify_custom(
        self,
        a1: _Atom,
        a2: _Atom,
        a3: _Atom,
        a4: _Atom,
        rn2: str,
        rn3: str,
    ) -> Union[DihedralType, int]:
        """
        Try user-provided lookup and residue-class routing.

        Atom-name quads in ``extra_lookup`` take priority over
        ``extra_residue_map`` routing.

        Parameters
        ----------
        a1, a2, a3, a4:
            Dihedral atoms in bond order.
        rn2, rn3:
            Residue names of atoms a2 and a3 (pre-computed by caller).

        Returns
        -------
        DihedralType or int
        """
        # Direct quad lookup (bidirectionalised at construction).
        val = self._extra_lookup.get(_names(a1, a2, a3, a4))
        if val is not None:
            return val

        # Residue-class routing: delegate to the appropriate sub-classifier.
        cls2 = self._extra_residue_map.get(rn2)
        cls3 = self._extra_residue_map.get(rn3)
        cls = cls2 or cls3
        if cls == "protein":
            return _classify_protein(a1, a2, a3, a4)
        if cls == "lipid":
            return _classify_lipid(a1, a2, a3, a4)
        if cls == "nucleic":
            return _classify_nucleic(a1, a2, a3, a4)

        return DihedralType.UNKNOWN
