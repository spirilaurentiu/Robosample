"""
DihedralType: integer enumeration of biochemical dihedral angle roles.
"""

from __future__ import annotations

from enum import IntEnum, auto


class DihedralType(IntEnum):
    """
    Integer enumeration identifying the biochemical role of a dihedral angle.

    Using ``IntEnum`` (rather than plain ``Enum``) allows direct integer
    comparison, numpy array storage, and fast range-based membership tests::

    User-defined types
    ------------------
    Extensions should use integer values >= ``DihedralType.USER_DEFINED`` (100).
    Assign them as plain ``int`` constants and pass them via the classifier's
    ``extra_lookup`` / ``extra_residue_map`` constructor arguments::

        MY_PTM_CHI = 100
        clf = DihedralClassifier(
            extra_lookup={("N", "CA", "CB", "CM"): MY_PTM_CHI},
            extra_residue_map={"MYMOD": "protein"},
        )

    ``classify()`` returns a ``DihedralType`` member for built-in values and a
    plain ``int`` for user-defined values.  Since ``DihedralType`` extends
    ``int``, all comparison and range operations work uniformly across both.
    """

    UNKNOWN = auto()

    # ------------------------------------------------------------------ protein

    PROTEIN_PHI = auto()
    """phi: C(i-1) - N(i) - CA(i) - C(i)"""

    PROTEIN_PSI = auto()
    """psi: N(i) - CA(i) - C(i) - N(i+1)"""

    PROTEIN_OMEGA = auto()
    """omega: CA(i) - C(i) - N(i+1) - CA(i+1)  (peptide-bond plane)"""

    PROTEIN_CHI_1 = auto()
    PROTEIN_CHI_2 = auto()
    PROTEIN_CHI_3 = auto()
    PROTEIN_CHI_4 = auto()
    PROTEIN_CHI_5 = auto()  # ARG only

    PROTEIN_RING_DIHEDRAL = auto()
    """
    Torsion angle whose four atoms all lie within a covalent ring.

    This label conveys that the dihedral is constrained by ring geometry rather
    than freely rotatable.  Covers:

    - PRO pyrrolidine ring:    CG  - CD  - N   - CA
    - HIS imidazole ring:      CG  - ND1 - CE1 - NE2
    - PHE / TYR benzene ring:  CG  - CD1 - CE1 - CZ
    - TRP pyrrole ring:        CG  - CD1 - NE1 - CE2
    - TRP benzene ring:        CE2 - CZ2 - CH2 - CZ3

    The CYX S-S bond is *not* included here.  It is a spanning-tree
    ring-closing cut handled by ``ring_closing.select_ring_closing_bond``.
    """

    # ------------------------------------------------------------------- lipid
    LIPID_ALPHA1 = auto()
    LIPID_ALPHA2 = auto()
    LIPID_ALPHA3 = auto()
    LIPID_ALPHA4 = auto()
    LIPID_ALPHA5 = auto()
    LIPID_ALPHA6 = auto()

    LIPID_THETA1 = auto()
    # LIPID_THETA2       = auto()
    # LIPID_THETA3       = auto()
    LIPID_THETA4 = auto()

    LIPID_SN1 = auto()
    LIPID_SN2 = auto()
    LIPID_AMIDE = auto()

    LIPID_GAMMA = auto()
    LIPID_BETA = auto()

    LIPID_D_ALKANE = auto()
    LIPID_D_ALKENE = auto()

    LIPID_CHL_OMEGA1 = auto()
    LIPID_CHL_OMEGA2 = auto()
    LIPID_CHL_OMEGA3 = auto()
    LIPID_CHL_OMEGA4 = auto()

    # ----------------------------------------------------------- nucleic acid
    NUCLEIC_ALPHA = auto()  # O3'(i-1) - P    - O5' - C5'
    NUCLEIC_BETA = auto()  # P         - O5' - C5' - C4'
    NUCLEIC_GAMMA = auto()  # O5'       - C5' - C4' - C3'
    NUCLEIC_DELTA = auto()  # C5'       - C4' - C3' - O3'
    NUCLEIC_EPSILON = auto()  # C4'       - C3' - O3' - P(i+1)
    NUCLEIC_ZETA = auto()  # C3'       - O3' - P(i+1) - O5'(i+1)

    NUCLEIC_CHI_PURINE = auto()  # O4' - C1' - N9  - C4
    NUCLEIC_CHI_PYRIMIDINE = auto()  # O4' - C1' - N1  - C2

    NUCLEIC_NU_0 = auto()  # C4' - O4' - C1' - C2'
    NUCLEIC_NU_1 = auto()  # O4' - C1' - C2' - C3'
    NUCLEIC_NU_2 = auto()  # C1' - C2' - C3' - C4'
    NUCLEIC_NU_3 = auto()  # C2' - C3' - C4' - O4'
    NUCLEIC_NU_4 = auto()  # C3' - C4' - O4' - C1'
    # -------------------------------------------------------------- glycan

    # --------------------------------------------------------- user-defined

    USER_DEFINED = auto()
    """
    First value in the user-defined range.

    Built-in members occupy 0 - 99.  Caller-registered custom types use
    integers >= 100.  Values returned by ``DihedralClassifier.classify()``
    in this range are plain ``int`` objects, not ``DihedralType`` members;
    test with ``dtype >= DihedralType.USER_DEFINED``.
    """

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return f"DihedralType.{self.name}"
