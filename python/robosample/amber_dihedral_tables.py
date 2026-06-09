"""
Lookup tables mapping 4-atom-name tuples to DihedralType values.

Each public ``*_LOOKUP`` dict is bidirectional: both the canonical forward
tuple and its reverse are stored at build time, so the classifier needs only
one ``dict.get()`` per dihedral regardless of atom ordering in the ParmEd
topology.

Build-time conflict detection raises ``ValueError`` immediately on import if
any two entries map the same 4-tuple to different types.

Residue-name frozensets (``PROTEIN_RESIDUES``, ``LIPID_RESIDUES``, etc.) are
used by the dispatcher to route each dihedral to the correct sub-classifier
without inspecting atom names.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

from amber_dihedral_types import DihedralType

# ---------------------------------------------------------------------------
# Internal types
# ---------------------------------------------------------------------------

_Quad = Tuple[str, str, str, str]
_Table = Dict[_Quad, DihedralType]
_Entries = List[Tuple[_Quad, DihedralType]]

# ---------------------------------------------------------------------------
# Build helpers
# ---------------------------------------------------------------------------


def _build_lookup(entries: _Entries) -> _Table:
    """
    Construct a bidirectional lookup table from a list of ``(quad, dtype)``
    pairs.

    Both the forward quad and its reverse are inserted.  A ``ValueError`` is
    raised at import time if any two entries (after expansion) map the same
    quad to *different* ``DihedralType`` values.

    Parameters
    ----------
    entries:
        ``((a1, a2, a3, a4), DihedralType)`` pairs in canonical forward order.

    Returns
    -------
    dict
        Flat ``{(a1, a2, a3, a4): DihedralType}`` mapping; both forward and
        reverse keys populated.

    Raises
    ------
    ValueError
        If the same 4-atom quad resolves to two different types.
    """
    table: _Table = {}
    for quad, dtype in entries:
        rev: _Quad = (quad[3], quad[2], quad[1], quad[0])
        for key in (quad, rev):
            existing = table.get(key)
            if existing is not None and existing is not dtype:
                raise ValueError(
                    f"Conflict: quad {key} maps to both "
                    f"{existing.name} and {dtype.name}"
                )
            table[key] = dtype
    return table


# ---------------------------------------------------------------------------
# Protein
# ---------------------------------------------------------------------------

#: Bidirectional lookup for all ff19SB protein dihedrals (backbone + side
#: chains).  Keys are 4-atom-name tuples in either order.
#:
#: Notable absences:
#: * CYX SG-SG: spanning-tree cut bond, not a torsional DOF; see ring_closing.
#: * MET methyl rotation (CG-SD-CE-HE): not a standard named angle.
#: * TYR O-H rotation (CZ-OH): not tracked in standard ff19SB analyses.
PROTEIN_DIHEDRAL_LOOKUP: _Table = _build_lookup([
    # ---- backbone -------------------------------------------------------
    (("C",   "N",   "CA",  "C"),   DihedralType.PROTEIN_PHI),
    (("N",   "CA",  "C",   "N"),   DihedralType.PROTEIN_PSI),
    # omega: canonical (CA-C-N-CA), carbonyl-O variant, ACE cap, NME cap
    (("CA",  "C",   "N",   "CA"),  DihedralType.PROTEIN_OMEGA),
    (("O",   "C",   "N",   "CA"),  DihedralType.PROTEIN_OMEGA),
    (("CH3", "C",   "N",   "CA"),  DihedralType.PROTEIN_OMEGA),
    (("CA",  "C",   "N",   "CH3"), DihedralType.PROTEIN_OMEGA),
    (("O",   "C",   "N",   "CH3"), DihedralType.PROTEIN_OMEGA),
    # ---- chi-1: shared across ALA-like side chains ----------------------
    # ARG, ASN, ASP, ASH, GLN, GLU, GLH, LEU, MET, PHE, PRO, TRP, TYR
    (("N",   "CA",  "CB",  "CG"),  DihedralType.PROTEIN_CHI_1),
    # CYS, CYX, CYM
    (("N",   "CA",  "CB",  "SG"),  DihedralType.PROTEIN_CHI_1),
    # ILE
    (("N",   "CA",  "CB",  "CG1"), DihedralType.PROTEIN_CHI_1),
    # SER
    (("N",   "CA",  "CB",  "OG"),  DihedralType.PROTEIN_CHI_1),
    # THR
    (("N",   "CA",  "CB",  "OG1"), DihedralType.PROTEIN_CHI_1),
    # VAL (CG1 already covered by ILE entry above)
    # ---- chi-2 ----------------------------------------------------------
    # ARG, LYS/LYN, GLN, GLU/GLH, PRO (CB-CG-CD)
    (("CA",  "CB",  "CG",  "CD"),  DihedralType.PROTEIN_CHI_2),
    # ASN, ASP, ASH: beta-carboxamide/carboxylate
    (("CA",  "CB",  "CG",  "OD1"), DihedralType.PROTEIN_CHI_2),
    # HIS (HID/HIE/HIP)
    (("CA",  "CB",  "CG",  "ND1"), DihedralType.PROTEIN_CHI_2),
    # ILE
    (("CA",  "CB",  "CG1", "CD1"), DihedralType.PROTEIN_CHI_2),
    # LEU, PHE, TRP, TYR (CG-CD1; LEU and PHE share the same tuple)
    (("CA",  "CB",  "CG",  "CD1"), DihedralType.PROTEIN_CHI_2),
    # MET
    (("CA",  "CB",  "CG",  "SD"),  DihedralType.PROTEIN_CHI_2),
    # ---- chi-3 ----------------------------------------------------------
    # ARG (CD-NE-CZ)
    (("CB",  "CG",  "CD",  "NE"),  DihedralType.PROTEIN_CHI_3),
    # GLN / GLU / GLH (CD-OE1)
    (("CB",  "CG",  "CD",  "OE1"), DihedralType.PROTEIN_CHI_3),
    # LYS / LYN (CD-CE)
    (("CB",  "CG",  "CD",  "CE"),  DihedralType.PROTEIN_CHI_3),
    # MET (SD-CE)
    (("CB",  "CG",  "SD",  "CE"),  DihedralType.PROTEIN_CHI_3),
    # ---- chi-4 ----------------------------------------------------------
    # ARG (NE-CZ)
    (("CG",  "CD",  "NE",  "CZ"),  DihedralType.PROTEIN_CHI_4),
    # LYS / LYN (CE-NZ)
    (("CG",  "CD",  "CE",  "NZ"),  DihedralType.PROTEIN_CHI_4),
    # ---- chi-5 (ARG only) -----------------------------------------------
    (("CD",  "NE",  "CZ",  "NH1"), DihedralType.PROTEIN_CHI_5),
    # ---- ring dihedrals -------------------------------------------------
    # PRO pyrrolidine: CD-N bond closes the ring; CG-CD-N-CA is the ring
    # dihedral that describes ring pucker.  Do NOT cut N-CA (backbone phi)
    # or CA-CB / CB-CG (chi-1/2).  Spanning-tree cut = N-CD bond.
    (("CG",  "CD",  "N",   "CA"),  DihedralType.PROTEIN_RING_DIHEDRAL),
    # HIS imidazole (aromatic 5-ring): CG-ND1-CE1-NE2
    (("CG",  "ND1", "CE1", "NE2"), DihedralType.PROTEIN_RING_DIHEDRAL),
    # PHE / TYR benzene (aromatic 6-ring): CG-CD1-CE1-CZ
    (("CG",  "CD1", "CE1", "CZ"),  DihedralType.PROTEIN_RING_DIHEDRAL),
    # TRP pyrrole ring (5-ring, fused): CG-CD1-NE1-CE2
    (("CG",  "CD1", "NE1", "CE2"), DihedralType.PROTEIN_RING_DIHEDRAL),
    # TRP benzene ring (6-ring, fused): CE2-CZ2-CH2-CZ3
    (("CE2", "CZ2", "CH2", "CZ3"), DihedralType.PROTEIN_RING_DIHEDRAL),
])

#: All ff19SB residue names recognised as protein by the dispatcher.
PROTEIN_RESIDUES: frozenset[str] = frozenset([
    "ALA", "ARG", "ASN", "ASP", "ASH", "AS4",
    "CYS", "CYX", "CYM",
    "GLN", "GLU", "GLH", "GL4", "GLY",
    "HID", "HIE", "HIP",
    "ILE", "LEU", "LYS", "LYN",
    "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "ACE", "NME",
])

# ---------------------------------------------------------------------------
# Lipid
# ---------------------------------------------------------------------------

LIPID_DIHEDRAL_LOOKUP: _Table = _build_lookup([
    # Alpha
    (("C2", "C3", "O31", "P31"),   DihedralType.LIPID_ALPHA1),
    (("C3", "O31", "P31", "O32"),  DihedralType.LIPID_ALPHA2),
    (("O31", "P31", "O32", "C31"), DihedralType.LIPID_ALPHA3),
    (("P31", "O32", "C31", "C32"), DihedralType.LIPID_ALPHA4),
    (("O32", "C31", "C32", "C33"), DihedralType.LIPID_ALPHA5),
    (("O32", "C31", "C32", "N31"), DihedralType.LIPID_ALPHA5),
    (("C31", "C32", "C33", "O36"), DihedralType.LIPID_ALPHA6),

    # Theta
    (("O31", "C3", "C2", "N11"),   DihedralType.LIPID_THETA1),
    (("O11", "C1", "C2", "O21"),   DihedralType.LIPID_THETA1),

    # (("O31", "C3", "C2", "C1"),    DihedralType.LIPID_THETA2),
    # (("O11", "C1", "C2", "C3"),    DihedralType.LIPID_THETA2),

    # (("C3", "C2", "C1", "O11"),    DihedralType.LIPID_THETA3),
    # (("C1", "C2", "C3", "O31"),    DihedralType.LIPID_THETA3),

    (("N11", "C2", "C1", "O11"),   DihedralType.LIPID_THETA4),
    (("O21", "C2", "C3", "O31"),   DihedralType.LIPID_THETA4),

    # Planarity / rigid
    (("C1", "O11", "C11", "O12"),  DihedralType.LIPID_SN1),
    (("C2", "O21", "C21", "O22"),  DihedralType.LIPID_SN2),
    (("C2", "N11", "C11", "O12"),  DihedralType.LIPID_AMIDE),

    # Beta / Gamma
    (("C3", "C2", "N11", "C11"),   DihedralType.LIPID_BETA),
    (("C2", "C1", "O11", "C11"),   DihedralType.LIPID_BETA),

    (("C3", "C2", "C1", "C12"),    DihedralType.LIPID_GAMMA),
    (("C3", "C2", "O21", "C21"),   DihedralType.LIPID_GAMMA),

    # Alkane chain torsions
    (("O11", "C11", "C12", "C13"),     DihedralType.LIPID_D_ALKANE),
    (("O21", "C21", "C12", "C13"),     DihedralType.LIPID_D_ALKANE),
    (("N11", "C11", "C12", "C13"),     DihedralType.LIPID_D_ALKANE),

    (("C11", "C12", "C13", "C14"),     DihedralType.LIPID_D_ALKANE),
    (("C21", "C12", "C13", "C14"),     DihedralType.LIPID_D_ALKANE),
    (("C2", "C1", "C12", "C13"),       DihedralType.LIPID_D_ALKANE),

    (("C12", "C13", "C14", "C15"),     DihedralType.LIPID_D_ALKANE),
    (("C13", "C14", "C15", "C16"),     DihedralType.LIPID_D_ALKANE),
    (("C14", "C15", "C16", "C17"),     DihedralType.LIPID_D_ALKANE),
    (("C15", "C16", "C17", "C18"),     DihedralType.LIPID_D_ALKANE),
    (("C16", "C17", "C18", "C19"),     DihedralType.LIPID_D_ALKANE),
    (("C17", "C18", "C19", "C110"),    DihedralType.LIPID_D_ALKANE),
    (("C18", "C19", "C110", "C111"),   DihedralType.LIPID_D_ALKANE),
    (("C19", "C110", "C111", "C112"),  DihedralType.LIPID_D_ALKANE),
    (("C110", "C111", "C112", "C113"), DihedralType.LIPID_D_ALKANE),
    (("C111", "C112", "C113", "C114"), DihedralType.LIPID_D_ALKANE),
    (("C112", "C113", "C114", "C115"), DihedralType.LIPID_D_ALKANE),
    (("C113", "C114", "C115", "C116"), DihedralType.LIPID_D_ALKANE),
    (("C114", "C115", "C116", "C117"), DihedralType.LIPID_D_ALKANE),
    (("C115", "C116", "C117", "C118"), DihedralType.LIPID_D_ALKANE),
    (("C116", "C117", "C118", "C119"), DihedralType.LIPID_D_ALKANE),
    (("C117", "C118", "C119", "C120"), DihedralType.LIPID_D_ALKANE),
    (("C118", "C119", "C120", "C121"), DihedralType.LIPID_D_ALKANE),
    (("C119", "C120", "C121", "C122"), DihedralType.LIPID_D_ALKANE),

    # Omega
    (("C13", "C17", "C20", "C22"), DihedralType.LIPID_CHL_OMEGA1),
    (("C17", "C20", "C22", "C23"), DihedralType.LIPID_CHL_OMEGA2),
    (("C20", "C22", "C23", "C24"), DihedralType.LIPID_CHL_OMEGA3),
    (("C22", "C23", "C24", "C25"), DihedralType.LIPID_CHL_OMEGA4),
])

#: AMBER LIPID21 residue names recognised by the dispatcher.
LIPID_RESIDUES: frozenset[str] = frozenset([
    "LAL", # C12:0 (lauroyl)
    "MY", # C14:0 (myristoyl)
    "PA", # C16:0 (palmitoyl)
    "ST", # C18:0 (stearoyl)
    "OL", # C18:1 (oleoyl)
    "AR", # C20:4 (arachidonoyl)
    "DHA", # C22:6 (docosahexaenoyl)
    "SA", # Sphingosine (SPM backbone)
    "PC", # Phosphatidylcholine
    "PE", # Phosphatidylethanolamine
    "PS", # Phosphatidylserine
    "PG", # Phosphatidylglycerol
    "PGR", # Phosphatidylglycerol stereoisomer (same torsions as PG)
    "PA", # Phosphatidic acid
    "PH", # Phosphatidyl phosphate
    "PH-", # Phosphatidyl hydrogen phosphate (PA with one less proton)
    "SPM", # Sphingomyelin
    "CHL", # Cholesterol
])

# ---------------------------------------------------------------------------
# Nucleic acid
# ---------------------------------------------------------------------------

#: Bidirectional lookup for AMBER OL3/OL15 nucleic-acid dihedrals.
NUCLEIC_DIHEDRAL_LOOKUP: _Table = _build_lookup([
    # ---- backbone -------------------------------------------------------
    # alpha:   O3'(i-1) - P(i) - O5'(i) - C5'(i)   [phi topology]
    (("O3'",  "P",    "O5'",  "C5'"),  DihedralType.NUCLEIC_ALPHA),
    # beta, gamma, delta: intra-residue
    (("P",    "O5'",  "C5'",  "C4'"),  DihedralType.NUCLEIC_BETA),
    (("O5'",  "C5'",  "C4'",  "C3'"),  DihedralType.NUCLEIC_GAMMA),
    (("C5'",  "C4'",  "C3'",  "O3'"),  DihedralType.NUCLEIC_DELTA),
    # epsilon: C4'(i) - C3'(i) - O3'(i) - P(i+1)   [psi topology]
    (("C4'",  "C3'",  "O3'",  "P"),    DihedralType.NUCLEIC_EPSILON),
    # zeta:    C3'(i) - O3'(i) - P(i+1) - O5'(i+1) [omega topology]
    (("C3'",  "O3'",  "P",    "O5'"),  DihedralType.NUCLEIC_ZETA),
    # ---- glycosidic chi -------------------------------------------------
    (("O4'",  "C1'",  "N9",   "C4"),   DihedralType.NUCLEIC_CHI_PURINE),
    (("O4'",  "C1'",  "N1",   "C2"),   DihedralType.NUCLEIC_CHI_PYRIMIDINE),
    # ---- sugar ring (nu angles) -----------------------------------------
    (("C4'",  "O4'",  "C1'",  "C2'"),  DihedralType.NUCLEIC_NU_0),
    (("O4'",  "C1'",  "C2'",  "C3'"),  DihedralType.NUCLEIC_NU_1),
    (("C1'",  "C2'",  "C3'",  "C4'"),  DihedralType.NUCLEIC_NU_2),
    (("C2'",  "C3'",  "C4'",  "O4'"),  DihedralType.NUCLEIC_NU_3),
    (("C3'",  "C4'",  "O4'",  "C1'"),  DihedralType.NUCLEIC_NU_4),
])

#: AMBER OL3 RNA and OL15 DNA residue names, including 5'- and 3'-terminal
#: variants.
NUCLEIC_RESIDUES: frozenset[str] = frozenset([
    # DNA
    "DA",  "DC",  "DG",  "DT",
    "DA3", "DA5", "DC3", "DC5", "DG3", "DG5", "DT3", "DT5",
    # RNA
    "A",   "C",   "G",   "U",
    "A3",  "A5",  "C3",  "C5",  "G3",  "G5",  "U3",  "U5",
])
