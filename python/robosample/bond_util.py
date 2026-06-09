"""
bond_classification.py
=======================
Heuristic classification of bonds in AMBER-family topology structures as
rigid (non-rotatable) or rotatable, based on atom-type strings.

Supported force fields
----------------------
  GAFF       Wang et al.,      J. Comput. Chem. 25, 1157-1174 (2004)
  GAFF2      He et al.,        J. Chem. Theory Comput. 16, 528-552 (2020)
  ff14SB     Maier et al.,     J. Chem. Theory Comput. 11, 3696-3713 (2015)
  ff19SB     Tian et al.,      J. Chem. Theory Comput. 16, 528-552 (2020)
  OL3        Zgarbova et al.,  J. Chem. Theory Comput. 7, 2886-2902 (2011)
  OL15       Zgarbova et al.,  J. Chem. Theory Comput. 11, 5723-5736 (2015)
  bsc1       Ivani et al.,     Nat. Methods 13, 55-58 (2016)
  LIPID17    Gould et al.,     J. Phys. Chem. B 122, 5573-5592 (2018)
  GLYCAM06j  Kirschner et al., J. Comput. Chem. 29, 622-655 (2008)

Atom-type string case conventions
----------------------------------
Atom types are compared case-sensitively (no .lower() normalisation).
The case conventions in standard AMBER prmtop files are:

  GAFF / GAFF2 / LIPID17   all lowercase   e.g. ca, ce, na, c2, os
  ff14SB / ff19SB           all uppercase   e.g. CA, CE, NA, C, N
  OL3 / OL15 / bsc1         uppercase, shared with protein types
  GLYCAM06j-1               mixed case      e.g. Cg, Os, Cp, Cy

Case-sensitive matching is essential to avoid type collisions that would
arise after case folding.  Known collisions that motivated this design:

  Lowercase "cp"  (GAFF)      -- bridgehead aromatic carbon (sp2)
  Mixed-case "Cp" (GLYCAM)    -- sp3 pyranose ring carbon (non-aromatic)

  Lowercase "os"  (GAFF)      -- ether / ester bridging oxygen (sp3)
  Mixed-case "Os" (GLYCAM)    -- glycosidic / ring oxygen (sp3)

Bond rigidity criteria
----------------------
A bond is classified as rigid if any of the following conditions holds,
tested in the order listed:

  0. Fast path: parmed.Bond.order >= 2 when an explicit bond order is
     available in the topology.  Bypasses all heuristics.

  1. Aromatic -- both atom types appear in AROMATIC_TYPES.

  2. Amide    -- the atom-type pair matches AMIDE_C x AMIDE_N,
                 capturing the ~15-20 kcal/mol C(=O)-N rotational barrier.

  3. Multiple -- both atom types appear in MULTIPLE_BOND_TYPES, indicating
                 sp or sp2 hybridisation consistent with double- or
                 triple-bond character.

The aromatic predicate does not perform ring-membership validation
(see Known Limitations).

Known Limitations
-----------------
  1. Biaryl inter-ring single bonds (e.g. biphenyl C-C, biindolyl C-C)
     connect two atoms that both carry aromatic type codes.  Without SSSR
     ring-membership validation, these bonds are over-classified as rigid.
     In practice this is rarely consequential for biological systems but
     affects drug-like small molecules with biaryl scaffolds.  Callers that
     require strict intra-ring validation should wrap is_rigid_bond() with
     an external ring check (e.g. via RDKit passed through parmed.rdkit).

  2. Phosphodiester P-O bonds are NOT classified as rigid.
     Non-bridging P-O bonds (atom type O2 in AMBER nucleic-acid FFs) have
     resonance character, but phosphorus geometry is approximately
     tetrahedral; the planarity constraint underlying sp2 rigidity does not
     apply.  The backbone P-OS bonds (bridging oxygens, alpha/zeta torsions)
     are correctly classified as rotatable because OS does not appear in any
     classification set.

  3. AMBER prmtop files do not encode bond order.  Criteria 1-3 therefore
     rely solely on atom-type strings.  Non-standard, custom, or
     parameterisation-file-specific types will not be recognised.

  4. GLYCAM06j-1 sp3 ring types (Cg, Cp, Cy) and oxygen types (Os, O3,
     O4, O5, O6) are intentionally absent from all classification sets.
     The glycosidic C1'-O-Cx linkage (involving Os) is correctly classified
     as rotatable.
"""

from __future__ import annotations

import parmed as pmd

# ==========================================================================
# Aromatic atom types
# ==========================================================================

# GAFF and GAFF2 aromatic types -- all lowercase.
# ce/cf are listed here because they appear in formally aromatic 5- and
# 6-membered ring systems in GAFF2, even though they also describe
# open-chain conjugated sp2 carbons.  The is_aromatic() predicate is
# always subject to the biaryl caveat documented above.
_AROMATIC_GAFF: frozenset[str] = frozenset(
    {
        # -- sp2 aromatic carbons ----------------------------------------------
        "ca",  # benzene-like aromatic C (6-membered rings)
        "cp",  # bridgehead aromatic C: ring-junction in fused systems (azulene)
        "cq",  # bridgehead aromatic C: 6-membered ring fused to 5-membered
        "cc",  # conjugated sp2 C in non-pure-aromatic 5-membered heterocycle
        "cd",  # conjugated sp2 C, alternating partner of cc
        "ce",  # conjugated sp2 C in open-chain or ring system (GAFF2 extended)
        "cf",  # conjugated sp2 C, alternating partner of ce
        # -- sp2 aromatic nitrogens --------------------------------------------
        "na",  # pyrrole-type N: lone pair in ring pi system, N-H bearing
        "nb",  # pyridine-type N: lone pair in sp2 plane (non-H-bearing)
        "nc",  # aromatic sp2 N in 5-membered ring
        "nd",  # aromatic sp2 N, alternating partner of nc
        "ne",  # GAFF2: aromatic N in 5-membered ring, additional variant
        "nf",  # GAFF2: aromatic N, alternating partner of ne
        # -- GAFF2: aromatic phosphorus and sulfur heteroatoms -----------------
        "pe",  # aromatic P: phosphole-type 5-membered ring
        "pf",  # aromatic P, alternating partner of pe
        "se",  # aromatic S: thiophene-type 5-membered ring
        "sf",  # aromatic S, alternating partner of se
    }
)

# ff14SB / ff19SB protein and OL3 / OL15 / bsc1 nucleic-acid base types --
# all uppercase.  Base carbons and nitrogens share the same type strings as
# the corresponding protein sidechain aromatic types.
_AROMATIC_PROTEIN_NUCLEIC: frozenset[str] = frozenset(
    {
        # -- sp2 aromatic carbons ----------------------------------------------
        "CA",  # Phe/Tyr/Trp/His ring C; nucleobase ring C (Ade, Gua, Cyt, Ura)
        "CB",  # ring-junction C: Trp C3a; Ade C4, Gua C4 (5/6-ring junction)
        "CC",  # His ring C (HIE and HID tautomers)
        "CD",  # His ring C, conjugated partner of CC
        "CK",  # 5-membered ring C: His C2; Ade C8, Gua C8
        "CM",  # pyrimidine 6-ring C: Cyt C5/C6, Ura C5/C6, Thy C5/C6
        "CN",  # Trp ring-junction C between 5- and 6-membered rings
        "CQ",  # 6-membered purine ring C: Ade C2, Gua C2
        "CR",  # protonated His (HIP) ring C
        "CV",  # epsilon-protonated His (HIE) ring C
        "CW",  # delta-protonated His (HID) ring C; Trp 5-membered ring C
        "C*",  # Trp ring-junction C (retained from legacy ff94 / ff99SB types)
        # -- aromatic / sp2 nitrogens ------------------------------------------
        "NA",  # pyrrole-type N: His N-H (HID, HIP), Trp N1; base N with H
        "NB",  # pyridine-type N: His lone-pair N (HIE, HID); purine ring N
        "NC",  # aromatic sp2 N in 6-ring: Cyt N3, Ade N1/N3, Gua N3/N7
        "N*",  # sp2 glycosidic N: Ade N9, Gua N9, Cyt N1, Ura N1, Thy N1
    }
)

AROMATIC_TYPES: frozenset[str] = _AROMATIC_GAFF | _AROMATIC_PROTEIN_NUCLEIC
"""
Union of aromatic atom types across GAFF, GAFF2, ff14SB, ff19SB,
OL3, OL15, and bsc1.

Atom-type strings are compared case-sensitively; GAFF/GAFF2 entries are
lowercase, ff14SB/ff19SB/nucleic-acid entries are uppercase.  LIPID17
aromatic types (where present, e.g. in cholesterol C5=C6) use GAFF-derived
types and are therefore covered by _AROMATIC_GAFF.  GLYCAM06j-1 does not
define aromatic types for standard carbohydrates.
"""


# ==========================================================================
# sp and sp2 atom types with multiple-bond (double or triple) character
# ==========================================================================

# GAFF and GAFF2 -- all lowercase.
# Notably absent: "os" (sp3 ether/ester bridging O, formerly included in
# error), "oh" (sp3 hydroxyl O), "o3" (sp3 O), "ss" (sp3 thioether S),
# "sh" (sp3 thiol S), and all sp3 carbon types (c3 etc.).
_MULTIPLE_GAFF: frozenset[str] = frozenset(
    {
        # -- sp2 / sp carbons --------------------------------------------------
        "c",  # sp2 carbonyl C: amides, esters, carboxylic acids, ketones
        "c1",  # sp C: internal alkynes (C#C)
        "c2",  # sp2 alkene C: isolated C=C (lipid double bonds, dehydroamino acids)
        "ce",  # conjugated sp2 C: dienes, enones, vinyl systems
        "cf",  # conjugated sp2 C, alternating partner of ce
        "cg",  # sp C: terminal alkynes (HC#C-)
        # -- sp2 oxygens -------------------------------------------------------
        "o",  # sp2 carbonyl O: C=O in amides, esters, ketones, aldehydes
        "o2",  # sp2 carboxylate / phosphate O: resonance-delocalized C-O(-), P=O
        # -- sp2 / sp nitrogens ------------------------------------------------
        "n",  # planar tertiary amide N: C(=O)-N (no H, e.g. N-methylamide)
        "nh",  # planar secondary amide N: C(=O)-NH (peptide bond in GAFF)
        "n2",  # sp2 imine / guanidinium N: C=N, Arg-like resonance
        "n1",  # sp nitrile N: C#N
        # -- sp2 sulfur --------------------------------------------------------
        "s2",  # sp2 thioketone / thioamide S: C=S
    }
)

# ff14SB / ff19SB protein and OL3 / OL15 / bsc1 nucleic acids -- uppercase.
_MULTIPLE_PROTEIN_NUCLEIC: frozenset[str] = frozenset(
    {
        # -- sp2 carbons -------------------------------------------------------
        "C",  # sp2 carbonyl C: protein backbone C=O; Asn/Gln/Asp/Glu side chains;
        #                 nucleobase carbonyls (Ura C2/C4, Gua C6, Cyt C2)
        # -- sp2 oxygens -------------------------------------------------------
        "O",  # sp2 carbonyl O: backbone and side-chain C=O
        "O2",  # sp2 carboxylate O: Asp/Glu side chains, C-terminus, nucleobase C=O
        # -- sp2 nitrogens -----------------------------------------------------
        "N",  # planar amide N: backbone and side-chain amide (Asn, Gln); also
        #                 in GLYCAM06j-1 N-acetyl groups (same uppercase type)
        "N2",  # sp2 guanidinium N: Arg delta/eta NH2 groups (C-N resonance in
        #                    the guanidinium group, partial double-bond char)
    }
)

# GLYCAM06j-1 -- mixed case.  The standard carbonyl / carboxylate types
# in GLYCAM (uronic acids, sialic acid, N-acetyl groups) use uppercase "C",
# "O", "O2", and "N", which are already covered by _MULTIPLE_PROTEIN_NUCLEIC.
# Only a GLYCAM-specific sp2 type that has no protein-FF counterpart is
# listed here.
_MULTIPLE_GLYCAM: frozenset[str] = frozenset(
    {
        "Cj",  # sp2 C in unsaturated modified sugar rings (e.g. 2,3-dehydro sugars)
    }
)

MULTIPLE_BOND_TYPES: frozenset[str] = (
    _MULTIPLE_GAFF | _MULTIPLE_PROTEIN_NUCLEIC | _MULTIPLE_GLYCAM
)
"""
Union of sp and sp2 atom types across GAFF, GAFF2, ff14SB, ff19SB,
OL3/OL15/bsc1, LIPID17, and GLYCAM06j-1 that predominantly participate in
double or triple bonds.

LIPID17 sp2 types are GAFF-derived (c, c2, o, o2, n) and thus covered by
_MULTIPLE_GAFF.  GLYCAM06j-1 standard carbonyl/carboxylate types ("C", "O",
"O2", "N") are uppercase and covered by _MULTIPLE_PROTEIN_NUCLEIC.
"""


# ==========================================================================
# Amide-bond atom-type pairs  C(=O)-N
# ==========================================================================

AMIDE_C: frozenset[str] = frozenset(
    {
        "c",  # GAFF / GAFF2 / LIPID17: sp2 carbonyl C
        "C",  # ff14SB / ff19SB / GLYCAM06j-1: sp2 carbonyl C
    }
)
"""
Carbonyl carbon types that participate in amide C(=O)-N bonds.

The amide predicate is applied as an explicit match for this high-barrier
bond class (~15-20 kcal/mol) in addition to the is_multiple_like() check,
providing a force-field-agnostic safeguard for edge cases where the
carbonyl carbon may have been assigned a non-standard type code.
"""

AMIDE_N: frozenset[str] = frozenset(
    {
        "n",  # GAFF / GAFF2: tertiary planar amide N (no H, e.g. N-methylamide)
        "nh",  # GAFF / GAFF2: secondary planar amide N (N-H, e.g. peptide bond)
        "N",  # ff14SB / ff19SB / GLYCAM06j-1: backbone and side-chain amide N
        #   (planar sp2; includes backbone NH and Asn/Gln side-chain NH2)
    }
)
"""
Amide nitrogen types that form C(=O)-N bonds.

Rotational barriers for representative amide bonds:
  backbone peptide bond C(=O)-N   : ~15-20 kcal/mol
  N-methylamide C(=O)-N(CH3)      : ~17 kcal/mol
  primary amide C(=O)-NH2         : ~12-15 kcal/mol

The Arg guanidinium nitrogens (type N2, uppercase) are NOT listed here
because they are better described as sp2 C=N resonance structures; they
are captured instead by is_multiple_like() via MULTIPLE_BOND_TYPES.
"""


# ==========================================================================
# Bond-level predicates
# ==========================================================================


def is_aromatic(a1: str, a2: str) -> bool:
    """
    Return True if both atom-type strings belong to AROMATIC_TYPES.

    Parameters
    ----------
    a1 : str
        Canonical (case-sensitive) atom-type string of the first atom.
    a2 : str
        Canonical (case-sensitive) atom-type string of the second atom.

    Returns
    -------
    bool

    Notes
    -----
    This predicate returns True for any pair of atoms that both carry
    aromatic type codes, including atoms in different aromatic rings
    connected by an inter-ring single bond (e.g. biphenyl, biaryl
    heterocycles).  Such bonds are flexible in solution
    (barrier ~1-3 kcal/mol) and represent false positives.  See the
    module-level Known Limitations section.  Callers requiring strict
    intra-ring validation must augment this predicate with an external
    ring-membership check.
    """
    return a1 in AROMATIC_TYPES and a2 in AROMATIC_TYPES


def is_amide(a1: str, a2: str) -> bool:
    """
    Return True if the atom-type pair matches the canonical amide pattern.

    An amide bond (C(=O)-N) is classified as rigid on the basis of its
    substantial rotational barrier (~15-20 kcal/mol), which arises from
    resonance delocalisation of the nitrogen lone pair into the carbonyl
    pi system, conferring partial double-bond character on the C-N bond.

    Parameters
    ----------
    a1 : str
        Canonical atom-type string of the first atom.
    a2 : str
        Canonical atom-type string of the second atom.

    Returns
    -------
    bool
        True for any ordering of an AMIDE_C type and an AMIDE_N type.
    """
    return (a1 in AMIDE_C and a2 in AMIDE_N) or (a2 in AMIDE_C and a1 in AMIDE_N)


def is_multiple_like(a1: str, a2: str) -> bool:
    """
    Return True if both atom types are consistent with sp or sp2 hybridisation.

    Both atoms must appear in MULTIPLE_BOND_TYPES.  This catch-all predicate
    covers multiple-bond character not handled by the aromatic or amide
    predicates: isolated alkenes and alkynes, carbonyls in non-amide contexts
    (esters, ketones, aldehydes), carboxylates, thioketones, nitriles,
    conjugated imines, and guanidinium C-N bonds.

    Parameters
    ----------
    a1 : str
        Canonical atom-type string of the first atom.
    a2 : str
        Canonical atom-type string of the second atom.

    Returns
    -------
    bool

    Notes
    -----
    Requiring both partners to carry sp/sp2 type codes is intentionally
    conservative.  A bond between an sp2 carbonyl carbon and an sp3
    substituent (e.g. C(=O)-CH3) is correctly NOT flagged, because the
    sp3 carbon type (c3, CT, Cg, etc.) is absent from MULTIPLE_BOND_TYPES.
    This matches the physical situation: the C-CH3 bond is freely rotatable
    even though the carbon on one side participates in a C=O double bond.
    """
    return a1 in MULTIPLE_BOND_TYPES and a2 in MULTIPLE_BOND_TYPES


def is_rigid_bond(parent_atom: pmd.Atom, child_atom: pmd.Atom) -> bool:
    """
    Determine whether the bond between two parmed.Atom objects is rigid.

    Applies three heuristic criteria in priority order after an optional
    fast-path bond-order check:

      0. Explicit bond order >= 2 (fast path).
         Used when parmed has resolved bond orders (e.g. from an attached
         RDKit molecule or SYBYL MOL2 input).  Bypasses all heuristics.

      1. Aromatic: both atoms carry types in AROMATIC_TYPES.
         Covers intra-ring bonds in all aromatic systems defined by the
         supported force fields.  Note the biaryl false-positive caveat
         in the module docstring.

      2. Amide: the type pair matches AMIDE_C x AMIDE_N.
         Explicit match for backbone and side-chain amide bonds,
         supplementing criterion 3 for robustness.

      3. Multiple-bond-like: both atoms carry types in MULTIPLE_BOND_TYPES.
         Covers C=O, C=C, C=N, C#N, C#C, C=S, and related bonds across
         all supported force fields.

    Parameters
    ----------
    parent_atom : parmed.Atom
        One endpoint of the bond.
    child_atom : parmed.Atom
        The other endpoint of the bond.

    Returns
    -------
    bool
        True if the bond is classified as rigid (non-rotatable).

    Raises
    ------
    AttributeError
        If either atom lacks a ``type`` attribute, indicating a malformed
        or incompletely initialised parmed.Atom object.

    Notes
    -----
    Approximate rotational barriers for representative bond classes:

      Bond class                   Example                Barrier (kcal/mol)
      ------------------------------------------------------------------
      Aromatic intra-ring C-C      benzene                ~30 (locked)
      Aromatic intra-ring C-N      pyridine C2-N1         ~25 (locked)
      Amide C(=O)-N                backbone peptide bond  15-20
      sp2 carbonyl C=O             ester, ketone          >40 (locked)
      sp2 alkene C=C               isolated double bond   ~60 (locked)
      sp2 imine C=N (guanidinium)  Arg Czeta-Neta         ~10-15
      Biaryl inter-ring C-C        biphenyl               ~1-3 (rotatable)
      sp3 C-C (ester C-OS)         ethyl ester            ~3-5 (rotatable)

    The last two entries are correctly classified as rotatable:
    - Biaryl: flagged as rigid by criterion 1 (known limitation).
    - Ester C-OS: the sp3 OS type is absent from all sets; not flagged.

    References
    ----------
    .. [1] Wang, J. et al. (2004) J. Comput. Chem. 25, 1157-1174.
    .. [2] He, X. et al. (2020) J. Chem. Theory Comput. 16, 528-552.
    .. [3] Maier, J. A. et al. (2015) J. Chem. Theory Comput. 11, 3696-3713.
    .. [4] Tian, C. et al. (2020) J. Chem. Theory Comput. 16, 528-552.
    .. [5] Zgarbova, M. et al. (2011) J. Chem. Theory Comput. 7, 2886-2902.
    .. [6] Zgarbova, M. et al. (2015) J. Chem. Theory Comput. 11, 5723-5736.
    .. [7] Ivani, I. et al. (2016) Nat. Methods 13, 55-58.
    .. [8] Gould, I. R. et al. (2018) J. Phys. Chem. B 122, 5573-5592.
    .. [9] Kirschner, K. N. et al. (2008) J. Comput. Chem. 29, 622-655.
    """
    # Fast path: use explicit bond order when available.
    bond = _find_bond(parent_atom, child_atom)
    if bond is not None and bond.order is not None and bond.order >= 2:
        return True

    # Atom types are compared case-sensitively; do NOT apply .lower().
    a1 = parent_atom.type
    a2 = child_atom.type

    if is_aromatic(a1, a2):
        return True

    if is_amide(a1, a2):
        return True

    if is_multiple_like(a1, a2):
        return True

    return False


# ==========================================================================
# Private helpers
# ==========================================================================


def _find_bond(a1: pmd.Atom, a2: pmd.Atom) -> "pmd.Bond | None":
    """
    Return the parmed.Bond connecting a1 and a2, or None if not found.

    Iterates over the bond list of a1 (complexity O(degree(a1))) and
    checks object identity for both endpoints.

    Parameters
    ----------
    a1 : parmed.Atom
    a2 : parmed.Atom

    Returns
    -------
    parmed.Bond or None
    """
    for bond in a1.bonds:
        if bond.atom1 is a2 or bond.atom2 is a2:
            return bond
    return None
