import parmed as pmd

AROMATIC_TYPES = {
    # GAFF aromatic carbons / heteroatoms
    "ca",
    "cp",
    "cq",
    "cc",
    "cd",
    "ce",
    "cf",
    "na",
    "nb",
    "nc",
    "nd",
    "ne",
    "nf",
    # ff19SB aromatic
    "c",
    "ca",
    "cb",
    "cc",
    "cd",
    "ce",
    "cf",
    "cw",
    "cr",
}

# sp / sp2 carbons and heteroatoms (very conservative)
MULTIPLE_BOND_TYPES = {
    # carbonyls, sp2/sp
    "c",
    "c1",
    "c2",
    "ce",
    "cf",
    "cg",
    "o",
    "o2",
    "os",
    "n",
    "n2",
    "n1",
}

# canonical AMBER amide pattern
AMIDE_C = {"c"}  # carbonyl carbon
AMIDE_N = {"n", "nh"}  # amide nitrogens


def is_aromatic(a1: str, a2: str) -> bool:
    return (a1 in AROMATIC_TYPES) and (a2 in AROMATIC_TYPES)


def is_amide(a1: str, a2: str) -> bool:
    return (a1 in AMIDE_C and a2 in AMIDE_N) or (a2 in AMIDE_C and a1 in AMIDE_N)


def is_multiple_like(a1: str, a2: str) -> bool:
    return a1 in MULTIPLE_BOND_TYPES and a2 in MULTIPLE_BOND_TYPES


def is_rigid_bond(parent_atom: pmd.Atom, child_atom: pmd.Atom) -> bool:
    a1_type = parent_atom.type.lower()
    a2_type = child_atom.type.lower()

    if is_aromatic(a1_type, a2_type):
        return True
    if is_amide(a1_type, a2_type):
        return True
    if is_multiple_like(a1_type, a2_type):
        return True

    return False
