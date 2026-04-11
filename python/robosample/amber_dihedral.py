from typing import Dict, List, Tuple

import parmed as pmd

PROTEIN_BACKBONE = {
    ("C", "N", "CA", "C"): "phi",
    ("N", "CA", "C", "N"): "psi",
    ("CA", "C", "N", "CA"): "omega",  # Standard peptide bond
    ("O", "C", "N", "CA"): "omega",  # Carbonyl Oxygen variant
    ("CH3", "C", "N", "CA"): "omega",  # N-terminal methylation (ACE)
    ("CA", "C", "N", "C"): "omega",  # C-terminal methylation (NME)
    ("O", "C", "N", "CH3"): "omega",  # Carbonyl Oxygen to NME cap
}

# For pretty images: https://emleddin.github.io/comp-chem-website/AMBERguide-AAs-DNA-RNA.html
PROTEIN_SIDECHAIN = {
    # N-terminus ACE
    "ACE": {},
    # C-terminus NME
    "NME": {},
    # Alanine does not have any side chain dihedrals since its definition requires four heavy atoms
    "ALA": {},
    # Arginine (ARG) has a terminal large, resonance-stabilized, planar structure with a diffuse positive charge guanidinium group
    # CZ is connected to two nitrogen atoms (NH1, NH2) which are connected to two hydrogens each (HH11, HH12 and HH21, HH22)
    # The actual guanidinium group is (HH11, HH12, NH1), CZ, (HH21, HH22, NH2) and chi5 allows rotation of this entire group around the NE-CZ bond
    "ARG": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "CD"),
        "chi3": ("CB", "CG", "CD", "NE"),
        "chi4": ("CG", "CD", "NE", "CZ"),
        "chi5": ("CD", "NE", "CZ", "NH1"),
    },
    # Asparagine (ASN)
    "ASN": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "OD1"),
    },
    # Aspartate (deprotonated, -1)
    "ASP": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "OD1"),
    },
    # Aspartic Acid (protonated, neutral)
    "ASH": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "OD1"),
    },
    # Cysteine (CYS) has a terminal thiol group (HG connected to SG)
    # If HG is lost, it can become CYX and form disulfide bonds or CYM and coordinate metals
    # IUPAC defines only chi1 (heavy atoms). HG rotation is a torsion, not a chi angle
    "CYS": {
        "chi1": ("N", "CA", "CB", "SG"),
    },
    # Cystine (Disulfide-bonded, neutral)
    # When bonded, the SG-SG bond creates a chi2 angle reaching into the partner residue
    "CYX": {
        "chi1": ("N", "CA", "CB", "SG"),
        "ring_closing": ("CB", "SG", "SG", "CB"),  # disulfide bridge
    },
    # Deprotonated Cysteine (Anionic or Metal-coordinated)
    # Metal coordination (e.g., Zn2+, Fe-S clusters) creates a rigid geometry,
    # but it is usually modeled as an external coordinate constraint, not a side-chain chi.
    "CYM": {
        "chi1": ("N", "CA", "CB", "SG"),
    },
    # Glutamine (GLN)
    # chi3 is defined by OE1 per IUPAC convention for amide groups.
    # Note: chi4 (hydrogen rotation) is omitted as per heavy-atom standards.
    "GLN": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "CD"),
        "chi3": ("CB", "CG", "CD", "OE1"),
    },
    # Glutamate (Deprotonated, -1)
    # chi3 is defined by OE1 per IUPAC convention for carboxylate symmetry.
    "GLU": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "CD"),
        "chi3": ("CB", "CG", "CD", "OE1"),
    },
    # Glutamic Acid (Protonated, neutral)
    # We treat the heavy-atom chi3 as the terminal dihedral.
    # Note: Does not include the rotation of the hydroxyl proton (HO-C-C-C).
    # Conventionally, OE1 is the carbonyl oxygen.
    "GLH": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "CD"),
        "chi3": ("CB", "CG", "CD", "OE1"),
    },
    # Glycine has no side chain
    "GLY": {},
    # Imidazole (HID/HIE/HIP) is aromatic and essentially planar
    "HIP": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "ND1"),
        "ring_closing_1": ("CG", "ND1", "CE1", "NE2"),
    },
    "HIE": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "ND1"),
        "ring_closing_1": ("CG", "ND1", "CE1", "NE2"),
    },
    "HID": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "ND1"),
        "ring_closing_1": ("CG", "ND1", "CE1", "NE2"),
    },
    # Isoleucine (ILE) is branched at CB and has a chiral center at the C3 (beta) position.
    "ILE": {
        "chi1": ("N", "CA", "CB", "CG1"),
        "chi2": ("CA", "CB", "CG1", "CD1"),
    },
    # Leucine (LEU)  is branched at the CG (gamma) carbon.
    "LEU": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "CD1"),
    },
    # Lysine (LYS, protonated, charge +1)
    # chi1-chi4 follow the heavy-atom chain to the terminal nitrogen.
    # IUPAC/PDB standards do not define chi5 (hydrogen rotation).
    "LYS": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "CD"),
        "chi3": ("CB", "CG", "CD", "CE"),
        "chi4": ("CG", "CD", "CE", "NZ"),
    },
    # Lysine (LYN, deprotonated, neutral)
    # Heavy-atom dihedrals remain the same as the protonated state.
    "LYN": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "CD"),
        "chi3": ("CB", "CG", "CD", "CE"),
        "chi4": ("CG", "CD", "CE", "NZ"),
    },
    # Methionine (MET) features a thioether linkage.
    "MET": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "SD"),
        "chi3": ("CB", "CG", "SD", "CE"),
    },
    # Phenylalanine (PHE)
    "PHE": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "CD1"),
        "ring_closing_1": ("CG", "CD1", "CE1", "CZ"),  # Same as TYR
    },
    # Proline (PRO)
    "PRO": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "CD"),
        # These dihedrals are unique to proline, they will not be picked up via other residues
        # I list them them here for completeness
        # 'chi3': ('CB', 'CG', 'CD', 'N'),
        "ring_closing": ("CG", "CD", "N", "CA"),
        # Topologies intentionally skip dihedrals that traverse the ring closure bond
        # Including it would create a dihedral that is linearly dependent on others already defined in the ring
        # Thus, this is absent and will never be detected
        # I include it here for completeness and to prove the point
        # 'pro-chi5': ('CD', 'N', 'CA', 'C'),
    },
    # Serine (SER)
    "SER": {
        "chi1": ("N", "CA", "CB", "OG"),
    },
    # Threonine (THR) is branched at CB and has a chiral center at the C3 (beta) position.
    "THR": {
        "chi1": ("N", "CA", "CB", "OG1"),
    },
    # Tryptophan (TRP)
    # The side chain is a rigid, planar bicyclic indole group.
    "TRP": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "CD1"),
        "ring_closing_1": ("CG", "CD1", "NE1", "CE2"),  # Pyrrole ring
        "ring_closing_2": ("CE2", "CZ2", "CH2", "CZ3"),  # Benzene ring
    },
    # Tyrosine (TYR)
    # The side chain is a rigid, planar aromatic ring with a hydroxyl group.
    "TYR": {
        "chi1": ("N", "CA", "CB", "CG"),
        "chi2": ("CA", "CB", "CG", "CD1"),
        "ring_closing_1": ("CG", "CD1", "CE1", "CZ"),  # Same as PHE
    },
    # Valine (VAL) is branched at the CB (beta) carbon with two methyl groups.
    "VAL": {
        "chi1": ("N", "CA", "CB", "CG1"),
    },
}

GLYCAN_PROTEIN_LINKAGE = {
    # O-Glycosylated Serine (OLS)
    "OLS": {
        "phi": ("CA", "CB", "OG", "C1"),  # C1 is the first carbon of the sugar
        "psi": ("CB", "OG", "C1", "C2"),
    },
    # O-Glycosylated Threonine (OLT)
    "OLT": {
        "phi": ("CA", "CB", "OG1", "C1"),
        "psi": ("CB", "OG1", "C1", "C2"),
    },
    # N-Glycosylated Asparagine (NLN)
    "NLN": {
        "phi": ("CB", "CG", "ND2", "C1"),
        "psi": ("CG", "ND2", "C1", "C2"),
    },
    # Hydroxyproline (HYP) - if glycosylated at OD1
    "HYP": {
        "phi": ("CB", "CG", "OD1", "C1"),
        "psi": ("CG", "OD1", "C1", "C2"),
    },
}

# Canonical pyranose ring torsions (tau1-tau6)
PYRANOSE_TAU_PATTERNS = {
    ("O5", "C1", "C2", "C3"): "tau1",
    ("C1", "C2", "C3", "C4"): "tau2",
    ("C2", "C3", "C4", "C5"): "tau3",
    ("C3", "C4", "C5", "O5"): "tau4",
    ("C4", "C5", "O5", "C1"): "tau5",
    ("C5", "O5", "C1", "C2"): "tau6",
}

FURANOSE_TAU_PATTERNS = {
    ("O4", "C1", "C2", "C3"): "tau0",
    ("C1", "C2", "C3", "C4"): "tau1",
    ("C2", "C3", "C4", "O4"): "tau2",
    ("C3", "C4", "O4", "C1"): "tau3",
    ("C4", "O4", "C1", "C2"): "tau4",
}

NUCLEIC_DIHEDRALS = {
    # Backbone torsions
    ("O3'", "P", "O5'", "C5'"): "alpha",
    ("P", "O5'", "C5'", "C4'"): "beta",
    ("O5'", "C5'", "C4'", "C3'"): "gamma",
    ("C5'", "C4'", "C3'", "O3'"): "delta",
    ("C4'", "C3'", "O3'", "P"): "epsilon",
    ("C3'", "O3'", "P", "O5'"): "zeta",
    # Glycosidic torsion (different base atoms for purines vs pyrimidines)
    ("O4'", "C1'", "N9", "C4"): "chi_purine",
    ("O4'", "C1'", "N1", "C2"): "chi_pyrimidine",
    # Sugar ring torsions (for pseudorotation analysis)
    ("C4'", "O4'", "C1'", "C2'"): "nu0",
    ("O4'", "C1'", "C2'", "C3'"): "nu1",
    ("C1'", "C2'", "C3'", "C4'"): "nu2",
    ("C2'", "C3'", "C4'", "O4'"): "nu3",
    ("C3'", "C4'", "O4'", "C1'"): "nu4",
}


def _is_dihedral_intraresidue(dihedral: pmd.topologyobjects.Dihedral) -> bool:
    resids = (
        dihedral.atom1.residue.idx,
        dihedral.atom2.residue.idx,
        dihedral.atom3.residue.idx,
        dihedral.atom4.residue.idx,
    )
    return len(set(resids)) == 1


def _extract_residue_indices(
    atoms: List[pmd.topologyobjects.Atom],
) -> Tuple[int, int, int, int]:
    resids = (
        atoms[0].residue.idx,
        atoms[1].residue.idx,
        atoms[2].residue.idx,
        atoms[3].residue.idx,
    )
    if not all(resid != -1 for resid in resids):
        raise ValueError(
            f"One or more atoms in the dihedral do not have a valid residue index: {resids}"
        )
    return resids


def _safe_insert(
    dictionary: Dict[Tuple[str, str, str, str], str],
    key: Tuple[str, str, str, str],
    value: str,
):
    if key in dictionary:
        if dictionary[key] != value:
            raise ValueError(
                f"Conflict in dihedral definitions for {value}: {dictionary[key]} vs {key}"
            )
    else:
        dictionary[key] = value


def _extract_atom_names(
    atoms: List[pmd.topologyobjects.Atom],
) -> Tuple[str, str, str, str]:
    return tuple(atom.name for atom in atoms)


def _has_nucleic_sugar(atoms: List[pmd.topologyobjects.Atom]) -> bool:
    for atom in atoms:
        if "'" in atom.name:
            return True
    return False


def _extract_carbon_index(name: str) -> int | None:
    if not name.startswith("C"):
        return None
    try:
        return int(name[1:])
    except ValueError:
        return None


def _extract_oxygen_index(name: str) -> int | None:
    if not name.startswith("O"):
        return None
    try:
        return int(name[1:])
    except ValueError:
        return None


class DihedralClassifier:
    def __init__(self):
        # Lookup table for dihedral definitions: keys are tuples of 4 atom types, values are dihedral labels
        # We need to flatten the protein backbone and side chain definitions into a single dictionary for efficient lookup
        self.protein_dihedral_definitions: Dict[Tuple[str, str, str, str], str] = {}

        # Protein backbone definitions
        for atom_types, dihedral_name in PROTEIN_BACKBONE.items():
            _safe_insert(self.protein_dihedral_definitions, atom_types, dihedral_name)

        # Protein side chain definitions
        for residue, dihedral_list in PROTEIN_SIDECHAIN.items():
            for dihedral_name, atom_types in dihedral_list.items():
                _safe_insert(
                    self.protein_dihedral_definitions, atom_types, dihedral_name
                )

        # Glycan-protein linkage definitions
        self.glycan_linkage_definitions: Dict[Tuple[str, str, str, str], str] = {}
        for residue, dihedral_list in GLYCAN_PROTEIN_LINKAGE.items():
            for dihedral_name, atom_types in dihedral_list.items():
                _safe_insert(self.glycan_linkage_definitions, atom_types, dihedral_name)

    def _classify_nucleic_acid(
        self, dihedral: pmd.topologyobjects.Dihedral
    ) -> str | None:
        # Check if the dihedral matches any of the defined nucleic acid dihedrals (backbone, glycosidic, or sugar)
        # We check both the forward and reverse order of the atoms to account for different orientations
        atoms = [dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4]
        label = NUCLEIC_DIHEDRALS.get(_extract_atom_names(atoms))
        if not label:
            atoms = list(reversed(atoms))
            label = NUCLEIC_DIHEDRALS.get(_extract_atom_names(atoms))

        return label

    def _classify_protein(
        self, gparent: pmd.Atom, parent: pmd.Atom, child: pmd.Atom, gchild: pmd.Atom
    ) -> str | None:

        # Check if the dihedral matches any of the defined protein dihedrals (backbone or side chain)
        # We check both the forward and reverse order of the atoms to account for different orientations
        atoms = (gparent, parent, child, gchild)
        label = self.protein_dihedral_definitions.get(_extract_atom_names(atoms))
        if not label:
            atoms = atoms[::-1]
            label = self.protein_dihedral_definitions.get(_extract_atom_names(atoms))

        resids = _extract_residue_indices(atoms)
        i = resids[1]

        # Phi: C(i-1) N(i) CA(i) C(i)
        if label == "phi":
            if (
                resids[0] == i - 1
                and resids[1] == i
                and resids[2] == i
                and resids[3] == i
            ):
                return label
            else:
                return None

        # Psi: N(i) CA(i) C(i) N(i+1)
        if label == "psi":
            if (
                resids[0] == i
                and resids[1] == i
                and resids[2] == i
                and resids[3] == i + 1
            ):
                return label
            else:
                return None

        # Omega: CA(i) C(i) N(i+1) CA(i+1)
        if label == "omega":
            if (
                resids[0] == i
                and resids[1] == i
                and resids[2] == i + 1
                and resids[3] == i + 1
            ):
                return label
            else:
                return None

        return label

    def _classify_glycan_linkage(
        self, dihedral: pmd.topologyobjects.Dihedral
    ) -> str | None:

        # Check if the dihedral matches any of the defined glycan-protein linkage dihedrals
        # We check both the forward and reverse order of the atoms to account for different orientations
        atoms = [dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4]
        label = self.glycan_linkage_definitions.get(_extract_atom_names(atoms))
        if not label:
            atoms = list(reversed(atoms))
            label = self.glycan_linkage_definitions.get(_extract_atom_names(atoms))

        return label

    def _classify_endocyclic_tau(
        self, dihedral: pmd.topologyobjects.Dihedral
    ) -> str | None:
        """
        Unified classifier for Pyranose (tau1-6) and Furanose (tau0-4).
        """
        if not _is_dihedral_intraresidue(dihedral):
            return None

        atoms = [dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4]
        atom_names = tuple(_extract_atom_names(atoms))
        rev_names = tuple(reversed(atom_names))

        # 1. Check Pyranose Patterns
        label = PYRANOSE_TAU_PATTERNS.get(atom_names) or PYRANOSE_TAU_PATTERNS.get(
            rev_names
        )
        if label:
            return f"pyranose-{label}"

        # 2. Check Furanose Patterns
        label = FURANOSE_TAU_PATTERNS.get(atom_names) or FURANOSE_TAU_PATTERNS.get(
            rev_names
        )
        if label:
            return f"furanose-{label}"

        return None

    def _is_glycosidic_phi(self, dihedral: pmd.topologyobjects.Dihedral) -> bool:
        """
        Pattern: O_ring(i) - C_anomeric(i) - L_link(j) - C_next(j)
        """

        def check(atoms: List[pmd.topologyobjects.Atom]) -> bool:

            # Phi spans across two residues: donor (sugar i) and acceptor (sugar j)
            resids = _extract_residue_indices(atoms)

            # First two atoms (O_ring, C_anomeric) belong to the donor (i)
            if not (resids[0] == resids[1]):
                return False

            # Last two atoms (L_link, C_next) belong to the acceptor (j)
            if not (resids[2] == resids[3]):
                return False

            # Donor and acceptor must be different residues and are not necessarily adjacent in sequence
            if resids[1] == resids[2]:
                return False

            # Now we try to match atom names
            names = _extract_atom_names(atoms)
            if _has_nucleic_sugar(atoms):
                return False

            # 1. Ring oxygen (O_ring)
            # Pyranose (O5), Furanose (O4), or Sialic Acid/Neu5Ac (O6)
            if names[0] not in ("O4", "O5", "O6"):
                return False

            # 2. Anomeric carbon (C_anomeric)
            # Aldoses use C1, ketoses (sialic acid, fructose) use C2
            if names[1] not in ("C1", "C2"):
                return False

            # Validation: Sialic acid specific check (O6 must pair with C2)
            if names[0] == "O6" and names[1] != "C2":
                return False

            # Validation: Standard pyranose/furanose (O5/O4 usually pair with C1)
            # Note: Fructose is an exception (O4-C2), so we allow the check to pass.

            # 3. Linkage atom (L_link)
            # Usually Oxygen (O), but can be Nitrogen (N) for N-linked glycans or Sulfur (S) for S-glycosides
            # GLYCAM names these by their linkage.
            # We check if it starts with the element symbol.
            if not any(names[2].startswith(e) for e in ("O", "N", "S")):
                return False

            # 4. Acceptor carbon (C_next)
            # Must be a carbon on the next sugar
            # In GLYCAM06j, these are C1 through C6 (or rarely C7-C9 in sialic chains)
            cx = _extract_carbon_index(names[3])
            if cx is None or not (1 <= cx <= 9):
                return False

            return True

        # Evaluate both directions
        fwd = [dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4]
        rev = list(reversed(fwd))
        return check(fwd) or check(rev)

    def _is_glycosidic_psi(self, dihedral: pmd.topologyobjects.Dihedral) -> bool:
        """
        Pattern: C_anomeric(i) | L_link(j) - C_next(j) - C_neighbor(j)
        """

        def check(atoms: List[pmd.topologyobjects.Atom]) -> bool:
            resids = _extract_residue_indices(atoms)
            names = _extract_atom_names(atoms)
            if _has_nucleic_sugar(atoms):
                return False

            if resids[0] == resids[1]:
                return False
            if not (resids[1] == resids[2] == resids[3]):
                return False

            # 1. Anomeric Carbon (Donor i)
            if names[0] not in ("C1", "C2"):
                return False

            # 2. Linkage Atom (Acceptor j)
            if not any(names[1].startswith(e) for e in ("O", "N", "S")):
                return False

            # 3. Acceptor Carbon (Acceptor j)
            cx = _extract_carbon_index(names[2])
            if cx is None:
                return False

            # 4. Neighbor Carbon (Acceptor j)
            # Usually Cx-1 or Cx+1. We just ensure it's a ring/backbone carbon.
            nx = _extract_carbon_index(names[3])
            if nx is None or not (1 <= nx <= 9):
                return False

            # Validation: Acceptor Carbon and Linkage Oxygen should match index
            # e.g., if names[1] is 'O4', names[2] should be 'C4'
            ox_idx = _extract_oxygen_index(names[1])
            if ox_idx is not None and ox_idx != cx:
                return False

            return True

        # Evaluate both directions
        fwd = [dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4]
        rev = list(reversed(fwd))
        return check(fwd) or check(rev)

    def _is_glycosidic_omega(self, dihedral: pmd.topologyobjects.Dihedral) -> bool:
        """
        Refined detector for true glycosidic omega (hinge) torsions.
        Pattern: O_link(j) - C_exo(j) - C_ring(j) - O_ring(j)
        """

        def check(atoms) -> bool:
            # 1. Intra-residue check: All atoms must be on the 'Acceptor' residue
            if not all(a.residue.idx == atoms[0].residue.idx for a in atoms):
                return False

            # 2. Linkage Bridge check: Atom1 must connect to a DIFFERENT residue
            linkage_atom = atoms[0]
            if not any(
                p.residue.idx != linkage_atom.residue.idx
                for p in linkage_atom.bond_partners
            ):
                return False

            names = [a.name for a in atoms]
            if _has_nucleic_sugar(atoms):
                return False

            # 3. Structural Pattern:
            # Atom 0: Linkage Oxygen (O5-O9 for Sialic/Higher sugars, O6 for Hexose)
            if not any(names[0].startswith(o) for o in ("O5", "O6", "O7", "O8", "O9")):
                return False

            # Atom 3: Ring Oxygen (O4, O5, or O6)
            if names[3] not in ("O4", "O5", "O6"):
                return False

            # Atom 1 & 2: Carbon sequence (Exocyclic -> Ring Neighbor)
            # Typically C6 -> C5 (Hexose) or C8 -> C7 (Sialic Acid)
            c_exo = _extract_carbon_index(names[1])
            c_ring = _extract_carbon_index(names[2])

            if c_exo is None or c_ring is None:
                return False

            # Ensure they are adjacent carbons (the exocyclic swing)
            return abs(c_exo - c_ring) == 1

        # Evaluate both directions
        fwd = [dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4]
        rev = list(reversed(fwd))
        return check(fwd) or check(rev)

    def _classify_exocyclic(self, dihedral: pmd.topologyobjects.Dihedral) -> str | None:
        """
        Classifies non-backbone glycan dihedrals (hydroxyls, anomeric substituents,
        and functional groups like sulfates/phosphates).
        """

        def check(atoms, names) -> str | None:
            a1, a2, a3, a4 = names
            if _has_nucleic_sugar(atoms):
                return None

            # 1. Chi (Hydroxyl/Exocyclic): ends in a pendant Oxygen
            # We exclude the ring oxygens (O4/O5) to avoid catching ring puckers.
            if a4.startswith("O") and a4 not in ("O4", "O5"):
                # Ensure the path is Carbon-heavy (C-C-C-O)
                c_indices = [_extract_carbon_index(n) for n in names[:3]]

                if all(idx is not None for idx in c_indices):
                    # Check if the oxygen is terminal (not part of a bridge)
                    # This prevents 'omega' from being double-counted as 'chi'
                    is_bridge = any(
                        p.residue.idx != atoms[3].residue.idx
                        for p in atoms[3].bond_partners
                    )
                    if not is_bridge:
                        return "glycosidic-chi-exocyclic"

            # 2. Anomeric Substituents (C1-O1-X)
            if (
                a2 == "C1"
                and a1 in ("O5", "O4")
                and a3.startswith("O")
                and a3 not in ("O4", "O5")
            ):
                return "glycosidic-chi-anomeric"

            # 3. Sulfates/Phosphates
            if _extract_carbon_index(a1) is not None and a2.startswith("O"):
                if "S" in a3:
                    return "glycosidic-chi-sulfate"
                if "P" in a3:
                    return "glycosidic-chi-phosphate"

            return None

        # Exocyclic dihedrals are intra-residue
        if not _is_dihedral_intraresidue(dihedral):
            return None

        # Check forward
        atoms = [dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4]
        result = check(atoms, _extract_atom_names(atoms))
        if result:
            return result

        # Check reverse
        atoms = list(reversed(atoms))
        result = check(atoms, _extract_atom_names(atoms))

        return result

    def _classify_substituent(
        self, dihedral: pmd.topologyobjects.Dihedral
    ) -> str | None:
        """
        Classifies substituent torsions (N-acetyl, O-acetyl, glycerol chains).
        Works for both Pyranoses and Furanoses.
        """
        if not _is_dihedral_intraresidue(dihedral):
            return None

        atoms = [dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4]
        names = [at.name for at in atoms]
        if _has_nucleic_sugar(atoms):
            return None

        def check_substituent(at_list, nm_list):
            a1, a2, a3, a4 = nm_list
            at2, at3 = at_list[1], at_list[2]

            # 1. Acetyl/Amide Groups (N-linked or O-linked)
            # Covers GlcNAc (N2-C2N) and O-acetylation (O3-C3A)
            # Bond 2-3 is the linkage (N-C or O-C)
            if (at2.element_name in ["N", "O"]) and at3.element_name == "C":
                neighbors = [at.element_name for at in at3.bond_partners]

                # Check if atom3 is a carbonyl carbon (has a double-bonded/carbonyl Oxygen)
                if "O" in neighbors:
                    # If a4 is the Carbonyl Oxygen -> Amide/Ester plane
                    if at_list[3].element_name == "O":
                        return (
                            "glycan-chi-subst-amide"
                            if at2.element_name == "N"
                            else "glycan-chi-subst-ester"
                        )

                    # If a4 is the Methyl Carbon (CME, C7, etc) -> Acetyl rotation
                    if at_list[3].element_name == "C":
                        return "glycan-chi-acetyl"

            # Sialic Acid / Polyols (Glycerol side chains)
            # Pattern: C6-C7-C8-C9 (common in Neu5Ac)
            if all(at.element_name == "C" for at in at_list):
                # If these atoms are part of the 'tail' and not the ring
                # (Checking if atoms are NOT in ring usually requires a ring-finder,
                # but we can infer via atom names/residue types)
                if any(name in ["C7", "C8", "C9"] for name in nm_list):
                    return "glycan-chi-exocyclic-tail"

            return None

        # Check forward and backward
        res = check_substituent(atoms, names)
        if res:
            return res
        return check_substituent(list(reversed(atoms)), list(reversed(names)))

    def classify(
        self, gparent: pmd.Atom, parent: pmd.Atom, child: pmd.Atom, gchild: pmd.Atom
    ) -> str | None:

        return self._classify_protein(gparent, parent, child, gchild)

        # # Extract atom types (not atom names) to dispatch to the appropriate classifier
        # atom_types = [dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, dihedral.atom4.type]

        # # Handle protein backbone and side chains
        # if all(t in atomtypes.AMBER_FF19SB_ATOM_TYPES for t in atom_types):
        #     label = self._classify_protein(dihedral)
        #     if label:
        #         return 'protein-' + label
        #     else:
        #         return 'protein-other'

        # # Handle lipids
        # if all(t in atomtypes.AMBER_LIPID_21_ATOM_TYPES for t in atom_types):
        #     return 'lipid'

        # # Try to classify as nucleic acid first
        # # Some glycosidic dihedrals may superficially resemble nucleic acid torsions, so we check for nucleic acid residues first to avoid misclassification.
        # NUCLEIC_RESIDUES = {
        #     "DA", "DC", "DG", "DU", "DA3", "DA5", "DC3", "DC5", "DG3", "DG5", "DU3", "DU5", # DNA
        #     "A", "C", "G", "U", "A3", "A5", "C3", "C5", "G3", "G5", "U3", "U5", # RNA
        # }
        # residue_names = [dihedral.atom1.residue.name, dihedral.atom2.residue.name, dihedral.atom3.residue.name, dihedral.atom4.residue.name]
        # if all(res in NUCLEIC_RESIDUES for res in residue_names):
        #     label = self._classify_nucleic_acid(dihedral)
        #     label = 'nucleic-' + label if label else None
        #     return label

        # # Finally, check for glycan linkages and other glycan-specific dihedrals
        # label = self._classify_glycan_linkage(dihedral)
        # if label:
        #     return 'glycan-linkage-' + label

        # # Try glycosidic phi
        # if self._is_glycosidic_phi(dihedral):
        #     return 'glycosidic-phi'

        # # Try glycosidic psi
        # if self._is_glycosidic_psi(dihedral):
        #     return 'glycosidic-psi'

        # # Try glycosidic omega (hinge)
        # if self._is_glycosidic_omega(dihedral):
        #     return 'glycosidic-omega'

        # # If it's not a backbone or linkage dihedral, it might be a pyranose or furanose ring tau or an exocyclic chi
        # # We check these before general substituents to capture specific glycan features.
        # label = self._classify_endocyclic_tau(dihedral)
        # if label:
        #     return label

        # # Check for exocyclic chi angles (hydroxyls, anomeric substituents, sulfates/phosphates)
        # label = self._classify_exocyclic(dihedral)
        # if label:
        #     return label

        # # Check for other substituent dihedrals (N-acetyl, O-acetyl, glycerol tails)
        # label = self._classify_substituent(dihedral)
        # if label:
        #     return label

        return None
