#!/usr/bin/env python
# Called by the plugin as a subprocess. Has nothing to do with PyMOL.

import sys
import json
import parmed as pmd

# adjust import to wherever your classifier lives
from amber_dihedral_classifier import DihedralClassifier, DihedralType

def main():
    print(sys.argv)
    pdb = sys.argv[1]
    structure = pmd.load_file(pdb)
    classifier = DihedralClassifier()

    for i, dih in enumerate(structure.dihedrals):

        atom_string = " - ".join(
            f"{a.residue.name}{a.residue.number}:{a.name}"
            for a in dih.atoms
        )

        try:
            result = classifier.classify(dih)

            # pretty printing
            if isinstance(result, DihedralType):
                label = result.name
            else:
                label = str(result)

        except Exception as e:
            label = f"ERROR({e})"

        print(f"{i:6d} | {label:20s} | {atom_string}")

if __name__ == "__main__":
    main()