"""Label a per-body reaction/force CSV by TM segment.

The reporter CSV identifies each rigid body only by `atom_idx` (a representative
atom). This maps atom_idx -> receptor residue (via the prmtop) -> segment label
(the same SEGMENTS the run_*.py scripts use), and writes `<csv>.labeled.csv` with
a leading `segment` column so you can tell which TM / microswitch / H8 body each
row belongs to.

Keep SEGMENTS below in sync with the run_<PDB>.*.py scripts if you edit them there.

Usage (from repo root):
    python3 examples/febs/label_reactions.py 3SN6.lig.0.reactions.csv examples/febs/3SN6.lig.nanodisc.prmtop
"""

import sys

import parmed as pmd

# MODEL/prmtop numbering; must match the run_<PDB>.*.py SEGMENTS.
SEGMENTS = {
    "3SN6": [  # human beta2AR
        ("TM1", 1, 41), ("TM2_EC", 42, 49), ("Na_D2.50", 50, 51), ("TM2_IC", 52, 74),
        ("TM3", 75, 100), ("DRY_ionic", 101, 120), ("TM4", 121, 168), ("TM5", 169, 237),
        ("TM6_EC", 238, 255), ("CWxP_toggle", 256, 259), ("TM6_cyto", 260, 275),
        ("TM7", 276, 292), ("NPxxY", 293, 297), ("H8", 298, 312),
    ],
    "7JJO": [  # turkey beta1AR
        ("TM1", 1, 39), ("TM2_EC", 40, 47), ("Na_D2.50", 48, 49), ("TM2_IC", 50, 72),
        ("TM3", 73, 98), ("DRY_ionic", 99, 120), ("TM4", 121, 165), ("TM5", 166, 250),
        ("TM6_EC", 251, 262), ("CWxP_toggle", 263, 266), ("TM6_cyto", 267, 284),
        ("TM7", 285, 299), ("NPxxY", 300, 304), ("H8", 305, 318),
    ],
    "7DH5": [  # dog beta3AR
        ("TM1", 1, 38), ("TM2_EC", 39, 47), ("Na_D2.50", 48, 49), ("TM2_IC", 50, 73),
        ("TM3", 74, 98), ("DRY_ionic", 99, 118), ("TM4", 119, 167), ("TM5", 168, 259),
        ("TM6_EC", 260, 268), ("CWxP_toggle", 269, 272), ("TM6_cyto", 273, 286),
        ("TM7", 287, 306), ("NPxxY", 307, 311), ("H8", 312, 325),
    ],
}


def receptor_atom_res(prmtop_path):
    """atom_idx (global prmtop index) -> residue number, for the molecule-0 receptor."""
    parm = pmd.load_file(prmtop_path)
    out = {}
    started = False
    for res in parm.residues:
        is_aa = len(res.name.strip()) == 3 and res.name.strip() not in ("ACE", "NME", "NHE")
        if is_aa:
            started = True
            for atom in res.atoms:
                out[atom.idx] = res.number
        elif started:
            break  # receptor chain ended
    return out


def seg_label(resid, segments):
    for label, lo, hi in segments:
        if lo <= resid <= hi:
            return label
    return "?"


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: label_reactions.py <reactions.csv> <prmtop>")
    csv_path, prmtop_path = sys.argv[1], sys.argv[2]

    tag = next((t for t in SEGMENTS if t in prmtop_path or t in csv_path), None)
    if tag is None:
        sys.exit("could not detect receptor (3SN6 / 7JJO / 7DH5) from the file names")
    segments = SEGMENTS[tag]
    atom_res = receptor_atom_res(prmtop_path)

    out_path = csv_path.rsplit(".", 1)[0] + ".labeled.csv"
    n = 0
    with open(csv_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.lstrip().startswith("#"):
                fout.write("# segment," + line.lstrip()[1:].lstrip() + "\n")
                continue
            fields = line.split(",")
            resid = atom_res.get(int(fields[3]))  # atom_idx is column 4 (0-based 3)
            label = seg_label(resid, segments) if resid is not None else "?"
            fout.write(label + "," + line + "\n")
            n += 1
    print(f"[{tag}] wrote {out_path} ({n} rows labeled)")


if __name__ == "__main__":
    main()
