"""
betaAR_pair_common.py  -  within-receptor agonist-bound vs apo reaction-force suite
====================================================================================
Shared engine behind analyse_3SN6.py / analyse_7JJO.py / analyse_7DH5.py.

Recreates the eight-figure analyse_ffar1.py poster suite (numbering preserved) as a
bound-vs-apo comparison for ONE β-adrenergic receptor, from the per-body spatial
reaction-force CSV (`<PDB>.<state>.<replica>.reactions.csv`) + paired DCD.

Spec: docs/specs/gpcr-world-design/30-ligand-vs-apo-reaction-analysis.md

Correctness constraints enforced here (see spec §3):
  * Receptor sits on a FREE root -> it tumbles in Ground. Only rotation-invariant
    scalars (|F|, |τ|, arm, cos∠(F,τ)) and per-frame-membrane-normal projections
    (F_axial/F_lat/T_spin/T_rock) are compared. Raw Ground components never are.
  * Frames are accepted MC/REX ROUNDS, not physical time -> distributional
    comparison + effect size; no causal-lag claim (Fig 5 headline is lag-0).
  * A frame where ANY body has |F| > 1e5 (native units) is dropped whole.
  * reactions.csv has no u/uDot -> u is reconstructed as the flexed boundary
    dihedral (φ, ψ at proline) of each segment's first residue, from the DCD.

Units are RAW native reporter units: force kJ·mol⁻¹·nm⁻¹, torque kJ·mol⁻¹.
Never cross-compare absolute numbers with the analyse_ffar1.py (kcal/Å) pipeline.
"""

import argparse
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal, stats

# ----------------------------------------------------------------------------
# Atom→segment table (inlined so this module has NO cross-file import and can be
# copied/relocated freely). MUST stay in sync with examples/febs/label_reactions.py
# and the run_<PDB>.*.py SEGMENTS (model/prmtop residue numbering).
# ----------------------------------------------------------------------------
SEGMENTS = {
    "3SN6": [  # human β2AR
        ("TM1", 1, 41), ("TM2_EC", 42, 49), ("Na_D2.50", 50, 51), ("TM2_IC", 52, 74),
        ("TM3", 75, 100), ("DRY_ionic", 101, 120), ("TM4", 121, 168), ("TM5", 169, 237),
        ("TM6_EC", 238, 255), ("CWxP_toggle", 256, 259), ("TM6_cyto", 260, 275),
        ("TM7", 276, 292), ("NPxxY", 293, 297), ("H8", 298, 312),
    ],
    "7JJO": [  # turkey β1AR
        ("TM1", 1, 39), ("TM2_EC", 40, 47), ("Na_D2.50", 48, 49), ("TM2_IC", 50, 72),
        ("TM3", 73, 98), ("DRY_ionic", 99, 120), ("TM4", 121, 165), ("TM5", 166, 250),
        ("TM6_EC", 251, 262), ("CWxP_toggle", 263, 266), ("TM6_cyto", 267, 284),
        ("TM7", 285, 299), ("NPxxY", 300, 304), ("H8", 305, 318),
    ],
    "7DH5": [  # dog β3AR
        ("TM1", 1, 38), ("TM2_EC", 39, 47), ("Na_D2.50", 48, 49), ("TM2_IC", 50, 73),
        ("TM3", 74, 98), ("DRY_ionic", 99, 118), ("TM4", 119, 167), ("TM5", 168, 259),
        ("TM6_EC", 260, 268), ("CWxP_toggle", 269, 272), ("TM6_cyto", 273, 286),
        ("TM7", 287, 306), ("NPxxY", 307, 311), ("H8", 312, 325),
    ],
}


def receptor_atom_res(prmtop_path):
    """atom_idx (global prmtop index) -> residue number, for the molecule-0 receptor."""
    import parmed as pmd
    parm = pmd.load_file(prmtop_path)
    out, started = {}, False
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


# ----------------------------------------------------------------------------
# Data-file location (configurable; set by run() from --datadir / --febsdir).
# DATADIR holds the *.reactions.csv and *.dcd; FEBSDIR holds the *.nanodisc.prmtop.
# Defaults auto-detect so the scripts work whether they sit in the repo tree or a
# relocated download/ folder.
# ----------------------------------------------------------------------------
DATADIR = "."
FEBSDIR = None  # None → search candidates in _febs_candidates()


def _febs_candidates():
    here = os.path.dirname(os.path.abspath(__file__))
    return [
        FEBSDIR,
        os.path.join(DATADIR, "examples", "febs"),
        os.path.join(here, "..", "..", "examples", "febs"),
        os.path.join(here, "examples", "febs"),
        DATADIR,   # prmtops sitting next to the CSVs/DCDs
        here,      # prmtops sitting next to the scripts
    ]


def prmtop_path(receptor, state):
    fn = f"{receptor}.{state}.nanodisc.prmtop"
    for d in _febs_candidates():
        if d and os.path.exists(os.path.join(d, fn)):
            return os.path.join(d, fn)
    return os.path.join(FEBSDIR or "examples/febs", fn)  # non-existent → clear error


def data_path(receptor, state, replica, suffix):
    return os.path.join(DATADIR, f"{receptor}.{state}.{replica}.{suffix}")

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
STATES = [("lig", "agonist-bound"), ("noLig", "apo")]
STATE_COLOR = {"noLig": "#9e9e9e", "lig": "#d1495b"}  # apo grey, bound red
STATE_NAME = {"noLig": "apo", "lig": "bound"}
RECEPTOR_LABEL = {"3SN6": "β2-AR · 3SN6", "7DH5": "β3-AR · 7DH5", "7JJO": "β1-AR · 7JJO"}

EXPLODE_CUTOFF = 1.0e5  # |F| above this = accepted clashing MC move; excise whole frame

# Split flex sub-segments merged into the canonical seven TM helices + H8 for
# Cα-axis geometry (tilt / centroid / interhelical distance). Spec §4 GEO_BODIES.
TM_MERGE = {
    "TM1": ["TM1"], "TM2": ["TM2_EC", "Na_D2.50", "TM2_IC"], "TM3": ["TM3", "DRY_ionic"],
    "TM4": ["TM4"], "TM5": ["TM5"], "TM6": ["TM6_EC", "CWxP_toggle", "TM6_cyto"],
    "TM7": ["TM7", "NPxxY"], "H8": ["H8"],
}
TM7_ORDER = ["TM1", "TM2", "TM3", "TM4", "TM5", "TM6", "TM7", "H8"]
ACT_PAIRS = [("TM3", "TM6"), ("TM5", "TM6"), ("TM6", "TM7"),
             ("TM3", "TM7"), ("TM1", "TM7"), ("TM2", "TM4")]

# (col, pretty, unit)
METRICS = [
    ("Fmag", "|F|", "kJ·mol⁻¹·nm⁻¹"),
    ("Tmag", "|τ|", "kJ·mol⁻¹"),
    ("arm", "moment arm |τ|/|F|", "nm"),
    ("F_axial", "F axial (∥ n̂)", "kJ·mol⁻¹·nm⁻¹"),
    ("F_lat", "F lateral (⊥ n̂)", "kJ·mol⁻¹·nm⁻¹"),
    ("T_spin", "τ spin (about n̂)", "kJ·mol⁻¹"),
    ("T_rock", "τ rock (⊥ n̂)", "kJ·mol⁻¹"),
    ("cosFT", "cos∠(F,τ)", "–"),
]
RAW_COLS = ["frame", "replica", "body_idx", "atom_idx", "fx", "fy", "fz", "tx", "ty", "tz"]
N_BOOT = 400
RNG = np.random.default_rng(0)

plt.rcParams.update(
    {
        "font.size": 12, "axes.titlesize": 13, "axes.labelsize": 12,
        "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 9,
        "figure.dpi": 200, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5, "lines.linewidth": 1.4,
    }
)


# ----------------------------------------------------------------------------
# Robust statistics (spec §7) — self-contained, autocorrelation-aware
# ----------------------------------------------------------------------------
def cliffs_delta(a, b):
    """δ = P(a>b) − P(a<b). Positive ⇒ a stochastically larger than b."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return np.nan
    bs = np.sort(b)
    gt = np.searchsorted(bs, a, side="left").sum()   # #(b < a)
    ge = np.searchsorted(bs, a, side="right").sum()   # #(b <= a)
    lt = n * m - ge                                    # #(b > a)
    return (gt - lt) / (n * m)


def _block_resample(x, L):
    n = len(x)
    nb = int(np.ceil(n / L))
    starts = RNG.integers(0, n - L + 1, size=nb) if n > L else np.zeros(nb, int)
    idx = (starts[:, None] + np.arange(L)[None, :]).ravel()[:n]
    return x[idx]


def cliffs_ci(a, b, nboot=N_BOOT):
    """Moving-block-bootstrap 95% CI on Cliff's δ (partial MC-autocorr correction)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) < 5 or len(b) < 5:
        return (np.nan, np.nan)
    La = max(1, int(round(len(a) ** 0.5)))
    Lb = max(1, int(round(len(b) ** 0.5)))
    d = np.empty(nboot)
    for i in range(nboot):
        d[i] = cliffs_delta(_block_resample(a, La), _block_resample(b, Lb))
    return tuple(np.percentile(d, [2.5, 97.5]))


def mwu_p(a, b):
    try:
        return stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
    except ValueError:
        return np.nan


def _stars(p):
    if not np.isfinite(p):
        return ""
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else ""


def circ_mean(theta):
    """Circular mean of angles (rad)."""
    theta = np.asarray(theta, float)
    return float(np.arctan2(np.sin(theta).mean(), np.cos(theta).mean())) if len(theta) else np.nan


def circ_std(theta):
    """Circular standard deviation (rad): sqrt(-2 ln R)."""
    theta = np.asarray(theta, float)
    if len(theta) == 0:
        return np.nan
    R = np.hypot(np.cos(theta).mean(), np.sin(theta).mean())
    R = min(max(R, 1e-12), 1.0)
    return float(np.sqrt(-2.0 * np.log(R)))


def circ_diff(theta):
    """Consecutive circular differences wrapped to (−π, π]."""
    d = np.diff(np.asarray(theta, float))
    return (d + np.pi) % (2 * np.pi) - np.pi


# ----------------------------------------------------------------------------
# Geometry primitives
# ----------------------------------------------------------------------------
def _dihedral(p0, p1, p2, p3):
    """Signed dihedral (rad) for stacks of points, each (n_frames, 3)."""
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1n = b1 / np.linalg.norm(b1, axis=1, keepdims=True)
    v = b0 - (b0 * b1n).sum(1, keepdims=True) * b1n
    w = b2 - (b2 * b1n).sum(1, keepdims=True) * b1n
    x = (v * w).sum(1)
    y = (np.cross(b1n, v) * w).sum(1)
    return np.arctan2(y, x)


def _pca_axis(pos):
    """Long axis (unit vector) of a Cα cloud via top eigenvector of covariance."""
    c = pos - pos.mean(0)
    return np.linalg.eigh(np.cov(c.T))[1][:, -1]


# ----------------------------------------------------------------------------
# Per-receptor topology maps
# ----------------------------------------------------------------------------
def _seg_of_body():
    """segment label -> canonical TM (via TM_MERGE); '?' if unmapped."""
    out = {}
    for tm, subs in TM_MERGE.items():
        for s in subs:
            out[s] = tm
    return out


def _tm_ranges(receptor):
    """canonical TM -> (lo, hi) residue span (merged segments)."""
    seg = {lbl: (lo, hi) for lbl, lo, hi in SEGMENTS[receptor]}
    out = {}
    for tm, subs in TM_MERGE.items():
        los = [seg[s][0] for s in subs if s in seg]
        his = [seg[s][1] for s in subs if s in seg]
        if los:
            out[tm] = (min(los), max(his))
    return out


def _atom_index(prmtop_path):
    """Receptor (molecule-0) maps: res.number->{name:idx} and P31 indices."""
    import parmed as pmd
    parm = pmd.load_file(prmtop_path)
    name_idx, started = {}, False
    for res in parm.residues:
        is_aa = len(res.name.strip()) == 3 and res.name.strip() not in ("ACE", "NME", "NHE")
        if is_aa:
            started = True
            name_idx[res.number] = {a.name: a.idx for a in res.atoms}
        elif started:
            break  # receptor chain ended
    p31 = np.array([a.idx for a in parm.atoms if a.name == "P31"])
    return name_idx, p31


# ----------------------------------------------------------------------------
# Geometry from the paired DCD (per-frame normal, tilt, centroids, u)
# ----------------------------------------------------------------------------
def compute_geometry(receptor, state, replica, stride):
    """
    Returns a dict of frame-indexed DataFrames, or None if the DCD is absent:
        normal : frame, nx, ny, nz
        tilt   : frame, tm, tilt          (deg, vs per-frame n̂; GEO_BODIES)
        cent   : frame, tm, cx, cy, cz
        u      : frame, body, u, uDot      (rad; DIH_BODIES, reconstructed dihedral)
    Frame index matches the reaction CSV 1:1 (spec §3.4); stride subsamples both.
    """
    try:
        import mdtraj as md
    except ImportError:
        print("  MDTraj not found — geometry (Figs 2/4/5/6/7/8) skipped.")
        return None

    dcd = data_path(receptor, state, replica, "dcd")
    prm = prmtop_path(receptor, state)
    if not os.path.exists(dcd):
        print(f"  {receptor}.{state}: no DCD ({dcd}) — geometry skipped.")
        return None

    name_idx, p31 = _atom_index(prm)
    ca_of_res = {n: d["CA"] for n, d in name_idx.items() if "CA" in d}

    # merged-TM Cα sets (need ≥4 Cα)
    tmr = _tm_ranges(receptor)
    tm_ca = {tm: np.array([ca_of_res[n] for n in range(lo, hi + 1) if n in ca_of_res])
             for tm, (lo, hi) in tmr.items()}
    tm_ca = {tm: idx for tm, idx in tm_ca.items() if len(idx) >= 4}

    # boundary-dihedral atom picks per DIH body (segment[1:]; first residue = seg lo)
    segs = SEGMENTS[receptor]
    resname = None
    import parmed as pmd
    parm = pmd.load_file(prm)
    resname = {r.number: r.name.strip() for r in parm.residues}
    dih_atoms = {}  # body label -> (kind, (i0,i1,i2,i3))
    for label, lo, _hi in segs[1:]:
        is_pro = resname.get(lo) == "PRO"
        if is_pro:  # ψ = N(n),CA(n),C(n),N(n+1)
            need = [(lo, "N"), (lo, "CA"), (lo, "C"), (lo + 1, "N")]
            kind = "psi"
        else:       # φ = C(n-1),N(n),CA(n),C(n)
            need = [(lo - 1, "C"), (lo, "N"), (lo, "CA"), (lo, "C")]
            kind = "phi"
        try:
            idx = tuple(name_idx[rn][an] for rn, an in need)
            dih_atoms[label] = (kind, idx)
        except KeyError:
            pass  # boundary at a chain terminus / missing atom — no u for this body

    # stride=s keeps DCD frames 0, s, 2s, … ; their TRUE frame numbers are i*stride,
    # which is what the reaction CSV is keyed on (load_forces subsamples to the same
    # set), so the force↔geometry join in §3.4 stays 1:1 for any stride.
    t = md.load(dcd, top=md.load_topology(prm), stride=stride)
    xyz = t.xyz * 10.0  # nm → Å
    nframes = t.n_frames
    fnum = np.arange(nframes) * stride  # true frame numbers of the kept frames

    # per-frame membrane normal
    up = None
    norm_rows, tilt_rows, cent_rows = [], [], []
    for f in range(nframes):
        P = xyz[f, p31]
        if up is None:
            up = P[:, 2] > np.median(P[:, 2])
        n = P[up].mean(0) - P[~up].mean(0)
        n /= np.linalg.norm(n)
        norm_rows.append({"frame": int(fnum[f]), "nx": n[0], "ny": n[1], "nz": n[2]})
        for tm, idx in tm_ca.items():
            pos = xyz[f, idx]
            cen = pos.mean(0)
            axis = _pca_axis(pos)
            cos = abs(float(np.dot(axis, n)))
            tilt_rows.append({"frame": int(fnum[f]), "tm": tm,
                              "tilt": np.degrees(np.arccos(np.clip(cos, 0, 1)))})
            cent_rows.append({"frame": int(fnum[f]), "tm": tm,
                              "cx": cen[0], "cy": cen[1], "cz": cen[2]})

    # reconstructed u / uDot per DIH body (vectorized over frames)
    u_rows = []
    for label, (kind, idx) in dih_atoms.items():
        p0, p1, p2, p3 = (xyz[:, idx[k]] for k in range(4))
        u = _dihedral(p0, p1, p2, p3)              # (nframes,)
        udot = np.concatenate([[0.0], circ_diff(u)])  # per-sample circular change
        for f in range(nframes):
            u_rows.append({"frame": int(fnum[f]), "body": label,
                           "u": float(u[f]), "uDot": float(udot[f])})

    return {
        "normal": pd.DataFrame(norm_rows),
        "tilt": pd.DataFrame(tilt_rows),
        "cent": pd.DataFrame(cent_rows),
        "u": pd.DataFrame(u_rows),
        "n_frames": nframes,
    }


# ----------------------------------------------------------------------------
# Force CSV load + membrane-frame decomposition
# ----------------------------------------------------------------------------
def _label_map(receptor, state):
    prm = prmtop_path(receptor, state)
    atom_res = receptor_atom_res(prm)
    segs = SEGMENTS[receptor]
    return {a: seg_label(r, segs) for a, r in atom_res.items()}


def load_forces(receptor, state, replica, geo, stride):
    """Load reactions.csv, label bodies, excise exploded frames, decompose vs n̂.
    Returns (df, n_bad, n_tot). df has: frame, segment, Fmag, Tmag, arm, cosFT,
    F_axial, F_lat, T_spin, T_rock. `stride` subsamples the CSV to the SAME frames
    the DCD was strided to (frame % stride == 0), so force↔geometry stays aligned."""
    csv = data_path(receptor, state, replica, "reactions.csv")
    df = pd.read_csv(csv, comment="#", names=RAW_COLS)
    # tolerate an accidental non-comment header row
    df = df[pd.to_numeric(df["frame"], errors="coerce").notna()].copy()
    for c in RAW_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["frame", "atom_idx"])
    df["frame"] = df["frame"].astype(int)
    if stride > 1:
        df = df[df.frame % stride == 0].copy()  # match the strided DCD frame set

    amap = _label_map(receptor, state)
    df["segment"] = df.atom_idx.astype(int).map(amap).fillna("?")
    df = df[df.segment != "?"].copy()  # drop belt/ground/lipid bodies

    df["Fmag"] = np.sqrt(df.fx**2 + df.fy**2 + df.fz**2)
    df["Tmag"] = np.sqrt(df.tx**2 + df.ty**2 + df.tz**2)

    # whole-frame explode excision (spec §3.3)
    bad = df.loc[df.Fmag > EXPLODE_CUTOFF, "frame"].unique()
    n_bad, n_tot = len(bad), df.frame.nunique()
    df = df[~df.frame.isin(bad)].copy()

    # rotation-invariant scalars
    df["arm"] = df.Tmag / df.Fmag
    dotFT = df.fx * df.tx + df.fy * df.ty + df.fz * df.tz
    df["cosFT"] = dotFT / (df.Fmag * df.Tmag)

    # membrane-frame projection against the per-frame n̂ (fallback [0,0,1])
    if geo is not None and len(geo["normal"]):
        df = df.merge(geo["normal"], on="frame", how="left")
        miss = df[["nx", "ny", "nz"]].isna().any(axis=1).sum()
        if miss:
            df[["nx", "ny", "nz"]] = df[["nx", "ny", "nz"]].fillna({"nx": 0.0, "ny": 0.0, "nz": 1.0})
        if miss:
            print(f"  {receptor}.{state}: {miss} force rows had no matching DCD normal → z fallback.")
    else:
        df["nx"], df["ny"], df["nz"] = 0.0, 0.0, 1.0
        print(f"  {receptor}.{state}: no per-frame normal → membrane frame uses fixed n̂=[0,0,1].")

    fdotn = df.fx * df.nx + df.fy * df.ny + df.fz * df.nz
    df["F_axial"] = fdotn.abs()
    df["F_lat"] = np.sqrt((df.fx - fdotn * df.nx) ** 2 + (df.fy - fdotn * df.ny) ** 2
                          + (df.fz - fdotn * df.nz) ** 2)
    tdotn = df.tx * df.nx + df.ty * df.ny + df.tz * df.nz
    df["T_spin"] = tdotn.abs()
    df["T_rock"] = np.sqrt((df.tx - tdotn * df.nx) ** 2 + (df.ty - tdotn * df.ny) ** 2
                           + (df.tz - tdotn * df.nz) ** 2)

    df["state"] = state
    print(f"  {receptor:5s} {state:6s}: {n_tot - n_bad:4d} clean frames "
          f"({n_bad} exploded, {100 * n_bad / max(n_tot, 1):.1f}%)")
    return df, n_bad, n_tot


# ----------------------------------------------------------------------------
# Assemble both states for one receptor
# ----------------------------------------------------------------------------
def load_receptor(receptor, replica, stride):
    data, expl = {}, []
    for state, _ in STATES:
        geo = compute_geometry(receptor, state, replica, stride)
        force, nb, nt = load_forces(receptor, state, replica, geo, stride)
        data[state] = {"force": force, "geo": geo}
        expl.append({"state": state, "exploded": nb, "total": nt})
    return data, pd.DataFrame(expl)


# ----------------------------------------------------------------------------
# Statistics table (spec §7)
# ----------------------------------------------------------------------------
def compute_stats(data, bodies, do_ci=("Fmag", "Tmag", "arm")):
    rows = []
    for body in bodies:
        lig = data["lig"]["force"]
        apo = data["noLig"]["force"]
        for col, _, _ in METRICS:
            a = lig[lig.segment == body][col].dropna().values
            b = apo[apo.segment == body][col].dropna().values
            if len(a) < 5 or len(b) < 5:
                continue
            d = cliffs_delta(a, b)
            p = mwu_p(a, b)
            lo, hi = cliffs_ci(a, b) if col in do_ci else (np.nan, np.nan)
            rows.append(dict(body=body, metric=col,
                             median_lig=np.median(a), median_apo=np.median(b),
                             n_lig=len(a), n_apo=len(b),
                             cliffs_delta=d, mwu_p=p, ci_lo=lo, ci_hi=hi))
    # Explicit columns so an EMPTY result (e.g. <5 frames/state) still exposes
    # .body/.metric for downstream selection instead of an attributeless frame.
    cols = ["body", "metric", "median_lig", "median_apo", "n_lig", "n_apo",
            "cliffs_delta", "mwu_p", "ci_lo", "ci_hi"]
    return pd.DataFrame(rows, columns=cols)


def _med_sem(v):
    v = np.asarray(v, float)
    if len(v) == 0:
        return np.nan, 0.0
    return np.median(v), 1.253 * np.std(v) / np.sqrt(len(v))  # SEM of the median


# ============================================================================
# FIGURES  (numbering preserved from analyse_ffar1.py; each is bound vs apo)
# ============================================================================
def _bodyvals(data, state, body, col):
    d = data[state]["force"]
    return d[d.segment == body][col].dropna().values


def fig1_force_decomposition(receptor, data, bodies, out):
    """Fig 1 — median F_lat / F_axial / |τ| per body, apo vs bound."""
    show = [("F_lat", "Lateral |F| (⊥ n̂)"), ("F_axial", "Axial |F| (∥ n̂)"), ("Tmag", "Torque |τ|")]
    fig, axes = plt.subplots(len(show), 1, figsize=(max(11, 0.9 * len(bodies)), 11), sharex=True)
    x, w = np.arange(len(bodies)), 0.38
    for ax, (col, ylab) in zip(axes, show):
        for k, (state, _) in enumerate(STATES):
            med, sem = zip(*[_med_sem(_bodyvals(data, state, b, col)) for b in bodies])
            ax.bar(x + (k - 0.5) * w, med, w, yerr=sem, capsize=3,
                   color=STATE_COLOR[state], label=STATE_NAME[state], alpha=0.9)
        ax.set_ylabel(ylab)
        ax.legend(loc="upper right")
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(bodies, rotation=35, ha="right")
    fig.suptitle(f"Fig 1 — Per-body force/torque decomposition, apo vs bound  ·  {RECEPTOR_LABEL[receptor]}\n"
                 "median ± SEM over MC rounds; membrane frame vs per-frame n̂", fontsize=13, y=1.0)
    plt.tight_layout(); plt.savefig(out, bbox_inches="tight"); plt.close()
    print(f"Saved {out}")


def fig2_helix_tilt(receptor, data, out):
    """Fig 2 — per-helix tilt vs membrane normal, apo vs bound (violins)."""
    tilts = {s: (data[s]["geo"]["tilt"] if data[s]["geo"] else None) for s, _ in STATES}
    if any(t is None for t in tilts.values()):
        print("Fig 2 skipped (no geometry)."); return
    tms = [tm for tm in TM7_ORDER if tm in set(tilts["lig"].tm) & set(tilts["noLig"].tm)]
    if not tms:
        print("Fig 2 skipped (no shared helices)."); return
    fig, ax = plt.subplots(figsize=(max(10, 1.3 * len(tms)), 5.5))
    w = 0.34
    for k, (state, _) in enumerate(STATES):
        data_sets = [tilts[state][tilts[state].tm == tm].tilt.values for tm in tms]
        pos = np.arange(len(tms)) + (k - 0.5) * w
        parts = ax.violinplot(data_sets, positions=pos, widths=w, showmedians=True)
        for pc in parts["bodies"]:
            pc.set_facecolor(STATE_COLOR[state]); pc.set_alpha(0.6)
        for key in ("cbars", "cmins", "cmaxes", "cmedians"):
            if key in parts:
                parts[key].set_color(STATE_COLOR[state])
        ax.plot([], [], color=STATE_COLOR[state], label=STATE_NAME[state])
    ax.set_xticks(np.arange(len(tms))); ax.set_xticklabels(tms)
    ax.set_ylabel("Helix tilt vs membrane normal (°)")
    ax.legend(loc="upper right")
    ax.set_title(f"Fig 2 — Per-helix tilt distribution, apo vs bound  ·  {RECEPTOR_LABEL[receptor]}\n"
                 "(rotation-invariant; absolute azimuth omitted under Free-root tumbling)", fontsize=12)
    plt.tight_layout(); plt.savefig(out, bbox_inches="tight"); plt.close()
    print(f"Saved {out}")


def _heat_pair(fig, gs, row, data, col, bodies, title, cmap, diverging, index_key):
    """Two side-by-side (body × frame) heatmaps (apo | bound), shared color scale."""
    mats = {}
    for state, _ in STATES:
        if index_key == "segment":
            d = data[state]["force"]
            mat = d.pivot_table(index="segment", columns="frame", values=col, aggfunc="mean").reindex(bodies)
        else:  # u / uDot from geo
            g = data[state]["geo"]
            if g is None:
                mats[state] = None; continue
            d = g["u"]
            mat = d.pivot_table(index="body", columns="frame", values=col).reindex(bodies)
        mats[state] = mat
    if any(m is None or m.size == 0 for m in mats.values()):
        return None
    allvals = np.concatenate([m.values.ravel() for m in mats.values()])
    allvals = allvals[np.isfinite(allvals)]
    if diverging:
        lim = np.nanpercentile(np.abs(allvals), 98) if len(allvals) else 1.0
        vmin, vmax = -lim, lim
    else:
        vmin, vmax = 0.0, (np.nanpercentile(allvals, 98) if len(allvals) else 1.0)
    im = None
    for j, (state, _) in enumerate(STATES):
        ax = fig.add_subplot(gs[row, j])
        im = ax.imshow(mats[state].values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                       interpolation="nearest")
        ax.set_yticks(range(len(bodies)))
        ax.set_yticklabels(bodies if j == 0 else [], fontsize=7)
        ax.set_title(f"{title} — {STATE_NAME[state]}", fontsize=9)
        ax.set_xlabel("MC round")
    fig.colorbar(im, ax=fig.axes[-2:], fraction=0.02, pad=0.02)
    return True


def fig3_heatmaps(receptor, data, bodies, dih_bodies, out):
    """Fig 3 — (body × frame) heatmaps, apo | bound side by side, per channel."""
    chans = [("F_lat", "|F_lat|", "YlOrRd", False, "segment"),
             ("F_axial", "|F_axial|", "YlOrRd", False, "segment"),
             ("Tmag", "|τ|", "YlOrRd", False, "segment"),
             ("u", "u (rad)", "RdBu_r", True, "u"),
             ("uDot", "u̇ (rad/round)", "RdBu_r", True, "u")]
    fig = plt.figure(figsize=(14, 3.0 * len(chans)))
    gs = fig.add_gridspec(len(chans), 2, hspace=0.5, wspace=0.06)
    drew = False
    for row, (col, title, cmap, div, key) in enumerate(chans):
        bset = dih_bodies if key == "u" else bodies
        ok = _heat_pair(fig, gs, row, data, col, bset, title, cmap, div, key)
        drew = drew or bool(ok)
    if not drew:
        plt.close(); print("Fig 3 skipped (no data)."); return
    fig.suptitle(f"Fig 3 — Mechanical & coordinate heatmaps, apo | bound  ·  {RECEPTOR_LABEL[receptor]}\n"
                 "shared color scale per channel; frames are MC rounds, not time", fontsize=13, y=1.0)
    plt.savefig(out, bbox_inches="tight"); plt.close()
    print(f"Saved {out}")


def fig4_u_distributions(receptor, data, dih_bodies, out):
    """Fig 4 — u (circular) and |u̇| distributions per body, apo vs bound."""
    if any(data[s]["geo"] is None for s, _ in STATES):
        print("Fig 4 skipped (no geometry)."); return
    bodies = [b for b in dih_bodies
              if all(b in set(data[s]["geo"]["u"].body) for s, _ in STATES)]
    if not bodies:
        print("Fig 4 skipped (no u bodies)."); return
    fig, axes = plt.subplots(2, 1, figsize=(max(11, 0.9 * len(bodies)), 9), sharex=True)
    x, w = np.arange(len(bodies)), 0.34
    # (a) circular mean of u with circular-std whiskers
    for k, (state, _) in enumerate(STATES):
        g = data[state]["geo"]["u"]
        cm = [np.degrees(circ_mean(g[g.body == b].u.values)) for b in bodies]
        cs = [np.degrees(circ_std(g[g.body == b].u.values)) for b in bodies]
        axes[0].errorbar(x + (k - 0.5) * w, cm, yerr=cs, fmt="o", capsize=3,
                         color=STATE_COLOR[state], label=STATE_NAME[state])
    axes[0].set_ylabel("u  (circular mean ± circ-std, °)")
    axes[0].legend(loc="upper right")
    # (b) median |u̇|
    for k, (state, _) in enumerate(STATES):
        g = data[state]["geo"]["u"]
        med = [np.median(np.abs(g[g.body == b].uDot.values)) if (g.body == b).any() else np.nan
               for b in bodies]
        axes[1].bar(x + (k - 0.5) * w, np.degrees(med), w, color=STATE_COLOR[state],
                    label=STATE_NAME[state], alpha=0.9)
    axes[1].set_ylabel("median |u̇|  (°/round)")
    axes[1].set_xticks(x); axes[1].set_xticklabels(bodies, rotation=35, ha="right")
    fig.suptitle(f"Fig 4 — Boundary coordinate u and |u̇| per body, apo vs bound  ·  {RECEPTOR_LABEL[receptor]}\n"
                 "u reconstructed as the flexed φ/ψ dihedral; u̇ is per-sample change", fontsize=13, y=1.0)
    plt.tight_layout(); plt.savefig(out, bbox_inches="tight"); plt.close()
    print(f"Saved {out}")


def fig5_xcorr(receptor, data, dih_bodies, out, max_lag=20):
    """Fig 5 — |F|↔u association per body; headline lag-0 r, correlogram exploratory."""
    if any(data[s]["geo"] is None for s, _ in STATES):
        print("Fig 5 skipped (no geometry)."); return
    bodies = [b for b in dih_bodies
              if all(b in set(data[s]["geo"]["u"].body) for s, _ in STATES)]
    if not bodies:
        print("Fig 5 skipped (no u bodies)."); return
    nrow = len(bodies)
    fig, axes = plt.subplots(nrow, 1, figsize=(11, max(2.0 * nrow, 3)), sharex=True, squeeze=False)
    axes = axes[:, 0]
    for ax, body in zip(axes, bodies):
        for state, _ in STATES:
            fdf = data[state]["force"]
            udf = data[state]["geo"]["u"]
            f = fdf[fdf.segment == body].groupby("frame").Fmag.mean()
            u = udf[udf.body == body].set_index("frame").u
            common = f.index.intersection(u.index)
            if len(common) < 5:
                continue
            fa = f.loc[common].values; ua = u.loc[common].values
            fz = fa - fa.mean(); uz = ua - ua.mean()
            n = len(fz)
            lags = signal.correlation_lags(n, n, mode="full")
            xc = signal.correlate(fz, uz, mode="full")
            denom = np.std(fz) * np.std(uz) * n
            xc = xc / denom if denom > 0 else xc
            m = np.abs(lags) <= max_lag
            ax.plot(lags[m], xc[m], color=STATE_COLOR[state], label=STATE_NAME[state])
            r0 = xc[lags == 0]
            ax.plot(0, r0[0] if len(r0) else 0, "o", color=STATE_COLOR[state], ms=5)
        ax.axhline(0, color="k", lw=0.5)
        ax.axvline(0, color="gray", lw=0.6, ls="--")
        ax.set_ylabel(body, rotation=0, ha="right", labelpad=40, fontsize=9)
        ax.set_ylim(-1.05, 1.05)
    axes[0].legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("lag (MC rounds — sampling sequence, NOT physical time)")
    fig.suptitle(f"Fig 5 — |F| ↔ u association per body, apo vs bound  ·  {RECEPTOR_LABEL[receptor]}\n"
                 "● = lag-0 correlation (the interpretable readout); curve is exploratory", fontsize=12, y=1.0)
    plt.tight_layout(); plt.savefig(out, bbox_inches="tight"); plt.close()
    print(f"Saved {out}")


def fig6_force_vs_tilt(receptor, data, out):
    """Fig 6 — merged-TM |F| vs tilt scatter per helix, apo vs bound, with r + slope."""
    if any(data[s]["geo"] is None for s, _ in STATES):
        print("Fig 6 skipped (no geometry)."); return
    seg2tm = _seg_of_body()
    tms = [tm for tm in TM7_ORDER
           if all(tm in set(data[s]["geo"]["tilt"].tm) for s, _ in STATES)]
    if not tms:
        print("Fig 6 skipped (no shared helices)."); return
    ncol = 4; nrow = int(np.ceil(len(tms) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.5 * ncol, 4 * nrow), squeeze=False)
    axes = axes.flatten()
    rows = []  # dumped to <PDB>_fig6_stats.csv so r/slope are quotable, not PNG-only
    for i, tm in enumerate(tms):
        ax = axes[i]
        for state, _ in STATES:
            fdf = data[state]["force"].copy()
            fdf["tm"] = fdf.segment.map(seg2tm)
            f = fdf[fdf.tm == tm].groupby("frame").Fmag.mean()
            til = data[state]["geo"]["tilt"]
            t = til[til.tm == tm].set_index("frame").tilt
            common = f.index.intersection(t.index)
            if len(common) < 5:
                continue
            xf, yt = f.loc[common].values, t.loc[common].values
            r, p = stats.pearsonr(xf, yt)
            slope = float(np.polyfit(xf, yt, 1)[0])  # °tilt per unit |F| (compliance)
            rows.append(dict(tm=tm, state=STATE_NAME[state], r=r, p=p, slope=slope,
                             n=len(common), F_median=float(np.median(xf)),
                             tilt_median=float(np.median(yt))))
            ax.scatter(xf, yt, s=10, alpha=0.35, color=STATE_COLOR[state],
                       label=f"{STATE_NAME[state]}  r={r:+.2f}{_stars(p)}")
            xx = np.linspace(xf.min(), xf.max(), 50)
            ax.plot(xx, np.polyval(np.polyfit(xf, yt, 1), xx), color=STATE_COLOR[state], lw=1.2, ls="--")
        ax.set_xlabel("|F| (merged TM, per round)"); ax.set_ylabel("tilt (°)")
        ax.set_title(tm); ax.legend(fontsize=7)
    for j in range(len(tms), len(axes)):
        axes[j].set_visible(False)
    fig.suptitle(f"Fig 6 — Mean |F| vs helix tilt, apo vs bound  ·  {RECEPTOR_LABEL[receptor]}\n"
                 "per-frame within each state (equal-round; mechanical–geometric coupling)", fontsize=13, y=1.0)
    plt.tight_layout(); plt.savefig(out, bbox_inches="tight"); plt.close()
    print(f"Saved {out}")
    if rows:
        csv = out.rsplit(".", 1)[0] + "_stats.csv"
        pd.DataFrame(rows).to_csv(csv, index=False)
        print(f"Saved {csv}")


def fig7_interhelical(receptor, data, out):
    """Fig 7 — interhelical Cα-centroid distances per activation pair, apo vs bound."""
    if any(data[s]["geo"] is None for s, _ in STATES):
        print("Fig 7 skipped (no geometry)."); return
    dist = {s: {} for s, _ in STATES}
    for state, _ in STATES:
        cent = data[state]["geo"]["cent"]
        piv = cent.pivot(index="frame", columns="tm", values=["cx", "cy", "cz"])
        for a, b in ACT_PAIRS:
            try:
                c1 = piv.loc[:, (["cx", "cy", "cz"], a)].values
                c2 = piv.loc[:, (["cx", "cy", "cz"], b)].values
            except KeyError:
                continue
            dist[state][(a, b)] = np.linalg.norm(c1 - c2, axis=1)
    pairs = [p for p in ACT_PAIRS if all(p in dist[s] for s, _ in STATES)]
    if not pairs:
        print("Fig 7 skipped (no shared pairs)."); return
    fig, ax = plt.subplots(figsize=(max(10, 1.5 * len(pairs)), 6))
    w = 0.34
    for k, (state, _) in enumerate(STATES):
        sets = [dist[state][p] for p in pairs]
        pos = np.arange(len(pairs)) + (k - 0.5) * w
        parts = ax.violinplot(sets, positions=pos, widths=w, showmedians=True)
        for pc in parts["bodies"]:
            pc.set_facecolor(STATE_COLOR[state]); pc.set_alpha(0.6)
        for key in ("cbars", "cmins", "cmaxes", "cmedians"):
            if key in parts:
                parts[key].set_color(STATE_COLOR[state])
        ax.plot([], [], color=STATE_COLOR[state], label=STATE_NAME[state])
    for i, p in enumerate(pairs):
        dmed = np.median(dist["lig"][p]) - np.median(dist["noLig"][p])
        d = cliffs_delta(dist["lig"][p], dist["noLig"][p])
        ax.text(i, ax.get_ylim()[1], f"Δμ={dmed:+.1f}Å\nδ={d:+.2f}", ha="center", va="top", fontsize=8)
    ax.set_xticks(range(len(pairs)))
    ax.set_xticklabels([f"{a}-{b}" for a, b in pairs], rotation=25, ha="right")
    ax.set_ylabel("Cα-centroid distance (Å)")
    ax.legend(loc="upper right")
    ax.set_title(f"Fig 7 — Interhelical activation distances, apo vs bound  ·  {RECEPTOR_LABEL[receptor]}\n"
                 "Δμ = median(bound) − median(apo); δ = Cliff's effect size", fontsize=12)
    plt.tight_layout(); plt.savefig(out, bbox_inches="tight"); plt.close()
    print(f"Saved {out}")


def fig8_u_dispersion(receptor, data, dih_bodies, out):
    """Fig 8 — circular dispersion of u per body, apo vs bound (sampling spread)."""
    if any(data[s]["geo"] is None for s, _ in STATES):
        print("Fig 8 skipped (no geometry)."); return
    bodies = [b for b in dih_bodies
              if all(b in set(data[s]["geo"]["u"].body) for s, _ in STATES)]
    if not bodies:
        print("Fig 8 skipped (no u bodies)."); return
    fig, ax = plt.subplots(figsize=(max(11, 0.9 * len(bodies)), 5))
    x, w = np.arange(len(bodies)), 0.38
    for k, (state, _) in enumerate(STATES):
        g = data[state]["geo"]["u"]
        cs = [np.degrees(circ_std(g[g.body == b].u.values)) for b in bodies]
        ax.bar(x + (k - 0.5) * w, cs, w, color=STATE_COLOR[state], label=STATE_NAME[state], alpha=0.9)
    ax.set_xticks(x); ax.set_xticklabels(bodies, rotation=35, ha="right")
    ax.set_ylabel("circular std of u  (°)")
    ax.legend(loc="upper right")
    ax.set_title(f"Fig 8 — Conformational spread of the flexed coordinate u, apo vs bound  ·  "
                 f"{RECEPTOR_LABEL[receptor]}\n(sampler-appropriate recast of internal-coordinate RMSD)", fontsize=12)
    plt.tight_layout(); plt.savefig(out, bbox_inches="tight"); plt.close()
    print(f"Saved {out}")


# ----------------------------------------------------------------------------
# Interpretations
# ----------------------------------------------------------------------------
def write_interpretations(receptor, stat, expl, bodies, path):
    L = [f"{RECEPTOR_LABEL[receptor]} — reaction-force analysis, agonist-bound vs apo",
         "=" * 68,
         "Units: |F| kJ·mol⁻¹·nm⁻¹, |τ| kJ·mol⁻¹, arm nm.  δ = Cliff's (bound−apo); δ>0 ⇒ rises on binding.",
         "Frames are accepted MC/REX rounds (NOT physical time); comparison is distributional.",
         "p-values are UNCORRECTED for MC autocorrelation (δ CIs use a moving-block bootstrap).",
         "Membrane frame projected on the per-frame phosphate normal (fixed [0,0,1] fallback if no DCD).",
         "",
         "Exploded rounds (|F|>1e5) excised (whole frame):"]
    for _, r in expl.iterrows():
        L.append(f"  {STATE_NAME[r.state]:5s}: {r.exploded}/{r.total} "
                 f"({100 * r.exploded / max(r.total, 1):.1f}%)")
    L.append("")
    L.append("Per-body effect size δ (bound−apo), * p<0.05 ** p<0.01 *** p<0.001:")
    for col, name, _ in METRICS:
        L.append(f"\n[{name}]")
        for body in bodies:
            row = stat[(stat.body == body) & (stat.metric == col)]
            if row.empty:
                L.append(f"  {body:12s} NA"); continue
            r = row.iloc[0]
            ci = (f"  CI[{r.ci_lo:+.2f},{r.ci_hi:+.2f}]"
                  if np.isfinite(r.ci_lo) else "")
            L.append(f"  {body:12s} δ={r.cliffs_delta:+.2f}{_stars(r.mwu_p)}"
                     f"  (med bound={r.median_lig:.2f}, apo={r.median_apo:.2f}){ci}")
    with open(path, "w") as f:
        f.write("\n".join(L))
    print(f"Saved {path}")


# ----------------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------------
def run(receptor):
    if receptor not in SEGMENTS:
        sys.exit(f"unknown receptor {receptor!r}; expected one of {sorted(SEGMENTS)}")
    ap = argparse.ArgumentParser(description=f"Bound-vs-apo reaction analysis: {receptor}")
    ap.add_argument("--replica", type=int, default=0, help="replica index in the CSV/DCD filenames")
    ap.add_argument("--stride", type=int, default=1,
                    help="subsample stride applied to BOTH the DCD and the reactions.csv "
                         "(keeps frames 0, stride, 2·stride, …); keeps force↔geometry aligned")
    ap.add_argument("--datadir", default=".",
                    help="directory holding the *.reactions.csv and *.dcd (default: cwd)")
    ap.add_argument("--febsdir", default=None,
                    help="directory holding the *.nanodisc.prmtop (default: auto-detect)")
    ap.add_argument("--outdir", default=".", help="directory for figures/tables")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    global DATADIR, FEBSDIR
    DATADIR = args.datadir
    FEBSDIR = args.febsdir

    def op(name):
        return os.path.join(args.outdir, f"{receptor}_{name}")

    print(f"\n=== {RECEPTOR_LABEL[receptor]}  (replica {args.replica}, stride {args.stride}) ===")
    print("Loading reaction CSVs + DCD geometry ...")
    data, expl = load_receptor(receptor, args.replica, args.stride)

    label_order = [lbl for lbl, _, _ in SEGMENTS[receptor]]
    common = set(data["lig"]["force"].segment.unique()) & set(data["noLig"]["force"].segment.unique())
    force_bodies = [b for b in label_order if b in common]
    dih_bodies = [b for b in label_order[1:] if b in common]  # SEGMENTS[1:], reported
    print(f"Force bodies (both states): {force_bodies}")

    print("Computing robust statistics (bootstrap CIs) ...")
    stat = compute_stats(data, force_bodies)
    stat.to_csv(op("stats.csv"), index=False)
    print(f"Saved {op('stats.csv')}")
    if stat.empty:
        nfr = {s: data[s]["force"].frame.nunique() for s, _ in STATES}
        print(f"  WARNING: no statistics computed — need ≥5 clean frames per state, "
              f"have {nfr}. Lower --stride (currently subsampling too aggressively). "
              f"Figures still render but are not statistically meaningful.")

    fig1_force_decomposition(receptor, data, force_bodies, op("fig1.png"))
    fig2_helix_tilt(receptor, data, op("fig2.png"))
    fig3_heatmaps(receptor, data, force_bodies, dih_bodies, op("fig3.png"))
    fig4_u_distributions(receptor, data, dih_bodies, op("fig4.png"))
    fig5_xcorr(receptor, data, dih_bodies, op("fig5.png"))
    fig6_force_vs_tilt(receptor, data, op("fig6.png"))
    fig7_interhelical(receptor, data, op("fig7.png"))
    fig8_u_dispersion(receptor, data, dih_bodies, op("fig8.png"))

    write_interpretations(receptor, stat, expl, force_bodies, op("interpretations.txt"))
    print(f"\nDone. {receptor}_fig1..fig8 + {receptor}_stats.csv + {receptor}_interpretations.txt")
