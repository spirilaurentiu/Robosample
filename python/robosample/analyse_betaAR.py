"""
analyse_betaAR.py  -  β-adrenergic GPCR reaction-force poster pipeline
======================================================================
Ligand (agonist-bound) vs apo comparison of per-body spatial reaction
forces (F, tau about the body origin, in Ground) for three β-ARs:

    3SN6  β2-AR (human)   7DH5  β3-AR (dog)   7JJO  β1-AR (turkey)

Design (deliberately fully matched — see docs/reaction_force_analysis_report.md):
  * Only the bodies REPORTED IN ALL SIX files are analysed, so every cell of the
    3 receptor × 2 state × N body grid is populated. That common set is
        TM7, TM6_cyto, NPxxY, H8
    — TM helices AND the NPxxY microswitch AND helix-8, not microswitches only.
    (TM2_EC / TM3 / TM6_EC are dropped because 7DH5 never reports them;
     Na_D2.50 is dropped because it is apo-only — the Na+ body collapses on binding.)
  * Exploded MC rounds (|F| > 1e5, accepted clashing moves) are excised.
  * Units are RAW native reporter units: force kJ·mol⁻¹·nm⁻¹, torque kJ·mol⁻¹.
    Never cross-compare absolute numbers with the ffar1 (kcal/Å) pipeline.

Metrics per body (all comparable bound-vs-apo despite the Free-root receptor
tumbling, because they are rotation-invariant or projected onto the membrane
normal, which the Tier-2 check found is ≈ [0,0,1] for every flat nanodisc):
  rotation-invariant : |F|, |τ|, moment arm |τ|/|F|, cos∠(F,τ)
  membrane-frame     : F_axial=|Fz|, F_lat=√(Fx²+Fy²)   (insertion vs in-plane shear)
                       τ_spin =|τz|, τ_rock=√(τx²+τy²)   (twist about normal vs rocking)

Stats: median per state; Cliff's δ (lig − apo, the load-bearing effect size);
Mann-Whitney U p; moving-block-bootstrap 95% CI on δ (partial MC-autocorrelation
correction). δ > 0 ⇒ load rises on agonist binding.

Run (repo root, robo_cpu env):
    python3 python/robosample/analyse_betaAR.py
Outputs: fig1..fig6 PNGs + betaAR_stats.csv + betaAR_interpretations.txt
"""

import os
import sys
import warnings

warnings.filterwarnings("ignore")

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy import stats

# label_reactions.py holds the authoritative atom→segment SEGMENTS table.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "febs"))
from label_reactions import SEGMENTS, receptor_atom_res, seg_label  # noqa: E402

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RECEPTORS = ["3SN6", "7DH5", "7JJO"]
STATES = [("lig", "agonist-bound"), ("noLig", "apo")]
RECEPTOR_LABEL = {"3SN6": "β2-AR · 3SN6", "7DH5": "β3-AR · 7DH5", "7JJO": "β1-AR · 7JJO"}

# Bodies reported in ALL six files (strict intersection). EC→IC-ish order.
COMMON_BODIES = ["TM7", "TM6_cyto", "NPxxY", "H8"]
# Extracellular pocket bodies — reported by 3SN6 & 7JJO only (7DH5 lacks them).
# Analysed as a SUPPLEMENTARY 2-receptor panel: the strongest binding signal lives here.
EC_BODIES = ["TM2_EC", "TM3", "TM6_EC"]
EC_RECEPTORS = ["3SN6", "7JJO"]
BODY_KIND = {"TM7": "TM helix", "TM6_cyto": "TM6 cytoplasmic", "NPxxY": "microswitch", "H8": "helix 8",
             "TM2_EC": "TM2 extracellular", "TM3": "TM helix", "TM6_EC": "TM6 extracellular"}

EXPLODE_CUTOFF = 1.0e5  # |F| above this = accepted clashing MC move; excise

# (col, pretty, unit, is_diverging_for_signed) — signed only for cosFT
METRICS = [
    ("Fmag", "|F|", "kJ·mol⁻¹·nm⁻¹"),
    ("Tmag", "|τ|", "kJ·mol⁻¹"),
    ("arm", "moment arm |τ|/|F|", "nm"),
    ("F_axial", "F axial (∥ normal)", "kJ·mol⁻¹·nm⁻¹"),
    ("F_lat", "F lateral (in-plane)", "kJ·mol⁻¹·nm⁻¹"),
    ("T_spin", "τ spin (about normal)", "kJ·mol⁻¹"),
    ("T_rock", "τ rock (in-plane)", "kJ·mol⁻¹"),
    ("cosFT", "cos∠(F,τ)", "–"),
]
METRIC_LABEL = {m[0]: m[1] for m in METRICS}
METRIC_UNIT = {m[0]: m[2] for m in METRICS}

STATE_COLOR = {"noLig": "#9e9e9e", "lig": "#d1495b"}  # apo grey, bound red
STATE_NAME = {"noLig": "apo", "lig": "bound"}
BODY_COLOR = {"TM7": "#377eb8", "TM6_cyto": "#4daf4a", "NPxxY": "#984ea3", "H8": "#ff7f00",
              "TM2_EC": "#e41a1c", "TM3": "#00b3b3", "TM6_EC": "#a65628"}
N_BOOT = 400
RNG = np.random.default_rng(0)

plt.rcParams.update(
    {
        "font.size": 12,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "figure.dpi": 200,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
        "lines.linewidth": 1.4,
    }
)


# ----------------------------------------------------------------------------
# Load + label + derive metrics
# ----------------------------------------------------------------------------
RAW_COLS = ["frame", "replica", "body_idx", "atom_idx", "fx", "fy", "fz", "tx", "ty", "tz"]


def _label_map(receptor, state):
    """atom_idx -> segment label, cached via the .labeled.csv if present."""
    prmtop = os.path.join(ROOT, "examples", "febs", f"{receptor}.{state}.nanodisc.prmtop")
    atom_res = receptor_atom_res(prmtop)
    segs = SEGMENTS[receptor]
    return {a: seg_label(r, segs) for a, r in atom_res.items()}


def load_one(receptor, state):
    csv = os.path.join(ROOT, f"{receptor}.{state}.0.reactions.csv")
    df = pd.read_csv(csv, comment="#", names=RAW_COLS)
    amap = _label_map(receptor, state)
    df["segment"] = df.atom_idx.map(amap).fillna("?")
    df = df[df.segment != "?"].copy()  # keep all receptor bodies; drop belt/ground

    df["Fmag"] = np.sqrt(df.fx**2 + df.fy**2 + df.fz**2)
    df["Tmag"] = np.sqrt(df.tx**2 + df.ty**2 + df.tz**2)

    # exploded-round bookkeeping: drop the WHOLE frame if ANY reported body blew up,
    # so every surviving frame is a clean, fully-populated snapshot (shared by the
    # common-body and extracellular analyses alike).
    bad_frames = df.loc[df.Fmag > EXPLODE_CUTOFF, "frame"].unique()
    n_bad = len(bad_frames)
    n_tot = df.frame.nunique()
    df = df[~df.frame.isin(bad_frames)].copy()

    # rotation-invariant
    df["arm"] = df.Tmag / df.Fmag
    dotFT = df.fx * df.tx + df.fy * df.ty + df.fz * df.tz
    df["cosFT"] = dotFT / (df.Fmag * df.Tmag)
    # membrane-frame (fixed n̂ = [0,0,1]; Tier-2: nanodiscs flat)
    df["F_axial"] = df.fz.abs()
    df["F_lat"] = np.sqrt(df.fx**2 + df.fy**2)
    df["T_spin"] = df.tz.abs()
    df["T_rock"] = np.sqrt(df.tx**2 + df.ty**2)

    df["receptor"] = receptor
    df["state"] = state
    print(f"  {receptor:5s} {state:6s}: {n_tot - n_bad:4d} clean frames "
          f"({n_bad} exploded, {100 * n_bad / max(n_tot, 1):.1f}%)")
    return df, n_bad, n_tot


def load_all():
    frames, expl = [], []
    for r in RECEPTORS:
        for s, _ in STATES:
            df, nb, nt = load_one(r, s)
            frames.append(df)
            expl.append({"receptor": r, "state": s, "exploded": nb, "total": nt})
    return pd.concat(frames, ignore_index=True), pd.DataFrame(expl)


# ----------------------------------------------------------------------------
# Robust statistics
# ----------------------------------------------------------------------------
def cliffs_delta(a, b):
    """δ = P(a>b) − P(a<b). Positive ⇒ a stochastically larger than b."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return np.nan
    bs = np.sort(b)
    gt = np.searchsorted(bs, a, side="left").sum()        # #(b < a)
    ge = np.searchsorted(bs, a, side="right").sum()        # #(b <= a)
    lt = n * m - ge                                        # #(b > a)
    return (gt - lt) / (n * m)


def _block_resample(x, L):
    n = len(x)
    nb = int(np.ceil(n / L))
    starts = RNG.integers(0, n - L + 1, size=nb) if n > L else np.zeros(nb, int)
    idx = (starts[:, None] + np.arange(L)[None, :]).ravel()[:n]
    return x[idx]


def cliffs_ci(a, b, nboot=N_BOOT):
    """Moving-block-bootstrap 95% CI on Cliff's δ (partial autocorr correction)."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) < 5 or len(b) < 5:
        return (np.nan, np.nan)
    La = max(1, int(round(len(a) ** 0.5)))
    Lb = max(1, int(round(len(b) ** 0.5)))
    d = np.empty(nboot)
    for i in range(nboot):
        d[i] = cliffs_delta(_block_resample(a, La), _block_resample(b, Lb))
    return tuple(np.percentile(d, [2.5, 97.5]))


def compute_stats(long, receptors=RECEPTORS, bodies=COMMON_BODIES, do_ci=("Tmag", "arm", "Fmag")):
    """Per (receptor, body, metric): medians, δ(lig−apo), MWU p, δ CI."""
    rows = []
    for r in receptors:
        for body in bodies:
            lig = long[(long.receptor == r) & (long.state == "lig") & (long.segment == body)]
            apo = long[(long.receptor == r) & (long.state == "noLig") & (long.segment == body)]
            for col, _, _ in METRICS:
                a = lig[col].dropna().values
                b = apo[col].dropna().values
                if len(a) < 5 or len(b) < 5:
                    continue
                delta = cliffs_delta(a, b)
                try:
                    p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
                except ValueError:
                    p = np.nan
                lo, hi = cliffs_ci(a, b) if col in do_ci else (np.nan, np.nan)
                rows.append(
                    dict(receptor=r, body=body, kind=BODY_KIND[body], metric=col,
                         median_lig=np.median(a), median_apo=np.median(b),
                         n_lig=len(a), n_apo=len(b),
                         cliffs_delta=delta, mwu_p=p, ci_lo=lo, ci_hi=hi)
                )
    return pd.DataFrame(rows)


def _piv(stat, r, value, bodies=COMMON_BODIES):
    """metric × body matrix of a stat column for one receptor."""
    sub = stat[stat.receptor == r]
    return sub.pivot(index="metric", columns="body", values=value).reindex(
        index=[m[0] for m in METRICS], columns=bodies
    )


def _draw_heatmap(ax, stat, r, bodies, vmax=0.65, first_col=True):
    """Cliff's-δ heatmap panel (metrics × bodies) for one receptor; returns the image."""
    dmat = _piv(stat, r, "cliffs_delta", bodies)
    pmat = _piv(stat, r, "mwu_p", bodies)
    im = ax.imshow(dmat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(bodies)))
    ax.set_xticklabels(bodies, rotation=30, ha="right")
    ax.set_yticks(range(len(METRICS)))
    ax.set_yticklabels([m[1] for m in METRICS] if first_col else [])
    ax.set_title(RECEPTOR_LABEL[r], fontsize=13)
    for i in range(dmat.shape[0]):
        for j in range(dmat.shape[1]):
            d = dmat.values[i, j]
            if not np.isfinite(d):
                continue
            ax.text(j, i, f"{d:+.2f}{_stars(pmat.values[i, j])}", ha="center", va="center",
                    fontsize=8, color="white" if abs(d) > 0.4 else "black")
    return im


def _stars(p):
    if not np.isfinite(p):
        return ""
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else ""


# ----------------------------------------------------------------------------
# Fig 1 — Fingerprint: median |F|, |τ|, moment arm, apo vs bound, per body
# ----------------------------------------------------------------------------
def fig1_fingerprint(long, out="fig1_fingerprint.png"):
    show = [("Fmag", "|F|  (kJ·mol⁻¹·nm⁻¹)"),
            ("Tmag", "|τ|  (kJ·mol⁻¹)"),
            ("arm", "moment arm |τ|/|F|  (nm)")]
    fig, axes = plt.subplots(len(show), len(RECEPTORS), figsize=(15, 10), sharex=True)
    x = np.arange(len(COMMON_BODIES))
    w = 0.38
    for ci, r in enumerate(RECEPTORS):
        for ri, (col, ylab) in enumerate(show):
            ax = axes[ri, ci]
            for k, (state, _) in enumerate(STATES):
                med, sem = [], []
                for body in COMMON_BODIES:
                    v = long[(long.receptor == r) & (long.state == state)
                             & (long.segment == body)][col].dropna().values
                    med.append(np.median(v) if len(v) else np.nan)
                    # SEM of the median ≈ 1.253·std/√n
                    sem.append(1.253 * np.std(v) / np.sqrt(len(v)) if len(v) else 0)
                ax.bar(x + (k - 0.5) * w, med, w, yerr=sem, capsize=3,
                       color=STATE_COLOR[state], label=STATE_NAME[state], alpha=0.9)
            if ri == 0:
                ax.set_title(RECEPTOR_LABEL[r], fontsize=13)
            if ci == 0:
                ax.set_ylabel(ylab, fontsize=11)
            ax.set_xticks(x)
            ax.set_xticklabels(COMMON_BODIES, rotation=30, ha="right")
            if ri == 0 and ci == len(RECEPTORS) - 1:
                ax.legend(loc="upper right", frameon=True)
    fig.suptitle("Fig 1 — Per-body reaction-force fingerprint: agonist-bound vs apo\n"
                 "(median ± SEM over MC rounds; bodies common to all 6 runs)",
                 fontsize=15, y=1.0)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


# ----------------------------------------------------------------------------
# Fig 2 — Effect-size redistribution map: Cliff's δ (bound − apo)
# ----------------------------------------------------------------------------
def fig2_heatmap(stat, out="fig2_heatmap.png"):
    fig, axes = plt.subplots(1, len(RECEPTORS), figsize=(16, 6), sharey=True)
    for ci, r in enumerate(RECEPTORS):
        im = _draw_heatmap(axes[ci], stat, r, COMMON_BODIES, first_col=(ci == 0))
    cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label("Cliff's δ  (bound − apo);  + = load rises on binding")
    fig.suptitle("Fig 2 — Force/torque redistribution map: effect size of agonist binding\n"
                 "(rows = mechanical channel, cols = body; * p<0.05  ** p<0.01  *** p<0.001, "
                 "uncorrected for MC autocorrelation)", fontsize=14, y=1.02)
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


# ----------------------------------------------------------------------------
# Fig 3 — Force↔torque conversion (moment-arm dumbbell) + load-share
# ----------------------------------------------------------------------------
def fig3_microswitch(long, out="fig3_microswitch.png"):
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), gridspec_kw={"width_ratios": [1.15, 1]})

    # (a) moment-arm dumbbell: apo→bound per body per receptor
    ax = axes[0]
    yl, yp = [], []
    y = 0.0
    for r in RECEPTORS:
        for body in COMMON_BODIES:
            apo = long[(long.receptor == r) & (long.state == "noLig") & (long.segment == body)]["arm"].median()
            lig = long[(long.receptor == r) & (long.state == "lig") & (long.segment == body)]["arm"].median()
            ax.plot([apo, lig], [y, y], color="#bbbbbb", lw=2, zorder=1)
            ax.scatter(apo, y, color=STATE_COLOR["noLig"], s=55, zorder=2)
            ax.scatter(lig, y, color=STATE_COLOR["lig"], s=55, zorder=3)
            ax.annotate("", xy=(lig, y), xytext=(apo, y),
                        arrowprops=dict(arrowstyle="->", color="#555", lw=1.2), zorder=2)
            yl.append(f"{r} · {body}")
            yp.append(y)
            y += 1
        y += 0.6
    ax.set_yticks(yp)
    ax.set_yticklabels(yl, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("median moment arm |τ|/|F|  (nm)")
    ax.set_title("Fig 3a — Force→torque conversion\n(apo ● → bound ●; right = more rotational load)")
    ax.scatter([], [], color=STATE_COLOR["noLig"], label="apo")
    ax.scatter([], [], color=STATE_COLOR["lig"], label="bound")
    ax.legend(loc="lower right")

    # (b) share of total intracellular |τ| carried by each body, apo vs bound
    ax = axes[1]
    xr = np.arange(len(RECEPTORS))
    w = 0.38
    for k, (state, _) in enumerate(STATES):
        bottoms = np.zeros(len(RECEPTORS))
        shares = {}
        for r in RECEPTORS:
            tot = sum(long[(long.receptor == r) & (long.state == state) & (long.segment == b)]["Tmag"].median()
                      for b in COMMON_BODIES)
            shares[r] = [long[(long.receptor == r) & (long.state == state) & (long.segment == b)]["Tmag"].median() / tot
                         for b in COMMON_BODIES]
        for bi, body in enumerate(COMMON_BODIES):
            vals = np.array([shares[r][bi] for r in RECEPTORS])
            ax.bar(xr + (k - 0.5) * w, vals, w, bottom=bottoms,
                   color=BODY_COLOR[body], edgecolor="white",
                   label=body if k == 0 else None, alpha=0.9 if state == "lig" else 0.55)
            bottoms += vals
    ax.set_xticks(xr)
    ax.set_xticklabels([RECEPTOR_LABEL[r] for r in RECEPTORS], rotation=15, ha="right")
    ax.set_ylabel("share of Σ median |τ|  (left bar = apo, right = bound)")
    ax.set_title("Fig 3b — Torque load-share redistribution\n(does binding move load between bodies?)")
    ax.legend(loc="upper center", ncol=4, fontsize=8, bbox_to_anchor=(0.5, -0.12))
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


# ----------------------------------------------------------------------------
# Fig 4 — Membrane-frame decomposition: axial/lateral force, spin/rock torque
# ----------------------------------------------------------------------------
def fig4_decomp(long, out="fig4_decomp.png"):
    fig, axes = plt.subplots(2, len(RECEPTORS), figsize=(15, 9), sharex=True)
    x = np.arange(len(COMMON_BODIES))
    w = 0.38
    rows = [(("F_axial", "F_lat"), "Force  (kJ·mol⁻¹·nm⁻¹)", "axial ∥n̂", "lateral ⊥n̂"),
            (("T_spin", "T_rock"), "Torque  (kJ·mol⁻¹)", "spin about n̂", "rock ⊥n̂")]
    for ci, r in enumerate(RECEPTORS):
        for ri, ((ca, cb), ylab, la, lb) in enumerate(rows):
            ax = axes[ri, ci]
            for k, (state, _) in enumerate(STATES):
                ma = [long[(long.receptor == r) & (long.state == state) & (long.segment == b)][ca].median()
                      for b in COMMON_BODIES]
                mb = [long[(long.receptor == r) & (long.state == state) & (long.segment == b)][cb].median()
                      for b in COMMON_BODIES]
                base = "//" if state == "lig" else None
                off = (k - 0.5) * w
                ax.bar(x + off, ma, w, color="#1f78b4", alpha=0.9 if state == "lig" else 0.5,
                       hatch=base, edgecolor="white")
                ax.bar(x + off, mb, w, bottom=ma, color="#e08214", alpha=0.9 if state == "lig" else 0.5,
                       hatch=base, edgecolor="white")
            if ri == 0:
                ax.set_title(RECEPTOR_LABEL[r], fontsize=13)
            if ci == 0:
                ax.set_ylabel(ylab)
            ax.set_xticks(x)
            ax.set_xticklabels(COMMON_BODIES, rotation=30, ha="right")
            if ci == len(RECEPTORS) - 1:
                from matplotlib.patches import Patch
                ax.legend(handles=[Patch(facecolor="#1f78b4", label=la),
                                   Patch(facecolor="#e08214", label=lb),
                                   Patch(facecolor="#888", alpha=0.5, label="apo (solid)"),
                                   Patch(facecolor="#888", hatch="//", label="bound (hatch)")],
                          fontsize=8, loc="upper right")
    fig.suptitle("Fig 4 — Membrane-frame decomposition (fixed n̂=[0,0,1]; flat nanodisc)\n"
                 "stacked bars: axial+lateral force / spin+rock torque, apo (solid) vs bound (hatched)",
                 fontsize=14, y=1.0)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


# ----------------------------------------------------------------------------
# Fig 5 — Consistency forest plot: δ(|τ|) and δ(arm) with block-boot CI
# ----------------------------------------------------------------------------
def fig5_forest(stat, out="fig5_forest.png"):
    metrics = [("Tmag", "δ  |τ|"), ("arm", "δ  moment arm")]
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharey=True)
    labels = [f"{r} · {b}" for r in RECEPTORS for b in COMMON_BODIES]
    ypos = np.arange(len(labels))[::-1]
    for ax, (col, title) in zip(axes, metrics):
        for y, r in zip(np.split(ypos, len(RECEPTORS)), RECEPTORS):
            for yi, body in zip(y, COMMON_BODIES):
                row = stat[(stat.receptor == r) & (stat.body == body) & (stat.metric == col)]
                if row.empty:
                    continue
                d = row.cliffs_delta.values[0]
                lo, hi = row.ci_lo.values[0], row.ci_hi.values[0]
                sig = np.isfinite(lo) and (lo > 0 or hi < 0)
                col_c = STATE_COLOR["lig"] if d > 0 else STATE_COLOR["noLig"]
                if np.isfinite(lo):
                    ax.plot([lo, hi], [yi, yi], color=col_c, lw=2, alpha=0.8)
                ax.scatter(d, yi, color=col_c, s=70 if sig else 40,
                           edgecolor="black" if sig else "none", zorder=3)
        ax.axvline(0, color="k", lw=1, ls="--")
        ax.set_xlim(-0.8, 0.8)
        ax.set_title(title)
        ax.set_xlabel("Cliff's δ (bound − apo)")
    axes[0].set_yticks(ypos)
    axes[0].set_yticklabels(labels, fontsize=9)
    fig.suptitle("Fig 5 — Effect-size consistency across receptors "
                 "(moving-block bootstrap 95% CI)\nfilled black-edged = CI excludes 0",
                 fontsize=14, y=1.0)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


# ----------------------------------------------------------------------------
# Fig 6 — Distribution ECDFs of |τ| + stability panel
# ----------------------------------------------------------------------------
def fig6_ecdf(long, expl, out="fig6_ecdf.png"):
    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(len(COMMON_BODIES), len(RECEPTORS) + 1,
                          width_ratios=[1, 1, 1, 0.8])
    for bi, body in enumerate(COMMON_BODIES):
        for ci, r in enumerate(RECEPTORS):
            ax = fig.add_subplot(gs[bi, ci])
            for state, _ in STATES:
                v = np.sort(long[(long.receptor == r) & (long.state == state)
                                 & (long.segment == body)]["Tmag"].dropna().values)
                if len(v):
                    ax.plot(v, np.linspace(0, 1, len(v)), color=STATE_COLOR[state],
                            label=STATE_NAME[state])
            if bi == 0:
                ax.set_title(RECEPTOR_LABEL[r], fontsize=12)
            if ci == 0:
                ax.set_ylabel(f"{body}\nECDF", fontsize=10)
            if bi == len(COMMON_BODIES) - 1:
                ax.set_xlabel("|τ|")
            if bi == 0 and ci == 0:
                ax.legend(fontsize=8, loc="lower right")
    # stability panel (right column, spanning)
    axs = fig.add_subplot(gs[:, -1])
    ys, labs, cols = [], [], []
    y = 0
    for r in RECEPTORS:
        for state, _ in STATES:
            row = expl[(expl.receptor == r) & (expl.state == state)].iloc[0]
            pct = 100 * row.exploded / max(row.total, 1)
            ys.append(pct)
            labs.append(f"{r}·{STATE_NAME[state]}")
            cols.append(STATE_COLOR[state])
            y += 1
    axs.barh(range(len(ys)), ys, color=cols)
    axs.set_yticks(range(len(ys)))
    axs.set_yticklabels(labs, fontsize=9)
    axs.invert_yaxis()
    axs.set_xlabel("% exploded rounds")
    axs.set_title("Stability", fontsize=12)
    fig.suptitle("Fig 6 — |τ| distributions (ECDF, apo vs bound) and per-run stability",
                 fontsize=14, y=1.0)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


# ----------------------------------------------------------------------------
# Fig 7 — SUPPLEMENTARY: extracellular pocket bodies (3SN6 & 7JJO only)
# ----------------------------------------------------------------------------
def fig7_extracellular(long, stat_ec, out="fig7_extracellular.png"):
    """Extracellular pocket bodies (TM2_EC, TM3, TM6_EC) — the strongest binding
    signal, but only 3SN6 and 7JJO report them, so this is a 2-receptor panel."""
    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.1], height_ratios=[1, 1], hspace=0.45, wspace=0.35)

    # top row: δ heatmaps for the two receptors
    for ci, r in enumerate(EC_RECEPTORS):
        ax = fig.add_subplot(gs[0, ci])
        im = _draw_heatmap(ax, stat_ec, r, EC_BODIES, first_col=(ci == 0))
    cbar = fig.colorbar(im, ax=fig.axes[:2], fraction=0.02, pad=0.02)
    cbar.set_label("Cliff's δ (bound − apo)")

    # top-right: |τ| bound vs apo bars, both receptors
    axb = fig.add_subplot(gs[0, 2])
    x = np.arange(len(EC_BODIES))
    w = 0.2
    for k, r in enumerate(EC_RECEPTORS):
        for j, (state, _) in enumerate(STATES):
            med = [long[(long.receptor == r) & (long.state == state) & (long.segment == b)]["Tmag"].median()
                   for b in EC_BODIES]
            off = (k - 0.5) * 2 * w + (j - 0.5) * w
            axb.bar(x + off, med, w, color=STATE_COLOR[state],
                    alpha=0.9 if r == "3SN6" else 0.55,
                    label=f"{r} {STATE_NAME[state]}")
    axb.set_xticks(x)
    axb.set_xticklabels(EC_BODIES, rotation=20, ha="right")
    axb.yaxis.set_label_position("right")
    axb.yaxis.tick_right()
    axb.set_ylabel("median |τ|  (kJ·mol⁻¹)")
    axb.set_title("Rotational load |τ|")
    axb.legend(fontsize=7, ncol=2)

    # bottom row spanning: moment-arm dumbbell for EC bodies
    axd = fig.add_subplot(gs[1, :])
    yl, yp = [], []
    y = 0.0
    for r in EC_RECEPTORS:
        for body in EC_BODIES:
            apo = long[(long.receptor == r) & (long.state == "noLig") & (long.segment == body)]["arm"].median()
            lig = long[(long.receptor == r) & (long.state == "lig") & (long.segment == body)]["arm"].median()
            axd.plot([apo, lig], [y, y], color="#bbbbbb", lw=2, zorder=1)
            axd.scatter(apo, y, color=STATE_COLOR["noLig"], s=60, zorder=2)
            axd.scatter(lig, y, color=STATE_COLOR["lig"], s=60, zorder=3)
            axd.annotate("", xy=(lig, y), xytext=(apo, y),
                         arrowprops=dict(arrowstyle="->", color="#555", lw=1.4), zorder=2)
            yl.append(f"{r} · {body}")
            yp.append(y)
            y += 1
        y += 0.6
    axd.set_yticks(yp)
    axd.set_yticklabels(yl, fontsize=10)
    axd.invert_yaxis()
    axd.set_xlabel("median moment arm |τ|/|F|  (nm);  apo ● → bound ●, right = more rotational")
    axd.set_title("Force→torque conversion at the extracellular pocket")
    axd.scatter([], [], color=STATE_COLOR["noLig"], label="apo")
    axd.scatter([], [], color=STATE_COLOR["lig"], label="bound")
    axd.legend(loc="lower right")

    fig.suptitle("Fig 7 (supplementary) — Extracellular pocket bodies: the strongest binding signal\n"
                 "TM2_EC / TM3 / TM6_EC, reported only by 3SN6 (β2) and 7JJO (β1); 7DH5 lacks them",
                 fontsize=14, y=0.98)
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


# ----------------------------------------------------------------------------
# Fig 8 — Collective loading modes: PCA of per-body force/torque channels
# ----------------------------------------------------------------------------
# Rotation-invariant / membrane-frame channels only, so the covariance reflects
# concerted LOADING across the bundle, not global receptor tumbling.
PCA_CHANNELS = ["F_axial", "F_lat", "T_spin", "T_rock"]


def _feature_matrix(long, r):
    parts = []
    for state, _ in STATES:
        sub = long[(long.receptor == r) & (long.state == state)]
        mats = []
        for ch in PCA_CHANNELS:
            piv = sub.pivot_table(index="frame", columns="segment", values=ch).reindex(columns=COMMON_BODIES)
            piv.columns = [f"{b}:{ch}" for b in COMMON_BODIES]
            mats.append(piv)
        M = pd.concat(mats, axis=1).dropna()
        M["__state"] = state
        parts.append(M)
    full = pd.concat(parts, axis=0)
    return full.drop(columns="__state"), full["__state"].values


def fig8_pca(long, out="fig8_pca.png"):
    fig, axes = plt.subplots(2, len(RECEPTORS), figsize=(15, 9),
                             gridspec_kw={"height_ratios": [1.3, 1]})
    for ci, r in enumerate(RECEPTORS):
        X, states = _feature_matrix(long, r)
        Xv = X.values.astype(float)
        sd = Xv.std(0)
        sd[sd == 0] = 1.0
        Z = (Xv - Xv.mean(0)) / sd
        U, S, Vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
        scores = U * S
        var = S**2 / (S**2).sum()

        ax = axes[0, ci]
        for state, _ in STATES:
            m = states == state
            ax.scatter(scores[m, 0], scores[m, 1], s=12, alpha=0.5,
                       color=STATE_COLOR[state], label=STATE_NAME[state])
        d1 = cliffs_delta(scores[states == "lig", 0], scores[states == "noLig", 0])
        ax.set_title(f"{RECEPTOR_LABEL[r]}\nPC1 {100*var[0]:.0f}%  PC2 {100*var[1]:.0f}%  "
                     f"| δ(PC1)={d1:+.2f}", fontsize=11)
        ax.set_xlabel("PC1 score")
        if ci == 0:
            ax.set_ylabel("PC2 score")
            ax.legend(fontsize=9)

        # PC1 loading map: channels (rows) × bodies (cols)
        axl = axes[1, ci]
        load = Vt[0].reshape(len(PCA_CHANNELS), len(COMMON_BODIES))
        vmax = np.abs(load).max()
        im = axl.imshow(load, cmap="PuOr_r", vmin=-vmax, vmax=vmax, aspect="auto")
        axl.set_xticks(range(len(COMMON_BODIES)))
        axl.set_xticklabels(COMMON_BODIES, rotation=30, ha="right")
        axl.set_yticks(range(len(PCA_CHANNELS)))
        axl.set_yticklabels(PCA_CHANNELS if ci == 0 else [])
        axl.set_title("PC1 loading (dominant collective mode)", fontsize=10)
        for i in range(load.shape[0]):
            for j in range(load.shape[1]):
                axl.text(j, i, f"{load[i,j]:+.2f}", ha="center", va="center",
                         fontsize=7, color="white" if abs(load[i, j]) > 0.6 * vmax else "black")
    fig.suptitle("Fig 8 — Collective loading modes: PCA of per-body force/torque channels across MC rounds\n"
                 "(membrane-frame channels; scatter = frames colored by state; "
                 "δ(PC1) = separation of bound vs apo along the top mode)", fontsize=13, y=1.0)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


# ----------------------------------------------------------------------------
# Tier-2 geometry (PRELIMINARY, descriptive) — from the sparse existing DCDs
# ----------------------------------------------------------------------------
# Merge the flex sub-segments into the 7 TM helices + H8 for Cα-axis geometry.
TM_MERGE = {
    "TM1": ["TM1"], "TM2": ["TM2_EC", "Na_D2.50", "TM2_IC"], "TM3": ["TM3", "DRY_ionic"],
    "TM4": ["TM4"], "TM5": ["TM5"], "TM6": ["TM6_EC", "CWxP_toggle", "TM6_cyto"],
    "TM7": ["TM7", "NPxxY"], "H8": ["H8"],
}
TM7_ORDER = ["TM1", "TM2", "TM3", "TM4", "TM5", "TM6", "TM7"]
ACT_PAIRS = [("TM3", "TM6"), ("TM5", "TM6"), ("TM6", "TM7"), ("TM3", "TM7"), ("TM1", "TM7"), ("TM2", "TM4")]


def _tm_ranges(receptor):
    seg = dict((lbl, (lo, hi)) for lbl, lo, hi in SEGMENTS[receptor])
    out = {}
    for tm, subs in TM_MERGE.items():
        los = [seg[s][0] for s in subs if s in seg]
        his = [seg[s][1] for s in subs if s in seg]
        if los:
            out[tm] = (min(los), max(his))
    return out


def compute_geometry(receptor, state):
    import mdtraj as md
    import parmed as pmd

    dcd = os.path.join(ROOT, f"{receptor}.{state}.0.dcd")
    prm = os.path.join(ROOT, "examples", "febs", f"{receptor}.{state}.nanodisc.prmtop")
    if not os.path.exists(dcd):
        return None
    parm = pmd.load_file(prm)
    # receptor CA: res.number -> ca atom idx (matches SEGMENTS numbering)
    ca_of_res, started = {}, False
    for res in parm.residues:
        is_aa = len(res.name.strip()) == 3 and res.name.strip() not in ("ACE", "NME", "NHE")
        if is_aa:
            started = True
            for a in res.atoms:
                if a.name == "CA":
                    ca_of_res[res.number] = a.idx
        elif started:
            break
    p31 = np.array([a.idx for a in parm.atoms if a.name == "P31"])
    tmr = _tm_ranges(receptor)
    tm_ca = {tm: np.array([ca_of_res[n] for n in range(lo, hi + 1) if n in ca_of_res])
             for tm, (lo, hi) in tmr.items()}
    tm_ca = {tm: idx for tm, idx in tm_ca.items() if len(idx) >= 4}

    t = md.load(dcd, top=md.load_topology(prm))
    xyz = t.xyz * 10.0  # nm → Å
    up = None
    tilt = {tm: [] for tm in tm_ca}
    dist = {p: [] for p in ACT_PAIRS}
    for f in range(t.n_frames):
        P = xyz[f, p31]
        if up is None:
            up = P[:, 2] > np.median(P[:, 2])
        n = P[up].mean(0) - P[~up].mean(0)
        n /= np.linalg.norm(n)
        cen = {}
        for tm, idx in tm_ca.items():
            pos = xyz[f, idx]
            cen[tm] = pos.mean(0)
            c = pos - pos.mean(0)
            axis = np.linalg.eigh(np.cov(c.T))[1][:, -1]
            cos = abs(np.dot(axis, n))
            tilt[tm].append(np.degrees(np.arccos(np.clip(cos, 0, 1))))
        for a, b in ACT_PAIRS:
            if a in cen and b in cen:
                dist[(a, b)].append(float(np.linalg.norm(cen[a] - cen[b])))
    return {
        "normal": n.tolist(),
        "n_frames": t.n_frames,
        "tilt": {tm: float(np.mean(v)) for tm, v in tilt.items() if v},
        "dist": {p: float(np.mean(v)) for p, v in dist.items() if v},
    }


def fig9_geometry(geo, out="fig9_geometry.png"):
    fig, axes = plt.subplots(2, len(RECEPTORS), figsize=(15, 9))
    for ci, r in enumerate(RECEPTORS):
        # (top) TM tilt, apo vs bound
        ax = axes[0, ci]
        tms = [tm for tm in TM7_ORDER if geo.get((r, "lig")) and tm in geo[(r, "lig")]["tilt"]]
        x = np.arange(len(tms))
        w = 0.38
        for k, (state, _) in enumerate(STATES):
            g = geo.get((r, state))
            vals = [g["tilt"].get(tm, np.nan) for tm in tms] if g else [np.nan] * len(tms)
            ax.bar(x + (k - 0.5) * w, vals, w, color=STATE_COLOR[state], label=STATE_NAME[state], alpha=0.9)
        ax.set_xticks(x)
        ax.set_xticklabels(tms, rotation=30, ha="right")
        nfr = geo.get((r, "lig"), {}).get("n_frames", "?"), geo.get((r, "noLig"), {}).get("n_frames", "?")
        ax.set_title(f"{RECEPTOR_LABEL[r]}  (n={nfr[0]}/{nfr[1]} frames)", fontsize=11)
        if ci == 0:
            ax.set_ylabel("helix tilt vs normal (°)")
            ax.legend(fontsize=9)
        # (bottom) interhelical activation distances
        ax = axes[1, ci]
        pairs = [f"{a}-{b}" for a, b in ACT_PAIRS]
        x = np.arange(len(pairs))
        for k, (state, _) in enumerate(STATES):
            g = geo.get((r, state))
            vals = [g["dist"].get(p, np.nan) for p in ACT_PAIRS] if g else [np.nan] * len(pairs)
            ax.bar(x + (k - 0.5) * w, vals, w, color=STATE_COLOR[state], label=STATE_NAME[state], alpha=0.9)
        ax.set_xticks(x)
        ax.set_xticklabels(pairs, rotation=35, ha="right")
        if ci == 0:
            ax.set_ylabel("Cα-centroid distance (Å)")
    fig.suptitle("Fig 9 (PRELIMINARY, descriptive — 6–10 DCD frames, no statistics) — helix geometry\n"
                 "top: per-helix tilt vs membrane normal; bottom: interhelical activation distances, apo vs bound",
                 fontsize=13, y=1.0)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


# ----------------------------------------------------------------------------
# Poster assembly — composite the panels into one page
# ----------------------------------------------------------------------------
def make_poster(out="poster_betaAR.png"):
    """Composite panels into a claim-structured poster mirroring the abstract's
    four Results claims, with an honest scope footer."""
    from PIL import Image, ImageDraw, ImageFont
    import matplotlib.font_manager as fm

    def _font(sz, bold=True):
        try:
            return ImageFont.truetype(fm.findfont(f"DejaVu Sans{':bold' if bold else ''}"), sz)
        except Exception:
            return ImageFont.load_default()

    def _img(name, w):
        p = os.path.join(ROOT, name)
        if not os.path.exists(p):
            return None
        im = Image.open(p).convert("RGB")
        return im.resize((w, round(im.height * w / im.width)), Image.LANCZOS)

    M, GAP = 60, 45
    COLW = 1500
    FULLW = 2 * COLW + GAP
    W = FULLW + 2 * M
    HEAD, SEC, FOOT = 210, 90, 470

    # layout program: rows rendered top→down
    #   ("sec", text) | ("full", name) | ("pair", left, right|None)
    prog = [
        ("sec", "① Redistribution of membrane forces & torques on agonist binding"),
        ("full", "fig2_heatmap.png"),
        ("pair", "fig1_fingerprint.png", "fig4_decomp.png"),
        ("pair", "fig3_microswitch.png", "fig7_extracellular.png"),
        ("sec", "②–③ Pivot regions & binding-site geometry  (PRELIMINARY — 6–10 DCD frames, descriptive)"),
        ("full", "fig9_geometry.png"),
        ("sec", "④ Subtype-specific differences & collective loading modes"),
        ("pair", "fig5_forest.png", "fig8_pca.png"),
        ("sec", "Quality control"),
        ("pair", "fig6_ecdf.png", None),
    ]

    # pre-scale & measure
    rows, body_h = [], 0
    for row in prog:
        if row[0] == "sec":
            rows.append(("sec", row[1], SEC))
            body_h += SEC + GAP
        elif row[0] == "full":
            im = _img(row[1], FULLW)
            if im:
                rows.append(("full", im, im.height))
                body_h += im.height + GAP
        else:
            li, ri = _img(row[1], COLW), _img(row[2], COLW) if row[2] else None
            h = max(li.height if li else 0, ri.height if ri else 0)
            rows.append(("pair", (li, ri), h))
            body_h += h + GAP

    H = HEAD + body_h + FOOT
    canvas = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(canvas)

    d.text((M, 40), "Mechanical drivers of transmembrane-helix conformation in three β-adrenergic receptors",
           fill="black", font=_font(50))
    d.text((M, 108), "Per-body spatial reaction forces (F, τ) from Robosample blocked-Gibbs / HMC  ·  "
           "β2-AR 3SN6 · β3-AR 7DH5 · β1-AR 7JJO  ·  apo vs agonist  ·  lipid nanodisc, implicit solvent",
           fill="#555555", font=_font(28, bold=False))

    y = HEAD
    for kind, payload, h in rows:
        if kind == "sec":
            d.rectangle([M, y, M + FULLW, y + h], fill="#2b3a55")
            d.text((M + 24, y + h // 2 - 22), payload, fill="white", font=_font(38))
        elif kind == "full":
            canvas.paste(payload, (M, y))
        else:
            li, ri = payload
            if li:
                canvas.paste(li, (M, y))
            if ri:
                canvas.paste(ri, (M + COLW + GAP, y))
        y += h + GAP

    # scope / caveats footer
    d.rectangle([M, y, M + FULLW, y + FOOT - GAP], outline="#2b3a55", width=4)
    fx, fy = M + 28, y + 22
    d.text((fx, fy), "Scope, methods & caveats", fill="#2b3a55", font=_font(34))
    lines = [
        "•  Signal that replicates: extracellular TM2_EC gains large rotational load on binding (δ|τ| +0.56 β2 / +0.63 β1, mostly in-plane rocking) — the",
        "    strongest, most reproducible effect. Intracellular bodies (TM7, TM6_cyto, NPxxY, H8) redistribute load only weakly and not sign-consistently.",
        "•  Units: raw native reporter units (force kJ·mol⁻¹·nm⁻¹, torque kJ·mol⁻¹). Effect size = Cliff's δ (bound−apo); CIs = moving-block bootstrap.",
        "    Exploded MC rounds (|F|>1e5, accepted clashing moves) excised. Membrane frame uses fixed n̂=[0,0,1] (nanodiscs verified flat).",
        "•  Data scope, NOT yet in this poster:  ligand set is apo vs one AGONIST per receptor (no antagonist / ligand-class comparison);  family is",
        "    β1/β2/β3 only (no α);  claims ②③ use sparse legacy DCDs (6–10 frames, descriptive, no statistics);  interhelical distance is an",
        "    accessibility PROXY, not pocket volume/SASA.  Force-force PCA (④) approximates collective modes from the force side, not coordinate space.",
        "•  Next: cadence-matched DCD re-run (per-round) unlocks quantitative geometry + force↔geometry coupling; antagonist runs unlock ligand-class bias.",
    ]
    for i, ln in enumerate(lines):
        d.text((fx, fy + 52 + i * 40), ln, fill="#222222", font=_font(25, bold=False))

    canvas.save(os.path.join(ROOT, out))
    print(f"Saved {out}  ({W}×{H}px)")


# ----------------------------------------------------------------------------
# Interpretations
# ----------------------------------------------------------------------------
def write_interpretations(stat, expl, path="betaAR_interpretations.txt"):
    L = []
    L.append("β-adrenergic GPCR reaction-force analysis — bound vs apo")
    L.append("=" * 62)
    L.append("Bodies analysed (common to all 6 runs): " + ", ".join(COMMON_BODIES))
    L.append("Units: |F| kJ·mol⁻¹·nm⁻¹, |τ| kJ·mol⁻¹, arm nm. δ = Cliff's (bound−apo).")
    L.append("Exploded rounds (|F|>1e5) excised per run:")
    for _, row in expl.iterrows():
        L.append(f"  {row.receptor} {STATE_NAME[row.state]:5s}: "
                 f"{row.exploded}/{row.total} ({100*row.exploded/max(row.total,1):.1f}%)")
    L.append("")
    L.append("Cross-receptor consistency of δ (sign in all 3 receptors ⇒ reproducible):")
    for col, name, _ in METRICS:
        L.append(f"\n[{name}]")
        for body in COMMON_BODIES:
            ds = stat[(stat.body == body) & (stat.metric == col)].set_index("receptor")["cliffs_delta"]
            ds = ds.reindex(RECEPTORS)
            signs = {np.sign(v) for v in ds.dropna()}
            tag = "  CONSISTENT" if len(signs) == 1 and ds.notna().all() else ""
            vals = "  ".join(f"{r}:{ds[r]:+.2f}" if np.isfinite(ds[r]) else f"{r}:NA" for r in RECEPTORS)
            L.append(f"  {body:10s} {vals}{tag}")
    with open(path, "w") as f:
        f.write("\n".join(L))
    print(f"Saved {path}")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    print("Loading + labeling 6 reaction CSVs ...")
    long, expl = load_all()
    print(f"Total clean rows: {len(long)}  bodies: {sorted(long.segment.unique())}")

    print("Computing robust statistics (this includes bootstrap CIs) ...")
    stat = compute_stats(long)
    stat_ec = compute_stats(long, receptors=EC_RECEPTORS, bodies=EC_BODIES)
    pd.concat([stat, stat_ec], ignore_index=True).to_csv("betaAR_stats.csv", index=False)
    print("Saved betaAR_stats.csv")

    fig1_fingerprint(long)
    fig2_heatmap(stat)
    fig3_microswitch(long)
    fig4_decomp(long)
    fig5_forest(stat)
    fig6_ecdf(long, expl)
    fig7_extracellular(long, stat_ec)
    fig8_pca(long)

    print("Computing Tier-2 geometry from sparse DCDs (descriptive) ...")
    geo = {}
    for r in RECEPTORS:
        for s, _ in STATES:
            try:
                g = compute_geometry(r, s)
                if g:
                    geo[(r, s)] = g
                    print(f"  {r} {s}: {g['n_frames']} frames, normal={np.round(g['normal'], 2)}")
            except Exception as e:
                print(f"  {r} {s}: geometry skipped ({e})")
    if geo:
        fig9_geometry(geo)
    else:
        print("  No geometry computed — Fig 9 skipped.")

    write_interpretations(stat, expl)
    make_poster()
    print("\nDone. fig1..fig9 + poster_betaAR.png + betaAR_stats.csv + betaAR_interpretations.txt")
