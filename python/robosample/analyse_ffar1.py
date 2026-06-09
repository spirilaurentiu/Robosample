"""
analyze_worldD_poster.py  -  FFAR1 World-D  |  poster-quality figures
======================================================================
Fixes vs previous version
--------------------------
  * Rotation angle: unwrapped with np.unwrap (no ±180° jumps)
  * Time series: aggregated per TM (7 lines, not 15) + heatmap option
  * Force decomposition: cleaner layout + written interpretation
  * All figures: poster-grade fonts, DPI=200, tight layout
  * RMSD in internal coordinates (u values) added
  * Interpretations printed to  interpretations.txt

Run
---
  python analyze_worldD_poster.py [vectors.dat] [ffar1.prmtop] [ffar1.dcd] [stride]
"""

import sys
import warnings

warnings.filterwarnings("ignore")

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal, stats

# -- global poster style -------------------------------------------------------
plt.rcParams.update(
    {
        "font.size": 13,
        "axes.titlesize": 14,
        "axes.labelsize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 10,
        "figure.dpi": 200,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
        "lines.linewidth": 1.4,
    }
)

# -- identity ------------------------------------------------------------------
DIAG_RESID = [14, 37, 41, 62, 86, 110, 130, 146, 182, 209, 223, 238, 247, 257, 290]
HELIX_LABEL = [
    "TM1-N",
    "TM1-C",
    "TM2-N",
    "TM2-C",
    "TM3-N",
    "TM3-C",
    "TM4-N",
    "TM4-C",
    "TM5-N",
    "TM5-C",
    "TM6-N",
    "TM6-M",
    "TM6-C",
    "TM7-N",
    "TM7-C",
]
BOUNDARY_MAP = {
    "TM1": ["TM1-N", "TM1-C"],
    "TM2": ["TM2-N", "TM2-C"],
    "TM3": ["TM3-N", "TM3-C"],
    "TM4": ["TM4-N", "TM4-C"],
    "TM5": ["TM5-N", "TM5-C"],
    "TM6": ["TM6-N", "TM6-M", "TM6-C"],
    "TM7": ["TM7-N", "TM7-C"],
}
HELIX_RANGES = {
    "TM1": (14, 37),
    "TM2": (41, 65),
    "TM3": (72, 110),
    "TM4": (130, 147),
    "TM5": (182, 210),
    "TM6": (223, 248),
    "TM7": (257, 290),
}
TM_LIST = list(BOUNDARY_MAP.keys())
TM_COLOR = {
    "TM1": "#e41a1c",
    "TM2": "#377eb8",
    "TM3": "#4daf4a",
    "TM4": "#984ea3",
    "TM5": "#ff7f00",
    "TM6": "#a65628",
    "TM7": "#f781bf",
}


def _c(tm):
    return TM_COLOR.get(tm.split("-")[0], "gray")


# -- load / aggregate ----------------------------------------------------------
def load_vectors(path):
    df = pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        names=[
            "frame",
            "inboard_idx",
            "outboard_idx",
            "fx",
            "fy",
            "fz",
            "tx",
            "ty",
            "tz",
            "u",
            "uDot",
        ],
    )

    raw_unique = df.frame.nunique()
    rows_per_block = (df.frame == df.frame.iloc[0]).sum()

    # -- Auto-detect whether "frame" col is really a pair index ---------------
    # Signature: small number of unique values (≈ n_pairs) that repeat
    # throughout the file, rather than a monotonically advancing frame counter.
    expected_frames = len(df) // rows_per_block
    if raw_unique <= rows_per_block:
        print(
            f"  ⚠  'frame' column has only {raw_unique} unique values "
            f"(looks like pair index, not frame number)."
        )
        print(
            f"  → Reconstructing frame numbers: "
            f"{expected_frames} frames × {rows_per_block} pairs"
        )
        df["pair_idx"] = df["frame"]
        df["frame"] = np.repeat(np.arange(expected_frames), rows_per_block)
    else:
        print(f"  ✓  frame column looks correct: {raw_unique} unique frames")

    # -- Helix label mapping ---------------------------------------------------
    idx_order = df[df.frame == df.frame.iloc[0]]["inboard_idx"].tolist()
    df["helix"] = df["inboard_idx"].map(
        {idx: HELIX_LABEL[i] for i, idx in enumerate(idx_order[: len(HELIX_LABEL)])}
    )

    null_helix = df.helix.isna().sum()
    if null_helix > 0:
        print(
            f"  ⚠  {null_helix} rows have no helix label "
            f"({null_helix / len(df) * 100:.1f}% of data) — check HELIX_LABEL vs inboard_idx values."
        )

    df["F_mag"] = np.sqrt(df.fx**2 + df.fy**2 + df.fz**2)
    df["T_mag"] = np.sqrt(df.tx**2 + df.ty**2 + df.tz**2)
    return df


def decompose_wrt_normal(df, geo_df):
    """
    Decompose force (and torque) vectors into components parallel and
    perpendicular to the per-frame membrane normal stored in geo_df.

    F_axial  = |F · n̂|          — along membrane normal (compression/extension)
    F_lat    = |F - F_axial·n̂| — in membrane plane (shear)

    Same decomposition applied to torque:
    T_axial  = |τ · n̂|          — spinning around membrane normal
    T_lat    = |τ - T_axial·n̂| — rocking / tilting torque
    """
    # geo_df has one row per (frame, helix); membrane normal is the same for all
    # helices in a frame, so just take one row per frame
    normals = (
        geo_df.groupby("frame")[["norm_x", "norm_y", "norm_z"]]
        .first()
        .reset_index()
        .rename(columns={"norm_x": "nx", "norm_y": "ny", "norm_z": "nz"})
    )

    df = df.merge(normals, on="frame", how="left")

    missing = df[["nx", "ny", "nz"]].isna().any(axis=1).sum()
    if missing > 0:
        print(
            f"  Warning: {missing} rows have no matching membrane normal "
            f"(frames in vectors.dat not in geo_df). Falling back to Z=[0,0,1]."
        )
        df[["nx", "ny", "nz"]] = df[["nx", "ny", "nz"]].fillna(
            {"nx": 0.0, "ny": 0.0, "nz": 1.0}
        )

    # -- Force decomposition ---------------------------------------------------
    # scalar projection onto normal: F · n̂
    df["F_axial_signed"] = df.fx * df.nx + df.fy * df.ny + df.fz * df.nz
    df["F_axial"] = df["F_axial_signed"].abs()

    # lateral component: F - (F·n̂)n̂, then take magnitude
    df["_flx"] = df.fx - df.F_axial_signed * df.nx
    df["_fly"] = df.fy - df.F_axial_signed * df.ny
    df["_flz"] = df.fz - df.F_axial_signed * df.nz
    df["F_lat"] = np.sqrt(df._flx**2 + df._fly**2 + df._flz**2)

    # -- Torque decomposition (same logic) ------------------------------------
    df["T_axial_signed"] = df.tx * df.nx + df.ty * df.ny + df.tz * df.nz
    df["T_axial"] = df["T_axial_signed"].abs()

    df["_tlx"] = df.tx - df.T_axial_signed * df.nx
    df["_tly"] = df.ty - df.T_axial_signed * df.ny
    df["_tlz"] = df.tz - df.T_axial_signed * df.nz
    df["T_lat"] = np.sqrt(df._tlx**2 + df._tly**2 + df._tlz**2)

    # drop scratch columns
    df = df.drop(
        columns=[
            "_flx",
            "_fly",
            "_flz",
            "_tlx",
            "_tly",
            "_tlz",
            "F_axial_signed",
            "T_axial_signed",
        ]
    )
    return df


def aggregate_by_tm(df):
    """Mean of all boundary residues per TM helix → one value per (frame, TM)."""
    rows = []
    for tm, labels in BOUNDARY_MAP.items():
        sub = (
            df[df.helix.isin(labels)]
            .groupby("frame")[["F_mag", "T_mag", "F_lat", "F_axial", "u", "uDot"]]
            .mean()
            .reset_index()
        )
        sub["helix"] = tm
        rows.append(sub)
    return pd.concat(rows, ignore_index=True)


# -----------------------------------------------------------------------------
# Fig 1 - Force & torque decomposition (bar chart)
# -----------------------------------------------------------------------------
def fig1_force_decomposition(df, tm_df, out="fig1_force_decomposition.png"):
    grp = tm_df.groupby("helix")[["F_lat", "F_axial", "T_mag"]]
    means = grp.mean().reindex(TM_LIST)
    stds = grp.std().reindex(TM_LIST)
    n = grp.size().reindex(TM_LIST)
    sems = stds.div(np.sqrt(n), axis=0)

    x, w = np.arange(len(TM_LIST)), 0.25
    fig, ax = plt.subplots(figsize=(11, 5))

    for i, (col, label, color) in enumerate(
        [
            ("F_lat", "Lateral |F|  (in-plane)", "steelblue"),
            ("F_axial", "Axial |F|  (membrane normal)", "darkorange"),
            ("T_mag", "Torque |τ|", "forestgreen"),
        ]
    ):
        offset = (i - 1) * w
        ax.bar(
            x + offset,
            means[col],
            w,
            yerr=sems[col],
            capsize=4,
            label=label,
            color=color,
            alpha=0.88,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(TM_LIST)
    ax.set_ylabel("Mean magnitude  (kcal mol⁻¹ Å⁻¹)")
    ax.set_title(
        "Fig 1 - Mechanical loading at each transmembrane helix\n"
        "(force decomposed relative to membrane normal)"
    )
    ax.legend(loc="upper right", fontsize=9)
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"Saved {out}")
    return means, stds


# -----------------------------------------------------------------------------
# Fig 2 - Helix tilt + unwrapped rotation
# -----------------------------------------------------------------------------
def _delta_tilt(series):
    """Tilt relative to frame 0. No unwrapping needed — tilt has no periodicity."""
    vals = series.values.copy()
    return vals - vals[0]


def _delta_rotation(series):
    """
    Express rotation relative to frame 0, then unwrap.
    Result: accumulated rotation change in degrees, starting at 0.
    """
    vals = series.values.copy()
    delta = vals - vals[0]
    delta = (delta + 180) % 360 - 180
    return np.rad2deg(np.unwrap(np.deg2rad(delta)))


def fig2_helix_tilt(geo_df, out="fig2_helix_tilt.png"):
    fig, axes = plt.subplots(2, 1, figsize=(30, 10), sharex=True)

    for tm in TM_LIST:
        sub = geo_df[geo_df.helix == tm].sort_values("frame")
        c = TM_COLOR[tm]
        axes[0].plot(sub.frame, _delta_tilt(sub.tilt), color=c, label=tm)
        axes[1].plot(sub.frame, _delta_rotation(sub.rotation), color=c, label=tm)

    axes[0].set_ylabel("ΔTilt  (°, relative to frame 0)")
    axes[0].set_xlabel("Frame")
    axes[0].set_title("Fig 2a - Relative helix tilt vs membrane normal")
    axes[0].axhline(0, color="k", lw=0.5, ls="--")  # reference line

    axes[1].set_ylabel("ΔRotation  (°, relative to frame 0)")
    axes[1].set_xlabel("Frame")
    axes[1].set_title("Fig 2b - Relative helix rotation in membrane plane")
    axes[1].axhline(0, color="k", lw=0.5, ls="--")

    for ax in axes:
        ax.legend(ncol=4, fontsize=9)

    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"Saved {out}")


# -----------------------------------------------------------------------------
# Fig 3 - Force, torque, u, uDot heatmaps (4 panels, shared x-axis)
# -----------------------------------------------------------------------------
def fig3_heatmap_combined(tm_df, out="fig3_heatmap_combined.png"):
    """
    Five stacked heatmaps sharing a common frame axis:
      3a  |F| lateral   sequential, vmin=0
      3b  |F| axial     sequential, vmin=0
      3c  |τ|           sequential, vmin=0
      3d  u             diverging, symmetric about 0
      3e  u̇             diverging, symmetric about 0
    """
    CMAP_SEQ = "YlOrRd"
    CMAP_DIV = "RdBu_r"

    panels = [
        (
            "F_lat",
            "Fig 3a - Lateral force (in membrane plane) per helix",
            "|F_lat|  (kcal mol⁻¹ Å⁻¹)",
            CMAP_SEQ,
            False,
        ),
        (
            "F_axial",
            "Fig 3b - Axial force (along membrane normal) per helix",
            "|F_axial|  (kcal mol⁻¹ Å⁻¹)",
            CMAP_SEQ,
            False,
        ),
        (
            "T_mag",
            "Fig 3c - Torque magnitude |τ| per helix",
            "|τ|  (kcal mol⁻¹)",
            CMAP_SEQ,
            False,
        ),
        (
            "u",
            "Fig 3d - Generalised coordinate u (φ dihedral, mean per TM)",
            "u  (rad)",
            CMAP_DIV,
            True,
        ),
        (
            "uDot",
            "Fig 3e - Generalised velocity u̇ (mean per TM)",
            "u̇  (rad ps⁻¹)",
            CMAP_DIV,
            True,
        ),
    ]

    fig, axes = plt.subplots(5, 1, figsize=(14, 15), sharex=True)

    for ax, (col, title, cbar_label, cmap, diverging) in zip(axes, panels):
        mat = tm_df.pivot(index="helix", columns="frame", values=col).reindex(TM_LIST)
        vals = mat.values

        if diverging:
            lim = np.nanpercentile(np.abs(vals), 98)
            vmin, vmax = -lim, lim
        else:
            vmin = 0.0
            vmax = np.nanpercentile(vals, 98)

        im = ax.imshow(
            vals,
            aspect="auto",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
        )
        ax.set_yticks(range(len(TM_LIST)))
        ax.set_yticklabels(TM_LIST, fontsize=8)
        ax.set_title(title, fontsize=10)
        plt.colorbar(im, ax=ax, label=cbar_label, fraction=0.025, pad=0.02)

    axes[-1].set_xlabel("Frame")
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"Saved {out}")


# -----------------------------------------------------------------------------
# Fig 5 - Cross-correlation |F| ↔ u  (summary heatmap of peak lag + r)
# -----------------------------------------------------------------------------
def fig5_xcorr_summary(
    tm_df,
    out_detail="fig5a_xcorr_detail.png",
    out_summary="fig5b_xcorr_summary.png",
    max_lag=40,
):
    """
    Panel a: full cross-correlogram per TM (7 panels, readable).
    Panel b: summary heatmap showing peak lag and peak r for poster.
    """
    ci = 1.96
    lags_all, xcorr_all = {}, {}
    peak_lag, peak_r = {}, {}

    for tm in TM_LIST:
        sub = tm_df[tm_df.helix == tm].sort_values("frame")
        f_arr = sub.F_mag.values
        u_arr = sub.u.values
        f_z = f_arr - f_arr.mean()
        u_z = u_arr - u_arr.mean()
        n = len(f_z)
        lags = signal.correlation_lags(n, n, mode="full")
        xc = signal.correlate(f_z, u_z, mode="full")
        denom = np.std(f_z) * np.std(u_z) * n
        xc = xc / denom if denom > 0 else xc
        mask = np.abs(lags) <= max_lag
        lags_all[tm] = lags[mask]
        xcorr_all[tm] = xc[mask]
        pk_idx = np.argmax(np.abs(xc[mask]))
        peak_lag[tm] = int(lags[mask][pk_idx])
        peak_r[tm] = float(xc[mask][pk_idx])

    # (a) detailed per-TM correlograms
    fig, axes = plt.subplots(
        len(TM_LIST), 1, figsize=(11, 2.2 * len(TM_LIST)), sharex=True
    )
    for ax, tm in zip(axes, TM_LIST):
        n = tm_df[tm_df.helix == tm].shape[0]
        ax.bar(lags_all[tm], xcorr_all[tm], width=1, color=TM_COLOR[tm], alpha=0.75)
        ax.axhline(0, color="k", lw=0.6)
        for sign in [1, -1]:
            ax.axhline(sign * ci / np.sqrt(n), color="gray", lw=0.8, ls="--")
        ax.set_ylabel(tm, fontsize=10, rotation=0, labelpad=48)
        ax.set_ylim(-1.05, 1.05)
    axes[-1].set_xlabel("Lag (frames)  — positive: |F| leads u")
    fig.suptitle(
        "Fig 5a - Cross-correlation  |F| ↔ u  per TM helix", y=1.002, fontsize=13
    )
    plt.tight_layout()
    plt.savefig(out_detail, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_detail}")

    # (b) summary heatmap: rows = TM, cols = [peak lag, peak r]
    summary = pd.DataFrame({"Peak lag (frames)": peak_lag, "Peak r": peak_r}).T[TM_LIST]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3))
    for ax, row_lbl in zip(axes, ["Peak lag (frames)", "Peak r"]):
        data = summary.loc[row_lbl].values.astype(float).reshape(1, -1)
        cmap = "coolwarm" if "r" in row_lbl else "PiYG"
        vmax = max(abs(data.min()), abs(data.max()))
        im = ax.imshow(data, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(len(TM_LIST)))
        ax.set_xticklabels(TM_LIST)
        ax.set_yticks([])
        ax.set_title(row_lbl)
        plt.colorbar(im, ax=ax, fraction=0.06, pad=0.04)
        for j, v in enumerate(data[0]):
            ax.text(j, 0, f"{v:.2f}", ha="center", va="center", fontsize=10, color="k")
    fig.suptitle("Fig 5b - Cross-correlation summary: peak lag and peak r", fontsize=13)
    plt.tight_layout()
    plt.savefig(out_summary, dpi=200)
    plt.close()
    print(f"Saved {out_summary}")
    return peak_lag, peak_r


# -----------------------------------------------------------------------------
# Fig 6 - Force vs tilt scatter
# -----------------------------------------------------------------------------
def fig6_force_vs_tilt(tm_df, geo_df, out="fig6_force_vs_tilt.png"):
    ncols = 4
    nrows = 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4 * nrows))
    axes = axes.flatten()
    results = {}
    for i, (tm, _) in enumerate(BOUNDARY_MAP.items()):
        ax = axes[i]
        mean_F = tm_df[tm_df.helix == tm].set_index("frame")["F_mag"]
        tilt_s = geo_df[geo_df.helix == tm].set_index("frame")["tilt"]
        common = mean_F.index.intersection(tilt_s.index)
        if len(common) < 5:
            ax.set_title(f"{tm}  (n<5)")
            continue
        x, y = mean_F[common].values, tilt_s[common].values
        r, p = stats.pearsonr(x, y)
        results[tm] = (r, p, len(common))
        star = (
            "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
        )
        ax.scatter(x, y, s=10, alpha=0.35, color=TM_COLOR[tm])
        xx = np.linspace(x.min(), x.max(), 100)
        ax.plot(xx, np.polyval(np.polyfit(x, y, 1), xx), "k--", lw=1.2)
        ax.set_xlabel("|F|  mean")
        ax.set_ylabel("Tilt  (°)")
        ax.set_title(f"{tm}   r={r:+.2f} {star}")
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    fig.suptitle(
        "Fig 6 - Mean boundary force vs helix tilt  (mechanical-geometric coupling)",
        fontsize=13,
    )
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"Saved {out}")
    return results


# -----------------------------------------------------------------------------
# Fig 7 - Interhelical distances
# -----------------------------------------------------------------------------
PAIRS = [
    ("TM3", "TM5"),
    ("TM3", "TM6"),
    ("TM5", "TM6"),
    ("TM1", "TM7"),
    ("TM2", "TM4"),
    ("TM6", "TM7"),
]


def fig7_interhelical(geo_df, out="fig7_interhelical.png"):
    COLOR = "steelblue"

    try:
        pivoted = geo_df.pivot(
            index="frame", columns="helix", values=["cx", "cy", "cz"]
        )
    except Exception as e:
        print(f"Pivot failed: {e}")
        return

    frames = sorted(geo_df.frame.unique())
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
    axes = axes.flatten()
    stats_rows = []

    for ax, (h1, h2) in zip(axes, PAIRS):
        try:
            c1 = pivoted.loc[:, (["cx", "cy", "cz"], h1)].values
            c2 = pivoted.loc[:, (["cx", "cy", "cz"], h2)].values
        except KeyError:
            ax.set_visible(False)
            continue

        d = np.linalg.norm(c1 - c2, axis=1)
        ax.plot(frames, d, color=COLOR, lw=1.2)
        ax.axhline(d.mean(), color="k", ls="--", lw=0.8, label=f"μ={d.mean():.1f} Å")
        ax.fill_between(
            frames,
            d.mean() - d.std(),
            d.mean() + d.std(),
            alpha=0.15,
            color=COLOR,
        )
        ax.set_title(f"{h1}-{h2}")
        ax.set_ylabel("Distance  (Å)")
        ax.legend(fontsize=8)
        stats_rows.append(
            {
                "pair": f"{h1}-{h2}",
                "mean_A": d.mean(),
                "std_A": d.std(),
                "min_A": d.min(),
                "max_A": d.max(),
            }
        )

    axes[-1].set_xlabel("Frame")
    axes[-2].set_xlabel("Frame")
    fig.suptitle("Fig 5 - Inter-helix centroid distances", fontsize=13)
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"Saved {out}")
    return pd.DataFrame(stats_rows)


# -----------------------------------------------------------------------------
# Fig 8 - RMSD in internal coordinates (u)
# -----------------------------------------------------------------------------
def fig8_internal_rmsd(df, out="fig8_internal_rmsd.png"):
    """
    For each TM: RMSD of u relative to frame-0 value.
    This is the internal-coordinate equivalent of Cartesian RMSD —
    it measures how much the hinge dihedral has drifted from the start.
    """
    u0 = df[df.frame == df.frame.min()].set_index("helix")["u"]

    fig, ax = plt.subplots(figsize=(11, 5))
    for tm, labels in BOUNDARY_MAP.items():
        sub = df[df.helix.isin(labels)]
        rmsd_per_frame = sub.groupby("frame").apply(
            lambda g: np.sqrt(
                (
                    (g.set_index("helix")["u"] - u0.reindex(g.helix.values).values) ** 2
                ).mean()
            )
        )
        ax.plot(
            rmsd_per_frame.index, rmsd_per_frame.values, label=tm, color=TM_COLOR[tm]
        )

    ax.set_xlabel("Frame")
    ax.set_ylabel("RMSD  u  (rad)")
    ax.set_title("Fig 8 - Internal-coordinate RMSD  (φ dihedral drift per TM helix)")
    ax.legend(ncol=4)
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"Saved {out}")


# -----------------------------------------------------------------------------
# Statistics summary
# -----------------------------------------------------------------------------
def compute_statistics(tm_df):
    cols = ["F_mag", "F_lat", "F_axial", "T_mag", "u", "uDot"]
    rows = []
    for tm in TM_LIST:
        sub = tm_df[tm_df.helix == tm]
        row = {"helix": tm}
        for c in cols:
            row[f"{c}_mean"] = sub[c].mean()
            row[f"{c}_std"] = sub[c].std()
            row[f"{c}_skew"] = float(stats.skew(sub[c].dropna()))
        rows.append(row)
    sdf = pd.DataFrame(rows).set_index("helix")
    sdf.to_csv("helix_statistics.csv")
    return sdf


# -----------------------------------------------------------------------------
# Helix geometry (MDTraj) — dynamic membrane normal from phosphate COMs
# -----------------------------------------------------------------------------


def _membrane_normal(xyz_all, p_indices, p_upper_mask):
    """
    Compute membrane normal from phosphate atom positions.

    Parameters
    ----------
    xyz_all       : (N_atoms, 3) array — full frame coordinates (Å)
    p_indices     : 1-D array of phosphorus atom indices in xyz_all
    p_upper_mask  : boolean mask (len = len(p_indices)) — True for upper leaflet P

    Returns
    -------
    normal : (3,) unit vector pointing from lower → upper leaflet
    midplane_z : float — Z coordinate of bilayer midplane (Å)
    """
    p_xyz = xyz_all[p_indices]
    upper_com = p_xyz[p_upper_mask].mean(0)
    lower_com = p_xyz[~p_upper_mask].mean(0)
    normal = upper_com - lower_com
    normal /= np.linalg.norm(normal)
    midplane_z = (upper_com[2] + lower_com[2]) / 2.0
    return normal, midplane_z


def _helix_axis(pos):
    """
    PCA on Cα positions → helix long axis (unit vector).
    Sign is fixed so axis always points toward +Z convention
    (will be re-signed vs membrane normal later).
    """
    c = pos.mean(0)
    _, vecs = np.linalg.eigh(np.cov((pos - c).T))
    axis = vecs[:, -1]  # largest eigenvalue = long axis
    if axis[2] < 0:
        axis = -axis
    return axis


def _tilt_vs_normal(axis, normal):
    """
    Tilt angle (degrees) between helix axis and membrane normal.
    Uses absolute dot product so the result is always in [0°, 90°].
    """
    cos_theta = np.dot(axis, normal)
    # ensure axis points same general direction as normal for consistency
    if cos_theta < 0:
        axis = -axis
        cos_theta = -cos_theta
    tilt = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
    return tilt, axis  # return (possibly flipped) axis for rotation calc


def _helix_rotation_vs_normal(axis, normal):
    """
    Azimuthal rotation of helix axis around the membrane normal.

    Projects axis onto the plane perpendicular to the membrane normal,
    then measures the angle of that projection. This is meaningful even
    when the membrane normal is not exactly Z.

    Returns angle in degrees, range (-180, 180].
    """
    # Component of axis in the membrane plane
    proj = axis - np.dot(axis, normal) * normal
    norm_proj = np.linalg.norm(proj)
    if norm_proj < 1e-6:
        return 0.0  # axis is perfectly parallel to normal — angle undefined
    proj /= norm_proj

    # We need a consistent reference direction in the membrane plane.
    # Use the component of global X that is perpendicular to the membrane normal.
    x_ref = np.array([1.0, 0.0, 0.0])
    x_ref = x_ref - np.dot(x_ref, normal) * normal
    x_ref_norm = np.linalg.norm(x_ref)
    if x_ref_norm < 1e-6:
        # normal is parallel to X; fall back to Y
        x_ref = np.array([0.0, 1.0, 0.0])
        x_ref = x_ref - np.dot(x_ref, normal) * normal
        x_ref /= np.linalg.norm(x_ref)
    else:
        x_ref /= x_ref_norm

    y_ref = np.cross(normal, x_ref)  # completes right-handed frame in membrane plane

    rotation = np.degrees(np.arctan2(np.dot(proj, y_ref), np.dot(proj, x_ref)))
    return rotation


def _assign_leaflets_first_frame(chunk0, p_indices_local, midplane_z_guess=None):
    """
    Assign phosphorus atoms to upper/lower leaflet using the first frame.
    Returns a boolean mask (True = upper leaflet).
    If midplane_z_guess is None, uses the median P z-coordinate.
    """
    p_xyz = chunk0.xyz[0][p_indices_local] * 10.0  # nm → Å
    if midplane_z_guess is None:
        midplane_z_guess = np.median(p_xyz[:, 2])
    upper_mask = p_xyz[:, 2] > midplane_z_guess
    n_up, n_lo = upper_mask.sum(), (~upper_mask).sum()
    print(f"  Leaflet assignment: {n_up} upper P, {n_lo} lower P atoms")
    if min(n_up, n_lo) == 0:
        raise ValueError(
            "All phosphates in one leaflet — check residue name or system centering."
        )
    return upper_mask


def compute_helix_geometry(
    top="ffar1.prmtop",
    traj="ffar1.dcd",
    stride=1,
    lipid_resnames=("POPC", "POPE", "CHOL", "CHL1", "DPPC", "PC"),  # extend as needed
):
    try:
        import mdtraj as md
    except ImportError:
        print("MDTraj not found — pip install mdtraj")
        return None

    # -- Load topology --------------------------------------------------------
    try:
        top_obj = md.load_topology(top)
    except Exception as e:
        print(f"Topology load failed: {e}")
        return None

    # -- Collect Cα indices per TM helix --------------------------------------
    ca_indices = {}
    for tm, (r0, r1) in HELIX_RANGES.items():
        idx = [
            a.index
            for a in top_obj.atoms
            if a.name == "CA" and r0 <= a.residue.resSeq <= r1
        ]
        if len(idx) >= 4:
            ca_indices[tm] = np.array(idx)
        else:
            print(f"  Warning: {tm} has only {len(idx)} Cα — skipped.")

    # -- Collect phosphorus indices --------------------------------------------
    p_indices = np.array(
        [
            a.index
            for a in top_obj.atoms
            if a.name == "P31" and a.residue.name in lipid_resnames
        ]
    )
    if len(p_indices) == 0:
        raise ValueError(
            f"No phosphorus atoms found for residues {lipid_resnames}. "
            "Check lipid_resnames or use atom name 'P' in your force field."
        )
    print(f"Found {len(p_indices)} phosphorus atoms for membrane normal.")

    # -- Atom indices to load (Cα + P) ----------------------------------------
    all_ca = np.unique(np.concatenate(list(ca_indices.values())))
    all_load = np.unique(np.concatenate([all_ca, p_indices]))

    # Local indices within the reduced atom set
    ca_local = {tm: np.searchsorted(all_load, ca_indices[tm]) for tm in ca_indices}
    p_local = np.searchsorted(all_load, p_indices)

    # -- Leaflet assignment from first chunk -----------------------------------
    first_chunk = next(
        iter(
            md.iterload(
                traj, top=top_obj, stride=stride, atom_indices=all_load, chunk=1
            )
        )
    )
    upper_mask = _assign_leaflets_first_frame(first_chunk, p_local)

    # -- Main trajectory loop --------------------------------------------------
    records, frame_counter = [], 0

    for chunk in md.iterload(traj, top=top_obj, stride=stride, atom_indices=all_load):
        for f in range(chunk.n_frames):
            xyz = chunk.xyz[f] * 10.0  # nm → Å

            # Dynamic membrane normal this frame
            normal, midplane_z = _membrane_normal(xyz, p_local, upper_mask)

            for tm, loc in ca_local.items():
                pos = xyz[loc]
                raw_axis = _helix_axis(pos)
                tilt, axis = _tilt_vs_normal(raw_axis, normal)
                rotation = _helix_rotation_vs_normal(axis, normal)

                records.append(
                    {
                        "frame": frame_counter + f,
                        "helix": tm,
                        "tilt": tilt,
                        "rotation": rotation,
                        "cx": pos.mean(0)[0],
                        "cy": pos.mean(0)[1],
                        "cz": pos.mean(0)[2],
                        "ax": axis[0],
                        "ay": axis[1],
                        "az": axis[2],
                        "norm_x": normal[0],  # membrane normal this frame
                        "norm_y": normal[1],
                        "norm_z": normal[2],
                        "midplane_z": midplane_z,
                    }
                )

        frame_counter += chunk.n_frames

    geo_df = pd.DataFrame(records)
    geo_df.to_csv("helix_geometry.csv", index=False)
    print(
        f"Saved helix_geometry.csv ({geo_df.frame.nunique()} frames, "
        f"{geo_df.helix.nunique()} helices)"
    )
    return geo_df


# -----------------------------------------------------------------------------
# Write interpretations + data file
# -----------------------------------------------------------------------------

# Section explanations (static) -----------------------------------------------
_SECTION_TEXT = {
    "header": """FFAR1 / GPCR World-D Analysis — Interpretations + Time Series Data
Generated by analyze_worldD_poster.py
All numerical values rounded to 2 decimal places.
Columns in data tables separated by tabs.
This file is intended to be read by an AI for automatic interpretation.
""",
    "fig1": """Fig 1 - Force & torque decomposition
=====================================
What is plotted
  Three grouped bars per TM helix (mean +/- 1 SD over all frames):
    Lateral |F|  - force in the membrane plane (xy). Drives in-plane translation/tilt.
    Axial   |Fz| - force along membrane normal (z). Drives insertion/extraction.
    Torque  |t|  - rotational moment around the inboard bond. Drives helix spin.
  Large SD = dynamic fluctuating load. Small SD = steady persistent load.
  High lateral + low torque = helix pushed sideways by lipid/neighbour.
  High torque = twisting drive (relevant for TM5/TM6 activation rotation).
""",
    "fig2": """Fig 2 - Helix tilt and unwrapped azimuthal rotation
=====================================================
What is plotted
  Tilt     = angle between helix principal axis (Ca PCA) and membrane normal z (degrees).
  Rotation = azimuthal angle of principal axis projected onto xy, unwrapped (degrees).
  Increasing tilt over time = helix leaning away from normal (activation signature).
  TM5/TM6 activation typically shows +5 to +15 degrees tilt increase.
  Monotone rotation drift = sustained spin. Oscillation = restoring torque at preferred register.
  Large TM6 rotation is associated with DRY-motif exposure and G-protein coupling.
""",
    "fig3": """Fig 3 - Force / torque time series (per TM helix, per frame)
=============================================================
What is plotted
  |F| = total force magnitude (kcal/mol/Ang). |T| = torque magnitude (kcal/mol).
  Values are mean over boundary residues of each TM helix.
  Persistent high values = sustained structural constraint.
  Transient spikes = conformational events (lipid rearrangement, water ingress, transition).
  Correlated spikes across multiple helices = collective bundle rearrangement.
""",
    "fig4": """Fig 4 - Generalised coordinate u and velocity uDot (per TM helix, per frame)
=============================================================================
What is plotted
  u    = mean phi dihedral at helix boundaries (radians).
  uDot = generalised velocity (radians/ps).
  Persistent non-zero u = helix settled in a preferred dihedral state.
  High |uDot| = dynamically active pivot point.
  Rapid sign alternation in uDot = high-frequency local oscillation (not rigid-body motion).
""",
    "fig5": """Fig 5 - Cross-correlation |F| vs u per TM helix
================================================
What is plotted
  Normalised cross-correlation between |F(t)| and u(t). Lag in frames.
  Positive lag peak: |F| leads u  -> force is the cause, dihedral is the effect (force-driven motion).
  Negative lag peak: u  leads |F| -> geometry changes first, force builds after (strain loading).
  Zero lag peak: simultaneous coupling (underdamped oscillator).
  95% CI for white noise shown as reference.
""",
    "fig6": """Fig 6 - Mean boundary force vs helix tilt (scatter)
====================================================
What is plotted
  X = mean |F| at TM boundary residues. Y = helix tilt (degrees). One point per frame.
  Pearson r and p-value per TM helix.
  Positive r (p<0.05): force drives tilt (mechanical-geometric coupling).
  Negative r: high force accompanies smaller tilt (restoring/membrane resistance force).
  Slope = compliance: delta-tilt per unit force.
  n.s.: force and tilt decoupled at this timescale.
""",
    "fig7": """Fig 7 - Inter-helix centroid distances
=======================================
What is plotted
  Euclidean distance between Ca centroids of selected TM helix pairs (Angstrom).
  TM3-TM5, TM3-TM6, TM5-TM6: binding pocket triangle. Closure = agonist-stabilised state.
  TM1-TM7: extracellular cap, structural integrity marker.
  TM2-TM4: intracellular coupling face.
  TM6-TM7: NPxxY-motif region, activation gate.
  Decreasing TM5-TM6 = binding pocket closure.
  Increasing TM3-TM6 = TM6 outward swing (G-protein coupling geometry).
  Delta_mean > 1.5 Ang with reduced SD = significant ligand effect.
""",
    "fig8": """Fig 8 - Internal-coordinate RMSD (phi dihedral drift per TM helix)
===================================================================
What is plotted
  RMSD of u relative to frame 0, averaged over boundary residues per TM helix (radians).
  Differences wrapped to (-pi, pi] before squaring to handle circularity.
  Rising RMSD = progressive departure from starting geometry (slow conformational transition).
  Flat RMSD   = structurally locked helix.
  Plateau     = helix reached a new minimum (converged).
  Still-rising at end = insufficient sampling for that helix.
  Preferred over Cartesian RMSD here because it requires no alignment and isolates
  torsional DOF — exactly what Robosample samples in World D.
""",
    "stats_note": """Descriptive Statistics Table
=============================
Columns: mean, std, skew for each of F_mag, F_lat, F_axial, T_mag, u, uDot.
Values are aggregated (mean of boundary residues) per TM helix per frame, then statistics
computed over all frames.
""",
    "xcorr_note": """Cross-correlation Summary Table
================================
peak_lag: lag (frames) at which |cross-correlation| is maximised.
  Positive = |F| leads u (force causes dihedral change).
  Negative = u leads |F| (geometry change precedes force buildup).
peak_r: normalised cross-correlation value at peak_lag.
  |peak_r| > 1.96/sqrt(N) is significant at 95% level.
""",
}

# Data serialisation helpers ---------------------------------------------------


def _ts_block(df_pivot, fmt=".2f"):
    """Convert a frame x helix pivot to a tab-separated string."""
    lines = ["frame\t" + "\t".join(str(c) for c in df_pivot.columns)]
    for idx, row in df_pivot.iterrows():
        vals = "\t".join(format(v, fmt) if pd.notna(v) else "nan" for v in row)
        lines.append(f"{idx}\t{vals}")
    return "\n".join(lines)


def _stat_block(stats_df, fmt=".2f"):
    lines = ["helix\t" + "\t".join(stats_df.columns)]
    for idx, row in stats_df.iterrows():
        vals = "\t".join(format(v, fmt) for v in row)
        lines.append(f"{idx}\t{vals}")
    return "\n".join(lines)


def _xcorr_block(peak_lag, peak_r, fmt=".2f"):
    helices = list(peak_lag.keys())
    lines = [
        "metric\t" + "\t".join(helices),
        "peak_lag\t" + "\t".join(format(peak_lag[h], ".0f") for h in helices),
        "peak_r\t" + "\t".join(format(peak_r[h], fmt) for h in helices),
    ]
    return "\n".join(lines)


def _dist_block(geo_df, pairs, fmt=".2f"):
    """frame x pair distance table."""
    try:
        pivoted = geo_df.pivot(
            index="frame", columns="helix", values=["cx", "cy", "cz"]
        )
    except Exception:
        return "(distance data unavailable)"
    frames = sorted(geo_df.frame.unique())
    header = "frame\t" + "\t".join(f"{h1}-{h2}" for h1, h2 in pairs)
    rows = [header]
    for f in frames:
        vals = []
        for h1, h2 in pairs:
            try:
                c1 = pivoted.loc[f, (["cx", "cy", "cz"], h1)].values
                c2 = pivoted.loc[f, (["cx", "cy", "cz"], h2)].values
                vals.append(format(float(np.linalg.norm(c1 - c2)), fmt))
            except Exception:
                vals.append("nan")
        rows.append(f"{f}\t" + "\t".join(vals))
    return "\n".join(rows)


def _rmsd_block(df, fmt=".2f"):
    u0 = df[df.frame == df.frame.min()].set_index("helix")["u"]
    records = {}
    for tm, labels in BOUNDARY_MAP.items():
        sub = df[df.helix.isin(labels)]
        rmsd = sub.groupby("frame").apply(
            lambda g: float(
                np.sqrt(
                    (
                        (
                            (
                                g.set_index("helix")["u"]
                                - u0.reindex(g.helix.values).values
                                + np.pi
                            )
                            % (2 * np.pi)
                            - np.pi
                        )
                        ** 2
                    ).mean()
                )
            )
        )
        records[tm] = rmsd
    rmsd_df = pd.DataFrame(records)
    return _ts_block(rmsd_df, fmt=fmt)


# Main writer ------------------------------------------------------------------


def write_interpretations(
    path="interpretations.txt",
    tm_df=None,
    df=None,
    geo_df=None,
    stats_df=None,
    peak_lag=None,
    peak_r=None,
):
    """
    Write section explanations + actual numerical time series data to path.
    All numbers rounded to 2 decimal places.
    """
    S = _SECTION_TEXT
    lines = [S["header"]]

    # -- Fig 1 --
    lines += [S["fig1"]]
    if tm_df is not None:
        lines += ["DATA (mean per TM per frame, tab-separated):"]
        for col, label in [
            ("F_lat", "F_lat"),
            ("F_axial", "F_axial"),
            ("T_mag", "T_mag"),
        ]:
            piv = tm_df.pivot(index="frame", columns="helix", values=col)[TM_LIST]
            lines += [f"  {label}:", _ts_block(piv), ""]

    # -- Fig 2 --
    lines += [S["fig2"]]
    if geo_df is not None:
        for col, label in [
            ("tilt", "tilt_deg"),
            ("rotation", "rotation_deg_unwrapped"),
        ]:
            sub = geo_df.copy()
            if col == "rotation":
                sub[col] = sub.groupby("helix")["rotation"].transform(
                    lambda s: np.rad2deg(np.unwrap(np.deg2rad(s.values)))
                )
            piv = sub.pivot(index="frame", columns="helix", values=col)[TM_LIST]
            lines += [f"  {label}:", _ts_block(piv), ""]

    # -- Fig 3 --
    lines += [S["fig3"]]
    if tm_df is not None:
        for col, label in [("F_mag", "|F|"), ("T_mag", "|T|")]:
            piv = tm_df.pivot(index="frame", columns="helix", values=col)[TM_LIST]
            lines += [f"  {label}:", _ts_block(piv), ""]

    # -- Fig 4 --
    lines += [S["fig4"]]
    if tm_df is not None:
        for col, label in [("u", "u_rad"), ("uDot", "uDot_rad_per_ps")]:
            piv = tm_df.pivot(index="frame", columns="helix", values=col)[TM_LIST]
            lines += [f"  {label}:", _ts_block(piv), ""]

    # -- Fig 5 --
    lines += [S["fig5"]]
    if peak_lag is not None and peak_r is not None:
        lines += [S["xcorr_note"], _xcorr_block(peak_lag, peak_r), ""]

    # -- Fig 6 --
    lines += [S["fig6"]]
    if tm_df is not None and geo_df is not None:
        lines += ["  DATA (mean |F| per TM per frame):"]
        piv = tm_df.pivot(index="frame", columns="helix", values="F_mag")[TM_LIST]
        lines += [_ts_block(piv), ""]
        lines += ["  DATA (tilt per TM per frame):"]
        piv2 = geo_df.pivot(index="frame", columns="helix", values="tilt")[TM_LIST]
        lines += [_ts_block(piv2), ""]

    # -- Fig 7 --
    lines += [S["fig7"]]
    if geo_df is not None:
        lines += ["  DATA (centroid distances in Ang, tab-separated):"]
        lines += [_dist_block(geo_df, PAIRS), ""]

    # -- Fig 8 --
    lines += [S["fig8"]]
    if df is not None:
        lines += ["  DATA (internal-coordinate RMSD, radians):"]
        lines += [_rmsd_block(df), ""]

    # -- Statistics --
    lines += [S["stats_note"]]
    if stats_df is not None:
        lines += [_stat_block(stats_df.round(2)), ""]

    with open(path, "w") as f:
        f.write("\n".join(lines))
    print(f"Saved {path}")


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    vec_path = "vmd/vectors.dat"
    top_path = "examples/ffar1.prmtop"
    traj_path = "vmd/ffar1_6000.repl0.dcd"
    stride = int(sys.argv[4]) if len(sys.argv) > 4 else 1

    print(f"\nLoading trajectory ({top_path} + {traj_path}, stride={stride}) ...")
    geo_df = compute_helix_geometry(top_path, traj_path, stride=stride)

    print(f"Loading {vec_path} ...")
    df = load_vectors(vec_path)
    df = decompose_wrt_normal(df, geo_df)
    tm_df = aggregate_by_tm(df)
    print(f"  {df.frame.nunique()} frames, {df.helix.nunique()} boundary residues")

    stats_df = compute_statistics(tm_df)

    fig1_force_decomposition(df, tm_df)
    fig3_heatmap_combined(tm_df)
    fig5_xcorr_summary(tm_df)
    fig8_internal_rmsd(df)

    fig2_helix_tilt(geo_df)
    dist_stats = fig7_interhelical(geo_df)
    if dist_stats is not None:
        dist_stats.to_csv("interhelical_stats.csv", index=False)
    fig6_force_vs_tilt(tm_df, geo_df)

    # cross-correlation data needed by writer
    peak_lag_data, peak_r_data = {}, {}
    for tm in TM_LIST:
        sub = tm_df[tm_df.helix == tm].sort_values("frame")
        f_arr = sub.F_mag.values
        u_arr = sub.u.values
        f_z = f_arr - f_arr.mean()
        u_z = u_arr - u_arr.mean()
        n = len(f_z)
        from scipy import signal as _sig

        lags = _sig.correlation_lags(n, n, mode="full")
        xc = _sig.correlate(f_z, u_z, mode="full")
        denom = np.std(f_z) * np.std(u_z) * n
        xc = xc / denom if denom > 0 else xc
        pk = np.argmax(np.abs(xc))
        peak_lag_data[tm] = int(lags[pk])
        peak_r_data[tm] = float(xc[pk])

    write_interpretations(
        path="interpretations.txt",
        tm_df=tm_df,
        df=df,
        geo_df=geo_df,
        stats_df=stats_df,
        peak_lag=peak_lag_data,
        peak_r=peak_r_data,
    )
    print("\nDone. Figures: fig1-fig8 + interpretations.txt")
