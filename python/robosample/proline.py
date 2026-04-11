import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from MDAnalysis.lib.distances import calc_dihedrals

# ========= USER INPUT =========
topology = "/home/victor/AllFrames/ensembleCluster1All.pdb"
# ==============================

u = mda.Universe(topology)


def cremer_pople_5ring(coords):
    """
    Compute Cremer-Pople puckering parameters (q2, q3) for 5-membered ring.
    coords: (5,3) array ordered sequentially around ring.
    """
    # center
    coords = coords - coords.mean(axis=0)

    # normal vector via SVD
    _, _, vh = np.linalg.svd(coords)
    normal = vh[2]

    # project onto normal
    z = np.dot(coords, normal)

    N = 5
    q2 = np.sqrt(2 / 5) * np.sum(z[k] * np.cos(4 * np.pi * k / N) for k in range(N))
    q3 = np.sqrt(2 / 5) * np.sum(z[k] * np.sin(4 * np.pi * k / N) for k in range(N))

    return q2, q3


prolines = u.select_atoms("resname PRO").residues

results = {}

for res in prolines:
    res_id = res.resid
    results[res_id] = {"omega": [], "q2": [], "q3": []}

    # ring atom order: N-CA-CB-CG-CD
    ring_atoms = res.atoms.select_atoms("name N CA CB CG CD")

    # omega atoms: C(i-1), N(i), CA(i), C(i)
    try:
        prev_res = res.universe.residues[res.ix - 1]
        omega_atoms = u.select_atoms(
            f"resid {prev_res.resid} and name C or (resid {res.resid} and name N CA C)"
        )
        if len(omega_atoms) != 4:
            continue
    except Exception:
        continue

    for ts in u.trajectory:
        # omega
        omega = calc_dihedrals(
            omega_atoms.positions[0:1],
            omega_atoms.positions[1:2],
            omega_atoms.positions[2:3],
            omega_atoms.positions[3:4],
        )[0]
        results[res_id]["omega"].append(np.degrees(omega))

        # puckering
        ring_coords = ring_atoms.positions
        if ring_coords.shape[0] == 5:
            q2, q3 = cremer_pople_5ring(ring_coords)
            results[res_id]["q2"].append(q2)
            results[res_id]["q3"].append(q3)

# ===== Plotting =====


def wrap_angle_deg(angle):
    return (angle + 180) % 360 - 180


for resid, data in results.items():
    omega = wrap_angle_deg(np.array(data["omega"]))
    q2 = np.array(data["q2"])
    q3 = np.array(data["q3"])

    cis_mask = np.abs(omega) < 30
    trans_mask = np.abs(np.abs(omega) - 180) < 30

    Q = np.sqrt(q2**2 + q3**2)
    phi = wrap_angle_deg(np.degrees(np.arctan2(q3, q2)))

    fig = plt.figure(figsize=(10, 10))
    gs = fig.add_gridspec(4, 1)

    # compute fractions
    trans_mask = np.abs(np.abs(omega) - 180) < 30
    cis_mask = np.abs(omega) < 30
    transition_mask = ~(cis_mask | trans_mask)

    cis_frac = np.mean(cis_mask) * 100
    trans_frac = np.mean(trans_mask) * 100
    transit_frac = np.mean(transition_mask) * 100

    # ===============================
    # 1) Omega with legend for cis/trans/transition
    # ===============================
    ax0 = fig.add_subplot(gs[0])

    # shaded cis/trans regions
    ax0.axhspan(-30, 30, color="red", alpha=0.15)
    ax0.axhspan(150, 180, color="blue", alpha=0.15)
    ax0.axhspan(-180, -150, color="blue", alpha=0.15)

    # omega trace
    ax0.plot(omega, color="black", linewidth=1)

    ax0.set_ylabel("Degrees")
    ax0.set_ylim(-180, 180)

    # add fractions to title
    ax0.set_title(
        f"Omega angles (Cis: {cis_frac:.2f}, Trans: {trans_frac:.2f}, Transition: {transit_frac:.2f})"
    )

    # create legend for the colored regions
    legend_elements = [
        Patch(facecolor="red", alpha=0.15, label="Cis"),
        Patch(facecolor="blue", alpha=0.15, label="Trans"),
        Patch(facecolor="none", edgecolor="none", label="Transition"),
    ]

    ax0.legend(handles=legend_elements, loc="upper right", frameon=True)

    # ===============================
    # 2) Puckering amplitude Q with legend
    # ===============================
    ax1 = fig.add_subplot(gs[1])
    # ax1.plot(Q, color="black", linewidth=1)
    ax1.scatter(range(len(Q)), Q, color="black", s=1)

    # Shaded bands
    ax1.axhspan(0, 0.05, color="gray", alpha=0.15)
    ax1.axhspan(0.05, 0.15, color="green", alpha=0.15)
    ax1.axhspan(0.15, max(Q) * 1.05, color="orange", alpha=0.15)

    ax1.set_ylabel("Å")
    ax1.set_title("Puckering amplitude Q")

    # Create legend entries for interpretation
    legend_elements = [
        Line2D([0], [0], color="orange", lw=4, alpha=0.5, label="Strong puckering"),
        Line2D([0], [0], color="green", lw=4, alpha=0.5, label="Moderate puckering"),
        Line2D([0], [0], color="gray", lw=4, alpha=0.5, label="Planar"),
    ]
    ax1.legend(handles=legend_elements, loc="upper right", frameon=True)

    # ===============================
    # 3) Pseudorotation phase phi
    # ===============================
    ax2 = fig.add_subplot(gs[2])

    # Shade Cγ-endo / Cγ-exo regions
    ax2.axhspan(-90, 0, color="purple", alpha=0.15)
    ax2.axhspan(90, 180, color="orange", alpha=0.15)
    ax2.axhspan(-180, -90, color="orange", alpha=0.15)

    # plot phi
    # ax2.plot(phi, color="black", linewidth=1)
    ax2.scatter(range(len(phi)), phi, color="black", s=1)

    ax2.set_ylabel("Degrees")
    ax2.set_ylim(-180, 180)
    ax2.set_title("Pseudorotation phase φ")

    # Add legend for interpretation
    legend_elements = [
        Patch(facecolor="purple", alpha=0.15, label="Cγ-endo region"),
        Patch(facecolor="orange", alpha=0.15, label="Cγ-exo region"),
    ]
    ax2.legend(handles=legend_elements, loc="upper right", frameon=True)

    # ===============================
    # 4) Polar representation (Cis/Trans + Cremer–Pople)
    # ===============================
    ax3 = fig.add_subplot(gs[3], projection="polar")

    # Map omega to colors
    colors = np.full_like(omega, fill_value="gray", dtype=object)
    colors[np.abs(omega) < 30] = "red"  # cis
    colors[np.abs(np.abs(omega) - 180) < 30] = "blue"  # trans
    # transition (gray) is default

    # Scatter in polar coordinates: φ = angle, Q = radius
    ax3.scatter(np.radians(phi), Q, s=1, c=colors, alpha=0.7)

    # Optional: highlight Cremer–Pople regions
    # Cγ-endo ~ φ -90° → 0°, Cγ-exo ~ φ 90° → 180° / -180° → -90° (shaded sectors)
    ax3.fill_between(
        np.radians(np.linspace(-90, 0, 100)), 0, max(Q) * 1.1, color="purple", alpha=0.1
    )
    ax3.fill_between(
        np.radians(
            np.concatenate([np.linspace(90, 180, 100), np.linspace(-180, -90, 100)])
        ),
        0,
        max(Q) * 1.1,
        color="orange",
        alpha=0.1,
    )

    ax3.set_title(
        "Polar view: φ = pseudorotation angle, Q = puckering amplitude\nColor = ω (cis/trans)"
    )

    # Build legend
    legend_elements = [
        Patch(facecolor="red", label="Cis ω (<30°)", alpha=0.5),
        Patch(facecolor="blue", label="Trans ω (~180°)", alpha=0.5),
        Patch(facecolor="gray", label="Transition ω", alpha=0.5),
        Patch(facecolor="purple", alpha=0.1, label="Cγ-endo φ region"),
        Patch(facecolor="orange", alpha=0.1, label="Cγ-exo φ region"),
    ]
    ax3.legend(
        handles=legend_elements,
        loc="upper left",  # location is relative to bbox
        bbox_to_anchor=(1.1, 1.1),  # just outside axes on the right
        borderaxespad=0.0,  # minimal padding
        frameon=True,
    )

    fig.suptitle(f"Proline {resid}")
    plt.tight_layout()
    plt.show()
