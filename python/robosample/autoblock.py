import argparse
import warnings

import matplotlib.pyplot as plt
import MDAnalysis as mda
import mdtraj as md
import numpy as np
import parmed as pmd
from deeptime.decomposition import TICA
from deeptime.markov import TransitionCountEstimator
from deeptime.markov.msm import MaximumLikelihoodMSM
from MDAnalysis.analysis import dihedrals as mda_dihedrals
from scipy.ndimage import label as nd_label
from scipy.stats import entropy as sp_entropy
from sklearn.cluster import MiniBatchKMeans

import robosample

# python3 python/robosample/autoblock.py 1apq examples/1APQ.prmtop examples/1APQ.rst7 6000 0 20000 1

# Create the parser
parser = argparse.ArgumentParser(description="Process PDB code and seed.")

# Add the arguments
parser.add_argument("name", type=str, help="Name of the simulation.")
parser.add_argument("prmtop", type=str, help="Relative path to the .prmtop file.")
parser.add_argument("inpcrd", type=str, help="Relative path to the .inpcrd file.")
parser.add_argument("seed", type=int, help="The seed.")
parser.add_argument("equil_steps", type=int, help="The number of equilibration steps.")
parser.add_argument("prod_steps", type=int, help="The number of production steps.")
parser.add_argument("write_freq", type=int, help="CSV and DCD write frequency.")

# Parse the arguments
args = parser.parse_args()

# Mean first passage time to cross an energy barrier
# 3 kcal/mol - sub-picosecond to picosecond transitions (modest barrier, ~5KbT)
# 6 kcal/mol - tens to hundreds of picoseconds (moderate barrier, ~10KbT)
# 10 kcal/mol - nanoseconds or longer (high barrier , ~16KbT)
# 2 ps of MD is enough to explore shallow wells, but not to cross deep barriers without enhanced sampling (e.g., HMC, replica exchange)
TIMESTEP_CARTESIAN = 0.001
MDSTEPS_CARTESIAN = 2000

# create robosample context
context = robosample.Context(
    name=args.name,
    seed=args.seed,
    prmtop=args.prmtop,
    inpcrd=args.inpcrd,
    write_freq=args.write_freq,
    testing=False,
)

# # Add cartesian world (will integrate with OpenMM)
# context.addCartesianWorld().add_sampler(
#     timeStep=TIMESTEP_CARTESIAN,
#     mdSteps=MDSTEPS_CARTESIAN,
#     boostMDSteps=MDSTEPS_CARTESIAN,
#     acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
#     use_nuts=False,
# )

# # Add replicas (geometric temperature ladder)
# context.initialize([300])
# context.run_rex(args.equil_steps, args.prod_steps, args.write_freq, True, False)

TIMESTEP_TD = 0.02
MDSTEPS_TD = 1024
dcd_file = "1apq_6000.repl0.20ksteps.dcd"

# windows = context.detect_transition_windows(dcd_file=dcd_file)
# matrix = context.compute_differential_correlation(dcd_file=dcd_file, result=windows)

# plt.imshow(matrix, cmap="viridis")
# plt.colorbar()
# plt.title("Differential Correlation Matrix")
# plt.xlabel("Dihedral Index")
# plt.ylabel("Dihedral Index")
# plt.show()


#######################################################3


#  1. Dihedral extraction


def extract_dihedrals(universe, dihedral_atom_groups):
    """
    Returns
    -------
    angles_rad : (T, N)   raw angles in radians
    X          : (T, 2N)  sin/cos torus embedding
    """
    ags = [
        mda.AtomGroup(universe.atoms[[gp, p, c, gc]])
        for gp, p, c, gc in dihedral_atom_groups
    ]
    angles_deg = mda_dihedrals.Dihedral(ags).run().angles  # (T, N)
    angles_rad = np.deg2rad(angles_deg).astype(np.float64)
    X = np.concatenate([np.sin(angles_rad), np.cos(angles_rad)], axis=1)
    return angles_rad, X


#  2. TICA


def run_tica(X, lagtime, n_components=10):
    """
    lagtime heuristic: ~10 % of the slowest process you care about.
    For loop rearrangements at 20k frames with 5-10 transitions,
    start with lagtime = 50-200 frames and scan.
    """
    tica = TICA(lagtime=lagtime, dim=n_components)
    Y = tica.fit_transform(X)  # deeptime returns ndarray directly
    return Y, tica


#  3. Micro-state clustering


def cluster_tica(Y, n_clusters=200, random_state=42):
    """
    KMeans is required here (not HDBSCAN) because:
      • every frame must receive a label  no -1 noise points
      • the resulting dtraj must be contiguous for MSM lag counting
      100-300 microstates is the standard range for MSM estimation
    """
    km = MiniBatchKMeans(n_clusters=n_clusters, n_init=5, random_state=random_state)
    dtraj = km.fit_predict(Y).astype(int)
    return dtraj, km


#  4. MSM with connected-set restriction


def build_msm(dtraj, lagtime):
    counts = (
        TransitionCountEstimator(lagtime=lagtime, count_mode="effective")
        .fit(dtraj)
        .fetch_model()
        .submodel_largest(directed=False)
    )
    msm = MaximumLikelihoodMSM(reversible=True).fit(counts).fetch_model()

    #  CRITICAL: use the MSM's own state list, not counts.states
    # MaximumLikelihoodMSM(reversible=True) may prune further for detailed
    # balance; msm.count_model.states reflects the final active set.
    msm_states = msm.count_model.states  # shape (n_msm_states,)
    state_map = {s: i for i, s in enumerate(msm_states)}
    active_dtraj = np.array([state_map.get(s, -1) for s in dtraj])

    print(f"  submodel_largest : {counts.n_states} microstates")
    print(f"  after MSM fit    : {msm.n_states} microstates")
    print(f"  frames pruned    : {(active_dtraj == -1).sum()}")

    return msm, counts, active_dtraj


#  5. Spectral gap -> number of metastable states


def find_n_metastable(msm, max_n=20, plot=True):
    """
    The eigenvalue spectrum has a large gap between eigenvalues k and k+1
    when there are k metastable basins.  We scan |Δλ_i| for i ≥ 1
    (skip the trivial i=0 gap from λ₀=1).
    """
    eigs = np.sort(np.real(msm.eigenvalues()))[::-1]
    eigs = eigs[: min(max_n + 2, len(eigs))]

    gaps = np.abs(np.diff(eigs))  # all positive (eigenvalues descend)
    # gaps[0] is the drop from λ₀≈1 -> always large, skip it
    best_idx = int(np.argmax(gaps[1:])) + 1  # index in gaps[]
    n_metastable = best_idx + 1  # #eigenvalues above the gap

    if plot:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        x = np.arange(len(eigs))
        axes[0].plot(x, eigs, "o-")
        axes[0].axvline(
            best_idx, color="red", ls="--", label=f"gap -> {n_metastable} states"
        )
        axes[0].set_xlabel("index")
        axes[0].set_ylabel("λ")
        axes[0].set_title("MSM eigenvalue spectrum")
        axes[0].legend()

        axes[1].bar(range(len(gaps)), gaps)
        axes[1].axvline(
            best_idx, color="red", ls="--", label=f"largest gap at {best_idx}"
        )
        axes[1].set_xlabel("gap index")
        axes[1].set_ylabel("|Δλ|")
        axes[1].set_title("Spectral gaps")
        axes[1].legend()
        plt.tight_layout()
        plt.show()

    print(f"[spectral gap] n_metastable = {n_metastable}")
    print(f"  eigenvalues: {np.round(eigs[: n_metastable + 2], 4)}")
    return n_metastable, eigs, gaps


#  6. PCCA+ metastable decomposition


def pcca_coarse_grain(msm, active_dtraj, n_metastable):
    """
    deeptime API: msm.pcca(n) -> PCCAModel
      pcca.memberships  : (n_active_microstates, n_metastable)
      pcca.assignments  : (n_active_microstates,)  hard argmax

    We map every frame through active_dtraj to get per-frame memberships.
    Frames outside the active set (active_dtraj == -1) get NaN rows.
    """
    pcca = msm.pcca(n_metastable)

    T = len(active_dtraj)
    memberships = np.full((T, n_metastable), np.nan)
    hard_traj = np.full(T, -1, dtype=int)

    active_mask = active_dtraj >= 0
    memberships[active_mask] = pcca.memberships[active_dtraj[active_mask]]
    hard_traj[active_mask] = np.argmax(memberships[active_mask], axis=1)

    return memberships, hard_traj, pcca


#  7. Transition detection via membership entropy


def detect_transitions(memberships, hard_traj, entropy_threshold=None, flank_window=5):
    """
    Core idea: a frame is "in transition" when its PCCA membership vector
    is spread across multiple states (high Shannon entropy).
    This is more robust than watching for hard-state changes, which flicker.

    entropy_threshold : normalised entropy in [0,1] above which a frame is
        labelled transitioning.  None -> auto as mean + 1.5 sigma of the trace.

    Returns
    -------
    events      : list of dicts  {start, end, from_state, to_state,
                                   frames, peak_entropy}
    norm_ent    : (T,) normalised entropy trace  (0 = in-basin, 1 = uniform)
    """
    n_states = memberships.shape[1]
    max_ent = np.log(n_states)

    # Fill inactive frames with uniform -> high entropy (won't mask real transitions)
    mem = memberships.copy()
    nan_rows = np.any(np.isnan(mem), axis=1)
    mem[nan_rows] = 1.0 / n_states
    mem = np.clip(mem, 1e-12, 1.0)
    mem /= mem.sum(axis=1, keepdims=True)

    raw_ent = sp_entropy(mem.T)  # (T,)   scipy expects (k, n)
    norm_ent = raw_ent / max_ent  # [0,1]

    if entropy_threshold is None:
        mu, sigma = norm_ent.mean(), norm_ent.std()
        entropy_threshold = min(mu + 1.5 * sigma, 0.85)
        print(
            f"[auto threshold] entropy_threshold = {entropy_threshold:.3f}  "
            f"(μ={mu:.3f}, sigma={sigma:.3f})"
        )

    trans_mask = norm_ent > entropy_threshold
    labeled, n_blobs = nd_label(trans_mask)  # contiguous high-entropy blobs

    events = []
    for blob_id in range(1, n_blobs + 1):
        blob_frames = np.where(labeled == blob_id)[0]
        t_start, t_end = int(blob_frames[0]), int(blob_frames[-1])

        # Vote for flanking states (majority of flank_window frames each side)
        pre_states = [
            hard_traj[f]
            for f in range(max(0, t_start - flank_window), t_start)
            if hard_traj[f] >= 0
        ]
        post_states = [
            hard_traj[f]
            for f in range(t_end + 1, min(len(hard_traj), t_end + flank_window + 1))
            if hard_traj[f] >= 0
        ]

        if not pre_states or not post_states:
            continue

        state_before = int(np.bincount(pre_states).argmax())
        state_after = int(np.bincount(post_states).argmax())

        if state_before == state_after:
            continue  # entropy spike but no net state change -> skip

        events.append(
            {
                "start": t_start,
                "end": t_end,
                "from_state": state_before,
                "to_state": state_after,
                "frames": blob_frames.tolist(),
                "peak_entropy": float(norm_ent[blob_frames].max()),
            }
        )

    print(f"[transitions] {len(events)} events detected")
    for ev in events:
        print(
            f"  frames {ev['start']:>6}-{ev['end']:>6}  "
            f"state {ev['from_state']} -> {ev['to_state']}  "
            f"peak_H = {ev['peak_entropy']:.3f}"
        )

    return events, norm_ent


#  8. Transition-only correlation matrix


def transition_correlation_matrix(angles_rad, events, window=10):
    """
    Compute NxN circular dihedral correlation restricted to transition frames.

    angles_rad : (T, N) raw angles in radians  ← use this, not sin/cos,
                 so the correlation matrix is (NxN) not (2Nx2N)
    window     : extra frames to include on each side of each transition blob

    Circular correlation via mean-centred sin features is a reasonable
    approximation; for a full von Mises implementation swap in pingouin or
    pycircstat if needed.
    """
    T, N = angles_rad.shape

    frame_set = set()
    for ev in events:
        for f in range(ev["start"] - window, ev["end"] + window + 1):
            if 0 <= f < T:
                frame_set.add(f)

    if not frame_set:
        warnings.warn("No transition frames  falling back to global matrix.")
        frame_idx = np.arange(T)
    else:
        frame_idx = np.array(sorted(frame_set))

    pct = 100.0 * len(frame_idx) / T
    print(f"[correlation] {len(frame_idx)} / {T} frames ({pct:.1f} %)")

    # Circular correlation: centre the sin-projection
    sel = angles_rad[frame_idx]  # (F, N)
    s = np.sin(sel)
    s -= s.mean(axis=0)
    corr = np.corrcoef(s.T)  # (N, N)
    return corr, frame_idx


#  8b. Per-transition correlation matrices


def plot_per_transition_corr(angles_rad, events, window=10):
    """
    Compute and plot a separate NxN circular correlation matrix
    for each detected transition window.
    """
    T, N = angles_rad.shape

    if not events:
        print("[per-transition] no events to plot.")
        return []

    # n_events = len(events)

    # # Lay out subplots in a grid  at most 3 columns
    # ncols = min(n_events, 3)
    # nrows = (n_events + ncols - 1) // ncols
    # fig, axes = plt.subplots(
    #     nrows,
    #     ncols,
    #     figsize=(5 * ncols, 4.5 * nrows),
    #     squeeze=False,
    # )

    corr_matrices = []
    for idx, ev in enumerate(events):
        frame_set = set()
        for f in range(ev["start"] - window, ev["end"] + window + 1):
            if 0 <= f < T:
                frame_set.add(f)
        frame_idx = np.array(sorted(frame_set))

        sel = angles_rad[frame_idx]  # (F, N)
        s = np.sin(sel)
        s -= s.mean(axis=0)
        corr = np.corrcoef(s.T)  # (N, N)
        corr_matrices.append(corr)

    #     ax = axes[idx // ncols][idx % ncols]
    #     im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
    #     plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    #     ax.set_title(
    #         f"Transition {idx + 1}:  {ev['from_state']} -> {ev['to_state']}\n"
    #         f"frames {ev['start']}-{ev['end']}  "
    #         f"({len(frame_idx)} frames,  peak H={ev['peak_entropy']:.2f})",
    #         fontsize=8,
    #     )
    #     ax.set_xlabel("Dihedral index", fontsize=7)
    #     ax.set_ylabel("Dihedral index", fontsize=7)
    #     ax.tick_params(labelsize=6)

    # # Hide any leftover empty axes
    # for empty in range(n_events, nrows * ncols):
    #     axes[empty // ncols][empty % ncols].set_visible(False)

    # fig.suptitle("Per-transition dihedral correlation matrices", fontsize=11, y=1.01)
    # plt.tight_layout()
    # plt.show()

    return corr_matrices


#  9. Plots


def plot_entropy_trace(norm_ent, events, hard_traj):
    fig, axes = plt.subplots(2, 1, figsize=(14, 6), sharex=True)

    axes[0].plot(norm_ent, lw=0.5, color="steelblue")
    for ev in events:
        axes[0].axvspan(
            ev["start"],
            ev["end"],
            alpha=0.3,
            color="tomato",
            label=f"{ev['from_state']}->{ev['to_state']}",
        )
    axes[0].set_ylabel("Normalised entropy")
    axes[0].set_title("PCCA membership entropy  (red bands = detected transitions)")

    axes[1].plot(hard_traj, lw=0.4, color="darkorange")
    for ev in events:
        axes[1].axvspan(ev["start"], ev["end"], alpha=0.3, color="tomato")
    axes[1].set_ylabel("Metastable state")
    axes[1].set_xlabel("Frame")
    axes[1].set_title("Hard metastable state assignment")
    plt.tight_layout()
    plt.show()


def plot_corr(corr, title="Dihedral correlation  transition frames only"):
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
    plt.colorbar(im, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Dihedral index")
    ax.set_ylabel("Dihedral index")
    plt.tight_layout()
    plt.show()


#  10. Full pipeline


def run_pipeline(
    universe,
    dihedral_atom_groups,
    tica_lagtime=50,  # tune: 10-200 depending on trajectory length
    tica_dim=10,
    n_microstates=200,  # tune: sqrt(T) is a rough heuristic
    msm_lagtime=50,  # should match or be multiple of tica_lagtime
    max_metastable=15,
    entropy_threshold=None,  # None = auto
    window=10,
):
    print("=== Step 1: dihedrals ===")
    angles_rad, X = extract_dihedrals(universe, dihedral_atom_groups)
    T, N = angles_rad.shape
    print(f"  {T} frames x {N} dihedrals")

    print("=== Step 2: TICA ===")
    Y, tica = run_tica(X, lagtime=tica_lagtime, n_components=tica_dim)

    print("=== Step 3: KMeans microstates ===")
    dtraj, km = cluster_tica(Y, n_clusters=n_microstates)

    print("=== Step 4: MSM ===")
    msm, counts, active_dtraj = build_msm(dtraj, lagtime=msm_lagtime)
    print(f"  {counts.n_states} microstates in largest connected set")

    print("=== Step 5: spectral gap ===")
    n_meta, eigs, gaps = find_n_metastable(msm, max_n=max_metastable, plot=False)

    print("=== Step 6: PCCA+ ===")
    memberships, hard_traj, pcca = pcca_coarse_grain(msm, active_dtraj, n_meta)

    print("=== Step 7: transition detection ===")
    events, norm_ent = detect_transitions(
        memberships,
        hard_traj,
        entropy_threshold=entropy_threshold,
        flank_window=window,
    )

    print("=== Step 8: transition correlation matrix ===")
    corr, trans_frames = transition_correlation_matrix(
        angles_rad, events, window=window
    )
    return [corr]

    # plot_entropy_trace(norm_ent, events, hard_traj)
    # plot_corr(corr)  # aggregate (unchanged)

    # per_corrs = plot_per_transition_corr(  # ← new
    #     angles_rad, events, window=window
    # )
    # return per_corrs

    # return dict(
    #     angles_rad=angles_rad,
    #     X=X,
    #     Y=Y,
    #     dtraj=dtraj,
    #     msm=msm,
    #     n_metastable=n_meta,
    #     memberships=memberships,
    #     hard_traj=hard_traj,
    #     events=events,
    #     norm_ent=norm_ent,
    #     corr=corr,
    #     trans_frames=trans_frames,
    #     per_transition_corr=per_corrs,  # ← new
    # )


universe = mda.Universe(args.prmtop, dcd_file)
# run_pipeline(universe, context.standard_dihedral_atom_groups)


def scale_to_range(values, n):
    if not values:
        return []

    vmin = min(values)
    vmax = max(values)

    # avoid division by zero if all values are equal
    if vmax == vmin:
        return [1 for _ in values]

    scaled = [1 + (v - vmin) * (n - 1) / (vmax - vmin) for v in values]

    return [int(round(x)) for x in scaled]


#######################################################3

for cycle in range(10):
    # (
    #     (strong_blocks, strong_blocks_correlation),
    #     (weak_blocks, weak_blocks_correlation),
    #     rogue_blocks,
    # ) = context.build_gibbs_blocks_from_trajectory(dcd_file)

    # all_blocks = []
    # min_correlation = (
    #     min(strong_blocks_correlation + weak_blocks_correlation)
    #     if strong_blocks_correlation + weak_blocks_correlation
    #     else 0
    # )

    # expansion_factor = 1
    # for block, corr in zip(strong_blocks, strong_blocks_correlation):
    #     num_times = round((corr / min_correlation) * expansion_factor)
    #     all_blocks.extend([block] * num_times)

    # randomize = True
    # if randomize:
    #     random.shuffle(all_blocks)

    # Continue from last frame of the previous simulation
    traj = md.load(dcd_file, top=args.prmtop)
    last = traj[-1]

    parm = pmd.load_file(args.prmtop)
    parm.coordinates = last.xyz[0] * 10.0  # nm to Angstrom
    if last.unitcell_lengths is not None:
        parm.box = list(last.unitcell_lengths[0] * 10.0) + list(last.unitcell_angles[0])
    parm.save("last_frame.rst7", format="rst7", overwrite=True)

    # Create a new context for the next cycle of sampling
    context = robosample.Context(
        name=args.name + f"_cycle{cycle}",
        seed=args.seed,
        prmtop=args.prmtop,
        inpcrd=args.inpcrd,
        write_freq=args.write_freq,
        testing=False,
    )

    # Add cartesian world (will integrate with OpenMM)
    context.addCartesianWorld().add_sampler(
        timeStep=0.001,
        mdSteps=500,
        boostMDSteps=500,
        acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
        use_nuts=False,
    )

    correlation_matrices = run_pipeline(universe, context.standard_dihedral_atom_groups)
    for matrix in correlation_matrices:
        strong_blocks, weak_blocks, rogue_blocks, modularity = (
            context.chose_correlated_bonds(matrix)
        )
        blocks = strong_blocks + weak_blocks

        intra_block_corr = [np.mean(matrix[np.ix_(block, block)]) for block in blocks]
        min_intra_block_corr = min(intra_block_corr)

        print("=== Gibbs Blocks (high intra-correlation) ===")
        for i, block in enumerate(strong_blocks):
            mean_corr = np.mean(matrix[np.ix_(block, block)])
            print(
                f"  Block {i}: {len(block)} dihedrals | "
                f"mean intra-corr: {mean_corr:.3f} | "
                f"modularity contrib: {modularity.get(i, float('nan')):.4f}"
            )
        print("=== Weakly Correlated Blocks ===")
        for i, block in enumerate(weak_blocks):
            mean_corr = np.mean(matrix[np.ix_(block, block)])
            print(
                f"  Block {i}: {len(block)} dihedrals | "
                f"mean intra-corr: {mean_corr:.3f} | "
                f"modularity contrib: {modularity.get(i, float('nan')):.4f}"
            )
        print("=== Rogue Dihedrals (low correlation) ===")
        print(f"  {len(rogue_blocks)} dihedrals")

        num_oversamples = scale_to_range(intra_block_corr, 2)

        expanded_blocks = []
        for block, num in zip(blocks, num_oversamples):
            expanded_blocks.extend([block] * num)

        expanded_blocks.append(np.array([rogue_blocks]))

        # random.shuffle(expanded_blocks)
        for block in expanded_blocks:
            # block = np.array([rogue_blocks])
            block = block.reshape(-1)
            bonds = context.standard_dihedral_bonds.loc[block]

            sele = context.build_flexibilities(bonds)
            context.addTorsionalWorld(sele).add_sampler(
                timeStep=0.025,
                mdSteps=2000,  # ignored if using NUTS
                boostMDSteps=2000,  # ignored if using NUTS
                acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
                use_nuts=True,
            )

            # break
        break

    context.initialize([300])

    import time

    start_time = time.time()
    context.run_rex(0, 100, args.write_freq, True)
    end_time = time.time()
    duration = end_time - start_time
    print(f"Rex run time: {duration} seconds")

    # import openmm as mm
    # from openmm import app, unit

    # prmtop = app.AmberPrmtopFile(args.prmtop)
    # inpcrd = app.AmberInpcrdFile(args.inpcrd)
    # system = prmtop.createSystem(
    #     nonbondedMethod=app.CutoffNonPeriodic,
    #     nonbondedCutoff=1.2 * unit.nanometer,
    #     implicitSolvent=app.OBC2,
    #     constraints=app.HBonds,
    # )

    # # Integrator
    # timestep = 0.002 * unit.picoseconds  # 2 fs
    # integrator = mm.LangevinIntegrator(
    #     300 * unit.kelvin, 1.0 / unit.picosecond, timestep
    # )

    # # Platform (optional)
    # platform = mm.Platform.getPlatformByName("CUDA")

    # # Simulation
    # simulation = app.Simulation(prmtop.topology, system, integrator, platform)
    # simulation.context.setPositions(inpcrd.positions)

    # if inpcrd.velocities is not None:
    #     simulation.context.setVelocities(inpcrd.velocities)
    # else:
    #     simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

    # # Reporters
    # steps_per_frame = 1000
    # output_dcd = f"{args.name}_cycle{cycle}_openmm.dcd"
    # simulation.reporters.append(app.DCDReporter(output_dcd, steps_per_frame))

    # simulation.runForClockTime(time=duration * unit.second)

    break
