#region Imports
import numpy as np
import mdtraj as md
import math
import matplotlib
import matplotlib.pyplot as plt
from batana import *
import seaborn as sns
import scipy
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import scipy.spatial.distance as ssd
import networkx as nx
#endregion # Imports

#region Functions
def share_mid_bond_dihedrals(dih1, dih2):
    """ Check if two dihedrals share the middle bond
    :param dih1: dihedral 1 indexes
    :param dih2: dihedral 2 indexes
    :return: True if share middle bond, False otherwise
    """
    # Extract middle bonds (b, c) from each dihedral
    bond1 = (dih1[1], dih1[2])
    bond2 = (dih2[1], dih2[2])

    # Compare bonds, ignoring direction
    return bond1 == bond2 or bond1 == bond2[::-1]
#
#endregion # Functions

#region Parse arguments 
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--prmtop', default=None, 
    help='Prmtop file.')
parser.add_argument('--inpcrd', default=None, 
    help='Inpcrd file.')
parser.add_argument('--dcd', default=None, 
    help='Trajectory file.')

parser.add_argument('--stride', default=1, type=int,
    help='Stride for the read lines.')
parser.add_argument('--analyze', default=[], nargs='+', 
    help='Decide what to analyze')
args = parser.parse_args()
#endregion # Parse arguments

# -----------------------------------------------------------------------------
#                            Main function
#region -----------------------------------------------------------------------
def main(prmtop, dcd):

    bat = BAT(dcd, prmtop)
    bat.calcBATIndexes()
    bat.calcBAT()

    nframes = bat.bos.shape[0]
    nofBos = bat.bos.shape[1]
    nofAngs = bat.angs.shape[1]
    nofDihs = bat.dihs.shape[1]
    print("Number of frames:", nframes)
    print("Bonds shape:", np.array(bat.bos).shape)
    print("Angles shape:", np.array(bat.angs).shape)
    print("Dihedrals shape:", np.array(bat.dihs).shape)
    
    dihs = bat.dihs
    nofDihs = dihs.shape[1]
    bats = BATStats()
    top = bat.mdtrajObj.topology

    nframes_check = min(2, nframes)
    nofBos_check = min(100, nofBos)
    nofAngs_check = min(100, nofAngs)
    nofDihs_check = min(100, nofDihs)

    #region Debug print
    # print("Bond indeces:", end=' ')
    # for boIx in range(nofBos_check):
    #     print(bat.boIxs[boIx], end=' ')
    # print()
    # print("Angle indeces:", end=' ')
    # for angIx in range(nofAngs_check):
    #     print(bat.angIxs[angIx], end=' ')
    # print()
    # print("Dihedral indeces:")
    # for dihIx in range(nofDihs_check):
    #     print(dihIx,":",  bat.dihIxs[dihIx][0], bat.dihIxs[dihIx][1], bat.dihIxs[dihIx][2], bat.dihIxs[dihIx][3], end=' ')
    #     for frameIx in range(nframes_check):
    #         print(bat.dihs[frameIx][dihIx], end=' ')
    #     print()
    # plt.figure()
    # plt.plot(range(bat.dihIxs.shape[0]), bat.dihIxs[:,1], label='Dihedral 1')
    # plt.show()
    # exit()
    # for frameIx in range(nframes_check):
    #     print("bonds:", end=' ')
    #     for boIx in range(nofBos_check):
    #         print(bat.bos[frameIx][boIx], end=' ')
    #     print()
    #     print("angles:", end=' ')
    #     for angIx in range(nofAngs_check):
    #         print(np.degrees(bat.angs[frameIx][angIx]), end=' ')
    #     print()
    #     print("dihedrals:", end=' ')
    #     for dihIx in range(nofDihs_check):
    #         print(np.degrees(bat.dihs[frameIx][dihIx]), end=' ')
    #     print()
    # exit(0)
    #endregion Debug print

    kept_indices = []
    seen_bonds = set()

    for i in range(nofDihs):
        # --- Filter 1: Terminal Hydrogen ---
        dih = bat.dihIxs[i]
        if top.atom(dih[0]).element.symbol == 'H' or top.atom(dih[3]).element.symbol == 'H':
            continue
            
        # --- Filter 2: One dihedral per unique middle bond ---
        # This prevents having 3 different dihedrals for the same rotating bond
        mid_bond = tuple(sorted((dih[1], dih[2])))
        if mid_bond in seen_bonds:
            continue
        
        # If it passes both, keep it
        kept_indices.append(i)
        seen_bonds.add(mid_bond)

    print(f"Reduced from {nofDihs} to {len(kept_indices)} representative dihedrals.")

    # Slice the main dihedrals array to get only the representative columns
    # dihs has shape (nframes, nofDihs)
    reduced_data = dihs[:, kept_indices] 
    n_reduced = len(kept_indices)

    # Compute All-vs-All Correlation for this subset
    reduced_corr_matrix = np.zeros((n_reduced, n_reduced))

    for i in range(n_reduced):
        for j in range(i, n_reduced):
            # Calculate correlation using your bats tool
            c = bats.dihedralsCorrelation(reduced_data[:, i], reduced_data[:, j])
            reduced_corr_matrix[i, j] = c
            reduced_corr_matrix[j, i] = c # Symmetric matrix

    # Create labels for the axes (e.g., "ALA15 C-N-CA-C")

    labels = []
    for idx in kept_indices:
        atoms = bat.dihIxs[idx]
        res = top.atom(atoms[0]).residue
        resName = res.name
        resSeq = res.resSeq + 1  # mdtraj is 0-based
        names = [top.atom(a).name for a in atoms]
        labels.append(f"{resName}{resSeq} {'-'.join(names)}")

    for i in range(n_reduced):
        print(f"{labels[i]}:", end=' ')
        for j in range(i, n_reduced):
            if i < j:
                if reduced_corr_matrix[i, j] > 0.8:
                    print(f"{labels[j]}: {reduced_corr_matrix[i, j]:.4f}", end=' ')
                    #print(f"{reduced_corr_matrix[i, j]:.4f}", end=' ')
        print()

    # 5. Plot
    plt.figure(figsize=(12, 10))
    sns.heatmap(reduced_corr_matrix, 
                xticklabels=labels, 
                yticklabels=labels, 
                cmap='viridis', 
                vmin=0, vmax=1)

    plt.title("Reduced All-vs-All Dihedral Correlation")
    plt.xticks(rotation=90, fontsize=8)
    plt.yticks(fontsize=8)
    plt.tight_layout()
    plt.show()
    
    # Clustering Threshold
    min_corr = 0.79
    dist_threshold = 1.0 - min_corr

    clustering_method = 'graph_based'  # Options: 'hierarchical' or 'graph_based'

    if clustering_method == 'graph_based':

        # 1. Create a graph object
        G = nx.Graph()

        # 2. Add nodes for every dihedral in your reduced set
        # We use the index relative to kept_indices (0 to n_reduced-1)
        for i in range(len(kept_indices)):
            G.add_node(i)

        # 3. Add edges between any two dihedrals that are highly correlated
        for i in range(len(kept_indices)):
            for j in range(i + 1, len(kept_indices)):
                if np.abs(reduced_corr_matrix[i, j]) > min_corr:
                    G.add_edge(i, j)

        # 4. Extract "Connected Components"
        # This is exactly what you asked for: if A-B and B-C, then (A, B, C)
        clusters = list(nx.connected_components(G))

        # Print the results
        print(f"--- Graph-Based Clusters (r > {min_corr}) ---")
        for count, cluster_nodes in enumerate(clusters):
            if len(cluster_nodes) > 1: # Skip dihedrals that aren't correlated with anything
                print(f"\nCluster {count + 1}:")
                for node_idx in cluster_nodes:

                    orig_idx = kept_indices[node_idx] # IMPORTANT: map back to original dihedral index
                    
                    atoms = bat.dihIxs[orig_idx]
                    res = top.atom(atoms[0]).residue
                    resName = res.name
                    resSeq = res.resSeq + 1  # mdtraj is 0-based
                    names = "-".join([top.atom(a).name for a in atoms])
                    print(f"  - [{orig_idx}] {resName}{resSeq} {names}", end=' ')
                    print(atoms[0], atoms[1], atoms[2], atoms[3])

        plt.figure(figsize=(10, 10))
        pos = nx.spring_layout(G) # Calculates positions for nodes
        nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=500, font_size=8)
        plt.title(f"Dihedral Correlation Network (r > {min_corr})")
        plt.show()

    pass

#endregion

#region Main
if __name__=="__main__":

    if not args.prmtop or not args.dcd:
        parser.error("Both --prmtop and --dcd are required.")   

    main(args.prmtop, args.dcd)
#endregion
