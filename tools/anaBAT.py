#region Imports
import numpy as np
import mdtraj as md
import math
import matplotlib
import matplotlib.pyplot as plt
from batana import *
#endregion

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
#endregion
# -----------------------------------------------------------------------------
#                            Main function
#region -----------------------------------------------------------------------
def main(prmtop, dcd):

    bat = BAT(dcd, prmtop)
    bat.calcBATIndexes()
    bat.calcBAT()

    dihs = bat.dihs
    nofFrames = dihs.shape[0]
    nofDihs = dihs.shape[1]

    bats = BATStats()

    # Get important dihedrals correlations
    for ix in range(nofDihs):
        for jx in range(nofDihs):
            if ix <= jx:
                corr = bats.dihedralsCorrelation(dihs[:, ix], dihs[:, jx])
                if corr > 0.8:
                    print(ix, jx, corr)
    exit()

    # Get all-vs-all dihedral correlations
    dihsCorrsAva = np.empty((nofDihs, nofDihs)) * np.nan
    for ix in range(nofDihs):
        for jx in range(nofDihs):
            if ix <= jx:
                dihsCorrsAva[ix, jx] = bats.dihedralsCorrelation(dihs[:, ix], dihs[:, jx])
            if jx > 10: break # COMMENT THIS
        if ix > 10: break # COMMENT THIS

    # Fill the opposite half
    for ix in range(nofDihs):
        for jx in range(nofDihs):
            if ix <= jx:
                dihsCorrsAva[jx, ix] = dihsCorrsAva[ix, jx]

    # Print dihedrals correlations
    for ix in range(nofDihs):
        for jx in range(nofDihs):
            print(dihsCorrsAva[ix, jx], end = ' ')
            if jx > 10: break # COMMENT THIS
        print()
        if ix > 10: break # COMMENT THIS

#endregion

#region Main
if __name__=="__main__":

    if not args.prmtop or not args.dcd:
        parser.error("Both --prmtop and --dcd are required.")   

    main(args.prmtop, args.dcd)
#endregion
