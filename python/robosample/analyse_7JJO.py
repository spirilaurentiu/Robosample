"""Agonist-bound vs apo reaction-force analysis for 7JJO (turkey β1-AR).

Recreates the eight-figure analyse_ffar1.py suite as a within-receptor bound-vs-apo
comparison. Reads 7JJO.{lig,noLig}.<replica>.reactions.csv + the paired DCDs.
Spec: docs/specs/gpcr-world-design/30-ligand-vs-apo-reaction-analysis.md

Run from the repo root:
    python3 python/robosample/analyse_7JJO.py [--replica 0] [--stride 1] [--outdir .]
"""

from betaAR_pair_common import run

if __name__ == "__main__":
    run("7JJO")
