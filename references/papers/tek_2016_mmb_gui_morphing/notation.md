# Notation — Tek et al. 2016 (MMB-GUI morphing)

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| A | initial crystal structure of a morph | atomic coordinate set | endpoint of trajectory |
| B | experimentally observed intermediate structure | atomic coordinate set | target of recapitulation (not used to drive the morph) |
| C | final crystal structure of a morph | atomic coordinate set | endpoint of trajectory |
| i | interpolated (morphed) structure along the trajectory | atomic coordinate set (per frame) | index over morph frames |
| rmsd(XY) | Cα root-mean-square deviation between X and Y | Å | computed with Chimera, Cα atoms only |
| sRMSD | domain-based RMSD: align on one domain, compute RMSD on the other | Å | emphasizes large-scale rearrangement |
| improvement | morph quality metric | percent (%) | higher is better; see equations.md |
| F | threading (spring) force constant for gappedThreading/threading springs | MMB internal force-constant units (dimensionless setting) | e.g. 30 or 60; larger = stronger pull toward target |
| physics zone | region (radius around flexible residues) where the MD force field is active | Å (radius) | e.g. 10 Å around all flexible residues |
| flexibility zone | region granted internal-coordinate degrees of freedom (hinges, interfaces, active sites) | — | rigid outside = single kinematic body |
| convergence criterion | stop when energy difference between consecutive frames < threshold for N frames | kJ/mol | threshold 50 kJ/mol over 5 consecutive frames (ribosome runs) |

Conventions / notes:
- Multiscale kinematics: rigid regions are treated as a single rigid body in internal coordinates
  (MMB is built on the Simbody internal-coordinate engine; see Flores et al. 2011).
- Springs connect corresponding atoms between initial and (rigid) final structures via a gapped
  sequence alignment (SeqAn); these are the same springs used for homology modeling and rigid alignment.
- RMSD/sRMSD are structural comparison metrics only; they do not drive the dynamics (springs do).
