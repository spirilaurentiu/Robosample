# Equations — Tek et al. 2016 (MMB-GUI morphing)

This is a software/application paper; it contains a single implementable formula: the morphing
*improvement* benchmark metric (adopted from Weiss & Levitt 2009). Everything else is protocol/prose.

<!-- eq:improvement -->
$$ \text{improvement} = \frac{\min[\text{rmsd}(AB),\, \text{rmsd}(CB)] - \min_i[\text{rmsd}(iB)]}{\min[\text{rmsd}(AB),\, \text{rmsd}(CB)]} \times 100\% $$
- **what:** Percentage by which the best interpolated morph structure improves on the trivial endpoint
  approximations of the known intermediate. It measures how much closer the closest-approach morph frame
  comes to the intermediate B than either endpoint (A or C) does.
- **symbols:** A - initial crystal structure; B - experimentally observed intermediate crystal structure;
  C - final crystal structure; i - an interpolated (morphed) structure along the trajectory;
  rmsd(XY) - Cα root-mean-square deviation between structures X and Y (Å), computed after alignment
  (Chimera); min over i - taken over all morph frames (closest approach to B); the numerator's first
  term is the better (smaller) endpoint-vs-B RMSD.
<!-- CHECK: the paper's rendered numerator reads "min[rmsd(AB),rmsd(CB)] - min[rmsd(iB)]"; the second
     min is over the morph-frame index i (closest approach). Reconstructed as min_i rmsd(iB) from the
     surrounding text ("interpolated structure ... at closest approach"). -->
