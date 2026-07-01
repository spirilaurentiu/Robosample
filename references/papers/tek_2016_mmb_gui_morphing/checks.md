# Checks / fixtures — Tek et al. 2016 (MMB-GUI morphing)

## Improvement-score benchmark (Table 1, Weiss & Levitt dataset)

For each protein: A = initial, B = intermediate, C = final. `improvement` computed per `equations.md`
from the RMSD columns. RMSD in Å.

| Protein | rmsd(A,B) | rmsd(B,C) | min rmsd(morph, B) | improvement | MMB rank | best benchmark improvement (method) |
|---|---|---|---|---|---|---|
| 5′-Nucleotidase | 5.42 | 4.72 | 1.42 | 70% | 1 | 68% (FATCAT) |
| Ca²⁺-ATPase | 13.75 | 10.10 | 8.67 | 14% | 1 (a) | 12% (Climber) |
| Myosin | 16.56 | 12.01 | 3.84 | 68% | 1 | 58% (NOMAD-Ref) |
| RBP (ribose-binding protein) | 2.22 | 4.20 | 0.77 | 65% | 2 | 68% (FATCAT) |
| RNase III | 7.26 | 13.15 | 3.96 | 45% | 2 (a) | 50% (NOMAD-Ref) |

Fixtures for the improvement formula:
- given rmsd(AB)=16.56, rmsd(CB)=12.01, min_i rmsd(iB)=3.84 -> expect improvement ≈ (min(16.56,12.01) − 3.84)/min(16.56,12.01) × 100% = (12.01 − 3.84)/12.01 × 100% ≈ 68%.
- given rmsd(AB)=2.22, rmsd(CB)=4.20, min_i rmsd(iB)=0.77 -> expect improvement = (min(2.22,4.20) − 0.77)/min(2.22,4.20) × 100% = (2.22 − 0.77)/2.22 × 100% ≈ 65%.
- given rmsd(AB)=5.42, rmsd(CB)=4.72, min_i rmsd(iB)=1.42 -> expect (4.72 − 1.42)/4.72 × 100% ≈ 70%.
- given rmsd(AB)=13.75, rmsd(CB)=10.10, min_i rmsd(iB)=8.67 -> expect (10.10 − 8.67)/10.10 × 100% ≈ 14%.
- given rmsd(AB)=7.26, rmsd(CB)=13.15, min_i rmsd(iB)=3.96 -> expect (7.26 − 3.96)/7.26 × 100% ≈ 45%.

(These confirm the improvement metric uses the *smaller* of the two endpoint-vs-B RMSDs as the denominator.)

Myosin PDB IDs (from Table 1 caption): A = 1QVI, B = 1KK7, C = 1KK8.
RBP example (Figure 2): initial 1BA2, final 2DRI, intermediate 1URP; min RMSD vs intermediate 0.77 Å.

## Ca²⁺-ATPase alternate flexibility scheme (Results text)
- Uniform domain-hinge-bending protocol (HingeMaster hinges): minimum RMSD vs final = 5.13 Å.
- Flexibilizing residues 42–47, 57–59, 80–84, 112–114, 122–126: RMSD vs final = 4.15 Å,
  but RMSD vs intermediate increased to 9.09 Å.
- Weiss & Levitt's Climber with extreme parameters: 16% (Ca²⁺-ATPase), 48% (RNase III) improvement.
- NOMAD-Ref: 50% on RNase III (but erratic trajectory).

## Ribosome translocation morph (large-system parameters)
- No-protein (except EF-G) morphs: 144 791 to 155 378 moving atoms; 20–43 min to converge on a laptop.
- Full-ribosome (states 4→5, most proteins) system: 249 313 atoms; ~45 min to converge.
- Convergence criterion: energy difference between consecutive frames < 50 kJ/mol for 5 consecutive frames.
- Physics zone radius (Weiss & Levitt benchmark morphs): 10 Å around all flexible residues.
- Threading force constant F: P-site tRNA needed F increased 30 → 60 to pass the gate (S7 contact);
  with F = 30 the tRNA became stuck in the gate.
- Base structure: 3.3 Å crystal structure from T. thermophilus, PDB 2WDG and 2WDI.
- Predictive test: morph state 1 (2WDG/2WDI) → state 3 (3J5W/3J5X), omitting state 2 (3J5T/3J5U);
  intermediate frame #20 most closely resembles the omitted A/P tRNA state 2.

## Structural/geometric observations (approximate, morph-derived)
- Intersubunit rotation: small subunit rotates clockwise up to 12° relative to large subunit.
- A→P translocation: A-tRNA CCA 3′-end travels up to ~10 Å; elbow up to ~40 Å.
- P→E translocation: P-tRNA CCA end ~40 Å; elbow ~60 Å.
- Anticodon-loop tip (residue 34) sweep: ~14 Å (A→P), ~17 Å (P→E).
- Gate width: ~13 Å (non-/partially rotated) widening to > 20 Å with head swivel.
- ASL 'untwisting': ASL rotated by nearly 30° relative to pE/E state.
- Claimed multiscale speedup: up to ~2000-fold vs conventional simulation.

## Hardware baseline
- Benchmark morphs: single laptop CPU core, MacBook Air, Intel Core i5 @ 2.8 GHz, OS X 10.8.5;
  convergence in a few minutes.
