# Notation - Stein, Rice & Brünger 1997 (Torsion-Angle MD)

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| E | total hybrid energy | kcal/mol | Eq. 1 |
| E_chem | chemical (force-field) energy | kcal/mol | E_geom + w_vdw E_vdw |
| E_geom | ideal-geometry energy (bonds, angles, planarity, chirality) | kcal/mol | |
| E_nmr | NMR-restraint energy | kcal/mol | Eq. 2 |
| E_NOE | NOE distance-restraint energy | kcal/mol | flat-bottom + soft asymptote, Eq. 6 |
| E_dihedral (E_cdih) | dihedral-angle restraint energy | kcal/mol | flat-bottom harmonic, Eq. 8 |
| E_vdw | van der Waals energy | kcal/mol | Eq. 5 (repel) during protocol; Eq. 4 (LJ) for final analysis |
| w_NOE | weight on NOE term | dimensionless | 150 (50 in minimization) |
| w_dihedral | weight on dihedral term | dimensionless | 100 (5 for nucleic acids; 300 in minimization) |
| w_vdw | weight on vdW term | dimensionless | 0.1 -> 1.0 across protocol |
| epsilon | Lennard-Jones well depth | kcal/mol | per atom-pair |
| sigma | Lennard-Jones distance parameter | Angstrom | per atom-pair |
| R | interatomic / spin-pair distance | Angstrom | |
| Delta | NOE distance violation | Angstrom | 0 inside [d_lower, d_upper], Eq. 7 |
| d_lower | lower distance bound (NOE) | Angstrom | |
| d_upper | upper distance bound (NOE) | Angstrom | switch to soft asymptote at d_upper + 0.5 |
| a, b | soft-asymptote continuity constants | kcal/mol (a), kcal/mol*Angstrom^softexp (b) | chosen for differentiability at R = d_upper + 0.5 |
| softexp | soft-asymptote exponent | dimensionless | |
| phi | dihedral angle | degrees | violations judged in degrees |
| phi_lower, phi_upper | dihedral bounds | degrees | |
| m_i | mass of atom i | atomic mass units | |
| r_{i,u} | Cartesian coordinate u of atom i | Angstrom | u in {x,y,z} |
| t | time | ps | |
| beta_i | Berendsen coupling force constant | (force/velocity) | Eq. 10 |
| T | instantaneous system temperature | K | |
| T_0 | bath (target) temperature | K | 50,000 K high-T (20,000 K nucleic acid) down to 300 K |
| v_i (vec) | velocity of atom i | Angstrom/ps | |
| r_i, r_j | centers of mass of bodies i, j | Angstrom, inertial frame | vector |
| r_ij | COM offset r_j - r_i | Angstrom | vector |
| h_ij | bond connecting bodies i and j | Angstrom | fixed length; vector |
| |h_ij| | bond length | Angstrom | constant (rigid) |
| h_hat_ij | unit vector along bond | dimensionless | h_ij / |h_ij| |
| s_ij | vector i-COM -> bond endpoint on body i | Angstrom | vector |
| s_ji | vector j-COM -> bond endpoint on body j | Angstrom | vector |
| q_ij | relative torsion angle about bond h_ij | rad | single rotational DOF per joint |
| q̇_ij | relative torsion-angle rate | rad/ps | |
| omega_i, omega_j | angular velocities of bodies i, j | rad/ps | vector |

## Conventions and notes

- Reduced/soft potential: during structure calculation the vdW term is a purely
  repulsive quartic (Eq. 5), NOT Lennard-Jones. LJ (Eq. 4) is used only for the
  final energetic analysis of accepted structures.
- vdW contact scaled by 0.8 (van der Waals radius scale) times 2^(1/6)*sigma.
- Solvent neglected: electrostatics excluded; hydrogen bonds modeled as pseudo-NOEs.
- Restraint energies are flat-bottomed: zero inside the experimental bounds,
  penalizing only violations. NOE has a soft asymptote beyond d_upper + 0.5 A to
  limit forces from large violations; dihedral restraint is plain harmonic outside
  its bounds.
- Torsion-angle MD imposes one rotational DOF per interbody bond (Eq. 11); rigid
  bond length |h_ij|. Closed (ring) bonding networks are handled approximately by
  allowing one bond in each ring to vibrate.
- Simulated annealing via Berendsen weak temperature coupling (Eq. 10).
