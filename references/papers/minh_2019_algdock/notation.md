# Notation - Minh 2019 (AlGDock)

Convention: reduced units throughout. Reduced potential $u = U/(k_B T)$; reduced free energies $f$ in units of $k_B T$. $\beta = (k_B T)^{-1}$. Standard = reduced / $\beta$.

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $B(r_R)$ | binding potential of mean force | energy | eq:1; function of receptor conf |
| $\beta$ | inverse temperature | $(k_B T)^{-1}$ | $\beta=(k_B T)^{-1}$ |
| $\beta_T$ | target inverse temp | $(k_B\cdot300\,\text{K})^{-1}$ | $\beta_T^{-1}=k_B(300\text{ K})$ |
| $\beta_H$ | high inverse temp | $(k_B\cdot600\,\text{K})^{-1}$ | $\beta_H^{-1}=k_B(600\text{ K})$ |
| $T_T$ | target temperature | 300 K | milestone D/E |
| $T_H$ | high temperature | 600 K | milestone C |
| $T(\alpha)$ | temperature ramp | K | $=(T_T-T_H)\alpha+T_H$; $T(0)=T_H$, $T(1)=T_T$ |
| $r_{RL}$ | complex internal coordinates | vector | excludes overall trans/rot |
| $r_R$ | receptor internal coordinates | vector | rigid receptor |
| $r_L$ | ligand internal coordinates | vector | flexible ligand |
| $\xi$ | relative translation + rotation of species | 6 DOF | |
| $I(\xi)$ | bound indicator function | {0,1} | 1 bound, 0 unbound |
| $J(\xi)$ | Jacobian of coordinate transform | scalar | Cartesian -> $(r_L,\xi)$ |
| $U(\cdot)$ | potential energy in solvent | energy | MM + implicit solvent |
| $U_S$ | sampling force-field energy | energy | MMTK; ligand-only solvation |
| $U_T$ | target force-field energy | energy | OpenMM; full-complex solvation |
| $u_\lambda(x)$ | reduced potential | $k_B T$ | $=U_\lambda/(k_B T_\lambda)$ |
| $u_I(d)$ | flat-bottom restraint (reduced) | $k_B T$ | eq:2 |
| $d$ | ligand COM distance to site center | nm | |
| $d_0$ | binding-site radius | 6.0 Å = 0.6 nm | |
| $k$ | restraint spring constant | 10000 kJ/(mol nm$^2$) | eq:2 |
| $\Psi_g$ | grid interaction energy | kJ/mol | $=\Psi_{PBSA}+\Psi_{vdW}$ (eq:3) |
| $\Psi_{PBSA}$ | electrostatic grid energy | kJ/mol | trilinear interp of PB potential |
| $\Psi_{vdW}$ | vdW grid energy | kJ/mol | transformed grid interp |
| $\Psi_{sg}$ | soft-grid interaction energy | kJ/mol | tanh-capped; $v_{max}\tanh(v_o/v_{max})$ |
| $v_{max}$ | soft LJ repulsive cap | 10.0 kJ mol$^{-1/2}$ | electrostatic cap = 10x min ratio |
| $\alpha$ | progress variable (states CD) | [0,1] | 0 = milestone C, 1 = milestone D |
| $\alpha_{sg}(\alpha)$ | soft-grid scaling | dimensionless | $-(2\alpha-1)^2+1$ |
| $\alpha_g(\alpha)$ | unperturbed-grid scaling | dimensionless | sigmoid-gated parabola (eq:4) |
| $p_{acc}$ | replica-exchange acceptance prob | [0,1] | eq:5 |
| $\langle p_{acc}\rangle$ | mean exchange rate | [0,1] | target 0.4-0.99; insert state if < 0.4 |
| $\bar{p}_{acc}$ | estimator for $\langle p_{acc}\rangle$ | [0,1] | sample mean of eq:5 |
| $K$ | number of thermodynamic states in a direction | int | exchanges up to $\min(5,K)$ apart |
| $N_{states}$ | number of thermodynamic states | int | system-specific |
| $\mathcal{L}$ | thermodynamic length | dimensionless | eq:6 |
| $g(\gamma)_{ij}$ | Fisher-information metric | matrix | eq:7 |
| $l_\lambda(x)$ | normalized log probability | $k_B T$ | $=-u_\lambda-\ln Z_\lambda$ |
| $Z_\lambda$ | partition function of state $\lambda$ | | $=\int e^{-u_\lambda}dx$ |
| $\lambda, \lambda^i$ | thermodynamic parameter vector / component | | temperature, grid scalings |
| $s$ | thermodynamic speed | adjustable | eq:8 |
| $s_{bc}$ | thermodynamic speed for states BC | 20.0 | eq:9 |
| $s_{cd}$ | thermodynamic speed for states CD | 0.2 | eq:10 |
| $\sigma_\lambda[\cdot]$ | std dev of quantity in state $\lambda$ | | |
| $f_{XY}$ | reduced free energy diff between milestones X,Y | $k_B T$ | |
| $f'_{CD}$ | reduced free energy (eq:11) | $k_B T$ | avoids $U(r_R)$ |
| $f_{EE_p}$ | reduced free energy of pose $p$ | $k_B T$ | eq:14 |
| $f_{AE,p}$ | pose-specific BPMF | $k_B T$ | eq:15 |
| $w_c$ | configuration reweighting factor D->E | unnormalized | eq:12/13 |

## Force fields and models
- Proteins/ions: AMBER ff14SB. Other molecules: GAFF2 with AM1BCC charges.
- Implicit solvent (sampling): generalized Born/surface area model II (OBC) from Onufriev et al., adapted from OpenMM; ligand-only solvation.
- PB grid: APBS 1.4, linear Poisson-Boltzmann, protein dielectric 2.0, solvent dielectric 80.0, solvent radius 1.4 Å, 300 K, fine spacing 0.5 Å.
- vdW grid transform: inverse-transformation power 4 for repulsive; no transform for attractive.

## Milestones
- A: target force field, full complex solvation. E: equivalent to A (both endpoints of pathways).
- B: sampling force field. C: high temperature (600 K), grids off. D: target temperature (300 K), grids fully on.
- States XY = intermediate states between milestones X and Y. Pathways Desolvated / Full differ in implicit solvent treatment across BD.
