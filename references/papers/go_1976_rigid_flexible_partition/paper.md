# On the Use of Classical Statistical Mechanics in the Treatment of Polymer Chain Conformation

Nobuhiro Gō and Harold A. Scheraga. *Macromolecules* **9**(4), 535-542 (1976).

## Abstract

Two different treatments of the degrees of freedom of bond stretching and bond angle
bending in chain polymers by classical statistical mechanics lead to different and
nonequivalent expressions of the partition functions. If we fix the bond lengths and
bond angles at the outset and treat them as constraints (the **classical rigid model**),
the partition function is given by an integral of $(\det \mathbf{G})^{-1/2}\exp[-\beta F(Q)]$
over the space of the dihedral angles $Q$, where the elements of the matrix $\mathbf{G}$
are the coefficients in the quadratic expression (in generalized momenta conjugate to
$Q$) for the kinetic energy of the chain, and $F(Q)$ is the conformational energy. If we
conceptually allow bond lengths and bond angles to vary under an infinitely strong
potential (the **classical flexible model**) and integrate the Boltzmann factor over the
momenta conjugate to the Cartesian coordinates, we obtain the partition function as an
integral of $\exp[-\beta F(Q)]$ over the dihedral angles $Q$ (no metric weight). The origin
of the difference lies in the different treatments of the vibrational motions of bond
lengths and bond angles. A quantum-mechanically correct partition function for these
vibrations is derived, and the approximations that reduce it to each classical form are
examined. The **classical rigid model** follows by (a) the ground-state approximation for
all bond-stretch/bend vibrations (neglecting excited vibrational states) and (b)
neglecting the conformational dependence of the zero-point energy. The **classical
flexible model** follows by treating all these vibrations classically. A quantitative
analysis shows that, of the two nonequivalent classical treatments, the **classical
flexible model (eq 1) is the better approximation** than the classical rigid model (eq 2).

## Background

Equilibrium conformational properties have been calculated from a classical partition
function of the form (eq 1)

$$ Z = (\text{constant}) \int \exp[-\beta F(Q)]\, dQ $$

where $\beta=1/kT$, $Q$ is a set of dihedral angles, and $F(Q)$ is the conformational (free)
energy / potential of mean force (it includes the free energy of solvation as a function
of $Q$). It was suggested that one might instead need (eq 2)

$$ Z = (\text{constant}) \int \left[\frac{1}{\det \mathbf{G}}\right]^{1/2}\exp[-\beta F(Q)]\, dQ $$

derived by fixing bond lengths and bond angles as constraints. The present paper explores
the origin of the difference and shows that eq 1 is more accurate than eq 2, and that eq 1
avoids the difficulty of the intractable conformation-dependent $\det\mathbf{G}$ factor.

The paper's earlier work (ref 6) erroneously supported eq 2; this paper corrects that,
in agreement with Flory (ref 7). The approximations reducing the QM expression to eq 2 are
introduced in three stages (approximations A, B, C of ref 6): (A) ground-state
approximation for hard-variable vibrations, retaining conformation-dependent zero-point
energy; (B) neglect of the conformational dependence of that zero-point energy; (C) neglect
of the dependence of the hard-variable averages on the soft variables (this is eq 2).

## I. Nonequivalence of the Two Classical Models

**A. Classical flexible model.** Independent variables are bond lengths, bond angles,
dihedral angles, and 6 external (overall translation + rotation) variables. Hard variables
$Q'=(q_1',\dots,q_l')$ = bond lengths + bond angles. Soft variables $Q=(q_1,\dots,q_m)$ =
dihedral + external. Vibrations of the hard variables are governed by a harmonic force
field, giving the conformational energy (eq 3)

$$ F(Q,Q') = F_0(Q) + \frac{1}{2}\sum_{i,j=1}^{l} f_{ij}''(q_i'-q_{i0}')(q_j'-q_{j0}') $$

The strain-free values $q_{i0}'$ and the force constants $f_{ij}''$ are approximated as
Q-independent (this is the same approximation used to obtain approximation C in ref 6, and
is justified because $f_{ij}''$ are set by covalent-bond force fields insensitive to
conformation). The kinetic energy is eq 4 (Cartesian) or eq 5 (internal), with the full
$3n\times3n$ mass metric $\mathcal{H}$ block-partitioned into soft/hard blocks
$\mathbf{H}^0,\mathbf{H}',\mathbf{H}''$. The flexible Hamiltonian is $F(Q,Q')+K_f$.

**B. Classical rigid model.** Hard variables kept fixed: set $\dot Q'=0$ in eq 5, giving
$K_r=\tfrac12\dot Q^+\mathbf{H}^0\dot Q$ (eq 6). Conformational energy is simply $F_0(Q)$.

**C. Partition function, flexible model.** Expressing $K_f$ in Cartesian momenta
$p_{k\alpha}=m_k\dot x_{k\alpha}$ and integrating over momentum space gives eq 7; the
prefactor is Q-independent. Changing variables Cartesian -> internal introduces the Jacobian
$D$ (Appendix B), which is independent of dihedral angles and factorizes into external- and
hard-variable factors. Assigning minimum-energy hard values $Q_0'$ and integrating out the
external and hard variables yields eq 8 (= eq 1): a flat integral over the internal soft
variables. The bond-stretch and bend degrees of freedom do not appear explicitly; they are
varied only *conceptually*.

**D. Partition function, rigid model.** Rewrite $K_r=\tfrac12 P_r^+\mathbf{G}P_r$ (eq 9)
with $P_r=\mathbf{H}^0\dot Q$ and $\mathbf{G}=(\mathbf{H}^0)^{-1}$ given by the Schur
complement (eq 10) of the inverse full metric $\mathcal{G}=\mathcal{H}^{-1}$ (eq 11). The
phase-space partition function (eq 12), after Gaussian momentum integration and integrating
out the external variables, yields eq 13 (= eq 2), which carries the extra factor
$(\det\mathbf{G})^{-1/2}$. This factor is a function of the internal soft variables $Q$ and
CANNOT be pulled out of the integral - this is the intractability Flory pointed out.

**E. Comparison of the two models.** Using the Jacobian identity (eq 14)

$$ D = \left[\det\mathbf{G}\,\det\mathbf{G}''\left(\prod_{k=1}^n m_k\right)^3\right]^{-1/2} $$

and the Wilson GF relation $\det(\mathbf{F}''\mathbf{G}'')=\prod_i(2\pi\nu_i)^2$ (eq 15),
$Z_f$ (eq 8) can be rewritten as eq 16:

$$ Z_f = \left(\frac{kT}{2\pi\hbar^2}\right)^{m/2}(8\pi^2 V)\int\left[\prod_{i=1}^l\frac{kT}{2\pi\hbar\nu_i}\right]\left[\frac{1}{\det\mathbf{G}}\right]^{1/2}\exp[-\beta F_0(Q)]\,dQ $$

Since eq 8 and eq 16 are the same quantity, the product of the first two integrand factors
in eq 16 is Q-independent, though each factor separately is not. **The flexible model equals
the rigid model plus a set of $l$ classical harmonic oscillators with conformation-dependent
frequencies $\nu_i$.** If the frequencies' conformation dependence were neglected, the two
models would be equivalent. Comparing eq 16 with eq 13, the rigid model is obtained by
replacing $\prod_i(kT/2\pi\hbar\nu_i)$ by unity. In the limit of infinitely strong hard
potentials, $\prod_i(kT/2\pi\hbar\nu_i)\to0$ (flexible) while the rigid model keeps it at
unity - a superficial paradox arising from the deficiency of classical statistical mechanics
for high-frequency motions, not from flexibility vs rigidity per se.

## II. Approximations Involved in the Two Classical Models

Because hard-variable frequencies are high vs $kT/h$, the correct treatment is quantum.
Replacing the classical oscillator factor $(2\pi\hbar\nu_i/kT)^{-1}$ in eq 16 by the QM one
$[2\sinh(\pi\hbar\nu_i/kT)]^{-1}$ gives the quantum-correct partition function eq 17. (Soft
variables may be kept classical - a very good approximation.)

$Z_f$ is recovered from $Z_{\rm QM}$ by the classical limit
$[2\sinh(\pi\hbar\nu_i/kT)]^{-1}\to(2\pi\hbar\nu_i/kT)^{-1}$, valid for small $\nu_i$. The
error is < 20% of the statistical weight for $2\pi\hbar\nu_i\le2.1kT$ (up to 440 cm$^{-1}$ at
room T). $Z_r$ is recovered by the high-frequency / ground-state replacement
$[2\sinh]^{-1}\to\exp(-\pi\hbar\nu_i/kT)$, valid with < 20% error for $2\pi\hbar\nu_i\ge1.6kT$
(down to 330 cm$^{-1}$), PLUS neglecting the conformation dependence of the zero-point
energy. Both introduce error, worst for softened low-frequency coupled modes (Hagler &
Lifson) where the ground-state approximation is invalid.

Define the ratios of the QM vibrational partition function to the corresponding factors in
$Z_r$ (unity) and $Z_f$ (classical oscillator), eqs 18-19:

$$ \Gamma_r = \prod_{i=1}^l\left(\frac{1}{2\sinh(x_i/2)}\right),\qquad \Gamma_f = \prod_{i=1}^l\left(\frac{x_i}{2\sinh(x_i/2)}\right),\qquad x_i=\frac{2\pi\hbar\nu_i}{kT} $$

Neither is conformation-independent; the smaller the conformation dependence, the better the
classical model ($\Gamma_f$ equals Flory's $\Gamma$). Relative to a reference conformation
(eqs 20-23),

$$ \ln(\Gamma_r/\Gamma_{r_0})=\sum_i g_r(x_{i0})\Delta x_i,\qquad \ln(\Gamma_f/\Gamma_{f_0})=\sum_i g_f(x_{i0})\Delta x_i $$
$$ g_r(x_{i0})=-\tfrac12\coth\tfrac{x_{i0}}{2},\qquad g_f(x_{i0})=\frac{1}{x_{i0}}-\tfrac12\coth\tfrac{x_{i0}}{2} $$

Both $g_r,g_f<0$; both $|g|\to\tfrac12$ as $x\to\infty$; for all finite $x$, $|g_f|<|g_r|$.
As $x\to0$, $|g_r|\to\infty$ (rigid model fails for softened modes) while $|g_f|\to0$
(flexible model error vanishes). Hence the distribution of $\ln(\Gamma_r/\Gamma_{r_0})$ is
broader than that of $\ln(\Gamma_f/\Gamma_{f_0})$: **the classical flexible model is the
better approximation.** Since most frequencies are high ($x\ge2.0$), the rigid model is not
totally unrealistic, but the flexible model is superior.

## III. Discussion

The earlier erroneous conclusion (ref 6) favoring eq 2 came from inadequate analysis of the
approximations. The rigid model's two approximations (ground-state approximation - unwarranted
for softened modes - and neglect of zero-point-energy conformation dependence) together
introduce more error than the classical treatment of high-frequency vibrations in eq 1.
Therefore **eq 1 is a better approximation than eq 2.** The flexible treatment is especially
good for the softened coupled motions that are the weakest point of the rigid model. In an
absolute sense eq 1 remains acceptable because high-frequency vibrations are localized and
insensitive to conformation ($\Delta x_i\to0$ as $x_{i0}\to\infty$), while $g_f(x_{i0})\to0$
as $x_{i0}\to0$, so appreciable error would arise only in an intermediate frequency range,
for which there is presently no evidence.

Another key approximation is neglecting the dependence of $Q_0'$ and $f_{ij}''$ on $Q$ in
eq 3; strictly $Q_0'$ should be $Q_0'-\mathbf{F}''^{-1}\mathbf{f}'$ (Appendix D), Q-dependent
through $\mathbf{f}'$. Retaining it would require varying bond lengths/angles with
conformation - intractable for large molecules. Note: assuming $f_{ij}''$ is Q-independent is
NOT the same as assuming the vibrational frequencies are Q-independent; the frequency
dependence enters through $\mathbf{G}''(Q)$ (eq A-8).

**Solvent.** Once the classical treatment is accepted, solvent enters through
$F(Q)=U(Q)+V(Q)$ (intramolecular energy + free energy of solvation), i.e. the potential of
mean force; the polymer is then treated as if in vacuum under this potential.

**Conformational entropy.** For a stable conformation with small fluctuations, minimize the
conformational energy, then compute the entropy from the second-derivative matrix
$\mathbf{F}$ at the minimum as $-\tfrac12 R\ln\det\mathbf{F}$.

## Appendices (derivation-heavy; see equations.md "Derivations")

- **Appendix A** — derivation of eq 14: explicit mass-weighted forms of the metric blocks
  $\mathbf{H}^0,\mathbf{H}',\mathbf{H}''$ and $\mathbf{G}^0,\mathbf{G}',\mathbf{G}''$ from
  Cartesian-internal Jacobians; leads to $D=[\det\mathbf{G}\det\mathbf{G}''(\prod_k m_k)^3]^{-1/2}$.
- **Appendix B** — proof that the Cartesian->internal Jacobian is independent of dihedral
  angles: stepwise change of variables yields $D=\prod_k d_k^{-2}\cdot\prod_k\sin\tau_k$
  (bond lengths $d_k$, bond angles $\tau_k$), each atom step contributing $d_k^{-2}\sin\tau_k$;
  used to derive eq 8.
- **Appendix C** (Historical, non-equilibrium) — with all masses equal $m_0$,
  $\mathbf{H}^0=m_0 g$ with $g$ the geometric metric tensor of the constrained space, giving
  the Kirkwood form $Z_r\propto\int(\det g)^{1/2}\exp[-\beta F_0]dQ$; shows Kirkwood (1949)
  used the classical rigid model for diffusion in constrained space.
- **Appendix D** — re-examines the separation of hard/soft motions: rewriting the flexible
  Hamiltonian decouples it into hard-oscillator + soft-variable parts (with the soft KE
  matrix changing from $\mathbf{G}^0$ to $\mathbf{G}$), justifying replacing classical by QM
  oscillator factors (eq 16 -> eq 17); also gives the corrected interpretation
  $Q_0'\to Q_0'-\mathbf{F}''^{-1}\mathbf{f}'$.

> Note: the Appendix A/B/D display equations are severely OCR-corrupted in the source
> conversion (Greek letters, subscripts, and integral/sum signs lost). The clean statements
> above are reconstructed from the main-text prose and standard Wilson GF-matrix theory.
