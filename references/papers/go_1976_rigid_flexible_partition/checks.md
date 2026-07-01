# Checks / fixtures - Gō & Scheraga 1976

This is a theory paper; the concrete numbers are asymptotics of the per-mode
sensitivity functions and frequency thresholds. They double as unit tests for any
implementation of the $(\det\mathbf{G})^{-1/2}$ metric weight (Fixman-type factor) and
the quantum/classical oscillator corrections.

## Per-mode sensitivity functions (eqs 22, 23)

Define $x = 2\pi\hbar\nu/kT = h\nu/kT$.

- $g_r(x) = -\tfrac12\coth(x/2)$.
- $g_f(x) = 1/x - \tfrac12\coth(x/2)$.

Fixtures:
- Given $x\to\infty$: expect $|g_r|\to \tfrac12$ and $|g_f|\to \tfrac12$ (both converge to 1/2 at high frequency).
- Given $x\to 0$: expect $|g_r|\to\infty$ (rigid model diverges / fails for softened modes) and $|g_f|\to 0$ (flexible model error vanishes for softened modes).
- Given any finite $x>0$: expect $|g_f| < |g_r|$ (flexible model strictly better). Both $g_r,g_f<0$.
  - Sanity fixture at $x=2.0$: $g_r = -\tfrac12\coth(1.0) \approx -0.6516$; $g_f = 1/2 - 0.6516 = -0.1516$. So $|g_f|\approx0.152 < |g_r|\approx0.652$.

## Classical-limit accuracy thresholds (statistical-weight error < 20%)

- Given the classical-limit replacement $[2\sinh(\pi\hbar\nu/kT)]^{-1}\to(2\pi\hbar\nu/kT)^{-1}$ (i.e. QM -> classical oscillator, eq 16 route): error < 20% of statistical weight for $2\pi\hbar\nu$ up to $2.1\,kT$, i.e. $x\le 2.1$.
  - At room temperature this corresponds to frequencies up to **440 cm$^{-1}$**.
- Given the ground-state replacement $[2\sinh(\pi\hbar\nu/kT)]^{-1}\to\exp(-\pi\hbar\nu/kT)$ (high-frequency / rigid-model route): error < 20% for $2\pi\hbar\nu$ down to $1.6\,kT$, i.e. $x\ge 1.6$.
  - At room temperature this corresponds to frequencies down to **330 cm$^{-1}$**.
- Most hard-variable vibrational frequencies lie $\ge 330$ cm$^{-1}$; the relevant infrared range quoted is **400-3000 cm$^{-1}$**.

## Model-equivalence limit

- Given the conformation dependence of the frequencies $\nu_i$ is neglected: expect flexible model (eq 16) and rigid model (eq 13) to become identical.
- Given all atom masses equal to $m_0$: expect $\mathbf{H}^0 = m_0 g$ (metric tensor) so $Z_r \propto \int(\det g)^{1/2}\exp[-\beta F_0]dQ$ (Kirkwood form, eq C-3). Note $\det\mathbf{G} \propto 1/\det g$, so the $(\det\mathbf{G})^{-1/2}$ weight equals $(\det g)^{+1/2}$ up to a mass constant.

## Jacobian identity (eq B-12, sanity for coordinate transform code)

- Given transform from bond-difference Cartesian vectors to (bond lengths $d_k$, bond angles $\tau_k$, dihedrals): expect Jacobian $D = \prod_{k=2}^{n} d_k^{-2} \cdot \prod_{k=3}^{n}\sin\tau_k$, **independent of all dihedral angles**. (Each per-atom step contributes $d_k^{-2}\sin\tau_k$.)

## Conformational entropy (Discussion)

- Given a stable conformation at an energy minimum with second-derivative (Hessian) matrix $\mathbf{F}$ of the conformational energy: expect conformational entropy $S = -\tfrac12 R\ln\det\mathbf{F}$ (harmonic approximation).

## Relative-stability correction (ref 22/23)

- Given relative stabilities computed with the $(\det\mathbf{G})^{-1/2}$ (rigid) weight vs without it: for the two lowest-energy conformations D' and E' of the cyclic molecule of ref 23, expect NO change in ordering because $\det\mathbf{G}$ is nearly equal for those conformations. Correct relative stabilities should use conformational free energies WITHOUT the $(RT/2)\ln\det\mathbf{G}$ term.
