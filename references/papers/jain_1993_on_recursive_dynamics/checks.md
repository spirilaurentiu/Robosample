# Checks / fixtures — Jain, Vaidehi, Rodriguez 1993

Concrete numbers stated by the paper, as regression fixtures.

## Complexity / flop-count fixtures (Section 5.4)

- **O(N) worst case:** given a serial chain, no point-mass clusters, only
  single-dof hinges, expect the recursive algorithm cost ≈ **$500\,\mathcal{N}$**
  floating-point operations. Point masses, multi-dof hinges, or branches reduce
  this.
- **O(N^3) reference:** given a serial chain, single-dof hinges, no point masses,
  expect the conventional (form-$\mathcal{M}$-and-solve) cost ≈
  **$\mathcal{N}^3/3 + 19\mathcal{N}^2 + 350\mathcal{N}$** flops.
- **Speedup fixture:** given a polypeptide with each residue a rigid cluster and
  2 bending dof between neighbors, so $\mathcal{N}\approx 2\times(\text{residues})$;
  for **400 residues** ($\mathcal{N}\approx 800$) expect the $O(\mathcal{N}^3)$
  cost ≈ **450×** the $O(\mathcal{N})$ cost.
  - Sanity cross-check at $\mathcal{N}=800$: $O(\mathcal{N}^3)$ ≈
    $800^3/3 + 19\cdot800^2 + 350\cdot800 \approx 1.71\times10^8 + 1.216\times10^7
    + 2.8\times10^5 \approx 1.836\times10^8$; $O(\mathcal{N})\approx 500\cdot800 =
    4.0\times10^5$; ratio ≈ **459×** — consistent with the stated ~450×.
    <!-- CHECK: ratio recomputed here from the two flop formulas; paper only states "~450". -->
- **General claim:** internal-variable models reduce MD cost by ~1 order of
  magnitude vs Cartesian; the $O(\mathcal{N})$ algorithm adds a further large
  speedup that grows with molecule size.

## Structural invariants an implementation must satisfy

- **Mass-matrix factorization consistency:** the Newton-Euler factorization
  $\mathcal{M}=H\phi M\phi^*H^*$ (eq:4.15a) and the innovations factorization
  $\mathcal{M}=[I+H\phi K]D[I+H\phi K]^*$ (eq:5.6) must produce the **same**
  $\mathcal{M}$. Test: for a random tree of clusters, assemble $\mathcal{M}$ both
  ways and require equality to machine tolerance.
- **Inverse identity:** $[I+H\phi K][I-H\psi K] = I$ (eq:5.7). Test the product
  equals identity.
- **Mass-matrix inverse:** $\mathcal{M}\,\mathcal{M}^{-1}=I$ with
  $\mathcal{M}^{-1}=[I-H\psi K]^*D^{-1}[I-H\psi K]$ (eq:5.8).
- **Solver equivalence:** the O(N) recursion (eq:5.11a+eq:5.11b) must return the
  same $\ddot\theta$ as directly solving $\mathcal{M}\ddot\theta = T-\mathcal{C}$
  (eq:4.14) via a dense linear solve, for the same inputs $T,\hat f_c,a,b$. This
  is the primary regression test for a port.
- **Riccati positivity:** each $P(k)$ (eq:5.1) is symmetric PD; each $D(k)$ is
  symmetric PD and invertible.
- **Spatial transform form:** $\phi(\mathcal{O}_x,\mathcal{O}_y)=[[I_3,\tilde
  l],[0_3,I_3]]$ (eq:A.2) has unit determinant and $\phi^{-1}$ = same with
  $-\tilde l$.
- **Zero-extent atom cluster:** a single atom is a cluster with zero extent and
  zero rotational inertia ($\mathcal{I}=0$, $p=0$ in eq:A.4).

No benchmark energy/trajectory tables are given (results deferred to a
"forthcoming publication").
