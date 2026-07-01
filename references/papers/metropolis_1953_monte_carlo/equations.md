# Equations - Metropolis et al. 1953

<!-- eq:1 -->
$$ E = \frac{1}{2} \sum_{\substack{i=1 \\ i \neq j}}^{N} \sum_{j=1}^{N} V(d_{ij}). $$
- **what:** Total configurational potential energy: half the sum over all ordered pairs of the pair potential evaluated at minimum-image distance.
- **symbols:** E - potential energy of the configuration; N - number of particles; V - pair potential; $d_{ij}$ - minimum-image distance between particles i and j; factor 1/2 corrects double counting of ordered pairs.

<!-- eq:2 -->
$$ \bar{F} = \left[ \int F \exp(-E/kT)\, d^{2N} p\, d^{2N} q \right] \Big/ \left[ \int \exp(-E/kT)\, d^{2N} p\, d^{2N} q \right]. $$
- **what:** Canonical-ensemble average of an observable F over the 4N-dimensional phase space; momentum integrals cancel, leaving a 2N-dimensional configuration-space average.
- **symbols:** $\bar{F}$ - ensemble average of F; F - observable; E - potential energy (Eq. 1); k - Boltzmann constant; T - temperature; $d^{2N}p\, d^{2N}q$ - phase-space volume element (2N momenta, 2N positions in 2D).

<!-- eq:3 -->
$$ \begin{array}{l}
X \to X + \alpha \xi_1 \\
Y \to Y + \alpha \xi_2,
\end{array} $$
- **what:** Trial-move proposal: displace a particle's x,y coordinates by uniform random offsets, giving a uniform proposal over a square of side $2\alpha$ centered on the current position (symmetric proposal -> $P_{rs}=P_{sr}$).
- **symbols:** X, Y - particle coordinates; $\alpha$ - maximum allowed displacement (step size); $\xi_1, \xi_2$ - independent random numbers uniform in $[-1, 1]$.

<!-- eq:4 -->
$$ \bar{F} = (1/M) \sum_{j=1}^{M} F_j, $$
- **what:** Monte Carlo estimator: sample mean of F over M configurations in the Markov chain. Every attempted move (accepted OR rejected) contributes one sample; a rejected move re-counts the current configuration.
- **symbols:** $\bar{F}$ - estimator of the ensemble average; M - number of moves/samples; $F_j$ - value of F after the jth move carried out per the full accept/reject prescription.

## Metropolis acceptance rule (from prose, not a numbered equation - the core of the paper)

<!-- eq:accept -->
$$ P_{\text{accept}} = \begin{cases} 1 & \Delta E \le 0 \\ \exp(-\Delta E/kT) & \Delta E > 0 \end{cases} $$
- **what:** Acceptance probability of a trial move. Draw $\xi_3 \sim U(0,1)$; accept if $\xi_3 < \exp(-\Delta E/kT)$ (always accept when $\Delta E<0$), else reject and keep the old position. Symmetric proposal makes this satisfy detailed balance w.r.t. the Boltzmann distribution.
- **symbols:** $\Delta E$ - energy change caused by the move ($E_{\text{new}} - E_{\text{old}}$); k - Boltzmann constant; T - temperature; $\xi_3$ - random number uniform in $[0,1]$.

<!-- eq:nu-canonical -->
$$ \nu_r \propto \exp(-E_r/kT). $$
- **what:** Target stationary distribution the chain converges to: number of ensemble members in state r is proportional to the Boltzmann factor of that state's energy.
- **symbols:** $\nu_r$ - number of ensemble systems in state r; $E_r$ - energy of state r; k - Boltzmann constant; T - temperature.

## Rigid-sphere equation of state

<!-- eq:7 -->
$$ \left\langle \sum_{i} \mathbf{X}_{i}^{(\text{tot})} \cdot \mathbf{r}_{i} \right\rangle_{\text{AV}} = 2PA + \left\langle \sum_{i} \mathbf{X}_{i}^{(\text{int})} \cdot \mathbf{r}_{i} \right\rangle_{\text{AV}} = 2E_{\text{kin}}. $$
- **what:** Clausius virial theorem in 2D relating total force virial, pressure-area work, internal-force virial, and kinetic energy.
- **symbols:** $\mathbf{X}_i^{(\text{tot})}$ - total force on particle i; $\mathbf{X}_i^{(\text{int})}$ - internal force on particle i; $\mathbf{r}_i$ - position of particle i; P - pressure; A - area (2D "volume"); $E_{\text{kin}}$ - total kinetic energy; AV - ensemble average.

<!-- eq:8 -->
$$ \bar{F}_{i} = m\bar{v}^{2}\pi d_{0}\bar{n}. $$
- **what:** Average collisional force on a central rigid disk from surrounding particles (isotropic Maxwellian velocity distribution).
- **symbols:** $\bar{F}_i$ - average force on central particle; m - particle mass; $\bar{v}^2$ - mean-square speed; $d_0$ - collision diameter (twice actual sphere radius); $\bar{n}$ - average number density of other particles at the surface.

<!-- eq:10 -->
$$ PA = E_{\text{kin}}(1 + \pi d_0^2 \bar{n}/2) \equiv NkT(1 + \pi d_0^2 \bar{n}/2). $$
- **what:** Equation of state for 2D rigid spheres: pressure-area product in terms of the contact density $\bar{n}$. Determining $\bar{n}(A)$ fixes the EOS.
- **symbols:** P - pressure; A - area; $E_{\text{kin}} = NkT$ - kinetic energy; N - particle count; k - Boltzmann constant; T - temperature; $d_0$ - collision diameter; $\bar{n}$ - contact/surface number density.

<!-- eq:11a -->
$$ d_0 = d(1 - 2^{\nu - 8}), \qquad d = (1/14). $$
- **what:** Forbidden (collision) distance parameterized by auxiliary variable $\nu$ at fixed lattice spacing d.
- **symbols:** $d_0$ - forbidden distance / collision diameter; d - nearest-neighbor spacing of the initial trigonal lattice = 1/14; $\nu$ - auxiliary parameter varied 0 to 7.

<!-- eq:11b -->
$$ (A/A_0) = 1/(3^{\frac{1}{2}} d_0^2 N/2) = 1/[0.98974329(1-2^{\nu-8})^2]. $$
- **what:** Ratio of actual area to close-packed area as a function of $\nu$. Used to map computed $\bar{n}$ onto the EOS curve.
- **symbols:** A - actual area (= 1, fixed); $A_0$ - close-packed area; $d_0$ - collision diameter (Eq. 11a); N - 224; $\nu$ - auxiliary parameter.

<!-- eq:deltaA2 -->
$$ \Delta A^2 = (K^2 - 1)\pi d_0^2/64. $$
- **what:** Equal-area bin width for the radial distribution histogram between $\pi d_0^2$ and $K^2\pi d_0^2$, split into 64 zones.
- **symbols:** $\Delta A^2$ - area of each histogram zone; K - outer-radius multiplier (>1, chosen per $\nu$); $d_0$ - collision diameter.

<!-- eq:12 -->
$$ (m-1)\Delta A^2 + \pi d_0^2 < \pi r^2 \le m\Delta A^2 + \pi d_0^2. $$
- **what:** Binning criterion: pair at separation r counts into histogram bin m if $\pi r^2$ falls in this equal-area zone. Counts $N_m$ form the radial distribution function.
- **symbols:** m - bin index (1..64); $\Delta A^2$ - zone area (Eq. deltaA2); $d_0$ - collision diameter; r - pair separation distance.

## Virial coefficient expansion

<!-- eq:13 -->
$$ (PA/NkT)-1 = C_1(A_0/A) + C_2(A_0/A)^2 + C_3(A_0/A)^3 + C_4(A_0/A)^4 + O(A_0/A)^5. $$
- **what:** Virial (density) expansion of the compressibility factor for the 2D rigid-sphere gas in powers of $(A_0/A)$ (proportional to density).
- **symbols:** P - pressure; A - area; $A_0$ - close-packed area; N, k, T - as above; $C_1..C_4$ - virial coefficients (Eq. 13-coeffs).

<!-- eq:13-coeffs -->
$$ C_1 = \pi/3^{\frac{1}{2}}, \qquad C_2 = 4\pi^2 A_{3,3}/9, $$
$$ C_3 = \pi^3 (6A_{4,5} - 3A_{4,4} - A_{4,6})/3^{\frac{3}{2}}, $$
$$ C_{4} = (8\pi^{3}/135) \cdot [12A_{5,5} - 60A_{5,6}' - 10A_{5,6}'' + 30A_{5,7}' + 60A_{5,7}'' + 10A_{5,7}''' - 30A_{5,8}' - 15A_{5,8}'' + 10A_{5,9} - A_{5,10}]. $$
- **what:** Virial coefficients in terms of cluster (irreducible area) integrals $A_{i,k}$.
- **symbols:** $C_1..C_4$ - virial coefficients; $A_{i,k}$ - cluster integral of i particles with k overlap bonds; primes distinguish topologically distinct bond arrangements.

<!-- eq:f-def -->
$$ f(r_{ij}) = 1 \text{ if } r_{ij} < d, \qquad f(r_{ij}) = 0 \text{ if } r_{ij} > d, $$
- **what:** Overlap (Mayer-like) bond indicator for hard disks: 1 when two particles overlap (separation less than diameter d), else 0.
- **symbols:** $f(r_{ij})$ - bond indicator; $r_{ij}$ - separation of particles i,j; d - hard-disk diameter (here the overlap threshold).

<!-- eq:A33 -->
$$ A_{3,3} = \frac{1}{\pi^2 d^4} \int \cdots \int dx_1 dx_2 dx_3\, dy_1 dy_2 dy_3\, (f_{12} f_{23} f_{31}). $$
- **what:** Example cluster integral: normalized configuration-space volume of 3 mutually overlapping disks (triangle of bonds).
- **symbols:** $A_{3,3}$ - three-particle three-bond cluster integral; d - diameter; $x_i,y_i$ - particle coordinates; $f_{ij}$ - bond indicators (Eq. f-def).

<!-- eq:14 -->
$$ (PA/NkT) - 1 = 1.813799(A_0/A) + 2.57269(A_0/A)^2 + 3.179(A_0/A)^3 + 3.38(A_0/A)^4 + O(A_0/A)^5. $$
- **what:** Final numerical four-term virial EOS for 2D rigid spheres (curve C).
- **symbols:** coefficients $C_1=1.813799$, $C_2=2.57269$, $C_3=3.179$, $C_4=3.38$; $A_0/A$ - reduced density variable.
