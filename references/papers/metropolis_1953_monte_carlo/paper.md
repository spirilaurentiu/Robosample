# Equation of State Calculations by Fast Computing Machines

**Authors:** Nicholas Metropolis, Arianna W. Rosenbluth, Marshall N. Rosenbluth, Augusta H. Teller (Los Alamos Scientific Laboratory), Edward Teller (University of Chicago).
**Venue:** J. Chem. Phys. 21, 1087 (1953). DOI: 10.1063/1.1699114.

## Abstract

A general method, suitable for fast computing machines, for investigating such
properties as equations of state for substances consisting of interacting
individual molecules is described. The method consists of a modified Monte Carlo
integration over configuration space. Results for the two-dimensional
rigid-sphere system have been obtained on the Los Alamos MANIAC and are presented
here. These results are compared to the free volume equation of state and to a
four-term virial coefficient expansion.

## I. Introduction

The purpose of this paper is to describe a general method, suitable for fast
electronic computing machines, of calculating the properties of any substance
which may be considered as composed of interacting individual molecules.
Classical statistics is assumed, only two-body forces are considered, and the
potential field of a molecule is assumed spherically symmetric. These are the
usual assumptions made in theories of liquids. Subject to the above assumptions,
the method is not restricted to any range of temperature or density. This paper
also presents results of a preliminary two-dimensional calculation for the
rigid-sphere system.

## II. The General Method for an Arbitrary Potential Between the Particles

To reduce the problem to a feasible size for numerical work we consider only a
finite number of particles. This number `N` may be as high as several hundred.
Our system consists of a square containing `N` particles. To minimize surface
effects we suppose the complete substance to be periodic, consisting of many such
squares, each square containing `N` particles in the same configuration. Thus we
define $d_{AB}$, the minimum distance between particles A and B, as the shortest
distance between A and any of the periodic images of B, of which there is one in
each square. If the potential falls off rapidly with distance, at most one of the
distances AB can make a substantial contribution; hence we need consider only the
minimum distance $d_{AB}$. (Two-dimensional nomenclature is used; the extension to
three dimensions is obvious.)

If we know the positions of the `N` particles in the square, we can calculate the
potential energy of the system:

<!-- eq:1 -->
$$ E = \frac{1}{2} \sum_{\substack{i=1 \\ i \neq j}}^{N} \sum_{j=1}^{N} V(d_{ij}). $$

Here `V` is the potential between molecules, and $d_{ij}$ is the minimum distance
between particles `i` and `j` as defined above.

To calculate properties we use the canonical ensemble. The equilibrium value of
any quantity of interest `F` is

<!-- eq:2 -->
$$ \bar{F} = \left[ \int F \exp(-E/kT)\, d^{2N} p\, d^{2N} q \right] \Big/ \left[ \int \exp(-E/kT)\, d^{2N} p\, d^{2N} q \right]. $$

$(d^{2N} p\, d^{2N} q)$ is a volume element in the 4N-dimensional phase space.
Since forces between particles are velocity-independent, the momentum integrals
separate off, and we need integrate only over the 2N-dimensional configuration
space. It is impractical to carry out a several-hundred-dimensional integral by
usual numerical methods, so we resort to the Monte Carlo method: integrate over a
random sampling of points instead of over a regular array of points.

The most naive method would be to place each of the `N` particles at a random
position (defining a random point in the 2N-dimensional configuration space),
calculate `E` via Eq. (1), and give the configuration weight $\exp(-E/kT)$. This
is impractical for close-packed configurations, since with high probability one
chooses a configuration where $\exp(-E/kT)$ is very small. So we use a modified
Monte Carlo scheme: instead of choosing configurations randomly and weighting
them with $\exp(-E/kT)$, we **choose configurations with probability
$\exp(-E/kT)$ and weight them evenly.**

We proceed as follows. Place the `N` particles in any configuration, e.g. a
regular lattice. Then move each particle in succession according to:

<!-- eq:3 -->
$$ \begin{array}{l}
X \to X + \alpha \xi_1 \\
Y \to Y + \alpha \xi_2,
\end{array} $$

where $\alpha$ is the maximum allowed displacement (arbitrary for this argument),
and $\xi_1, \xi_2$ are random numbers between $-1$ and $1$. After a move, the
particle is equally likely to be anywhere within a square of side $2\alpha$
centered about its original position. (By periodicity, if the move would put the
particle outside the square, it re-enters from the opposite side.)

We then calculate the change in energy $\Delta E$ caused by the move.

- If $\Delta E < 0$ (move brings the system to lower energy), we allow the move
  and put the particle in its new position.
- If $\Delta E > 0$, we allow the move with probability $\exp(-\Delta E/kT)$:
  take a random number $\xi_3$ between 0 and 1; if $\xi_3 < \exp(-\Delta E/kT)$
  move the particle to its new position; if $\xi_3 > \exp(-\Delta E/kT)$ return
  it to its old position.

Then, whether the move has been allowed or not (whether we are in a new
configuration or the original one), we consider that we are in a new
configuration for the purpose of taking averages:

<!-- eq:4 -->
$$ \bar{F} = (1/M) \sum_{j=1}^{M} F_j, $$

where $F_j$ is the value of property `F` after the `j`th move is carried out
according to the complete prescription above. Having attempted to move a particle
we proceed similarly with the next one.

### Derivation (not implemented): proof that the scheme samples $\exp(-E/kT)$

Since a particle is allowed to move to any point within a square of side
$2\alpha$ with finite probability, a large enough number of moves lets it reach
any point in the complete square. Since this holds for all particles, we may
reach any point in configuration space. Hence the method is **ergodic**. (In
practice one need not make enough moves for a particle to diffuse evenly, since
configuration space is symmetric with respect to interchange of particles.)

Consider a very large ensemble of systems. Suppose there are a finite number of
states of the system, and that $\nu_r$ is the number of systems of the ensemble
in state `r` (a state = a given point in configuration space). We must prove that
after many moves the ensemble tends to a distribution

<!-- eq:nu-canonical -->
$$ \nu_r \propto \exp(-E_r/kT). $$

Make a move in all systems of the ensemble. Let the a priori probability that the
move carries a system from state `r` to state `s` be $P_{rs}$ (the probability
before discriminating on $\exp(-\Delta E/kT)$). First, $P_{rs} = P_{sr}$, since a
particle is equally likely to be moved anywhere within a square of side $2\alpha$
centered about its original position: if states `r` and `s` differ only by the
position of the moved particle and these positions are within each other's
squares, the transition probabilities are equal; otherwise they are zero.

Assume $E_r > E_s$. The number of systems moving from `r` to `s` is $\nu_r P_{rs}$
(all moves to lower energy are allowed). The number moving from `s` to `r` is
$\nu_s P_{sr} \exp(-(E_r - E_s)/kT)$ (weighted by the exponential factor). Thus
the net number of systems moving from `s` to `r` is

<!-- eq:5 -->
$$ P_{rs}\left(\nu_s \exp(-(E_r - E_s)/kT) - \nu_r\right). $$

So between any two states `r` and `s`, if

<!-- eq:6 -->
$$ (\nu_r/\nu_s) > \left[ \exp(-E_r/kT)/\exp(-E_s/kT) \right], $$

then on the average more systems move from state `r` to state `s`. Combined with
ergodicity, these facts mean the ensemble must approach the canonical
distribution. It is clear from the derivation that after a forbidden move we must
count again the initial configuration; not doing so would unjustifiably reduce
the number in state `s` relative to `r`.

The argument does not specify how rapidly the canonical distribution is
approached. The maximum displacement $\alpha$ must be chosen with care: if too
large, most moves are forbidden; if too small, the configuration changes too
slowly. In either case it takes longer to reach equilibrium.

For the rigid-sphere case the game of chance on $\exp(-\Delta E/kT)$ is
unnecessary since $\Delta E$ is either zero or infinity. Particles are moved one
at a time per Eq. (3); if a sphere overlaps another after a move, it is returned
to its original position.

**Note on random numbers:** the random numbers used were generated by the middle
square process: if $\xi_n$ is an `m`-digit random number, a new random number
$\xi_{n+1}$ is the middle `m` digits of the complete 2m-digit square of $\xi_n$.

## III. Specialization to Rigid Spheres in Two Dimensions

### A. The Equation of State

The virial theorem of Clausius gives an equation of state in terms of $\bar{n}$,
the average density of other particles at the surface of a particle. Let
$X_i^{(\text{tot})}$ and $X_i^{(\text{int})}$ represent the total and internal
force, respectively, acting on particle `i` at position $\mathbf{r}_i$. The virial
theorem can be written

<!-- eq:7 -->
$$ \left\langle \sum_{i} \mathbf{X}_{i}^{(\text{tot})} \cdot \mathbf{r}_{i} \right\rangle_{\text{AV}} = 2PA + \left\langle \sum_{i} \mathbf{X}_{i}^{(\text{int})} \cdot \mathbf{r}_{i} \right\rangle_{\text{AV}} = 2E_{\text{kin}}. $$

Here `P` is the pressure, `A` the area, and $E_{\text{kin}}$ the total kinetic
energy $E_{\text{kin}} = Nm\bar{v}^2/2$ of the system of `N` particles.

### Derivation (not implemented): virial contribution of rigid-sphere collisions

Consider the collisions of the spheres as represented by those of a particle of
radius $d_0$ (twice the radius of the actual spheres) surrounded by $\bar{n}$
point particles per unit area. Surrounding particles in an area of
$2\pi d_0 v \cos\phi\, \Delta t$, traveling with velocity `v` at an angle $\phi$
with the radius vector, collide with the central particle provided
$|\phi| < \pi/2$. Assuming elastic recoil, they each exert an average force during
time $\Delta t$ on the central particle of $2mv\cos\phi/\Delta t$.

All $\phi$'s are equally probable, since for any velocity-independent potential
the velocity distribution is Maxwellian, hence isotropic. The total force acting
on the central particle, averaged over $\phi$, over time, and over velocity, is

<!-- eq:8 -->
$$ \bar{F}_{i} = m\bar{v}^{2}\pi d_{0}\bar{n}. $$

The internal-force sum is

<!-- eq:internal-sum -->
$$ \left\langle \sum_{i} \mathbf{X}_{i}^{(\text{int})} \cdot \mathbf{r}_{i} \right\rangle_{\text{AV}} = -\frac{1}{2} \sum_{i} \left\{ \sum_{j \neq i} \mathbf{r}_{ij} F_{ij} \right\}, $$

with $F_{ij}$ the magnitude of the force between two particles and $r_{ij}$ the
distance between them. Since $r_{ij} = d_0$ and $\sum_j F_{ij}$ is given by
Eq. (8), we have

<!-- eq:9 -->
$$ \left\langle \sum_{i} \mathbf{X}_{i}^{(\text{int})} \cdot \mathbf{r}_{i} \right\rangle_{\text{AV}} = -(Nm\bar{v}^{2}/2)\pi d_{0}^{2}\bar{n}. $$

Substituting (9) into (7) and replacing $(N/2)m\bar{v}^2$ by $E_{\text{kin}}$
gives finally

<!-- eq:10 -->
$$ PA = E_{\text{kin}}(1 + \pi d_0^2 \bar{n}/2) \equiv NkT(1 + \pi d_0^2 \bar{n}/2). $$

This shows that a determination of the single quantity $\bar{n}$ (per Eq. (4)) as
a function of the area `A` suffices to determine the equation of state for rigid
spheres.

### B. The Actual Calculation of $\bar{n}$

We set up the calculation on a system composed of $N=224$ particles
$(i=0, 1, \dots, 223)$ placed inside a square of unit side and unit area. The
particles were arranged initially in a trigonal lattice of fourteen particles per
row by sixteen particles per column, alternate rows displaced relative to each
other. This gives each particle six nearest neighbors at approximately equal
distances of $d = 1/14$.

Instead of performing the calculation for various areas `A` and a fixed distance
$d_0$, we solve the equivalent problem of leaving $A=1$ fixed and changing $d_0$.
We denote by $A_0$ the area the particles occupy in close-packed arrangement. For
numerical convenience we defined an auxiliary parameter $\nu$, varied from zero to
seven, in terms of which the ratio $(A/A_0)$ and the forbidden distance $d_0$ are
defined as:

<!-- eq:11a -->
$$ d_0 = d(1 - 2^{\nu - 8}), \qquad d = (1/14), $$

<!-- eq:11b -->
$$ (A/A_0) = 1/(3^{\frac{1}{2}} d_0^2 N/2) = 1/[0.98974329(1-2^{\nu-8})^2]. $$

The unit cell is a parallelogram with interior angle 60°, side $d_0$, and altitude
$3^{\frac{1}{2}} d_0/2$ in the close-packed system.

Every configuration reached was analyzed in terms of a radial distribution
function $N(r^2)$. We chose a $K > 1$ for each $\nu$ and divided the area between
$\pi d_0^2$ and $K^2 \pi d_0^2$ into sixty-four zones of equal area $\Delta A^2$:

<!-- eq:deltaA2 -->
$$ \Delta A^2 = (K^2 - 1)\pi d_0^2/64. $$

We then had the machine calculate for each configuration the number of pairs of
particles $N_m$ $(m = 1, 2, \dots, 64)$ separated by distances `r` which satisfy

<!-- eq:12 -->
$$ (m-1)\Delta A^2 + \pi d_0^2 < \pi r^2 \le m\Delta A^2 + \pi d_0^2. $$

The $N_m$ were averaged over successive configurations per Eq. (4), and after
every sixteen cycles (a cycle = moving every particle once) were extrapolated
back to $r^2 = d_0^2$ to obtain $N_{\frac{1}{2}}$. This $N_{\frac{1}{2}}$ differs
from $\bar{n}$ in Eq. (10) by a constant factor depending on `N` and `K`.

The quantity `K` was chosen for each $\nu$ to give reasonable statistics for the
$N_m$. The maximum displacement $\alpha$ of Eq. (3) was set to $(d - d_0)$. About
half the moves in a cycle were forbidden by this choice, and the initial approach
to equilibrium from the regular lattice was fairly rapid.

## IV. Numerical Results for Rigid Spheres in Two Dimensions

We first ran for something less than sixteen cycles to remove the effects of the
initial regular configuration on the averages. Then about forty-eight to
sixty-four cycles were run at $\nu = 2, 4, 5, 5.5, 6, 6.25, 6.5,$ and $7$. A
smaller amount of data was obtained at $\nu = 0, 1,$ and $3$. The time per cycle
on the Los Alamos MANIAC is approximately three minutes, and a given point on the
pressure curve was obtained in four to five hours of running.

Deviation from the free volume theory begins with a fairly sudden break at
$\nu = 6\ (A/A_0 \sim 1.8)$.

For $\nu = 5$, least-square fits with a straight line to the first sixteen $N_m$
values gave extrapolated values $N_{\frac{1}{2}}^{(1)} = 6367$,
$N_{\frac{1}{2}}^{(2)} = 6160$, and $N_{\frac{1}{2}}^{(3)} = 6377$; the average was
used in constructing PA/NkT. In general, least-square fits of the first sixteen to
twenty $N_m$'s by a parabola, or a straight line where suitable, were made. The
resultant value of (PA/NkT)-1 is $64\bar{N}_{\frac{1}{2}}/[N^2(K^2-1)]$. The
average error was about 3 percent (root-mean-square deviations for the three or
four $N_{\frac{1}{2}}$ values).

Table I gives the results in numerical form (see `checks.md`).

## V. The Virial Coefficient Expansion

One can show that

<!-- eq:13 -->
$$ (PA/NkT)-1 = C_1(A_0/A) + C_2(A_0/A)^2 + C_3(A_0/A)^3 + C_4(A_0/A)^4 + O(A_0/A)^5, $$

with coefficients

<!-- eq:13-coeffs -->
$$ C_1 = \pi/3^{\frac{1}{2}}, \qquad C_2 = 4\pi^2 A_{3,3}/9, $$

$$ C_3 = \pi^3 (6A_{4,5} - 3A_{4,4} - A_{4,6})/3^{\frac{3}{2}}, $$

$$ C_{4} = (8\pi^{3}/135) \cdot [12A_{5,5} - 60A_{5,6}' - 10A_{5,6}'' + 30A_{5,7}' + 60A_{5,7}'' + 10A_{5,7}''' - 30A_{5,8}' - 15A_{5,8}'' + 10A_{5,9} - A_{5,10}]. $$

The coefficients $A_{i,k}$ are cluster integrals over configuration space of `i`
particles with `k` bonds between them. A bond is established if the two particles
overlap. The cluster integral is the volume of configuration space for which the
appropriate bonds are established. If `k` bonds can be distributed over the `i`
particles in two or more different ways without destroying irreducibility, the
separate cases are distinguished by primes.

Define $f(r_{ij})$ by

<!-- eq:f-def -->
$$ f(r_{ij}) = 1 \text{ if } r_{ij} < d, \qquad f(r_{ij}) = 0 \text{ if } r_{ij} > d, $$

then, for example,

<!-- eq:A33 -->
$$ A_{3,3} = \frac{1}{\pi^2 d^4} \int \cdots \int dx_1 dx_2 dx_3\, dy_1 dy_2 dy_3\, (f_{12} f_{23} f_{31}). $$

The coefficients $A_{3,3}$, $A_{4,4}$, and $A_{4,5}$ were calculated
algebraically, the remainder numerically by Monte Carlo integration. For example
for $A_{5,5}$: particle 1 was placed at the origin, and particles 2, 3, 4, 5 were
put down at random subject to $f_{12} = f_{23} = f_{34} = f_{15} = 1$. The number
of trials for which $f_{45} = 1$, divided by the total number of trials, is just
$A_{5,5}$.

The data on $A_{4,6}$ is reliable:

<!-- eq:A46-ratio -->
$$ A_{4,6}/A_{4,4} = 0.752\ (\pm 0.002). $$

Because of the relatively large positive and negative terms in $C_4$ of Eq. (13),
$C_4$ (a small difference) is less accurate:

<!-- eq:C4-value -->
$$ C_4 = 8\pi^3 (0.585)/135 \quad (\pm \sim 5\ \text{percent}). $$

Our final formula is

<!-- eq:14 -->
$$ (PA/NkT) - 1 = 1.813799(A_0/A) + 2.57269(A_0/A)^2 + 3.179(A_0/A)^3 + 3.38(A_0/A)^4 + O(A_0/A)^5. $$

This formula agrees very well with our calculated equation of state for
$(A/A_0) > 2.5$.

## VI. Conclusion

The method of Monte Carlo integrations over configuration space seems a feasible
approach to statistical-mechanical problems not yet analytically soluble. At least
for a single-phase system a sample of several hundred particles seems sufficient.
Runs made with 56 particles and with 224 particles agreed within statistical
error. For a computing time of a few hours it seems possible to obtain the
pressure for a given volume and temperature to an accuracy of a few percent.

For two-dimensional rigid spheres our results agree with the free volume
approximation for $A/A_0 < 1.8$ and with a five-term virial expansion for
$A/A_0 > 2.5$. There is no indication of a phase transition.
