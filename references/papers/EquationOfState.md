![](_page_0_Picture_0.jpeg)

### Equation of State Calculations by Fast Computing Machines

[Nicholas Metropolis](http://jcp.aip.org/search?sortby=newestdate&q=&searchzone=2&searchtype=searchin&faceted=faceted&key=AIP_ALL&possible1=Nicholas Metropolis&possible1zone=author&alias=&displayid=AIP&ver=pdfcov), [Arianna W. Rosenbluth,](http://jcp.aip.org/search?sortby=newestdate&q=&searchzone=2&searchtype=searchin&faceted=faceted&key=AIP_ALL&possible1=Arianna W. Rosenbluth&possible1zone=author&alias=&displayid=AIP&ver=pdfcov) [Marshall N. Rosenbluth,](http://jcp.aip.org/search?sortby=newestdate&q=&searchzone=2&searchtype=searchin&faceted=faceted&key=AIP_ALL&possible1=Marshall N. Rosenbluth&possible1zone=author&alias=&displayid=AIP&ver=pdfcov) [Augusta H. Teller,](http://jcp.aip.org/search?sortby=newestdate&q=&searchzone=2&searchtype=searchin&faceted=faceted&key=AIP_ALL&possible1=Augusta H. Teller&possible1zone=author&alias=&displayid=AIP&ver=pdfcov) and [Edward Teller](http://jcp.aip.org/search?sortby=newestdate&q=&searchzone=2&searchtype=searchin&faceted=faceted&key=AIP_ALL&possible1=Edward Teller&possible1zone=author&alias=&displayid=AIP&ver=pdfcov)

Citation: [J. Chem. Phys. 2](http://jcp.aip.org/?ver=pdfcov)1, 1087 (1953); doi: 10.1063/1.1699114

View online: [http://dx.doi.org/10.1063/1.1699114](http://link.aip.org/link/doi/10.1063/1.1699114?ver=pdfcov)

View Table of Contents: [http://jcp.aip.org/resource/1/JCPSA6/v21/i6](http://jcp.aip.org/resource/1/JCPSA6/v21/i6?ver=pdfcov)

Published by the [AIP Publishing LLC.](http://www.aip.org/?ver=pdfcov)

### Additional information on J. Chem. Phys.

Journal Homepage: [http://jcp.aip.org/](http://jcp.aip.org/?ver=pdfcov)

Journal Information: [http://jcp.aip.org/about/about\\_the\\_journal](http://jcp.aip.org/about/about_the_journal?ver=pdfcov) Top downloads: [http://jcp.aip.org/features/most\\_downloaded](http://jcp.aip.org/features/most_downloaded?ver=pdfcov)

Information for Authors: [http://jcp.aip.org/authors](http://jcp.aip.org/authors?ver=pdfcov)

![](_page_0_Picture_12.jpeg)

instead, only water molecules with different amounts of excitation energy. These may follow any of three paths:

- (a) The excitation energy is lost without dissociation into radicals (by collision, or possibly radiation, as in aromatic hydrocarbons).
- (b) The molecules dissociate, but the resulting radicals recombine without escaping from the liquid cage.
- (c) The molecules dissociate and escape from the cage. In this case we would not expect them to move more than a few molecular diameters through the dense medium before being thermalized.

In accordance with the notation introduced by Burton, Magee, and Samuel,<sup>22</sup> the molecules following

paths (a) and (b) can be designated  $H_2O^*$  and those following path (c) can be designated  $H_2O^{\dagger}$ . It seems reasonable to assume for the purpose of these calculations that the ionized  $H_2O$  molecules will become the  $H_2O^{\dagger}$  molecules, but this is not likely to be a complete correspondence.

In conclusion we would like to emphasize that the qualitative result of this section is not critically dependent on the exact values of the physical parameters used. However, this treatment is classical, and a correct treatment must be wave mechanical; therefore the result of this section cannot be taken as an a priori theoretical prediction. The success of the radical diffusion model given above lends some plausibility to the occurrence of electron capture as described by this crude calculation. Further work is clearly needed.

THE JOURNAL OF CHEMICAL PHYSICS

VOLUME 21, NUMBER 6

JUNE, 1953

#### Equation of State Calculations by Fast Computing Machines

NICHOLAS METROPOLIS, ARIANNA W. ROSENBLUTH, MARSHALL N. ROSENBLUTH, AND AUGUSTA H. TELLER,

Los Alamos Scientific Laboratory, Los Alamos, New Mexico

ANI

EDWARD TELLER,\* Department of Physics, University of Chicago, Chicago, Illinois (Received March 6, 1953)

A general method, suitable for fast computing machines, for investigating such properties as equations of state for substances consisting of interacting individual molecules is described. The method consists of a modified Monte Carlo integration over configuration space. Results for the two-dimensional rigid-sphere system have been obtained on the Los Alamos MANIAC and are presented here. These results are compared to the free volume equation of state and to a four-term virial coefficient expansion.

### I. INTRODUCTION

HE purpose of this paper is to describe a general method, suitable for fast electronic computing machines, of calculating the properties of any substance which may be considered as composed of interacting individual molecules. Classical statistics is assumed, only two-body forces are considered, and the potential field of a molecule is assumed spherically symmetric. These are the usual assumptions made in theories of liquids. Subject to the above assumptions, the method is not restricted to any range of temperature or density. This paper will also present results of a preliminary twodimensional calculation for the rigid-sphere system. Work on the two-dimensional case with a Lennard-Jones potential is in progress and will be reported in a later paper. Also, the problem in three dimensions is being investigated.

# II. THE GENERAL METHOD FOR AN ARBITRARY POTENTIAL BETWEEN THE PARTICLES

In order to reduce the problem to a feasible size for numerical work, we can, of course, consider only a finite number of particles. This number N may be as high as several hundred. Our system consists of a square† containing N particles. In order to minimize the surface effects we suppose the complete substance to be periodic, consisting of many such squares, each square containing N particles in the same configuration. Thus we define  $d_{AB}$ , the minimum distance between particles A and B, as the shortest distance between A and any of the particles B, of which there is one in each of the squares which comprise the complete substance. If we have a potential which falls off rapidly with distance, there will be at most one of the distances AB which can make a substantial contribution; hence we need consider only the minimum distance  $d_{AB}$ .

<sup>&</sup>lt;sup>22</sup> Burton, Magee, and Samuel, J. Chem. Phys. 20, 760 (1952).

<sup>\*</sup> Now at the Radiation Laboratory of the University of California, Livermore, California.

<sup>†</sup> We will use the two-dimensional nomenclature here since it is easier to visualize. The extension to three dimensions is obvious.

Our method in this respect is similar to the cell method except that our cells contain several hundred particles instead of one. One would think that such a sample would be quite adequate for describing any onephase system. We do find, however, that in two-phase systems the surface between the phases makes quite a perturbation. Also, statistical fluctuations may be sizable.

If we know the positions of the N particles in the square, we can easily calculate, for example, the potential energy of the system,

$$E = \frac{1}{2} \sum_{\substack{i=1 \ i \neq j}}^{N} \sum_{j=1}^{N} V(d_{ij}). \tag{1}$$

(Here V is the potential between molecules, and  $d_{ij}$  is the minimum distance between particles i and j as defined above.)

In order to calculate the properties of our system we use the canonical ensemble. So, to calculate the equilibrium value of any quantity of interest F,

$$\bar{F} = \left[ \int F \exp(-E/kT) d^{2N} p d^{2N} q \right] / \left[ \int \exp(-E/kT) d^{2N} p d^{2N} q \right], \quad (2)$$

where  $(d^{2n} p d^{2n} q)$  is a volume element in the 4N-dimensional phase space. Moreover, since forces between particles are velocity-independent, the momentum integrals may be separated off, and we need perform only the integration over the 2N-dimensional configuration space. It is evidently impractical to carry out a several hundred-dimensional integral by the usual numerical methods, so we resort to the Monte Carlo method. The Monte Carlo method for many-dimensional integrals consists simply of integrating over a random sampling of points instead of over a regular array of points.

Thus the most naive method of carrying out the integration would be to put each of the N particles at a random position in the square (this defines a random point in the 2N-dimensional configuration space), then calculate the energy of the system according to Eq. (1), and give this configuration a weight  $\exp(-E/kT)$ . This method, however, is not practical for close-packed configurations, since with high probability we choose a configuration where  $\exp(-E/kT)$  is very small; hence a configuration of very low weight. So the method we employ is actually a modified Monte Carlo scheme, where, instead of choosing configurations randomly, then weighting them with  $\exp(-E/kT)$ , we choose

configurations with a probability  $\exp(-E/kT)$  and weight them evenly.

This we do as follows: We place the N particles in any configuration, for example, in a regular lattice. Then we move each of the particles in succession according to the following prescription:

$$\begin{array}{l}
X \to X + \alpha \xi_1 \\
Y \to Y + \alpha \xi_2,
\end{array} \tag{3}$$

where  $\alpha$  is the maximum allowed displacement, which for the sake of this argument is arbitrary, and  $\xi_1$  and  $\xi_2$ are random numbers between (-1) and 1. Then, after we move a particle, it is equally likely to be anywhere within a square of side  $2\alpha$  centered about its original position. (In accord with the periodicity assumption. if the indicated move would put the particle outside the square, this only means that it re-enters the square from the opposite side.)

We then calculate the change in energy of the system  $\Delta E$ , which is caused by the move. If  $\Delta E < 0$ , i.e., if the move would bring the system to a state of lower energy, we allow the move and put the particle in its new position. If  $\Delta E > 0$ , we allow the move with probability  $\exp(-\Delta E/kT)$ ; i.e., we take a random number  $\xi_3$  between 0 and 1, and if  $\xi_3 < \exp(-\Delta E/kT)$ , we move the particle to its new position. If  $\xi_3$  $> \exp(-\Delta E/kT)$ , we return it to its old position. Then, whether the move has been allowed or not, i.e., whether we are in a different configuration or in the original configuration, we consider that we are in a new configuration for the purpose of taking our averages. So

$$\bar{F} = (1/M) \sum_{j=1}^{M} F_j,$$
 (4)

where  $F_j$  is the value of the property F of the system after the jth move is carried out according to the complete prescription above. Having attempted to move a particle we proceed similarly with the next one.

We now prove that the method outlined above does choose configurations with a probability  $\exp(-E/kT)$ . Since a particle is allowed to move to any point within a square of side  $2\alpha$  with a finite probability, it is clear that a large enough number of moves will enable it to reach any point in the complete square. || Since this is true of all particles, we may reach any point in configuration space. Hence, the method is ergodic.

Next consider a very large ensemble of systems. Suppose for simplicity that there are only a finite number of states¶ of the system, and that  $\nu_r$  is the number of

<sup>‡</sup> This method has been proposed independently by J. E. Mayer and by S. Ulam. Mayer suggested the method as a tool to deal with the problem of the liquid state, while Ulam proposed it as a procedure of general usefulness. B. Alder, J. Kirkwood, S. Frankel, and V. Lewinson discussed an application very similar to ours.

<sup>§</sup> It might be mentioned that the random numbers that we used were generated by the middle square process. That is, if  $\xi^u$ is an m digit random number, then a new random number  $\xi_{n+1}$ is given as the middle m digits of the complete 2m digit square of  $\xi_n$ .

In practice it is, of course, not necessary to make enough moves to allow a particle to diffuse evenly throughout the system since configuration space is symmetric with respect to interchange of particles.

¶ A state here means a given point in configuration space.

systems of the ensemble in state r. What we must prove is that after many moves the ensemble tends to a distribution

$$\nu_r \propto \exp(-E_r/kT)$$
.

Now let us make a move in all the systems of our ensemble. Let the a priori probability that the move will carry a system in state r to state s be  $P_{rs}$ . [By the a priori probability we mean the probability before discriminating on  $\exp(-\Delta E/kT)$ .] First, it is clear that  $P_{rs} = P_{sr}$ , since according to the way our game is played a particle is equally likely to be moved anywhere within a square of side  $2\alpha$  centered about its original position. Thus, if states r and s differ from each other only by the position of the particle moved and if these positions are within each other's squares, the transition probabilities are equal; otherwise they are zero. Assume  $E_r > E_s$ . Then the number of systems moving from state r to state s will be simply  $\nu_r P_{rs}$ , since all moves to a state of lower energy are allowed. The number moving from s to r will be  $\nu_s P_{sr} \exp(-(E_r - E_s)/kT)$ , since here we must weigh by the exponential factor. Thus the net number of systems moving from s to r is

$$P_{rs}(\nu_s \exp(-(E_r - E_s)/kT) - \nu_r).$$
 (5)

So we see that between any two states r and s, if

$$(\nu_r/\nu_s) > \lceil \exp(-E_r/kT)/\exp(-E_s/kT) \rceil,$$
 (6)

on the average more systems move from state r to state s. We have seen already that the method is ergodic; i.e., that any state can be reached from any other, albeit in several moves. These two facts mean that our ensemble must approach the canonical distribution. It is, incidentally, clear from the above derivation that after a forbidden move we must count again the initial configuration. Not to do this would correspond in the above case to removing from the ensemble those systems which tried to move from s to r and were forbidden. This would unjustifiably reduce the number in state s relative to r.

The above argument does not, of course, specify how rapidly the canonical distribution is approached. It may be mentioned in this connection that the maximum displacement  $\alpha$  must be chosen with some care; if too large, most moves will be forbidden, and if too small, the configuration will not change enough. In either case it will then take longer to come to equilibrium.

For the rigid-sphere case, the game of chance on  $\exp(-\Delta E/kT)$  is, of course, not necessary since  $\Delta E$  is either zero or infinity. The particles are moved, one at a time, according to Eq. (3). If a sphere, after such a move, happens to overlap another sphere, we return it to its original position.

# III. SPECIALIZATION TO RIGID SPHERES IN TWO DIMENSIONS

#### A. The Equation of State

The virial theorem of Clausius can be used to give an equation of state in terms of  $\bar{n}$ , the average den-

![](_page_3_Picture_14.jpeg)

Fig. 1. Collisions of rigid spheres.

sity of other particles at the surface of a particle. Let  $X_i^{\text{(tot)}}$  and  $X_i^{\text{(int)}}$  represent the total and the internal force, respectively, acting on particle i, at a position  $\mathbf{r}_i$ . Then the virial theorem can be written

$$\langle \sum_{i} \mathbf{X}_{i}^{(\text{tot})} \cdot \mathbf{r}_{i} \rangle_{\text{AV}} = 2PA + \langle \sum_{i} \mathbf{X}_{i}^{(\text{int})} \cdot \mathbf{r}_{i} \rangle_{\text{AV}} = 2E_{\text{kin}}.$$
 (7)

Here P is the pressure, A the area, and  $E_{kin}$  the total kinetic energy,

$$E_{\rm kin} = Nm\bar{v}^2/2$$

of the system of N particles.

Consider the collisions of the spheres for convenience as represented by those of a particle of radius  $d_0$ , twice the radius of the actual spheres, surrounded by  $\bar{n}$  point particles per unit area. Those surrounding particles in an area of  $2\pi d_0 v \cos\phi \Delta t$ , traveling with velocity v at an angle  $\phi$  with the radius vector, collide with the central particle provided  $|\phi| < \pi/2$ . (See Fig. 1.) Assuming elastic recoil, they each exert an average force during the time  $\Delta t$  on the central particle of

#### $2mv\cos\phi/\Delta t$ .

One can see that all  $\phi$ 's are equally probable, since for any velocity-independent potential between particles the velocity distribution will just be Maxwellian, hence isotropic. The total force acting on the central particle, averaged over  $\phi$ , over time, and over velocity, is

$$\bar{F}_{i} = m\bar{v}^{2}\pi d_{0}\bar{n}. \tag{8}$$
The sum
$$\langle \sum_{i} \mathbf{X}_{i}^{(\mathrm{int})} \cdot \mathbf{r}_{i} \rangle_{\mathsf{AV}}$$
\nis
$$-\frac{1}{2} \sum_{i} \{ \sum_{j} \mathbf{r}_{ij} F_{ij} \},$$

$$\sum_{i \neq i} \mathbf{r}_{ij} F_{ij} \},$$

with  $F_{ij}$  the magnitude of the force between two particles and  $r_{ij}$  the distance between them. We see that

![](_page_4_Figure_1.jpeg)

Fig. 2. Initial trigonal lattice.

 $r_{ij} = d_0$  and  $\sum_{j} F_{ij}$  is given by Eq. (8), so we have

$$\langle \sum_{i} \mathbf{X}_{i}^{(\text{int})} \cdot \mathbf{r}_{i} \rangle_{Av} = -(Nm\bar{v}^{2}/2)\pi d_{0}^{2}\bar{n}. \tag{9}$$

Substitution of (9) into (7) and replacement of  $(N/2)m\bar{v}^2$  by  $E_{\rm kin}$  gives finally

$$PA = E_{\text{kin}}(1 + \pi d_0^2 \bar{n}/2) \equiv NkT(1 + \pi d_0^2 \bar{n}/2).$$
 (10)

This equation shows that a determination of the one quantity  $\bar{n}$ , according to Eq. (4) as a function of A, the area, is sufficient to determine the equation of state for the rigid spheres.

#### B. The Actual Calculation of $\overline{n}$

We set up the calculation on a system composed of N=224 particles  $(i=0,\,1\cdots 223)$  placed inside a square of unit side and unit area. The particles were arranged initially in a trigonal lattice of fourteen particles per row by sixteen particles per column, alternate rows being displaced relative to each other as shown in Fig. 2. This arrangement gives each particle six nearest neighbors at approximately equal distances of d=1/14 from it.

Instead of performing the calculation for various areas A and for a fixed distance  $d_0$ , we shall solve the equivalent problem of leaving A=1 fixed and changing  $d_0$ . We denote by  $A_0$  the area the particles occupy in close-packed arrangement (see Fig. 3). For numerical convenience we defined an auxiliary parameter  $\nu$ , which we varied from zero to seven, and in terms of which the ratio  $(A/A_0)$  and the forbidden distance  $d_0$  are defined as follows:

$$d_0 = d(1 - 2^{\nu - 8}), \quad d = (1/14),$$
 (11a)

$$(A/A_0) = 1/(3\frac{1}{2}d_0^2N/2) = 1/0.98974329(1-2^{\nu-8})^2$$
. (11b)

![](_page_4_Figure_13.jpeg)

Fig. 3. The close-packed arrangement for determining  $A_0$ .

The unit cell is a parallelogram with interior angle 60°, side  $d_0$ , and altitude  $3^{\frac{1}{2}}d_0/2$  in the close-packed system.

Every configuration reached by proceeding according to the method of the preceding section was analyzed in terms of a radial distribution function  $N(r^2)$ . We chose a K>1 for each  $\nu$  and divided the area between  $\pi d_0^2$  and  $K^2\pi d_0^2$  into sixty-four zones of equal area  $\Delta A^2$ ,

$$\Delta A^2 = (K^2 - 1)\pi d_0^2/64$$
.

We then had the machine calculate for each configuration the number of pairs of particles  $N_m$   $(m=1, 2, \cdots 64)$ separated by distances r which satisfy

$$(m-1)\Delta A^2 + \pi d_0^2 < \pi r^2 \le m\Delta A^2 + \pi d_0^2. \tag{12}$$

The  $N_m$  were averaged over successive configurations according to Eq. (4), and after every sixteen cycles (a cycle consists of moving every particle once) were extrapolated back to  $r^2=d_0^2$  to obtain  $N_{\frac{1}{2}}$ . This  $N_{\frac{1}{2}}$  differs from  $\bar{n}$  in Eq. (10) by a constant factor depending on N and K.

The quantity K was chosen for each  $\nu$  to give reasonable statistics for the  $N_m$ . It would, of course, have been possible by choosing fairly large K's, with perhaps a larger number of zones, to obtain  $N(r^2)$  at large distances. The oscillatory behavior of  $N(r^2)$  at large distances is of some interest. However, the time per cycle goes up fairly rapidly with K and with the number of zones in the distance analysis. For this reason only the behavior of  $N(r^2)$  in the neighborhood of  $d_0^2$  was investigated.

The maximum displacement  $\alpha$  of Eq. (3) was set to  $(d-d_0)$ . About half the moves in a cycle were forbidden by this choice, and the initial approach to equilibrium from the regular lattice was fairly rapid.

![](_page_5_Figure_2.jpeg)

FIG. 4. A plot of (PA/NkT)-1 versus  $(A/A_0)-1$ . Curve A (solid line) gives the results of this paper. Curves B and C (dashed and dot-dashed lines) give the results of the free volume theory and of the first four virial coefficients, respectively.

## IV. NUMERICAL RESULTS FOR RIGID SPHERES IN TWO DIMENSIONS

We first ran for something less than sixteen cycles in order to get rid of the effects of the initial regular configuration on the averages. Then about forty-eight to sixty-four cycles were run at

$$\nu = 2, 4, 5, 5.5, 6, 6.25, 6.5,$$
 and 7.

Also, a smaller amount of data was obtained at  $\nu=0$ , 1, and 3. The time per cycle on the Los Alamos MANIAC is approximately three minutes, and a given point on the pressure curve was obtained in four to five hours of running. Figure 4 shows (PA/NkT)-1 versus  $(A/A_0)-1$  on a log-log scale from our results (curve A), compared to the free volume equation of Wood¹ (curve B) and to the curve given by the first four virial coefficients (curve C). The last two virial coefficients were obtained by straightforward Monte Carlo integration on the MANIAC (see Sec. V). It is seen that the agreement between curves A and B at small areas and between curves A and C at large areas is good. Deviation from the free volume theory begins with a fairly sudden break at  $\nu=6(A/A_0 \sim 1.8)$ .

A sample plot of the radial distribution function for  $\nu=5$  is given in Fig. 5. The various types of points represent values after sixteen, thirty-two, and forty-eight cycles. For  $\nu=5$ , least-square fits with a straight line to the first sixteen  $N_m$  values were made, giving extrapolated values of  $N_{\frac{1}{2}}{}^{(1)}=6367$ ,  $N_{\frac{1}{2}}{}^{(2)}=6160$ , and  $N_{\frac{1}{2}}{}^{(3)}=6377$ . The average of these three was used in constructing PA/NkT. In general, least-square fits of the first sixteen to twenty  $N_m$ 's by means of a parabola, or, where it seemed suitable, a straight line, were made.

![](_page_5_Figure_10.jpeg)

Fig. 5. The radial distribution function  $N_m$  for  $\nu=5$ ,  $(A/A_0)=1.31966$ , K=1.5. The average of the extrapolated values of  $N_{\frac{1}{2}}$  in  $\overline{N}_{\frac{1}{2}}=6301$ . The resultant value of (PA/NkT)-1 is  $64\overline{N}_{\frac{1}{2}}/N^2(K^2-1)$  or 6.43. Values after 16 cycles,  $\bullet$ ; after 32,  $\times$ ; and after 48,  $\odot$ .

The errors indicated in Fig. 4 are the root-mean-square deviations for the three or four  $N_{\frac{1}{2}}$  values. Our average error seemed to be about 3 percent.

Table I gives the results of our calculations in numerical form. The columns are  $\nu$ ,  $A/A_0$ , (PA/NkT)-1, and, for comparison purposes, (PA/NkT-1) for the free volume theory and for the first four coefficients in the virial coefficient expansion, in that order, and finally  $PA_0/NkT$  from our results.

#### V. THE VIRIAL COEFFICIENT EXPANSION

One can show<sup>2</sup> that

$$(PA/NkT)-1=C_1(A_0/A)+C_2(A_0/A)^2$$
  
+ $C_3(A_0/A)^3+C_4(A_0/A)^4+0(A_0/A)^5$ ,

$$C_1 = \pi/3^{\frac{1}{2}}, \quad C_2 = 4\pi^2 A_{3,3}/9,$$

$$C_3 = \pi^3 (6A_{4.5} - 3A_{4.4} - A_{4.6})/3^{\frac{3}{2}}$$

$$C_{4} = (8\pi^{3}/135) \cdot [12A_{5,5} - 60A_{5,6}' - 10A_{5,6}'' + 30A_{5,7}' + 60A_{5,7}'' + 10A_{5,7}''' - 30A_{5,8}' - 15A_{5,8}'' + 10A_{5,9} - A_{5,10}].$$
(13)

Table I. Results of this calculation for  $(PA/NkT)-1=X_1$  compared to the free volume theory  $(X_2)$  and the four-term virial expansion  $(X_3)$ . Also  $(PA_0/NkT)$  from our calculations.

| ν    | $(A/A_0)$ | $X_1$  | $X_2$ | $X_{2}$ | $(PA_0/NkT)$ |
|------|-----------|--------|-------|---------|--------------|
| 2    | 1.04269   | 49.17  | 47.35 | 9.77    | 48.11        |
| 4    | 1.14957   | 13.95  | 13.85 | 7.55    | 13.01        |
| 5    | 1.31966   | 6.43   | 6.72  | 5.35    | 5.63         |
| 5.5  | 1.4909    | 4.41   | 4.53  | 4.02    | 3.63         |
| 6    | 1.7962    | 2,929  | 2.939 | 2.680   | 2.187        |
| 6.25 | 2.04616   | 2.186  | 2.323 | 2.065   | 1.557        |
| 6.5  | 2.41751   | 1.486  | 1.802 | 1.514   | 1.028        |
| 7    | 4.04145   | 0.6766 | 0.990 | 0.667   | 0.4149       |

<sup>&</sup>lt;sup>2</sup> J. E. Mayer and M. G. Mayer, Statistical Mechanics (John Wiley and Sons, Inc., New York, 1940), pp. 277-291.

<sup>&</sup>lt;sup>1</sup> William W. Wood, J. Chem. Phys. 20, 1334 (1952).

![](_page_6_Picture_1.jpeg)

Fig. 6. Schematic diagrams for the various area integrals.

The coefficients  $A_{i,k}$  are cluster integrals over configuration space of i particles, with k bonds between them. In our problem a bond is established if the two particles overlap. The cluster integral is the volume of configuration space for which the appropriate bonds are established. If k bonds can be distributed over the i particles in two or more different ways without destroying the irreducibility of the integrals, the separate cases are distinguished by primes. For example,  $A_{33}$  is given schematically by the diagram

and mathematically as follows: if we define  $f(r_{ij})$  by

$$f(r_{ij}) = 1$$
 if  $r_{ij} < d$ ,  
 $f(r_{ij}) = 0$  if  $r_{ij} > d$ ,

then

$$A_{3,3} = \frac{1}{\pi^2 d^4} \int \cdots \int dx_1 dx_2 dx_3 dy_1 dy_2 dy_3 (f_{12} f_{23} f_{31}).$$

The schematics for the remaining integrals are indicated in Fig. 6.

The coefficients  $A_{3,3}$ ,  $A_{4,4}$ , and  $A_{4,5}$  were calculated algebraically, the remainder numerically by Monte Carlo integration. That is, for  $A_{5,5}$  for example, particle 1 was placed at the origin, and particles 2, 3, 4, and 5

were put down at random, subject to  $f_{12}=f_{23}=f_{34}=f_{15}=1$ . The number of trials for which  $f_{45}=1$ , divided by the total number of trials, is just  $A_{5,5}$ .

The data on  $A_{4,6}$  is quite reliable. We obtained

$$A_{4,6}/A_{4,4} = 0.752(\pm 0.002)$$
.

However, because of the relatively large positive and negative terms in  $C_4$  of Eq. (13), the coefficient  $C_4$ , being a small difference, is less accurate. We obtained

$$C_4 = 8\pi^3 (0.585)/135$$
 ( $\pm \sim 5$  percent).

Our final formula is

$$(PA/NkT) - 1 = 1.813799(A_0/A) + 2.57269(A_0/A)^2 + 3.179(A_0/A)^3 + 3.38(A_0/A)^4 + 0(A_0/A)^5.$$
(14)

This formula is plotted in curve C of Fig. 4 and tabulated for some values of  $(A/A_0)$  in column 5 of Table I. It is seen in Fig. 4 that the curves agrees very well with our calculated equation of state for  $(A/A_0) > 2.5$ . In this region both the possible error in our last virial coefficients and the contribution of succeeding terms in the expansion are quite small (less than our probable statistical error) so that the virial expansion should be accurate.

#### VI. CONCLUSION

The method of Monte Carlo integrations over configuration space seems to be a feasible approach to statistical mechanical problems which are as yet not analytically soluble. At least for a single-phase system a sample of several hundred particles seems sufficient. In the case of two-dimensional rigid spheres, runs made with 56 particles and with 224 particles agreed within statistical error. For a computing time of a few hours with presently available electronic computers, it seems possible to obtain the pressure for a given volume and temperature to an accuracy of a few percent.

In the case of two-dimensional rigid spheres our results are in agreement with the free volume approximation for  $A/A_0 < 1.8$  and with a five-term virial expansion for  $A/A_0 > 2.5$ . There is no indication of a phase transition.

Work is now in progress for a system of particles with Lennard-Jones type interactions and for three-dimensional rigid spheres.