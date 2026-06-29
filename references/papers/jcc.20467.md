![](_page_0_Picture_0.jpeg)

# **Quantum Mechanical Calculation of the Effects of Stiff and Rigid Constraints in the Conformational Equilibrium of the Alanine Dipeptide**

# **PABLO ECHENIQUE,1,2 IVÁN CALVO,1,2 J. L. ALONSO1,2**

*1Departamento de Física Teórica, Facultad de Ciencias, Universidad de Zaragoza, Pedro Cerbuna 12, 50009 Zaragoza, Spain 2Instituto de Biocomputación y Física de los Sistemas Complejos (BIFI), Edificio Cervantes, Corona de Aragón 42, 50009 Zaragoza, Spain*

*Received 25 January 2006; Revised 2 March 2006; Accepted 15 March 2006 DOI 10.1002/jcc.20467 Published online 9 August 2006 in Wiley InterScience (www.interscience.wiley.com).*

**Abstract:** If constraints are imposed on a macromolecule, two inequivalent classical models may be used: the stiff and the rigid one. This work studies the effects of such constraints on the conformational equilibrium distribution (CED) of the model dipeptide HCO-l-Ala-NH2 *without any simplifying assumption*. We use *ab initio* quantum mechanics calculations including electron correlation at the MP2 level to describe the system, and we measure the conformational dependence of all the correcting terms to the naive CED based in the potential energy surface that appear when the constraints are considered. These terms are related to mass-metric tensors determinants and also occur in the Fixman's compensating potential. We show that some of the corrections are non-negligible if one is interested in the whole Ramachandran space. On the other hand, if only the energetically lower region, containing the principal secondary structure elements, is assumed to be relevant, then, all correcting terms may be neglected up to peptides of considerable length. This is the first time, as far as we know, that the analysis of the conformational dependence of these correcting terms is performed in a relevant biomolecule with a realistic potential energy function.

© 2006 Wiley Periodicals, Inc. J Comput Chem 27: 1733–1747, 2006

**Key words:** constraints; alanine dipeptide; Fixman; *ab initio*; mass metric tensors

# **Introduction**

In computer simulations of large complex systems, such as macromolecules and, specially, proteins,1–6 one of the main bottlenecks to design efficient algorithms is the necessity to sample an astronomically large conformational space.3, 7 In addition, being the typical timescales of the different movements in a wide range, demandingly small timesteps must be used in molecular dynamics simulations in order to properly account for the fastest modes, which lie in the femtosecond range. However, most of the biological interesting behavior (allosteric transitions, protein folding, enzymatic catalysis) is related to the slowest conformational changes, which occur in the timescale of milliseconds or even seconds.4, 8–11 Fortunately, the fastest modes are also the most energetic ones and are rarely activated at room temperature. Therefore, in order to alleviate the computational problems and also simplify the images used to think about these elusive systems, one may naturally consider the reduction of the number of degrees of freedom describing macromolecules via the imposition of constraints.<sup>12</sup>

How to study the conformational equilibrium of these constrained systems has been an object of much debate.13–17 Two different classical models exist in the literature, which are conceptually13–16, 18, 19 and practically6, 13, 20–24 inequivalent. In the *classical rigid* model, the constraints are assumed to be *exact* and all the velocities that are orthogonal to the hypersurface defined by them vanish. In the *classical stiff* ∗ model, on the other hand, the constraints are assumed to be *approximate* and they are implemented by a steep potential that drives the system to the constrained hypersurface. In this case, the orthogonal velocities are activated and may act as "heat containers."

*Correspondence to:* P. Echenique; e-mail: pnique@unizar.es

Contract/grant sponsor: Aragón Government

Contract/grant sponsor: MEC (Spain) (to P.E., I.C.); contract/grant numbers: FIS2004-05073, FPA2003-02948

Contract/grant sponsor: MYCT (Spain); contract/grant number: BFM2003- 08532

<sup>∗</sup>Some authors use the word *flexible* to refer to this model.15, 21, 22, 25 We, however, prefer to term it *stiff* <sup>18</sup> and keep the name *flexible* to refer to the case in which no constraints are imposed.

In this work, we do not address the question of which model is a better approximation of physical reality. Although, in the literature, it is commonly assumed (often implicitly) that the classical stiff model should be taken as a reference, 6,9,16,19,20,22,26 we believe that this opinion is much influenced by the use of popular classical force fields<sup>6,27-37</sup> (which are stiff by construction) and by the goal of reproducing their results at a lower computational cost, i.e., using rigid molecular dynamics simulations. 4,5,8,9,14,19,21–23,25,26,38–42 In our opinion, the question whether the rigid or the stiff model should be used to approximate the real quantum mechanical statistics of an arbitrary organic molecule has not been satisfactorily answered yet. For discussions about the topic, see ref. 13–15, 17, 18 and 43–45. In this work, we adopt the cautious position that any of the two models may be useful in certain cases or for certain purposes and we study them both on equal footing. Our concern is, then, to study the effects that either way of imposing constraints causes in the conformational equilibrium of macromolecules.

In the Born-Oppenheimer approximation  $^{46}$  customarily used in quantum mechanics and in the majority of the classical force fields, the relevant degrees of freedom are the Euclidean (also called *Cartesian* by some authors) 3n coordinates of the n nuclei. However, it is frequent to define a different set of coordinates in which the overall translation and rotation of the system are distinguished and the remaining 3n-6 degrees of freedom are chosen (according to different prescriptions as *internal coordinates*, which are simple geometrical parameters (typically consisting of bond lengths, bond angles, and dihedral angles) that describe the internal structure of the system.  $^{47}$ 

In macromolecules, the natural constraints are those derived from the relative rigidity of the internal covalent structure of groups of atoms that share a common center (and also from the rigidity of rotation around double or triple bonds) compared to the energetically "cheaper" rotation around single bonds. In internal coordinates, these chemical constraints may be directly implemented by asking that some conveniently selected hard coordinates (normally, bond lengths, bond angles, and some dihedrals) have constant values or values that depend on the remaining soft coordinates (see ref. 15 for a definition). In Euclidean coordinates, on the other hand, the expression of the constraints is more cumbersome and complicated procedures<sup>25,26,40,48–50</sup> must be used at each timestep to implement them in molecular dynamics simulations. This is why, in the classical stiff model, as well as in the rigid one, it is common to use internal coordinates and they are also the choice throughout this work.

In the equilibrium statistical mechanics of both the stiff and rigid models, the marginal probability density in the coordinate part of the phase space in these internal coordinates is not proportional to the naive  $\exp[-\beta V_{\Sigma}(q^i)]$ , where  $V_{\Sigma}(q^i)$  denotes the potential energy on the constrained hypersurface<sup>†</sup>. Instead, some correcting terms that come from different sources must be added to the potential energy  $V_{\Sigma}(q^i)$ .  $^{13,15,18,19,39,51,52}$  These terms involve determinants of massmetric tensors and also of the Hessian matrix of the constraining part of the potential (see sec. 2). If Monte Carlo simulations in the coordinate space are to be performed  $^{5,53-57}$  and the probability densities

that correspond to any of these two models sampled, the corrections should be included or, otherwise, showed to be negligible.

Additionally, the three different correcting terms are involved in the definition of the so-called Fixman's compensating potential, <sup>16</sup> which is frequently used to reproduce the stiff equilibrium distribution using rigid molecular dynamics simulations. <sup>9,14,19,21–23,38,39,42,51</sup>

Customarily in the literature, some of these corrections to the potential energy are assumed to be independent of the conformation and thus dropped from the basic expressions. Also, subtly entangled to the assumptions underlying many classical results, a second type of approximation is made that consists of assuming that the equilibrium values of the hard coordinates do not depend on the soft coordinates

In this work, we measure the conformational dependence of *all correcting terms* and of the Fixman's compensating potential in the model dipeptide HCO-L-Ala- $NH_2$  without any simplifying assumption. The potential energy function is considered to be the effective Born-Oppenheimer potential for the nuclei derived from *ab initio* quantum mechanical calculations including electron correlation at the MP2 level. We also repeat the calculations, with the same basis set (6-31++G(d,p)) and at the Hartree-Fock level of the theory in order to investigate if this less demanding method without electron correlation may be used in further studies. It is *the first time*, as far as we are aware, that this type of study is performed in a relevant biomolecule with a realistic potential energy function.

In the next section, we introduce the notation to be used and derive the statistical mechanics formulae of the rigid and stiff models in the general case. In the Methods section, we describe the computational methods used and we summarize the factorization of the external coordinates presented in ref. 58. The Results section is devoted to the presentation and discussion of the assessment of the approximation that consists of neglecting the different corrections to the potential energy in the model dipeptide HCO-L-Ala-NH<sub>2</sub>, without any simplifying assumption, which is the central aim of this work, and is the following section conclusions are summarized. Finally, in the appendix, we discuss the use of the different approximations in the literature and we give a precise definition of exactly and approximately separable hard and soft coordinates that will shed some light on the relation between the different types of simplifications aforementioned.

### Theory

First of all, it is convenient to introduce certain notational conventions that will be used extensively in the rest of the work:

- The system under scrutiny will be a set of n mass points termed *atoms*. The Euclidean coordinates of the atom  $\alpha$  in a set of axes fixed in space are denoted by  $\vec{x}_{\alpha}$ . The subscript  $\alpha$  runs from 1 to n.
- The curvilinear coordinates suitable to describe the system will be denoted by  $q^{\mu}$ ,  $\mu = 1, ..., 3n$  and the set of Euclidean coordinates by  $x^{\mu}$  when no explicit reference to the atoms index needs to be made. We shall often use N := 3n for the total number of degrees of freedom.
- The coordinates  $q^{\mu}$  are split into  $(q^A, q^a)$ , a = 7, ..., N. The first six are termed *external coordinates* and are denoted by  $q^A$ . They

 $<sup>^\</sup>dagger$ By  $q^i$ , we denote the soft internal coordinates of the system. See the next section and the Appendix for a precise definition.

Table 1. Definition of the Indices Used.

| Indices                        | Range              | Number        | Description               |
|--------------------------------|--------------------|---------------|---------------------------|
| $\alpha, \beta, \gamma, \dots$ | $1,\ldots,n$       | n             | Atoms                     |
| $\mu, \nu, \rho, \dots$        | $1, \ldots, N$     | N = 3n        | All coordinates           |
| $A, B, C, \dots$               | 1,,6               | 6             | External coordinates      |
| $a, b, c, \dots$               | $7,\ldots,N$       | N-6           | Internal coordinates      |
| $i, j, k, \dots$               | $7, \ldots, M + 6$ | M             | Soft internal coordinates |
| $I,J,K,\dots$                  | $M+7,\ldots,N$     | L = N - M - 6 | Hard internal coordinates |
| $u, v, w, \dots$               | $1, \dots, M + 6$  | M + 6         | All soft coordinates      |

describe the overall position and orientation of the system with respect to a frame fixed in space (see ref. 58 for further details). The coordinates  $q^a$  are said *internal coordinates* and determine the positions of the atoms in the frame fixed in the system. They parameterize what we shall call the *internal subspace* or *conformational space*, denoted by  $\mathcal{I}$  and the coordinates  $q^A$  parameterize the *external subspace*, denoted by  $\mathcal{E}$ .

- The general set-up of the problem may be described as follows. Instead of us being interested on the conformational equilibrium of the system in the external subspace  $\mathcal{E}$  plus the whole internal subspace  $\mathcal{I}$  (i.e., the *whole space*, denoted by  $\mathcal{E} \times \mathcal{I}$ ), we wish to find the probability density on a hypersurface  $\Sigma \subset \mathcal{I}$  of dimension M (plus the external subspace  $\mathcal{E}$ ), i.e., on  $\mathcal{E} \times \Sigma$ .
- In typical internal coordinates  $q^a$ , normally consisting of bond lengths, bond angles, and dihedral angles (see ref. 59 and references therein), the hypersurface  $\Sigma$  is described via L:=N-M-6 constraints:

$$q^{I} = f^{I}(q^{i}) \quad I = M + 7, \dots, N,$$
 (1)

where the  $q^a$  are split into  $q^a \equiv (q^i, q^I)$ , and the  $q^i$ , i = 7, ..., M+6, which parameterize  $\Sigma$ , are called *internal soft coordinates*, whereas the  $q^I$  are termed *hard coordinates*. The external coordinates  $q^A$ , together with the  $q^i$ , form the whole set of *soft coordinates*, denoted by  $q^u \equiv (q^A, q^i)$ , u = 1, ..., M+6.

In Table 1, a summary of the indices used is given.

## Classical Stiff Model

In the classical stiff model, the constraints in eq. (1) are implemented by imposing an strong energy penalization when the internal conformation of the system, described by  $q^a$ , departs from the constrained hypersurface  $\Sigma$ . To ensure this, we must have that the potential energy function in  $\mathcal I$  satisfies certain conditions. First, we write the potential  $V(q^a)$  as follows<sup>‡</sup>:

$$V(q^{i}, q^{I}) = \underbrace{V(q^{i}, f^{I}(q^{i}))}_{V_{\Sigma}(q^{i})} + \underbrace{\left[V(q^{i}, q^{I}) - V(q^{i}, f^{I}(q^{i}))\right]}_{V_{c}(q^{i}, q^{I})}. \tag{2}$$

Next, we impose the following conditions on the *constraining* potential  $V_c(q^i, q^I)$  defined above:

- (i) That  $V_c(q^i, f^l(q^i)) \le V_c(q^i, q^l) \, \forall q^i, q^l$ , i.e., that  $\Sigma$  be the global minimum of  $V_c$  (and, henceforth, a local one too) with respect to variations of the hard coordinates.
- (ii) That, for small variations  $\Delta q^I$  on the hard coordinates (i.e., for changes  $\Delta q^I$  considered as physically irrelevant), the associated changes in  $V_c(q^i,q^I)$  are much larger than the thermal energy RT.

The advantages of this formulation, much similar to that on ref. 15, are many. First, it sets a convenient framework for the derivation of the statistical mechanics formulae of the classical stiff model relating it to the fully flexible model in the whole space  $\mathcal{E} \times \mathcal{I}$ . Second, it clearly separates the potential energy on  $\Sigma$  from the part that is responsible of implementing the constraints. Third, contrarily to the formulation based on delta functions, <sup>51</sup> it allows to clearly understand the necessity of including the correcting term associated to the determinant of the Hessian of  $V_c$  (see the derivation that follows). Finally, and more importantly for us, it provides a direct prescription for calculating  $V_{\Sigma}(q^i)$  and  $\Sigma$  (the Potential Energy Surface (PES), frequently used in quantum chemistry calculations  $^{60-64}$ ) via geometry optimization at fixed values of the soft coordinates.

We also remark that, in order to satisfy point (ii) above and to allow the derivation of the different correcting terms that follows and the validity of the final expressions, the hard coordinates  $q^I$  must be indeed hard; however, the *soft coordinates*  $q^i$  do not have to be soft (in the sense that they produce energetic changes much smaller than RT when varied). They may be interesting for some other reason and hence voluntarily picked to describe the system studied, without altering the formulae presented in this section. Despite this qualifications, the terms *soft* and *hard* will be kept in this work for consistence with most of the existing literature,  $^{15, 18, 52, 65, 66}$  although, in some cases, the labels *important* and *unimportant* (for  $q^i$  and  $q^I$  respectively), proposed by Karplus and Kushick,  $^{67}$  may be more appropriate.

In the case of the model dipeptide HCO-L-Ala-NH<sub>2</sub> investigated in this work, for example, the barriers in the Ramachandran angles  $\phi$  and  $\psi$  may be as large as  $\sim$ 40 RT; however, the study of small dipeptides is normally aimed to the design of effective potentials for polypeptides, <sup>68–70</sup> where long-range interactions in the sequence may compensate these local energy penalizations. This and the fact that the Ramachandran angles are the relevant degrees of freedom to describe the conformation of the backbone of these systems make it convenient to choose them as *soft coordinates*  $q^i$  despite the fact that they may be energetically hard in the case of the dipeptide HCO-L-Ala-NH<sub>2</sub>. As remarked above, this does not affect the calculations

Now, because of the condition (ii) above, the statistical weights of the conformations that lie far away from the constrained hypersurface  $\Sigma$  are negligible and, therefore, it suffices to describe the system in the vicinity of the equilibrium values of the  $q^I$ . In this region, for each value of the internal soft coordinates  $q^i$ , we may expand  $V_c(q^i, q^I)$  in eq. (2) up to second-order in the hard

<sup>&</sup>lt;sup>‡</sup>Note that we have simply added and subtracted from the total potential energy  $V(q^i, q^I) \equiv V(q^a)$  of the system the same quantity,  $V(q^i, f^I(q^i))$ .

coordinates around  $\Sigma$  (i.e., around  $q^I=f^I(q^i)$ ) and drop the higher order terms:

$$V_{c}(q^{i}, q^{I}) \simeq V_{c}(q^{i}, f^{I}(q^{i})) + \left[\frac{\partial V_{c}}{\partial q^{J}}\right]_{\Sigma} (q^{J} - f^{J}(q^{i}))$$

$$+ \frac{1}{2} \underbrace{\left[\frac{\partial^{2} V_{c}}{\partial q^{J} \partial q^{K}}\right]_{\Sigma}}_{\mathcal{H}_{JK}(q^{i})} (q^{J} - f^{J}(q^{i})) (q^{K} - f^{K}(q^{i})), \quad (3)$$

where the subindex  $\Sigma$  indicates evaluation on the constrained hypersurface and a more compact notation,  $\mathcal{H}(q^i)$ , has been introduced for the Hessian matrix of  $V_c$  with respect to the hard variables evaluated on  $\Sigma$ . Also, the Einstein's sum convention is assumed on repeated indices.

In this expression, the zeroth order term  $V_c(q^i, f^I(q^i))$  is zero by definition of  $V_c$  (see eq. (2)) and the linear term is also zero, because of the condition (i) above. Hence, the first non-zero term of the expansion in eq. (3) is the second order one. Using this, together with eq. (2), we may write the *stiff Hamiltonian* 

$$H_{s}(q^{\mu}, p_{\mu}) := \frac{1}{2} p_{\nu} G^{\nu\rho}(q^{u}, q^{I}) p_{\rho} + V_{\Sigma}(q^{i}) + \frac{k}{2} \mathcal{H}_{JK}(q^{i}) (q^{J} - f^{J}(q^{i})) (q^{K} - f^{K}(q^{i})), \quad (4)$$

the mass-metric tensor  $G_{\nu\rho}$  being

$$G_{\nu\rho}(q^{u}, q^{I}) := \sum_{\sigma=1}^{N} \frac{\partial x^{\sigma}(q^{\mu})}{\partial q^{\nu}} m_{\sigma} \frac{\partial x^{\sigma}(q^{\mu})}{\partial q^{\rho}}$$
 (5)

and  $G^{\nu\rho}$  its inverse, defined by

$$G^{\nu\sigma}(q^u, q^I) G_{\sigma\rho}(q^u, q^I) = \delta^{\nu}_{\rho}, \tag{6}$$

where  $\delta^{\nu}_{\rho}$  denotes the Kronecker's delta.

Therefore, the *stiff partition function* of the system is<sup>§</sup>

$$Z_{\rm s} = \frac{\alpha_{QM}}{h^N} \int dq^{\mu} dp_{\mu} \exp\left[-\beta H_{\rm s}(q^{\mu}, p_{\mu})\right], \tag{7}$$

where h is Planck's constant, we denote  $\beta := 1/RT$  (per mole energy units are used throughout the article, so RT is preferred over  $k_BT$ ) and  $\alpha_{QM}$  is a combinatorial number that accounts for quantum indistinguishability and that must be specified in each particular case (e.g., for a gas of N indistinguishable particles,  $\alpha_{QM} = 1/N!$ ).

Now, using the condition (ii) again, the  $q^I$  appearing in the massmetric tensor G in  $H_s$  (in eq. (7)) can be approximately evaluated at their equilibrium values  $f^I(q^i)$ , yielding, for the stiff partition function.

$$Z_{s} = \frac{\alpha_{QM}}{h^{N}} \int dq^{u} dq^{I} dp_{\mu} \exp \left[ -\beta \left( \frac{1}{2} p_{v} G^{v\rho}(q^{u}, f^{I}(q^{i})) p_{\rho} + V_{\Sigma}(q^{i}) + \frac{1}{2} \mathcal{H}_{JK}(q^{i}) (q^{J} - f^{J}(q^{i})) (q^{K} - f^{K}(q^{i})) \right) \right].$$
(8)

If we now integrate over the hard coordinates  $q^I$ , we have

$$Z_{s} = \left(\frac{2\pi}{\beta}\right)^{\frac{L}{2}} \frac{\alpha_{QM}}{h^{N}} \int dq^{u} dp_{\mu} \exp\left[-\beta \left(\frac{1}{2}p_{\nu}G^{\nu\rho}(q^{u}, f^{I}(q^{i}))p_{\rho}\right) + V_{\Sigma}(q^{i}) + T\frac{R}{2}\ln[\det\mathcal{H}(q^{i})]\right]$$
(9)

where the part of the result of the Gaussian integral consisting of  $\det^{-1/2}\mathcal{H}$  has been taken to the exponent.

Note that the Hessian matrix  $\mathcal{H}_{JK}$  involves only derivatives with respect to the hard coordinates (see eq. (3)), so that the minimization protocol embodied in eq. (2) (which is identical to the procedure followed in quantum chemistry for computing the PES along reaction coordinates) guarantees that  $\mathcal{H}_{JK}$  is positive defined and, hence, det  $\mathcal{H}$  is positive, allowing to take its logarithm as in the previous expression. The fact that it is only this "partial Hessian" that makes sense in the computation of equilibrium properties along soft (or reaction) coordinates, has been recently pointed out in ref. 72.

It is also frequent to integrate over the momenta in the partition function. Doing this in eq. (9) and taking the determinant of the mass-metric tensor that shows up<sup>||</sup> to the exponent, we may write the partition function as an integral only on the coordinates:

$$Z_{s} = \chi_{s}(T) \int dq^{u} \exp\left[-\beta \left(V_{\Sigma}(q^{i}) + T\frac{R}{2}\ln[\det \mathcal{H}(q^{i})]\right) - T\frac{R}{2}\ln\left[\det G(q^{u}, f^{I}(q^{i}))\right]\right)\right], \quad (10)$$

where the multiplicative factor that depends on T has been defined as follows:

$$\chi_{s}(T) := \left(\frac{2\pi}{\beta}\right)^{\frac{N+L}{2}} \frac{\alpha_{QM}}{h^{N}}.$$
 (11)

If the exponent in eq. (10) is seen as a free energy, then,  $V_{\Sigma}(q^i)$  may be regarded as the internal energy and the two conformation-dependent correcting terms that are added to it as effective entropies (which is compatible with their being linear in RT). The second one comes only from the desire to write the marginal probabilities in the

<sup>§</sup> No Jacobian appears in the integral measure because  $q^{\mu}$  and  $p_{\mu}$  are obtained from the Euclidean coordinates via a canonical transformation.<sup>71</sup>

Note that, by G, we denote the matrix that corresponds to the mass-metric tensor with two covariant indices  $G_{\mu\nu}$ . The same convention has been followed for the Hessian matrix  $\mathcal{H}$  in eq. (9) and for the reduced mass-metric tensor g in eq. (21).

coordinate space (i.e., averaging the momenta) and may be called a *kinetic entropy*, <sup>17</sup> the first term, on the other hand, is truly an entropic term that comes from the averaging out of certain degrees of freedom and it is reminiscent of the *conformational* or *configurational entropies* appearing in quasiharmonic analysis. <sup>6,67,73</sup>

In this spirit, we define

$$F_{s}(q^{u}) := V_{\Sigma}(q^{i}) - T(S_{s}^{c}(q^{i}) + S_{s}^{k}(q^{u})),$$
 (12a)

$$S_{\rm s}^c(q^i) := -\frac{R}{2} \ln[\det \mathcal{H}(q^i)], \tag{12b}$$

$$S_{\rm s}^{\rm k}(q^u) := \frac{R}{2} \ln \left[ \det G(q^u, f^I(q^i)) \right].$$
 (12c)

In such a way that the *stiff equilibrium probability* in the soft subspace  $\mathcal{E} \times \Sigma$  is given by

$$P_{s}(q^{u}) = \frac{\exp[-\beta F_{s}(q^{u})]}{Z'_{s}}, \quad \text{with} \quad Z'_{s} := \int dq^{u} \exp[-\beta F_{s}(q^{u})].$$
(13)

Now, it is worth remarking that, although the kinetic entropy  $S_s^k$  depends on the external coordinates  $q^A$ , we have recently shown<sup>58</sup> that the determinant of the mass-metric tensor G may be written, for any molecule, general internal coordinates and arbitrary constraints, as a product of two functions: one depending only on the external coordinates, and the other only on the internal ones  $q^a$ . Hence the externals-dependent factor in eq. (12c) may be integrated out independently to yield an effective free energy and a probability density  $P_s$  that depend only on the soft internals  $q^i$  (see sec. Factorization of the External Coordinates).

#### Classical Rigid Model

If the relations in eq. (1) are considered to hold *exactly* and are treated as holonomic constraints, the Hamiltonian function that describes the classical mechanics in the subspace  $(\mathcal{E} \times \Sigma) \subset (\mathcal{E} \times \mathcal{I})$ , spanned by the coordinates  $g^u$ , may be written as follows:

$$H_{\rm r}(q^u, \eta_u) := \frac{1}{2} \eta_{\nu} g^{\nu w}(q^u) \eta_w + V_{\Sigma}(q^i), \tag{14}$$

where the reduced mass-metric tensor  $g_{vw}(q^u)$  in  $\mathcal{E} \times \Sigma$ , that appears in the kinetic energy, is (see what follows)

$$g_{vw}(q^{u}) = G_{vw}(q^{u}, f^{I}(q^{i})) + \frac{\partial f^{J}(q^{i})}{\partial q^{v}} G_{JK}(q^{u}, f^{I}(q^{i})) \frac{\partial f^{K}(q^{i})}{\partial q^{w}}$$

$$+ G_{vK}(q^{u}, f^{I}(q^{i})) \frac{\partial f^{K}(q^{i})}{\partial q^{w}} + \frac{\partial f^{J}(q^{i})}{\partial q^{v}} G_{Jw}(q^{u}, f^{I}(q^{i}))$$

$$:= \frac{\partial \tilde{f}^{\mu}}{\partial q^{\nu}} G_{\mu\nu}(q^{u}, f^{I}(q^{i})) \frac{\partial \tilde{f}^{\nu}}{\partial q^{w}}, \quad (15)$$

and  $g^{vw}(q^u)$  is defined to be its inverse in the sense of eq. (6). Also, the notation

$$\tilde{f}^{\mu} := \begin{cases} q^{u} & \text{if } u := \mu = 1, \dots, M + 6 \\ f^{I}(q^{i}) & \text{if } I := \mu = M + 7, \dots, N \end{cases}$$
 (16)

has been introduced for convenience.

Note that eq. (15) may derived from the unconstrained Hamiltonian in  $(\mathcal{E} \times \mathcal{I})$ ,

$$H(q^{\mu}, p_{\mu}) := \frac{1}{2} p_{\nu} G^{\nu \rho}(q^{\mu}) p_{\rho} + V(q^{a}), \tag{17}$$

using the constraints in eq. (1), together with its time derivatives (denoted by an overdot: as in  $\dot{A}$ )

$$\dot{q}^I := \frac{\partial f^I(q^i)}{\partial q^j} \dot{q}^j \tag{18}$$

and defining the momenta  $\eta_{\nu}$  as

$$\eta_{\nu} := g_{\nu w}(q^{u}) \, \dot{q}^{w} = g_{\nu w}(q^{u}) \, G^{w \mu}(q^{u}, f^{I}(q^{i})) \, p_{\mu}. \tag{19}$$

Hence, the rigid partition function is

$$Z_{\rm r} = \frac{\alpha_{QM}}{h^{M+6}} \int dq^u \, d\eta_u \exp\left[-\beta \left(\frac{1}{2}\eta_{\nu}g^{\nu\nu}(q^u)\eta_{\nu} + V_{\Sigma}(q^i)\right)\right]. \tag{20}$$

Integrating over the momenta, we obtain the marginal probability density in the coordinate space analogous to eq. (10):

$$Z_{\rm r} = \chi_{\rm r}(T) \int dq^u \exp\left[-\beta \left(V_{\Sigma}(q^i) - T\frac{R}{2}\ln[\det g(q^u)]\right)\right], \tag{21}$$

where

$$\chi_{\rm r}(T) := \left(\frac{2\pi}{\beta}\right)^{\frac{M+6}{2}} \frac{\alpha_{QM}}{h^{\frac{M+6}{2}}}.$$
(22)

Repeating the analogy with free energies and entropies in the last paragraphs of the previous subsection, we define

$$F_{\rm r}(q^u) := V_{\Sigma}(q^i) - TS_{\rm r}^{\rm k}(q^u),$$
 (23a)

$$S_{\mathbf{r}}^{\mathbf{k}}(q^{u}) := \frac{R}{2} \ln[\det g(q^{u})], \tag{23b}$$

being the *rigid equilibrium probability* in the soft subspace  $\mathcal{E} \times \Sigma$ 

$$P_{r}(q^{u}) = \frac{\exp[-\beta F_{r}(q^{u})]}{Z'_{r}}, \quad \text{with} \quad Z'_{r} := \int dq^{u} \exp[-\beta F_{r}(q^{u})].$$
(24)

**Table 2.** Equilibrium Probability Densities and Correcting Terms to the Potential Energy  $V_{\Sigma}(q^i)$  in the Classical Stiff and Rigid Models of Constraints.

| Classical stiff model                                                            | Classical rigid model                                                        |
|----------------------------------------------------------------------------------|------------------------------------------------------------------------------|
| $P_{\rm s}(q^{\rm u}) = \frac{\exp[-\beta F_{\rm s}(q^{\rm u})]}{Z_{\rm s}'}$    | $P_{\rm r}(q^u) = \frac{\exp[-\beta F_{\rm r}(q^u)]}{Z_{\rm r}'}$            |
| $F_{s}(q^{u}) := V_{\Sigma}(q^{i}) - T(S_{s}^{c}(q^{i}) + S_{s}^{k}(q^{u}))$     | $F_{\mathbf{r}}(q^u) := V_{\Sigma}(q^i) - TS_{\mathbf{r}}^{\mathbf{k}}(q^u)$ |
| $S_{s}^{k}(q^{u}) := \frac{R}{2} \ln \left[ \det G(q^{u}, f^{I}(q^{i})) \right]$ | $S_{\mathbf{r}}^{\mathbf{k}}(q^u) := \frac{R}{2} \ln[\det g(q^u)]$           |
| $S_{\rm s}^c(q^i) := -\frac{R}{2} \ln[\det \mathcal{H}(q^i)]$                    |                                                                              |

As in the case of G, we have shown in ref. 58 that the determinant of the reduced mass-metric tensor g may be written, for any molecule, general internal coordinates and arbitrary constraints, as a product of two functions: one depending only on the external coordinates, and the other only on the internal ones  $q^i$ . Hence the externals-dependent factor in det  $g(q^u)$  may be integrated out independently to yield a free energy and a probability density  $P_r$  that depend only on the soft internals  $q^i$  (see the next section).

To end this subsection, we remark that it is frequent in the literature  $^{9,18,19,21-23,38,42,51,57}$  to define the so-called *Fixman's compensating potential*  $^{16}$  as the difference between  $F_s(q^u)$ , in eq. (12), and  $F_r(q^u)$ , defined above, i.e.,

$$V_{\mathrm{F}}(q^{u}) := TS_{\mathrm{r}}^{\mathrm{k}}(q^{u}) - TS_{\mathrm{s}}^{c}(q^{i}) - TS_{\mathrm{s}}^{\mathrm{k}}(q^{u})$$

$$= \frac{RT}{2} \ln \left[ \frac{\det G(q^{u})}{\det \mathcal{H}(q^{i}) \det g(q^{u})} \right]. \tag{25}$$

Hence, performing rigid molecular dynamics simulations, which would yield an equilibrium distribution proportional to  $\exp[-\beta F_{\rm r}(q^u)]$ , and adding  $V_{\rm F}(q^u)$  to the potential energy  $V_{\Sigma}(q^i)$ , one can reproduce instead the stiff probability density  $P_{\rm s} \propto \exp[-\beta F_{\rm s}(q^u)]$ . It is allows to obtain at a lower computational cost (due to the timescale problems discussed in the introduction) equilibrium averages that otherwise must be extracted from expensive fully flexible whole-space simulations. In fact, it seems that this particular application of the theoretical tools herein described, and not the search for the correct probability density to sample in Monte Carlo simulations, was what prompted the interest in the study of mass-metric tensors effects.

Finally, in Table 2, we summarize the equilibrium probability densities and the different correcting terms derived in this section.

# Methods

# Factorization of the External Coordinates

In the recent work,  $^{58}$  we have shown that the determinant of the mass-metric tensor G in eq. (12c) can be written as follows if

the SASMIC<sup>59</sup> coordinates for general branched molecules are used:

$$\det G = \left(\prod_{\alpha=1}^{n} m_{\alpha}^{3}\right) \sin^{2} \theta \left(\prod_{\alpha=2}^{n} r_{\alpha}^{4}\right) \left(\prod_{\alpha=3}^{n} \sin^{2} \theta_{\alpha}\right), \tag{26}$$

where the  $r_{\alpha}$  are bond lengths and the  $\theta_{\alpha}$  bond angles.

Note that this expression, whose validity was proved for the more particular case of serial polymers by  $G\bar{o}$  and Scheraga<sup>15</sup> and, before, by Volkenstein,<sup>74</sup> does not explicitly depend on the dihedral angles. However, it may depend on them via the hard coordinates if the constraints in the form presented in eq. (1) are used.

The term depending on the masses of the atoms in the expression above may be dropped from eq. (12c), because it does not depend on the conformation, and the only part of det G that depend on the external coordinates,  $\sin^2 \theta$ , may be integrated out in eq. (10) ( $\theta$  is one of the externals  $q^A$  that describe the overall orientation of the molecule; see ref. 58 for further details). Hence, the kinetic entropy due to the mass-metric tensor G in the stiff case, may be written, up to additive constants, as

$$S_{\rm s}^{\rm k}(q^i) = \frac{R}{2} \left[ \sum_{\alpha=2}^n \ln(r_{\alpha}^4) + \sum_{\alpha=3}^n \ln(\sin^2 \theta_{\alpha}) \right],$$
 (27)

where the individual contributions of each degree of freedom have been factorized

Also in ref. 58, we have shown that the determinant of the reduced mass-metric tensor g in eq. (23b) can be written as follows:

$$\det g = \sin^2 \theta \, \det g_2(q^i), \tag{28}$$

being the matrix  $g_2$ 

$$g_{2} = \begin{pmatrix} m_{\text{tot}} I^{(3)} & m_{\text{tot}} v(\vec{R}) & \cdots & m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{j}} & \cdots \\ m_{\text{tot}} v^{T}(\vec{R}) & \mathcal{J} & \cdots & \sum_{\alpha} m_{\alpha} \frac{\partial \vec{x}_{\alpha}'}{\partial q^{j}} \times \vec{x}_{\alpha}' & \cdots \\ \vdots & \vdots & & \vdots & & \vdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \sum_{\alpha} m_{\alpha} \left( \frac{\partial \vec{x}_{\alpha}'}{\partial q^{i}} \times \vec{x}_{\alpha}' \right)^{\mathsf{T}} & \cdots & \sum_{\alpha} m_{\alpha} \frac{\partial \vec{x}_{\alpha}'^{T}}{\partial q^{i}} \frac{\partial \vec{x}_{\alpha}'}{\partial q^{j}} & \cdots \\ \vdots & & \vdots & & \vdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \sum_{\alpha} m_{\alpha} \left( \frac{\partial \vec{x}_{\alpha}'}{\partial q^{i}} \times \vec{x}_{\alpha}' \right)^{\mathsf{T}} & \cdots & \sum_{\alpha} m_{\alpha} \frac{\partial \vec{x}_{\alpha}'^{T}}{\partial q^{i}} \frac{\partial \vec{x}_{\alpha}'}{\partial q^{j}} & \cdots \\ \vdots & & \vdots & & \vdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \sum_{\alpha} m_{\alpha} \left( \frac{\partial \vec{x}_{\alpha}'}{\partial q^{i}} \times \vec{x}_{\alpha}' \right)^{\mathsf{T}} & \cdots & \sum_{\alpha} m_{\alpha} \frac{\partial \vec{R}}{\partial q^{i}} \frac{\partial \vec{R}}{\partial q^{i}} & \cdots \\ \vdots & & \vdots & & \vdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \sum_{\alpha} m_{\alpha} \left( \frac{\partial \vec{R}}{\partial q^{i}} \times \vec{x}_{\alpha}' \right)^{\mathsf{T}} & \cdots \\ \vdots & & \vdots & & \vdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \sum_{\alpha} m_{\alpha} \left( \frac{\partial \vec{R}}{\partial q^{i}} \times \vec{R}_{\alpha}' \right)^{\mathsf{T}} & \cdots \\ \vdots & & \vdots & \ddots & \vdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \sum_{\alpha} m_{\alpha} \left( \frac{\partial \vec{R}}{\partial q^{i}} \times \vec{R}_{\alpha}' \right)^{\mathsf{T}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \sum_{\alpha} m_{\alpha} \left( \frac{\partial \vec{R}}{\partial q^{i}} \times \vec{R}_{\alpha}' \right)^{\mathsf{T}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \sum_{\alpha} m_{\alpha} \left( \frac{\partial \vec{R}}{\partial q^{i}} \times \vec{R}_{\alpha}' \right)^{\mathsf{T}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \sum_{\alpha} m_{\alpha} \left( \frac{\partial \vec{R}}{\partial q^{i}} \times \vec{R}_{\alpha}' \right)^{\mathsf{T}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \sum_{\alpha} m_{\alpha} \left( \frac{\partial \vec{R}}{\partial q^{i}} \times \vec{R}_{\alpha}' \right)^{\mathsf{T}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \sum_{\alpha} m_{\alpha} \left( \frac{\partial \vec{R}}{\partial q^{i}} \times \vec{R}_{\alpha}' \right)^{\mathsf{T}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \cdots \\ m_{\text{tot}} \frac{\partial \vec{R}}{\partial q^{i}} & \cdots \\ m_{\text{tot}}$$

where the superindex T indicates matrix transposition,  $I^{(3)}$  denotes the  $3 \times 3$  identity matrix and  $\vec{x}'_{\alpha}$  is the position of atom  $\alpha$  in the reference frame fixed in the system (the "primed" reference frame).

Additionally, we denote the *total mass* of the system by  $m_{\text{tot}} := \sum_{\alpha} m_{\alpha}$ , the position of the *center of mass* of the system in the primed reference frame by  $\vec{R} := m_{\text{tot}}^{-1} \sum_{\alpha} m_{\alpha} \vec{x}'_{\alpha}$  and the *inertia tensor* of

the system, also in the primed reference frame, by

$$\mathcal{J} := \begin{pmatrix} \sum_{\alpha} m_{\alpha} ((x_{\alpha}^{\prime 2})^{2} + (x_{\alpha}^{\prime 3})^{2}) & -\sum_{\alpha} m_{\alpha} x_{\alpha}^{\prime 1} x_{\alpha}^{\prime 2} & -\sum_{\alpha} m_{\alpha} x_{\alpha}^{\prime 1} x_{\alpha}^{\prime 3} \\ -\sum_{\alpha} m_{\alpha} x_{\alpha}^{\prime 1} x_{\alpha}^{\prime 2} & \sum_{\alpha} m_{\alpha} ((x_{\alpha}^{\prime 1})^{2} + (x_{\alpha}^{\prime 3})^{2}) & -\sum_{\alpha} m_{\alpha} x_{\alpha}^{\prime 2} x_{\alpha}^{\prime 3} \\ -\sum_{\alpha} m_{\alpha} x_{\alpha}^{\prime 1} x_{\alpha}^{\prime 3} & -\sum_{\alpha} m_{\alpha} x_{\alpha}^{\prime 2} x_{\alpha}^{\prime 3} & \sum_{\alpha} m_{\alpha} ((x_{\alpha}^{\prime 1})^{2} + (x_{\alpha}^{\prime 2})^{2}) \end{pmatrix}.$$

$$(30)$$

The matrix  $v(\vec{R})$  is defined as:

$$v(\vec{R}) := \begin{pmatrix} 0 & -R^3 & R^2 \\ R^3 & 0 & -R^1 \\ -R^2 & R^1 & 0 \end{pmatrix}$$
 (31)

and  $\times$  denotes the usual vector cross product.

Then, since  $\sin^2 \theta$  may be integrated out in eq. (21), we can write, omitting additive constants, the kinetic entropy associated to the reduced mass-metric tensor g depending only on the soft internals  $q^i$ :

$$S_{\rm r}^{\rm k}(q^i) = \frac{R}{2} \ln[\det g_2(q^i)].$$
 (32)

Finally, one may note that, since  $\sin^2 \theta$  divides out in the second line of eq. (25) or, otherwise stated, eqs. (27) and (32) may be introduced in the first line, then the Fixman's potential is independent of the external coordinates as well.

#### Computational Methods

In the particular molecule treated in this work (the model dipeptide HCO-L-Ala-NH<sub>2</sub> in Fig. 1), the formulae in the preceding sections must be used with M=2, being the internal soft coordinates  $q^i \equiv (\phi, \psi)$  the typical Ramachandran angles<sup>75</sup> (see Table 3), the total number of coordinates N=48 and the number of hard internals L=40.

**Table 3.** SASMIC Internal Coordinates (ref. 59) in Z-Matrix Form of the Protected Dipeptide HCO-L-Ala-NH<sub>2</sub>.

| Atom name       | Bond length | Bond angle | Dihedral angle            |  |
|-----------------|-------------|------------|---------------------------|--|
| H <sub>1</sub>  |             |            |                           |  |
| $C_2$           | (2,1)       |            |                           |  |
| $N_3$           | (3,2)       | (3,2,1)    |                           |  |
| $O_4$           | (4,2)       | (4,2,1)    | (4,2,1,3)                 |  |
| $C_5$           | (5,3)       | (5,3,2)    | $\omega_0 := (5,3,2,1)$   |  |
| $H_6$           | (6,3)       | (6,3,2)    | (6,3,2,5)                 |  |
| C <sub>7</sub>  | (7,5)       | (7,5,3)    | $\phi := (7,5,3,2)$       |  |
| C <sub>8</sub>  | (8,5)       | (8,5,3)    | (8,5,3,7)                 |  |
| H <sub>9</sub>  | (9,5)       | (9,5,3)    | (9,5,3,7)                 |  |
| $H_{10}$        | (10,8)      | (10,8,5)   | $\chi := (10,8,5,3)$      |  |
| H <sub>11</sub> | (11,8)      | (11,8,5)   | (11,8,5,10)               |  |
| $H_{12}$        | (12,8)      | (12,8,5)   | (12,8,5,10)               |  |
| $N_{13}$        | (13,7)      | (13,7,5)   | $\psi := (13,7,5,3)$      |  |
| $O_{14}$        | (14,7)      | (14,7,5)   | (14,7,5,13)               |  |
| H <sub>15</sub> | (15,13)     | (15,13,7)  | $\omega_1 := (15,13,7,5)$ |  |
| H <sub>16</sub> | (16,13)     | (16,13,7)  | (16,13,7,15)              |  |

Principal dihedrals are indicated in bold face and their typical biochemical name is given.

![](_page_6_Picture_15.jpeg)

Figure 1. Atom numeration of the protected dipeptide HCO-L-Ala-NH<sub>2</sub>.

Regarding the side chain angle  $\chi$ , it has been argued elsewhere<sup>59</sup> that it is soft with the same right as the angles  $\phi$  and  $\psi$ , i.e., the barriers that hinder the rotation on this dihedral are comparable to the ones existing in the Ramachandran surface. However, the height of these barriers is sufficient ( $\sim$ 6-12 RT, see ref. 59) for the condition (ii) in Classical Stiff Model section to hold and, therefore, its inclusion in the set of hard coordinates is convenient due to its unimportant character (see discussion in the earlier section). Moreover, to describe the behaviour associated to  $\chi$  with a probability density different from a Gaussian distribution (i.e., its potential energy different from an harmonic oscillator), for example with the tools used in the field of circular statistics, <sup>76–78</sup> would severely complicate the derivation of the classical stiff model without adding any conceptual insight to the problem. In addition, although  $\chi$  is a periodic coordinate with threefold symmetry, the considerable height of the barriers between consecutive minima allows to make the quadratic assumption in eq. (3) at each equivalent valley and permits the approximation of the integral on  $\chi$  by three times a Gaussian integral. The multiplicative factor 3 simply adds a temperature- and conformation-independent reference to the configurational entropy  $S_s^c$  in eq. (12b).

The same considerations are applied to the dihedral angles,  $\omega_0$  and  $\omega_1$  (see Table 3), that describe the rotation around the peptide bond, and the quadratic approximation described above can also be used, since the heights of the rotation barriers around these degrees of freedom are even larger than the ones in the case of  $\chi$ .

The *ab initio* quantum mechanical calculations have been done with the package GAMESS<sup>79</sup> under Linux and in 3.20 GHz PIV machines. The coordinates used for the HCO-L-Ala-NH<sub>2</sub> dipeptide in the GAMESS input files and the ones used to generate them with automatic Perl scripts are the SASMIC coordinates introduced in ref. 59. They are presented in Table 3, indicating the name of the conventional dihedral angles (see also Fig. 1 for reference). To perform the energy optimizations, however, they have been converted to Delocalized Coordinates<sup>80</sup> in order to accelerate convergence.

First, we have calculated the typical potential energy surface (PES) in a regular  $12 \times 12$  grid of the bidimensional space spanned by the Ramachandran angles  $\phi$  and  $\psi$ , with both angles ranging from  $-165^{\circ}$  to  $165^{\circ}$  in steps of  $30^{\circ}$ . This has been done by running constrained energy optimizations at the MP2/6-31++G(d,p) level of the theory, freezing the two Ramachandran angles at each value of the grid, starting from geometries previously optimized at a lower

level of the theory and setting the gradient convergence criterium to OPTTOL =  $10^{-5}$  and the self-consistent Hartree-Fock convergence criterium to CONV=  $10^{-6}$ .

The results of these calculations (which took  $\sim$ 100 days of CPU time) are 144 conformations that define  $\Sigma$  and the values of  $V_{\Sigma}(\phi, \psi)$  at these points (the PES itself).

Then, at each optimized point of  $\Sigma$ , we have calculated the Hessian matrix in the coordinates of Table 3 removing the rows and columns corresponding to the soft angles  $\phi$  and  $\psi$ , the result being the matrix  $\mathcal{H}(\phi,\psi)$  in eq. (12b). This has been done, again, at the MP2/6-31++G(d,p) level of the theory, taking  $\sim$ 140 days of CPU time.

Equations (27) and (32) in sec. Factorization of the External Coordinates have been used to calculate the kinetic entropy terms associated to the determinants of the mass-metric tensors G and g, respectively. The quantities in eq. (27), being simply internal coordinates, have been directly extracted from the GAMESS output files via automated Perl scripts. On the other hand, in order to calculate the matrix  $g_2$  in eq. (29) that appears in the kinetic entropy of the classical rigid model, the Euclidean coordinates  $\vec{x}'_{\alpha}$  of the 16 atoms in the reference frame fixed in the system, as well as their derivatives with respect to  $q^i \equiv (\phi, \psi)$ , must be computed. For this, two additional 12 × 12 grids as the one described above have been computed; one of them displaced  $2^{\circ}$  in the positive  $\phi$ -direction and the other one displaced  $2^{\circ}$  in the positive  $\psi$ -direction. This has been done, again, at the MP2/6-31++G(d,p) level of the theory, starting from the optimized structures found in the computation of the PES described above and taking ~75 days of CPU time each grid. Using the values of the positions  $\vec{x}'_{\alpha}$  in these two new grids and also in the original one, the derivatives of these quantities with respect to the angles  $\phi$  and  $\psi$ , appearing in  $g_2$ , have been numerically obtained as finite differences.

The three calculations have been repeated for six special points in the Ramachandran space that correspond to important elements of secondary structure (see the next section), the total CPU time needed for computing all correcting terms at these points has been

**Table 4.** Maximum Variation<sup>a</sup>, Average<sup>b</sup> and Standard Deviation<sup>c</sup> in the  $12 \times 12$  Grid Defined in the Ramachandran Space of the Protected Dipeptide HCO-L-Ala-NH<sub>2</sub> for the Three Energy Surfaces,  $V_{\Sigma}$ ,  $F_s$ , and  $F_r$ , the Three Correcting Terms,  $-TS_s^k$ ,  $-TS_s^c$ , and  $-TS_r^k$  and the Fixman's Compensating Potential  $V_F$ .

|                             | MP2/6-31++G(d,p) |       |       | HF/6-31++G(d,p) |       |       |
|-----------------------------|------------------|-------|-------|-----------------|-------|-------|
|                             | Max.a            | Ave.b | Std.c | Max.a           | Ave.b | Std.c |
| $V_{\Sigma}$                | 21.64            | 6.76  | 3.88  | 23.62           | 6.92  | 4.35  |
| $\overline{F_{\mathrm{s}}}$ | 21.43            | 6.47  | 3.93  | 23.78           | 7.17  | 4.38  |
| $F_{\rm r}$                 | 21.09            | 6.46  | 3.82  | 23.09           | 6.76  | 4.31  |
| $-TS_s^k$                   | 0.24             | 0.09  | 0.05  | 0.23            | 0.09  | 0.04  |
| $-TS_s^c$                   | 1.67             | 0.98  | 0.32  | 1.34            | 0.63  | 0.30  |
| $-TS_{r}^{k}$               | 0.81             | 0.37  | 0.12  | 0.75            | 0.38  | 0.12  |
| $V_{\rm F}$                 | 1.68             | 0.89  | 0.30  | 1.35            | 0.55  | 0.27  |

The results at both MP2/6-31++G(d,p) and HF/6-31++G(d,p) levels of the theory are presented and all the functions have been referenced to zero in the grid. The units used are kcal/mol.

![](_page_7_Figure_10.jpeg)

**Figure 2.** PES of the model dipeptide HCO-L-Ala-NH<sub>2</sub>, computed at the MP2/6-31++G(d,p) level of the theory. The surface has been referenced to zero and smoothed with bicubic splines for visual convenience. The units in the *z*-axis are kcal/mol.

 $\sim$ 16 days. A total of  $\sim$ 406 days of CPU time has been needed to perform the whole study at the MP2/6-31++G(d,p) level of the theory.

Finally, we have repeated all the calculations at the HF/6-31++G(d,p) level of the theory in order to investigate if this less demanding method ( $\sim$ 10 days for the PES,  $\sim$ 8 days for the Hessians,  $\sim$ 10 days for each displaced grid,  $\sim$ 2 days for the special secondary structure points, being a total of  $\sim$ 40 days of CPU time) may be used instead of MP2 in further studies.

# Results

In Table 4, the maximum variation, the average and the standard deviation in the  $12 \times 12$  grid defined in the Ramachandran space of the protected dipeptide HCO-L-Ala-NH<sub>2</sub> are shown for the three energy surfaces,  $V_{\Sigma}$ ,  $F_{\rm s}$ , and  $F_{\rm r}$  (see eqs. (12) and (23)), for the three correcting terms,  $-TS_{\rm s}^k$ ,  $-TS_{\rm s}^c$ , and  $-TS_{\rm r}^k$  and for the Fixman's compensating potential  $V_{\rm F}$  (see eq. (25)). All the functions have been referenced to zero in the grid.

In Fig. 2, the PES  $V_{\Sigma}$ , at the MP2/6-31++G(d,p) level of the theory, is depicted with the reference set to zero for visual convenience. Neither the surfaces defined by  $F_s$  and  $F_r$  at the MP2/6-31++G(d,p) level of the theory nor the three energy surfaces  $V_{\Sigma}$ ,  $F_s$ , and  $F_r$  at HF/6-31++G(d,p) are shown graphically since they are visually very similar to the surface in Fig. 2.

In Fig. 3, the three correcting terms,  $-TS_s^k$ ,  $-TS_s^c$ , and  $-TS_r^k$  and the Fixman's compensating potential  $V_F$ , at the MP2/6-31++G(d,p) level of the theory, are depicted with the reference set to zero. The analogous surfaces at the HF/6-31++G(d,p) level of the theory are visually very similar to the ones in Fig. 3 and have been therefore omitted

From the results presented, one may conclude that, although the conformational dependence of the correcting terms  $-TS_{\rm s}^{\rm k}$ ,  $-TS_{\rm s}^{\rm c}$ , and  $-TS_{\rm r}^{\rm k}$  is more than an order of magnitude smaller than the conformational dependence of the PES  $V_{\Sigma}$  in the worst case, if

<sup>¶</sup> At the level of the theory used in the calculations, the minimum of  $V_{\Sigma}(\phi,\psi)$  in the grid is -416.0733418995 hartree.

![](_page_8_Figure_2.jpeg)

**Figure 3.** Ramachandran plots of the correcting terms appearing in eqs. (12) and (23), together with the Fixman's compensating potential defined in eq. (25), computed at the MP2/6-31++G(d,p) level of the theory in the model dipeptide HCO-L-Ala-NH<sub>2</sub>. The surfaces have been referenced to zero and smoothed with bicubic splines for visual convenience. The units in the *z*-axes are kcal/mol.

chemical accuracy (typically defined in the field of *ab initio* quantum chemistry as 1 kcal/mol<sup>81</sup>) is sought, they may be relevant. In fact, they are of the order of magnitude of the differences between the energy surfaces  $V_{\Sigma}$ ,  $F_s$ , and  $F_r$  calculated at MP2/6-31++G(d,p) and the ones calculated at HF/6-31++G(d,p).

For the same reasons, we may conclude that, if *ab initio* derived potentials are used to carry out Molecular Dynamics simulations of peptides, the Fixman's compensating potential  $V_{\rm F}$  should be included. Finally, regarding the relative importance of the different correcting terms  $-TS_{\rm s}^{\rm k}$ ,  $-TS_{\rm s}^{\rm c}$ , and  $-TS_{\rm r}^{\rm k}$ , the results in Table 4 suggest that the less important one is the kinetic entropy  $-TS_{\rm s}^{\rm k}$  of the stiff case (related to the determinant of the mass-metric tensor G) and that the most important one is the one related to the determinant of the Hessian matrix  $\mathcal H$  of the constraining part of the potential, i.e., the conformational entropy  $-TS_{\rm s}^{\rm k}$ . The first conclusion is in agreement with the approximations typically made in the literature, the second one, however, is not (see the Appendix).

Now, although the relative sizes of the conformational dependence of the different terms may be indicative of their importance, the degree of correlation among the surfaces is also relevant (see Table 5). Hence, in order to arrive to more precise conclusions, we re-examine here the results using a physically meaningful criterium to compare potential energy functions that has been introduced in ref. 82. The *distance*, denoted by  $d_{12}$ , between any two different potential energy functions,  $V_1$  and  $V_2$ , is an statistical quantity that, from a working set of conformations (in this case, the 144 points of the grid), measures the typical error that one makes in the *energy differences* if  $V_2$  is used instead of  $V_1$ , admitting a linear rescaling.

In Table 6, which contains the central results of this work, the distances between some of the energy surfaces that play a role in the problem are shown. We present the result in units of RT (at 300° K,

**Table 5.** Correlation Between the Different Correcting Terms Involved in the Study of the Constrained Equilibrium of the Protected Dipeptide HCO-L-Ala-NH<sub>2</sub>.

| $V_1{}^a$                                                         |                  | $V_2{}^{\mathrm{b}}$                                                                            | $r_{12}^{\mathrm{c}}$ |
|-------------------------------------------------------------------|------------------|-------------------------------------------------------------------------------------------------|-----------------------|
|                                                                   | MP2/6-           | 31++G(d,p)                                                                                      |                       |
| $V_{\Sigma}$                                                      | vs.              | $-TS_{\rm s}^{\rm c}$                                                                           | 0.1572                |
| $V_{\Sigma}$                                                      | vs.              | $-TS_{\mathrm{s}}^{\mathrm{c}} \ -TS_{\mathrm{s}}^{\mathrm{k}} \ -TS_{\mathrm{r}}^{\mathrm{k}}$ | -0.0008               |
| $V_{\Sigma}$                                                      | vs.              | $-TS_{\mathrm{r}}^{\mathrm{k}}$                                                                 | -0.3831               |
| $V_{\Sigma}$                                                      | vs.              | $V_{\rm F}$                                                                                     | 0.3334                |
|                                                                   | HF/6-3           | 31++G(d,p)                                                                                      |                       |
| $V_{\Sigma}$                                                      | vs.              | $-TS_{\rm s}^{\rm c}$                                                                           | 0.0682                |
| $V_{\Sigma}$                                                      | vs.              | $-TS_{\rm s}^{\rm k}$                                                                           | 0.0897                |
| $V_{\Sigma}$                                                      | vs.              | $-TS_{\mathrm{r}}^{\mathrm{k}}$                                                                 | -0.3544               |
| $V_{\Sigma}$                                                      | vs.              | $V_{\rm F}$                                                                                     | 0.2404                |
|                                                                   | MP2/6-31++G(d,t) | o) vs. HF/6-31++G(d,                                                                            | p)                    |
| $-TS_s^c$                                                         | vs.              | $-TS_s^c$                                                                                       | 0.9136                |
| $-TS_{\rm s}^{\rm c} \ -TS_{\rm s}^{\rm k} \ -TS_{\rm r}^{\rm k}$ | vs.              |                                                                                                 | 0.9808                |
| $-TS_{\rm r}^{\rm k}$                                             | vs.              | $ \begin{array}{l} -TS_{\rm s}^{\rm k} \\ -TS_{\rm r}^{\rm k} \end{array} $                     | 0.9316                |
| $V_{\rm F}$                                                       | vs.              | $V_{\rm F}$                                                                                     | 0.9217                |

<sup>&</sup>lt;sup>a</sup>Reference potential energy.

<sup>&</sup>lt;sup>b</sup>Approximated potential energy.

<sup>&</sup>lt;sup>c</sup>Pearson's correlation coefficient.

**Table 6.** Comparison of Different Energy Surfaces Involved in the Study of the Constrained Equilibrium of the Protected Dipeptide HCO-L-Ala-NH<sub>2</sub>. All quantities are dimensionless, except for  $d_{12}$ , which is given in units of the thermal energy RT at 300° K.

| Corr.a                | $V_1^{b}$                              | $V_2^{\rm c}$                     | $d_{12}^{\rm d}$ | $N_{\rm res}^{\ e}$ | $b_{12}^{\mathrm{f}}$ | $r_{12}^{g}$ |  |
|-----------------------|----------------------------------------|-----------------------------------|------------------|---------------------|-----------------------|--------------|--|
|                       |                                        | MP2/6-3                           | 1++G(d,p)        |                     |                       |              |  |
| $-TS_s^k - TS_s^c$    | $F_{\mathrm{s}}$                       | $V_{\Sigma}$                      | 0.74  RT         | 1.82                | 0.98                  | 0.9967       |  |
| $-TS_s^c$             | $F_{\mathrm{s}}$                       | $V_{\Sigma} - TS_{\rm s}^{\rm k}$ | 0.74~RT          | 1.83                | 0.98                  | 0.9967       |  |
| $-TS_{s}^{k}$         | $F_{\mathrm{s}}$                       | $V_{\Sigma} - TS_{\rm s}^{\rm c}$ | 0.11  RT         | 80.45               | 1.00                  | 0.9999       |  |
| $-TS_{\rm r}^{\rm k}$ | $F_{\rm r}$                            | $V_{\Sigma}$                      | 0.29  RT         | 11.62               | 1.01                  | 0.9995       |  |
| $V_{\rm F}$           | $F_{\rm s}$                            | $F_{\rm r}$                       | 0.67~RT          | 2.24                | 0.97                  | 0.9972       |  |
|                       |                                        | HF/6-31                           | ++G(d,p)         |                     |                       |              |  |
| $-TS_s^k - TS_s^c$    | $F_{\mathrm{s}}$                       | $V_{\Sigma}$                      | 0.73~RT          | 1.90                | 0.99                  | 0.9975       |  |
| $-TS_{s}^{c}$         | $F_{\mathrm{s}}$                       | $V_{\Sigma} - TS_{\rm s}^{\rm k}$ | 0.71  RT         | 2.00                | 0.99                  | 0.9976       |  |
| $-TS_{s}^{k}$         | $F_{\mathrm{s}}$                       | $V_{\Sigma} - TS_{\rm s}^{\rm c}$ | 0.10  RT         | 90.99               | 1.00                  | 0.9999       |  |
| $-TS_{\rm r}^{\rm k}$ | $F_{\rm r}$                            | $V_{\Sigma}$                      | 0.26RT           | 14.83               | 1.01                  | 0.9997       |  |
| $V_{\rm F}$           | $F_{\rm s}$                            | $F_{\rm r}$                       | 0.61~RT          | 2.69                | 0.98                  | 0.9982       |  |
|                       | MP2/6-31++G(d,p) vs. $HF/6-31++G(d,p)$ |                                   |                  |                     |                       |              |  |
|                       | $V_{\Sigma}$                           | $V_{\Sigma}$                      | 1.25 RT          | 0.64                | 1.12                  | 0.9925       |  |
|                       | $F_{\rm s}$                            | $F_{\rm s}$                       | 1.18 RT          | 0.72                | 1.11                  | 0.9934       |  |
|                       | $F_{\rm r}$                            | $F_{\rm r}$                       | 1.18 <i>RT</i>   | 0.72                | 1.12                  | 0.9932       |  |

<sup>&</sup>lt;sup>a</sup>Correcting term whose importance is measured in the corresponding row. <sup>b</sup>Reference potential energy  $V_1$  (the "correct" one, the one containing the correcting term).

where  $RT \simeq 0.6$  kcal/mol) because it has been argued in ref. 82 that, if the distance between two different approximations of the energy of the same system is less than RT, one may safely substitute one by the other without altering the relevant physical properties. Moreover, if one assumes that the effective energies compared will be used to construct a polypeptide potential and that it will be designed as simply the sum of mono-residue ones (making each term suitably depend on different pairs of Ramachandran angles), then, the number  $N_{\rm res}$  of residues up to which one may go keeping the distance between the two approximations of the N-residue potential below RT is (see eq. (23) in ref. 82):

$$N_{\rm res} = \left(\frac{RT}{d_{12}}\right)^2. \tag{33}$$

This number is also shown in Table 6, together with the slope  $b_{12}$  of the linear rescaling between  $V_1$  and  $V_2$  and the Pearson's correlation coefficient, <sup>83</sup> denoted by  $r_{12}$ .

The results at both MP2/6-31++G(d,p) and HF/6-31++G(d,p) levels of the theory are presented. The first three rows in each of the first two blocks are related to the classical stiff model, the next row to the classical rigid model, and the last one in each block to the comparison between the two models. The third block in the table is

associated to the comparison between the two different levels of the theory used.

The  $F_s$  vs.  $V_{\Sigma}$  row (in the first two blocks) assess the importance of the two correcting terms,  $-TS_s^k$  and  $-TS_s^c$ , in the stiff case. The result  $d_{12} = 0.74RT$  indicates that, for the alanine dipeptide,  $V_{\Sigma}$  may be used as an approximation of  $F_{\rm s}$  with caution if accurate results are sought. In fact, the low value of  $N_{\rm res} = 1.82 < 2$ shows that, if we wanted to describe a 2-residue peptide omitting the stiff correcting terms, we would typically make an error greater than the thermal noise in the energy differences. The next two rows investigate the effect of each one of the individual correcting terms. The conclusion that can be extracted from them (as the relative sizes in Table 4 already suggested) is that the conformational entropy associated to the determinant of the Hessian matrix  $\mathcal{H}$  is much more relevant than the correcting term  $-TS_s^k$ , related to the mass-metric tensor G, allowing to drop the latter up to  $\sim$ 80 residues (according to MP2/6-31++G(d,p) calculations). As has been already remarked, this second conclusion is in agreement with the approximations frequently done in the literature; however, it turns out that the importance of the Hessian-related term has been persistently underestimated (see the Appendix for a discussion).

The  $F_{\rm r}$  vs.  $V_{\Sigma}$  row, in turn, shows the data associated to the kinetic entropy term  $-TS_{\rm r}^{\rm k}$ , which is related to the determinant of the reduced mass-metric tensor g in the classical rigid model. From the results there  $(d_{12}=0.29RT$  and  $N_{\rm res}=11.62$  at the MP2/6-31++G(d,p) level), we can conclude that the only correction term in the rigid case is less important than the ones in the stiff case and that  $V_{\Sigma}$  may be used as an approximation of  $F_{\rm r}$  for oligopeptides of up to  $\sim$ 12 residues.

The last row in each of the first two blocks in Table 6 is related to the interesting question in molecular dynamics of whether or not one should include the Fixman's compensating potential  $V_F$  (see eq. (25)) in rigid simulations in order to obtain the stiff equilibrium distribution,  $\exp(-\beta F_s)$ , instead of the rigid one,  $\exp(-\beta F_r)$ . This question is equivalent to asking whether or not  $F_r$  is a good approximation of  $F_s$ . From the results in the table, we can conclude that the Fixman's potential is relevant for peptides of more than 2 residues and its omission may cause an error greater than the thermal noise in the energy differences.

The appreciable sizes of the different correcting terms, shown in Table 4, together with their low correlation with the PES  $V_{\Sigma}$ , presented in the first two blocks of Table 5, explain their considerable relevance discussed in the preceding paragraphs.

Moreover, from the comparison of the MP2/6-31++G(d,p) and the HF/6-31++G(d,p) blocks, one can tell that the study herein performed may well have been done at the lower level of the theory (if we had known) with a tenth of the computational effort (see the Methods section). This fact, explained by the high correlation, presented in the third block of Table 5, between the correcting terms calculated at the two levels, is very relevant for further studies on more complicated dipeptides or longer chains and it indicates that the differences in size between the different correcting terms at MP2/6-31++G(d,p) and HF/6-31++G(d,p), which are presented in Table 4, are mostly due to a harmless linear scaling effect similar to the well-known empirical scale factor frequently used in *ab initio* vibrational analysis.  $^{84-86}$  This view is supported by the data in the third block of Table 6, related to the comparison between the energy

<sup>&</sup>lt;sup>c</sup>Approximated potential energy  $V_2$  (i.e.  $V_1$  minus the correcting term in column a).

<sup>&</sup>lt;sup>d</sup>Statistical distance between  $V_1$  and  $V_2$  (see ref. 82).

<sup>&</sup>lt;sup>e</sup>Maximum number of residues in a polypeptide potential up to which the correcting term in column *a* may be omitted.

<sup>&</sup>lt;sup>f</sup>Slope of the linear rescaling between  $V_1$  and  $V_2$ .

<sup>&</sup>lt;sup>g</sup>Pearson's correlation coefficient.

**Table 7.** Ramachandran Angles (in Degrees) of Some Important Secondary Structure Elements in Polypeptides.

|                             | φ    | ψ   |
|-----------------------------|------|-----|
| α-helix                     | -57  |     |
| 3 <sub>10</sub> -helix      | -49  | -26 |
| $\pi$ -helix                | -57  | -70 |
| polyproline II              | -79  | 149 |
| parallel $\beta$ -sheet     | -119 | 113 |
| antiparallel $\beta$ -sheet | -139 | 135 |

Data taken from ref. 88.

surfaces calculated at MP2/6-31++G(d,p) and HF/6-31++G(d,p), where the slopes  $b_{12}$  are consistently larger than unity.

A last conclusion that may be extracted from the block labeled "MP2/6-31++G(d,p) vs. HF/6-31++G(d,p)" in Table 6 is that the typical error in the energy differences (given by the distances  $d_{12}$ ) produced when one reduces the level of the theory from MP2/6-31++G(d,p) to HF/6-31++G(d,p) is comparable (less than twice) to the error made if the most important correcting terms of the classical constrained models studied in this work are dropped. This is a useful hint for researchers interested in the conformational analysis of peptides with quantum chemistry methods  $^{60-64,70,87}$  and also to those whose aim is the design and parametrization of classical force fields from ab initio quantum mechanical calculations.  $^{68-70}$ 

Finally, in order to enrich and qualify the analysis, a new *working* set of conformations, different from the 144 points of the grid in the Ramachandran space, have been selected and the whole study has been repeated on them. These new conformations are six important

secondary structure elements which form repetitive patterns stabilized by hydrogen bonds in polypeptides. Their conventional names and the corresponding values of the  $\phi$  and  $\psi$  angles have been taken from ref. 88 and are shown in Table 7.

In Fig. 4, the relative energies of these conformations are shown for the three relevant potentials,  $V_{\Sigma}$ ,  $F_{\rm s}$ , and  $F_{\rm r}$ , at both MP2/6-31++G(d,p) and HF/6-31++G(d,p) levels of the theory. Since the antiparallel  $\beta$ -sheet is the structure with the minimum energy in all the cases, it has been set as the reference and the rest of energies in the figure should be regarded as relative to it.

The meaningful assessment, using the statistical distance described above, of the typical error made in the energy differences has been also performed on this new working set of conformations. The results are presented in Table 8.

The distances between the free energies,  $F_s$  and  $F_r$ , and their corresponding approximations obtained dropping the correcting entropies,  $-TS_s^k$ ,  $-TS_s^c$ , and  $-TS_r^k$ , or the Fixman's compensating potential  $V_F$ , in the first two blocks of the table, are consistently smaller than the ones found in the study of the grid defined in the whole Ramachandran space (cf. Table 6). And so are the distances between the three relevant potentials,  $V_{\Sigma}$ ,  $F_s$ , and  $F_r$ , calculated at the MP2/6-31++G(d,p) and HF/6-31++G(d,p) levels of the theory.

Although the distance  $d_{12}$  used is a statistical quantity and, therefore, one must be cautious when working with such a small set of conformations (of size six, in this case), the conclusion drawn from this second part of the study is that, if one is interested only in the "lower region" of the Ramachandran surface, where the typical secondary structure elements lie, then, one may safely neglect the conformational dependence of the different correcting terms appearing in the study of the constrained equilibrium of peptides. At least, up to oligopeptides (polyalanines) of  $\sim$ 10 residues in the worst case

![](_page_10_Figure_13.jpeg)

**Figure 4.** Relative energies of some important elements of secondary structure for the three potentials  $V_{\Sigma}$ ,  $F_s$ , and  $F_r$ , in the model dipeptide HCO-L-Ala-NH<sub>2</sub> and at both MP2/6-31++G(d,p) and HF/6-31++G(d,p) levels of the theory. The energy of the antiparallel β-sheet has been taken as reference. The units are kcal/mol.

**Table 8.** Comparison of Different Approximations to the Energies of Some Important Elements of Secondary Structure (see Table 7) in the Study of the Constrained Equilibrium of the Protected Dipeptide HCO-L-Ala-NH<sub>2</sub>.

| Corr.a                | $V_1^{b}$                              | $V_2^{\rm c}$                     | $d_{12}^{\mathrm{d}}$ | $N_{\rm res}^{}$ | $b_{12}^{\mathrm{f}}$ | $r_{12}{}^{g}$ |  |  |
|-----------------------|----------------------------------------|-----------------------------------|-----------------------|------------------|-----------------------|----------------|--|--|
|                       |                                        | MP2/6-3                           | 31++G(d,p)            | )                |                       |                |  |  |
| $-TS_s^k - TS_s^c$    | $F_{\mathrm{s}}$                       | $V_{\Sigma}$                      | 0.22RT                | 19.72            | 0.99                  | 0.9990         |  |  |
| $-TS_{\rm s}^{\rm c}$ | $F_{\mathrm{s}}$                       | $V_{\Sigma}-TS_{\rm s}^{\rm k}$   | 0.26 RT               | 14.07            | 0.98                  | 0.9985         |  |  |
| $-TS_{s}^{k}$         | $F_{\mathrm{s}}$                       | $V_{\Sigma} - TS_{\rm s}^{\rm c}$ | 0.06~RT               | 298.13           | 1.01                  | 0.9999         |  |  |
| $-TS_{\rm r}^{\rm k}$ | $F_{\rm r}$                            | $V_{\Sigma}$                      | 0.20RT                | 25.64            | 0.99                  | 0.9992         |  |  |
| $V_{\rm F}$           | $F_{\rm s}$                            | $F_{\rm r}$                       | 0.34~RT               | 8.73             | 0.99                  | 0.9977         |  |  |
|                       |                                        | HF/6-3                            | 1++G(d,p)             |                  |                       |                |  |  |
| $-TS_s^k - TS_s^c$    | $F_{\mathrm{s}}$                       | $V_{\Sigma}$                      | 0.14  RT              | 47.94            | 1.00                  | 0.9997         |  |  |
| $-TS_{\rm s}^{\rm c}$ | $F_{\mathrm{s}}$                       | $V_{\Sigma}-TS_{\rm s}^{\rm k}$   | 0.15  RT              | 46.12            | 1.00                  | 0.9997         |  |  |
| $-TS_{s}^{k}$         | $F_{\mathrm{s}}$                       | $V_{\Sigma} - TS_{\rm s}^{\rm c}$ | 0.05~RT               | 380.30           | 1.00                  | 0.9999         |  |  |
| $-TS_{\rm r}^{\rm k}$ | $F_{\rm r}$                            | $V_{\Sigma}$                      | 0.15  RT              | 41.85            | 0.99                  | 0.9997         |  |  |
| $V_{\rm F}$           | $F_{\rm s}$                            | $F_{\rm r}$                       | 0.18RT                | 30.12            | 1.01                  | 0.9996         |  |  |
|                       | MP2/6-31++G(d,p) vs. $HF/6-31++G(d,p)$ |                                   |                       |                  |                       |                |  |  |
|                       | $V_{\Sigma}$                           | $V_{\Sigma}$                      | 0.77  RT              | 1.68             | 1.28                  | 0.9929         |  |  |
|                       | $F_{\mathrm{s}}$                       | $F_{\mathrm{s}}$                  | 0.77  RT              | 1.69             | 1.26                  | 0.9928         |  |  |
|                       | $F_{\rm r}$                            | $F_{\rm r}$                       | 0.71 <i>RT</i>        | 1.96             | 1.28                  | 0.9939         |  |  |

See the footnotes of Table 6 for an explanation of the keys in the different columns

(the neglect of the Fixman's compensating potential  $V_F$  in the  $F_s$  vs.  $F_r$  comparison at MP2/6-31++G(d,p)).

This difference between the two working set of conformations may be explained looking at one of the ways of expressing the statistical distance used (see eq. (12a) in ref. 82):

$$d_{12} = \sqrt{2}\,\sigma_2 \big(1 - r_{12}^2\big)^{1/2},\tag{34}$$

where  $r_{12}$  is the Pearson's correlation coefficient between the potential energies denoted by  $V_1$  and  $V_2$  and  $\sigma_2$  is the standard deviation in the values of  $V_2$  on the relevant working set of conformations.

This last quantity,  $\sigma_2$ , is responsible for the differences between Tables 5 and 8, since the set of conformations comprised by the six secondary structure elements in Table 7 spans a smaller energy range than the whole PES in Fig. 2 (or  $F_s$ , or  $F_r$ , which have very similar variations). Accordingly, the dispersion in the energy values is smaller:  $\sigma_2 \simeq 2 \, \text{kcal/mol}$  in the case of the secondary structure elements and  $\sigma_2 \simeq 4 \, \text{kcal/mol}$  for the grid in the whole Ramachandran space (see Table 4). Since the correlation coefficient in both cases are of similar magnitude, the differences in  $\sigma_2$  produce a smaller distance  $d_{12}$  for the second set of conformations studied, i.e., a smaller typical error made in the energy differences when omitting the correcting terms derived from the consideration of constraints.

To end this section, we remark that although this "lower region" of the Ramachandran space contains the most relevant secondary structure elements (which are also the most commonly found in experimentally resolved native structures of proteins<sup>89–92</sup>) and may be the only region explored in the dynamical or thermodynamical study of small peptides, if the aim is the design of effective potentials for computer simulation of polypeptides, <sup>68–70</sup> then, some caution is recommended, since long-range interactions in the sequence may temporarily compensate local energy penalizations and the

higher regions of the energy surfaces studied could be important in transition states or in some relevant dynamical paths of the system.

In the following section, the many results discussed in the preceding paragraphs are summarized.

#### Conclusions

In this work, the theory of classical constrained equilibrium has been collected for the stiff and rigid models. The pertinent correcting terms, which may be regarded as effective entropies, as well as the Fixman's compensating potential, have been derived and theoretically discussed (see eqs. (12), (23), and (25), together with the formulae in the Methods section). In addition, the common approximation of considering that, for typical internals, the equilibrium values of the hard coordinates do not depend on the soft ones has also been discussed and related to the rest of simplifications. The treatment of both assumptions in the literature is thoroughly reviewed and discussed in the Appendix.

In the central part of the work (Results section), the relevance of the different correcting terms has been assessed in the case of the model dipeptide HCO-L-Ala-NH<sub>2</sub>, with quantum mechanical calculations including electron correlation. Also, the possibility of performing analogous studies at the less demanding Hartree-Fock level of the the theory has been investigated. The results found are summarized in the following points:

- In Monte Carlo simulations of the classical stiff model at room temperature, the effective entropy −TS<sub>s</sub><sup>k</sup>, associated to the determinant of the mass-metric tensor G, may be neglected for peptides of up to ~80 residues. Its maximum variation in the Ramachandran space is 0.24 kcal/mol.
- In Monte Carlo simulations of the classical stiff model at room temperature, the effective entropy  $-TS_s^c$ , associated to the determinant of the Hessian  $\mathcal{H}$  of the constraining part of the potential, should be included for peptides of more than 2 residues. Its maximum variation in the Ramachandran space is 1.67 kcal/mol.
- In Monte Carlo simulations of the classical rigid model at room temperature, the effective entropy  $-TS_{\rm r}^{\rm k}$ , associated to the determinant of the reduced mass-metric tensor g, may be neglected for peptides of up to  $\sim$ 12 residues. Its maximum variation in the Ramachandran space is 0.81 kcal/mol.
- In rigid Molecular Dynamics simulations intended to yield the stiff equilibrium distribution at room temperature, the Fixman's compensating potential V<sub>F</sub> should be included for peptides of more than 2 residues. Its maximum variation in the Ramachandran space is 1.68 kcal/mol.
- If the assumption that only the more stable region of the Ramachandran space, where the principal elements of secondary structure lie, is relevant, then, the importance of the correcting terms decreases and the limiting number of residues in a polypeptide potential up to which they may be omitted is approximately four times larger in each of the previous points.
- In both cases (i.e., either if the whole Ramachandran space is considered relevant or only the lower region), the errors made if the most important correcting terms are neglected are of the same order of magnitude as the errors due to a decrease in the level of theory from MP2/6-31++G(d,p) to HF/6-31++G(d,p).

• The whole study of the relevance of the different correcting terms (or future analogous investigations) may be performed at the HF/6-31++G(d,p) level of the theory, yielding very similar results to the ones obtained at MP2/6-31++G(d,p) and using a tenth of the computational effort.

To end this discussion, some qualifications should be made. On one hand, the conclusions above refer to the case in which a classical potential directly extracted from the quantum mechanical (Born-Oppenheimer) one is used; for the considerably simpler force fields typically used for macromolecular simulations, the study should be repeated and different results may be obtained. On the other hand, the investigation performed in this work has been done in one of the simplest dipeptides; both its isolated character and the relatively small size of its side chain play a role in the results obtained. Hence, for bulkier residues included in polypeptides, these conclusions should be approached with caution and much interesting work remains to be done.

## **Appendix**

Many approximations may be done to simplify the calculation of the different correcting terms introduced in the previous subsections. The most frequently found in the literature are the following three:

- (i) To neglect the conformational dependence of  $\det G$ .
- (ii) To neglect the conformational dependence of det  $\mathcal{H}$ .
- (iii) To assume that the hard coordinates are constant, i.e, that the f<sup>I</sup>(q<sup>I</sup>) in eq. (1) do not depend on the soft coordinates q<sup>I</sup>.

The conformational dependence of det g is customarily regarded as important, since it was shown to be non-negligible even for simple systems some decades ago<sup>13,21–23</sup> (normally in an indirect way, while studying the influence of the Fixman's compensating potential in eq. (25); see discussion below). With this same aim, Patriciu et al.<sup>20</sup> have very recently measured the conformational dependence of det g for a serial polymer with fixed bond lengths and bond angles (in the approximation (iii)), showing that it is non-negligible and suggesting that it may be so also for more general systems.

Note that, if approximations (i) and (ii) are assumed, then the Fixman's potential depends only on det g. In fact, whereas in the general case the Fixman's compensating potential can not be simplified beyond the expression in eq. (25), if one assumes approximation (iii), then the reduced mass-metric tensor g turns out to be the subblock of G with soft indices and, in this case, the quotient det G/ det g has been shown to be equal to 1/ det g0 by Fixman, where g1 denotes the sub-block of g2.

$$h^{IJ}(q^{\mu}) := \sum_{\sigma=1}^{N} \frac{\partial q^{I}}{\partial x^{\sigma}} \frac{1}{m_{\sigma}} \frac{\partial q^{J}}{\partial x^{\sigma}}.$$
 (A1)

This result has been extensively used in the literature,  $^{21-23,38,42,57}$  since each of the internal coordinates  $q^a$  typically used in macromolecular simulations only involves a small number of atoms, thus rendering the matrix h above sparse and allowing for efficient algorithms to be used in order to find its determinant.

Now, although det g is customarily regarded as important, the conformational variations of det G are almost unanimously neglected (approximation (i)) in the literature  $^{15,55}$  and may only be said to be indirectly included in h by the authors that use the expression above.  $^{20-23,38,57}$  This is mainly due to the fact, reported by  $G\bar{o}$  and Scheraga  $^{15}$  and, before, by Volkenstein,  $^{74}$  that det G in a serial polymer may be expressed as in eq. (26), being independent of the dihedral angles (which are customarily taken as the soft coordinates). If one also assumes approximation (iii), which, as will be discussed later, is very common, then det G is a constant for every conformation of the molecule.

Probably because of computational considerations, but also sometimes to the use of a formulation of the stiff case based on delta functions,  $^{51}$  the conformational dependence of det  ${\cal H}$  is almost unanimously neglected (approximation (ii)) in the literature.  $^{15,16,19,20,41,55,65,66}$  Only a few authors include this term in different stages of the reasoning,  $^{13-15,18,19,26,39}$  most of them only to argue later that it is negligible.

Although for some simple *ad hoc* designed potentials that lack long-range terms, <sup>21,22,57</sup> the aforementioned simplifying assumptions and the ones that will be discussed in the following paragraphs may be exactly fulfilled, in the case of the potential energies used in force fields for macromolecular simulation, <sup>6,27–37</sup> they are not. The typical energy function in this case has the form

$$\begin{split} V_{\mathrm{ff}}(q^{a}) &:= \frac{1}{2} \sum_{\alpha=1}^{N_{r}} K_{r_{\alpha}} (r_{\alpha} - r_{\alpha}^{0})^{2} + \frac{1}{2} \sum_{\alpha=1}^{N_{\theta}} K_{\theta_{\alpha}} (\theta_{\alpha} - \theta_{\alpha}^{0})^{2} \\ &+ V_{\mathrm{ff}}^{\mathrm{tors}}(\phi_{\alpha}) + V_{\mathrm{ff}}^{\mathrm{long-range}}(q^{a}), \quad (A2) \end{split}$$

where  $r_{\alpha}$  are bond lengths,  $\theta_{\alpha}$  are bond angles,  $\phi_{\alpha}$  are dihedral angles, and, for the sake of simplicity, no harmonic terms have been assumed for out-of-plane angles or for hard dihedrals (such as the peptide bond  $\omega$ ).  $N_r$  is the number of bond lengths,  $N_{\theta}$  the number of bond angles, and the quantities  $K_{r_{\alpha}}$ ,  $K_{\theta_{\alpha}}$ ,  $r_{\alpha}^{0}$ , and  $\theta_{\alpha}^{0}$  are constants. The term denoted by  $V_{\rm ff}^{\rm tors}(\phi_{\alpha})$  is a commonly included torsional potential that depends only on the dihedral angles  $\phi_{\alpha}$  and  $V_{\rm ff}^{\rm long-range}(q^{a})$  normally comprises long-range interactions such as Coulomb or van der Waals; hence, it depends on the atomic positions  $\vec{x}_{\alpha}'$  which, in turn, depend on all the internal coordinates  $q^{a}$ .

One of the reasons given for neglecting det  $\mathcal{H}$ , when classical force fields are used with potential energy functions such as the one in eq. (A2), is that the harmonic constraining terms dominate over the rest of interactions and, since the constants appearing on these terms (the  $K_{r_{\alpha}}$ ,  $K_{\theta_{\alpha}}$  in eq. (A2)) are independent of the conformation by construction, so is det  $\mathcal{H}$ . <sup>15,19,26</sup> Here, we analyze a more realistic quantum-mechanical potential and these considerations are not applicable; however, they also should be checked in the case of classical force fields, since, for a potential energy such as the one in eq. (A2)), the quantities  $K_{r_{\alpha}}$  and  $K_{\theta_{\alpha}}$  are finite and the long-range terms will also affect the Hessian at each point of the constrained hypersurface  $\Sigma$ , rendering its determinant conformation-dependent.

For the same reason, even in classical force fields, the equilibrium values of the hard coordinates are not the constant quantities  $r_{\alpha}^{0}$  and  $\theta_{\alpha}^{0}$  in eq. (A2) but some functions  $f^{I}(q^{i})$  of the soft coordinates (see eq. (1)). This fact, recognized by some authors, <sup>15,25,45,93</sup> provokes that, if one chooses to assume approximation (iii) and

the constants  $r_{\alpha}^0$  and  $\theta_{\alpha}^0$  appearing in eq. (A2) are designated as the equilibrium values, the potential energy in  $\Sigma$  may be heavily distorted, the cause being simply that the long-range interactions between atoms separated by three covalent bonds are not fully relaxed. This effect is probably larger if bond angles, and not only bond lengths, are also constrained, which may partially explain the different dynamical behavior found in ref. 6 when comparing these types of constraints in molecular dynamics simulations. In quantum mechanical calculations of small dipeptides, on the other hand, the fact that the bond lengths and bond angles depend on the Ramachandran angles  $(\phi, \psi)$  has been pointed out by Schäffer et al. Herefore, approximation (iii), which is very common in the literature,  $^{6,13,15,16,18-20,26,38-42,51,52,55,65,66,95,96}$  should be critically analyzed in each particular case.

Apart from the typical internal coordinates  $q^a$  used until now, in terms of which the constrained hypersurface  $\Sigma$  is described by the relations  $q^I = f^I(q^i)$  in eq. (1), with  $I = M + 7, \ldots, N$ , one may define a different set  $Q^a$  such that, on  $\Sigma$ , the corresponding hard coordinates are arbitrary constants  $Q^I = C^I$  (the external coordinates  $q^A$  and  $Q^A$  are irrelevant for this part of the discussion). To do this, for example, let

$$Q^{i} := q^{i}$$
  $i = 7, ..., M + 6$  and  $Q^{I} := q^{I} - f^{I}(q^{i}) + C^{I}$   $I = M + 7, ..., N$ . (A3)

Well then, while the relation between bond lengths, bond angles, and dihedral angles (the typical  $q^{a59}$ ) and the Euclidean coordinates is straightforward and simple, the expression of the transformation functions  $Q^a(x^{\mu})$  needs the knowledge of the  $f^I$ , which must be calculated numerically in most real cases. This drastically reduce the practical use of the  $Q^a$ ; however, it is also true that they are conceptually appealing, since they have a property that closely match our intuition about what the soft and hard coordinates should be (namely, that the hard coordinates  $Q^I$  are constant on the relevant hypersurface  $\Sigma$ ), and this is why we term them exactly separable hard and soft coordinates. Now, we must also point out that, although the real internal coordinates  $q^a$  do not have this property, they are usually close to it. The customary labeling of soft and hard coordinates in the literature is based on this circumstance. Somehow, the dihedral angles are the "softest" of the internal coordinates, i.e., the ones that "vary the most" when the system visits different regions of the hypersurface  $\Sigma$ , and this is why we term the real q<sup>a</sup> approximately separable hard and soft coordinates, considering approximation (iii) as a useful reference case.

To sum up, the three simplifying assumptions (i), (ii), and (iii) in the beginning of this section should be regarded as approximations in the case of classical force fields, as well as in the case of the more realistic quantum-mechanical potential investigated in this work, and they should be critically assessed in the systems of interest. Here, while studying the model dipeptide HCO-L-Ala-NH<sub>2</sub>, no simplifying assumptions of this type have been made.

# Acknowledgments

We would like to thank F. Falceto and V. Laliena for illuminating discussions and also the reviewers of the manuscript for much useful

suggestions. The numerical calculations have been performed at the BIFI computing facilities. We thank I. Campos, for the invaluable CPU time and the efficiency at solving the problems encountered.

### References

- Alonso, J. L.; Chass, G. A.; Csizmadia, I. G.; Echenique, P.; Tarancón, A. Meeting on Fundamental Physics 'Alberto Galindo'. In Álvarez-Estrada, R. F.; Dobado, A.; Fernández, L. A.; Martín-Delgado, M. A.; Muñox Sudupe, A. Eds.; Aula Documental: Madrid, 2004.
- 2. Dobson, C. M. Nature 2003, 426, 884.
- 3. Dill, K. A. Prot Sci 1999, 8, 1166.
- 4. He, S.; Scheraga, H. A. J Chem Phys 1998, 108, 287.
- Abagyan, R. A.; Totrov, M. M.; Kuznetsov, D. A. J Comp Chem 1994, 15, 488.
- 6. Van Gunsteren, W. F.; Karplus, M. Macromolecules 1982, 15, 1528.
- Levinthal, C. Mossbauer Spectroscopy in Biological Systems; In DeBrunner, J. T. P.; Munck, E. Eds.; University of Illinois Press: Illinois, 1969.
- Chun, H. M.; Padilla, C. E.; Chin, D. N.; Watanabe, M.; Karlov, V. I.; Alper, H. E.; Soosaar, K.; Blair, K. B.; Becker, O. M.; Caves, L. S. D.; Nagle, R.; Haney, D. N.; Farmer, B. L. J Comp Chem 2000, 21, 159.
- 9. Reich, S. Physica D 2000, 118, 210.
- 10. Reich, S. J Comput Phys 1999, 151, 49.
- Schlick, T.; Barth, E.; Mandziuk, M. Annu Rev Biophys Biomol Struct 1997, 26, 181.
- 12. Van Kampen, N. G.; Lodder, J. J. Am J Phys 1984, 52, 419.
- 13. Rallison, J. M. J Fluid Mech 1979, 93, 251.
- 14. Helfand, E. J Chem Phys 1979, 71, 5000.
- 15. Gō, N.; Scheraga, H. A. Macromolecules 1976, 9, 535.
- 16. Fixman, M. Proc Natl Acad Sci USA 1974, 71, 3050.
- 17. Gō, N.; Scheraga, H. A. J Chem Phys 1969, 51, 4751.
- 18. Morse, D. C. Adv Chem Phys 2004, 128, 65.
- 19. Den Otter, W. K.; Briels, W. J. Mol Phys 2000, 98, 773.
- Patriciu, A.; Chirikjian, G. S.; Pappu, R. V. J Chem Phys 2004, 121, 12708.
- 21. Perchak, D.; Skolnick, J.; Yaris, R. Macromolecules 1985, 18, 519.
- 22. Pear, M. R.; Weiner, J. H. J Chem Phys 1979, 71, 212.
- 23. Chandler, D.; Berne, B. J. J Chem Phys 1979, 71, 5386.
- 24. Gottlieb, M.; Bird, R. B. J Chem Phys 1976, 65, 2467.
- 25. Zhou, J.; Reich, S.; Brooks, B. R. J Chem Phys 2000, 111, 7919.
- Berendsen, H. J. C.; Van Gunsteren, W. F. The Physics of Superionic Conductors and Electrode Materials, In Perram, J. W., Ed.; NATO ASI Series B92. Plenum Press: New York, 1983; pp. 221–240.
- MacKerell, A. D., Jr.; Brooks, B.; Brooks, C. L. III; Nilsson, L.; Roux, B.; Won, Y.; Karplus, M. The Encyclopedia of Computational Chemistry, In Schleyer, P. V. R.; Schreiner, P. R.; Allinger, N. L.; Clark, T.; Gasteiger, J.; Kollman, P.; Schaefer III, H. F. Eds.; Wiley: Chichester, 1998; pp. 217–277.
- Brooks, B. R.; Bruccoleri, R. E.; Olafson, B. D.; States, D. J.; Swaminathan, S.; Karplus, M. J Comp Chem 1983, 4, 187.
- Cornell, W. D.; Cieplak, P.; Bayly, C. I.; Gould, I. R. Jr.; Merz, K. M.; Ferguson, D. M.; Spellmeyer, D. C.; Fox, T.; Caldwell, J. W.; Kollman, P. A. J Am Chem Soc 1995, 117, 5179.
- Pearlman, D. A.; Case, D. A.; Caldwell, J. W.; Ross, W. R.; Cheatham,
   T. E. III; DeBolt, S.; Ferguson, D.; Seibel, G.; Kollman, P. Comp Phys
   Commun 1995, 91, 1.
- 31. Jorgensen, W. L.; Tirado-Rives, J. J Am Chem Soc 1988, 110, 1657.
- Jorgensen, W. L.; Maxwell, D. S.; Tirado-Rives, J. J Am Chem Soc 1996, 118, 11225.
- 33. Halgren, T. A. J Comp Chem 1996, 17, 490.
- 34. Halgren, T. A. J Comp Chem 1996, 17, 520.

- 35. Halgren, T. A. J Comp Chem 1996, 17, 553.
- 36. Halgren, T. A. J Comp Chem 1996, 17, 587.
- 37. Halgren, T. A. J Comp Chem 1996, 17, 616.
- 38. Pasquali, M.; Morse, D. C. J Chem Phys 2002, 116, 1834.
- 39. Den Otter, W. K.; Briels, W. J. J Chem Phys 1998, 109, 4139.
- 40. Ciccotti, G.; Ryckaert, J. P. Comput Phys Rep 1986, 4, 345.
- 41. Berendsen, H. J. C.; Van Gunsteren, W. F. In Molecular Liquids— Dynamics and Interactions; Barnes, A. J.; Orville-Thomas, W. J.; Yarwood, J. Eds.; Reidel: Dordrecht, Holland, 1984; pp. 475–500.
- 42. Fixman, M. J Chem Phys 1978, 69, 1527.
- 43. Álvarez-Estrada, R. F.; Calvo, G. F. Mol Phys 2002, 100, 2957.
- 44. Álvarez-Estrada, R. F. Macromol Theory Simul 2000, 9, 83.
- 45. Hess, B.; Saint-Martin, H.; Berendsen, H. J. C. J Chem Phys 2002, 116, 9602.
- 46. Born, M.; Oppenheimer, J. R. Ann Phys Leipzig 1927, 84, 457.
- 47. Wilson, E. B. Jr.; Decius, J. C.; Cross, P. C. Molecular Vibrations: The Theory of Infrared and Raman Vibrational Spectra; Dover: New York, 1980.
- 48. Barth, E.; Kuczera, K.; Leimkuhler, B.; Skeel, R. D. J Comp Chem 1995, 16, 1192.
- 49. Andersen, H. C. J Comput Phys 1983, 52, 24.
- 50. Ryckaert, J. P.; Ciccotti, G.; Berendsen, H. J. C. J Comput Phys 1977, 23, 327.
- 51. Schlitter, J.; Klän, M. Mol Phys 2003, 101, 3439.
- 52. Van Gunsteren, W. F. In Computer Simulations of Biomolecular Systems, Van Gunsteren, W. F.; Weiner, P. K. Eds.; Escom Science: Netherlands, 1989.
- 53. Dinner, A. R. J Comp Chem 2000, 21, 1132.
- 54. Schofield, J.; Ratner, M. A. J Chem Phys 1998, 109, 9177.
- 55. Pertsin, A. J.; Hahn, J.; Grossmann, H. P. J Comp Chem 1994, 15, 1121.
- 56. Knapp, E. W.; Irgens-Defregger, A. J Fluid Mech 1993, 14, 19.
- 57. Almarza, N. G.; Enciso, E.; Alonso, J.; Bermejo, F. J.; Álvarez, M. Mol Phys 1990, 70, 485.
- 58. Echenique, P.; Calvo, I. J Comp Chem (accepted).
- 59. Echenique, P.; Alonso, J. L. J Comp Chem 2006, 27, 1076.
- 60. Láng, A.; Csizmadia, I. G.; Perczel, A. Proteins: Struct Funct Bioinf 2005, 58, 571.
- 61. Perczel, A.; Farkas, O.; Jakli, I.; Topol, I. A.; Csizmadia, I. G. J Comp Chem 2003, 24, 1026.
- 62. Vargas, R.; Garza, J.; Hay, B. P.; Dixon, D. A. J Phys Chem A 2002, 106, 3213.
- 63. Yu, C.-H.; Norman, M. A.; Schäfer, L.; Ramek, M.; Peeters, A.; van Alsenoy, C. J Mol Struct 2001, 567–568, 361.
- 64. Császár, A. G.; Perczel, A. Prog Biophys Mol Biol 1999, 71, 243.
- 65. Allen, M. P.; Tildesley, D. J. Computer Simulation of Liquids; Clarendon: Oxford, 2005.
- 66. Frenkel, D.; Smit, B. Understanding Molecular Simulations: From Algorithms to Applications, 2nd ed.; Academic Press: Orlando FL, 2002.
- 67. Karplus, M.; Kushick, J. N. Macromolecules 1981, 14, 325.

- 68. MacKerell, A. R. Jr.; Feig, M.; Brooks, C. L. III. J Comp Chem 2004, 25, 1400.
- 69. Bordner, A. J.; Cavasotto, C. N.; Abagyan, R. A. J Phys Chem B 2003, 107, 9601.
- 70. Beachy, M.; Chasman, D.; Murphy, R.; Halgren, T.; Friesner, R. J Am Chem Soc 1997, 119, 5908.
- 71. Arnold, V. I. Mathematical Methods of Classical Mechanics (Graduate Texts in Mathematics); 2nd ed.; Springer: New York, 1989.
- 72. Viskolcz, B.; Fejer, S. N.; Csizmadia, I. G. J Phys Chem A. (in press). ASAP article http://pubs.acs.org/cgi-bin/abstract.cgi/jpcafh/asap/abs/ jp058219t.html.
- 73. Andricioaei, I.; Karplus, M. J Chem Phys 2001, 115, 6289.
- 74. Volkenstein, M. V. Configurational Statistical of Polymeric Chains; Interscience: New York, 1959.
- 75. Ramachandran, G. N.; Ramakrishnan, C. J Mol Biol 1963, 7, 95.
- 76. Hnizdo, V.; Fedorowicz, A.; Singh, H.; Demchuk, E. J Comp Chem 2003, 24, 1172.
- 77. Demchuk, E.; Singh, H. Mol Phys 2001, 99, 627.
- 78. Mardia, K. V.; Jupp, P. E. Directional Statistics; Wiley: Chichester, 2000.
- 79. Schmidt, M. W.; Baldridge, K. K.; Boatz, J. A.; Elbert, S. T.; Gordon, M. S.; Jensen, H. J.; Koseki, S.; Matsunaga, N.; Nguyen, K. A.; Su, S.; Windus, T. L.; Dupuis, M.; Montgomery, J. A. J Comp Chem 1993, 14, 1347.
- 80. Baker, J.; Kessi, A.; Delley, B. J Chem Phys 1996, 105, 192.
- 81. Daykov, I. P.; Arias, T. A.; Engeness, T. D. Phys Rev Lett 2003, 90, 216402.
- 82. Alonso, J. L.; Echenique, P. J Comp Chem 2006, 27, 238.
- 83. Dobson, J. D. Applied Multivariate Data Analysis, Vol. I; Springer-Verlag: New York, 1991.
- 84. Levine, I. N. Quantum Chemistry; 5th ed.; Prentice Hall: Upper Saddle River, 1999.
- 85. Halls, M. D.; Velkovski, J.; Schlegel, H. B. Theo Chem Acc 2001, 105, 413.
- 86. Scott, A. P.; Radom, L. J Phys Chem 1996, 100, 16502.
- 87. Elstner, M.; Jalkanen, K. J.; Knapp-Mohammady, M.; Suhai, S. Chem Phys 2001, 263, 203.
- 88. Lesk, A. M. Introduction to Protein Architecture; Oxford University Press: Oxford, 2001.
- 89. Chakrabarti, P.; Pal, D. Prog Biophys Mol Biol 2001, 76, 1.
- 90. Berman, H. M.; Westbrook, J.; Feng, Z.; Gilliland, G.; Bhat, T. N.; Weissig, H.; Shindyalov, I. N.; Bourne, P. E. Nucleic Acids Research 2000, 28, 235.
- 91. Gunasekaran, K.; Ramakrishnan, C.; Balaram, P. J Mol Biol 1996, 264, 191.
- 92. Creighton, T. E. Proteins: Structures and Molecular Properties; 2nd ed.; W. H. Freeman: New York, 1992.
- 93. Chen, J.; Im, W.; Brooks, C. L. III. J Comp Chem 2005, 26, 1565.
- 94. Schäfer, L.; Ming, C. J Mol Struct 1995, 333, 201.
- 95. Mazars, M. J Phys A: Math Gen 1998, 31, 1949.
- 96. Mazars, M. Phys Rev E 1996, 53, 6297.