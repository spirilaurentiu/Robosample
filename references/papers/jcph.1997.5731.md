# Compensating Mass Matrix Potential for Constrained Molecular Dynamics

Abhinandan Jain<sup>1</sup>

Jet Propulsion Laboratory/California Institute of Technology, 4800 Oak Grove Drive, Pasadena, California 91109

Received January 1997

Rigid internal constraints are used in molecular models to speed up molecular dynamics (MD) simulations. It is well recognized that statistical averages from such constrained MD simulations differ by a metric tensor-dependent term from similar averages computed using conventional unconstrained MD simulations. Fixman proposed augmenting the standard potential with a compensating term which depends on the metric tensor to nullify the effects of this bias term. However, in the absence of tractable algorithms to compute this compensating tensor potential and its gradient its use has been impractical. This paper derives a new algorithm for computing the compensating potential, as well as its gradient for tree topology molecular systems. The algorithm is quite straightforward and is an extension of the spatial operators based  $O(\mathcal{N})$  algorithm that has been recently proposed for constrained dynamics. Indeed, the compensating potential is closely related and computed from the articulated body inertia quantities available from this  $O(\mathcal{N})$ algorithm. © 1997 Academic Press

#### 1. INTRODUCTION

Rigid internal constraints are often used in molecular models to eliminate high frequency modes and to enable the use of large numerical integration time steps necessary for speeding up molecular dynamics (MD) simulations [1, 2]. There has been considerable debate regarding the relationship between the statistical averages obtained from such constrained MD simulations and those obtained from conventional unconstrained Cartesian model MD simulations [3–7]. In particular, Fixman [3] pointed out that ensemble averages obtained from MD simulation using the constrained and unconstrained models will differ due to the presence of a metric tensor dependent term in the partition function for constrained molecular models. and this has been verified in simulations. Fixman proposed the augmentation of the standard potential by a compensating metric tensor potential in constrained MD simulations to compensate for the effects of the bias term. Several researchers [8-11] have verified the efficacy of this method for simple molecular systems. The prohibitive complexity of computing the metric tensor potential and its gradient has been a major hurdle on the use of these compensating potentials more generally in constrained MD simulations.

This paper analyzes the structure of the compensating potential and its gradient and develops substantially simpler expressions for them for tree topology molecular models. These expressions are used to derive computational algorithms for use in constrained MD simulations. The algorithms are straightforward extensions of the recently proposed spatial operators based  $O(\mathcal{N})$  algorithm for constrained MD simulations [2]. Indeed, the compensating potential and its gradient are readily computable from the articulated body inertia quantities available from this  $O(\mathcal{N})$  algorithm.

### 2. ENSEMBLE AVERAGES IN CONSTRAINED DYNAMICS SIMULATIONS

The partition function  $\mathcal{K}(\mathcal{T})$  for an n degree of freedom *unconstrained* Cartesian molecular model is given by the expression

$$\mathcal{Z}(\mathcal{T}) = \frac{1}{h^n} \int_{-\infty}^{\infty} \cdots \int_{-\infty}^{\infty} \int_{-\alpha_n}^{\gamma_n} \cdots \int_{-\alpha_1}^{\gamma_1} \exp\left[-\left(\frac{1}{2}p^* \mathcal{M}^{-1}p + \mathcal{V}\right) \middle/ k\mathcal{T}\right]$$

$$\times dp_1 \cdots dp_n dq_1 \cdots dq_n$$
(2.1)

where  $\mathscr{T}$  denotes the temperature, and  $q_i$  and  $p_i$  are the configuration and momentum coordinates,  $\mathscr{V}$  is the standard potential energy,  $\mathscr{M} \in \mathbb{R}^{n \times n}$  denotes the system **mass matrix** (or *metric tensor*). The system kinetic energy is given by  $\frac{1}{2}p^*\mathscr{M}^{-1}p$ , and the limits of integration  $\alpha_i$  and  $\gamma_i$  are determined by the geometry of the problem [12]. For conventional unconstrained Cartesian dynamics, the mass matrix  $\mathscr{M}$  is constant (and diagonal) and does not depend upon the system configuration. As a consequence, the ensemble average of a function f(q) is given by

$$\langle f(q) \rangle = \frac{2\pi k \mathcal{T}^{n/2}}{h^n} \int_{-\alpha_n}^{\gamma_n} \cdots \int_{-\alpha_1}^{\gamma_1} f(q) \exp[-\mathcal{V}/k\mathcal{T}] \times dq_1 \cdots dq^n / \mathcal{X}(\mathcal{T}),$$
(2.2)

<sup>&</sup>lt;sup>1</sup> E-mail: Abhinandan.Jain@jpl.nasa.gov.

where the momentum variables have been integrated over and eliminated.

In *constrained* dynamics models, the molecular system is modeled as a collection of rigid **clusters** coupled together by articulable **hinges** [2]. For such constrained dynamics models, the mass matrix  $\mathcal{M}(\boldsymbol{\theta}) \in \mathbb{R}^{\mathcal{N} \times \mathcal{N}}$  is a function of the internal configuration coordinates vector  $\boldsymbol{\theta}$  and  $\mathcal{N}$  the number of degrees of freedom for the constrained model. The kinetic energy K.E. is now given by the expression

K.E. = 
$$\frac{1}{2}\beta^* \mathcal{M}(\boldsymbol{\theta})\beta = \frac{1}{2}p^* \mathcal{M}^{-1}(\boldsymbol{\theta})p,$$
 (2.3)

where  $\beta$  denotes the internal velocity coordinates vector, and where the conjugate momenta vector p is given by the expression

$$p = \mathcal{M}(\boldsymbol{\theta})\beta. \tag{2.4}$$

The partition function  $\mathcal{K}(\mathcal{T})$  for the system is given by the expression

$$\mathcal{Z}'(\mathcal{T}) = \frac{1}{h^{\mathcal{N}}} \int_{-\infty}^{\infty} \cdots \int_{-\infty}^{\infty} \int_{-\alpha_{\mathcal{N}}}^{\gamma_{\mathcal{N}}} \cdots \int_{-\alpha_{1}}^{\gamma_{1}} \exp\left[-\left(\frac{1}{2}p^{*}\mathcal{M}^{-1}(\boldsymbol{\theta})p + \mathcal{V}\right) \middle/ k\mathcal{T}\right]$$

$$\times dp'_{1} \cdots dp'_{\mathcal{N}} d\boldsymbol{\theta}_{1} \cdots d\boldsymbol{\theta}_{\mathcal{N}}.$$

$$(2.5)$$

Using a diagonalizing transformation on the momentum coordinates of the form

$$p' = \mathcal{M}^{-1/2}(\boldsymbol{\theta})p \tag{2.6}$$

we can integrate Eq. (2.5) over the new momentum coordinates to get the following expression for the partition function

$$\mathcal{Z}(\mathcal{T}) = \frac{2\pi k \mathcal{T}^{\mathcal{X}/2}}{h^{\mathcal{X}}} \int_{-\alpha_{\mathcal{X}}}^{\gamma_{\mathcal{X}}} \cdots \int_{-\alpha_{1}}^{\gamma_{1}} \det\{\mathcal{M}^{1/2}(\boldsymbol{\theta})\}$$

$$\exp[-\mathcal{V}/k\mathcal{T}] \times d\boldsymbol{\theta}_{1} \cdots d\boldsymbol{\theta}_{K}.$$
(2.7)

This implies that the ensemble average of a function  $f(\theta)$  over the configuration space defined by the internal coordinates is given by

$$\langle f(\boldsymbol{\theta}) \rangle' = \frac{2\pi k \mathcal{T}^{N/2}}{h^{N}} \int_{-\alpha_{N}}^{\gamma_{N}} \cdots \int_{-\alpha_{1}}^{\gamma_{1}} \det\{\mathcal{M}^{1/2}(\boldsymbol{\theta})\} f(\boldsymbol{\theta})$$

$$\exp[-\mathcal{V}/k\mathcal{T}] \times d\boldsymbol{\theta}_{1} \cdots d\boldsymbol{\theta}_{N}/\mathcal{X}'.$$
(2.8)

In contrast with the unconstrained Cartesian dynamics expression in Eq. (2.2), the integrand in the partition function integral in Eq. (2.5) involves the determinant of the mass matrix. This additional term introduces a bias in statistical averages computed using constrained dynamics models versus those obtained using Cartesian models [3].

A method for bridging this gap in statistical estimates was proposed by Fixman [3]. He suggested replacing the potential function  $\mathcal{V}(\boldsymbol{\theta})$  by the modified potential function  $\mathcal{V}'(\boldsymbol{\theta})$ 

$$\mathcal{V}'(\boldsymbol{\theta}) \stackrel{\Delta}{=} \mathcal{V}(\boldsymbol{\theta}) + \mathcal{V}_c(\boldsymbol{\theta}), \text{ where } \mathcal{V}_c(\boldsymbol{\theta}) \stackrel{\Delta}{=} \frac{1}{2} \ln \det \{ \mathcal{M}(\boldsymbol{\theta}) \}$$
(2.9)

in constrained MD simulations. It is easy to see that replacing  $\mathscr{V}(\theta)$  in Eq. (2.7) and Eq. (2.8) by  $\mathscr{V}'(\theta)$  eliminates the metric tensor from the partition functions and the expression for the ensemble average. The extra potential term  $\mathscr{V}_c(\theta)$  is referred to as the **compensating mass matrix potential** or the *metric-tensor potential* [7], since it effectively compensates for the mass matrix determinant term in Eq. (2.8).

With the use of  $\mathcal{V}'(\boldsymbol{\theta})$ , the remaining differences between statistical averages computed using constrained and unconstrained MD simulations are due to the coarser sampling of the conformational space during the averaging process in constrained MD simulations. The appropriate use of constraints is important for ensuring that the loss in fidelity is within acceptable limits for the simulation experiment at hand. The use of the compensating mass matrix potential helps ensure that at least the systematic bias term from the metric tensor does not contribute to the averaging errors.

The implication of using  $\mathcal{V}'(\boldsymbol{\theta})$  is that now its gradient must be used for computing the forces during constrained MD simulations. The overall hinge torque vector T' is defined as the gradient of  $\mathcal{V}'(\boldsymbol{\theta})$  and is given by

$$T' = \nabla_{\theta} \mathcal{V}'(\boldsymbol{\theta}) = T + T_c, \text{ where}$$

$$T \stackrel{\Delta}{=} \nabla_{\theta} \mathcal{V}(\boldsymbol{\theta}), T_c \stackrel{\Delta}{=} \nabla_{\theta} \mathcal{V}_c(\boldsymbol{\theta}).$$
(2.10)

 $T_c$  represents the compensating hinge torque arising from the compensating mass-matrix potential  $\mathcal{V}_c(\boldsymbol{\theta})$  and must be used in addition to the standard torque term T during constrained MD simulations. Its kth element,  $T_c(k)$ , for the kth hinge is given by

$$T_c(k) = \frac{\partial \mathcal{V}_c(\boldsymbol{\theta})}{\partial \boldsymbol{\theta}(k)} = \frac{1}{2} \frac{\partial \ln \det{\{\mathcal{M}(\boldsymbol{\theta})\}}}{\partial \boldsymbol{\theta}(k)}.$$
 (2.11)

The need for the use of the compensating potential has been verified via computer simulations for several small systems [7, 8, 10, 11, 13, 14]. There is general consensus that the compensating potential is less important for rigid constraints on the bond stretching degrees of freedom while they are a significant factor for rigid constraints involving bond angles. It was found [7, 8] that the use of the compensating potential in constrained MD simulations for *n*-butane other systems effectively bridged the gap in the number of dihedral transitions using constrained and unconstrained MD simulations. Go and Scheraga [5, 6] examined the relative merits of Cartesian dynamics models and constrained dynamics models and concluded that the former was the more correct model for molecular dynamics simulations. However, Refs. [8, 9] argued that the constrained MD simulations produce similar statistical estimates as the Cartesian case if the compensating potential is included.

Despite the accepted importance of using the compensating potential during constrained dynamics simulations, it is rarely used in practice. The primary factor has been the lack of a tractable method for computing the compensating torque,  $T_c$ , for all but simple molecular systems [13]. A clever method for computing  $\mathcal{V}_c$  was proposed by Fixman [4]. However, this method is also limited to moderatively sized molecular systems and does not include a procedure for computing the compensating torque  $T_c$ .

In the following sections we take a closer look at the expression for  $T_c$  in Eq. (2.11) and derive simpler expressions for it using techniques from the spatial operator alge**bra** [2]. These expressions lead to a simple method for computing this compensating torque for arbitrary treetopology molecular systems. This method is an extension of the spatial operators based algorithm for constrained dynamics simulations [2]. This algorithm provides a highly efficient  $O(\mathcal{N})$  recursive method for solving the equations of motion for constrained systems without using iterative procedures such as in the SHAKE algorithm [15]. The equations of motion for the system are solved exactly and the complexity of the algorithm is  $O(\mathcal{N})$ , i.e., the computational cost grows only linearly with the number of unconstrained degrees of freedom in the system. The algorithm—which has been implemented as the NEIMO (Newton–Euler inverse mass operator method) software package [2, 16, 17] is based upon closed-form spatial operator expressions for the factorization and inversion of the mass matrix.

Expression for the Compensating Torque  $T_c(k)$ 

In general, if g(X) is a scalar function of a matrix  $X \in \mathbb{R}^{m \times n}$ , then its derivative with respect to a variable y is given by

$$\frac{\partial g(X)}{\partial y} = \sum_{i=1}^{m} \sum_{j=1}^{n} \frac{\partial g(X)}{\partial X(i,j)} \frac{\partial X(i,j)}{\partial (y)}$$

$$= \operatorname{Trace} \left\{ \frac{\partial g}{\partial X} * \frac{\partial X}{\partial y} \right\}, \tag{2.12}$$

where  $\partial g(X)/\partial X(i, j)$  and  $\partial X(i, j)/\partial y$  are  $m \times n$  matrices whose elements are defined as

$$\frac{\partial g}{\partial X}(i,j) \triangleq \frac{\partial g(X)}{\partial X(i,j)}, \qquad \frac{\partial X}{\partial y}(i,j) \triangleq \frac{\partial X(i,j)}{\partial y}.$$

For the scalar function  $g(X) \stackrel{\Delta}{=} \ln \det\{X\}$ , it is a well established fact [18] that

$$\frac{\partial g(X)}{\partial X} = \frac{\partial \ln \det\{X\}}{\partial X} = \{X^*\}^{-1}.$$
 (2.13)

Using Eq. (2.13) and Eq. (2.12) in Eq. (2.11) leads to the following expression for  $T_c(k)$ :

$$T_{c}(k) = \frac{1}{2}\operatorname{Trace}\left\{\mathcal{M}^{-1}(\boldsymbol{\theta})\frac{\partial\mathcal{M}(\boldsymbol{\theta})}{\partial\boldsymbol{\theta}(k)}\right\}$$
$$= \frac{1}{2}\operatorname{Trace}\left\{\mathcal{M}^{-1}(\boldsymbol{\theta})\mathcal{M}_{\boldsymbol{\theta}(k)}(\boldsymbol{\theta})\right\}.$$
 (2.14)

We have used the notational shorthand  $\mathcal{M}_{\theta(k)}(\theta)$  in place of  $\partial \mathcal{M}(\theta)/\partial \theta(k)$  in Eq. (2.14). In the following sections we use spatial operator expressions for  $\mathcal{M}(\theta)$  to further simplify Eq. (2.14).

#### 3. SPATIAL OPERATOR FORM OF MASS MATRIX

Now we briefly review the derivation of the spatial operator equations of motion and refer the reader to Ref. [2] for more detailed discussion on the notation and the concepts. For notational simplicity we limit out initial discussion to an n-cluster system with serial chain structure and single degree of freedom rotational hinges between clusters. The number of degrees of freedom for this system is  $\mathcal{N} = n + 5$ . Later we discuss the steps involved in extending the derivations and algorithms to general tree-topology systems with multiple degrees of freedom hinges.

As shown in Ref. [2], the Newton–Euler recursive equations of motion for the whole system have the form:

for 
$$k = n \cdots 1$$
  
 $V(k) = \phi^*(k+1, k)V(k+1) + H^*(k)\theta(k)$   
 $\alpha(k) = \phi^*(k+1, k)\alpha(k+1) + H^*(k)\theta(k) + a(k)$   
\nend loop  

$$f(0) = 0$$
  
for  $k = 1 \cdots n$   
 $f(k) = \phi(k, k-1)f(k-1) + M(k)\alpha(k) + b(k) + \hat{f}_c(k)$   
 $T(k) = H(k)f(k)$   
\nend loop

 $V(n + 1) = 0, \alpha(n + 1) = 0$ 

The introduction of **spatial operators** allows the expression of the equations of motion in the following more concise form:

$$V = \phi^* \mathbf{H}^* \dot{\boldsymbol{\theta}}$$

$$\alpha = \phi^* (\mathbf{H}^* \ddot{\boldsymbol{\theta}} + a)$$

$$f = \phi(\mathbf{M}\alpha + b + \hat{f}_c)$$

$$= \phi \mathbf{M} \phi^* \mathbf{H}^* \ddot{\boldsymbol{\theta}} + \phi(\mathbf{M} \phi^* a + b + \hat{f}_c)$$

$$T = \mathbf{H} f = \mathbf{H} \phi \mathbf{M} \phi^* \mathbf{H}^* \ddot{\boldsymbol{\theta}} + \mathbf{H} \phi(\mathbf{M} \phi^* a + b + \hat{f}_c).$$
(3.2)

In particular, the equations of motion have the form

$$T = \mathcal{M}(\boldsymbol{\theta})\ddot{\boldsymbol{\theta}} + \mathcal{C}(\boldsymbol{\theta}, \dot{\boldsymbol{\theta}}), \tag{3.3}$$

where

$$\mathcal{M}(\boldsymbol{\theta}) \stackrel{\Delta}{=} \boldsymbol{H} \phi \boldsymbol{M} \phi^* \boldsymbol{H}^* \in \mathbb{R}^{\mathcal{N} \times \mathcal{N}}$$
 (3.4a)

$$\mathscr{C}(\boldsymbol{\theta}, \dot{\boldsymbol{\theta}}) \stackrel{\Delta}{=} \boldsymbol{H} \phi(\boldsymbol{M} \phi^* a + b + \hat{f}_c) \in \mathbb{R}^{\mathcal{N}}. \tag{3.4b}$$

Here,  $\mathcal{M}(\boldsymbol{\theta})$  is the **mass matrix** of the serial chain and  $\mathcal{C}(\boldsymbol{\theta}, \boldsymbol{\theta})$  is the vector of Coriolis, centrifugal, gyroscopic, and Cartesian forces. Note that  $\mathcal{M}$  and  $\mathcal{C}$  are nonlinear functions of  $\boldsymbol{\theta}$  and  $\boldsymbol{\dot{\theta}}$ . The factorization in Eq. (3.4) of the mass matrix  $\mathcal{M}$  is referred to as the **Newton–Euler operator factorization** [19] because it is equivalent to the recursive Newton–Euler inverse dynamics algorithm in Eq. (3.1).

The solution of the equations of motion in Eq. (3.3) for the accelerations vector  $\ddot{\boldsymbol{\theta}}$  is used by the numerical integrator to propagate the state of the system during molecular dynamics simulations. However, Eq. (3.3) represents only a conceptual statement of the dynamics problem since  $\mathcal{M}$  and  $\mathcal{C}$  are not explicitly available. The conventional approach for computing the accelerations  $\ddot{\boldsymbol{\theta}}$  consists of first computing both  $\mathcal{M}$  and  $\mathcal{C}$  and solving the linear matrix equation for the vector  $\ddot{\boldsymbol{\theta}}$ . In general,  $\mathcal{M}$  is fully

populated and, as a result, the computational cost of solving the equations of motion using this method grows cubically with the number of degrees of freedom in the system; i.e., this method is of  $O(\mathcal{N}^3)$  computational complexity. The computational advantage of larger integration timesteps using constrained dynamics [20] can be lost due to the large computational costs for large molecules from the  $O(\mathcal{N}^3)$  dependency.

In the next section we review an alternative recursive algorithm for computing the vector of generalized accelerations  $\boldsymbol{\theta}$  without having to explicitly compute the mass matrix [2]. The complexity of this method is only  $O(\mathcal{N})$ ; i.e., its computational cost grows only linearly with the number of degrees of freedom in the model.

#### Innovations Factorization of the Mass Matrix

The  $O(\mathcal{N})$  spatial algebra algorithm for solving the equations of motion in Eq. (3.3) described in [2] depends on the following key results that give explicit analytical operator expressions for the square factorization and inversion of the mass matrix.

LEMMA 3.1.

$$\mathcal{M} = [I + \mathbf{H}\phi \mathbf{K}]\mathbf{D}[I + \mathbf{H}\phi \mathbf{K}]^*$$
 (3.5a)

$$[I + \mathbf{H}\phi\mathbf{K}]^{-1} = [I - \mathbf{H}\psi\mathbf{K}]$$
 (3.5b)

$$\mathcal{M}^{-1} = [I - \mathbf{H}\psi \mathbf{K}] * \mathbf{D}^{-1}[I - \mathbf{H}\psi \mathbf{K}]. \tag{3.5c}$$

*Proof.* See Ref. [2].

The new square factorization described in Eq. (3.5a) is also referred to as the Innovations Operator Factorization of the mass matrix and is an alternative to the factorization in Eq. (3.4). The factor  $[I + H\phi K] \in \mathbb{R}^{\mathcal{N} \times \mathcal{N}}$  is square, block lower triangular, and nonsingular, while D is a block diagonal matrix. This factorization provides a closed form operator expression for the block LDL\* decomposition of M. The following lemma gives the closed form operator expression for the inverse of the factor  $[I + H\phi K]$ . Once again, the factor  $[I - H\psi K]$  is square, block lower triangular, and nonsingular and so Eq. (3.5c) provides a closed form expression for the block LDL\* decomposition of  $\mathcal{M}^{-1}$ . The spatial operators  $\phi$ , K, and D embedded in these factorizations are based on spatially recursive filtering and smoothing algorithms [19, 21, 22]. The following Riccati equation for the articulated body inertia P is a key element of these filtering and smoothing algorithms.

ALGORITHM 3.1. The articulated body inertia quantities  $P(\cdot)$ ,  $D(\cdot)$ ,  $G(\cdot)$ ,  $K(\cdot)$ ,  $\tau(\cdot)$ ,  $\overline{\tau}(\cdot)$ ,  $P^+(\cdot)$ , and  $\psi(\cdot, \cdot)$  are computed by the following recursive procedure:

$$P^{+}(0) = 0$$
for  $k = 1 \cdots \mathcal{N}$ 

$$P(k) = \phi(k, k - 1)P^{+}(k - 1)\phi^{*}(k, k - 1) + M(k)$$

$$D(k) = H(k)P(k)H^{*}(k)$$

$$G(k) = P(k)H^{*}(k)D^{-1}(k)$$

$$K(k + 1, k) = \phi(k + 1, k)G(k)$$

$$\tau(k) = G(k)H(k)$$

$$\bar{\tau}(k) = I - \tau(k)$$

$$P^{+}(k) = \bar{\tau}(k)P(k)$$

$$\psi(k + 1, k) = \phi(k + 1, k)\bar{\tau}(k)$$
The series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the series of the

#### end loop

Algorithm 3.1 is by now the classical [19, 23] Riccati equation of Kalman filtering. Its solution P(k) is the articulated body inertia [19, 24] of the part of the system outboard of hinge k. The operator **P** is a block-diagonal  $6n \times$ 6n matrix with its kth diagonal element being  $P(k) \in \mathbb{R}^{6\times 6}$ . Define also

$$D = HPH^* \in \mathbb{R}^{N \times N}$$

$$G = PH^*D^{-1} \in \mathbb{R}^{6n \times N}$$

$$K = \mathcal{E}_{\phi}G \in \mathbb{R}^{6n \times n}$$

$$\overline{\tau} = I - GH \in \mathbb{R}^{6n \times 6n}$$

$$\mathcal{E}_{\psi} = \mathcal{E}_{\phi}\overline{\tau} \in \mathbb{R}^{6n \times 6n}$$

$$\psi = (I - \mathcal{E}_{\psi})^{-1} \in \mathbb{R}^{6n \times 6n}.$$
(3.7)

The operators D, G, and  $\overline{\tau}$  are all block diagonal. The operators K and  $\mathcal{E}_{\psi}$  are not block-diagonal, but their only nonzero block elements are K(k, k-1)'s and  $\psi(k, k-1)$ 1)'s, respectively, along the first subdiagonal. The following lemma describes the operator expression for the generalized accelerations  $\ddot{\theta}$  in terms of the hinge forces T and Cartesian spatial forces  $\hat{f}_c$ .

LEMMA 3.2.

$$\ddot{\boldsymbol{\theta}} = [I - \boldsymbol{H}\psi\boldsymbol{K}]^*\boldsymbol{D}^{-1}[T - \boldsymbol{H}\psi\{\boldsymbol{K}T + \boldsymbol{P}a + b + \hat{f}_c\}] - \boldsymbol{K}^*\psi^*a$$
(3.8)

*Proof.* See Ref. [2]. ■

 $O(\mathcal{N})$  Algorithm for Solving the Equations of Motion

The recursive implementation of Eq. (3.8) leads to the following  $O(\mathcal{N})$  computational algorithm for the accelerations.  $\ddot{\boldsymbol{\theta}}$ :

$$z^{+}(0) = 0$$

$$\mathbf{for} \, k = \mathbf{1} \cdots \mathbf{n}$$

$$z(k) = \phi(k, k - 1)z^{+}(k - 1)$$

$$+ \mathbf{P}(k)a(k) + b(k) + \hat{f}_{c}(k)$$

$$\varepsilon(k) = T(k) - \mathbf{H}(k)z(k)$$

$$\mathbf{\nu}(k = \mathbf{D}^{-1}(k)\varepsilon(k)$$

$$z^{+}(k) = z(k) + \mathbf{G}(k)\varepsilon(k)$$
\nend loop

end loop

$$\alpha(n+1) = 0$$

$$\mathbf{for} \, \mathbf{k} = \mathbf{n} \cdots \mathbf{1}$$

$$\alpha^{+}(k) = \phi^{*}(k+1, k)\alpha(k+1)$$

$$\ddot{\boldsymbol{\theta}}(k) = \boldsymbol{\nu}(k) - \boldsymbol{G}^{*}(k)\alpha^{+}(k)$$

$$\alpha(k) = \alpha^{+}(k) + \boldsymbol{H}^{*}(k)\ddot{\boldsymbol{\theta}}(k) + a(k)$$

$$(3.10)$$

end loop

This algorithm does not require either the explicit computation of the mass matrix  $\mathcal{M}$ , nor the numerical solution of the matrix equation (3.3). The steps in the above algorithm can be summarized as follows:

- 1. The first step is a recursion from the base to the tip to compute the orientation, location, and spatial velocities, V(k), and the Coriolis and gyroscopic terms a(k) and b(k)for each of the clusters using the first base-to-tip recursion in Eq. (3.1).
- 2. Next follows a recursion from the tip towards the base as defined by Eq. (3.6) to compute the P(k)'s etc.
- 3. The recursion in Eq. (3.9) from the tip to the base is used next to compute the residual forces z(k) etc. This recursion can be combined with the tip to base recursion in the previous step to obtain a single tip to base recursion sequence.
- 4. Finally, the base to tip recursion described by Eq. (3.10) computes the  $\ddot{\theta}(k)$  accelerations for all the clusters.

The computational cost of this algorithm depends only linearly on the number of clusters. The structure of this algorithm closely resembles those found in Kalman filtering and smoothing theory [21, 25].

#### 4. COMPENSATING MASS MATRIX TORQUE $T_c(i)$

Using Eq. (2.14) as a starting point, we now develop an expression for the compensating torque  $T_c(i)$  that is simple to compute. While we have already seen an expression for the mass matrix inverse,  $\mathcal{M}^{-1}(\boldsymbol{\theta})$ , we need an expression for the derivative of the mass matrix with respect to hinge coordinates.

Spatial Operator Expression for  $\mathcal{M}_{\theta}$ .

**LEMMA 4.1.** 

$$\mathscr{M}_{\boldsymbol{\theta}_{i}} = \boldsymbol{H} \boldsymbol{\phi} [\mathbb{H}_{\delta}^{i} \boldsymbol{\phi} \boldsymbol{M} - \boldsymbol{M} \boldsymbol{\phi}^{*} \mathbb{H}_{\delta}^{i}] \boldsymbol{\phi}^{*} \boldsymbol{H}^{*}. \tag{4.11}$$

Proof. See Ref. [26].

The matrix  $\mathbb{H}_{\delta}^{i}$  a new quantity in this result.  $\mathbb{H}_{\delta}^{i}$  is the  $6n \times 6n$  matrix whose elements are all zero, except for a single  $6 \times 6$  block  $\mathbb{H}(i)$  at the ith location on the diagonal. The index i corresponds to the joint-angle  $\theta_{i}$  with respect to which the sensitivity  $\mathcal{M}_{\theta_{i}}$  is being taken. The nonzero block-diagonal element  $\mathbb{H}(i) \in \mathbb{R}^{6 \times 6}$  is obtained as follows from the hinge rotational axis unit vector h(i):

$$\mathbb{H}(i) = \begin{pmatrix} \tilde{h}(i) & 0\\ 0 & \tilde{h}(i) \end{pmatrix}. \tag{4.12}$$

The notation  $\tilde{v}$  above denotes the 3  $\times$  3 cross-product tensor matrix associated with a 3-vector

$$v = \begin{pmatrix} x \\ y \\ z \end{pmatrix}$$

and is defined as

$$\tilde{v} \stackrel{\Delta}{=} \begin{pmatrix} 0 & -z & y \\ z & 0 & -x \\ -y & x & 0 \end{pmatrix}.$$

The formula in Eq. (4.11) is closed-form, in the sense that it contains an explicit analytical expression for the mass matrix sensitivity in terms of the operators  $\phi$ , M, and H appearing in the mass matrix itself. That the formula is closed-form is of extreme importance, because it implies that the mass matrix derivatives can be easily computed using operations and spatially recursive algorithms similar to those used to compute the mass matrix itself.

Spatial Operator Expressions for  $T_c(i)$ 

The following lemma uses Lemma 4.11 to develop a new expression for  $T_c(i)$ .

LEMMA 4.2.

$$T_{c}(i) = \operatorname{Trace} \{ \mathbf{P} \mathbf{\Omega} \mathbb{H}_{\delta}^{i} \},$$

$$where \mathbf{\Omega} \stackrel{\Delta}{=} \psi^{*} \mathbf{H}^{*} \mathbf{D}^{-1} \mathbf{H} \psi \in \mathbb{R}^{6n \times 6n}.$$

$$(4.13)$$

*Proof.* Two spatial operator identities that we need for the proof are given by the following equations. Their derivation can be found in the appendix of Ref. [27]:

$$[I - H\psi K]H\phi = H\psi \tag{4.14a}$$

$$\phi M\Omega = (\phi - \psi) + P\Omega. \tag{4.14b}$$

The derivation of Eq. (4.13) goes as follows:

$$T_{c}(i) \stackrel{2.11}{=} \frac{1}{2}\operatorname{Trace}\{\mathcal{M}^{-1}(\boldsymbol{\theta})\mathcal{M}_{\boldsymbol{\theta}(k)}(\boldsymbol{\theta})\}$$

$$\stackrel{3.5c}{=} \frac{4.11}{2}\operatorname{Trace}\{[I - \boldsymbol{H}\psi\boldsymbol{K}]*\boldsymbol{D}^{-1}[I - \boldsymbol{H}\psi\boldsymbol{K}]$$

$$\boldsymbol{H}\phi[\mathbb{H}^{i}_{\delta}\phi\boldsymbol{M} - \boldsymbol{M}\phi*\mathbb{H}^{i}_{\delta}]\phi*\boldsymbol{H}^{*}\}$$

$$\stackrel{4.14a}{=} \operatorname{Trace}\{[I - \boldsymbol{H}\psi\boldsymbol{K}]*\boldsymbol{D}^{-1}\boldsymbol{H}\psi\mathbb{H}^{i}_{\delta}\phi\boldsymbol{M}\phi*\boldsymbol{H}^{*}\}$$

$$= \operatorname{Trace}\{\phi*\boldsymbol{H}*[I - \boldsymbol{H}\psi\boldsymbol{K}]*\boldsymbol{D}^{-1}\boldsymbol{H}\psi\mathbb{H}^{i}_{\delta}\phi\boldsymbol{M}\}$$

$$\stackrel{4.14a}{=} \operatorname{Trace}\{\psi*\boldsymbol{H}*\boldsymbol{D}^{-1}\boldsymbol{H}\psi\mathbb{H}^{i}_{\delta}\phi\boldsymbol{M}\}$$

$$\stackrel{4.13}{=} \operatorname{Trace}\{\boldsymbol{\Omega}\mathbb{H}^{i}_{\delta}\phi\boldsymbol{M}\}$$

$$= \operatorname{Trace}\{\boldsymbol{\phi}\boldsymbol{M}\boldsymbol{\Omega}\mathbb{H}^{i}_{\delta}\}$$

$$\stackrel{4.14b}{=} \operatorname{Trace}\{(\phi - \psi + \boldsymbol{P}\boldsymbol{\Omega})\mathbb{H}^{i}_{\delta}\}$$

$$= \operatorname{Trace}\{\boldsymbol{P}\boldsymbol{\Omega}\mathbb{H}^{i}_{\delta}\}.$$

In the above steps, we have used the fact that  $\operatorname{Trace}\{AB\} = \operatorname{Trace}\{BA\}$ . Also, the last step used the fact that  $\operatorname{Trace}\{(\phi - \psi)\mathbb{H}^i_{\delta}\} = 0$ . This is true because  $(\phi - \psi)$  is strictly lower triangular and  $\mathbb{H}^i_{\delta}$  is block diagonal.

It has been shown in Ref. [19] that  $\Omega$  can be decomposed as

$$\mathbf{\Omega} = \mathbf{Y} + \tilde{\psi} * \mathbf{Y} + \mathbf{Y} \tilde{\psi}, \tag{4.15}$$

where  $\tilde{\psi} \stackrel{\Delta}{=} \psi - I$ , and the diagonal elements  $\mathbf{Y}(k, k) \in \mathbb{R}^{6 \times 6}$  of the block diagonal matrix  $\mathbf{Y} \in \mathbb{R}^{6n \times 6n}$  are defined via the recursion:

$$\mathbf{Y}(n+1) = 0$$

$$\mathbf{for} \ k = n \cdots 1$$

$$\mathbf{Y}(k) = \psi^*(k+1, k)\mathbf{Y}(k+1)\psi(k+1, k) \quad (4.16)$$

$$+ \mathbf{H}^*(k)\mathbf{D}^{-1}(k)\mathbf{H}(k)$$

end loop

This allows us to further simplify the expression for  $T_c(i)$  as described in the following lemma.

LEMMA 4.3.

$$T_c(i) = \operatorname{Trace}\{\mathbf{P}(i)\mathbf{Y}(i)\mathbb{H}(i)\}. \tag{4.17}$$

*Proof.* Using Eq. (4.15) in Eq. (4.13) it follows that

$$T_c(i) = \operatorname{Trace}\{P\Omega \mathbb{H}_{\delta}^i\} = \operatorname{Trace}\{P(Y + \tilde{\psi}^*Y + Y\tilde{\psi})\mathbb{H}_{\delta}^i\}$$
  
=  $\operatorname{Trace}\{PY \mathbb{H}_{\delta}^i\} = \operatorname{Trace}\{P(i)Y(i)\mathbb{H}(i)\}.$ 

Here we used the fact that Trace  $\{PY\tilde{\psi}\mathbb{H}^i_{\delta}\}=0$  because  $\tilde{\psi}$  is strictly lower triangular.

A final simplification step for  $T_c(i)$  is described in the following lemma.

Lemma 4.4. Let the partitioned form of the  $6 \times 6$  matrix P(i)Y(i) be given by

$$\mathbf{P}(i)\mathbf{Y}(i) = \begin{pmatrix} Q_{11} & Q_{12} \\ Q_{21} & Q_{22} \end{pmatrix},$$

where  $Q_{ij} \in \mathbb{R}^{3\times 3}$ . Then

$$T_c(i) = -h^*(i)\mathcal{F}[Q_{11} + Q_{22}],$$
 (4.18)

where the mapping  $\mathscr{F}[\cdot]:\mathbb{R}^{3\times 3}\to\mathbb{R}^3$  is defined via the relation:

$$v = \mathcal{F}[A]$$
 if and only if  $\tilde{v} = A - A^*$ 

Proof. We have

$$\mathbf{P}(i)\mathbf{Y}(i)\mathbb{H}(i) = \begin{pmatrix} Q_{11}\tilde{h}(i) & Q_{22}\tilde{h}(i) \\ Q_{21}\tilde{h}(i) & Q_{22}\tilde{h}(i) \end{pmatrix}$$
(4.19)

Therefore.

Trace 
$$\{P(i)Y(i)H(i)\} = \text{Trace}\{(Q_{11} + Q_{22})\tilde{h}(i)\}.$$

Using the easily established identity that  $\operatorname{Trace}\{A\tilde{v}\}=-v^*\mathscr{F}[A]$  for an arbitrary matrix  $A\in\mathbb{R}^{3\times 3}$  and vector  $v\in\mathbb{R}^3$  in Eq. (4.19) establishes the result.

Notice that this expression for  $T_c(i)$  is vastly simpler compared with the original expression involving the mass matrix inverse and its sensitivity. Furthermore, the expression involves articulated body inertia quantities many of whom are available from the computational steps for solving the equations of motion. It is easy to verify that the

the compensating torque vector for the six base cluster degrees of freedom is zero.

## O(N) Constrained Dynamics Algorithm with Compensating Potential

The overall algorithm for solving the equations of motion with the compensating mass matrix potential is defined by the following steps:

- 1. Carry out the first base-to-tip part of the recursion in Eq. (3.1) to compute the V(k), a(k) and b(k) terms for all the clusters.
- 2. Carry out the tip-to-base recursion in Eq. (3.6) to compute all the articulated body inertia quantities such as **P**, **D** etc.
- 3. Carry out the base-to-tip recursion in Eq. (4.16) to compute  $\mathbf{Y}(k)$  for all the links. Also compute the compensating torque  $T_c(k)$  for each link simultaneously using Eq. (4.18) and T'(k) using Eq. (2.10).
- 4. Carry out the recursions in Eq. (3.9) and Eq. (3.10) to solve for the  $\ddot{\theta}(k)$  hinge accelerations with T(k) replaced with T'(k).

The only significant change from the constrained dynamics algorithm in [2] is the additional Step 3 for the compensating torque. This additional step is also of  $O(\mathcal{N})$  computational complexity, and hence, the overall computational cost of the constrained dynamics algorithm remains  $O(\mathcal{N})$ . Moreover this method adds only marginal cost since it makes use of the articulated body inertia quantities available from the computations for the regular solution for the internal coordinate accelerations. This algorithm makes possible the easy incorporation of the compensating potential into potential into constrained molecular dynamics simulations involving tree-topology molecular systems. This algorithm has been implemented as a part of the **NEIMO** software package.

#### Extensions

Reference [28] describes an alternative method for computing the  $\mathbf{Y}(k)$  terms using dual articulated body inertias instead of Eq. (4.16). This method offers advantages for parallel implementation since these computations can be done concurrently with Step 2 instead of sequentially following it.

The derivations and algorithmic descriptions in the previous sections have focused upon serial-chain molecules with one degree of freedom rotational hinges. This was done for notational simplicity. The extension of the algorithm to tree-topology molecules is identical to the extension discussed in Ref. [2]. In this case the system has multiple tips and a single designated base cluster. To summarize, all the tip-to-base and base-to-tip recursions are replaced

by tips-to-base and base-to-tips recursions. At each hinge with branches, the recursions proceed through "scatter" and "gather" steps. To handle multiple degree of freedom hinges, it is simply a matter of recognizing that such hinges can be modeled as a sequence of single degrees of freedom hinges interconnected by pseudo-clusters of zero mass and inertia. With these modifications, the algorithm described above extends to general tree-topology molecular systems.

Spatial Operator Expression for  $\mathcal{V}_c$ 

While explicit knowledge of the compensating potential  $\mathcal{V}_c(\boldsymbol{\theta})$  is not required during dynamics simulations, the nature and sensitivity of its dependence on the configuration coordinates can be used to study its effect on the statistical averages obtained from constrained dynamics simulations. The following lemma gives a closed form expression for the compensating potential  $\mathcal{V}_c$ .

Lemma 4.5.

$$\mathcal{V}_c(\boldsymbol{\theta}) = \frac{1}{2} \sum_{i=1}^n \ln \det \{ \boldsymbol{D}(i) \}.$$
 (4.20)

*Proof.* From the innovations factorization in Eq. (3.5a), it follows that

$$\det\{\mathcal{M}\} = \det\{I + \mathbf{H}\phi\mathbf{K}\}^2 \det\{\mathbf{D}\}$$

However, since  $[I + H\phi K]$  is lower triangular with identity matrix blocks along its diagonal, its determinant is 1. Hence

$$\det\{\mathcal{M}\} = \det\{\boldsymbol{D}\}.$$

The result of the lemma follows simply by taking the logarithm of the above identity. ■

The expression in Eq. (4.20) allows the easy computation of  $\mathcal{V}_c$  by making use of Algorithm 3.1 for computing the required  $\mathbf{D}(k)$  quantities.

#### 5. CONCLUSIONS

Internal rigid constraints, together with internal coordinates, are often used to speed up molecular dynamics computations. However, it has been recognized that the use of such constrained dynamics introduces systematic biases into the computation of statistical averages, and a compensating mass matrix potential is required to offset these biases. The lack of a tractable method for computing these compensating terms has rendered their use impractical to date.

This paper derives analytical expressions and algorithms for including the compensating mass matrix potential and its gradient into constrained MD simulations for general tree-topology molecular systems. The algorithm is an extension of the previously described  $O(\mathcal{N})$  algorithms for internal coordinate molecular dynamics simulations. The computational complexity of the new algorithm remains  $O(\mathcal{N})$ . Extension of this compensating potential algorithm to closed-topology molecular systems is the subject of ongoing research.

#### **ACKNOWLEDGMENTS**

The author thanks Dr. N. Vaidehi for reviewing the paper and for many helpful suggestions. The research described in this paper has been performed partially at the Jet Propulsion Laboratory, California Institute of Technology, under contract with the National Aeronautics and Space Administration and with support from NSF Grant ACS 92-17368.

#### REFERENCES

- A. Mazur and R. Abagyan, New methodology for computer-aided modeling of biomolecular structure and dynamics: 1. Non-cyclic structures, J. Biomol. Struct. & Dyn. 6(4), 815 (1989).
- A. Jain, N. Vaidehi, and G. Rodriguez, A fast recursive algorithm for molecular dynamics simulations, J. Comput. Phys. 106, 258 (1993).
- 3. M. Fixman, Classical statistical mechanics of constraints: A theorm and application to polymers, *Proc. Nat. Acad. Sci.* **71**, 3050 (1974).
- M. Fixman, Simulation of polymer dynamics. I. General theory, J. Chem. Phys. 69, 1527 (1978).
- N. Go and H. Scheraga, Analysis of the contribution of internal vibrations to the statistical weights of equilibrium conformations of macromolecules, J. Chem. Phys. 51, 4751 (1969).
- N. Go and H. Scheraga, On the use of classical statistical mechanics in the treatment of polymer chain conformation, *Macromolecules* 9, 535 (1976).
- W. van Gunsteren and M. Karplus, Effects of constraints on the dynamics of macro-molecules, *Macromolecules* 15, 1528 (1982).
- 8. M. Pear and J. Weiner, Brownian dynamics study of a polymer chain of linked rigid bodies, *J. Chem. Phys.* **71**, 212 (1979).
- E. Helfand, Flexible vs rigid constraints in statistical mechanics, J. Chem. Phys. 71, 5000 (1979).
- D. Chandler and B. Berne, Comment on the role of constraints on the conformational structure of *n*-butane in lid solventseiqu, *J. Chem. Phys.* 71, 5386 (1979).
- 11. W. van Gunsteren, Constrained dynamics of flexible of molecules, *Molecular Physics* **40**(4), 1015 (1980).
- A. Munster, Statistical Thermodynamics (Springer-Verlag, New York/Berlin, 1969).
- H. Berendsen and W. van Gunsteren, The Physics of Superionic Conductors and Electrode Materials (Plenum, New York, 1980), p. 221.
- M. Allen and D. Tildsley, Computer Simulation of Liquids. Clarendon Press. 1987.
- J. Mcammon and S. Harvey, *Dynamics of Proteins and Nucleic Acids*. Cambridge University Press, Cambridge, 1987.
- 16. A. M. Mathiowetz, A. Jain, N. Karasawa, and W. A. Goddard III, Protein simulations using techniques suitable for very large systems: The cell multipole method for nonbond interactions and the Newton– Euler inverse mass operator method for internal coordinate dynamics, in *Proteins: Structure Function Genetics* 20, 227 (1994).
- 17. A. Vaidehi, A. Jain, and W. Goddard, Constant temperature con-

- strained molecular dynamics: The Newton–Euler inverse mass operator method, *J. Phys. Chem.* **100**(25), 10508 (1996).
- 18. A. Graham, Kronecker Products and Matrix Calculus: With Applications, Halsted, New York, 1981.
- 19. G. Rodriguez, K. Kreutz-Delgado, and A. Jain, A spatial operator algebra for manipulator modeling and control, *The International Journal of Robotics Research* **10**, 371 (1991).
- A. Mazur, V. Dorofeev, and R. Abagyan, Derivation and testing of explicit equations of motion for polymers described by internal coordinates, *Journal of Computational Physics* 92, 261 (1991).
- 21. G. Rodriguez, Kalman filtering, smoothing and recursive robot arm forward and inverse dynamics, *IEEE Journal of Robotics and Automation* **3**, 624 (1987).
- G. Rodriguez and K. Kreutz-Delgado, Spatial operator factorization and inversion of the manipulator mass matrix, *IEEE Transactions* on Robotics and Automation 8, 65 (1992).

- 23. T. Kailath, The Innovations Approach to Detection and Estimation Theory, *Proceedings of the IEEE* **58**, 680 (1970).
- R. Featherstone, The calculation of robot dynamics using articulatedp-body inertias, *The International Journal of Robotics Research* 2, 13 (1983).
- B. D. O. Anderson and J. B. Moore, *Optimal Filtering*. Prentice-Hall Inc., 1979.
- 26. A. Jain and G. Rodriguez, Diagonalized lagrangian robot dynamics, *IEEE Transactions on Robotics and Automation* **11,** 571 (1995).
- 27. A. Jain and G. Rodriguez, An analysis of the kinematics and dynamics of underactuated manipulators, *IEEE Transactions on robotics and automation* **9**, 411 (1993).
- A. Jain and G. Rodriguez, Base-invariant symmetric dynamics of free-flying space manipulators, *IEEE Transactions on Robotics and Automation* 11, 585 (1995).