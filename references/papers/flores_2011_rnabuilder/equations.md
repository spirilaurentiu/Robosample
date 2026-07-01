# Equations - Flores 2011, RNABuilder / Simbody base-pairing force field

Reduced/coarse-grained RNA modeling in internal coordinates. The implementable
core is the base-pairing force-torque potential (eqs 5-11). Units: ps, nm,
kJ/mol, Daltons (inherited from Molmodel).

<!-- eq:1 -->
$$ \dot{q}=N(q)\,u, \qquad M(q)\,\dot{u}=f(t,q,u) $$
- **what:** Unconstrained internal-coordinate multibody equations of motion (ODE form).
- **symbols:** q - generalized coordinates (nq-vector); u - generalized speeds (n-vector); N(q) - nq×n block-diagonal kinematic coupling matrix; M(q) - n×n composite system mass matrix; f(t,q,u) - n generalized forces (applied + Coriolis); t - time.

<!-- eq:2 -->
$$ \dot{q} = N(q)\,u, \qquad M(q)\,\dot{u} = f(t,q,u) - G^{T}\lambda, \qquad g(t,q) = 0 $$
- **what:** Constrained multibody EOM (index-3 DAE); adds m constraint equations and constraint forces.
- **symbols:** g(t,q) - m constraint equations; G = ∂g/∂q (m×nq); λ - Lagrange multipliers (unknown constraint forces, m-vector); other symbols as eq:1.

<!-- eq:3 -->
$$ {}^{A1}R^{B2} = {}^{A1}R^{G} \cdot {}^{G}R^{B2} = \left( {}^{G}R^{A1} \right)^{-1} \cdot {}^{G}R^{B2} $$
- **what:** Rotation aligning body frame B2 (residue 2) onto attachment frame A1 (residue 1); the residual misalignment rotation whose angle θ is penalized.
- **symbols:** ${}^{X}R^{Y}$ - rotation matrix (3×3, SO(3)) expressing frame Y in frame X; G - ground/global frame; A1 - attachment frame on residue 1; B1, B2 - body frames of residues 1, 2.

<!-- eq:4 -->
$$ {}^{G}R^{A1} = {}^{G}R^{B1} \cdot {}^{B1}R^{A1} $$
- **what:** Attachment-frame orientation in ground = body-frame orientation times the constant attachment-in-body rotation.
- **symbols:** ${}^{B1}R^{A1}$ - constant orientation of attachment frame A1 in residue 1's body frame (from the base-pairing model); ${}^{G}R^{B1}$ - body-frame orientation vs ground, function of q.

<!-- eq:5 -->
$$ U(r,\theta) = \left[ \frac{\theta^{2} \cdot \kappa}{2 \cdot k} + 1 \right] \cdot g(r, k, c) \cdot m, \qquad -\pi < \theta < \pi $$
- **what:** Base-pairing potential energy; couples the angular misalignment θ (from eq:3, Euler theorem) and the translational distance r between A1 and B2.
- **symbols:** U - potential energy (kJ/mol); r - distance between A1 and B2 (nm); θ - Euler rotation angle of ${}^{A1}R^{B2}$ (rad, in (-π,π)); κ - angular stiffness constant (per interaction type, typically positive); k - radial depth constant (per interaction type, typically negative); c - radial range (global, nm); m - global scaling factor (forceMultiplier); g - radial shape function (eq:6).

<!-- eq:6 -->
$$ g(r,k,c) = \begin{cases} -\dfrac{k \cdot r^{2}}{2 \cdot c^{2}} + \dfrac{3 \cdot k}{2}, & 0 \leq r < c \\[2mm] \dfrac{k \cdot c}{r}, & r \geq c \end{cases} $$
- **what:** Radial shape function; harmonic in r for r<c, inverse-r decay for r>=c. Continuous at r=c (both branches give k·c... note: short branch at r=c gives -k/2+3k/2 = k; long branch gives k). Value = k at inflection r=c.
- **symbols:** as eq:5. At r=c, g=k (the inflection-point / well-depth value).

<!-- eq:7 -->
$$ \vec{r} = r \cdot \hat{r} = \vec{x}_{B2} - \vec{x}_{A1} $$
- **what:** Displacement vector from attachment frame A1 origin to body frame B2 origin; defines r and unit direction r-hat.
- **symbols:** $\vec{x}_{A1}$, $\vec{x}_{B2}$ - positions of A1 and B2 origins (nm); $\hat{r}$ - unit vector along $\vec{r}$; r = |$\vec{r}$|.

<!-- eq:8 -->
$$ \vec{F} = -\vec{\nabla}U = -\frac{\partial U}{\partial\theta} \cdot \hat{\theta} - \frac{\partial U}{\partial r} \cdot \hat{r} = -\frac{\theta \cdot \kappa}{k} \cdot g(r, k, c) \cdot m \cdot \hat{\theta} - \left[\frac{\theta^{2} \cdot \kappa}{2 \cdot k} + 1\right] \cdot g'(r, k, c) \cdot m \cdot \hat{r} $$
<!-- CHECK: OCR wrote the second partial as ∂/∂τ; corrected to ∂/∂r from the explicit r-derivative result on the RHS. -->
- **what:** Force from the potential = negative gradient, split into angular (θ-hat, torque-like) and radial (r-hat) components.
- **symbols:** $\vec{F}$ - force; $\hat{\theta}$ - unit vector of the angular (rotation) direction; g' - radial derivative dg/dr (eq:9); other symbols as eq:5.

<!-- eq:9 -->
$$ g'(r,k,c) = \begin{cases} -\dfrac{k \cdot r}{c^{2}}, & r < c \\[2mm] -\dfrac{k \cdot c}{r^{2}}, & r \geq c \end{cases} $$
- **what:** Derivative dg/dr of the radial shape function (used in the radial force component).
- **symbols:** as eq:5/eq:6.

<!-- eq:10 -->
$$ \vec{f}_{A1} = \left[ \frac{\theta^{2} \cdot \kappa}{2 \cdot k} + 2 \right] \cdot g'(r, k, c) \cdot m \cdot \hat{r} = -\vec{f}_{B2} $$
<!-- CHECK: constant is "+2" here vs the "+1" radial coefficient in eq:8; transcribed as printed in the paper. Sign of overall term differs from eq:8's radial part ("we dropped the negative sign because of the sense of r-hat"). -->
- **what:** Translational (radial) force applied at attachment frame A1; equal and opposite force at B2. This is the net pairing pull between the two bases.
- **symbols:** $\vec{f}_{A1}$ - translational force on residue-1 attachment frame; $\vec{f}_{B2} = -\vec{f}_{A1}$ - reaction on residue 2; other symbols as eq:5.

<!-- eq:11 -->
$$ \vec{\tau}_{A1}^{*} = \frac{\theta \cdot \kappa \cdot m}{k} \cdot g(r, k, c) \cdot m \cdot \hat{\theta} + (\vec{x}_{A1} - \vec{x}_{O1}) \times \vec{f}_{A1} $$
$$ \vec{\tau}_{B2}^{*} = \frac{\theta \cdot \kappa \cdot m}{k} \cdot g(r, k, c) \cdot m \cdot \hat{\theta} - (\vec{x}_{B2} - \vec{x}_{O2}) \times \vec{f}_{A1} $$
<!-- CHECK: numerator shows κ·m and the trailing ·m, i.e. an m^2 factor; the angular-force term from eq:8 has a single m (θκ/k · g · m). One m is likely OCR duplication. Transcribed as printed. -->
- **what:** Adjusted torques applied at each base body: the angular (θ-hat) pairing torque plus the moment-arm correction, because the force is applied at the body origin O rather than at A1/B2.
- **symbols:** $\vec{\tau}_{A1}^{*}$, $\vec{\tau}_{B2}^{*}$ - adjusted torques on bodies 1, 2; $\vec{x}_{O1}$, $\vec{x}_{O2}$ - body-origin positions of bases 1, 2 (nm); $\times$ - cross product; other symbols as eq:5/eq:10.
