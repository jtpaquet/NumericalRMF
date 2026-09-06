

<!-- Page 1 -->

* [1] JULY 01 1999

A numerical study of rotating magnetic fields as a current drive for field reversed configurations \(\copyright\)

Richard D. Milroy



<!-- Page 2 -->

# A numerical study of rotating magnetic fields as a current drive for field reversed configurations

Richard D. Milroy

University of Washington, Redmond Plasma Physics Laboratory, Seattle, Washington 98195

17 February 1999; accepted 15 March 1999

###### Abstract

A fixed ion model has been developed to study the use of a Rotating Magnetic Field (RMF) as a current drive mechanism in a Field Reversed Configuration (FRC). This model is used to investigate the physics of RMF current drive in a parameter range of interest to two experiments at the University of Washington. Empirical expressions are found to characterize the critical RMF magnitude required for full penetration and the rate of RMF penetration. It is shown that in the presence of a strong anisotropic plasma resistivity, the direction and magnitude of the axial bias field can have a strong influence on the penetration of an RMF. Calculations that include the effects of realistic RMF antennae at finite radius are used to find the effects of coil spacing and positioning.

## I. INTRODUCTION
There are two experiments at the University of Washington designed to investigate the use of a Rotating Magnetic Field ~RMF! as an azimuthal current drive mechanism for a Field Reversed Configuration ~FRC!. In the star thruster experiment ~STX!,1 the RMF is being employed during the formation phase to reverse an initial forward bias field and form an FRC, using a methodology similar to previous Rotamak2 and FRC3,4 experiments. A major goal of this experiment is to investigate the use of RMF as an alternate formation method for FRCs. In the translation, confinement, and sustainment ~TCS!5 experiment, an RMF will be applied to an existing FRC. The plasma configuration will be formed conventionally and then translated into a sustainment chamber where an RMF will be applied. The primary goal of this experiment is to use an RMF to overcome the resistive loss of trapped flux, and demonstrate flux buildup in an existing FRC. To drive an azimuthal current in a cylindrical plasma, antennae are used to generate an external magnetic field of the form \( \bf B_\text{rmf} = B_\omega \cos(\omega t) \hat{\bf x} + B_\omega \sin(\omega t) \hat{\bf y} \). The frequency of the applied field is in the range \(\omega_{ci} \ll \omega \ll \omega_{ce} \), where \(\omega_{ci} = e B_\omega/m_i\) is the ion cyclotron frequency with respect to the rotating field strength, and \(\omega_{ce} = e B_\omega/m_e\) is the electron cyclotron frequency with respect to the rotating field strength. Under these conditions, the electrons are tied to the rotating field while the ions are unaffected, and an azimuthal electron current is generated. 

In the past, most experiments were designed with an RMF magnitude comparable to that of the induced axial field. Both of the University of Washington experiments are designed to demonstrate the induction and sustainment of axial magnetic fields significantly greater than the RMF current drive fields. To achieve this, higher temperatures ~lower plasma resistivity! must be reached. An objective of this paper is to extend previous theoretical6 and computational7,8 work into this new parameter range where the classical skin depth for the RMF is small compared to the plasma radius. A numerical model similar to those of Hugrass and Grimm,7 and Ohnishi et al.,8 has been developed and applied to this problem. Empirical relationships have been developed to predict the conditions required for RMF penetration and the rate of penetration. In addition, the effects of anisotropic plasma resistivity have been investigated. It has been found that while the imposed axial bias field has no effect when the assumed plasma resistivity is isotropic, it can have a large effect on the RMF penetration when the resistivity is aniso- tropic. Calculations that include the effects of RMF antennae at finite radius are used to find the effects of coil spacing and positioning.

A full three-dimensional magnetohydrodynamic ~MHD! code that includes the Hall term in its Ohm’s law is required to get a thorough understanding of RMF current drive in an FRC. This study is only a first step in that quest: however, this simplified model provides essential insight into some of the complex physical phenomena introduced when RMF fields are added to an FRC. The physical model employed in this study is discussed in Sec. II, along with a brief description of the numerical procedure. Section III presents the results of the study, fol-
lowed by a summary and conclusions in Sec. IV.

# II. PHYSICAL MODEL

A fixed ion computer model has been developed to study the interaction of an RMF with a cylindrical plasma. In this model, the plasma is represented as an infinitely long cylinder of radius R, surrounded by a vacuum. The model assumes the electrons move through a uniform density of stationary ions, and electron inertia is neglected.

## A. Model equations

With these assumptions, the relevant equations are

\[ \frac{\bf A}{\partial t} = \bf u \times \bf B - \eta \bf J - \frac{1}{en} [(\bf J \times \bf B - \nabla P_e )], \tag{1} \]

<!-- Page 3 -->

\[{\bf B} = \nabla\times{\bf A}, \tag{2}\] \[{\bf J} = \frac{1}{\mu_{o}}\nabla\times{\bf B}, \tag{3}\]

where \({\bf A}\) is the magnetic vector potential, \({\bf u}\) is the ion fluid velocity, \({\bf B}\) is the magnetic field strength, \({\bf J}\) is the electric current density, \(\eta\) is the plasma resistivity, and \(P_{e}\) is the electron pressure.

Equations (1) to (3) are solved numerically in the \(r-\theta\) plane assuming no axial variation (\(\partial/\partial z = 0\)). The ions are assumed to have a fixed background density \(n = N_{o}\), with zero velocity \({\bf u} = 0\). The electron pressure is assumed uniform (\(\nabla P_{e} = 0\)). It is convenient to Fourier expand the variables in \(\theta\) so that each quantity \({\bf Q}(\theta)\) can be expanded as

\[{\bf Q}(\theta) = {\bf Q}_{0} + \left(\sum_{n=1}^{N}\,{\bf Q}_{n}e^{in \theta} + {\rm c.c.}\right). \tag{4}\]

Boundary conditions at the outer radius, \(r = R\), must be applied to all three components of \({\bf A}\). The combination of boundary conditions on \(A_{r}\) and \(A_{\theta}\) are arranged to set \(B_{z}\) = \(B_{z0}\), where \(B_{z0}\) is the applied axial bias field independent of \(\theta\). The applied boundary conditions on \(A_{r}\) and \(A_{\theta}\) are

\[\frac{\partial A_{r,n-0}}{\partial r} = 0,\quad\frac{1}{r}\,\frac{\partial}{\partial r}\,(rA_{\phi,n-0}) = B_{z0},\quad A_{r,n-0}, \tag{5}\] \[= 0,\quad\frac{\partial}{\partial r}\,(rA_{\phi,n-0}) = 0.\]

The boundary conditions applied to \(A_{z}\) are more complex since they must accurately represent both the external rotating magnetic field and the effects of internal screening currents. It is assumed that there is no net axial current, so for the \(n = 0\) component the condition \(B_{\theta,n-0} = 0\), or equivalently, \((\partial/\partial r)A_{z,n-0} = 0\) is applied. Except for calculations with anisotropic resistivity, this boundary condition is not important since \(A_{z} = 0\)  for even harmonics, including \(n = 0\). For higher mode numbers the boundary conditions derived by Hugrass and Grimm [7] are applied. In the vacuum region around the FRC, \(A_{z}\) must satisfy the equation

\[\nabla^{2}A_{z} = 0,\quad r \geq R.\]

The solution to this equation can be expressed as

\[A_{z} = \sum_{n}\,A_{z,n}e^{in\theta},\,\,\,{\rm where}\,\,\,\,A_{z,n} =  \alpha_{n}(r/R)^{n} +\beta_{n}(r/R)^{-n},\]

where the first term (with the \(\alpha_{n}\) coefficient) accounts for the field produced by external currents, and the second term (with the \(\beta_{n}\) coefficient) accounts for the field produced by internal plasma currents. Since the magnetic field is continuous across the plasma boundary,

\[\frac{\partial A_{z,n}}{\partial r}\bigg|_{r=R} = \frac{n}{R}\,\alpha_{n } - \frac{n}{R}\,\beta_{n}\,.\]

The above two equations can be combined to yield [7]

\[A_{z,n} = 2\,\alpha_{n} - \left(\frac{R}{n}\right) \left.\frac{\partial A _{z,n}}{\partial r}\right|_{r=R}. \tag{6}\]

The numerical model has two options for representing the RMF. The first option, which is used for most of the calculations, assumes an idealized \(n = 1\)  single mode external field. The second option accurately accounts for discrete coils at finite radius. For the idealized single mode field, the magnitude of \(\alpha_{n}\) is adjusted to give the desired external RMF magnitude,

\[\alpha_{1} = \frac{RB_{\omega}}{2}e^{-i\omega t},\quad\alpha_{n\geq 2} =  0, \tag{7}\]

where \(B_{\omega}\) is the magnitude of the RMF in the vacuum region when the plasma is absent. This formulation correctly accounts for plasma screening currents, which tend to increase the vacuum field magnitude near the plasma boundary. In calculations that include the finite radius coil geometry of an experiment, the coefficients \(\alpha_{n\geq 2}\) are nonzero, as shown in the next section.

As an initial condition, \(B_{r} = 0\), \(B_{\theta} = 0\), and \(B_{z} = B_{z0}\) is applied throughout the plasma region. Except for calculations with an anisotropic resistivity, the solution is independent of \(B_{z0}\), and a value of \(0\) is used.

## B. Scaling

Hugrass [9] has shown that this system can be characterized by two important dimensionless parameters; \(\lambda = R/\delta\), where \(\delta = (2\,\eta/\mu_{o}\omega)^{1/2}\) is the classical skin depth, and by \(\gamma = \alpha_{ce}/\upsilon_{ei}\), where \(\upsilon_{ei} = \eta(ne^{2}/m_{e})\) is the electron-ion collision frequency. Expressed in terms of fundamental variables, these two dimensionless parameters are

\[\lambda = R\left(\frac{\mu_{o}\omega}{2\,\eta}\right)^{1/2}, \tag{8}\] \[\gamma = \frac{1}{e}\left(\frac{B_{\omega}}{n\,\eta}\right). \tag{9}\]

These parameters can be used to write Eqs. (1) to (3) in a dimensionless form. Each quantity \(q\) is transformed to a dimensionless variable \(\vec{q}\) using \(\vec{j} = q/q_{0}\). The scaling parameters are defined as \(t_{o} = 1/\omega\), \(r_{o} = R\), \(u_{o} = R\omega\), \(A_{o}\)\(= R\,B_{u}/\gamma_{u}\), \(B_{o} = A_{o}/R\), and \(J_{o} = (1/\mu_{o})B_{o}/R\). With these definitions, and neglecting the \(\nabla P_{e}\) term, Eqs. (1) to (3) can expressed as

\[\frac{\partial\vec{\Lambda}}{\partial\vec{t}} = (\vec{\bf u} \times \vec{\bf B})-\frac{1}{2\lambda^{2}}\,(\vec{ \bf J} + (\vec{\bf J} \times \vec{\bf B}))\,, \tag{10}\] \[\vec{\bf B} = \vec{\nabla} \times \vec{\Lambda}, \tag{11}\]

and

\[\vec{\bf J} = \vec{\nabla} \times \vec{\bf B}, \tag{12}\]

with a boundary condition

\[\vec{A}_{z1} + \frac{\partial\vec{A}_{z1}}{\partial\vec{r}} = \gamma e^{-i\vec{t}}. \tag{13}\]

If the RMF is fully penetrated and the electrons are rotating at an angular velocity of \(\omega\), the induced swing in the axial field between \(r = 0\) and \(r = R\) is

<!-- Page 4 -->

\[\Delta B_{z}=\frac{\mu_{o}ne\omega}{2}R^{2}.\]

In dimensionless units, this can be expressed as \(\Delta\vec{B}_{z}=\lambda^{2}\). This means that for an FRC with zero pressure at the separatrix, the axial field swings from \(+\lambda^{2}/2\) at \(r=R\), to \(-\lambda^{2}/2\) at \(r=0\).

Clearly, the problem can be completely defined through the specification of \(\lambda\) and \(\gamma\), where \(\lambda\) is a measure of the diffusivity, and \(\gamma\) is a measure of the RMF strength.

## C. Numerical method

The equations are solved numerically using a finite difference mesh in the \(r\)-direction and by Fourier expanding in the \(\theta\)-direction. This technique is found to be very efficient, especially in situations where few azimuthal mode numbers are required to accurately represent the model. For most of the problems studied in this paper, only the \(n=1\) mode needs to be included, since the calculated magnitudes of higher order modes are found to be small. However, an arbitrary number of modes can be included, and for some of the nonlinear calculations that include anisotropic resistivity or realistic coil geometry, up to 42 azimuthal modes are included.

Temporally, the equations are advanced using a predictor-corrector algorithm, which has two basic options. The first option is a straightforward explicit algorithm that can advance each time step very quickly, but with a very restrictive time-step stability condition. The second option is a semi-implicit algorithm that is unconditionally stable for any time-step size, but takes somewhat longer for each time step. Accuracy requirements demand a relatively small time step for the semi-implicit algorithm.

# III. Results

## A. Penetration condition and penetration rates

Hugrass and Grimm[7] first studied the time-dependent penetration of the RMF into a plasma column with a model very similar to that presented here. They showed that the magnitude of the external-driving field must exceed a certain critical value before much current is driven in the plasma column. Their measure of the current drive is the parameter \(\alpha\), which is defined as the ratio of the driven current divided by the current obtained if all of the electrons rotate synchronously with the RMF. From Ampere's law the actual driven current per unit length is given by

\[\overline{I_{\theta\theta}}=\frac{1}{\mu_{o}}(B_{z}(r=0)-B_{z}(r=R)).\]

If all the electrons were driven synchronously with the RMF, the current density per unit length would be

\[\overline{I_{\theta\theta}}=-1/2\pi e\omega R^{2}.\]

Thus, \(\alpha\) can be expressed as

\[\overline{\alpha}=\frac{2}{\mu_{o}ne\omega R^{2}}(B_{z}(r=R)-B_{z}(r=0)).\]

Hugrass and Grimm[7] found that the steady-state value of \(\alpha\) (expressed as \(\alpha_{x}\)) is quite small for values of \(\gamma\) below a threshold value, but when \(\gamma\) exceeds this threshold, \(\alpha_{x}\) jumps to a value near 1. Later, Hugrass[9] extended these results, showing the existence of nonunique solutions to the steady-state problem. The fundamental meaning of these results is that for values of \(\lambda\) greater than about 6, the RMF will not fully penetrate until \(\gamma\) (or equivalently \(B_{\omega}\)) exceeds a given threshold value, Nevertheless, once it has penetrated, \(B_{\omega}\) can be reduced significantly below this value before the RMF will be expelled and efficient current drive lost.

Figure 1 shows the calculated evolution of magnetic field lines into a plasma column with parameters relevant to the STX experiment (\(R=10\) cm, \(n=0.333\times 10^{20}\) m\({}^{-3}\), \(B_{\omega}\) = 100 G, \(\omega=2.2\times 10^{5}\) sec\({}^{-1}\), \(T_{e}=5\) eV, \(\lambda=11.07\), and \(\gamma\) = 16.6). Here time is labeled in the dimensionless unit \(T\) = (2 \(\pi\))/\(\omega\), or 2.856 \(\mu\)sec. The RMF is rotating in a counter-clockwise direction and the field lines are distorted as they attempt to drag the electrons with them. For this case, the steady normalized azimuthal current is found to be \(\alpha_{x}\) = 0.98, with full penetration achieved at \(t=40T\).

If this calculation is repeated with the RMF magnitude reduced by 10% to \(B_{\omega}=90\) G, or \(\gamma=14.9\), we find that complete penetration is never achieved and the steady normalized azimuthal current is only \(\alpha_{x}=0.42\). The equilibrium field lines for this case look very much like those at \(t=5\,T\), in Fig. 1. If, on the other hand, the RMF magnitude is increased by 10% to \(B_{\omega}=110\) G, or \(\gamma=18.2\), we find that complete penetration is achieved by \(t=27T\), compared to \(t\) = 40\(T\) for the case with \(B_{\omega}=100\) G. This is illustrated graphically in Fig. 2, which shows \(\alpha\) as a function of time during penetration for these three calculations. Clearly, the magnitude of the RMF must exceed a critical value before the field can fully penetrate the plasma column, and after this critical field is reached, the rate of penetration is a strong function of the field strength.

Figure 1: Evolution of magnetic field lines as the RMF penetrates a plasma column.

If the RMF is reduced below a critical value after it has fully penetrated, the force on the electrons will not be sufficient to overcome the electron-ion drag. Then the current



<!-- Page 5 -->

drive level will drop back to a level consistent with an RMF that is not fully penetrated. Figure 3 shows the magnetic field lines for the same calculation illustrated in Fig. 1, except at 
\(t=50T\), \(B_\omega\) is suddenly reduced from 100 to 70 G (\(\gamma\) is reduced to 11.6). As the RMF field magnitude drops in the plasma column, the rate of electron rotation slows. Since the magnetic field lines are tied to the electron fluid, field lines rotate more slowly at the column center than at the outside. This leads to an antiparallel configuration causing quick field annihilation, and a rapid reduction in current drive.

There are three key questions concerning RMF penetration that this study attempts to answer. First, what is the critical RMF magnitude required for full penetration, and how does it scale with \(\gamma\)? Second, what is the time required to achieve full penetration when the RMF magnitude is above its critical value? Finally, what efficiency can be achieved when the RMF magnitude is below the critical value and the current drive is restricted to a region near the separatrix?

A procedure similar to that of Hugrass9 is used to find the critical RMF magnitude. A calculation is started with small subcritical RMF magnitude. It is run long enough to establish a constant \(\alpha_s\) value for this field, after which the field is incremented to the next level where a new \(\alpha_s\) is established. This is done repeatedly until the field level reaches the critical value and full penetration is achieved.
After this, the field level can be decremented in steps until it
reaches a level where the RMF is expelled from the plasma,
establishing the critical value for field expulsion. Figure 4
shows \(\gamma_c\), which is the critical value of \(\gamma\) required to achieve
complete RMF penetration, plotted as a function of \(\lambda\). These
results agree with Hugrass’ observation9 that for values of \(\lambda\) less than about 6.5, the critical value of
\(\gamma\) required to achieve
penetration is the same as that for expulsion. It is found to be
approximately


\[\gamma_c=1.12\lambda, \qquad \text{for } \lambda\leq 6.5. \tag{14}\]

For larger values of \(\lambda\) it has been found that the following
empirical formula is a good fit to the numerical calculations


\[\gamma_c=1.12\lambda (1.0 + 0.12(\lambda-6.5)^{0.4}), \qquad \text{for } \lambda\gt 6.5. \tag{14}\]

The point at which the field is expelled from the plasma
remains the same \( \gamma_c=1.12\lambda \) for all values of \(\lambda\). For large
values of \(\lambda\), a larger RMF magnitude is required to achieve
complete penetration than is required to maintain it after pen-
etration is achieved.

An estimate of the time required for the RMF to penetrate \(\tau_P\) has been made by running several calculations with
various values of \(\gamma \g \gamma_c\). From the form of Eq. ~10!, it is
expected that \(\tau_P\) should scale as \(\lambda^2\). Since
\(\tau_P\) becomes infinite at \(\gamma=\gamma_c\), it is expected that
\(\tau_P\) should vary with a normalized form of \(\gamma\) that is a measure of the difference
between \(\gamma\) and \(\gamma_c\). This measure is defined as

\[\gamma_N \equiv \frac{\gamma-\gamma_c}{\gamma_c}. \tag{16}\]

Figure 5 shows the penetration time vs \(\gamma_N\) for three series of calculations with \(\lambda=100\), \(\lambda=11\), and \(\lambda=6.5\). The numerical
data are compared with the simple curve

\[\tau_P = \frac{\lambda^2}{2\sqrt{\gamma_N}}, \tag{17}\]

which is seen to be a reasonable approximation to the data
over a very wide range of
\(\gamma\) and \(\lambda\). Clearly, when designing
an experiment that requires full penetration, the design point
\(\gamma\) should be significantly greater than
\(\gamma_c\).

Next, we examine the steady-state level of current drive
that can be achieved when the RMF magnitude is below the
critical value. Several calculations have been made where \(\lambda\)
is held constant while
\(\gamma\) is slowly incremented in steps. The



<!-- Page 6 -->

time between each increment is long enough to allow the current drive parameter \(\alpha\) to reach its steady-state value \(\alpha_{s}\). The results of this series of calculations are summarized in Fig. 6, where curves for \(\lambda\) = 6.5, 11, 32, and 100 are shown. The dotted lines in this figure show an empirical fit to the numerical data,

\[\alpha_{s} = \frac{1.6}{\sqrt{\lambda}}\exp \left(-4 \left(\frac{\gamma_{ c}-\gamma}{\gamma_{c}}\right)^{8}\right) . \tag{18}\]

This fit is accurate to approximately 20% for \(.2 \in  \gamma/\gamma_{c}\) \(\approx .95\). Clearly, for large values of \(\lambda\), the current drive efficiency is very low for \(\gamma  < \gamma_{c}\).

## B. Anisotropic resistivity effects

In all of the calculations reported above, an isotropic resistivity has been assumed. For a hot plasma with classical resistivity, the ratio of perpendicular to parallel resistivity is \(\eta_{\perp}\) / \(\eta_{\parallel}\)= 1.96. In FRC experiments, the transport is observed to be consistent with an anomalous perpendicular resistivity larger than classical. However, in the absence of a toroidal magnetic field, the plasma current is perpendicular to the field and no measurement of the parallel component of resistivity has been made. For other toroidal plasmas it is found that the ratio of \(\eta_{\perp}\) / \(\eta_{\parallel}\) is considerably greater than 2.

To account for anisotropic resistivity in this study, Eq. (10) is modified to

\[\frac{\partial\mathbf{\overline{A}}}{\partial\mathbf{\overline{t}}} = ( \mathbf{\overline{u}} \times \mathbf{\overline{B}})-\frac{1}{2\lambda^{2}}  \left(\mathbf{J}_{\perp}+\frac{\eta_{\parallel}}{\eta_{\perp}} \mathbf{\overline{J}}_{\parallel} + (\mathbf{\overline{J}} \times  \mathbf{\overline{B}})\right) , \tag{18}\] 

where \(\mathbf{\overline{J}}_{\parallel} = (\mathbf{\overline{J}} \cdot \mathbf{ \overline{B}})/\mathbf{\overline{B}}^{2}\), \(\mathbf{\overline{J}}_{\perp} = \mathbf{\overline{J}} - \mathbf{\overline{ J}}_{\parallel}\), and \(\lambda =R(\mu_0 \omega/2\eta_\perp)^{1/2}\). Figure 7 shows the time history of the azimuthal current, for
a set of calculations with different forms of anisotropic re-
sistivity. All of the calculations have the same basic parameters (\( R=10~\text{cm}, B_\omega=120~\text{G}, \omega =2.2\times 10^6~\text{sec}^{-1}, T_e=5~\text{eV}, n=0.333 \times 10^{20}~\text{m}^{-3}, \lambda=11.07, \gamma=19.2, \) and \( \tilde{B}_{z0}=\lambda^2/2 \)) similar to those of Fig. 1, but with a
higher value of \(B_\omega\), and consequently \(\gamma \). Also, it should be
noted that since the resistivity is now dependent on the mag-
netic field direction, these calculations are now dependent on
the axial bias field. Initially \( B_{z} = B_{z0} \) uniformly throughout
the calculation region, and subsequently this value is applied
as a boundary condition at \(r=R\). All of the calculations in
Fig. 7 have the same value for \(\eta_\perp\), but \(\eta_\parallel\) varies from 1 to \(\frac{1}{4}\) of that value. For two of the curves in Fig. 7 the parallel
direction is assumed to be the \(z\)-direction, while the others
correctly employ the local direction of \(\bf{B}\) to specify the par-
allel and perpendicular directions. When the parallel direc-
tion is assumed to be in the \(z\-direction, decreasing\(\eta_\parallel\) leads
to a slower rate of penetration as expected,7 since the resis-
tive skin depth decreases with decreasing resistivity. How-
ever, it is found that when the local value of \(\bf{B}\) is used to
specify the parallel direction, the penetration rate does not
slow with decreasing \( \eta_\parallel\), and in fact increases.

To understand these results consider the schematic in
Fig. 8. Assume \(B_z\) is positive, pointing out of the page. This
field combined with the RMF is upward and pointing out of
the page. The RMF field ~shown as vertical lines in Fig. 8! is
rotating counterclockwise, dragging the electrons in the same
direction. If \(\bf{J}\) tends to follow the field lines, it flows into the
page on the right-hand side ~RHS!, and out of the page on
the left-hand side ~LHS!. This is the same direction as the
current in the RMF coils, and will therefore enhance the
rotating magnetic field. If \(B_z\) is negative ~pointing into the
page!, the effect is reversed and the induced axial current
will oppose the penetration of the RMF.

An axial bias field of \(\tilde{B}_{z0}=\lambda^2/2\) was assumed for the
calculations of Fig. 7. For these conditions and when the


<!-- Page 7 -->

RMF is fully penetrated, the induced current sustains an axial magnetic field profile with \(\tilde{B}_{z} = \lambda^{2}/2\) at the outer boundary, decreasing to \(- \lambda^{2}/2\) by \(r = 0\). Thus in the outer regions where \(B_{z}\) is positive, the anisotropic resistivity leads to the induction of axial currents which enhance the RMF. In the inner region where \(B_{z}\) is reversed, the induced axial currents oppose the penetration of the RMF.

Clearly, with the inclusion of anisotropic resistivity effects, the calculations become dependent on the magnitude and direction of the applied axial bias field. If a bias field of zero is applied (\(B_{z0}  = 0\)), the induced axial field is negative so the anisotropic resistivity induces axial currents which oppose the RMF throughout the entire plasma volume.

Figure 9 shows the steady normalized current as a function of \(\gamma\) for three different parameter sets. In all three curves \(\lambda = 11.07\). These curves are generated by running the code at a constant \(\gamma\) until \(\alpha\) reaches a steady-state value. The parameter \(\gamma\) is then changed by a small amount and the code is again run until \(\alpha\) reaches a new constant value. The curves are double-valued because the full penetration threshold value of \(\gamma\) is larger when \(\gamma\) is increasing than the threshold value where the RMF is expelled when \(\gamma\) is decreasing. The curve labeled \(\eta_\parallel  = \eta_\perp\) corresponds to the parameters of Fig. 1. As expected, the curve with \(B_{z0}  = 0\) shows that a much higher value of \(\gamma\) is required to achieve full penetration, and it displays much less hysteresis. The curve with \(\vec{B}_{z0}  = \lambda^{2}/2\) has a wider hysteresis profile, as a higher value of \(\gamma\) is required to achieve full penetration when \(\gamma\) is increasing, and \(\gamma\) drops to a lower value before the RMF is expelled when it is decreasing.

Figure 10 shows the steady-state radial profile of the RMF field for the three calculations of Fig. 9 at points on the curves with full penetration. The \(y\)-axis, which is labeled \(B_{\theta}/B_{\omega}\), is the normalized magnitude of the \(n = 1\) component of \(B_{\theta}\). The curve labeled \(\vec{B}_{z0}  = \lambda^{2}/2\) shows an amplification of \(B_{\theta}\) relative to the boundary value as it penetrates in the outer regions where \(B_{z}\) is positive. As it penetrates further into regions of negative \(B_{z}\), the magnitude of \(B_{\theta}\) decreases. On the other hand, the curve labeled \(\vec{B}_{z0}  = 0\) rapidly and monotonically decreases from the boundary value as it penetrates the plasma. In this case, the induced \(B_{z}\) is always negative, so the axial currents induced by the anisotropic resistivity always oppose the penetration of the RMF field. The curves do not begin (at the outer radial boundary) with a value of 1 as may be expected, because the applied boundary conditions (see Sec. II) correctly account for the effects of the induced internal plasma currents. When the internal currents screen (or oppose) the penetration of the RMF they amplify the RMF at the boundary. On the other hand, if the internal currents enhance the penetration of the RMF, they reduce the magnitude of the RMF at the boundary. This is also illustrated in Fig. 11, which shows the calculated steady-state magnetic field lines for the same two anisotropic resistivity calculations.

Figure 8: Schematic of plasma current in the presence of the rotating magnetic field.

Figure 11: Steady-state magnetic field lines for the anisotropic calculations illustrated in Fig. 10.

Figure 9: Steady-state current drive as a function of \(\gamma\).

Figures 10 and 11 show that when an FRC bias field (\(\tilde{B}_{z0}  = \lambda^{2}/2\)) is applied, a net amplification of the RMF is


<!-- Page 8 -->

realized. The amplification in the outer regions is significantly greater than the shielding of the inner region. The net effect is that the rotating field at \(r = 0\) is more than 20% greater than the external field far removed from the plasma, and more than 50% greater than the rotating field at the plasma boundary. It is interesting to note that Blevin and Thonemann[3] observed a 500% amplification of the rotating field in some of their pioneering experiments. They attributed the amplification to three-dimensional end effects but it is possible that anisotropic resistivity effects also contributed.

In Fig. 9, the curve with \(\vec{B}_{z0} = \lambda^{2}/2\) has a relatively flat region where \(\alpha_{s} \simeq \).5 for 15\(< \gamma < \)21. In this region, the character of the RMF penetration is markedly different from that of previous cases. It is found that RMF penetration proceeds normally at the start of the calculation, with the axial currents induced by the anisotropic resistivity enhancing the rate of penetration. However, after the RMF has penetrated about half way in, the magnetic field lines begin to tear and a closed field-line structure forms inside the plasma. This structure rotates with the RMF, but not as fast. Figure 12 shows the calculated evolution of induced closed field-line structure during one revolution of the RMF for a calculation with \(\gamma = 16.6\). (Here, the structure appears to rotate clockwise, but this is because in this paper the graphs are plotted in a rotating frame of reference that makes the external RMF appear stationary.) When the RMF is in this mode, it develops an \(\alpha_{z}\) of about 0.5. The electron velocity profile is approximately that of a rigid-rotor, but there is significant slippage at the outer radial boundary. Figure 13 shows the \(n\)\(= 0\) component of the \(\theta\)-current profile at the same four times that are illustrated in Fig. 12. This figure shows that the current profile varies significantly during a single rotation of the RMF field, and that there is significant slippage at the outer radial boundary. Presumably, the magnetic field structures that form inside the plasma force an almost rigid rotation of the electrons. When there is full penetration, the current profile is a straight line (rigid-rotor) as illustrated by the "Full Penetration" line in Fig. 13.

As shown in Fig. 9, when \(\gamma\) is increased beyond 21, the RMF has sufficient power to overcome this mode, leading to normal full penetration. It is interesting to compare the curve (iv) of Fig. 7 with the \(\eta_{ {}_{1}}/\,\eta_{ {}_{2}} = \frac{1}{4}\), \(B_{z0} = \lambda^{2}/2\) curve of Fig. 9. It has the same basic parameters and its value of \(\gamma = 19.2\) is below the threshold for penetration in Fig. 9, but here penetration is even faster than it was for the comparable isotropic case. The key difference is that this calculation started out with \(\gamma = 19.2\) and the mode structures illustrated in Fig. 12 never formed, while the other calculations started with a small value of \(\gamma\) which was then slowly increased. Once the structures have formed, a larger value of \(\gamma\) is required to achieve full penetration.

## C. Finite radius coil effects

The applied boundary conditions for all of the calculations reported above correspond to an ideal pure \(n\)\(= 1\) RMF. Hugrass[10] has studied the effects of realistic coil geometry, but using a simplified model where the electron motion is forced to have a rigid-rotor velocity profile. He showed that in addition to the desired \(n = 1\) RMF, coils at finite radius produce undesired odd spatial harmonics. These undesirable effects can all be minimized by moving the coils out to a larger radius, by using two conductor coil pairs, or by using a three-phase system to generate the RMF. Here, we extend these results to examine the influence of these added spatial harmonics on RMF penetration, and without the rigid electron motion constraint. In general, we find agreement with Hugrass' earlier work.

In a vacuum, the vector potential due to a set of infinitely long dipole coils arranged on a cylindrical surface of radius \(R_{c}\) can be expanded as[10]

\[A_{z} = \frac{\mu_{o}I}{4\,\pi} \sum_{n =  \text{odd}}\frac{1}{n}  \left[\frac{r}{R_{c}}\right]^{n}\cos(n(\theta -\varphi)), \tag{20}\]

where \(\varphi\) is the angular position of the coil set, and \(I\) is the current flowing in them. In the presence of \(I_{c}\) polyphase coil sets, each with current \(I_{i} = I\sin(\omega t - \xi_{i})\), the vector potential is

\[A_{z} = \frac{\mu_{o}I}{4\,\pi} \sum_{n = \text{odd}}\left[\frac{1}{n}  \left(\frac{r}{R_{c}}\right)^{n} \sum_{i=1}^{I_{c}}\,\left[\cos(n(\theta  -\varphi_{i}))\cdot\sin(\omega t - \xi_{i})\right]\right] . \tag{21}\]


<!-- Page 9 -->

Thus finite coil radius effects can be included in the calculation by specifying \(\alpha_{n}\) in Eq. (6) as

\[\alpha_{n} = \frac{\mu_{o}I}{8\,\pi}\bigg(\frac{r}{R_{c}}\bigg)^{n}\frac{1} {n}\sum_{i=1}^{I_{c}}\ \{\sin(\omega t-\xi_{i})e^{-in\varphi_{i}}\},\ \ n=\text{odd},\] \[\alpha_{n} = 0,\ \ \ n=\mbox{even}. \tag{22}\]

The simplest RMF coil arrangement consists of a two-phase system with a single coil pair per phase and each coil pair separated by \(\pi/2\) in space and phase. If \(\varphi_{1} = 0\) and \(\varphi_{2}\) = \(\pi/2\), and \(\xi_{1} = \pi\) and \(\xi_{2} = -\pi/2\), the vacuum vector potential can be expressed as [10]

\[A_{z} = -\frac{\mu_{o}I_{o}}{\pi}\sum_{n= \text{odd}}\bigg[\frac{1}{n} \bigg(\frac{r}{R_{c}}\bigg)^{n}\sin(\omega t+(-1)^{(n+1)/2}n\,\theta) \bigg].\]

The resulting field consists of the desired \(n = 1\) component rotating with an angular frequency \(\omega\), as well as odd harmonics which rotate at a frequency of \((-1)^{(n-1)/2}(\omega/n)\). It has been found that these have a harmful effect on the current drive and need to be minimized. Increasing the coil radius \(R_{c}\) is an obvious way to reduce the relative magnitude of the harmonics. However, this is not practical experimentally since the required energy to reach the same magnitude for the fundamental \(n = 1\) field scales with \(R_{c}^{2}\).

Using two coils for each phase of the antenna can also reduce the harmonics. For a two-phase system with two coils per phase, the vacuum vector potential can be expressed as [10]

\[A_{z} = -\frac{2\,\mu_{o}I_{o}}{\pi}\sum_{n= \text{odd}}\Bigg[\frac{\cos ( \frac{1}{2}n\,\alpha_{c})}{n}\bigg(\frac{r}{R_{c}}\bigg)^{n} \tag{24}\]

\[\times\sin(\omega t+(-1)^{(n+1)/2})n\,\theta)\Bigg],\]

where \(\alpha_{c}\) is the angular separation between two coils in the same phase. This expression shows that each of the higher order modes is reduced by the factor \(\cos((1/2)n\,\alpha_{c})/\cos((1/2)\,\alpha_{c})\) relative to the primary \(n = 1\) component. Choosing \(\alpha_{c} = 60^{\circ}\) eliminates the third harmonic, while choosing \(\alpha_{c} = 36^{\circ}\) eliminates the fifth harmonic. It will be shown that \(\alpha_{c} = 45^{\circ}\) is a good compromise, in agreement with previous studies. [10]

The harmonics can be reduced even further using a three-phase system with two coils per phase. Here the vacuum vector potential can be expressed as [10]

\[A_{z} = -\frac{3\,\mu_{o}I_{o}}{\pi}\sum_{n= \text{odd}}\Bigg[\frac{ \cos( \frac{1}{2}n\,\alpha_{c})}{n}\Bigg(\frac{r}{R_{c}}\Bigg)^{n} \tag{25}\]

\[\times\gamma_{n}\sin(\omega t+\beta_{n}n\,\theta)\Bigg],\]

where \(\gamma_{1} = 1\), \(\gamma_{3} = 0\), \(\gamma_{5} = 1\), and \(\gamma_{n+6} = \gamma_{n}\), and \(\beta_{1} = 1\), \(\beta_{5}\) = \(-1\), and \(\beta_{n+6} = \beta_{n}\). For this configuration the 3rd, 9th, and 15th harmonics are not present and the 5th, 11th, and 17th harmonics rotate counter to the primary. The fifth harmonic can be eliminated by setting the coil separation to \(36^{\circ}\), leaving the seventh as the first harmful harmonic.

Figure 14 shows the azimuthal current time history for three calculations with \(\lambda = 11.07\), and the simplest coil arrangement (two-phase with a single coil set per phase). The ratio \(k = r/R_{c}\) varies from 0.7 to 0.5 for these calculations. Equations (20) to (25) show this ratio is important because the boundary values of \({\bf A}\), and hence \({\bf B}\), for the higher order modes are attenuated by the factor \((r/R_{c})^{n}\). The three calculations with \(k = 0.7\) have \(\gamma = 19.9\), 24.9, and 33.2, respectively, all well above the critical value of \(\gamma_{c} = 15.1\). Full penetration is never achieved for the curve labeled (i), while the curve labeled (ii) with \(\gamma = 24.9\) and \(k = 0.7\) is very close to the critical \(\gamma_{c}\). Both this curve and curve (iii), which has a much larger \(\gamma = 33.2\), exhibit large oscillations in the driven current. For the curve labeled (iv), the coil radius has been increased, reducing \(k\) to 0.5, and \(\gamma = 19.9\) has been reduced back to the level of curve (i). Here the RMF penetration is relatively quick and the oscillation magnitude is reduced. The curve labeled (v) is for comparison and is from a calculation where ideal \(n = 1\) boundary conditions are applied (i.e., finite radius coil effects are not included). Experiments with a two-phase antenna system and a single coil per phase should have the antenna at a large enough radius to ensure that the parameter \(k\) is not greater than 0.5.

Figure 15: Magnetic field lines during \(\frac{1}{4}\) revolution of RMF from (iii) of Fig. 14 starting at \(t = 70.25T\).


Figure 15 shows the magnetic field lines at eight times during a quarter rotation of the RMF for the calculation corresponding to curve (iii) of Fig. 14. The plots are drawn in a


<!-- Page 10 -->

rotating frame of reference so the ideal \(n=1\) component of the RMF is vertical for all eight plots, and the coils appear to rotate. Clearly, the field lines inside the plasma region continue to oscillate after the RMF has fully penetrated. Furthermore, the electron velocity profile does not evolve to the stationary rigid-rotor profile observed in calculations that neglect effects of having coils at finite radius. Figure 16 shows the radial profile of the \(n=0\) component of \(J_{\theta}\) from the same calculation and times as Fig. 15. The velocity profile changes rapidly with relatively large excursions from the ideal rigid-rotor profile, with a periodicity of \(\frac{1}{4}\) of an RMF cycle. The profile oscillations are not confined to the edge of the plasma and the magnitudes do not dampen as you move in radially from the outer boundary where the driving force exists.

If two coils are used for each phase of the antenna, the harmonic content diving these oscillations is reduced[10], as shown by Eq. (24). Figure 17 shows the azimuthal current time history for three calculations with a coil separation specified by \(\alpha_{c}=36\), \(45\), and \(60^{\circ}\). The calculations at \(36\) and \(60^{\circ}\) suppress the \(n=5\) and \(n=3\) modes, respectively. The calculation with \(\alpha_{c}=45^{\circ}\) is a compromise where both the magnitudes of both \(n=3\) and \(n=5\) modes are partially suppressed by a factor of \(0.414\). All three calculations have \(\lambda\) = \(11.07\), \(\gamma=19.9\), and \(k=0.7\); the same parameters as curve (i) of Fig. 14, where complete penetration was never achieved. Clearly, employing two conductors per phase yields a dramatic performance enhancement, but the results are insensitive to the exact angle between conductors. The calculation with \(\alpha_{c}=45^{\circ}\) has the most rapid penetration, with the penetration time history approaching that of the ideal \(n=1\) calculation. However, the calculation with \(\alpha_{c}\) = \(60^{\circ}\) has the smallest oscillations in driven current.

Figure 18 shows the azimuthal current time history for a similar set of calculations, except the coils are assumed to be closer to the plasma, by increasing the parameter \(k\) from \(0.7\) to \(0.8\). As expected, this has a deleterious effect on the RMF performance. Only the calculation with \(\alpha_{c}=45^{\circ}\) achieves complete penetration, and for that calculation, the rate of penetration is much slower than it was for \(k=0.7\), as shown in Fig. 17. Experiments with a two-phase antenna and two coils per phase should have the antenna at a large enough radius to ensure that the parameter \(k\) is smaller than \(0.7\).

Finally, a set of calculations has been made to investigate the effects of employing a three-phase antenna that uses two coils per phase. Equation (25) shows that for this type of antenna, the \(3\)rd, \(9\)th, \(15\)th,... harmonics are absent, so the \(5\)th and \(7\)th are the first harmful harmonics. If the coil separation parameter \(\alpha_{c}\) is set to \(36^{\circ}\), the \(5\)th harmonic is eliminated, leaving the \(7\)th as the first harmful harmonic. If \(\alpha_{c}\) is set to \(30^{\circ}\), both the \(5\)th and \(7\)th harmonics are attenuated by a factor of approximately \(0.27\). Hugrass[10] found that for a three-phase system with two coils per phase with \(\alpha_{c}\) set to \(30^{\circ}\) and \(k=0.8\), higher harmonics had a negligible effect on the system. Figure 19 shows the azimuthal time history for a set of calculations with \(\lambda = 11.07\), \(\gamma = 19.9\), \(k = 0.8\) and 0.9, and \(\alpha_{c} = 30\) and \(36^{\circ}\). For \(k = 0.8\) the penetration rate is almost as fast as for the ideal \(n = 1\) calculation, in stark contrast to the similar two-phase calculations of Fig. 18. When \(k\) is increased to 0.9, however, the penetration rate slows considerably. The calculations with \(k = .9\) show a clear advantage of \(\alpha_{c} = 30^{\circ}\), over \(\alpha_{c} = 36^{\circ}\). It is better to attenuate both the 5th and 7th harmonics, rather than eliminate the 5th harmonic without attenuating the 7th. Figure 20 illustrates the evolution of the magnetic field lines during penetration for the calculation with \(k = 0.8\), and \(\alpha_{c} = 30^{\circ}\).

Figure 16: Radial profile of the \(n=0\) component of \(J_{\theta}\) for calculation illustrated in Fig. 15.

Figure 17: Azimuthal current time history for a two-phase antenna with two coils per phase, and \(\lambda=11.6\), \(\gamma=19.9\), and \(k=0.7\).

Figure 18: Azimuthal current time history for a two-phase antenna with two coils per phase, and \(\lambda=11.6\), \(\gamma=19.9\), and \(k=0.8\).

Figure 19: Azimuthal current time history for a three-phase antenna with two coils per phase, and \(\lambda=11.6\), and \(\gamma=19.9\).

It is noted that for calculations with larger values of \(k\), higher order harmonics become more important and must be included in the numerical solutions. For these calculations, all harmonics up to \(n = 21\) are included for \(k = 0.8\), and up to \(n = 42\) are included for \(k = 0.9\). Harmonics up to \(n = 10\) are included for calculations with \(k \approx 0.7\).

# IV. Summary and Conclusions

The penetration of an RMF into a plasma column has been studied with a numerical model. Several empirical formulas have been established by curve-fitting the results of multiple calculations. In particular these empirical relationships predict (1) the conditions required for RMF penetration, (2) the point at which the RMF field will be expelled from a plasma column, (3) the time it takes for an RMF field to penetrate, and (4) the current drive efficiency when the penetration condition is not met.

It was found that an anisotropic resistivity can have a strong effect on RMF current drive for an FRC, affecting both the penetration and final equilibrium configuration. The underlying physical mechanism was shown to result from an induced axial current. This current can be in the same direction (phase) as the current in the external RMF coils and thus enhance the RMF, or it can be in the opposite direction and retard the penetration of the RMF. The phase of the induced current depends on the direction of the local axial magnetic field \(B_{z}\), so RMF penetration becomes strongly dependent on the imposed bias on the axial magnetic field. With an axial bias field consistent with an FRC, the magnitude of the RMF is significantly amplified at the center (\(r = 0\)) as compared to its magnitude at the outer plasma boundary.

The effects of using a realistic coil set at a finite radius have been investigated. These results are in agreement with a previous study by Hugrass[10] based on an equilibrium model that assumed a rigid-rotor electron velocity profile: however, dynamic effects and the effects on penetration were also studied in this paper. It was found that the harmonics introduced by utilizing coils at finite radius can seriously impede the RMF penetration unless the coils are sufficiently far from the plasma. For a two-phase coil system that utilizes two coils per phase, the coil radius should be at least 1.4 times the plasma radius to minimize these effects. The optimum angle between the two coils in a single phase is found to be about \(45^{\circ}\). If a three-phase coil system with two coils per phase is employed, it is found that the coils can be as close as 1.25 times the plasma radius before the higher harmonics displays serious effects on RMF penetration. Here the optimal angle between two coils of the same phase is about \(30^{\circ}\).

# Acknowledgements.

 The author acknowledges many useful discussions with J. T. Slough, A. L. Hoffman, and S. A. Cohen during the course of this work. The present work has been supported by grants from the Office of Fusion Energy Sciences of the U.S. Department of Energy.

## References

* (1) K. E. Miller, J. T. Slough, and A. L. Hoffman, _Space Technology and Applications International Forum_, AIP Conf. Proc. No. 420, Part 3 (AIP, New York, 1998), p. 1352.
* (2) W. N. Hugrass, I. R. Jones, K. F. McKenna, M. G. R. Phillips, R. G. Storer, and H. Tuczek, Phys. Rev. Lett. **44**, 1676 (1980).
* (3) H. A. Blevin and P. C. Thonemann, Nucl. Fusion Suppl. 55 (1962).
* (4) A. J. Knight and I. R. Jones, Plasma Phys. Controlled Fusion **32**, 575 (1990).
* (5) A. L. Hoffman, Phys. Plasmas **5**, 979 (1998).
* (6) I. R. Jones and W. N. Hugrass, J. Plasma Phys. **26**, 441 (1981).
* (7) W. N. Hugrass and R. C. Grimm, J. Plasma Phys. **26**, 455 (1981).
* (8) M. Ohnishi, A. Ishida, Y. Yamamoto, and K. Yoshikawa, Trans. Fusion Technol. **27**, 391 (1995).
* (9) W. N. Hugrass, Aust. J. Phys. **38**, 157 (1985).
* (10) W. N. Hugrass, Aust. J. Phys. **39**, 513 (1986).

Figure 20: Evolution of RMF field lines during penetration for calculation with a three-phase antenna, two coils per phase, \(\lambda = 11.6\), and \(\gamma = 19.9\), \(k = 8\), and \(\alpha = 30^{\circ}\).