---
paper_id: "2000_Milroy"
title: "* [1]A magnetohydrodynamic model of rotating magnetic field current drive in a field-reversed configuration \\(\\copyright\\)"
authors: ""
year: "2000"
institution: ""
publisher: ""
doi: "10.1081/1.1290279"
keywords: ""
abstract: ""
num_references: "14"
num_citations: ""
source_file: "2000 Milroy - A magnetohydrodynamic model of rotating magnetic field current drive in a field-reversed configuration.pdf"
---



<!-- Page 1 -->

* [1]A magnetohydrodynamic model of rotating magnetic field current drive in a field-reversed configuration \(\copyright\)
* [2]Richard D. Mitroy

* [3]Crack for updates

* [4]Phys. Plasmas 7, 4135-4142 (2000)

[https://doi.org/10.1081/1.1290279](https://doi.org/10.1081/1.1290279)


<!-- Page 2 -->

A magnetohydrodynamic model of rotating magnetic field current drive in a field-reversed configuration

Richard D. Milroy

University of Washington, Redmond Plasma Physics Laboratory, Seattle, Washington 98195

24 April 2000; accepted 5 July 2000

A numerical model has been developed to study the use of a Rotating Magnetic Field (RMF) as an electron current drive mechanism for the formation and sustainment of a field-reversed configuration (FRC). Previous models assumed a fixed ion model, but here a full two-dimensional (_r_- \(\theta\)) magnetohydrodynamic model has been developed. The model has been applied to two classes of problems: (1) For the sustainment problem, a RMF is applied to a preexisting FRC. (2) For the formation problem, a RMF is applied to a plasma column with an initially uniform axial magnetic field and background plasma density. The RMF-induced current reverses this bias field, forming a FRC. The code employs an option to include some three-dimensional effects to satisfy the average \(\beta\) condition and equalize pressure and density between inner and outer field lines, when it is applied to sustainment simulations. (c 2000 American Institute of Physics.

[S1070-664X(00)03610-7]

## I Introduction

A transverse rotating magnetic field (RMF) has been used to generate toroidal currents in a plasma column.[1] This mechanism has been studied extensively in both the Rotamak[2] and the field-reversed configuration (FRC).[3, 4] Jones has recently compiled a review of this work.[5] The RMF field penetration and the resultant current drive have been studied analytically and numerically.[6, 7, 8, 9, 10]

Previous numeric studies have not included fluid flow or plasma profiles. In this paper, an existing code,9 which solved for the time evolution of the magnetic field in a plasma with stationary background ions, has been extended through the addition of the ion momentum and temperature equations of form a fully two-dimensional (r-theta) magneto-hydrodynamic ~MHD! model. In addition, the model has been adapted to optionally include some of the effects of flow parallel to the magnetic field lines ~in the z direction! and axial pressure balance for a FRC. In this paper we focus on calculations with parameters based on the results of the Star Thrust Experiment ~STX! at the University of Washington.4

Here, the numerical model is applied to two types of calculations. First, the formation of a FRC from an initial uniform background plasma and axial magnetic field is examined. In these calculations, the evolution of the configuration is followed as a RMF-induced toroidal current reverses the axial magnetic field, creating a FRC.

In a separate set of calculations, the sustainment of an existing FRC is examined. In this case, the calculation is initialized with a density, temperature, and a magnetic field profile corresponding to a FRC. Then the evolution profiles is calculated, as a new equilibrium is sought with RMF current drive. For these calculations, the (r-theta) MHD model incldes some three-dimensional (3-D) effects by equalizing pressure on the inner and outer field lines, and by satisfying the condition [11] \[ <\beta> = 1-0.5x_s^2\]. With a suitable choice of resistivity, equilibrium is predicted with only partial penetration of the RMF, as observed in the STX experiments [4]. In fact, the RMF only penetrates to near the magnetic field null, providing just enough current drive to make \[E_\theta\] zero there. A steady radial inward flow of the fluid provides the current drive on the inner field lines. It is assumed that the fluid is transported from the inner field lines back to the outer field lines by a steady flow parallel to the
magnetic field (primarily in the z direction).

A description of the MHD model used for these studies can be found in Sec. II. The results of the numerical calculations and a discussion of the relevant physics are given in Sec. III. Finally, the implications of these results and plans for future work are discussed in Sec. IV.

## II Computation Model

### A. The basic MHD model

A two-dimensional r-tehta MHD model has been developed to study the evolution of the RMF and plasma profiles. This model is an extension of a previous model9 that assumed a uniform density of stationary ions. The MHD equations for this model are

\[\frac{\partial\mathbf{A}}{\partial t}=\mathbf{u}\times\mathbf{B} -\,\eta\mathbf{J}-\,\frac{1}{en}\bigg[(\mathbf{J}\times\mathbf{B} )-\boldsymbol{\nabla}P_{e}\bigg], \tag{1}\]

\[\frac{\partial n}{\partial t}+\boldsymbol{\nabla}\cdot n\mathbf{u}=0, \tag{2}\]

\[Mn\bigg(\frac{\partial\mathbf{u}}{\partial t}+\mathbf{u}\cdot\boldsymbol {\nabla}\mathbf{u}\bigg)=(\mathbf{J}\times\mathbf{B})-\boldsymbol{ \nabla}\mathbf{P}-\boldsymbol{\nabla}\cdot\Pi, \tag{3}\]

\[ \frac{3}{2} n k_B T \left( \frac{\partial T}{\partial t} + \mathbf{u} \cdot \nabla T \right) = \eta \mathbf{J}^2 - P \nabla \cdot \mathbf{u} + \nabla \cdot (\kappa \nabla T) - \Pi : \nabla \mathbf{u} - R_a, \tag{4}\]

\[ \mathbf{B} = \nabla \times \mathbf{A}, \tag{5} \]


\[ \mathbf{J} = \frac{1}{\mu_0} \nabla \times \mathbf{B}, \tag{6} \]

\[ P=n k_B T, \tag{7}\]

\[ P_e = f_{TE} P, \tag{8}\]

\[ \pi = -\mu(\nabla \mathbf{u} + \nabla \mathbf{u}^T + \frac{2}{3}\nabla \cdot \mathbf{u}), tag{9} \]

\[ R_a = 10^{-31} f_{imp} n^2, \tag{10} \]



<!-- Page 3 -->


where \({\bf A}\) is the magnetic vector potential, \(n\) is the number density, \({\bf u}\) is the ion fluid velocity, \(M\) is the ion mass, \({\bf B}\) is the magnetic field strength, \({\bf J}\) is the electric current density, \(\vec{\eta}\) is the plasma resistivity, \(P\) is the total plasma pressure, \(P_{e}\) is the electron pressure, \(k_{\rm B}\) is Boltzmann's constant, \(\mu\) is the viscosity, \(R_{a}\) is the power loss due to radiation, and \(f_{\rm imp}\) is the assumed concentration of impurities in the plasma.

It is assumed that the plasma is confined within an infinitely long cylinder of radius \(R\), surrounded by a vacuum. Equations (1)-(9) are solved numerically in the \(r\)-\(\theta\) plane, assuming no axial variation (\(\partial/\partial z\)\(=\) 0). The variables are Fourier expanded in \(\theta\), and finite difference techniques are used to solve for the Fourier coefficients. The numerical techniques and boundary conditions for the magnetic field equations are discussed in an earlier paper,[9] however, the boundary condition for the \(n\) = 0 component of \(A_{\theta}\) (or \(B_{z}\)) has been modified to correspond to that of a flux conserved at a specified radius. In this paper, the previous model is extrapolated to include fluid flow and time varying plasma profiles, but the magnetic field solver is essentially unchanged.

A small artificial viscosity coefficient \(\mu\), is used in advancing Eq. (3) to provide numerical stability. The thermal conductivity coefficient \(\kappa\), in Eq. (4), is chosen to keep the temperature profile relatively flat, as observed in current experiments. (It would be incorrect to assume a classical \(\kappa_{e\perp}\) here because the RMF opens the field lines in the region that it has penetrated.)

While the numerical model is capable of including an arbitrary number of Fourier modes, only the \(n\) = 0 and \(n\) = 1 modes have been retained. It was shown earlier[9] that this is sufficient, since the magnitude of the higher-order modes is small. In addition, the \(n\) = 0 component of \(u_{\theta}\) is zeroed. This prevents a _spin-up_ of the ions that is not observed in experiments, presumably due to charge exchange[12]--a process not included in this model.

### B. Phenomenological inclusion of 3-D effects

Two physical effects that are very important to equilibrium FRCs are missing from an \(r\)-\(\theta\) plasma model. First, if it is confined by a flux-conserving coil, a FRCs length will adjust to produce an average \(\beta\) defined by[11] (\(\beta\)) = 1 - \(\frac{1}{2\chi_{s}^{2}}\), where \(x_{s}\) = \(r_{s}/R_{c}\), \(r_{s}\) is the separatrix radius, and \(R_{c}\) is the flux conserving coil radius. Second, since the plasma is free to flow parallel to the magnetic field lines, the plasma pressure is equalized between inner and outer field lines on an Alfven time scale. In addition, strong parallel thermal conductivity rapidly reduces temperature differences between inner and outer field lines. Thus, pressure, temperature, and density should rapidly equalize between inner and outer field lines.

The numerical model has optional switches; one adds a term that forces the solution to evolve toward satisfying the average \(\beta\) condition; the other adds terms to equalize the density, temperature, and pressure between inner and outer field lines.

A FRC maintains the average \(\beta\) condition by expanding or contracting axially. If \(\langle\beta\rangle\) is less than 1 - \(\frac{1}{2\chi_{s}^{2}}\), it will shrink axially. This has the effect of increasing the density and temperature (and pressure) at the axial midplane, and consequently increasing the \(\langle\beta\rangle\). Conversely, if \(\langle\beta\rangle\) is greater than 1 - \(\frac{1}{2}\chi_{s}^{2}\), the FRC will expand axially. To account for these effects in the code, Eq. (2) has been modified to read as

\[\frac{\partial n}{\partial t}+\vec{\nabla}\cdot n{\bf u}=\dot{n}_{\beta} +\dot{n}_{\parallel}\,, \tag{11}\]

where \(\dot{n}_{\beta}\) is a contribution from the change in FRC length that moves the solution toward satisfying the average \(\beta\) condition, and \(\dot{n}_{\parallel}\) is the contribution from flow along field lines to equalize pressure on inner and outer field lines. The term \(\dot{n}_{\beta}\) is defined as

\[\dot{n}_{\beta} = \frac{n}{\gamma\tau_{\beta}}\bigg(\frac{\langle\beta_{D}\rangle -\langle\beta\rangle}{\langle\beta\rangle}\bigg), \tag{12}\]

on the closed field lines (\(\psi\)\(<\) 0), and is zero on the open field lines (\(\psi\)\(\approx\) 0). \(\psi\) is the magnetic flux function defined as the \(n\) = 0 component of \(rA_{\theta}\). Here, \(\langle\beta_{D}\rangle\) = 1 - \(\frac{1}{2\chi_{s}^{2}}\) is the desired value of \(\langle\beta\rangle\), and \(\tau_{\beta}\) is an input parameter that defines the time scale over which the FRC length will change.

The term \(\dot{n}_{\parallel}\) is defined differently on the inner field lines than on the outer field lines. On the inner field lines it is defined as

\[\dot{n}_{\parallel} = \frac{1}{2\,\gamma\tau_{\parallel}}\bigg(\frac{n_{o}T_{o}-n _{i}T_{i}}{\frac{1}{2}(T_{i}+T_{o})}\bigg), \tag{13}\]

where the subscripts "\(i\)" and "\(o\)" designate whether the value is measured on an inner field line or on the corresponding outer field line, and \(\tau_{\parallel}\) is an input parameter defining the time scale over which the pressure difference is equalized. The rate of change of density on the outer field lines is defined as

\[\dot{n}_{\rm oil} = -\frac{|B_{zo}|}{|B_{z\parallel}|}\dot{n}_{\parallel}\,. \tag{14}\]

The sign difference between outer and inner field lines results from the fact that if this process is increasing density in one region, it must be taking it from the other. The magnetic field ratio term accounts for the fact that the cross-sectional area of a magnetic flux bundle is inversely proportional to the local field strength.

The temperature equation is modified as

\[\dot{T} = \cdots+(\gamma-1)\,\frac{T}{n}(\dot{n}_{\beta}+\dot{n}_{ \parallel})+\dot{T}_{T}, \tag{15}\]


<!-- Page 4 -->

where the first term is from adiabatic changes with density flow, and the second term, \(\dot{T}_{T}\), is from parallel thermal conduction. It is also defined differently on inner and outer field lines,

\[\dot{T}_{iT}=\frac{(T_{o}-T_{i})}{2\,\tau_{T}}\quad\mbox{and}\quad \dot{T}_{oT}=-\frac{n_{i}}{n_{o}}\frac{B_{o}}{B_{i}}\dot{T}_{iT}, \tag{16}\]

where \(\tau_{T}\) is the time scale over which the temperature difference is equalized.

## III Results

### FRC formation

A FRC can be formed through the application of a RMF to a cold plasma column with a positive embedded axial magnetic field.[4] The plasma is preionized prior to the application of the RMF. The RMF induces an azimuthal current in the plasma sufficient to reverse the magnetic field, thereby forming a FRC. In this section, the results of numerical calculations to simulate this process are reported.

These calculations are initialized with a uniform positive \(B_{z}\). The initial density is uniform, and the temperature is set to 1 eV. A plasma resistivity of the maximum of classical, or 25 \(\mu\)m\(\,\)m is employed. Cooling due to impurity radiation with \(f_{\rm imp}=0.03\) is assumed. A temperature clamp limiting temperature to 50 eV, a typical value inferred from the STX results[4] is also applied. Without this clamp, the simulation temperature rises to a significantly higher value. The reason for this discrepancy is not understood, and will be discussed further in Sec. IV.

The initial conditions for the _standard_ calculation reported in this section include a uniform \(n=5\times 10^{18}\,\)m\({}^{-3}\), and \(B_{z}=30\,\)G. The rotating magnetic field has a vacuum strength of \(B_{\omega}=20\,\)G, and rotates at a frequency of \(\omega_{\rm mrf}=2.2\,\)\(\times\,10^{6}\,\)rad per second. The stated strength of \(B_{\omega}=20\,\)G is the magnitude that the RMF would have if there was no plasma. In fact, the strength of the \(n=1\,\) component of \(B_{\theta}\) is greater than 20 G, due to the effects of screening currents in the plasma, as explained elsewhere.[9] No corrections are made to include 3-D effects. These effects cannot be included since the initial conditions do not correspond to a FRC, as required for the average \(\beta\) condition to apply.

The time history of several key quantities (\(B_{z}\),\(r_{s}\),\(r_{\Delta\phi}\),\(T_{\rm ave}\), \(\phi\)) for this calculation is shown in Fig. 1. The boundary condition conserving flux on the \(\theta\)-pinch coil (at 22.6 cm radius) causes the external axial magnetic field to rise from the initial 30 G bias to over 100 G by 200 \(\mu\)s. The induced reverse flux remains zero until a sufficient reverse field is induced to compensate for the forward bias field, but almost 0.2 mWb of flux is induced by 350 \(\mu\)s. The excluded flux radius is seen to rise rapidly to about 15 cm, even though there is insufficient reversal to form a separatrix until 175 \(\mu\)s. After this time, \(r_{\Delta\phi}\), and \(r_{s}\) are in close agreement. The temperature rises, almost linearly to its final clamped value of 50 eV in about 175 \(\mu\)s.

Radial profiles of several key variables are shown in Fig. 2. The axial magnetic field profile shows the RMF rapidly starts to penetrate the initially cold resistive plasma, but as the plasma heats and becomes more conductive, the penetration rate slows. Initial axial field reversal is achieved at 20 \(\mu\)s, but the field does not reverse at \(r=0\) until 350 \(\mu\)s. The

Figure 1: Time history of (a) external field, (b) trapped flux, (c) excluded flux radius and separatrix radius, and (d) temperature, for the _standard_ formation calculation.

Figure 2: Profiles of (a) axial magnetic field, (b) plasma density, (c) plasma pressure, and (d) \(\alpha\) at 10, 50, 150, and 300 \(\mu\)s, for the _standard_ formation calculation.


<!-- Page 5 -->

density plot shows that the plasma is initially compressed inward radially, but as it heats and the pressure increases it moves back out. The \(\alpha\)=\(u_{e,\theta}\)/(\(\omega_{\rm{mf}}r\)) curve shows the RMF penetrating radially into the plasma column, even as flow slowly moves the plasma out. The parameter \(\alpha\) is 1, where the electrons are rotating synchronously with the RMF. This plot shows that \(\alpha\) is almost exactly 1 in the outer region, where the RMF field magnitude is large, but then falls off rapidly. This calculation clearly illustrates the time scale and profiles that may be expected when RMF is utilized to induce a reverse field and form a FRC.

It is anticipated that 3-D effects play an important role in the formation of an equilibrium FRC, and could become important shortly after field reversal is achieved at about 50 \(\mu\)s. A full three-dimensional code is required to accurately represent these effects during the formation process, and is beyond the scope of this paper.

### Sustainment of an existing FRC

The calculations of the previous section show how the application of a RMF to a plasma column with an embedded axial bias field can lead to field reversal. However, these calculations cannot predict the formation of an equilibrium FRC since the 3-D effects described in Sec. II.2 are not applied. In this section, we examine the use of RMF to sustain a FRC. The calculations are initialized with a \(B_{z}\) and \(n\) profile of an equilibrium FRC. Furthermore, the FRC length is adjusted to maintain the average \(\beta\) condition and adjustments are made to maintain \(n\), \(T\), and \(P\) as a function of \(\psi\), as described in Sec. II.2. If RMF boundary conditions are not applied, the configuration decays monotonically, because the initial trapped flux is lost through resistive diffusion. However, if an appropriate RMF boundary condition is applied, this decay stops and flux growth is observed once the transverse field penetrates to the magnetic null.

The initial pressure profile is specified as \(P(u)\) = \(P_{o}e^{-a|a|b}\), where \(P_{o}\)= \(B_{oz}^{2}\)/2\(\mu_{0}\),\(B_{zo}\) is the initial external axial magnetic field strength, \(a\) and \(b\) are profile adjusting parameters, \(u\) = (\(r/r_{0}\))\({}^{2}\)\(-\)1, and \(r_{0}\) is the radius at the field null. The initial axial magnetic field is specified to produce pressure balance.

The initial conditions for the _standard_ numerical simulation of this section are a uniform temperature \(T_{eo}\) = 20 eV, external axial magnetic field \(B_{zo}\)= 80 G, and \(x_{s}\) = 0.8. The profile parameters are set to \(a\) = 2, and \(b\) = 2. The simulation region extends from \(r\)= 0 to \(r\) = 21.5 cm, and a flux conserving coil is at \(R_{c}\) = 22.6 cm. The RMF field ramps on in 1 \(\mu\)s, to \(B_{\rm{mf}}\)= 20 G, and has a frequency of \(\omega_{\rm{mf}}\) = 2.2\(\times\)10\({}^{6}\) s\({}^{-1}\). The relaxation parameters \(\tau_{\beta}\), \(\tau_{1}\), and \(\tau_{T}\) are all set to 1 \(\mu\)s. A uniform resistivity of 25 \(\mu\)\(\Omega\) m is employed, and cooling due to impurity radiation with \(f_{\rm{imp}}\) = 0.03 is assumed. Like the formation calculations, the plasma temperature is limited to values less than the experimentally observed 50 eV.

The time history of the external axial magnetic field \(B_{ze}\), the trapped flux \(\phi\), and the FRC length \(L_{\rm{FRC}}\) are plotted in Fig. 3. Note: the parameter \(L_{\rm{FRC}}\) is a length relative to the initial starting length. For the first 25 \(\mu\)s the trapped flux decays, displaying a flux lifetime of about 70 \(\mu\)s. At about this time the RMF reaches the magnetic null and the current drive generates a growth in flux, and after some initial oscillations, equilibrium is reached.

The radial profiles of some key parameters after equilibrium is reached are shown in Fig. 4. The magnetic field profile in Fig. 4(a) shows that the RMF barely penetrates to the magnetic null. The \(B_{z}\) field profile has a very distinctive flat region around the null--a feature observed in the STX experiments.[12]

Figure 4: Equilibrium profiles for (a) magnetic field, (b) radial velocity, (c) pressure, and (d) normalized azimuthal current and \(\alpha\), for the _standard_ sustainment calculation.

Figure 3: Time history of (a) external field, (b) trapped flux, and (c) FRC length for the _standard_ sustainment calculation.


<!-- Page 6 -->

It is found that the equilibrium radial velocity, shown in Fig. 4(b), is uniformly negative, or inward. This profile can be sustained with an unchanging density profile, due to the contributions of \(n_{\parallel}\) in Eq. (11). Particles are convected radially inward from the outer field lines to the inner field lines by \(u_{r}\). This leads to a small pressure imbalance, with the inner field line pressure higher than the outer field line pressure. This, in turn, leads to an axial flow toward the FRC ends on the inner field lines, and a corresponding flow toward the midplane on the outer field lines. Equations (11), (12), and (13) permit this behavior in the two-dimensional (\(r\)-\(\theta\)) MHD model.

This equilibrium flow pattern is an essential feature of RMF current drive for this plasma configuration. As shown in Fig. 4(a), the RMF field only penetrates to the magnetic null. To achieve equilibrium, \(E_{\theta}^{0}=u_{r}^{0}B_{z}^{0}+\eta J_{\theta}^{0}\)\(+\langle J_{z}^{1}B_{r}^{1}\rangle/(\epsilon n)\) must be zero, or equivalently the \(n=0\) component of the \(\theta\) component of the right-hand side of Eq. (1) must be zero. The superscript "0" or "1" is used to denote the \(n=0\) or \(n=1\) components, respectively. The \(\eta J_{\theta}^{0}\) term is negative everywhere, but it can be balanced by the \(n=0\) component of \((J_{z}^{1}\timesB_{r}^{1})\). However, since the RMF only penetrates to the null, both \(J_{z}^{1}\) and \(B_{r}^{1}\) vanish on the inner field lines. Here the \(\eta J_{\theta}^{0}\) term is balanced by the radial flow through the (\(u_{r}^{0}\timesB_{z}^{0}\)) term, which is positive on the inner field lines. This term is negative on the outer field lines, so a larger contribution from the \(n=0\) component of \((J_{z}^{1}\timesB_{r}^{1})\) is required there.

The equilibrium pressure profile, shown in Fig. 4(c), also has a distinct flat top, as required for pressure balance, since \(B_{z}\) has a flat region around the null. Figure 4(d) shows the normalized current profile, and the parameter \(\alpha\)\(=u_{e}/(\omega_{\rm{mfl}}\tau)\). This parameter is almost exactly 1 in the outer region where the RMF field magnitude is large, but falls to a relatively small value near the magnetic null. It rises again on the inner field lines, as it must, since pressure is maintained as a function of flux. The fluctuations on \(\alpha\) at small \(r\) reflect small numerical errors in equalizing pressure between inner and outer field lines. The algorithm used is accurate to first order.

A penetration condition is often written in terms of two dimensionless quantities \(\gamma\) and \(\lambda\), where, \(\gamma=\omega_{ee}/v_{ei},\omega_{ee}\) is the electron cyclotron frequency in the RMF field, \(\nu_{ei}\) is the electron-ion collision frequency, \(\lambda=R/\delta\), \(R\) is the plasma radius, and \(\delta=(2\,\eta/\mu_{e}\omega)\) is the classical skin depth. These two parameters can be expressed as

\[\gamma = \frac{1}{e}\left(\frac{B_{\omega}}{n\,\eta}\right), \tag{17}\] \[\lambda = r_{s}\left(\frac{\mu_{e}\omega}{2\,\eta}\right)^{1/2}. \tag{18}\]

Previous calculations showed [9] that full RMF penetration is expected when \(\gamma>\gamma_{c}\), where

\[\gamma_{c}=1.12\lambda[\,1.0+0.12(\lambda-6.5)^{0.4}],\quad\mbox{for } \lambda>6.5. \tag{19}\]

This condition was derived for a fixed ion model with a uniform plasma density profile. For the equilibrium parameters of this calculation, \(\lambda=47\), \(\gamma_{c}=80\), and \(\gamma=53\). These parameters are calculated using the peak density for \(n\), and a \(B_{\omega}\) corresponding to the vacuum value that would exist if the plasma were not present. The profile apparently adjusts itself so that the parameters almost allow for full penetration, but fall short of it. This is exactly what is required if the RMF is to penetrate to the null, but not beyond that point. It is relatively easy to understand the mechanism that causes this. With the fixed \(\eta\) used for this calculation, the parameter \(\lambda\) does not vary much. However, the parameter \(\gamma\) varies inversely with density. When the RMF penetrates beyond the field null, flux is generated and the FRC shrinks axially and expands radially. This compresses the flux between the separatrix and the coil increasing the external magnetic field. To maintain pressure balance, the density must also increase, which in turn decreases \(\gamma\). On the other hand, if the RMF does not reach the null, the FRCs poloidal flux decays and the FRC elongates and shrinks radially, which decreases the density. This leads to an increase in \(\gamma\), which continues until the RMF can penetrate to the null. The only balance that can be reached is when the RMF just barely reaches the null so that the flux neither increases nor decreases.

We now examine the process through which the equilibrium is established. Figures 5 and 6 show the magnetic field and current profiles early in the calculation at 20 and 40 \(\mu\)s, respectively. At 20 \(\mu\)s the RMF has penetrated to about 14

Figure 5: Profiles of (a) magnetic field, and (b) normalized azimuthal current and \(\alpha\), for the _standard_ sustainment calculation, at \(t=20\)\(\mu\)s.

Figure 6: Profiles of (a) magnetic field, and (b) normalized azimuthal current and \(\alpha\), for the _standard_ sustainment calculation, at \(t=40\)\(\mu\)s.


<!-- Page 7 -->

cm, and the \(B_{z}\) profile is steep near here due to the strong current drive. A flat region has developed just inside of where the current drive ends. Figure 5(b) shows that \(\alpha\)\(\approx\) 1 in the outer region, where the RMF is strong. Figure 6 shows that by 40 \(\mu\)s, the RMF current drive has penetrated all the way to 11 cm, but the \(\alpha\) profile is significantly below 1 inside of 16 cm. That raises an interesting question: Are the electrons _slipping_ relative to the RMF field, or has the rate of rotation of the RMF field slowed in this region. Figure 7 shows that it is the latter. This figure shows a projection of the magnetic field lines in the \(r\)-\(\theta\) plane, at four times during the penetration phase of the calculation. By 40 \(\mu\)s, the RMF has penetrated too far and cannot provide the torque required to overcome the friction between the electrons and ions, so the electron velocity slows. This causes the field lines to wind up around these slowing electrons leading to an anti-parallel configuration, which quickly decays from resistive effects. This phenomenon is found to repeat itself three times between 20 and 100 \(\mu\)s for this calculation. This figure vividly illustrates how the RMF can be expelled from the FRC when it penetrates too far.

The only simulation parameters that affect the final equilibrium, assuming equilibrium is achieved, are plasma temperature, resistivity \(\eta\), RMF field strength \(B_{\omega}\), RMF frequency \(\omega_{\rm RMF}\), and flux on the coil \(\phi_{c}\). Other parameters that specify the initial conditions determine the path to the equilibrium, but do not affect the equilibrium. The initial conditions must be chosen suitably, however, or equilibrium will not be found. For example, if the initial density is too low to support the \(B_{z}\) field gradient with all the electrons rotating at the RMF angular velocity, the configuration will quickly decay along a path that does not lead to FRC equilibrium.

The fact that our equilibrium \(B_{z}\) profile has a large flat region near the magnetic null means that many of the electrons are rotating slower than \(\omega_{\rm RMF}\). In other words, the plasma line density \(N\) exceeds the critical line density \(N_{c}\), where \(N_{c}\) is defined so that if all of the electrons corotate with the RMF, the resulting \(J_{\theta}\) would exactly support the \(B_{z}\) profile. If we further simplify the definition of \(N_{c}\) by assuming momentarily that the plasma is concentrated in a thin sheath at the null, Faraday's law can be used to derive \(N_{c}\) as

\[N_{c}=\frac{4\,\pi B_{z,\rm Ext}}{\mu_{0}q_{e}\omega}, \tag{20}\]

where we have used the fact that the change in axial field \(\Delta B_{z}\) is approximately \(2B_{z\epsilon}\). With this definition, a new parameter that characterizes the profile shape can be defined as \(\zeta=N_{c}/N\). For a RMF-driven equilibrium, \(\zeta\) must be less than 1 to support the axial field (otherwise some electrons would have to rotate faster than the RMF). For small \(\zeta\), the plasma line density is much larger than that required to support the \(B_{z}\) profile, which implies that there must be a large region with a low \(J_{\theta}\). This corresponds to a \(B_{z}\) profile with a large flat region near the null. On the other hand, as \(\zeta\) approaches 1, there is barely enough plasma to support the \(B_{z}\) field swing with all of the electrons participating, so the flattened region near the null must vanish. For the standard calculation in equilibrium, \(B_{z\epsilon}=0.130\,\mathrm{G}\), \(N_{c}=3.7\times10^{17}\,\mathrm{m}^{-1}\), and \(N=6.9\times10^{17}\,\mathrm{m}^{-1}\), so \(\zeta=0.55\).

The 50 eV temperature limit was chosen to agree with experimental observations with no further justification. However, the solution is quite sensitive to this parameter. Table 1 summarizes the equilibrium parameters for three calculations where the temperature limit \(T_{\rm Limit}\), was varied. The parameter \(\zeta\) varies from 0.3 to 0.7, with \(T_{\rm Limit}\) ranging from 20 to 90 (with the resistivity held constant). The dependence of profile shape on \(\zeta\) is illustrated in Fig. 8, which shows the normalized \(B_{z}\) profiles for these three calculations. As the temperature is increased, the density is forced lower (pressure balance) and the current drive becomes more efficient. This results in an increased \(J_{\theta}\), and consequently a larger \(B_{z}\) swing and a larger \(B_{z\epsilon}\). However, there is a limit--if the parameter \(\zeta\) increases beyond 1, there are not enough electrons to support the field swing and the configuration collapses. In fact, for a calculation like those summarized in Table 1, except with \(T_{\rm Limit}=100\,\mathrm{eV}\), the trapped flux increases for the first 150 \(\mu\)s, but then decays until field reversal is lost by about 350 \(\mu\)s.

It is interesting to note that the RMF must penetrate all the way to the field null if \(E_{\theta}\) is to be zero there, as required for unchanging trapped flux. At first glance, this may seem almost impossible when there is a large flat region near the null, since the RMF can only penetrate a classical skin depth \(\delta\), in the region of low current. This dilemma is resolved by

\begin{table}
\begin{tabular}{c c c c c c} \(T_{\rm Limit}\) & \(\gamma\) & \(\lambda\) & \(\gamma^{\prime}/\gamma_{c}\) & \(B_{z\epsilon}\) & \(x_{s}\) & \(\zeta\) \\ \hline
20 & 37 & 45 & 0.48 & 101 & 0.847 & 0.30 \\
50 & 53 & 47 & 0.66 & 134 & 0.887 & 0.55 \\
90 & 55 & 49 & 0.66 & 175 & 0.91 & 0.70 \\ \end{tabular}
\end{table}
Table 1: Predicted equilibrium parameters when \(T_{\rm Limit}\) is varied.

Figure 7: Magnetic field lines at four times during RMF penetration.


<!-- Page 8 -->

making most of the flat region have a small negative \(B_{z}\), thereby positioning the null near the outer edge of the flat region. This is clearly illustrated in Fig. 8, where a dashed vertical line, drawn at the expected position of the null (\(x_{s}/\sqrt{2}\)), is well inside of the position of the actual null. This positioning of the null requires that the pressure is not a perfect function of \(\psi\), as normally assumed for a FRC. This is allowed in the model. While Eqs. (13)\(-\)(15) continually work to equalize pressure between inner and outer field lines, the steady-state radial flow ensures that the pressure on the inner field lines remains higher.

## IV Discussion and Summary

A full MHD \(r\)-\(\theta\) model has been developed to study RMF current drive in a FRC. This is the first time that plasma motion and profiles have been self-consistently included in RMF current drive calculations. It was found that a steady-state RMF-driven equilibrium is predicted when some 3-D effects are included. Both an axial pressure balance (average \(\beta\)) condition and equalization of pressure between inner and outer field lines are included. Pressure equalization between inner and outer field lines is an essential component of the RMF current drive since it permits a steady flow velocity, which provides the current drive on the inner field lines. Steady-state solutions can be found without employing the average \(\beta\) condition, but only when the initial conditions are carefully specified. Inclusion of the average \(\beta\) condition leads to steady-state solution over a much wider set of initial conditions. However, it is acknowledged that there are many other 3-D effects that could play an important role in RMF current drive. For example, real-world RMF antennas have an azimuthal current flow at their ends, which makes the RMF field structure much more complex in this region, but it is beyond the scope of this model to address this effect.

It is also acknowledged that the average \(\beta\) condition is only approximately satisfied in the presence of a RMF. It was derived assuming that the external magnetic field is parallel to the confining flux observer, which is no longer strictly correct due to the presence of the transverse RMF field. However, since the pressure due to this component is small compared to the poloidal field pressure, the relationship is still a good approximation. In addition, there is some question about the validity of equalizing the pressure and temperature between inner and outer field lines. It has been shown that the transverse field associated with a fully penetrated RMF will open all the field lines.[13] With partial penetration, it is expected that the outer and inner field lines are no longer connected. However, the time-averaged magnetic field (which the ions see) is unaltered from that of a traditional FRC. It is beyond the scope of this paper for us to address these concerns.

A temperature clamp was employed to prevent the temperature from rising significantly above the temperatures observed in experiments. The need for this clamp is not understood at this time. It could be due to high thermal conduction losses not accounted for in the model. Because the RMF field opens the field lines, the electrons actually see a short path to the cold quartz walls. (The ions, on the other hand, do not really see the RMF component of the field since it oscillates many times while they execute one gyro-orbit.) The very low plasma density near the walls (which is also observed in experiments) should limit the loss rate through this channel, but no analysis of the thermal conduction loss has been made. The temperature could also be limited by radiation from impurities, however, the calculations already assume radiation cooling due to a 3% impurity content. This is believed to be in excess of the actual impurity content, since doping with a 0.5% CO\({}_{2}\) was found to roughly double the radiation level and cause a major degradation in the RMF current drive.[14] It has also been suggested[12] that the resistivity employed by the model is higher than it is in the STX experiments, leading to high heating rates. The resistivity used in the model is consistent with estimates from past FRC experiments and gives general agreement with the observed field strength and profile. Perhaps, in experiments, the resistivity is reduced by the stabilizing effect of the RMF--especially in the outer regions where the RMF field strength is larger. Reducing the resistivity in the calculations is found to increase the trapped flux level significantly above the observed levels, but other effects may be reducing the effective current drive in the experiments.

Clearly, it is desirable to apply a full 3-D model to this problem. The main difficulty is associated with the full inclusion of the Hall term, which is the essence of RMF current drive. The semi-implicit numerical algorithm employed in this study is unconditionally stable, but it is found that the time step must be very small for sufficient accuracy. In fact, typical calculations are run for 2\(\times\)10\({}^{6}\) time steps with 100 radial cells. It is expected that a full 3-D calculation would need the same spatial and temporal resolution, making it computationally expensive.

Detailed comparisons between numerical calculations and experimental results await further reduction of the experimental data. Hopefully, some of the outstanding questions discussed above will be illuminated during the course of this work.

Figure 8: Normalized equilibrium \(B_{z}\) profiles for three calculations, with \(T_{\rm{limit}}\)=20, 50, and 90 eV.


<!-- Page 9 -->

###### Acknowledgements.

 The author acknowledges many useful discussions with J. T. Slough and A. L. Hoffman during the course of this work. The present work has been supported by grants from the Office of Fusion Energy Sciences of the U.S. Department of Energy.

## References

* [1] H. A. Blevin and P. C. Thonemann, Nucl. Fusion Suppl. **55**, Pt. 1 (1962).
* [2] W. N. Hugrass, I. R. Jones, K. F. McKenna, M. G. R. Phillips, R. G. Storer, and H. Tuczek, Phys. Rev. Lett. **44**, 1676 (1980).
* [3] A. J. Knight and I. R. Jones, Plasma Phys. Controlled Fusion **32**, 575 (1980).
* [4] J. T. Slough and K. E. Miller, Phys. Plasmas **7**, 1945 (2000).
* [5] I. R. Jones, Phys. Plasmas **6**, 1950 (1999).
* [6] I. R. Jones and W. N. Hugrass, J. Plasma Phys. **26**, 441 (1981).
* [7] W. N. Hugrass and R. C. Grimm, J. Plasma Phys. **26**, 455 (1981).
* [8] M. Ohnishi, A. Ishida, Y. Yamamoto, and K. Yoshikawa, Trans. Fusion Technol. **27**, 391 (1995).
* [9] R. D. Miltroy, Phys. Plasmas **6**, 2771 (1999).
* [10] A. L. Hoffman, "Rotating magnetic field current drive of FRCs subject to equilibrium constraints" to be published in Nucl. Fusion.
* [11] W. T. Armstrong, R. K. Linford, J. Lipson, D. A. Platts, and E. G. Sherwood, Phys. Fluids **24**, 2068 (1981).
* [12] J. T. Slough (private communication, 2000).
* [13] S. A. Cohen and R. D. Miltroy, Phys. Plasmas **7**, 2539 (2000).
* [14] J. T. Slough and K. E. Miller, "Enhanced confinement and stability of a field-reversed configuration with rotating magnetic field current drive," to be published in Phys. Rev. Lett.