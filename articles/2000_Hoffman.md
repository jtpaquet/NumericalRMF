---
paper_id: "2000_Hoffman"
title: "## Appendix A Rotaing magnetic field current drive of FRCs subject to equilibrium constraints"
authors: ""
year: "2000"
institution: ""
publisher: ""
doi: ""
keywords: ""
abstract: "The standard analysis for rotating magnetic field (RMF) current production in simple fixed density plasma columns is extended to include the critical interaction between the RMF drive and an elongated field reversed configuration (FRC) equilibrium inside a flux observer. The standard analysis involves penetration of the RMF into a highly conducting plasma through production of synchronous rotation of the electron fluid due to Hall terms in Ohm's law. In the present article the RMF analysis is combined with two dimensional equilibrium constraints, which have a defining role in governing the RMF penetration. Very simple analytic models are developed, based on instantaneous RMF penetration into an edge layer, that illuminate the basic parameters necessary for current drive and flux buildup or sustainment to be successful. The model does not account for time dependent RMF penetration or transport consistent density profiles, but it clearly points out the basic conditions that must be satisfied, and how they scale, for RMF sustainment of FRCs to be effective."
num_references: "26"
num_citations: ""
source_file: "2000 Hoffman - Rotating magnetic field current drive of FRCs subject to equilibrium constraints.pdf"
---



<!-- Page 1 -->

## Appendix A Rotaing magnetic field current drive of FRCs subject to equilibrium constraints

This article has been downloaded from IOPscience. Please scroll down to see the full text article.

2000 Nucl. Fusion 40 1523

([http://iopscience.iop.org/0029-5515/40/8/310](http://iopscience.iop.org/0029-5515/40/8/310))

View the table of contents for this issue, or go to the journal homepage for more

Download details:

IP Address: 139.184.30.133

The article was downloaded on 12/01/2013 at 23:06

Please note that terms and conditions apply.


<!-- Page 2 -->

# Rotating magnetic field current drive

of FRCs subject to equilibrium constraints

A.L. Hoffman

Redmond Plasma Physics Laboratory, University of Washington,

Washington, DC, United States of America

###### Abstract

The standard analysis for rotating magnetic field (RMF) current production in simple fixed density plasma columns is extended to include the critical interaction between the RMF drive and an elongated field reversed configuration (FRC) equilibrium inside a flux observer. The standard analysis involves penetration of the RMF into a highly conducting plasma through production of synchronous rotation of the electron fluid due to Hall terms in Ohm's law. In the present article the RMF analysis is combined with two dimensional equilibrium constraints, which have a defining role in governing the RMF penetration. Very simple analytic models are developed, based on instantaneous RMF penetration into an edge layer, that illuminate the basic parameters necessary for current drive and flux buildup or sustainment to be successful. The model does not account for time dependent RMF penetration or transport consistent density profiles, but it clearly points out the basic conditions that must be satisfied, and how they scale, for RMF sustainment of FRCs to be effective.

## 1 Introduction

Rotating magnetic fields (RMFs) have been shown to be capable of driving azimuthal currents in a plasma column through a whole series of rotamak experiments [1, 2, 3, 4, 5, 6, 7]. Many numerical and analytical calculations have also been carried out to describe both the penetration of the RMF into a conducting plasma column [8, 9, 10] and the plasma and RMF conditions that are consistent with full penetration [11, 12, 13]. Both types of calculations show approximately the same condition necessary for penetration, namely that the relative (to friction) RMF drive parameter \(\gamma=\omega_{ce}/\nu_{ei}\) (where \(\omega_{ce}=eB_{\omega}/m_{e}\) is the electron cyclotron frequency in a rotating field of strength \(B_{\omega}\) and \(\nu_{ei}\) is the electron-ion collision frequency) exceeds the penetration parameter \(\lambda=a/\delta\) (where \(a\) is the column radius and \(\delta=(2\eta_{\parallel}/\mu_{0}\omega)^{1/2}\) is the classical penetration distance into a column with parallel resistivity \(\eta_{\parallel}\) of an RMF with frequency \(\omega\)). Penetration is made possible by the driven synchronous rotation of the electrons with frequency \(\omega_{r}\) such that the relative frequency seen by the electrons \(\varpi=\omega-\omega_{r}\) is a small quantity. For conditions where full penetration is achieved, it can be shown that \(\omega_{r}\approx\omega/(1+2/\gamma^{2})\) so that for large \(\gamma\), \(\varpi=(2/\gamma^{2})\omega\) and the effective penetration distance \(\delta^{*}=(2\eta_{\parallel}/\mu_{0}\varpi)^{1/2}\) is equal to 0.71\(\gamma\delta\)[12, 13]. Having the length \(\delta^{*}\) exceed \(a\) is the origin of the \(\gamma>\lambda\) requirement. Numerical calculations show the penetration process to be extremely non-linear, with very little penetration and current drive realized until \(\gamma\) exceeds \(\lambda\)[8]. Recent calculations by Milroy show that a value of \(\gamma_{c}=1.12\lambda[1+0.12(\lambda-6.5)^{0.4}]\) is needed for RMF penetration [10].

RMFs have been proposed for use in a method for building up and sustaining the flux in a field reversed configuration (FRC) [13]. Such flux buildup and sustainment have been demonstrated in rotamak experiments, but calculations have never been performed that try to self-consistently combine the two dimensional FRC equilibrium constraints with RMF drive. This is critical since in an FRC the total line current will be set by that needed to reverse the external field. In order for RMF drive to be successful, the current profile must be consistent with the conditions necessary for RMF penetration. The RMF will tend to drive the electron velocity synchronously, and thus somewhat change the current profile, if these changes are not too large. There are thus two requirements for successful RMF current drive of FRCs. The first is the aforementioned \(\gamma>\lambda\) overall force condition, and the second is that the RFM frequency be chosen such that if all the electrons between the field null and the separatrix rotate at this velocity, they will produce sufficient current to approximately reverse the external equilibrium field.

RMF drive in an FRC will result in radial expansion rather than simple current increase. In Australian rotamak experiments at Flinders University this expansion continues until the FRC separatrix encounters the plasma tube wall [6, 7]. In more recent experiments at the University of Washington, using a flux observer, FRC expansion is limited by


<!-- Page 3 -->

external field buildup, which allows the line current to increase [14]. Since the RMF can only penetrate when the electron azimuthal velocity is nearly synchronous with the RMF, successful RMF drive consistent with equilibrium must result in an annular current layer at the edge, whose thickness is determined by the external field, RMF frequency and electron density. In this article a model will be developed to describe this process self-consistently. The model will be an equilibrium one, and thus assume instantaneous radial pressure equilibration and instantaneous RMF penetration (where possible due to \(\gamma\) exceeding \(\lambda\)) to a depth governed by the equilibrium. Since this depth is generally shallow, and radial pressure equilibration is rapid, this model should accurately depict the overall dynamics.

The calculations are based on the analytical RMF penetration solutions obtained by Hugrass [15] for the case of a fixed electron rotational speed \(\omega_{r}\) throughout the plasma column. This simplification is justified in an equilibrium approach since the outer current layer must be near synchronous for the RMF to penetrate. Rather than allow the RMF drive to determine the current, having \(\gamma>\lambda\) will result in a positive azimuthal electric field, \(E_{\theta}\), which will increase the FRC flux. We will calculate this electric field in a simple manner by evaluating overall torques due to RMF drive and resistive friction. It is, of course, necessary that some RMF penetrate to the field null to result in a positive value of \(E_{\theta}\) there, but we will assume that if the penetration distance, using simple analytic expressions for density profiles, is close to the field null, enough flexibility exists in the detailed interaction to allow positive rates of flux increase. This phenomena is helped by the edge drive flattening the axial field profile, resulting in diminished resistive \(\eta j_{\theta}\) at the field null. The above treatment is somewhat crude, but it is the best that can be done without a full transport calculation. It will be seen to be extremely useful in determining the conditions necessary for RMF drive to proceed in FRCs, in deriving basic scaling laws for rates of flux increase and in explaining experimental results.

The calculations of Hugrass are used in Section 2 to derive analytic expressions, relevant to an FRC, for the RMF and resistive forces exerted on the electrons. These forces will be combined with FRC equilibrium constraints and representative density profiles in Section 3 to determine net torques and azimuthal voltages at the field null. The voltages and equilibrium model are used in Section 4 to study the temporal evolution of an FRC under various conditions. The results will be compared with recent experiments. Finally, a summary and criteria for RMF flux sustainment, as related to larger devices or reactors, is given in Section 5.

## 2 Derivation of RMF force and torque expressions

RMF current drive will be calculated for a one dimensional FRC (\(B_{z}=B_{z}(r)\)). Radial growth and a non-constant line density will later be allowed for, which presumes axial shrinkage or growth. Axial length changes are a standard feature of elongated FRC equilibria, but are not relevant to our one dimensional model. Electron inertia will be ignored, which implies that \(\nu_{ei}\gg\omega\). This is appropriate for most of the low electron temperature experiments that have been conducted, but electron inertia may be important at higher temperatures. The actual requirement is that \(\nu_{ei}\gg\varpi\), which will generally apply in the edge layer, so that electron inertia may not be important in any case. Solutions including the effect of electron inertia are given and discussed in the Appendix. Having \(\nu_{ei}\gg\varpi\) ensures that the electrons oscillate resistively in the RMF generated axial electric field. Their axial oscillation is thus in-phase with the RMF generated \(E_{z}\) and \(B_{r}\), resulting in a net \(\langle\nu_{ez}B_{r}\rangle\) azimuthal drive term. Ions have large inertia (\(\omega\gg\nu_{ie}\)) and will oscillate inductively (out of phase with \(B_{r}\)), and thus not be subject to the azimuthal drive. This has all been discussed extensively in previous work[12, 13].

Ignoring inertia, the electron equation of motion can be written as

\[\mathbf{E}+\mathbf{v}_{e}\times\mathbf{B}=\eta\mathbf{j}. \tag{1}\]

We will use the general definition \(\eta=m_{e}\nu/ne^{2}\). For the initial RMF penetration calculations we need only have a uniform value of parallel resistivity \(\eta_{\parallel}\), but later in calculating total RMF drive and resistive forces we will utilize appropriate density profiles. We will also assume that all the current is carried by electrons, but later show how the results are affected by some specified fraction of the azimuthal current being carried by ions. The gradient in electron pressure has also been ignored in this form of the generalized Ohm's law since we will only deal with the axial and azimuthal equations, and \(p_{e}\) is assumed to be only a function of \(r\). It is important to note that this equation implicitly includes the Hall term, which is critical to RMF current drive. When a


<!-- Page 4 -->

rotating magnetic field is applied with the vacuum form

\[B_{r}=B_{\omega}\cos(\omega t-\theta),B_{\theta}=B_{\omega}\sin(\omega t-\theta) \tag{2}\]

it has been shown that to first order there will be steady azimuthal and radial velocities and a steady azimuthal electric field \(E_{\theta}(r)\)[12]. The average (over an RMF cycle) \(\theta\) component of Eq. (1) can then be written as

\[E_{\theta}=\eta_{\perp}j_{\theta}+v_{er}B_{z}-\langle v_{ez}B_{r}\rangle. \tag{3}\]

Ignoring ion axial motion, the electron axial oscillation can be calculated from the axial component of Eq. (1),

\[E_{z}=\eta_{\parallel}j_{z}+v_{e\theta}B_{r}. \tag{4}\]

Cross-field radial flow can contribute to current drive, but is absent at the FRC field null (\(B_{z}=0\)) and thus cannot alone sustain the FRC flux.

Experimentally the RMF only penetrates at most to the field null so that current drive on the inner field lines can only result from an inward radial velocity there. The RMF on the outer field lines can reverse the normally outward diffusion there, resulting in an overall inward flow which can be accommodated in two dimensions by an axial flow from the inner to the outer field lines along the same flux surface. Since in this work we are concerned only with total torques and average values of \(E_{\theta}\), we will ignore this radial flow (setting \(v_{r}=0\)) and take into account the resistive torques from both the inner and outer currents in balances against the RMF produced torque. This will account, in an approximate manner, for the need for the RMF to be strong enough to drive the entire FRC current. It is assumed that in equilibrium a flow pattern will be set up so that \(E_{\theta}\) in Eq. (3) can be identically zero everywhere. Numerical calculations showing this process will be reported in a subsequent paper by Milroy [16].

Different parallel and perpendicular resistivities have been allowed to show the individual effects of each component. Equation (4) is non-linear since both \(v_{e\theta}\) and \(B_{r}\) are normally functions of radius. It has been solved numerically by mode analysis in Refs [8, 9, 10]. However, useful analytic solutions can be obtained for a specified \(v_{e\theta}=\omega_{r}r\), which makes the equation linear.

Following the method first used by Hugrass, Grimm and Jones [8, 11] we will utilize the vector potential \(\mathbf{A}\) (\(\mathbf{B}=\mathbf{\nabla}\times\mathbf{A}\)) for the radial and azimuthal components of the RMF. Thus \(\mathbf{A}=A_{z}\mathbf{e}_{z}\) and

\[B_{r}=\frac{1}{r}\frac{\partial A_{z}}{\partial\theta},\ B_{\theta}=-\frac{ \partial A_{z}}{\partial r} \tag{5}\]

Also, \(\mu_{0}\mathbf{j}=\mathbf{\nabla}\times\mathbf{B}=-\nabla^{2}\mathbf{A}\) (choosing \(\mathbf{A}\) such that \(\mathbf{\nabla}\cdot\mathbf{A}=0\)) so that

\[j_{z}=-\frac{1}{\mu_{0}}\nabla^{2}A_{z}. \tag{6}\]

Since \(\mathbf{\nabla}\times\mathbf{E}=-\partial\mathbf{B}/\partial t=-\partial(\mathbf{\nabla} \times\mathbf{A})\partial t\), \(\mathbf{E}=-\partial\mathbf{A}/\partial t\). Taking the time derivative of Eq. (4) we thus obtain

\[\frac{\partial E_{z}}{\partial t}+\omega_{r}\frac{\partial E_{z}}{\partial \theta}=\frac{\eta_{\parallel}}{\mu_{0}}\nabla^{2}E_{z}. \tag{7}\]

As might be expected, this has the form of a diffusion equation. We will solve it for a'steady state', after the RMF diffuses into the column, with \(E_{z}\) specified to have the form

\[E_{z}=\hat{E}_{z}(r)e^{i(\omega t-\theta)}. \tag{8}\]

Full multimode calculations show that this single mode is very dominant when an external field as given by Eq. (2) is specified [10]. \(\hat{E}_{z}\) would initially also be a function of \(t\), but we will only solve for the steady state result. Diffusion to the FRC field null will be rapid when the initial azimuthal currents are close to the final ones. \(\hat{E}_{z}\) will be complex in this notation, but henceforth we will drop the 'over-hat' notation and recognize that all oscillating quantities are complex and must be multiplied by \(e^{i(\omega t-\theta)}\). A simple result is then that \(E_{z}=-i\omega A_{z}\).

Equation (7) can now be written as

\[i(\omega-\omega_{r})E_{z}=\frac{\eta_{\parallel}}{\mu_{0}}\nabla^{2}E_{z}. \tag{9}\]

Using our previous definitions of \(\varpi\) and \(\delta^{*}\), calling \(k=\sqrt{2}/\delta^{*}\) and writing out the Laplacian operator in cylindrical co-ordinates, Eq. (9) takes the form of a Bessel equation,

\[r^{2}\frac{\partial^{2}E_{z}}{\partial r^{2}}+r\frac{\partial E_{z}}{\partial r }-(1+ik^{2}r^{2})E_{z}=0 \tag{10}\]

Since one of the coefficients is complex, this has the unusual solution

\[E_{z}(r)=\frac{2\omega}{\sqrt{ik}}\ B_{\omega}\frac{I_{1}(\sqrt{ik}r)}{I_{0}( \sqrt{ika})} \tag{11}\]

where the coefficient in front has been chosen to match to the external vacuum solution \(E_{z}=\omega rB_{\omega}-C/r\), where \(C\) is a constant and both \(E_{z}\) and \(\partial E_{z}/\partial r\)


<!-- Page 5 -->

must be continuous at \(r=a\). \(B_{r}\) and \(B_{\theta}\) can then be obtained from Eq. (5) making use of \(A_{z}=iE_{z}/\omega\):

\[B_{r}=\frac{2}{\sqrt{ikr}}\frac{I_{1}(\sqrt{ik}r)}{I_{0}(\sqrt{ik}a)}B_{\omega} \tag{12}\]

\[B_{\theta}=i\left(B_{r}-2\frac{I_{1}(\sqrt{i}kr)}{I_{0}(\sqrt{ik}a)}B_{\omega} \right). \tag{13}\]

These expressions, containing both amplitude and phase information, are useful in interpreting internal probe measurements of the RMF penetration profiles. The modified Bessel functions of the argument \(\sqrt{i}x\) can be written in terms of the her functions, \(I_{0}(\sqrt{i}x)=\mbox{ber}_{0}(x)+i\mbox{bei}_{0}(x)\) and \(iI_{1}(\sqrt{i}x)=\mbox{ber}_{1}(x)+i\mbox{bei}_{1}(x)\). However, the Bessel solutions are useful without such resolution into real and imaginary parts when taking limits for small and large \(kr\). These limits will be given later after expressions for forces and total torques have been derived.

The oscillating axial current is calculated from Eq. (4) as

\[j_{z}=\frac{E_{z}-\omega_{r}rB_{r}}{\eta_{\parallel}}=\frac{\varpi}{\omega} \frac{E_{z}}{\eta_{\parallel}}. \tag{14}\]

The \(j_{z}\) screening currents are reduced by near synchronous rotation, which is the factor that allows the RMF to penetrate a distance \(\delta^{*}=(\omega/\varpi)^{1/2}\delta\). The average \(\theta\) force \(F_{\theta}=\langle j_{z}B_{r}\rangle\) is now simply given by

\[F_{\theta}=\frac{\varpi}{\omega}\frac{\langle E_{z}^{2}\rangle}{\omega r\eta_ {\parallel}}=\varpi r\frac{\langle B_{r}^{2}\rangle}{\eta_{\parallel}} \tag{15}\]

since \(j_{z}\), \(E_{z}\) and \(B_{r}\) are all in-phase. Thus \(\langle E_{z}^{2}\rangle\) is simply \(E_{z}^{2}/2\). Using the magnitude of the expression for \(E_{z}\),

\[F_{\theta}=\frac{2B_{\omega}^{2}}{\mu_{0}r}\left|\frac{I_{1}(\sqrt{i}kr)}{I_{ 0}(\sqrt{i}ka)}\right|^{2}. \tag{16}\]

Using the her expansions for \(I_{0}\) and \(I_{1}\), the force can be integrated to find the total torque exerted on a plasma column by the RMF. The result is

\[T_{M}=\int_{0}^{a}2\pi r^{2}F_{\theta}\,dr=\frac{2\pi B_{\omega}^{2}}{\mu_{0} }a\delta^{*}f(ka) \tag{17}\]

where

\[f(ka) = \sqrt{2}\] \[\times\,\left(\frac{\mbox{ber}_{1}(ka)\mbox{bei}_{1}^{\prime}(ka) -\mbox{ber}_{1}^{\prime}(ka)\mbox{bei}_{1}(ka)}{\mbox{bei}_{0}^{2}(ka)+\mbox{bei }_{0}^{2}(ka)}\right).\]

This is the same expression as that derived by Hugrass in Ref. [15] except for the small misprint of leaving out a factor \(ka\) in his work. For large \(ka\) (small \(\delta^{*}\) and little penetration) \(f(ka)\) is approximately equal to unity (\(f(ka)\approx 1-0.7/ka\) for \(ka>3\)) and the total torque is proportional to the penetration distance, as might be expected. For small \(ka<1\), \(f(ka)\approx(\sqrt{2}/16)(ka)^{3}\). In this full penetration case the torque is proportional to \((1/\delta^{*})^{2}\) or to \(\varpi\), as also might be expected from a consideration of induction motors. No slip means no torque.

It is useful to write the torque in terms of a reference value

\[T_{0}=\frac{2\pi a^{2}B_{\omega}^{2}}{\mu_{0}} \tag{19}\]

\[T_{M}=T_{0}\left(\frac{\sqrt{2}f(ka)}{ka}\right) \tag{20}\]

so that a plot of \(\sqrt{2}f(ka)/ka\) versus \(ka\) shows at what value of \(ka\) and thus \(\varpi/\omega=(ka/\lambda)^{2}/2\) the torque peaks. This is done in Fig. 1. The torque is seen to peak at \(ka\approx 2.5\) or \(a/\delta^{*}\approx 2\). At this value \(\varpi/\omega=4/\lambda^{2}\), so that very close to synchronous rotation is required for large \(\lambda\). The torque decreases to about one third its peak value at \(ka=1\) (\(a/\delta^{*}=0.71\) and \(\varpi/\omega=0.5/\lambda^{2}\)) and at \(ka=10\) (\(a/\delta^{*}=7.1\) and \(\varpi/\omega=50/\lambda^{2}\)). For \(\lambda>50\) this is not a very wide range. Plots of the radial component of the rotating field, \(B_{r}\), and the azimuthal component, \(B_{\theta}\), are shown on Figs 2(a) and 2(b) for the above and other cases. For small \(a/\delta^{*}\) there is approximately full penetration and \(B_{r}\) and \(B_{\theta}\) have uniform profiles. The phases can also be calculated from Eqs (12) and (13), and they are also the same as in vacuum. For large \(a/\delta^{*}\) the RMF only penetrates approximately a distance \(\delta^{*}\) into the plasma column and the current drive force is localized there. It is for this reason that the full rigid rotor solutions are applicable to an edge current layer configuration where penetration does not occur beyond a certain depth (the current layer thickness). \(B_{\theta}\) is about a factor of two higher at the edge due to the axial image currents which hinder RMF penetration. \(B_{\theta}\) does not contribute to current drive, but it does add to the radial confinement force provided by the axial field. This effect has been seen in the University of Washington STX experiments [14] where the internal \(B_{z}\) field becomes higher than the external \(B_{z}\) field, and RMF radial force gradients could be an important stabilization factor for FRCs. It will be seen in the next section that effective RMF


<!-- Page 6 -->

current drive of FRCs always results in these strong edge gradients.

It is useful to look at the calculated fields and forces under the two limiting assumptions of small and large \(ka\) (small and large \(a/\delta^{*}\)). For small \(kr\), \(I_{0}(\sqrt{ik}r)=1+(\sqrt{i}kr/2)^{2}+...\approx 1\) and \(I_{1}(\sqrt{i}kr)=I_{0}^{\prime}(\sqrt{i}kr)\approx\sqrt{i}kr/2\). This yields \(E_{z}\approx\omega rB_{r}\), \(B_{r}\approx B_{\omega}\) and \(B_{\theta}\approx-iB_{\omega}\). In complex notation these are just the vacuum fields given by Eq. (2). The force \(F_{\theta}\approx(B_{\omega}^{2}/2\mu_{0})k^{2}r=(rB_{\omega}^{2}/2\eta_{0})\varpi\) and the torque \(T_{M}\approx(\pi B_{\omega}^{2}a^{4}/4\eta_{\parallel})\varpi=(\lambda^{2}/4) T_{0}(\varpi/\omega)\). If we calculate a torque due to electron-ion friction as

\[T_{\eta}=\int_{0}^{a}2\pi^{2}\nu_{\perp}nm_{e}\omega_{r}\,dr\] \[\qquad=\frac{\pi}{2}m_{e}n\nu_{\perp}a^{4}\omega_{r}=\left(\frac{ \lambda^{2}}{2\gamma^{2}}\right)T_{0}\frac{\omega_{r}}{\omega} \tag{21}\]

we can also write this in terms of \(T_{0}\) plus our \(\lambda\) and \(\gamma\) parameters. In that way we can see how the resistive and RMF torques scale relative to one another. We have changed to a definition of

\[\gamma=\frac{\omega_{ce}}{\sqrt{\nu_{\parallel}|\nu_{\perp}}} \tag{22}\]

to allow for different parallel and perpendicular resistivities. Equating \(T_{M}\) to \(T_{\eta}\) yields the relationship \(\varpi=2\omega_{r}/\gamma^{2}\), which is the same result derived from full penetration analysis assuming the necessary condition that \(\omega_{r}\approx\omega\) in this small-\(ka\) limit. This limit only applies when \(\gamma>2\lambda\), so that it does not violate the basic requirement that \(\gamma\) exceed \(\lambda\).

In the large-\(ka\) limit we can use the expansions \(I_{0}(x)\), \(I_{1}(x)\approx e^{x}/(2\pi x)^{1/2}\). This applies perfectly well for complex arguments, resulting in \(I_{1}(\sqrt{i}kr)/I_{0}(\sqrt{i}ka)\approx(a/r)^{1/2}e^{-\sqrt{i}k(a-r)}\). Using the fact that \(\sqrt{i}=(1+i)/\sqrt{2}\), the exponential can be written as \(e^{-i(a-r)/\delta^{*}}e^{-(a-r)/\delta^{*}}\). The following expressions can then be easily calculated:

\[E_{z}\approx \left\{\exp\left[-i\left(\frac{\pi}{4}+\frac{a-r}{\delta^{*}} \right)\right]\right\}\] \[\times\omega\sqrt{\frac{2a}{r}}\delta^{*}B_{\omega}\exp\left[- \left(\frac{a-r}{\delta^{*}}\right)\right]\] \[B_{r}=\frac{1}{\omega r}E_{z} \tag{23}\] \[B_{\theta}\approx \left\{\exp\left[-i\left(\frac{\pi}{2}+\frac{a-r}{\delta^{*}} \right)\right]\right\}\] \[\times 2\sqrt{\frac{a}{r}}B_{\omega}\exp\left[-\left(\frac{a-r}{ \delta^{*}}\right)\right]\] (24) \[F_{\theta}\approx \frac{2B_{\omega}^{2}}{\mu_{0}a}\left(\frac{a}{r}\right)^{2}\exp \left[-2\left(\frac{a-r}{\delta^{*}}\right)\right]\] (25) \[T_{M}\approx\left(1-e^{-2a/\delta^{*}}\right)\frac{2\pi B_{ \omega}^{2}}{\mu_{0}}\,a\delta^{*}. \tag{26}\]

Figure 1: Torque function dependence on \(ka\).

Figure 2: Calculated (a) \(B_{r}\) and (b) \(B_{\theta}\) profiles for various values of \(ka\).


<!-- Page 7 -->

### Hoffman

The quantities inside the curly brackets are just the phases. In this limit the torque is proportional to \(\delta^{*}\) as mentioned earlier. The force at the edge is just \(2B_{\omega}^{2}/\mu_{0}a\), independent of the plasma resistivity. The reason for this is that as the resistivity increases, the axial electron oscillation is more restricted, but the RMF penetration increases, increasing \(E_{z}\) and \(B_{r}\) and thus the driving force (in Eq. (15), \(B_{r}\propto\delta^{*}\propto\sqrt{\eta_{||}}\)). The radial force \(F_{r}=-\langle j_{z}B_{\theta}\rangle\) can also be calculated from these expressions.

If we equate the RMF drive force at the edge, \(r=a\), to the electron friction force \(F_{\eta}=\nu_{\perp}nm_{e}\omega_{r}a\), the condition that \(F_{\theta}\) exceed \(F_{\eta}\) can be written as \(\gamma^{2}\omega>\lambda^{2}\omega_{r}\). Thus, a necessary condition for the RMF to spin the edge electrons up to near synchronous velocity for large-\(\lambda\) columns is the same \(\gamma>\lambda\) condition derived for full penetration. This is the reason why numerical solutions of the full non-linear diffusion equation do not show any significant RFM penetration or current drive until this condition is met.

## 3 FRC equilibrium and flux drive model

We wish to compare the RMF produced torque \(T_{M}\) with the torque \(T_{\eta}\) due to resistivity (given by Eq. (21)), and then use this to estimate effective values of \(E_{\theta}\) at the field null. Intuitively, \(E_{\theta}\) will be positive when \(T_{M}\) exceeds \(T_{\eta}\) and vice versa. It is instructive to plot \(T_{M}/T_{0}\) as a function of \(\omega_{r}\) rather than \(ka\) as has been done in Fig. 1. Since \(ka=\lambda(2\varpi/\omega)^{1/2}\), the plots will look very different for different values of \(\lambda\). Such plots are shown in Fig. 3 for \(\lambda=5\) and \(\lambda=50\). The frictional torque, as given by Eq. (21), simply varies linearly as \(\omega_{r}/\omega\), with the relative magnitude scaling as \((\lambda/\gamma)^{2}\). We can then plot \(T_{\eta}/T_{0}\) on the same graph, realizing that its slope will depend on the value of \(\lambda/\gamma\). This has been done in Fig. 3 for a value of \(\gamma/\lambda=\sqrt{2}\). It can be seen that a \(T_{M}=T_{\eta}\) solution will always be possible for \(\gamma>1.15\lambda\) (\(T_{M}({\rm max})=0.38T_{0}\) and \(T_{\eta}\approx(\lambda^{2}/2\gamma^{2})T_{0}\) for \(\omega_{r}\approx\omega\)). This is very close to Milroy's numerical calculations of \(\gamma_{c}=1.12\lambda[1+0.12(\lambda-6.5)^{0.4}]\). The second term in his bracketted solution reflects the higher initial \(\gamma\) required for initial RMF penetration when \(\lambda\) exceeds 6.5, whereas the \(1.12\lambda\) reflects the minimum value, once the RMF has penetrated, to keep it from being expelled [10]. The analytical model used here is based only on the final equilibrium penetrated condition. The reason for the good applicability of the simple fixed \(\omega_{r}\) calculation is that the true non-linear numerical solution has very similar field profiles once the RMF has penetrated far enough to yield large torques.

It can be seen that for large \(\lambda\) and \(\gamma/\lambda\) not too far above unity, three solutions of \(T_{M}=T_{\eta}\) are possible. (The small (\(\omega_{r}/\omega\) crossing is not shown in Fig. 3.) The same feature is seen in the full numerical solutions. In the simple fixed column numerical solutions the intermediate-\(\omega_{r}\) solution is not a stable solution since a small increase in \(\omega_{r}\) causes a positive value of \(T_{M}-T_{\eta}\), resulting in a continued increase of \(\omega_{r}\). Likewise, a decrease in \(\omega_{r}\) would result in a negative value of \(T_{M}-T_{\eta}\) and a further decrease in electron rotation. However, for an FRC constrained by radial equilibrium requirements, the situation will be different. In that case a positive value of \(T_{M}-T_{\eta}\) will result in an increase in flux, which causes the FRC radius to increase and the average value of \(\omega_{r}\) to decrease (assuming \(n\,{\buildrel\sim\over{\sim}}\,B_{e}\)). The opposite holds true when \(T_{M}-T_{\eta}\) becomes negative. The intermediate solution of Fig. 3 thus becomes the stable solution. (Of course, as seen in Fig. 3 for large \(\lambda\) there is very little difference in \(\omega_{r}\) between the two large \(\omega_{r}/\omega\) solutions.)

We could simply use the rigid rotor solutions of Eq. (20) to calculate net torques in our equilibrium model. However, it is obvious that this will not duplicate the full numerical solutions, nor what is observed experimentally. The curves in Fig. 3 show that, for intermediate values of \(\omega_{r}/\omega\) and modest values of \(\gamma/\lambda\), the resistive torque will exceed the RMF torque even with \(\gamma\) exceeding Milroy's \(\gamma_{c}\). In the numerical solutions the RMF 'eats' its way inwards by sequentially spinning up the outer layers of electrons. An accurate analytic model must reflect this

Figure 3: Normalized RMF and frictional torque dependences on \(\omega_{r}/\omega\).


<!-- Page 8 -->

factor. The same phenomenon must occur experimentally or an FRC could not be built up from a simple plasma column unless \(\gamma\) was exceedingly large compared with \(\lambda\). In our model we will thus assume that if the \(\gamma>\lambda\) criterion is met, a nearly synchronous edge current will be built up until it is sufficient to produce an internal field of opposite sign and equal magnitude to the external field. At that point the FRC must expand rather than produce more current if \(T_{M}\) exceeds \(T_{\eta}\). Before making this adjustment, however, we will introduce a reasonable expression for FRC flux drive.

Choosing a realistic way of relating the calculated torques to a representative value of \(E_{\theta}\) at the field null of an FRC is the most problematic feature of this simple model. In reality, an edge current layer will produce a reduced field gradient and current at the field null and reduce the RMF drive requirements there. As mentioned earlier, the \(E_{\theta}(r)\) profile will be strongly influenced by radial flow, but modelling this accurately requires a full MHD code, which will soon be reported by Milroy [16]. For the present simple analytic model we will make the reasonable assumption that the FRC will adjust itself in such a manner as to produce positive \(E_{\theta}\) at the field null and increasing flux if the RMF produced torque exceeds the resistive torque. Ignoring radial flow, at any given location the azimuthal electric field in Eq. (3) can be written as

\[E_{\theta}=\frac{1}{ne}F_{\theta}-\frac{m_{e}}{e}\,\nu_{\perp}(v_{e\theta}-v_{ i\theta}). \tag{27}\]

At this point we will allow for an ion azimuthal velocity \(v_{i\theta}=\alpha v_{e\theta}\) with \(\alpha\) negative if the ions move in their diamagnetic direction. We will not address the long term question of ion spin-up in the electron diamagnetic direction, which has been discussed elsewhere [17]. Consistent with a rigid rotor rotation we will stipulate an azimuthal field distribution \(E_{\theta}(r)=E_{\theta R}(r/R)\), where \(R=r_{s}/\sqrt{2}\) is the location of the FRC field null for an FRC with separatrix radius \(r_{s}\approx a\). This distribution nearly duplicates the flux decay rates calculated for a non-driven rigid rotor FRC. Multiplying by \(r\) and integrating over the column area gives a value for \(E_{\theta R}\) in terms of the torques,

\[\int_{0}^{a}2\pi r^{2}E_{\theta}dr=\frac{\pi a^{3}}{\sqrt{2}}\,E_{\theta R}= \frac{1}{ne}\,(T_{M}-T_{\eta}). \tag{28}\]

Then the rate of change of FRC flux can be expressed as

\[\left.\frac{d\phi}{dt}\right|_{R}=2\pi RE_{\theta R}=\frac{2}{nea^{2}}(T_{M}- T_{\eta}). \tag{29}\]

Now we must calculate the magnitude of the edge current and find out how it affects the RMF drive and frictional torques. The current is assumed to flow almost synchronously in an edge layer of thickness \(a-y\). For simplicity we choose a parabolic density distribution \(n=n_{0}(a^{2}-r^{2})/(a^{2}-y^{2})\) in this layer. This distribution varies only slightly from the constant temperature rigid rotor density profile \(n=n_{0}\mbox{sech}[K(r/y)^{2}-1]\), where \(K=n_{0}e\omega y^{2}/(2B_{e}/\mu_{0})\approx y/(a-y)\) must be large for the density near the separatrix to be low. This has the advantage that both the density and temperature go to zero at \(r=a\). This parabolic profile is obviously artificial and not based on any transport analysis, but the chosen profile is only meant to be representative, and we are already assuming that the profile can adjust to provide consistent force balance with at least some RMF penetration to the field null if y is close to \(r_{s}/\sqrt{2}\).

The FRC is assumed to have a 2-D equilibrium, so there must also be current flowing on the corresponding inner field lines. (The axial magnetic field must be a function of \(u=1-(r/R)^{2}\).) In addition, the FRC must be allowed to adjust its length to diameter ratio to satisfy the average beta condition \(\langle\beta\rangle=1-\frac{1}{2}x_{s}^{2}\)[18]. \(x_{s}\equiv r_{s}/r_{c}\), \(r_{s}\) is the separatrix radius and \(r_{c}\) is the radius of a flux conserver. This requires an additional parameter in the FRC density profile beyond the radial position \(y\). We will accomplish this by letting the plasma radius \(a\) be less than the separatrix radius \(r_{s}\), so that the two profile parameters are \(y/a\) and \(a/r_{s}\). A sketch of the assumed density and axial magnetic field profiles is shown in Fig. 4. An FRC with low pressure near the separatrix must have such a thin current layer profile in order to satisfy the average beta condition.

All the electrons between \(y\) and \(a\) are assumed to be rotating synchronously with velocity \(v_{e\theta}\approx\omega r\), and there is no current between the field null \(R=r_{s}/\sqrt{2}\) and \(y\). The total line current between the field null and the outer edge is given by

\[I^{\prime}=\frac{1}{4}\,n_{0}e\,(1-\alpha)\,\omega_{r}(a^{2}-y^{2})=\frac{B_{ e}}{\mu_{0}}. \tag{30}\]

\(y\) must be greater than \(R\) for this model to apply, which means that the FRC density must be high enough to reverse the external field, \(B_{e}\), if there was full synchronous rotation. (\((a-y)\) is actually determined by the average beta condition, and must be small when \(\langle\beta\rangle\) is large, requiring an even larger \(n_{0}\).) Otherwise \(\omega_{r}\) would be larger than \(\omega\) and flux drive would be negative. Equation (30) could be used to determine \(y\) if the other quantities were known.


<!-- Page 9 -->

### Hoffman

Under the thin flux layer assumption, Eq. (25) will be used to calculate the total torque by stipulating that \(\delta^{*}=a-y\). Then

\[T_{M}=\int_{y}^{a}2\pi r^{2}F_{\theta}dr=0.865\left(1-\frac{y}{a}\right)T_{0}. \tag{31}\]

This has a maximum value of \(0.253T_{0}\) at \(y=R\approx a/\sqrt{2}\). Using Eq. (20) for the pure rigid rotor solution and Fig. 1 for \(\sqrt{2}f(ka)/ka\), this value is seen to occur for \(ka=4.8\), or for \(\delta^{*}=0.29a\), which is the same as \(a-R\) when \(a\approx r_{s}\). We thus make a continuous transition as \(\delta^{*}\) becomes equal to the null to separatrix distance, using Eq. (31) when \(ka\) exceeds \(4.8\) and Eq. (20) for \(ka\) less then \(4.8\). For realistic cases with large \(\lambda\), \(a-y\) will never be as large as \(0.29a\), and Eq. (31) will always apply.

One may ask why the torque should achieve this maximal value, which requires an optimum value of slip \(\varpi/\omega\). Fixed column numerical calculations show that when \(\gamma\) well exceeds \(\lambda\), RMF penetration will reach all the way to the centre of a column. This would produce too much current, so that the FRC will expand, the density will increase and the average ratio of \(\gamma/\lambda\) will decrease until expansion is halted. At that point the RMF and frictional torques will be in balance. An increase in \(\omega_{r}\) beyond that needed to produce the equilibrium \(a-y\) penetration would result in lowered torque and restoration of the original \(\omega_{r}\). A basic result of the full numerical calculations (and our simple edge force analysis) is that for \(\gamma\) exceeding \(\lambda\) the RMF will always 'eat' its way in as long as the diamagnetic electron current can adjust, consistent with the overall equilibrium, to be nearly synchronous. A basic assumption of this model is that the electron velocity adjusts in a stable manner to produce the torque given by Eq. (31).

Similar adjustments must be made for the resistive torque to account for the regions where the current flows (this is how we account for the fact that the RMF torque on the outer field lines must be sufficient to also drive current on the inner field lines through the aforementioned radial flow pattern). Recalculating the integral in Eq. (21), and stipulating that current on the outer flux surfaces is mirrored by current on the inner surfaces (at the same value of \(|u|\) with \(u=(r/R)^{2}-1\), the resistive torque becomes

\[T_{\eta}=\left(\frac{\lambda^{2}}{2\gamma^{2}}\right)\left[1+\left(\frac{r_{ s}}{a}\right)^{4}-2\left(\frac{y}{a}\frac{r_{s}}{a}\right)^{2}\right](1- \alpha)\frac{\omega_{r}}{\omega}\,T_{0}. \tag{32}\]

In this notation the collision frequencies, and thus \(\gamma\), are based on the average value of density \(\bar{n}=n_{0}/2\). We have also accounted for the ion component of azimuthal current. It should be noted that having the RMF drive only take place in the edge layer results in the RMF torque \(T_{M}\) decreasing more rapidly with decreases in \(a-y\) than the resistive torque \(T_{M}\), so that an even higher ratio of \(\gamma/\lambda\) will be required.

The calculated resistive torque depends explicitly on the chosen density profile, but is certainly sufficiently accurate within the framework of this model. The RMF drive torque depends on the density profile insomuch as the current layer thickness \(a-y\) is affected. The true density profile between \(y\) and \(r_{s}\) can only be calculated with a transport model but it will always have the feature of low density near the separatrix, as sketched in Fig. 4. This is expected due to the strong RMF drive there, and is seen in MHD calculations now being performed by Milroy [16]. The exact details will not have a large effect on the computed torques. For the purposes of this article, the general features are more important.

Now that a method exists for calculating the rate of FRC flux change, it is necessary to calculate a new FRC equilibrium. The equilibrium will be based on the profile given in Fig. 4. We need to calculate \(r_{s}\), \(a\) and \(y\) for a given flux conserver radius \(r_{c}\). The three determining factors will be the ratio of FRC flux to external flux (\(\pi r_{c}^{2}B_{0}\)), the width of the edge layer necessary to carry the reversal current as given by Eq. (30) and the FRC average beta condition, \(\langle\beta\rangle=1-x_{s}^{2}\). Using the representative density distribution \(n=n_{0}(a^{2}-r^{2})/(a^{2}-y^{2})\), with electrons rotating

Figure 4: FRC equilibrium profiles used in the model considered.


<!-- Page 10 -->

at frequency \(\omega_{r}\approx\omega\), the FRC magnetic field profile can be calculated as

\[B_{z}=\left[1-\left(\frac{a^{2}-r^{2}}{a^{2}-y^{2}}\right)^{2}\right]B_{e} \tag{33}\]

with

\[B_{e}=\tfrac{1}{4}\,\mu_{0}n_{0}e(1-\alpha)\omega(a^{2}-y^{2}). \tag{34}\]

\(B_{z}\) is zero between \(R\) and \(y\) and is \(B_{e}\) between \(a\) and \(r_{s}\) as sketched in Fig. 4. The FRC flux can be obtained by integration,

\[\phi=\int_{y}^{r_{s}}2\pi rB_{z}\,dr=\pi\left[(r_{s}^{2}-a^{2})+\frac{2}{3}(a^{ 2}-y^{2})\right]B_{e}. \tag{35}\]

The average beta can also be calculated as

\[\langle\beta\rangle =\frac{1}{0.5\pi r_{s}^{2}}\left[\pi(y^{2}-0.5r_{s}^{2})+\int_{y} ^{a}\left(1-\frac{B_{z}^{2}}{B_{e}^{2}}\right)2\pi r\,dr\right]\] \[=\frac{a^{2}}{r_{s}^{2}}\left[\left(1+\frac{y^{2}}{a^{2}}\right) -\frac{1}{15}\left(1-\frac{y^{2}}{a^{2}}\right)\right]-1. \tag{36}\]

Both these expressions give sharp boundary profile results if \(y=a\)[18].

By setting \(\langle\beta\rangle=1-\tfrac{1}{2}(r_{s}/r_{c})^{2}\), Eqs (34), (35) and (36) can be solved simultaneously to yield \(r_{s}\), \(y\) and \(a\). It is also necessary to use the flux conservation expression \(B_{e}=B_{0}/(1-x_{s}^{2})\) with the vacuum flux, represented by \(B_{0}\), specified. Defining a characteristic frequency \(\omega_{f}\) as

\[\omega_{f}=\frac{2B_{0}}{(1-\alpha)\bar{n}e\mu_{0}r_{c}^{2}} \tag{37}\]

and the normalized flux as

\[\bar{\phi}=\frac{\phi}{\pi r_{c}^{2}B_{0}} \tag{38}\]

a cubic equation for \(x_{s}^{2}\) can be derived,

\[x_{s}^{4}(1-x_{s}^{2})=4(1-x_{s}^{2})^{2}\,\bar{\phi}-\frac{8\omega_{f}}{15 \omega}. \tag{39}\]

Again, this reduces to the sharp boundary result for \(n\to\infty\), which is equivalent to all the azimuthal current flowing on an infinitely thin layer. Once \(x_{s}\) is found, the remaining profile parameters are easily derived,

\[\left(\frac{a}{r_{s}}\right)^{2}=1-\frac{1}{4}\,x_{s}^{2}+\frac{8}{15x_{s}^{2 }(1-x_{s}^{2})}\frac{\omega_{f}}{\omega} \tag{40}\]

\[\left(\frac{y}{a}\right)^{2}=1-\frac{(r_{s}/a)^{2}}{x_{s}^{2}(1-x_{s}^{2})} \frac{\omega_{f}}{\omega}. \tag{41}\]

### Rotating magnetic field current drive of FRCs

Since this model only applies when \(y>r_{s}/\sqrt{2}\), it is necessary that \(\omega\) be greater than \(\omega_{f}/\{x_{s}^{2}(1-x_{s}^{2})[(a/r_{s})^{2}-0.5]\}\), which is just a statement that there is a high enough RMF frequency for all the electrons, when moving synchronously, to carry the current necessary to reverse \(B_{e}\).

## 4 FRC evolution

Before presenting actual model calculations, it is useful to enumerate the RMF drive parameters for comparison with those of STX. Aside from the ratio of \(\gamma\) to \(\lambda\), \(\omega_{f}\) is the most important parameter. Taking classical values for \(\eta_{\parallel}=515/T_{e}^{3/2}\) (eV) \(\mu\Omega\,\mathrm{m}\) (\(\nu_{\parallel}=1450\times 10^{6}n(10^{20}\ \mathrm{m}^{-3})/T_{e}^{3/2}\) (eV) \(\mathrm{s}^{-1}\)) and basing \(\nu_{\perp}\) on a specified diffusivity ((\(\nu_{\perp}=3.55\times 10^{6}n(10^{20}\ \mathrm{m}^{-3})D_{\perp}(\mathrm{m}^{2}/ \mathrm{s})\mathrm{s}^{-1}\) (or optionally \(v_{\perp}=2v_{\parallel}\))),

\[\lambda=35a(\mathrm{m})\,T_{e}^{3/4}(\mathrm{eV})\,\omega^{1/2}(10^{6}\ \mathrm{s}^{-1}) \tag{42}\]

\[\gamma=0.25\frac{B_{\omega}(\mathrm{G})\,T_{e}^{3/4}\,(\mathrm{eV})}{D_{\perp} ^{1/2}(\mathrm{m}^{2}/\mathrm{s})\bar{n}(10^{20}\ \mathrm{m}^{2}/\mathrm{s})} \tag{43}\]

\[\omega_{f}=\frac{0.1B_{0}(\mathrm{T})}{\bar{n}(10^{20}\,\mathrm{m}^{-3})r_{c}^ {2}(\mathrm{m})(1-\alpha)}\times 10^{6}\,\mathrm{s}^{-1}. \tag{44}\]

It is useful to consider \(x_{s}\) as known, and examine what effect this key FRC parameter, through its determination of \(\langle\beta\rangle=1-\tfrac{1}{2}x_{s}^{2}\), has on RMF drive. Approximating \(a\) as \(r_{s}\), Eq. (41) can be written as \(1-(y/a)^{2}=\omega_{f}/x_{s}^{2}(1-x_{s}^{2})\omega\). Since the maximum value of \(1-(y/a)^{2}\) is \(\tfrac{1}{2}\), \(\omega_{f}/\omega\) must be less than \(0.5x_{s}^{2}(1-x_{s}^{2})\). This has a maximum possible value of \(0.125\) at \(x_{s}=0.707\). Equations (40) and (41) can be used to derive an expression for \(a/r_{s}\),

\[\left(\frac{a}{r_{s}}\right)^{2}=\frac{1-0.25x_{s}^{2}}{1-0.533[1-(y/a)^{2}]}. \tag{45}\]

The density profile shown in Fig. 4 with \(a\) less than \(r_{s}\) will only be valid for \(1-(y/a)^{2}<0.469x_{s}^{2}\). The reason for this latter inequality is that for small \(x_{s}\)\(\langle\beta\rangle\) is large, and there is not much flux inside the separatrix. Since \(1-y/a\) needs to be as large as possible for efficient flux drive, \(x_{s}\) should be very large. In fact, for our simple equilibrium model it will be exceedingly difficult for the RMF to increase the flux of an elongated FRC inside a flux conserver (as will be seen in Fig. 7) if \(x_{s}\) is smaller than about \(0.65\). At high beta there will be steep density gradients near the FRC edge, with large azimuthal currents there. Even if the RMF frequency were made


<!-- Page 11 -->

large enough to exceed the electron rotational speed there, it would not penetrate close to the field null. Flux drive of small-\(x_{s}\) FRCs can only occur through dynamic processes where a field null is generated at \(r>r_{s}/\sqrt{2}\) and flux is convected inwards. This process is seen in both STX internal field measurements and numerical calculations.

When \(a\) is calculated to be larger than \(r_{s}\) (as will occur at smaller \(x_{s}\)), we will use a somewhat different representative density profile, extending it fully to the separatrix. The flexibility necessary for satisfying both the RMF produced reversal current condition and the average beta condition will be provided by allowing \(\bar{n}\) in the expression for \(\omega_{f}\) to increase to as large as \(n_{0}\). In the following calculations \(a\) is then taken as \(r_{s}\). Equations (35) and (36) are then sufficient to determine \(x_{s}\),

\[x_{s}^{2}=\sqrt{\left(\frac{8}{5}\bar{\phi}\right)^{2}+\frac{16}{5}\bar{\phi} -\frac{8}{5}\bar{\phi}}. \tag{46}\]

The average beta condition then requires that \(1-(y/a)^{2}=0.469x_{s}^{2}\), which again points out the importance of having large \(x_{s}\). In order to satisfy the amount of driven current required to reverse the external field, we allow the average density to increase to \(f_{n}n_{0}\) so that, using the original definition of \(\omega_{f}\) with \(\bar{n}=0.5n_{0}\), Eq. (41) becomes

\[\frac{1}{2f_{n}}\frac{\omega_{f}}{\omega x_{s}^{2}(1-x_{s}^{2})}=1-\left(\frac {y}{a}\right)^{2}=0.469x_{s}^{2}. \tag{47}\]

This can be reduced to \(f_{n}=1.07\omega_{f}/x_{s}^{4}(1-x_{s}^{2})\omega\). In the model it is assumed that if \(\omega_{f}/\omega\) is too high, and \(f_{n}\) is calculated as greater than unity for a given \(x_{s}\), then even an average current layer density equal to the peak density \(n_{0}\) is insufficient to carry the equilibrium current if the electrons have a rotational speed equal to \(\omega\). In that case it is assumed that the electrons must rotate faster than \(\omega\) and RMF drive ceases (\(T_{M}\) becomes zero in Eq. (29) and \(\omega_{r}\) becomes correspondingly greater than \(\omega\) in Eq. (32). The maximum value of \(x_{s}^{4}(1-x_{s}^{2})\) is 0.148 at \(x_{s}=0.817\). Thus the maximum value of \(\omega_{f}/\omega\) which can result in current drive in our model is 0.14. This average beta related restriction is somewhat less restrictive than the \(\omega_{f}/\omega<0.125\) restriction mentioned earlier. Realistically, one should keep \(\omega_{f}/\omega\) below about 0.10 to allow for a range of \(x_{s}\), which has not quite been done for the STX experimental conditions listed in Table 1.

Calculations have been performed to model the recent STX experiment. The measured plasma conditions for both the STX experiment and two rotamak experiments are given in Table 1 for the final conditions after flux buildup. The listed values of \(\lambda\) are based on classical resistivity at the listed electron temperature. The listed values of \(\nu_{\perp}\) are calculated for one half the listed maximum density. For the low temperature rotamak cases the perpendicular resistivity is assumed to be classical, but for the STX conditions a value of \(\eta_{\perp}=25\;\mu\Omega\) m is assumed. In all cases \(\gamma\) exceeds \(\lambda\). The rotamak plasmas are spherical and have no flux conserver, so the analysis in this article does not apply. In both cases the separatrix is outside the glass vacuum vessel wall, and wall contact probably limits the size. Thus no value of \(\omega_{f}\) is given. However, we have still listed a value of \(a-y\) needed to carry the reversal current (assumed to be \(2B_{e}/\mu_{0}\) due to the spherical shape) at the widest point. This value is seen to agree well with the measured RMF penetration depth for both the rotamaks and the STX FRC. It is interesting that in all cases the RMF is measured to penetrate approximately to the field null. (During the formation dynamics a field null may be transiently generated at a larger radius than \(r_{s}/\sqrt{2}\), which is not accounted for in our quasi-steady model.)

Internal probe measurements of the axial confinement field and the rotating field are shown in Fig. 5 for the STX experiment. The RMF field for the experiments is started at 150 \(\mu\)s. Its \(B_{\theta}\) amplitude without plasma would be about 20 G everywhere, but this is seen to be doubled outside a plasma column by the axial shielding currents, as shown in Fig. 2. In the experiment a weak FRC is first formed by the standard (but very low voltage) field reversed theta pinch technique [18] and then allowed to decay away before the RMF was applied. This was an easy way to start with full ionization. An FRC was recreated by the RMF and its flux built up to about 0.35 mWb in 100 \(\mu\)s. During that time the RMF is seen to just penetrate to the field null at about \(r=15\) cm. After that time the FRC flux is seen to slowly decay. This is in contrast to the rotamak results, where the configuration is maintained in steady state for over 20 ms. The reasons for this are uncertain, but they could be due to the minimum conditions necessary for flux buildup and penetration outlined in this article not being met. The ratio of \(\omega_{f}/\omega=0.12\) listed on Table 1 is fairly high, so that a further increase in electron temperature could reduce the average electron density to a value which, when the electrons are rotating synchronously, is insufficient to provide for a full \(2B_{e}\) field change. The only possible equilibrium would then be the


<!-- Page 12 -->

\begin{table}
\begin{tabular}{l c c c} \hline Property & Flinders 10 L & Flinders 50 L & STX \\ & rotamak [6] & rotamak [7] & FRC [14] \\ \hline \(r_{s}\) (cm)/\(r_{s}\) (cm) & 15 & 25 & 20/25 \\ \(\omega\) (10\({}^{6}\) s\({}^{-1}\)) & 3.1 & 3.1 & 2.2 \\ \(B_{\omega}\) (G) & 100 & 23 & 20 \\ \(B_{e}\) (G) & 90 & 23 & 90 \\ \(T_{e}\) (eV) & \(\sim\)10 & \(\sim\)17 & \(\sim\)50 \\ \(\delta\) (cm) & 0.3 & 0.2 & 0.1 \\ \(\lambda=r_{s}\delta\) & 50 & 130 & 200 \\ \(n_{e0}\) (10\({}^{19}\) m\({}^{-3}\)) & 4 & 0.3 & 0.4 \\ \(\nu_{\perp}(10^{6}\) s\({}^{-1}\)) & 18 & 0.6 & 1.4 \\ \(\gamma=\omega_{ce}/\sqrt{\nu_{\parallel}\nu_{\perp}}\) & 140 & 950 & 730

non-field-reversed, high beta, theta pinch plasma which the column is observed to evolve towards. It should be mentioned that the RMF power was continually decaying in the STX experiments, and that continual flux buildup could be observed under different operating conditions.

The STX results are unique, compared with those of the spherical rotamaks, for several reasons. The plasmas are elongated and confined in a flux converter. They thus compress the external field, limiting their radial expansion (at peak flux the separatrix is slightly beyond the 20 cm inner quartz plasma tube radius, but no density is measured there), and they must satisfy the average beta condition. For the first time in an RMF driven FRC the external field is also significantly higher than the rotating field. This was also seen in rotamak ST experiments where a toroidal field presumably prevented unbridled expansion and wall contact [19].

A calculation designed to approximately mimic the STX results is shown in Fig. 6. The vacuum field \(B_{0}\) is chosen as 40 G so that \(B_{e}\) will increase to about the measured 90 G at large \(x_{s}\). The perpendicular diffusion coefficient was kept constant at a 20 m\({}^{2}\)/s value typical of that inferred from flux decay rates measured without RMF drive. Initial conditions were specified so that the final density and temperature would reach close to the measured values after flux buildup. The temperature and density in all the calculations are assumed to increase adiabatically and isentropically (\(n\propto B_{e}^{1.2}\), \(T_{e}\propto B_{e}^{0.8}\)). For the assumed fixed cross-field resistivity, the ratio of \(\gamma\) to \(\lambda\) is not dependent on electron temperature but decreases linearly with increasing density. The ions are cold in the STX experiments, and \(\alpha\) is taken as zero.

Figure 5: Internal probe measurements of equilibrium \(B_{z}\) and RMF \(B_{\theta}\) fields in the STX experiment. The red \(B_{\theta}\) traces are for a vacuum reference.


<!-- Page 13 -->

In order to start the flux buildup from a small radius in the calculations, the value of \(a\) is initially allowed to exceed \(r_{s}\). This is seen to occur for the first 20 \(\mu\)s in the calculations when the initial \(x_{s}\) is taken as 0.5 and there are too few electrons in the necessarily thin current layer (thin due to a high average beta) to carry the reversal current. Otherwise, the equilibrium model presented here would not allow flux buildup to occur under such conditions. (In the STX experiments, flux buildup begins from a large radius mirror confined high beta plasma extending nearly to the plasma tube wall.) From 20 to 60 \(\mu\)s, with \(x_{s}\) ranging from 0.63 to 0.77, \(a\) is taken as equal to \(r_{s}\) in the calculations, and the average density in the current layer falls from the peak density to one half of the peak density. The plotted values of \(\gamma\) are based on \(\bar{n}\), which is the reason why \(\gamma\) is increasing during this time.

It may be noted that the model distance \(y\) is larger than \(r_{s}/\sqrt{2}\), so that the model RMF penetration does not reach the quasi-steady field null. This will always be the case for this simple model, and the effect on flux buildup is only accounted for by the reduction in \(T_{M}\) given by Eq. (32). We can do no more with this simple analytic approach. If the calculated value of \(y\) is too far from \(r_{s}/\sqrt{2}\) (as will occur when the density is large, as reflected in a small value of \(\omega_{f}/\omega\), then it may not be possible to experimentally obtain a quasi-steady (field null at \(r_{s}/\sqrt{2}\)) solution, and either flux buildup will cease, or will only proceed dynamically. Such dynamic effects are seen in numerical solutions, with the internal flux increasing in an oscillating manner, but this is beyond the scope of the present quasi-steady analysis. Our analysis, based on overall torques, only shows the minimum conditions necessary for flux buildup to occur. However, these minimum conditions also appear sufficient in many cases, due to reduction of the current at the field null and some RMF penetration to that position, when compared with numerical calculations [16].

Flux buildup is fairly rapid for this small FRC, with a peak rate of about 0.0035 mWb/\(\mu\)s (3.5 Wb/s), agreeing roughly with the experimental measurements, which also must be limited somewhat by RMF diffusion times. Using \(\phi\approx 0.67\pi a^{2}[1-(y/a)]B_{e}\) and \(V_{\theta}=d\phi/dt\) given by Eq. (29), a useful numerical relation for the flux buildup time is

\[\begin{split}\tau_{bu}&=\frac{\phi}{V_{\theta}}\\ &=390\,\frac{a^{2}(\mathrm{m})B_{e}(\mathrm{T})\bar{n}(10^{20}\ \mathrm{m}^{-3})}{B_{\omega}^{2}(\mathrm{G})}\left(1-\frac{T_{\eta}}{T_{M}} \right)^{-1}\mathrm{s}.\end{split} \tag{48}\]

The model solution shows an equilibrium flux level being reached, which occurs at the time that \(T_{\eta}\) equals \(T_{M}\). There are many reasons, other than the ratio \(\omega_{f}/\omega\) becoming too large, why the experimental results could show a decaying flux. The equilibrium is due to a fine balance between \(T_{M}\) and \(T_{\eta}\), which can easily be disturbed by a change in density or external field. The effect of an increasing magnetic field is shown in the Fig. 7 calculations. In this case \(B_{0}\) is increased sinusoidally by a factor of 2 over 500 \(\mu\)s. \(x_{s}\) never becomes large enough to produce a low beta profile with a 'vacuum' region near the separatrix. The increasing density also causes the relative resistive torque, as represented by a decreasing value of \(\gamma\), to become too large. When \(x_{s}\) becomes too small, at about 375 \(\mu\)s, there are not enough electrons in the thin current layer to reverse the external field. In the model this means RMF flux drive ceases. The


<!-- Page 14 -->

FRC flux then decays at a rate characteristic of the assumed resistivity (\(D_{\perp}=\eta_{\perp}/\mu_{0}=20\) m\({}^{2}\)/s). The \(\sim\)30 \(\mu\)s decay timescale is in rough agreement with experimental results if the RMF is shut off early. The limitation on flux buildup in a rising external field (which is what is desired experimentally) can be rectified analytically by allowing the cross-field resistivity to decrease with rising temperature, or by stipulating a lower value of \(D_{\perp}\).

A final calculation is shown to elucidate what would be necessary to achieve about 100 times the flux level that was produced in STX. A flux observer radius of \(r_{c}=50\) cm and a final field of \(B_{e}=2.5\) kG result in \(\sim\)30 mW of flux at an \(x_{s}\) value of about 0.8. Various combinations of electron density and RMF frequency are possible to produce the necessary synchronous current. A value of \(n_{0}=2\times 10^{20}\) m\({}^{-3}\) has been chosen for the example calculation. Using the above \(x_{s}\) value the vacuum field \(B_{0}\) would be about 1 kG. The value of \(\omega_{f}\) for those conditions is \(0.04\times 10^{6}\) s\({}^{-1}\), so a reasonable RMF frequency is of the order of \(\omega=0.5\times 10^{6}\) s\({}^{-1}\). The ratio of \(\gamma\) over \(\lambda\), as given by Eqs (42) and (43) is

\[\frac{\gamma}{\lambda}=0.0071\] \[\times\frac{B_{\omega}(\mathrm{G})}{\bar{n}(10^{20}\ \mathrm{m}^{-3})D_{\perp}^{1/2}( \mathrm{m}^{2}/\mathrm{s})\omega^{1/2}(10^{6}\ \mathrm{s}^{-1})a(\mathrm{m})}. \tag{49}\]

Choosing \(\bar{n}=1.0\times 10^{20}\) m\({}^{-3}\), \(\omega=0.5\times 10^{6}\) s\({}^{-1}\) and an optimistic value of \(D_{\perp}=2\) m\({}^{2}\)/s (about the best measured for high density FRCs [20] and corresponding classically to an electron temperature of 55 eV) the calculations were performed for a value of \(B_{\omega}=100\) G, yielding a value of \(\gamma/\lambda=1.8\) at the final conditions, which is high enough to allow continued flux buildup or sustainment.

Results from the calculation are shown in Fig. 8 for a vacuum external field rising from 0.4 to 0.8 kG in 5 ms. Since the classical resistivity is low, \(D_{\perp}=\eta_{\perp}/\mu_{0}\) was kept constant at 2 m\({}^{2}\)/s throughout the calculation. If the low field results are indicative of high field behaviour, the success of such future experiments will depend primarily on the value of the anomalous resistivity. This key resistivity requirement is discussed further in Section 5.

The calculation shown in Fig. 8 was also repeated for a value of \(\alpha=-1\). Since the ions then carry half the current, the value of \(\omega_{f}\) is reduced by a factor of two and the RMF frequency \(\omega\) has to be cut in half. The net result is to have no effect on the flux buildup and FRC evolution. The value of \(T_{M}\) in Eq. (31) is unchanged. \(\lambda\) is reduced by a factor of \(1/\sqrt{2}\), as indicated by the dashed curve in Fig. 8, but the resistive torque given by Eq. (32) is unchanged. There would only be an effect if the RMF were nearly fully penetrated, and Eq. (20) were used for the RMF torque, but this will almost never occur within the framework of this edge model.

Ion azimuthal current has no effect on the flux buildup physics since the resistive torque that must be overcome is proportional to the difference between \(v_{e\theta}\) and \(v_{i\theta}\) (proportional to \(j_{\theta}\)). When making simple comparisons of \(\lambda\) with \(\gamma\), one should use an effective frequency of \((1-\alpha)\omega\) in computing an effective \(\lambda\). This value is unchanged by ion current contributions. The effect of ion current is important, however, in determining the appropriate RMF frequency to use. It can be beneficial in lowering the required frequency, although it has no effect on the required


<!-- Page 15 -->

\(B_{\omega}\) magnitude. The effect on frequency can be significant since the flux buildup and sustainment process is very sensitive to the correct value of \(\omega_{f}/\omega\). This is shown by the calculations illustrated in Fig. 9 of the final flux as a function of \(\omega\) reached in the calculations illustrated in Fig. 8. For too small an RMF frequency the FRC current will be such that the average electron azimuthal rotational velocity exceeds the RMF frequency, and RMF penetration will be impossible. The model calculations show a moderately quickly decreasing performance with too high an RMF frequency due to the current being carried in an ever narrowing sheath thickness. In reality the results may be far worse if the sheath is so thin that no RMF field can penetrate to the field null. Unfortunately, the model presented here cannot predict exactly when that may occur. It is safe to say that experiments should be carefully planned to optimize the \(\omega_{f}/\omega\) ratio.

## 5 Summary and conclusions

The modelling presented here has profound implications for FRC physics. The creation of a toroidal plasma, with a true vacuum boundary, is fairly unique. In addition to the flux being maintained, particle losses are also strongly restricted. Refuelling, originally thought necessary for long term operation of present experiments may not be necessary. The RMF also produces a \(\langle j_{z}B_{\theta}\rangle\) inward force which has a sharp gradient for the type of edge currents that are produced. This should be highly stabilizing, which may account for the observed stability of the cold, collisional, STX FRCs. The analytic results show the type of scaling which must be followed in extrapolating to larger devices. Parameters for the present STX experiment, a planned upgrade, the newly constructed TCS (Translation, Confinement and Sustainment), a proof of principle (POP) level device, and a reactor are shown in Table 2. The after flux buildup conditions were chosen to satisfy both \(\gamma>\lambda\) (for conceivable values of \(D_{\perp}\)) and \(\omega_{f}/\omega\approx 0.1\).

On the basis of \(\omega\propto\omega_{f}\) as given in Eq. (44) and \(\gamma/\lambda\) given by Eq. (49), the required RMF strength \(B_{\omega}\propto(D_{\perp}B_{e}n_{e})^{1/2}\). \(B_{\omega}\sim 100\) G is a realistic upper limit based on both technological and input power considerations, and the required values of \(\eta_{\perp}=\mu_{0}D_{\perp}\) can be found by setting \(\gamma/\lambda\approx 2\) in Table 2. For the STX cases the ions are assumed to be cold. They are assumed to be equilibrated with the electrons for the other cases, but the electrons are assumed to carry all the azimuthal current for all cases, so that \(\alpha=0\) in the previous equa

Figure 8: Calculation of flux buildup due to \(B_{\omega}=100\) G, \(\omega=0.5\times 10^{6}\) s\({}^{-1}\) RMF, \(n_{0}=0.3\times 10^{20}\) m\({}^{-3}\) and \(T_{e0}=80\) eV, with \(D_{\perp}=2\) m\({}^{2}\)/s. Compression is assumed with \(B_{0}\) increasing from 0.4 to 0.8 kG.

Figure 9: Effect of RMF frequency on final flux for the conditions of Fig. 8


<!-- Page 16 -->

**Article: Rotating magnetic field current drive of FRCs**


<!-- BLOCK FAILED: 2000_Hoffman_page016_block2 (see markdown_new\2000_Hoffman\2000_Hoffman_page016_block_fail_t34-b64_l0-r612.pdf) -->



<!-- BLOCK FAILED: 2000_Hoffman_page016_block3 (see markdown_new\2000_Hoffman\2000_Hoffman_page016_block_fail_t64-b217_l0-r612.pdf) -->


tions. The \(\gamma/\lambda\) values are based on an average density of one half the \(n_{e0}\) value listed. The FRC flux is based on a ratio of \(a/r_{c}=0.8\), \(x_{s}=0.83\) and \(y/a=0.85\), typical of the analytic calculations. The expression for flux, based on Eq. (34) then becomes \(\phi=0.24\pi r_{s}^{2}B_{e}\), which is very close to the rigid rotor value of \(0.31x_{s}\pi r_{s}^{2}B_{e}\) for the chosen value of \(x_{s}\). An important point, based on this analysis, is that the ratio of \(\omega/\omega_{f}\) be kept at a value close to 0.1 in order to be consistent with the FRC equilibrium. This, in turn, determines the RMF frequency \(\omega\), which must be smaller for larger devices since the current carrying radius is larger.

The column labelled STX/ug involves only a straightforward extension of the STX RMF power supply, and should be relatively straightforward. Varying the magnitude of the RMF needed to operate at a certain density will provide an indication of the value of the effective diffusivity. The POP parameters involve a hundred-fold increase in flux over the STX parameters, logarithmically halfway to the reactor flux level. If small scale experiments are fully successful, this could be an ultimate goal of the new TCS device. The analytic results show the transport scaling needed for such an experiment to be fully successful. The cross-field resistive diffusion coefficient must be about 4 m\({}^{2}\)/s in TCS to achieve the desired ratio of \(\gamma/\lambda\), and even lower for the POP device or reactor.

The reactor conditions listed in Table 2 are a relatively straightforward extension from the POP conditions, with the device size and field being increased by factors of about 5. Of course, producing the high electron temperature is a separate question. The ratio of \(B_{e}/n_{e}\) increases by a factor of 2, but \(r_{c}^{2}\) increases by a factor of 25, requiring a tenfold decrease in RMF frequency. Depending on the realized value of cross-field resistivity, the RMF may need to be strong. For the listed \(B_{\omega}\) of 100 G, a cross-field diffusivity under 1 m\({}^{2}\)/s is necessary. This is probably as high as is tolerable since a current of 2 MA/m is needed to reverse the FRC external field, and a 100 G RMF at the listed density produces an 8 V loop voltage, implying an RMF maintenance power of about 16 MW/m.

The 100 G RMF field strength for a 2.5 m radius device is also about at the limit of what is technologically feasible. The voltage required to drive the antenna current \(I_{a}\) scales as \(V_{a}=I_{a}\omega L_{a}\) where \(L_{a}\), the antenna impedance for a given length, scales as \(r_{c}^{2}\). Since the required antenna current scales as \(B_{\omega}r_{c}\), and the RMF frequency scales as \(B_{e}/nr_{c}^{2}\), the antenna voltage \(V_{a}\propto B_{\omega}(B_{e}/n)r_{c}\). Ten kilovolt voltages will be used on TCS, driving a parallel tuned circuit (the antenna is the inductance), so that about 250 kV would be required for the listed reactor conditions.

One way of alleviating the RMF power requirements would be to centrally fuel the FRC. This would tend to provide an overall \(\langle v_{r}B_{z}\rangle\) driving force to complement the RMF drive. The strong RMF drives produce a net torque on the plasma, which will eventually spin up the total fluid if not compensated. Fuelling alone, as discussed by Ohnishi and Ishida [17], cannot provide this compensation in a reactor since the required fuelling rate \(\nu_{n}\sim(m_{e}/m_{i})\nu_{\perp}\) is too large. If fuelling can reduce the required RMF torque, then neutral beam injection in the ion diamagnetic direction should be able to provide the necessary source of oppositely directed momentum.

The results given in this article are by no means definitive of RMF-FRC interactions. The FRC beta


<!-- Page 17 -->

condition can be changed by altering its shape and, for long timescales, the flux conserver properties will be modified by external field programming. Key questions remain as to how well the FRC profile can adjust to provide sufficient \(E_{\theta}\) at the field null, and how well the essentially synchronous electron current can be maintained as \(\lambda\) becomes very large. End effects relating to the two dimensional nature of the FRC have also not been studied. However, the analytic results presented here are consistent with experimental measurements on small devices with significant \(\lambda\) values, and should serve as a valuable scaling tool in planning future experiments.

## Appendix

If electron momentum is non-negligible, Eq. (1) is modified to be

\[m_{e}\frac{d\mathbf{v}_{e}}{dt}=-e(\mathbf{E}+\mathbf{v}_{e}\times\mathbf{B}-\eta\mathbf{j}) \tag{50}\]

where \(d/dt\) is the total derivative \(\partial/\partial t+\mathbf{v}\cdot\mathbf{\nabla}\). Thus, using the previous \(e^{i(\omega t-\theta)}\) notation, the \(z\) component of \(dv_{ez}/dt\) becomes \(i(\omega-v_{e\theta}/r)v_{ez}\). Following the same analysis as used in deriving Eq. (9), and again calling \(v_{e\theta}=r\omega_{r}\) and \(\varpi=\omega-\omega_{r}\), the diffusion equation for \(E_{z}\) is modified to

\[i\varpi E_{z}=\frac{\eta_{\perp}}{\mu_{0}}\left(1+i\frac{\varpi}{\nu_{\parallel }}\right)\nabla^{2}E_{z}. \tag{51}\]

A Bessel equation similar to Eq. (10) follows, except in this case \(k\) is complex. Putting \(\varepsilon=\varpi/\nu_{\parallel}\),

\[r^{2}\,\frac{\partial^{2}E_{z}}{\partial r^{2}}+r\frac{\partial E_{z}}{ \partial r}-(1+k^{2}r^{2})E_{z}=0 \tag{52}\]

where

\[k^{2}=i\,\frac{\mu_{0}}{\eta_{\parallel}}\frac{\varpi}{1+i\varepsilon}=i\left( \frac{\eta_{\parallel}}{\mu_{0}\varpi}+i\frac{m_{e}}{\mu_{e}ne^{2}}\right)^{-1}.\]

Recognizing \((m_{e}/\mu_{0}ne^{2})^{1/2}\) as \(\delta_{e}=c/\omega_{pe}\), the collisionless skin depth, \(k\) takes the form

\[k^{2}=\frac{i}{0.5\delta^{*2}+i\delta_{e}^{2}},\quad k=\frac{1}{ \delta_{T}}\exp\left[i\left(\frac{\pi}{4}-\frac{\theta_{\varepsilon}}{2} \right)\right] \tag{53}\]

where \(\delta_{T}\) is defined as \((0.5\delta^{*2}+\delta_{e}^{2})^{1/2}\) and \(\theta_{\varepsilon}=\tan^{-1}\varepsilon\). Now the solution for \(E_{z}\) is the same as given by Eq. (11), only modified by the different expression for \(k\)

\[E_{z}(r)=\frac{2\omega}{k}\,B_{\omega}\,\frac{I_{1}(kr)}{I_{0}(ka)}. \tag{54}\]

In the limit of large \(kr\), \(E_{z}\) is given by

\[E_{z} \approx \exp\left[-i\left(\chi+\frac{a-r}{\delta_{T}}\sin\chi\right)\right] \tag{55}\] \[\times 2\omega\sqrt{\frac{a}{r}}\,\delta_{T}\,\exp\left[-\left(\frac{a- r}{\delta_{T}}\right)\cos\chi\right]\]

where \(\chi=(\pi/4)-\theta_{\varepsilon}/2\). For \(\varepsilon=0\) this is the same as Eq. (23). The electron skin depth then plays no role since the electrons behave resistively. In the opposite inductive limit, where \(\varpi\gg\nu_{\parallel}\) the electron skin depth is dominant and the phase of \(E_{z}\) is different by \(45^{\circ}\). The ions display this type of behaviour, but due to electron-ion interactions it can be shown that RMF penetration is not limited to an ion skin depth. For moderate \(\varpi\sim\nu_{\parallel}\) the RMF penetration actually increases due to electron inertia, but the azimuthal force exerted on the electrons decreases since their oscillation becomes out of phase with \(B_{z}=E_{z}/\omega r\). This can be seen from the equation for axial current oscillation

\[j_{z}=(1+i\varepsilon)^{-1}\frac{\varpi}{\omega}\frac{E_{z}}{\eta_{\parallel}}. \tag{56}\]

The force \(F_{\theta}=\langle j_{z}B_{r}\rangle\) then becomes

\[F_{\theta}=\frac{1}{1+\varepsilon^{2}}\frac{\varpi}{\omega}\frac{\langle E_{z }^{2}\rangle}{\omega r\eta_{\parallel}}. \tag{57}\]

This is reduced by electron inertia from the value given by Eq. (16) when \(\varpi\) in not small compared with \(\nu_{\parallel}\).

## Acknowledgements

The author would like to acknowledge extensive conversations with Drs R. Milroy and J. Slough and Mr. K. Miller on results from numerical RMF calculations and the STX experiments, which have stimulated the development of this analytical model. He also appreciates the extremely careful reading, and resultant improvement, by one of the referees. This work was supported by USDOE Grant No. DE-FG03-96ER54376.

## References

* [1] Blevin, H.A., Thonemann, P.C., Nucl. Fusion: 1962 Suppl., Part 1 (1962) 55.
* [2] Hugrass, W.N., et al., Phys. Rev. Lett. **44** (1980) 1676.
* [3] Hugrass, W.N., Jones, I.R., Phillips, M.G.R., J. Plasma Phys. **26** (1981) 465.


<!-- Page 18 -->

* [4] Durance, G., Hogg, G.R., Tendys, J., Watterson, P.A., Plasma Phys. Control. Fusion **29** (1987) 227.
* [5] Zwi, H.R., Kuthi, A., Wong, A.Y., Wells, B., Phys. Fluids B **3** (1991) 126.
* [6] Donaldson, N., Euripides, P., Jones, I.R., Xu, S., Plasma Phys. Control. Fusion **37** (1995) 209.
* [7] Euripides, P.E., Jones, I.R., Deng, C., Nucl. Fusion **37** (1997) 1505.
* [8] Hugrass, W.N., Grimm, R.C., J. Plasma Phys. **26** (1981) 455.
* [9] Ohnishi, M., Ishida, A., Yamamoto, Y., Yoshikawa, K., Trans. Fusion Technol. **27** (1995) 391.
* [10] Milroy, R.D., Phys. Plasmas **6** (1999) 2771.
* [11] Jones, I.R., Hugrass, W.N., J. Plasma Phys. **26** (1981) 441.
* [12] Hugrass, W.N., J. Plasma Phys. **28** (1982) 369.
* [13] Hoffman, A.L., Phys. Plasmas **5** (1998) 979.
* [14] Slough, J.T., Miller, K.E., Phys. Plasmas **7** (2000) 1945.
* [15] Hugrass, W.N., Aust. J. Phys. (1986) 513.
* [16] Milroy, R.D., A Magnetohydrodynamic Model of RMF Current Drive in FRCs (in preparation).
* [17] Ohnishi, M., Ishida, A., Nucl. Fusion **36** (1996) 232.
* [18] Tuszewski, M., Nucl. Fusion **28** (1988) 2033.
* [19] Jones, I.R., Deng, C., El-Fayoumi, I.M., Euripides, P., Phys. Rev. Lett. **81** (1998) 2072.
* [20] Hoffman, A.L., Slough, Nucl. Fusion **33** (1993) 23.

(Manuscript received 22 November 1999

Final manuscript accepted 4 May 2000)

E-mail address of A.L. Hoffman:

hoffman@aa.washington.edu

Subject classification: H1, Ct; B0, Ct