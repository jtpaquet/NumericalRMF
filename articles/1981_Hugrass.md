---
paper_id: "1981_Hugrass"
title: "J. Plasma Physics_ (1981), _vol. 26_, _part 3_, _pp_. 455-464"
authors: ""
year: "1981"
institution: "School of Physical Sciences, The Flanders University of South Australia,"
publisher: ""
doi: ""
keywords: ""
abstract: ""
num_references: "7"
num_citations: ""
source_file: "1981 Hugrass - Numerical study of the generation of an azimuthal current in a plasma cylinder using a transverse rotating magnetic field.pdf"
---



<!-- Page 1 -->

J. Plasma Physics_ (1981), _vol. 26_, _part 3_, _pp_. 455-464

_Printed in Great Britain_

A numerical study of the generation of an azimuthal current in a plasma cylinder using a transverse rotating magnetic field

**By W. N. HUGRASS+ and R. C. GRIMM+**

School of Physical Sciences, The Flanders University of South Australia,

South Australia 5042, Australia

(Received 22 September 1980)

The generation of a steady azimuthal current in a cylindrical plasma column using a rotating magnetic field is numerically investigated. The mixed initial-boundary-value problem is solved using a finite difference method. It is shown that substantial azimuthal current can be driven provided that the amplitude of the rotating magnetic field is greater than a certain threshold value which depends on the plasma resistivity.

## 1 Introduction

In the preceding paper (Jones & Hugrass 1981), the technique, first investigated by Blevin & Thonemann (1962), of generating unidirectional azimuthal electron currents in a plasma cylinder by means of radio-frequency rotating magnetic fields is described and discussed. In particular, steady-state solutions are presented which show that, provided the amplitude and rotation frequency of the field are suitably chosen, the required penetration of the field into the plasma is not limited by the usual classical skin effect. These steady-state solutions are obtained for two limiting cases: the case where the resistive term \(\eta\mathbf{j}\) in Ohm's law is considered dominant and the case where the plasma behaviour is considered to be mainly governed by the nonlinear \(\mathbf{j}\times\mathbf{B}\) Hall term in Ohm's law. The purpose of this paper is to present the results of a numerical treatment of the corresponding initial-value problem. Such a study allows one to examine the question of accessibility of the steady state and, moreover, to incorporate arbitrary contributions from the \(\eta\mathbf{j}\) and \(\mathbf{j}\times\mathbf{B}\) terms. The modifications introduced by anisotropic resistivity and the presence of flux conserving rings are also examined.

## 2 The initial value problem

### Physical model

We consider an infinitely long cylindrical plasma column of radius \(r=a\) which is surrounded by a vacuum region enclosed by a set of flux conserving rings at \(r=b\). We assume that these rings conduct in the azimuthal directions only and


<!-- Page 2 -->

that they are sufficiently far from the plasma column that the system has no \(z\) variation (\(\partial/\partial z=0\)).

The plasma is assumed cold and stationary with uniform electron number densities. Initially the current density is zero but there may be a uniform externally applied axial magnetic field \(B_{a}\). At \(t=0\) we apply a transverse rotating magnetic field,

\[B_{r}=B_{t}(t)\cos{(\omega t-\theta)}, \tag{1}\]

\[\dot{B}_{\theta}=B_{t}(t)\sin{(\omega t-\theta)}, \tag{2}\]

where the amplitude rises in a time \(\tau_{r}\) characteristic of the RF source used to generate the rotating field, i.e.

\[B_{t}(t)=\left(1-e^{-t|\tau_{r}}\right)B_{\omega}, \tag{3}\]

\[t\geq 0.\]

We study the penetration of this field into the plasma column and the generation of a steady azimuthal current. If the frequency, \(\omega\), of the applied field satisfies the condition \(\omega\geq\omega_{el}\), then, to a good approximation, the ions are not affected by the rotating field and form a neutralizing background of positive charges. Assuming that the system is initially in equilibrium, we can neglect the ion motion and study the field penetration using the equation of motion for the electrons.

\[\mathbf{E}=\eta\mathbf{J}+\frac{1}{ne}\,\mathbf{J}\times\mathbf{B}+\frac{m}{ne ^{2}}\frac{\partial\mathbf{J}}{\partial t}, \tag{4}\]

together with Maxwell's equations

\[\nabla\times\mathbf{B}=\mu_{0}\mathbf{J}, \tag{5}\]

\[\nabla\times\mathbf{E}=-\frac{\partial\mathbf{B}}{\partial t} \tag{6}\]

and

\[\nabla\,.\,\mathbf{B}=0. \tag{7}\]

We assume a classical resistivity tensor (Braginskii 1965)

\[\eta=\eta_{\perp}\left(\mathbf{I}-\frac{\mathbf{B}\mathbf{B}}{B^{2}}\right)+ \eta_{\perp}\frac{\mathbf{B}\mathbf{B}}{B^{2}} \tag{8}\]

with

\[\eta_{\perp}=\gamma\eta_{\perp}=\gamma\frac{mv_{ei}}{ne^{2}}, \tag{9}\]

Here \(\gamma\) is a numerical coefficient which can be chosen phenomenologically in the calculations; in a strong magnetic field (\(\omega_{ce}>v_{ei}\)), \(\gamma=1\cdot 96\) (Brakinskii 1965). In a cold resistive plasma, typical of the majority of experimental measurements, \(\omega\ll\nu_{ei}\) and it is reasonable to neglect the electron inertia term in (4). At high temperatures or large rotating field frequencies \(\omega\lesssim\nu_{ei}\) and this term should be retained.


<!-- Page 3 -->

### The basic equations

The magnetic field can be expressed in the form

\[\mathbf{B}=\nabla\times\left(A_{s}\mathbf{e}_{z}\right)+B_{s}\mathbf{e}_{z} \tag{10}\]

where \(A_{s}\) is the \(z\) component of the vector potential and \(\mathbf{e}_{z}\) is the unit vector in


<!-- BLOCK FAILED: 1981_Hugrass_page003_block2 (see markdown_new\1981_Hugrass\1981_Hugrass_page003_block_fail_t126-b256_l0-r459.pdf) -->


Equations (12), (13) form a coupled set of nonlinear partial differential equations which, with the appropriate boundary conditions, describe the evolution of the magnetic field in the plasma from the initial conditions,

\[A_{g}(r,\theta,0)=0\quad\quad\text{for}\quad r\leq a, \tag{14}\] \[B_{g}(r,\theta,0)=B_{a}\quad\text{for}\quad r\leq b. \tag{15}\]

### The boundary conditions

Boundary conditions are required for \(A_{g}\) and \(B_{g}\) at the plasma vacuum interface, \(r=a.\) Since no currents flow in the vacuum region,

\[B_{g}(r,\theta,t)=B_{g}(a,t) \tag{16}\]

for \(a\leq r\leq b.\) Noting that \(B_{g}\) is continuous at \(r=a,\) the appropriate boundary condition for \(B_{g}\) follows from the conservation of axial flux inside the conducting rings.

\[B_{g}(a,t)=\frac{b^{2}B_{a}(0)}{\frac{1}{\alpha}\frac{1}{\alpha}-\frac{1}{ \alpha}\frac{1}{\alpha}}\left[{}^{a}rdr\left[{}^{2\pi}d\theta B_{g}(r,\theta,t).\right.\right. \tag{17}\]


<!-- BLOCK FAILED: 1981_Hugrass_page003_block4 (see markdown_new\1981_Hugrass\1981_Hugrass_page003_block_fail_t474-b542_l0-r459.pdf) -->


well as the externally applied rotating field. Straightforward calculation of the vector potential due to the internal currents leads to serious numerical difficulties because these currents are such as almost to cancel the imposed external field (during the linear phase of the calculation) and a solution of the equation

\[\nabla^{2}A_{z}=0,\ \ \ r\geqslant a, \tag{18}\]

is required. This solution can be expressed as the sum of harmonics

\[A_{z}=\Sigma_{m}A_{m}e^{in\theta}, \tag{19}\]


<!-- Page 4 -->


<!-- BLOCK FAILED: 1981_Hugrass_page004_block1 (see markdown_new\1981_Hugrass\1981_Hugrass_page004_block_fail_t0-b197_l0-r459.pdf) -->



<!-- BLOCK FAILED: 1981_Hugrass_page004_block2 (see markdown_new\1981_Hugrass\1981_Hugrass_page004_block_fail_t197-b390_l0-r459.pdf) -->


\[A_{0}=a\ln\frac{a}{R}\frac{\partial A_{0}}{\partial r}\bigg|_{r=a}. \tag{25}\]

For the particular form of external field studied in this paper (see (1), (2)),

\[A_{\rm ext}=aB_{t}(\cos\omega t+i\sin\omega t), \tag{26}\]

\[A_{\rm ext}=0,\ \ \ m\neq 1. \tag{27}\]

## 3 The numerical method

The coupled equations (11), (12) subject to the initial and boundary conditions (14), (15), (17), (24) and (25) are solved using a time centred finite difference method on a uniformly spaced polar mesh,

\[\begin{cases}A_{s}(t+\Delta t)\\ B_{s}(t+\Delta t)\end{cases}-\begin{cases}A_{s}(t)\\ B_{s}(t)\end{cases}=0\cdot 5\,\Delta tD\begin{cases}A_{s}(t)\\ B_{s}(t)\end{cases}+0\cdot 5\,\Delta tD\begin{cases}A_{s}(t+\Delta t)\\ B_{s}(t+\Delta t)\end{cases} \tag{28}\]

where the operator \(D\) is the space centred finite difference approximation of the right-hand side of (11), (12). A fully implicit solution of (28) is not possible because of the nonlinearity introduced by the Hall term in (4). The linear terms are treated semi-implicitly and the resulting system of algebraic equations is solved


<!-- Page 5 -->

using a successive over-relaxation scheme which is conveniently organized to include the iteration over the nonlinear part at interior points as well as the iteration involved in satisfying the boundary conditions. The first approximation in the iterative process is obtained by linear extrapolation from the previous time step.

The Fourier transforms in the boundary condition (24) are evaluated using a modified fast Fourier transform (FFT) algorithm which makes use of the symmetry properties of the discrete Fourier transforms (DFT) of real variables. The point \(r=0\) requires special treatment to remove the apparent singularity in (11), (12) when written explicitly in polar co-ordinates. This is performed by assuming that \(A_{z},B_{z}\) are constant over the first radial zone and using Gauss's theorem.

The calculations were carried out on a DEC-10 computer at Flinders University, a typical simulation required 1-2 cpu hours. Usually sufficient accuracy was attained with a \(16\times 16\) mesh in \((r,\theta)\), but for small values of \(\eta\) (i.e. small values of \(\delta\)), radial mesh refinement was required. Although at early times in the simulation (during the linear phase) a time step much larger than the Courant condition (Richtmyer & Morton 1967) could be tolerated, efficient convergence of the iteration over the nonlinear part, which becomes dominant during the later phase of the simulation, usually sets a practical limit to the time step of about 10 times the Courant condition at the origin.

## 4 Results

A series of numerical calculations was performed for a plasma column of radius \(a=0\cdot 02\) m and a uniform number density \(n_{e}=10^{20}\) m\({}^{-3}\). The externally applied rotating field has an angular frequency \(\omega=2\pi f=5\times 10^{6}\) sec\({}^{-1}\) and its amplitude rises to its steady-state values in an \(e\)-folding time, \(\tau_{\tau}=0\cdot 4\)\(\mu\)sec.

In the ideal case, the electron fluid rotates as a rigid body at the rotating field frequency. The corresponding azimuthal current per unit length is

\[I_{\theta p}=\pi a^{2}n_{e}ef. \tag{29}\]

The azimuthal current per unit length is obtained from the numerical calculation using Ampere's law,

\[I_{\theta}=\frac{1}{\mu_{0}}\left|B_{z}(a)-B_{z}(0)\right| \tag{30}\]

and the ratio \(\alpha(t)=I_{\theta}(t)/I_{\theta p}\) can be used as a quantitative measure of the penetration of the rotating field into the plasma; in particular we are interested in the steady-state value, \(\alpha_{s}\) and the transient time \(\tau\) defined by the relation

\[\alpha(\tau)=0\cdot 5\,\alpha_{s}. \tag{31}\]

### The threshold regime

Figure 1 shows a plot


<!-- Page 6 -->

is neglected (since \(\nu_{ei}\gg\omega\)), and the radius of the flux conserving rings is chosen much larger than \(a\). It is seen that \(\alpha_{s}\ll 1\) for very small values of \(B_{\omega}\). For each value of \(\eta\) there is a range of \(B_{\omega}\) in which the curve is very steep. This transition region is steeper and better defined for smaller values of \(\eta\) (i.e. large values of \(a/\delta\)). For large rotating field amplitudes, \(\alpha_{s}\) asymptotically approaches its maximum possible value, 1. It is important to note that, above a certain threshold value, only little gain can be achieved by increasing \(B_{\omega}\). Figure 2 shows the development of the field lines in a plasma cylinder of \(2\cdot 5\times 10^{-5}\)\(\Omega\) m resistivity (\(\delta=2\cdot 8\) mm and \(a/\delta=7\cdot 1\)) for rotating field amplitudes of 30 and 50 gauss. It is seen that the rotating magnetic field does not penetrate into the plasma column for \(B_{\omega}=30\) gauss whereas it does penetrate for \(B_{\omega}=50\) gauss. Note that the


<!-- Page 7 -->

[MISSING_PAGE_POST]


<!-- BLOCK FAILED: 1981_Hugrass_page007_block2 (see markdown_new\1981_Hugrass\1981_Hugrass_page007_block_fail_t36-b71_l0-r459.pdf) -->



<!-- BLOCK FAILED: 1981_Hugrass_page007_block3 (see markdown_new\1981_Hugrass\1981_Hugrass_page007_block_fail_t71-b107_l0-r459.pdf) -->



<!-- BLOCK FAILED: 1981_Hugrass_page007_block4 (see markdown_new\1981_Hugrass\1981_Hugrass_page007_block_fail_t107-b142_l0-r459.pdf) -->


threshold value is about 35 gauss in this case (see figure 1). The effect of anisotropy of the plasma resistivity was examined by studying the case \(\eta_{\perp}=0\cdot 5\,\eta_{\perp}\). It is seen from figure 3 that the general behaviour of the system is not changed. The threshold value of \(B_{u}\) is larger than the isotropic case of the same \(\eta_{\perp}\). This behaviour is expected, since the expression for the effective skin depth in this case is given by (Hugrass 1979)

\[\delta^{*}=\frac{\omega_{\rm ee}}{v_{ei_{\perp}}}\left(\frac{\eta_{\perp}}{ \mu_{0}\omega}\right)^{\frac{1}{2}} \tag{32}\]

assuming that \(B_{z}\gg B_{r}B_{\theta}\).


<!-- Page 8 -->


<!-- BLOCK FAILED: 1981_Hugrass_page008_block1 (see markdown_new\1981_Hugrass\1981_Hugrass_page008_block_fail_t0-b62_l0-r459.pdf) -->



<!-- BLOCK FAILED: 1981_Hugrass_page008_block2 (see markdown_new\1981_Hugrass\1981_Hugrass_page008_block_fail_t62-b89_l0-r459.pdf) -->


* [16]


<!-- BLOCK FAILED: 1981_Hugrass_page008_block4 (see markdown_new\1981_Hugrass\1981_Hugrass_page008_block_fail_t122-b155_l0-r459.pdf) -->


Figure 4: The transient time \(\tau\) as a function of the plasma conductivity \(\sigma=1/\eta\). The values of \(B_{\omega}\) (in gauss) are shown against the curves.

Figure 6: Radial distribution of the axial magnetic field, \(B_{z}\) for two values of the radius of the flux conserving rings, \(b\), shown against the curves. \(B_{{}_{bg}}=80\) gauss, \(\eta=5\times 10^{-6}\)\(\Omega\) m.


<!-- Page 9 -->

### The transient time

It is obvious from figure 2 that the rotating field does not penetrate instantaneously into the plasma. Figure 4 shows the transient time, \(\tau\), plotted against the conductivity for various values of the rotating field amplitude, \(B_{\omega}\). A scalar conductivity was assumed in these calculations and the radius of the conducting rings was chosen much larger than \(a\). As expected, the transient time is larger

for larger values of the conductivity. It appears that for sufficiently large rotating field amplitude, the transient time is approximately given by

\[\tau\approx\tfrac{1}{8}\mu_{0}a^{2}\sigma. \tag{33}\]

The transient time can be much larger than that predicted by (33) if the rotating field magnitude is not well above the threshold value corresponding to the

assumed plasma conductivity (see figure 1). On the other hand, smaller values of \(\tau\) can be achieved in the presence of flux conserving rings (figure 5). This is explained by the fact that less axial magnetic flux needs to be diffused into the plasma for smaller values of \(b\) (see figure 6).

## 5 Discussion

The generation of steady-state currents in confined plasmas would be of great importance to the development of fusion reactors. In a number of interesting confinement configurations, experiments with low-temperature plasmas(Hugrass, Jones & Phillips 1979; Hugrass _et al._ 1980) have shown that significant plasma currents can be driven by using an externally applied rotating magnetic field. An understanding of the physical processes involved in these experiments is essential in order to study the scaling of these schemes to the reactor regimes of interest.

In this paper we have studied a very simple physical model in which only electron motion is considered. In the absence of the Hall terms in (4), the external field would penetrate only a classical skin depth and currents would be generated only in the colder edge region of the plasma column. When the Hall terms are included the model equations become nonlinear and these must be solved numerically. The principal results of these numerical calculations show that, at sufficiently large applied magnetic field amplitudes, the rotating field does penetrate into the column considerably further than the classical skin depth, \(\delta=(2\eta/\mu_{0}\omega)^{\frac{1}{4}}\). This is caused by the downward Doppler shift in the RF field frequency in a frame of reference moving with the electron fluid. It has been shown (Hugrass 1979), and is further evidenced by these calculations, that the effective skin depth is given by

\[\delta*\simeq\frac{\omega_{\varepsilon\varepsilon}}{\nu_{\varepsilon \varepsilon}}(\eta/\mu_{0}\omega)^{\frac{1}{4}},\]

and thus substantial azimuthal current can be produced in a plasma cylinder of radius \(a\), provided \(a<\delta*\).

Although these calculations have identified the nonlinear physical mechanism


<!-- Page 10 -->

responsible for the successful penetration of the externally applied field, additional work is required to confidently scale the results to plasmas of thermonuclear interests. The most obvious extension of the fluid model considered here would be to include the motion of the ions, which is essential for simulation times longer than an ion sound speed. Since the enhanced current penetration shown by these calculations occurs on a time-scale longer than the ion sound speed, this improvement to the model remains an important area for future work.

We gratefully acknowledge the many useful and stimulating discussions with Drs I. R. Jones, R. G. Storer and Professor H. Blevin which made this work possible. One of us (R. C. G.) wishes to express his gratitude to Professor M. H. Brennan, Dr I. R. Jones and Dr R. G. Storer for their hospitality and kindness shown during the period of sabbatical leave from Princeton University, when part of this work was completed. One of the authors (W. N. H.) was the recipient of a Flanders University Research Scholarship during the course of this work.

## References

* [Blevin & Thonemann1962]Blevin, H. A. & Thonemann, P. C. 1962 _Nucl. Fusion Suppl._ part I, p. 55.
* [Braginskii1965]Braginskii, S. I. 1965 _Reviews of Plasma Physics_, vol. 1 (ed. M. A. Leontovich). Consultants Bureau.
* [Hugrass1979]Hugrass, W. N. 1979 Ph.D. Thesis, Flinders University.
* [Hugrass, Jones & Phillips1979]Hugrass, W. N., Jones, I. R. & Phillips, M. G. R. 1979 _Nucl. Fusion_, **19**, 1546.
* [Hugrass, Jones, McKenna, Phillips, Storer & Stoyer1980]Hugrass, W. N., Jones, I. R., McKenna, K. F., Phillips, M. G. R., Storer, R. G. & Stoyer, H. 1980 _Phys. Rev. Lett._**44**, 1676.
* [Jones & Hugrass1981]Jones, I. R. & Hugrass, W. N. 1981 _J. Plasma Phys._**26**, 441.
* [Richtmyer Morton1967]Richtmyer, R. D. & Morton, K. W. 1967 _Difference Method for Initial Value Problems_, 2nd ed. Interscience.