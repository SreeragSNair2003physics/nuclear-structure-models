# Nuclear Structure Models: From SEMF to Nilsson

This repository is a computational companion to my MSc project on nuclear structure. It follows a progression from the global properties described by the semi-empirical mass formula (SEMF) to single-particle shell levels and their evolution with nuclear deformation.

## 1. Binding energy: the SEMF

The SEMF estimates a nucleus’s total binding energy from its mass number \(A\), proton number \(Z\), and neutron number \(N=A-Z\). Its terms approximate the bulk attraction of the nuclear force, the surface correction, proton repulsion, neutron–proton imbalance, and pairing.

$$
B(A,Z)=a_vA-a_sA^{2/3}
-a_c\frac{Z(Z-1)}{A^{1/3}}
-a_a\frac{(A-2Z)^2}{A}
+\delta_p(A,Z)
$$

The pairing term used in the scripts is

$$
\delta_p =
\begin{cases}
+a_pA^{-3/4}, & Z,N\ \text{both even},\\
-a_pA^{-3/4}, & Z,N\ \text{both odd},\\
0, & A\ \text{odd}.
\end{cases}
$$

The calculation divides \(B\) by \(A\) to obtain binding energy per nucleon. For the data curve, it selects the largest CSV value of \(B/A\) at each \(A\), forming an upper envelope. For the SEMF curve, it searches over \(Z\) at each \(A\) and plots the largest predicted \(B/A\); this model curve starts at \(A=16\).

![Nuclear binding energy curve](figures/binding_energy_curve.png)

## 2. Beta stability and drip lines

For each fixed \(A\), the code searches over \(Z\) and selects the composition with the largest SEMF binding energy:

$$
Z_\beta(A)=\underset{1\leq Z<A}{\operatorname{arg\,max}}\;
B\bigl(Z,A-Z\bigr).
$$

The neutron and proton drip-line estimates use one-nucleon separation energies. In this code, a positive separation energy means that adding the nucleon increases the total binding energy:

$$
S_n(Z,N+1)=B(Z,N+1)-B(Z,N),
$$

$$
S_p(Z+1,N)=B(Z+1,N)-B(Z,N).
$$

The plotted boundary is where the corresponding separation energy first becomes non-positive after being positive. These are SEMF estimates; the formula does not include the detailed shell effects that influence real drip lines.

![SEMF beta stability and drip lines](figures/semf_stability_driplines.png)

## 3. Spherical 3D harmonic oscillator

The shell model describes nucleons occupying quantized single-particle orbitals. The 3D harmonic oscillator is a simple spherical potential for introducing those shells. Its energy depends only on the major-shell number \(N\):

$$
E_N=\left(N+\frac{3}{2}\right)\hbar\omega.
$$

Including the two spin states, the degeneracy of shell \(N\) is

$$
g_N=(N+1)(N+2).
$$

The script plots these shell energies and prints the cumulative particle counts after filling each shell. These are closures of the ideal harmonic oscillator; the model alone does not reproduce the full set of observed nuclear magic numbers.

![3D harmonic oscillator energy levels](figures/3d_ho_energy_levels.png)

## 4. Spin-orbit splitting in a shell-model diagram

The shell-model diagram adds an orbital correction and a spin-orbit term to the oscillator energy. The script uses the schematic expression

$$
E_{N\ell j}
=
N+\frac{3}{2}
+D\ell(\ell+1)
+C\langle\boldsymbol{\ell}\cdot\mathbf{s}\rangle,
$$

with chosen values \(C=-0.1\) and \(D=0.0225\), in oscillator-energy units. For spin \(s=\tfrac12\),

$$
\langle\boldsymbol{\ell}\cdot\mathbf{s}\rangle
=
\begin{cases}
\ell/2, & j=\ell+\tfrac12,\\
-(\ell+1)/2, & j=\ell-\tfrac12.
\end{cases}
$$

The spin-orbit term splits the two \(j\) values when \(\ell>0\); an \(s\) orbital does not split. This figure is a schematic illustration using the selected \(C,D\) values, rather than a parameter-fitted spectrum.

![Shell-model energy levels](figures/shell_model_energy_levels.png)

## 5. First-order axial deformation of the oscillator

An axial deformation distinguishes motion along the symmetry axis \(z\) from motion in the perpendicular plane. In the first-order approximation used here,

$$
\frac{E}{\hbar\omega_0}
=
N+\frac{3}{2}
-\frac{\delta}{3}(2n_z-n_\perp),
\qquad
n_\perp=n_x+n_y=N-n_z.
$$

The deformation parameter \(\delta\) shifts the spherical oscillator energies and lifts degeneracies. This calculation illustrates the deformation effect by itself; it does not yet include spin-orbit or Nilsson \(\ell^2\) terms.

![First-order axial deformation splitting of 3D HO levels](figures/anisotropic_ho_first_order.png)

## 6. Nilsson single-particle levels

The Nilsson model extends the spherical shell-model picture to deformed nuclei. The basis states in the script are written as \(|N,\ell,\Lambda,\Sigma\rangle\), where \(\Lambda\) is the orbital angular-momentum projection on the symmetry axis, \(\Sigma\) is the spin projection, and

$$
\Omega=\Lambda+\Sigma.
$$

The Hamiltonian implemented is

$$
H_{\rm Nilsson}
=
H_{\rm osc}
+h_\delta
-2\kappa_N\hbar\omega_0\,
\boldsymbol{\ell}\cdot\mathbf{s}
-\kappa_N\mu_N\hbar\omega_0
\left[
\ell^2-\frac{N(N+3)}{2}
\right],
$$

with the deformation term

$$
h_\delta
=
-\delta\hbar\omega_0
\frac{4}{3}\sqrt{\frac{\pi}{5}}\,
r^2Y_{20}.
$$

Here \(\kappa_N\) and \(\mu_N\) are the shell-dependent parameters listed in the script. The program constructs the Hamiltonian matrix in the stated basis, groups states by \(N\) and \(\Omega\), diagonalizes each block, and plots the resulting energies as \(\delta\) varies.

**Approximation:** the matrix elements of the deformation term are restricted to states with the same major-shell number \(N\). Couplings between different major shells are omitted.

![Neutron Nilsson diagram](figures/nilsson_diagram.png)

## 7. Shell gaps versus deformation

At each deformation, the code collects the calculated single-particle energies, counts each positive-\(\Omega\) level twice for its time-reversed partner, and sorts the levels. The gap after filling \(p\) states is

$$
\Delta_p(\delta)=E_{p+1}(\delta)-E_p(\delta).
$$

In this plot, \(p\) is the number of filled neutron states; it does not mean proton number. The plotted gaps are properties of this truncated model spectrum, not measured nuclear gaps.

![Neutron shell gaps versus deformation](figures/neutron_shell_gaps_vs_deformation.png)

## Data and references

The binding-energy CSV is described in the script as based on AME2020. The AME2020 mass tables and evaluation are available from the [Atomic Mass Data Center](https://amdc.impcas.ac.cn/web/masseval.html); the evaluation is described in [AME2020 Part II](https://doi.org/10.1088/1674-1137/abddaf).

The exact source for the shell-dependent \(\kappa_N,\mu_N\) values and the \(C,D\) values used in the shell-model diagram should also be cited here.
