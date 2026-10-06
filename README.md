# ARRAKIS

This project is licensed under the GNU GPL v3. See the LICENSE file for more information.

## Description:

The Arrakis solver is dedicated to the simulation of unsteady trans-critical flow with non-cohesive sediment transport in rivers. It is based on finite Volume methods for the Saint-Venant-Exner equations in 1D channel with variable rectangular cross sections. The main equations are as follows

$$ \partial_t h + \partial_x hu = - hu \dfrac{B'}{B}, $$
$$ \partial_t hu + \partial_x \left( hu^2 + \frac{gh^2}{2} \right) = - hu^2 \dfrac{B'}{B} + \partial_x \mathcal{D}_u -gh\partial_x z_b -gh J, $$
$$ \partial_t z_b + \xi \partial_x q_s = - \xi q_s \dfrac{B'}{B} + \partial_x \mathcal{D}_z - \dfrac{\xi}{\rho_s} S^{ED}, $$

where:

- $h$ is the water depth,
- $u$ is the water velocity,
- $z_b$ is the bottom elevation,
- $B$ is the channel width,
- $J=C_f u|u|/(2g R_h)$ is the friction slope,
- $C_f$ is the friction coefficient, determined by empirical formulas, 
- $R_h=hB/(B+2h)$ is the hydraulic radius,
- $q_s$ is the bedload solid flux, determined by empirical formulas,
- $\mathcal{D}_u$ and $\mathcal{D}_z$ account for the diffusion of velocity, and the diffusion of the mobile bed due to slope effect,
- $\mathcal{D}_u = h \nu_u \partial_x u$, where $\nu_u$ is the velocity diffusivity
in $m^2.s^{-1}$,
- $\mathcal{D}_z = \nu_z \left( |\partial_x z_b|-J_c \right)$, if $|\partial_x z_b| \geq J_c$ and $0$ else, with $\nu_z$ the slope diffusion coefficient,
- $S^{ED}$ is the erosion-deposition source term.

Suspended sediments are modeled via the following advection-diffusion equation:

$$ \partial_t hC + \partial_x huC = - hu C \dfrac{B'}{B} + \partial_x ( h\nu_c \partial_x C) + S^{ED}. $$

where:

- $S^{ED} = \alpha_s w_s (C_{eq} -C)$ in which
- $\alpha_s$ is a calibration coefficient,
- $w_s$ is the settling velocity,
- and $C_{eq}$ is the equilibrium concentration which is detemined by empirical formulas.

Bedload and suspended sediments can be activated separately or simultaneously.
If bedload is deactivated, the term $\partial_x q_s$ is neglected and the standard Saint-Venant (SV) equations are solved.
If bedload is activated, the coupled Saint-Venant-Exner (SVE) equations are solved as a fully coupled hyperbolic system where wave speeds are influenced by sediment transport processes.
Various formulas are implemented in Arrakis for the bedload solid flux $q_s$ and for the suspended equilibrium concentration $C_{eq}$ which can be set via the keywords BEDLOAD FORMULA and SUSPENSION FORMULA.

Various numerical schemes are implemented in Arrakis to solve the SV equations (Roe, HLL, ACU, Kinetic)
and the coupled SVE equations (Roe-SVE, HLL-SVE, ACU-SVE). Second order scheme is available for both SV and SVE equations. 
Second order in time is achieved by a modified Heun scheme and second order in space by linear piecewise state reconstructions and slope limiters following the classical MUSCL technique.
For tracer transport Arrakis also includes fifth order WENO reconstructions.

For more details on input parameters, please refer to the files

- **scripts/arrakis/arrakis.dico**: in which all available options are listed and,
- **scripts/arrakis/arrakis.yml**: which is an example of input parameter file.


## Quick installation:

* install gfortran and python
* install python prerequisites:

```
pip install -r requirements.txt
```

* for an installation with OpenTELEMAC, add the following lines in your ~/.bashrc (replace <path> by the path to your telemac-mascaret installation repertory):

```
export ARRAKIS_ROOT='<path>/telemac-mascaret/optionals/addons/arrakis'
export PYTHONPATH=$PYTHONPATH:<path>/telemac-mascaret/optionals/addons/arrakis/scripts
export PATH=$PATH:<path>/telemac-mascaret/optionals/addons/arrakis/scripts/arrakis
```

* for a standard installation, add the following lines in your ~/.bashrc:

```
export ARRAKIS_ROOT='<path>/arrakis'
export PYTHONPATH=$PYTHONPATH:<path>/arrakis/scripts
export PATH=$PATH:<path>/arrakis/scripts/arrakis
```

* compile the code:
```
mkdir build
cd build
cmake ..
make
```

## How to use:

* Run using a .yml input file

```
arrakis.py <my_case>.yml
```

* Use the quick plot option of arrakis 

```
plot_arrakis.py <my_case>.yml
```

* add -b to plot mobile sediment bed if bedload is activated

```
plot_arrakis.py <my_case>.yml -b
```

## Compilation options:

* compile the code in release mode : 
```
make clean
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

* compile the code in debug mode : 
```
make clean
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```

* compile the code in optim mode : 
```
make clean
cmake -DCMAKE_BUILD_TYPE=Optim ..
make
```

* compile the code in profiling mode : 
```
make clean
cmake -DCMAKE_BUILD_TYPE=Profiling ..
make
```

## Run the valdiation:

* Run full validation

```
validate_arrakis.py
```

## Integration checklist for developpers:

* Compile in normal and debug mode
* Clean the code until no warnings
* Run full validation in normal and debug mode
* If validation succeed, make a merge request

