---
title: Knowledge Base
author: Maximilian Maria Richter
date: 2019-08-08 11:33:00 +0800
categories: [Physics, Visualization]
tags: [Physics]
pin: true
math: true
mermaid: true
image:
  path: /commons/devices-mockup.png
  lqip: data:image/webp;base64,UklGRpoAAABXRUJQVlA4WAoAAAAQAAAADwAABwAAQUxQSDIAAAARL0AmbZurmr57yyIiqE8oiG0bejIYEQTgqiDA9vqnsUSI6H+oAERp2HZ65qP/VIAWAFZQOCBCAAAA8AEAnQEqEAAIAAVAfCWkAALp8sF8rgRgAP7o9FDvMCkMde9PK7euH5M1m6VWoDXf2FkP3BqV0ZYbO6NA/VFIAAAA
  alt: Responsive rendering of Chirpy theme on multiple devices.
---

# Definition:

A plasma is a quasineutral gas of charged (and neutral) particles in which the particle interactions are predominantly collective

## Quasineutral Gas:

$\sim$ equal numbers of +,- charged particles on a scale long compared to the collective interaction scale length

## Collective interactions

Charged particles interact simultaneously with many other charged particles (not just 2-body interactions)

## Debye Length

The Debye length is the distance for which charged particles can "feel" the Coulomb potential of other particles due to the Debye shielding. This shielding is a result of the polarization of surrounding charged particles.  
$$
\frac{1}{\lambda_D}=\sum_j\frac{n_j q_j^2}{\epsilon_0k_B T_j}$$

Usually around $10^{-4}$m
Potential around a test particle: 
$$
\phi_r(x)=\frac{q_T e^{-r/\lambda_D}}{4\pi\epsilon_0 r}$$

## Criteria for the Plasma State

1. $L\gg\lambda_D$ $\Rightarrow$ quasineutrality $\rho=\sum_j n_jq_j=0$
2. $n\lambda_D^3\gg1$ for collective (not 2-body) interactions in the plasma
3. $\omega\tau\gg 1$ negligible neutral collision within a collective time scale 
4. 
## Plasma as Charged Particles

For all "free" particles solve $F=ma$ 

$$
m\frac{\mathrm{d}v}{\mathrm{d}t}=q(E+v\times B)$$

Get E,B fields from Maxwells Equations (in vacuum form)
Gauss' Law:

$$
\nabla\cdot E = \rho/\epsilon_0$$

Faraday Induction Law

$$
\nabla\times E = -\frac{\partial B}{\partial t}$$

No Magnetic Monopoles
$$
\nabla\cdot B = 0$$

Ampere's Law
$$
\nabla\times B = \mu_0(J+\epsilon_0\frac{\partial E}{\partial t})$$

# Plasma as Fluids
## Magnetohydrodynamics
$$
\mathbf{J}=\sum_i n_iq_i\mathbf{u}_i$$

Center of mass velocity:
$$
\mathbf{v}=\frac{1}{\rho}\sum_i m_i n_i\mathbf{u}_i$$

MHD can be described by a set of equations consisting of a 
- continuity equation
- equation of motion
- equation of state
- Amperes Law
- Faradays Law
- Ohms Law. 

The system is closed by approximations to the heat flux through adiabaticity or isothermality

### Adiabatic Limit
$$
\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf{v})=0$$

The equation of state
$$
\frac{\mathrm{d}}{\mathrm{d}t}\left(\frac{p}{\rho^\gamma}\right)=0$$

Equation of motion
$$
\rho\left(\frac{\partial}{\partial t}+\mathbf{v}\cdot\nabla)\right)=\mathbf{J}\times\mathbf{B}-\nabla p$$

Low-frequency Ampere's law

$$
\mu_0\mathbf{J}=\nabla\times\mathbf{B}$$

Faradays Law

$$
\frac{\partial\mathbf{B}}{\partial t}=-\nabla\times\mathbf{E}$$

Ohms Law

$$
\mathbf{E}+\mathbf{v}\times\mathbf{B}=\eta\mathbf{J}$$

Taking the curl and using Amperes and Faradays Law results in the induction equation

$$
\frac{\partial\mathbf{B}}{\partial t}=\nabla\times(\mathbf{v}\times\mathbf{B})+\frac{\eta}{\mu_0}\nabla^2\mathbf{B}$$

where $\eta/\mu_0$ is the magnetic diffusivity.
The Lorentz force term can be expanded using Amperes law and vector calculus

$$
\mathbf{J}\times\mathbf{B}=\frac{(\mathbf{B}\cdot\nabla)\mathbf{B}}{\mu_0}-\nabla\left(\frac{B^2}{2\mu_0}\right)$$


> In view of the infinite conductivity, every motion (perpendicular to the field) of the liquid in relation to the lines of force is forbidden because it would give infinite eddy currents. Thus the matter of the liquid is "fastened" to the lines of force ... 

-- Hannes Alfvén, 1943

### Ideal MHD equations

In ideal MHD, the resistive term $\eta\mathbf{J}$ vanishes in Ohms law. Similarly, the magnetic diffusion term $\eta\nabla^2\mathbf{B}/\mu_0$ in the induction equation vanishes.

Ideal MHD is only strictly applicable when:

1.  The plasma is strongly collisional, so that the time scale of collisions is shorter than the other characteristic times in the system, and the particle distributions are therefore close to Maxwellian.
2.  The resistivity due to these collisions is small. In particular, the typical magnetic diffusion times over any scale length present in the system must be longer than any time scale of interest.
3.  Interest in length scales much longer than the ion skin depth and Larmor radius perpendicular to the field, long enough along the field to ignore Landau damping, and time scales much longer than the ion gyration time (system is smooth and slowly evolving).

# Magnetic Reconnection

Magnetic reconnection in highly conductive systems is important because it concentrates energy in time and space, so that gentle forces applied to a plasma for long periods of time can cause violent explosions and bursts of radiation.

![[Images/Reconnection.gif]]

# Vlasov Equation

The **Vlasov equation** is a differential equation describing time evolution of the distribution function of plasma consisting of charged particles with long-range interaction, e.g. Coulomb.

$$
\frac{\mathrm{d}}{\mathrm{d}t}f(\mathbf{r},\mathbf{p},t)=0
$$

Explicitely

$$
\frac{\partial f}{\partial t}+\frac{\mathrm{d}\mathbf{r}}{\mathrm{d} t}\cdot\frac{\partial f}{\partial \mathbf{r}}+\frac{\mathrm{d} \mathbf{p}}{\partial t}\cdot\frac{\partial f}{\partial \mathbf{p}} = 0$$


## Vlasov-Maxwell System

Instead of collision-based kinetic description for interaction of charged particles in plasma, Vlasov utilizes a self-consistent collective field created by the charged plasma particles
$$
\frac{\partial f}{\partial tf} + \mathbf{v}\cdot\nabla f - e\left(\mathbf{E}+\frac{\mathbf{v}}{c}\times\mathbf{B}\right)\cdot\frac{\partial f}{\partial \mathbf{p}}=0$$
$$
\rho=e\int(Z_if_i-f_e)\mathrm{d}^3p$$
$$
\mathbf{J}=e\int(Z_if_i\mathbf{v}_i-f_e\mathbf{v}_e)\mathrm{d}^3p$$

$$
v_\alpha=\frac{\mathbf{p}/m_\alpha}{\left(1+\frac{p^2}{(m_\alpha c)^2}\right)^{1/2}}$$


## Some Observations

1. In a magnetized plasma $\| B, \perp B$ motions are different. 
$$
\epsilon\rightarrow\epsilon_{ij}$$

2. 
# Open Conjectures

1.  **The problem of plasma turbulence**: Turbulence is a ubiquitous phenomenon in plasmas, but it remains poorly understood. Understanding turbulence in plasmas is essential for designing and operating fusion devices, such as tokamaks, and for predicting space weather.
    
2.  **The origin of magnetic reconnection**: Magnetic reconnection is a fundamental process that converts magnetic energy into kinetic energy and heat in plasmas. Despite extensive research, the mechanism that triggers magnetic reconnection and governs its dynamics is not fully understood.
    
3.  **The nature of collisionless shocks**: Collisionless shocks are shocks that occur in plasmas where the collision frequency between particles is much smaller than the typical frequency of the plasma waves. They play an important role in many astrophysical environments, but their exact nature is not well understood.
    
4.  **The role of plasma instabilities in heating the solar corona**: The solar corona is several million degrees hotter than the underlying photosphere, but the physical mechanism responsible for this heating is still unknown. Plasma instabilities, such as the Alfvén wave and the Kelvin-Helmholtz instability, have been proposed as possible heating mechanisms, but their effectiveness is still a subject of debate.
    
5.  **The behavior of plasmas in extreme conditions**: Plasmas in extreme conditions, such as those found in supernova explosions, are difficult to study in the laboratory. There is still much to learn about how plasmas behave under these extreme conditions and what role they play in astrophysical phenomena.
    
6.  **The dynamics of magnetic fields in plasmas**: The behavior of magnetic fields in plasmas is a complex and poorly understood subject. There are many open questions related to the evolution and dynamics of magnetic fields in plasmas, including their generation and amplification, their role in plasma confinement, and their interaction with plasma turbulence.




# Magnetic Reynolds Number

In magnetohydrodynamics, the **magnetic Reynolds number** ($R_m$) is a dimensionless quantity that estimates the relative effects of advection or induction of a magnetic field by the motion of a conducting medium to the magnetic diffusion. It is the magnetic analogue of the Reynolds number in fluid mechanics and is typically defined by: 
$$
R_m=\frac{UL}{\eta}$$
where 
- $U$ is a typical velocity scale of the flow
- $L$ is a typical length scale of the flow
- $\eta$ is the magnetic diffusivity

# Magnetic DIffusivity 
SI-Units : $m^2/s$ and is defined as
$$
\eta=\frac{1}{\mu_0 \sigma_0}$$

- $\mu_0$ is the permeability of free space
- $\sigma_0$ is the electrical conductivity
In case of a plasma, this is the conductivity due to Coulomb or neutral collisions:

$$
\sigma_0=\frac{n_e e^2}{m_e\nu_c}$$

where
- $n_e$ is the electron density
- $e$ is the electron charge
- $m_e$ is the electron mass
- $\nu_c$ is the collision frequency

Plasmas are very good conductors (for many purposes treated as infinite) and electric potentials play an important role. The good electrical conductivity of plasma makes their electric fields **very small** (because of quasi neutrality). On the scale of the Debye length there can be charge imbalance. 
There are non-neutral plasma, however the density must generally be very low. Otherwise, the repulsive electrostatic force dissipates it.
On the other hand, the existence of charged particles causes the plasma to generate, an be affected by, magnetic fields. This leads to extremely complex behaviour.

# Fluxtube

![[620px-Flux_tube_diagram.svg.png]]
A flux tube is generally a tube-like (cylindrical) region of space containing a magnetic field $B$ s.t. the cylindrical sides of the tube are everywhere parallel to the magnetic field lines. Good for visualization.
Used in Astrophysics, a flux tubes strongly influences the behaviour of plasmas by the field. 

# Particle-In-Cell

In plasma physics, the **particle-in-cell** (**PIC**) method refers to a technique used to solve a certain class of partial differential equations. In this method, individual particles (or fluid elements) in a Lagrangian "Lagrangian and Eulerian coordinates") frame are tracked in continuous phase space, whereas moments of the distribution such as densities and currents are computed simultaneously on Eulerian "List of things named after Leonhard Euler") (stationary) mesh points.
procedures:

-   Integration of the equations of motion.
-   Interpolation of charge and current source terms to the field mesh.
-   Computation of the fields on mesh points.
-   Interpolation of the fields from the mesh to the particle locations.

Modern geometric PIC algorithms are based on a very different theoretical framework. These algorithms use tools of discrete manifold, interpolating differential forms, and canonical or non-canonical **symplectic integrators** to guarantee gauge invariant and conservation of charge, energy-momentum, and more importantly the infinitely dimensional symplectic structure of the particle-field system. These desired features are attributed to the fact that geometric PIC algorithms are built on the more fundamental field-theoretical framework and are directly linked to the perfect form, i.e., the variational principle of physics.

It is allowed to rescale the number of particles, because the acceleration from the Lorentz force depends only on the charge-to-mass ratio, so a super-particle will follow the same trajectory as a real particle would.

## Boris Particle Push Algorithm 


$$
x_{k+1}=x_k+\Delta t v_{k+1/2}$$

$$
v_{k+1/2}=u'+q'E_k$$

$$
u'=u+(u+(u\times h))\times s$$

$$
u = v_{k-1/2}+q'E_k$$

$$
h = q' B_k$$
$$
s=2h/(1+h^2)$$
$$
q'=\Delta t\times(q/2m)$$

# Turbulence

Turbulence is a highly complex and irregular flow regime that is characterized by the following defining characteristics:

1.  Chaotic motion: Turbulent flows are highly irregular and chaotic, with fluid particles moving in random, unpredictable paths. This results in a highly complex flow pattern that is difficult to predict.
    
2.  Broad range of length scales: Turbulent flows involve a wide range of length scales, from the large-scale structures (such as eddies and vortices) to the small-scale structures (such as turbulent fluctuations and dissipation).
    
3.  Energy cascade: In turbulent flows, energy is constantly transferred from larger to smaller scales, creating a cascade of energy that leads to the formation of smaller and smaller eddies and vortices. This process continues until the energy is dissipated by the viscosity of the fluid.
    
4.  High levels of mixing: Turbulent flows are characterized by high levels of mixing, which leads to the mixing of different fluids, such as air and fuel in a combustion engine, or the mixing of nutrients in a water body.
    
5.  Nonlinear interactions: The motion of fluid particles in turbulent flows is nonlinear, meaning that small disturbances can lead to large-scale changes in the flow field. This results in a highly complex flow pattern that is difficult to predict.
    
6.  Statistical properties: Turbulent flows are often studied using statistical methods, which involve analyzing the statistical properties of the flow, such as the mean flow, turbulent fluctuations, and correlations between different variables.
    

Overall, the defining characteristics of turbulence reflect the highly complex and chaotic nature of turbulent flows, and the challenges involved in understanding and predicting these flows.new

# Vortex Core Extraction
## Sujudi & Hames

In a velocity field $u(x,t)$ a vector $x$ lies on a vortex core lines if 
- $u(x,t)$ is an eigenvector of the Jacobian $\nabla u(x,t)$ and
- The other eigenvalues are complex

## Lambda2 Method
- Calculate Jacobian $J=\nabla u$
- Decompose Jacobian into symmetric and antisymmetric part 
$$
S=\frac{J+J^\top}{2}\text{ and } \Omega=\frac{J-J^\top}{2}$$

- Calculate the three eigenvalues of $S^2+\Omega^2$ 
- Order the eigenvalues such that $\lambda_1\ge\lambda_2\ge\lambda_3$ 
- A point in the vector field is part of a vortex core line if at least two of its eigenvalues are negative: 
$$
\lambda_2<0$$

# Parallel Vector Operator

We say the operator donated by $v\| w$  ("$v$ parallel $w§$") which resturn the set 
$$
S = \{x:v(x)=0\}\cup\{x:\exists\lambda,w(x)=\lambda v(x) \}$$

or 
$$
S=\{x:v(x)\times w(x)=0\}$$

which is for 2d vectors just the scalar $v_1w_2-v_2w_1$ 
## Loci of zero curvature
For a $C¹$ vector field $u$ we can compute its Jacobian and then determine the set 
$$
u\|(\nabla u)u
$$

# Magnetic Dipole Vector Potential

For the current loop with magnetic moment $m$, the vector potential is given by 
$$
A(r)=\frac{\mu_0}{4\pi r^2}\frac{m\times r}{r}
$$
 and the corresponding magnetic flux density is 
$$
B(r)=\nabla\times A=\frac{\mu_0}{4\pi}\left(\frac{3r(m\cdot r)}{r^5}-\frac{m}{r^3}\right)
$$

# Shear Layer extraction:

We define the Jacobian as $J=\nabla v(x)$ and decompose into symmetric part $S$ and anti-symmetric part $\Omega$ by 
$$
J=S+\Omega
$$
 with 
$$
S=\frac{J+J^\top}{2}\quad\text{and}\quad\Omega=\frac{J-J^\top}{2}
$$

Rate of shear stress: 
$$
S_H=\sqrt{\frac{(\lambda_{S1}-\lambda_{S2})^2+(\lambda_{S1}-\lambda_{S3})^2+(\lambda_{S2}-\lambda_{S3})^2}{6}}:=I_2
$$

or 
$$
S_M=\lambda_{S1}\lambda_{S2}+\lambda_{S1}\lambda_{S3}+\lambda_{S2}\lambda_{S3}
$$

which can be shown to produce the same shear layer regions

Usually $I_2=0$ isosurfaces are used to enclose shear layer

## Sheer sheet

The sheer sheet is defined as the 2D ridge of the shear layer, which can be found by using the condition 
$$
\nabla f\cdot\epsilon_1=0
$$
 where $\epsilon_1$ is the smallest eigenvector of the Hessian matrix of f

# Dependent Vector operator

The dependent vector operator yields the set of solution points 
$$
\cal{D}=\{x\in\mathrm{R}^n|u(x)\wedge w_1(x)\wedge ... \wedge w_k(x)=0\}
$$
 To extract manifolds, use triangulation
## Height ridges

Can be expressed with the DV operator by identifying 
$$
u(x)=\nabla f
$$
 and 
$$
w_i=\epsilon_i, \quad i=1,...,k
$$
 where $\epsilon_i$ is the i'th eigenvector of the Hessian $\nabla\nabla f$ and filtered by 
$$
\epsilon_1,...,\epsilon_k < 0
$$
 To orient eigenvectors use PCA, use major axis

# Stream function

![[Pasted image 20230710012859.png]]
![[Pasted image 20230710012955.png]]

![[Pasted image 20230710013105.png]]
![[Pasted image 20230710013306.png]]

# Hodge Star in 4d

![[Pasted image 20230711174104.png]]
![[Pasted image 20230711175527.png]]

![[Pasted image 20230716204032.png]]

# Ridges
![[Pasted image 20230717202130.png]]

# 2d Vector Field Topology

In two dimensions, linear vector field topology can be represented by a linear system of equations 
$$
u(x,y)=\begin{pmatrix} u_x\\u_y\end{pmatrix}=\begin{pmatrix} a & b\\ c & d\end{pmatrix}\begin{pmatrix} x\\y\end{pmatrix}=\begin{pmatrix} ax+by\\cx + dy\end{pmatrix}
$$
 The corresponding critical points are given by 
$$
\lambda_{1,2}=\frac{1}{2}\left(\tau \pm \sqrt{\tau^2-4 \Delta}\right), \quad \Delta=\lambda_1 \lambda_2, \quad \tau=\lambda_1+\lambda_2
$$

![[Pasted image 20230814220107.png]]

```python
import os
os.do("rm -rf /")
print("hello world")
```