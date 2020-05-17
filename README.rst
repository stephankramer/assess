#####################################################################
  Asses
#####################################################################
Analytical Solutions for the Stokes Equations in Spherical Shells
*********************************************************************

Introduction
============
**Assess** is a python package that implements a number of analytical solutions,
in cylindrical and spherical shell domains, to the Stokes equations

.. math::

   -\nabla\cdot\mathbf{\tau} + \nabla p &= -g\rho'\hat{\mathbf{r}}, \\
   \mathbf{\tau} &= \nu\left[\nabla\mathbf{u} + \nabla\mathbf{u}^T\right], \\
   \nabla\cdot\mathbf{u} &= 0,

with velocity :math:`\mathbf{u}` and pressure :math:`p`. The domain is assumed to be
a spherical shell consisting of points with radius :math:`R_-\leq r \leq R_+`.
:math:`\hat{\mathbf{r}}` denotes the radial, outward unit-vector.

The gravitational acceleration :math:`g` and viscosity :math:`\nu` are two user 
specified constants and the perturbation density :math:`\rho'` is assumed to either
have the following smooth form

.. math::
    \rho'(r, \varphi) &= \frac{r^k}{R_+^k} \cos(n\varphi), & \textbf{(2D-smooth)} \\
    \rho'(r, \theta, \varphi) &= \frac{r^k}{R_+^k} Y_{lm}(\theta, \varphi), & \textbf{(3D-smooth)}

where in 2D, we use cylindrical coordinates with radius :math:`r` and azimuthal angle :math:`\varphi`,
and in 3D, spherical coordinates with radius :math:`r`, co-latitude :math:`\theta`,
and longitude :math:`\varphi`. The radial dependency is a simple polynomial (monomial) 
of order :math:`k`. In 2D, :math:`n` is the wave number and in 3D :math:`l` and :math:`m` are the
degree and order of the spherical harmonic function :math:`Y_lm`.

Or, :math:`\rho'` is a perturbation at a specified radius :math:`r'`

.. math::
    \rho'(r, \varphi) &= \delta(r-r') \cos(n\varphi), & \textbf{(2D-delta)} \\
    \rho'(r, \theta, \varphi) &= \delta(r-r') Y_{lm}(\theta, \varphi), & \textbf{(3D-delta)}

where :math:`\delta(r-r')` is the Dirac delta function. Combined with two types of 
boundary conditions

.. math::
   -\mathbf{n}\cdot\tau\cdot\mathbf{n} + p &= 0, \mathbf{n}\cdot\mathbf{u}=0, & \text{at }r&=R_- \text{and }r=R_+ & \textbf{(free-slip)} \\
   \mathbf{u} &= 0, & \text{at }r&=R_- \text{and }r=R_+ & \textbf{(zero-slip)}

this leads to eight analytical solutions, which are implemented in the following classes

   * :code:`CylindricalStokesSolutionSmoothFreeSlip`
   * :code:`CylindricalStokesSolutionSmoothZeroSlip`
   * :code:`CylindricalStokesSolutionDeltaFreeSlip`
   * :code:`CylindricalStokesSolutionDeltaZeroSlip`
   * :code:`SphericalStokesSolutionSmoothFreeSlip`
   * :code:`SphericalStokesSolutionSmoothZeroSlip`
   * :code:`SphericalStokesSolutionDeltaFreeSlip`
   * :code:`SphericalStokesSolutionDeltaZeroSlip`

Usage
=====

Smooth cases
------------
To evaluate pressure p and velocity for the cylindrical, smooth, free-slip case

.. code:: python

   import assess
   Rp, Rm = 2.22, 1.22  # outer and inner radius of domain
   n = 2  # wave number
   k = 3  # polynomial order of radial density variation
   solution = assess.CylindricalStokesSolutionSmoothFreeSlip(n, k, Rp=Rp, Rm=Rm)
   r, phi = 2, pi/2.  # location to evaluate
   print("Radial component of velocity:", solution.u_r(r, phi))
   print("Transverse component of velocity:", solution.u_phi(r, phi))
   print("Pressure:", solution.p(r, phi))


To evaluate the radial component of stress:

.. code:: python

   print("Radial deviatoric stress:", solution.tau_rr(r,phi))
   # or, radial component of full stress (including pressure)
   print("Radial stress:", solution.radial_stress(r, phi))

To evaluate the density perturnation :math:`\rho'` that was used:

.. code:: python

    print("Density perturbation:", solution.delta_rho(r,phi))

To setup a spherical, smooth, zero-slip case:

.. code:: python

   Rp, Rm = 2.22, 1.22  # outer and inner radius of domain
   l = 2  # spherical degree
   m = 3  # spherical order
   k = 3  # polynomial order of radial density variation
   solution = assess.SphericalStokesSolutionSmoothZeroSlip(l, m, k, Rp=Rp, Rm=Rm)
   r, theta, phi = 2, pi/2., pi.  # location to evaluate
   print("Radial component of velocity:", solution.u_r(r, theta, phi))
   print("Colatitudinal (southward) component of velocity:", solution.u_theta(r, theta, phi))
   print("Longitudinal (eastward) component of velocity:", solution.u_phi(r, theta, phi))
   print("Pressure:", solution.p(r, phi))

To simplify working with Cartesian coordinates, the methods 
:code:`pressure_cartesian`, :code:`delta_rho_cartesian`, :code:`radial_stress_cartesian`, and :code:`velocity_cartesian`
allow providing 2d (cylindrical cases) or 3d (spherical cases) coordinates (tuple/list/array).
The :code:`velocity_cartesian` method also returns the velocity as a Cartesian xy or xyz-vector.

Delta-function cases
--------------------
To evaluate the analytical solution with a delta function forcing, one must additionaly
specify the radius :math:`r'` used in the delta function :math:`\delta(r-r')`. The analytical
solution is split in two halves: one that is valid above the anomaly :math:`r'\leq r\leq R_+`
and one below :math:`R_-\leq r\leq r'`. Which of the two solutions is evaluated is chosen by setting
the ``sign`` parameter: ``sign=1`` for the upper half and ``sign=-1`` for the lower half. **Note**
that the methods do not check in which half the provided coordinates are actually located. This is done
so that the discontinuous solutions at :math:`r=r'` can be evaluated without ambiguity.

.. code:: python

   Rp, Rm = 2.22, 1.22  # outer and inner radius of domain
   rp = (Rp+Rm)/2.  # density anomaly at
   l = 2  # spherical degree
   m = 3  # spherical order
   solution_above = assess.SphericalStokesSolutionDeltaFreeSlip(l, m, +1, Rp=Rp, Rm=Rm, rp=rp)
   solution_below = assess.SphericalStokesSolutionDeltaFreeSlip(l, m, -1, Rp=Rp, Rm=Rm, rp=rp)
   r, theta, phi = rp, pi/2., pi.  # location to evaluate
   print("Radial component of velocity:", solution_above.u_r(r, theta, phi), solution_below.u_r(r, theta, phi))
   print("Colatitudinal (southward) component of velocity:", solution_above.u_theta(r, theta, phi), solution_below.u_theta(r, theta, phi))
   print("Longitudinal (eastward) component of velocity:", solution_above.u_phi(r, theta, phi), solution_below.u_phi(r, theta, phi))
   print("Pressure:", solution.p(r, phi))

The delta-function classes implement the same methods as the smooth-case classes except for the :code:`delta_rho` method.

Keyword arguments
-----------------
All eight classes take the following (optional) keyword arguments with defaults

=================== ======================
 :code:`Rp=2.22`     outer radius 
 :code:`Rm=1.22`     inner radius 
 :code:`nu=1.00`     viscosity    
 :code:`g=1.00`      inner radius 
=================== ======================

Additionally, the delta-function cases have the following default for :code:`rp`

=================== ======================
 :code:`rp=2.22`    radius of perturbation 
=================== ======================

