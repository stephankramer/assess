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
 :code:`rp=1.72`    radius of perturbation 
=================== ======================

