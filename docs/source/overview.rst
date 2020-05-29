Overview
==========

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
specified constants and the perturbation density :math:`\rho'` is assumed to **either**
have the following smooth form

.. math::
    \rho'(r, \varphi) &= \frac{r^k}{R_+^k} \cos(n\varphi), & \textbf{(2D-smooth)} \\
    \rho'(r, \theta, \varphi) &= \frac{r^k}{R_+^k} Y_{lm}(\theta, \varphi), & \textbf{(3D-smooth)}

where in 2D, we use cylindrical coordinates with radius :math:`r` and azimuthal angle :math:`\varphi`,
and in 3D, spherical coordinates with radius :math:`r`, co-latitude :math:`\theta`,
and longitude :math:`\varphi`. The radial dependency is a simple polynomial (monomial) 
of order :math:`k`. In 2D, :math:`n` is the wave number and in 3D :math:`l` and :math:`m` are the
degree and order of the spherical harmonic function :math:`Y_{lm}` (see :func:`assess.Y` for definition).

**Or**, :math:`\rho'` is a perturbation at a specified radius :math:`r'`

.. math::
    \rho'(r, \varphi) &= \delta(r-r') \cos(n\varphi), & \textbf{(2D-delta)} \\
    \rho'(r, \theta, \varphi) &= \delta(r-r') Y_{lm}(\theta, \varphi), & \textbf{(3D-delta)}

where :math:`\delta(r-r')` is the Dirac delta function. Combined with two types of 
boundary conditions

.. math::
   -\mathbf{n}\cdot\tau\cdot\mathbf{n} + p &= 0, \mathbf{n}\cdot\mathbf{u}=0, & \text{at }r&=R_- \text{and }r=R_+ & \textbf{(free-slip)} \\
   \mathbf{u} &= 0, & \text{at }r&=R_- \text{and }r=R_+ & \textbf{(zero-slip)}

this leads to eight analytical solutions, which are implemented in the following classes

   * :class:`assess.CylindricalStokesSolutionSmoothFreeSlip`
   * :class:`assess.CylindricalStokesSolutionSmoothZeroSlip`
   * :class:`assess.CylindricalStokesSolutionDeltaFreeSlip`
   * :class:`assess.CylindricalStokesSolutionDeltaZeroSlip`
   * :class:`assess.SphericalStokesSolutionSmoothFreeSlip`
   * :class:`assess.SphericalStokesSolutionSmoothZeroSlip`
   * :class:`assess.SphericalStokesSolutionDeltaFreeSlip`
   * :class:`assess.SphericalStokesSolutionDeltaZeroSlip`
