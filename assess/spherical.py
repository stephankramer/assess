"""Analytical solutions in spherical shell domains."""
from __future__ import division
from assess.cylindrical import AnalyticalStokesSolution
from assess.smooth import coefficients_sphere_smooth_fs, coefficients_sphere_smooth_ns
from assess.delta import coefficients_sphere_delta_fs, coefficients_sphere_delta_ns
import scipy
import numpy
from math import sqrt, atan2, cos, sin, tan, acos, pi

__all__ = ['SphericalStokesSolutionSmoothFreeSlip', 'SphericalStokesSolutionSmoothZeroSlip',
           'SphericalStokesSolutionDeltaFreeSlip', 'SphericalStokesSolutionDeltaZeroSlip']


# Legendre polynomials and their derivatives:

def Y(m, l, theta, phi):
    """Real part of spherical harmonic function Y^m_l.

    :param m: order of the harmonic
    :param l: degree of the harmonic
    :param theta: co-latitude in [0, pi]
    :param phi: longitude in [0, 2*pi]
    """
    # everywhere we take the real part of Y, corresponding to the cos(m phi) part of the solution
    # NOTE: theta and phi are swapped in scipy wrt our definitions
    return scipy.special.sph_harm(m, l, phi, theta).real


def dYdphi(m, l, theta, phi):
    """Colatitudinal derivative of spherical harmonic function Y^m_l.

    :param m: order of the harmonic
    :param l: degree of the harmonic
    :param theta: co-latitude in [0, pi]
    :param phi: longitude in [0, 2*pi]
    """
    # except in theta derivatives (of odd order)
    return -m * scipy.special.sph_harm(m, l, phi, theta).imag


def dYdtheta(m, l, theta, phi):
    """Longitudinal derivative of spherical harmonic function Y^m_l.

    :param m: order of the harmonic
    :param l: degree of the harmonic
    :param theta: co-latitude in [0, pi]
    :param phi: longitude in [0, 2*pi]
    """
    # this is from http://functions.wolfram.com/Polynomials/SphericalHarmonicY/20/01/01/
    # which can be derived from (x^2-1)d/dx P^ml = sqrt(1-x^2) P^(m+1)l -mx P^ml
    dydt = m/tan(theta) * Y(m, l, theta, phi)
    if m < l:
        # for m==l, Y^(m+1)_l=0
        # note we fiddle with phi to obtain the desired exp(im phi)
        # despite raising m to m+1
        dydt += sqrt((l-m)*(l+m+1)) * Y(m+1, l, theta, phi*m/(m+1))
    return dydt


def to_spherical(X):
    """Convert Cartesian 3D coordinates X to spherical r, theta, phi

    :param X: Cartesian 3D coordinates
    :returns r, theta, phi: radius, longitude, colatitude"""
    r = sqrt(X[0]**2+X[1]**2+X[2]**2)
    theta = acos(X[2]/r)
    phi = atan2(X[1], X[0])
    return r, theta, phi


class SphericalStokesSolution(AnalyticalStokesSolution):
    r"""Base class for solutions in spherical shell domains.

    This implements analytical solutions based on a poloidal function
    of the form.

    .. math::

        \mathcal{P}(r,\theta,\varphi) = \mathcal{P}_l(r)Y_lm(\theta, \varphi)

    so that the velocity solution takes the form:

    .. math::

        \mathbf{u} = \nabla\times\left(\mathbf{r}\times\nabla\mathcal{P}\right)

    The function :math:`\mathcal{P}_l` and its derivative should
    be implement in methods :meth:`Pl` and :meth:`dPldr` in the subclass.
    """
    def __init__(self, l, m, Rp=2.22, Rm=1.22, nu=1.0, g=1.0):
        super(SphericalStokesSolution, self).__init__(Rp=Rp, Rm=Rm, nu=nu, g=g)
        self.l = l
        self.m = m

    def Pl(self, r):
        """Abstract method to be implemented by subclass.

        :param r: radius
        """
        raise NotImplementedError("Should implement Pl() method")

    def dPldr(self, r):
        """Abstract method to be implemented by subclass.

        :param r: radius
        """
        raise NotImplementedError("Should implement dPldr() method")

    def u_theta(self, r, theta, phi):
        """Return colatitudinal (southward) component of velocity.

        :param r: radius
        :param theta: co-latitude in [0, pi]
        :param phi: longitude in [0, 2*pi]
        """
        # u_theta = -1/r d/dtheta d/dr (r P)
        #         = -1/r dP/dtheta - d/dtheta d/dr P
        #         = -(1/r Pl + dPl/dr) dY/dtheta
        return -(self.Pl(r)/r + self.dPldr(r)) * dYdtheta(self.m, self.l, theta, phi)

    def u_phi(self, r, theta, phi):
        """Return longitudinal (eastward) component of velocity.

        :param r: radius
        :param theta: co-latitude in [0, pi]
        :param phi: longitude in [0, 2*pi]
        """
        # u_phi = -1/(r sin(theta)) d/dphi d/dr (r P)
        #       = -1/(r sin(theta)) * (dP/dphi + r d/dphi d/dr P)
        #       = -1/sin(theta) * (Pl/r + dPl/dr) * dY/dphi
        return -(self.Pl(r)/r + self.dPldr(r)) / sin(theta) * dYdphi(self.m, self.l, theta, phi)

    def u_r(self, r, theta, phi):
        """Return radial component of velocity.

        :param r: radius
        :param theta: co-latitude in [0, pi]
        :param phi: longitude in [0, 2*pi]
        """
        # u_r = 1/r \Lambda^2 P = -1/r l(l+1) P
        return -self.l*(self.l+1)*self.Pl(r)*Y(self.m, self.l, theta, phi)/r

    def tau_rr(self, r, theta, phi):
        """Return radial component of deviatoric stress.

        :param r: radius
        :param theta: co-latitude in [0, pi]
        :param phi: longitude in [0, 2*pi]
        """
        # tau_rr = 2 nu d/dr (1/r \Lambda^2 P)
        #        = -2 nu l(l+1) [1/r dP/dr - 1/r^2 P]
        return -2*self.nu*self.l*(self.l+1)*(self.dPldr(r) - self.Pl(r)/r)*Y(self.m, self.l, theta, phi)/r

    def radial_stress(self, r, theta, phi):
        """Return radial component of stress.

        :param r: radius
        :param theta: co-latitude in [0, pi]
        :param phi: longitude in [0, 2*pi]
        """
        return self.tau_rr(r, theta, phi) - self.p(r, theta, phi)

    def radial_stress_cartesian(self, X):
        """Return radial component of stress at Cartesian location.

        :param X: 3D Cartesian location"""
        return self.radial_stress(*to_spherical(X))

    def velocity_cartesian(self, X):
        """Return Cartesian velocity at Cartesian location.

        :param X: 3D Cartesian location"""
        r, theta, phi = to_spherical(X)
        if theta < 1e-7*pi or theta > pi*(1.-1e-7):
            # workaround pole problem by averaging in 4 points near the pole
            dx = 1e-6*r
            return tuple(numpy.mean([self.velocity_cartesian((x, y, X[2])) for x, y in [[dx, dx], [-dx, dx], [dx, -dx], [-dx, -dx]]], axis=0))
        ur = self.u_r(r, theta, phi)
        uth = self.u_theta(r, theta, phi)
        uph = self.u_phi(r, theta, phi)
        costh = cos(theta)
        req = sqrt(X[0]**2+X[1]**2)
        return (X[0]/r*ur + X[0]/req*costh*uth - X[1]/req*uph,
                X[1]/r*ur + X[1]/req*costh*uth + X[0]/req*uph,
                X[2]/r*ur - sin(theta)*uth)

    def pressure_cartesian(self, X):
        """Return pressure solution at Cartesian location.

        :param X: 3D Cartesian location"""
        return self.p(*to_spherical(X))


class SphericalStokesSolutionDelta(SphericalStokesSolution):
    r"""Base class for solutions in spherical shell domains with delta(r-r') forcing

    This implements the analytical solution in one half (above or below r')
    of the domain which is based on a poloidal function

    .. math ::

        \mathcal{P}(r,\theta,\varphi) = \mathcal{P}_l(r)Y_lm(\theta, \varphi)

    and velocity

    .. math ::

        \mathbf{u} = \nabla\times\left(\mathbf{r}\times\nabla\mathcal{P}\right)

    where for biharmonic solutions, :math:`\mathcal{P}_l(r)` is determinded by four
    coefficients, A, B, C, and D.
    """
    def __init__(self, ABCD, l, m, Rp=2.22, Rm=1.22, nu=1.0, g=1.0):
        """:param ABCD: list or tuple of 4 floats, coefficients for the 4 independent biharmonic
           solutions.
        :param l: degree of the harmonic
        :param m: order of the harmonic
        :param Rp: outer radius
        :param Rm: inner radius
        :param nu: viscosity
        :param g: forcing strength"""
        super(SphericalStokesSolutionDelta, self).__init__(l, m, Rp=Rp, Rm=Rm, nu=nu, g=g)
        self.ABCD = ABCD
        A, B, C, D = self.ABCD
        self.G = -2*nu*(l+1)*(2*l+3)*C
        self.H = -2*nu*l*(2*l-1)*D

    def Pl(self, r):
        """Radial part of poloidal function

        :param r: radius"""
        A, B, C, D = self.ABCD
        l = self.l
        return A*r**l + B*r**(-l-1) + C*r**(l+2) + D*r**(-l+1)

    def dPldr(self, r):
        """Radial derivative of radial part of poloidal function

        :param r: radius"""
        A, B, C, D = self.ABCD
        l = self.l
        return l*A*r**(l-1) + (-l-1)*B*r**(-l-2) + (l+2)*C*r**(l+1) + (-l+1)*D*r**-l

    def p(self, r, theta, phi):
        """Pressure solution

        :param r: radius
        :param theta: co-latitude in [0, pi]
        :param phi: longitude in [0, 2*pi]
        """
        l, m = self.l, self.m
        return (self.G*r**l + self.H*r**(-l-1))*Y(m, l, theta, phi)


class SphericalStokesSolutionSmooth(SphericalStokesSolution):
    r"""Base class for solutions in spherical shell domains with r^k forcing

    This implements the analytical solution in one half (above or below r')
    of the domain which is based on a poloidal function

    .. math ::

        \mathcal{P}(r,\theta,\varphi) = \mathcal{P}_l(r)Y_lm(\theta, \varphi)

    and velocity

    .. math ::

        \mathbf{u} = \nabla\times\left(\mathbf{r}\times\nabla\mathcal{P}\right)

    where the solution :math:`\mathcal{P}_l(r)` is determinded by four
    coefficients, A, B, C, and D of four independent biharmonic solutions,
    and one coefficient E associated with the inhomogenous part.
    """
    def __init__(self, ABCDE, k, l, m, Rp=2.22, Rm=1.22, nu=1.0, g=1.0):
        """:param ABCDE: list or tuple of 5 floats, coefficients for the 4 biharmonic
           and one inhomogenous solution.
        :param l: degree of the harmonic
        :param m: order of the harmonic
        :param Rp: outer radius
        :param Rm: inner radius
        :param nu: viscosity
        :param g: forcing strength"""
        super(SphericalStokesSolutionSmooth, self).__init__(l, m, Rp=Rp, Rm=Rm, nu=nu, g=g)
        self.k = k
        self.ABCDE = ABCDE
        A, B, C, D, E = self.ABCDE
        self.G = -2*nu*(l+1)*(2*l+3)*C
        self.H = -2*nu*l*(2*l-1)*D
        self.K = -g*(k+2)/((k+1)*(k+2)-l*(l+1))/Rp**k

    def Pl(self, r):
        """Radial part of poloidal function

        :param r: radius"""
        A, B, C, D, E = self.ABCDE
        l, k = self.l, self.k
        return A*r**l + B*r**(-l-1) + C*r**(l+2) + D*r**(-l+1) + E*r**(k+3)

    def dPldr(self, r):
        """Radial derivative of radial part of poloidal function

        :param r: radius"""
        A, B, C, D, E = self.ABCDE
        l, k = self.l, self.k
        return l*A*r**(l-1) + (-l-1)*B*r**(-l-2) + (l+2)*C*r**(l+1) + (-l+1)*D*r**-l + (k+3)*E*r**(k+2)

    def p(self, r, theta, phi):
        """Pressure solution

        :param r: radius
        :param theta: co-latitude in [0, pi]
        :param phi: longitude in [0, 2*pi]
        """
        k, l, m = self.k, self.l, self.m
        return (self.G*r**l + self.H*r**(-l-1) + self.K*r**(k+1))*Y(m, l, theta, phi)

    def delta_rho(self, r, theta, phi):
        r"""Perturbation density :math:`\rho'` in forcing term: :math:`g\rho'\hat r`

        :param r: radius
        :param theta: co-latitude in [0, pi]
        :param phi: longitude in [0, 2*pi]
        """
        k, l, m = self.k, self.l, self.m
        return r**k * Y(m, l, theta, phi) / self.Rp**k

    def delta_rho_cartesian(self, X):
        r"""Perturbation density :math:`\rho'` in forcing term: :math:`g\rho'\hat r`

        :param X: 3D Cartesian coordinate"""
        return self.delta_rho(*to_spherical(X))


class SphericalStokesSolutionSmoothFreeSlip(SphericalStokesSolutionSmooth):
    """Analytical Solution in cylindrical domain with smooth r^k forcing and free-slip boundary conditions"""
    def __init__(self, l, m, k, Rp=2.22, Rm=1.22, nu=1.0, g=1.0):
        """:param l: degree of the harmonic
        :param m: order of the harmonic
        :param k: polynomial order of forcing
        :param Rp: outer radius
        :param Rm: inner radius
        :param nu: viscosity
        :param g: forcing strength"""
        ABCDE = coefficients_sphere_smooth_fs(Rp, Rm, k, l, g, nu)
        super(SphericalStokesSolutionSmoothFreeSlip, self).__init__(ABCDE, k, l, m, Rp=Rp, Rm=Rm, nu=nu, g=g)


class SphericalStokesSolutionSmoothZeroSlip(SphericalStokesSolutionSmooth):
    """Analytical Solution in cylindrical domain with smooth r^k forcing and zero-slip boundary conditions"""
    def __init__(self, l, m, k, Rp=2.22, Rm=1.22, nu=1.0, g=1.0):
        """:param l: degree of the harmonic
        :param m: order of the harmonic
        :param k: polynomial order of forcing
        :param Rp: outer radius
        :param Rm: inner radius
        :param nu: viscosity
        :param g: forcing strength"""
        ABCDE = coefficients_sphere_smooth_ns(Rp, Rm, k, l, g, nu)
        super(SphericalStokesSolutionSmoothZeroSlip, self).__init__(ABCDE, k, l, m, Rp=Rp, Rm=Rm, nu=nu, g=g)


class SphericalStokesSolutionDeltaFreeSlip(SphericalStokesSolutionDelta):
    """Analytical Solution in cylindrical domain with delta(r-r') forcing and free-slip boundary conditions"""
    def __init__(self, l, m, sign, Rp=2.22, Rm=1.22, rp=1.72, nu=1.0, g=1.0):
        r""":param l: degree of the harmonic
        :param m: order of the harmonic
        :param sign: +1 for upper half solution r'<r<Rp
                     -1 for lower half solution Rm<r<r'
        :param Rp: outer radius
        :param Rm: inner radius
        :param nu: viscosity
        :param g: forcing strength"""
        ABCD = coefficients_sphere_delta_fs(Rp, Rm, rp, l, g, nu, sign)
        super(SphericalStokesSolutionDeltaFreeSlip, self).__init__(ABCD, l, m, Rp=Rp, Rm=Rm, nu=nu, g=g)


class SphericalStokesSolutionDeltaZeroSlip(SphericalStokesSolutionDelta):
    """Analytical Solution in cylindrical domain with delta(r-r') forcing and zero-slip boundary conditions"""
    def __init__(self, l, m, sign, Rp=2.22, Rm=1.22, rp=1.72, nu=1.0, g=1.0):
        r""":param l: degree of the harmonic
        :param m: order of the harmonic
        :param sign: +1 for upper half solution r'<r<Rp
                     -1 for lower half solution Rm<r<r'
        :param Rp: outer radius
        :param Rm: inner radius
        :param nu: viscosity
        :param g: forcing strength"""
        ABCD = coefficients_sphere_delta_ns(Rp, Rm, rp, l, g, nu, sign)
        super(SphericalStokesSolutionDeltaZeroSlip, self).__init__(ABCD, l, m, Rp=Rp, Rm=Rm, nu=nu, g=g)
