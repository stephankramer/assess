"""Analytical solutions in cylindrical domains."""
from __future__ import division
from assess.smooth import coefficients_cylinder_smooth_fs, coefficients_cylinder_smooth_ns
from assess.delta import coefficients_cylinder_delta_fs, coefficients_cylinder_delta_ns
from math import sqrt, atan2, cos, sin


__all__ = ['CylindricalStokesSolutionSmoothFreeSlip', 'CylindricalStokesSolutionSmoothZeroSlip',
           'CylindricalStokesSolutionDeltaFreeSlip', 'CylindricalStokesSolutionDeltaZeroSlip']


class AnalyticalStokesSolution(object):
    """Base class for solutions in sperical or cylindrical shell domains."""
    def __init__(self, Rp=2.22, Rm=1.22, nu=1.0, g=1.0):
        """
        :param Rp: outer radius
        :param Rm: inner radius
        :param g: magnitude of source term
        :param nu: viscosity
        """


class CylindricalStokesSolution(AnalyticalStokesSolution):
    r"""Base class for solutions in cylindrical shell domains.

    This implements analytical solutions based on a streamfunction
    of the form.

    .. math::

        \psi(r,\varphi) = \psi_r(r)\sin(n\varphi)

    where :math:`\psi_r` should be defined in a method :meth:`psi_r`
    and its derivative in :meth:`dpsi_rdr`.
    """
    def __init__(self, n, Rp=2.22, Rm=1.22, nu=1.0, g=1.0):
        """Initialize basic parameters of analytical solution in cylindrical domain."""
        super(CylindricalStokesSolution, self).__init__(Rp=Rp, Rm=Rm, nu=nu, g=g)
        self.n = n

    def psi_r(self, r):
        """Abstract method to be implemented by subclass.

        :param r: radius
        """
        raise NotImplementedError("Should implement psi_r() method")

    def dpsi_rdr(self, r):
        """Abstract method to be implemented by subclass.

        :param r: radius"""
        raise NotImplementedError("Should implement dpsi_rdr() method")

    def u_r(self, r, phi):
        """Return radial component of velocity.

        :param r: radius
        :param phi: angle with x-axis."""
        dpsi_dphi = self.n*cos(self.n*phi)*self.psi_r(r)
        return -dpsi_dphi/r

    def u_phi(self, r, phi):
        """Return tangential component of velocity.

        :param r: radius
        :param phi: angle with x-axis."""
        dpsi_dr = sin(self.n*phi)*self.dpsi_rdr(r)
        return dpsi_dr

    def tau_rr(self, r, phi):
        """Return radial component of deviatoric stress.

        :param r: radius
        :param phi: angle with x-axis."""
        dpsi_dphi = self.n*cos(self.n*phi)*self.psi_r(r)
        dpsi_drdphi = self.n*cos(self.n*phi)*self.dpsi_rdr(r)
        return 2*self.nu*(dpsi_dphi/r**2 - dpsi_drdphi/r)

    def radial_stress(self, r, phi):
        """Return radial component of stress.

        :param r: radius
        :param phi: angle with x-axis."""
        return self.tau_rr(r, phi) - self.p(r, phi)

    def radial_stress_cartesian(self, X):
        """Return radial component of stress at Cartesian location.

        :param X: 2D Cartesian location"""
        r = sqrt(X[0]**2+X[1]**2)
        phi = atan2(X[1], X[0])
        return self.radial_stress(r, phi)

    def velocity_cartesian(self, X):
        """Return Cartesian velocity at Cartesian location.

        :param X: 2D Cartesian location"""
        r = sqrt(X[0]**2+X[1]**2)
        phi = atan2(X[1], X[0])
        ur = self.u_r(r, phi)
        ut = self.u_phi(r, phi)
        return [ur*X[0]/r - ut*X[1]/r, ur*X[1]/r + ut*X[0]/r]

    def pressure_cartesian(self, X):
        """Return pressure solution at Cartesian location.

        :param X: 2D Cartesian location"""
        r = sqrt(X[0]**2+X[1]**2)
        phi = atan2(X[1], X[0])
        return self.p(r, phi)


class CylindricalStokesSolutionDelta(CylindricalStokesSolution):
    r"""Base class for solutions in cylindrical shell domains with delta(r-r') forcing

    This implements the analytical solution in one half (above or below r')
    of the domain which is based on a biharmonic streamfunction

    .. math::

        \psi(r,\varphi) = \psi_r(r)\sin(n\varphi)

    determinded by 4 coefficients, A, B, C, and D.
    """
    def __init__(self, ABCD, n, Rp=2.22, Rm=1.22, nu=1.0, g=1.0):
        """:param ABCD: list or tuple of 4 floats, coefficients for the 4 biharmonic
           solutions of the streamfunction
        :param n: wave number
        :param Rp: outer radius
        :param Rm: inner radius
        :param nu: viscosity
        :param g: forcing strength"""
        super(CylindricalStokesSolutionDelta, self).__init__(n, Rp=Rp, Rm=Rm, nu=nu, g=g)
        self.ABCD = ABCD
        A, B, C, D = self.ABCD
        self.G = -4*self.nu*C*(n+1)
        self.H = -4*self.nu*D*(n-1)

    def psi_r(self, r):
        """Radial part of streamfunction

        :param r: radius"""
        A, B, C, D = self.ABCD
        n = self.n
        return A*r**n+B*r**(-n)+C*r**(n+2)+D*r**(-n+2)

    def dpsi_rdr(self, r):
        """Radial derivative of radial part of streamfunction

        :param r: radius"""
        A, B, C, D = self.ABCD
        n = self.n
        return A*n*r**(n-1) + B*-n*r**(-n-1) + C*(n+2)*r**(n+1) + D*(-n+2)*r**(-n+1)

    def p(self, r, phi):
        """Pressure solution

        :param r: radius
        :param phi: angle with x-axis"""
        n = self.n
        return (self.G*r**n + self.H*r**(-n))*cos(n*phi)


class CylindricalStokesSolutionSmooth(CylindricalStokesSolution):
    r"""Base class for solutions in cylindrical shell domains with r^k forcing

    The analytical solution is based on a streamfunction

    .. math:

        \psi(r,\varphi) = \psi_r(r)\sin(n\varphi)

    determinded by 4 coefficients, A, B, C, and D corresponding to 4 independent
    biharmonic (i.e. homogenous) solutions and a fifth coefficient E that
    corresponds to the inhomogeneous part.
    """
    def __init__(self, ABCDE, k, n, Rp=2.22, Rm=1.22, nu=1.0, g=1.0):
        """:param ABCDE: list or tuple of 5 floats, coefficients for the 4 biharmonic
           and one inhomogenous solutions of the streamfunction
        :param k: polynomial order of forcing
        :param n: wave number
        :param Rp: outer radius
        :param Rm: inner radius
        :param nu: viscosity
        :param g: forcing strength"""
        super(CylindricalStokesSolutionSmooth, self).__init__(n, Rp=Rp, Rm=Rm, nu=nu, g=g)
        self.k = k
        self.ABCDE = ABCDE
        A, B, C, D, E = self.ABCDE
        self.G = -4*self.nu*C*(n+1)
        self.H = -4*self.nu*D*(n-1)
        self.F = -g*(k + 1)*Rp**(-k)/((k+1)**2-n**2)

    def psi_r(self, r):
        """Radial part of streamfunction

        :param r: radius"""
        A, B, C, D, E = self.ABCDE
        n, k = self.n, self.k
        return A*r**n+B*r**(-n)+C*r**(n+2)+D*r**(-n+2)+E*r**(k+3)

    def dpsi_rdr(self, r):
        """Radial derivative of radial part of streamfunction

        :param r: radius"""
        A, B, C, D, E = self.ABCDE
        n, k = self.n, self.k
        return A*n*r**(n-1) + B*-n*r**(-n-1) + C*(n+2)*r**(n+1) + D*(-n+2)*r**(-n+1)+E*(k+3)*r**(k+2)

    def p(self, r, phi):
        """Pressure solution

        :param r: radius
        :param phi: angle with x-axis"""
        n, k = self.n, self.k
        return (self.G*r**n + self.H*r**(-n) + self.F*r**(k+1))*cos(n*phi)

    def delta_rho(self, r, phi):
        r"""Perturbation density :math:`\rho'` in forcing term: :math:`g\rho'\hat r`

        :param r: radius
        :param phi: angle with x-axis"""
        n, k = self.n, self.k
        return r**k * cos(n*phi) / self.Rp**k

    def delta_rho_cartesian(self, X):
        r"""Perturbation density :math:`\rho'` in forcing term: :math:`g\rho'\hat r`

        :param X: 2D Cartesian coordinate"""
        r = sqrt(X[0]**2+X[1]**2)
        phi = atan2(X[1], X[0])
        return self.delta_rho(r, phi)


class CylindricalStokesSolutionSmoothFreeSlip(CylindricalStokesSolutionSmooth):
    """Analytical Solution in cylindrical domain with smooth r^k forcing and free-slip boundary conditions"""
    def __init__(self, n, k, Rp=2.22, Rm=1.22, nu=1.0, g=1.0):
        """:param n: wave number
        :param k: polynomial order of forcing
        :param Rp: outer radius
        :param Rm: inner radius
        :param nu: viscosity
        :param g: forcing strength"""
        ABCDE = coefficients_cylinder_smooth_fs(Rp, Rm, k, n, g, nu)
        super(CylindricalStokesSolutionSmoothFreeSlip, self).__init__(ABCDE, n, k, Rp=Rp, Rm=Rm, nu=nu, g=g)


class CylindricalStokesSolutionSmoothZeroSlip(CylindricalStokesSolutionSmooth):
    """Analytical Solution in cylindrical domain with smooth r^k forcing and zero-slip boundary conditions"""
    def __init__(self, n, k, Rp=2.22, Rm=1.22, nu=1.0, g=1.0):
        """:param n: wave number
        :param k: polynomial order of forcing
        :param Rp: outer radius
        :param Rm: inner radius
        :param nu: viscosity
        :param g: forcing strength"""
        ABCDE = coefficients_cylinder_smooth_ns(Rp, Rm, k, n, g, nu)
        super(CylindricalStokesSolutionSmoothZeroSlip, self).__init__(ABCDE, n, k, Rp=Rp, Rm=Rm, nu=nu, g=g)


class CylindricalStokesSolutionDeltaFreeSlip(CylindricalStokesSolutionDelta):
    """Analytical Solution in cylindrical domain with delta(r-r') forcing and free-slip boundary conditions"""
    def __init__(self, n, sign, Rp=2.22, Rm=1.22, rp=1.72, nu=1.0, g=1.0):
        r""":param n: wave number
        :param sign: +1 for upper half solution r'<r<Rp
                     -1 for lower half solution Rm<r<r'
        :param Rp: outer radius
        :param Rm: inner radius
        :param nu: viscosity
        :param g: forcing strength"""
        ABCD = coefficients_cylinder_delta_fs(Rp, Rm, rp, n, g, nu, sign)
        super(CylindricalStokesSolutionDeltaFreeSlip, self).__init__(ABCD, n, Rp=Rp, Rm=Rm, nu=nu, g=g)


class CylindricalStokesSolutionDeltaZeroSlip(CylindricalStokesSolutionDelta):
    """Analytical Solution in cylindrical domain with delta(r-r') forcing and zero-slip boundary conditions"""
    def __init__(self, n, sign, Rp=2.22, Rm=1.22, rp=1.72, nu=1.0, g=1.0):
        r""":param n: wave number
        :param sign: +1 for upper half solution r'<r<Rp
                     -1 for lower half solution Rm<r<r'
        :param Rp: outer radius
        :param Rm: inner radius
        :param nu: viscosity
        :param g: forcing strength"""
        ABCD = coefficients_cylinder_delta_ns(Rp, Rm, rp, n, g, nu, sign)
        super(CylindricalStokesSolutionDeltaZeroSlip, self).__init__(ABCD, n, Rp=Rp, Rm=Rm, nu=nu, g=g)
