import pytest
import assess
from math import pi, cos, sin
from numpy.testing import assert_allclose

#######################
# the actual tests
######################


class _TestNoNormalFlow:
    """Test zero radial flow at top and bottom"""
    def test_no_normal_flow_top(self, case_top, angle):
        assert_allclose(case_top.u_r(case_top.Rp, *angle), 0, atol=1e-12)

    def test_no_normal_flow_bottom(self, case_bottom, angle):
        assert_allclose(case_bottom.u_r(case_bottom.Rm, *angle), 0, atol=1e-12)


class _TestFreeSlip(_TestNoNormalFlow):
    """Test zero shear stress in phi direction (and no normal from above)"""
    def test_free_slip_top(self, case_top, angle):
        assert_allclose(case_top.tau_rphi(case_top.Rp, *angle), 0, atol=1e-12)

    def test_free_slip_bottom(self, case_bottom, angle):
        assert_allclose(case_bottom.tau_rphi(case_bottom.Rm, *angle), 0, atol=1e-12)


class _TestSphericalFreeSlip(_TestFreeSlip):
    """Test zero shear stress in theta direction (and phi direction above, and no-normal flow)"""
    def test_free_slip_top_theta(self, case_top, angle):
        assert_allclose(case_top.tau_rtheta(case_top.Rp, *angle), 0, atol=1e-12)

    def test_free_slip_bottom_theta(self, case_bottom, angle):
        assert_allclose(case_bottom.tau_rtheta(case_bottom.Rm, *angle), 0, atol=1e-12)


class _TestZeroSlip(_TestNoNormalFlow):
    """Test zero velocity at top and bottom. Also seperately tests no-normal
    flow, so we test both u_r() and u_cartesian()"""
    def test_free_slip_top(self, case_top, top_X):
        assert_allclose(case_top.velocity_cartesian(top_X), 0, atol=1e-12)

    def test_free_slip_bottom(self, case_bottom, bottom_X):
        assert_allclose(case_bottom.velocity_cartesian(bottom_X), 0, atol=1e-12)


class _TestInterfaceConditions:
    """Test interface condition in r and phi directions."""
    def test_interface_conditions(self, case_top, case_bottom, angle, case_kwargs):
        rp = case_kwargs.get('rp', 1.72)
        assert_allclose(case_top.u_r(rp, *angle), case_bottom.u_r(rp, *angle))
        assert_allclose(case_top.u_phi(rp, *angle), case_bottom.u_phi(rp, *angle))
        assert_allclose(case_top.tau_rr(rp, *angle), case_bottom.tau_rr(rp, *angle))
        assert_allclose(case_top.tau_rphi(rp, *angle), case_bottom.tau_rphi(rp, *angle))


class _TestCylindricalInterfaceConditions(_TestInterfaceConditions):
    """Test interface condition in r and phi directions, including normal stress jump."""
    def test_interface_normal_stress(self, case_top, case_bottom, angle, case_kwargs):
        rp = case_kwargs.get('rp', 1.72)
        n, g = case_top.n, case_top.g
        phi = angle[0]
        assert_allclose(case_top.radial_stress(rp, *angle) - case_bottom.radial_stress(rp, *angle), g*cos(n*phi))


class _TestSphericalInterfaceConditions(_TestInterfaceConditions):
    """Test interface condition in r, phi and theta directions, including normal stress jump."""
    def test_interface_conditions_theta(self, case_top, case_bottom, angle, case_kwargs):
        rp = case_kwargs.get('rp', 1.72)
        assert_allclose(case_top.u_theta(rp, *angle), case_bottom.u_theta(rp, *angle))
        assert_allclose(case_top.tau_rr(rp, *angle), case_bottom.tau_rr(rp, *angle))
        assert_allclose(case_top.tau_rtheta(rp, *angle), case_bottom.tau_rtheta(rp, *angle))

    def test_interface_normal_stress(self, case_top, case_bottom, angle, case_kwargs):
        rp = case_kwargs.get('rp', 1.72)
        l, m, g = case_top.l, case_top.m, case_top.g
        assert_allclose(case_top.radial_stress(rp, *angle) - case_bottom.radial_stress(rp, *angle), g*assess.Y(l, m, *angle))


########################################################
# the fixtures that form the parameters for the tests
########################################################


class CylindricalFixtures:
    @pytest.fixture(params=[2, 3, 4, 5, 8])
    def n(self, request):
        return request.param

    @pytest.fixture(params=[0, pi/2, -pi/3])
    def angle(self, request):
        return [request.param]

    @pytest.fixture
    def top_X(self, case_top, angle):
        return [case_top.Rp * cos(angle[0]), case_top.Rp * sin(angle[0])]

    @pytest.fixture
    def bottom_X(self, case_bottom, angle):
        return [case_bottom.Rm * cos(angle[0]), case_bottom.Rm * sin(angle[0])]


class SphericalFixtures:
    @pytest.fixture(params=[2, 3, 4, 5, 8])
    def l(self, request):  # NOQA: E743
        return request.param

    @pytest.fixture(params=[-1, 0, 1])
    def m(self, request, l):
        return request.param*l

    @pytest.fixture(params=[[0.01, 0], [pi/2, 0], [pi-0.01, 0], [pi/2, 1.0], [1.5, 4.0]])
    def angle(self, request):
        return request.param

    @pytest.fixture
    def top_X(self, case_top, angle):
        return assess.from_spherical(case_top.Rp, *angle)

    @pytest.fixture
    def bottom_X(self, case_bottom, angle):
        return assess.from_spherical(case_bottom.Rm, *angle)


class SmoothFixtures:
    @pytest.fixture()
    def case_top(self, case):
        return case

    case_bottom = case_top

    @pytest.fixture(params=[
        {},
        {'Rm': 1.0, 'Rp': 2.0, 'nu': 3.0, 'g': 4.0},
        {'Rm': 0.9, 'Rp': 1.9, 'g': 9.81}])
    def case_kwargs(self, request):
        return request.param


class DeltaFixtures:
    @pytest.fixture(params=[
        {},
        {'rp': 2.},
        {'Rm': 1.0, 'rp': 1.2, 'Rp': 2.0, 'nu': 3.0, 'g': 4.0},
        {'Rm': 0.9, 'Rp': 1.9, 'g': 9.81}])
    def case_kwargs(self, request):
        return request.param


class SmoothCylindricalFixtures(CylindricalFixtures, SmoothFixtures):
    @pytest.fixture(params=[-7, -1, 0, 1, 2, 9])
    def k(self, request, n):
        k = request.param
        if (k+1)**2 == n**2 or (k+3)**2 == n**2:
            # NOTE: return None for not implemented combination
            return
        else:
            return k


class SmoothSphericalFixtures(SphericalFixtures, SmoothFixtures):
    @pytest.fixture(params=[-7, -1, 0, 1, 2, 9])
    def k(self, request, l):
        k = request.param
        if (k+1)*(k+2) == l*(l+1) or (k+3)*(k+4) == l*(l+1):
            # NOTE: return None for not implemented combination
            return
        else:
            return k


########################################################
# the 8 different cases
########################################################


class TestCylindricalStokesSolutionSmoothFreeSlip(SmoothCylindricalFixtures, _TestFreeSlip):
    @pytest.fixture()
    def case(self, n, k, case_kwargs):
        if k is None:
            pytest.xfail()
        return assess.CylindricalStokesSolutionSmoothFreeSlip(n, k, **case_kwargs)


class TestCylindricalStokesSolutionSmoothZeroSlip(SmoothCylindricalFixtures, _TestZeroSlip):
    @pytest.fixture()
    def case(self, n, k, case_kwargs):
        if k is None:
            pytest.xfail()
        return assess.CylindricalStokesSolutionSmoothZeroSlip(n, k, **case_kwargs)


class TestSphericalStokesSolutionSmoothFreeSlip(SmoothSphericalFixtures, _TestSphericalFreeSlip):
    @pytest.fixture()
    def case(self, l, m, k, case_kwargs):
        if k is None:
            pytest.xfail()
        return assess.SphericalStokesSolutionSmoothFreeSlip(l, m, k, **case_kwargs)


class TestSphericalStokesSolutionSmoothZeroSlip(SmoothSphericalFixtures, _TestZeroSlip):
    @pytest.fixture()
    def case(self, l, m, k, case_kwargs):
        if k is None:
            pytest.xfail()
        return assess.SphericalStokesSolutionSmoothZeroSlip(l, m, k, **case_kwargs)


class TestCylindricalStokesSolutionDeltaFreeSlip(CylindricalFixtures, DeltaFixtures, _TestFreeSlip, _TestCylindricalInterfaceConditions):
    @pytest.fixture()
    def case_top(self, n, case_kwargs):
        return assess.CylindricalStokesSolutionDeltaFreeSlip(n, +1, **case_kwargs)

    @pytest.fixture()
    def case_bottom(self, n, case_kwargs):
        return assess.CylindricalStokesSolutionDeltaFreeSlip(n, -1, **case_kwargs)


class TestCylindricalStokesSolutionDeltaZeroSlip(CylindricalFixtures, DeltaFixtures, _TestZeroSlip, _TestCylindricalInterfaceConditions):
    @pytest.fixture()
    def case_top(self, n, case_kwargs):
        return assess.CylindricalStokesSolutionDeltaZeroSlip(n, +1, **case_kwargs)

    @pytest.fixture()
    def case_bottom(self, n, case_kwargs):
        return assess.CylindricalStokesSolutionDeltaZeroSlip(n, -1, **case_kwargs)


class TestSphericalStokesSolutionDeltaFreeSlip(SphericalFixtures, DeltaFixtures, _TestSphericalFreeSlip, _TestSphericalInterfaceConditions):
    @pytest.fixture()
    def case_top(self, l, m, case_kwargs):
        return assess.SphericalStokesSolutionDeltaFreeSlip(l, m, +1, **case_kwargs)

    @pytest.fixture()
    def case_bottom(self, l, m, case_kwargs):
        return assess.SphericalStokesSolutionDeltaFreeSlip(l, m, -1, **case_kwargs)


class TestSphericalStokesSolutionDeltaZeroSlip(SphericalFixtures, DeltaFixtures, _TestZeroSlip, _TestSphericalInterfaceConditions):
    @pytest.fixture()
    def case_top(self, l, m, case_kwargs):
        return assess.SphericalStokesSolutionDeltaZeroSlip(l, m, +1, **case_kwargs)

    @pytest.fixture()
    def case_bottom(self, l, m, case_kwargs):
        return assess.SphericalStokesSolutionDeltaZeroSlip(l, m, -1, **case_kwargs)
