import assess
import math
import numpy as np
import scipy

# Y(), dYdphi(), and dYdtheta() depend on scipy sph_harm_y
# test whether the conventions correspond to what I assumed:
# * phi is longitude [0, 2pi]
# * theta is co-latitude [0,pi]
# As that went wrong when changing from sph_harm -> sph_harm_y
# sph_harm_y uses now uses the same coordinate convention as assess
# (where sph_harm didn't)
# Further we use:
# * l for degree (n in scipy)
# * m for order  (same in scipy


def assumed_normalisation(l, m):
    # based on scipy docstring
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm_y.html
    return np.sqrt((2*l+1)/(4*math.pi))*np.sqrt(math.factorial(l-m)/math.factorial(l+m))


def test_Y():
    for phi in np.linspace(0, 2*math.pi, 10):
        for theta in np.linspace(0, math.pi, 10):
            for l in range(5):
                for m in range(0, l):
                    N = assumed_normalisation(l, m)
                    # NOTE: lpmv takes order, degree (i.e. m, l)
                    P = scipy.special.lpmv(m, l, math.cos(theta))
                    expected = N * P * math.cos(m*phi)
                    Y = assess.Y(l, m, theta, phi)
                    np.testing.assert_allclose(Y, expected, atol=1e-15)
