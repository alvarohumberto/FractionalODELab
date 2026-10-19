import numpy as np
import pytest
from math import gamma
from fractional_ode_lab import (pece, trapezoidal, lubich, vo_pece, caputo_l1,
                                caputo_strong_fdm, residual, mittag_leffler)


def rate(E):
    return np.log2(E[-2] / E[-1])


def bench(a):
    f = lambda t, y: (40320 / gamma(9 - a) * t ** (8 - a) - 3 * gamma(5 + a / 2) / gamma(5 - a / 2) * t ** (4 - a / 2)
                      + 9 / 4 * gamma(a + 1) + (1.5 * t ** (a / 2) - t ** 4) ** 3 - abs(y) ** 1.5)
    ex = lambda t: t ** 8 - 3 * t ** (4 + a / 2) + 9 / 4 * t ** a
    return f, ex


@pytest.mark.parametrize("a", [0.3, 0.5, 0.8, 1.0, 1.5, 1.8])
def test_pece_mittag_leffler(a):
    y0 = [1.0, 0.0] if a > 1 else 1.0
    E = [abs(pece(lambda t, y: -y, a, y0, 2.0, N)[1][-1] - mittag_leffler(a, -2.0 ** a)) for N in (200, 400)]
    assert E[-1] < 2e-5 and rate(E) > min(2, 1 + a) - 0.25


@pytest.mark.parametrize("a", [0.25, 0.5, 0.75, 1.25, 1.5, 1.75])
def test_benchmark_second_order(a):
    f, ex = bench(a)
    y0 = 0.0 if a < 1 else [0.0, 0.0]
    for solver in (lambda N: pece(f, a, y0, 1.0, N, corrector_iterations=3),
                   lambda N: trapezoidal(f, a, y0, 1.0, N),
                   lambda N: lubich(f, a, y0, 1.0, N)):
        E = [abs(solver(N)[1][-1] - ex(1.0)) for N in (160, 320)]
        assert E[-1] < 2e-4 and rate(E) > 1.6


def test_system_and_vo():
    fLV = lambda t, y: np.array([y[0] - 0.5 * y[0] * y[1], -y[1] + 0.5 * y[0] * y[1]])
    _, yf = pece(fLV, 0.8, np.array([2.0, 1.0]), 5.0, 3200)
    _, y = pece(fLV, 0.8, np.array([2.0, 1.0]), 5.0, 400)
    assert np.max(np.abs(y[-1] - yf[-1])) < 1e-3
    _, y1 = vo_pece(lambda t, y: -y, lambda t: 0.7, 1.0, 2.0, 200)
    _, y2 = pece(lambda t, y: -y, 0.7, 1.0, 2.0, 200)
    assert np.max(np.abs(y1 - y2)) == 0.0


@pytest.mark.parametrize("a", [0.3, 0.5, 0.8, 1.3, 1.7])
def test_strong_fdm_second_order_smooth(a):
    ex = sum(1.0 / gamma(k + 1 - a) for k in range(int(np.ceil(a)), 60))   # D^a e^t at t = 1
    d0 = 1.0
    E = [abs(caputo_strong_fdm(np.exp, a, 1.0 / N, N, d0=d0)[-1] - ex) for N in (80, 160)]
    assert rate(E) > 1.9


def test_l1_order():
    a = 0.5; ex = gamma(4.5) / gamma(4.5 - a)
    E = [abs(caputo_l1((np.arange(N + 1) / N) ** 3.5, a, 1.0 / N)[-1] - ex) for N in (80, 160)]
    assert abs(rate(E) - (2 - a)) < 0.1


def test_residual_decreases():
    a = 0.5
    f = lambda t, y: gamma(4) / gamma(4 - a) * t ** (3 - a) + 1 + t ** 3 - y
    R = [np.max(np.abs(residual(f, a, *pece(f, a, 1.0, 1.0, N, corrector_iterations=3), d0=0.0))) for N in (80, 160)]
    assert R[1] < R[0] / 3
