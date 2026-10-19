"""FractionalODELab core: Caputo fractional ODE solvers and Caputo-derivative operators.

Problem:  D_C^alpha y(t) = f(t, y(t)),  y^(k)(0) = y0[k],  k = 0..ceil(alpha)-1,
y scalar or vector, 0 < alpha < 2.

Solvers
  pece         fractional Adams-Bashforth-Moulton predictor-corrector (Diethelm-Ford-Freed)
  trapezoidal  implicit product-integration trapezoidal rule (PECE corrector solved by Newton)
  lubich       implicit Lubich convolution quadrature (BDF1 / BDF2 generating functions,
               optional starting weights)
  vo_pece      variable-order PECE, alpha = alpha(t) in (0, 1)
Operators
  caputo_l1          L1 scheme on nodal values, 0 < alpha < 1
  caputo_strong_fdm  Strong-FDM scheme of A. H. Salas, any non-integer alpha > 0
  residual           discrete residual D^alpha y - f(t, y) of a computed solution
"""
from math import gamma, ceil, factorial, comb, fsum, lgamma, exp, log
import numpy as np

__all__ = ["taylor_part", "pece", "trapezoidal", "lubich_weights", "lubich", "vo_pece",
           "caputo_l1", "caputo_strong_fdm", "residual", "mittag_leffler"]


# ----------------------------------------------------------------- helpers
def _as_ic(y0, alpha):
    """Initial data: for alpha <= 1, y0 = y(0); for alpha > 1, y0 = [y(0), y'(0), ...]
    (a single value is padded with zero derivatives)."""
    m = max(1, ceil(alpha))
    if m == 1 or np.isscalar(y0) or isinstance(y0, np.ndarray):
        ic = [np.atleast_1d(np.asarray(y0, dtype=float))]
    else:
        ic = [np.atleast_1d(np.asarray(v, dtype=float)) for v in y0]
    while len(ic) < m:
        ic.append(np.zeros_like(ic[0]))
    return ic[:m]


def taylor_part(t, y0, alpha):
    """T(t) = sum_{k < ceil(alpha)} y^(k)(0) t^k / k!, returned with shape (len(t), d)."""
    ic = _as_ic(y0, alpha)
    t = np.atleast_1d(np.asarray(t, dtype=float))
    return sum(np.outer(t ** k / factorial(k), ic[k]) for k in range(len(ic)))


def _wrap(f, d):
    def F(t, y):
        return np.atleast_1d(np.asarray(f(t, y if d > 1 else y[0]), dtype=float))
    return F


def _newton(G, x0, tol, maxit):
    x = np.array(x0, dtype=float)
    for _ in range(maxit):
        g = G(x); m = x.size
        J = np.empty((m, m))
        for i in range(m):
            e = np.zeros(m); e[i] = 1e-7 * (1 + abs(x[i]))
            J[:, i] = (G(x + e) - g) / e[i]
        dx = np.linalg.solve(J, -g); x = x + dx
        if np.max(np.abs(dx)) < tol * (1 + np.max(np.abs(x))):
            break
    return x


# ----------------------------------------------------------------- PECE / trapezoidal
def pece(f, alpha, y0, T, N, corrector_iterations=1, implicit=False, tol=1e-13, maxit=50):
    """Fractional Adams-Bashforth-Moulton method on t_n = n T/N.  Returns (t, y).
    corrector_iterations = M gives P(EC)^M E; implicit=True solves the corrector exactly.
    Error O(h^min(2, 1+alpha)) for smooth f (Diethelm, Ford & Freed 2004)."""
    h = T / N
    t = h * np.arange(N + 1)
    Tp = taylor_part(t, y0, alpha)
    d = Tp.shape[1]
    F = _wrap(f, d)
    y = np.empty((N + 1, d)); fv = np.empty((N + 1, d))
    y[0] = Tp[0]; fv[0] = F(t[0], y[0])
    cP = h ** alpha / gamma(alpha + 1); cC = h ** alpha / gamma(alpha + 2)
    k = np.arange(N + 2, dtype=float)
    kA = k ** alpha; kA1 = k ** (alpha + 1)
    for n in range(N):
        j = np.arange(n + 1)
        b = kA[n + 1 - j] - kA[n - j]
        a = kA1[n - j + 2] + kA1[n - j] - 2 * kA1[n - j + 1]
        a[0] = n ** (alpha + 1) - (n - alpha) * (n + 1) ** alpha
        yp = Tp[n + 1] + cP * (b @ fv[: n + 1])
        hist = Tp[n + 1] + cC * (a @ fv[: n + 1])
        if implicit:
            yp = _newton(lambda x: x - hist - cC * F(t[n + 1], x), yp, tol, maxit)
        else:
            for _ in range(corrector_iterations):
                yp = hist + cC * F(t[n + 1], yp)
        y[n + 1] = yp; fv[n + 1] = F(t[n + 1], yp)
    return t, (y[:, 0] if d == 1 else y)


def trapezoidal(f, alpha, y0, T, N, tol=1e-13, maxit=50):
    """Implicit fractional product-integration trapezoidal rule (stiff problems)."""
    return pece(f, alpha, y0, T, N, implicit=True, tol=tol, maxit=maxit)


# ----------------------------------------------------------------- Lubich CQ
def lubich_weights(alpha, N, order=2):
    """Weights w_0..w_N of omega(z) = delta(z)^alpha, delta = 1-z (order 1) or (1-z)(3-z)/2 (order 2)."""
    g = np.empty(N + 1); g[0] = 1.0
    for j in range(1, N + 1):
        g[j] = g[j - 1] * (1 - (alpha + 1) / j)
    if order == 1:
        return g
    q = np.empty(N + 1); q[0] = 1.0
    for j in range(1, N + 1):
        q[j] = q[j - 1] * (1 - (alpha + 1) / j) / 3
    return (1.5 ** alpha) * np.convolve(g, q)[: N + 1]


def _starting_exponents(alpha, order, max_count=4):
    qs = sorted({round(k + l * alpha, 12) for k in range(order + 1) for l in range(int(order / alpha) + 2)
                 if 0 < k + l * alpha < order - 1e-12})
    return qs[:max_count]


def lubich(f, alpha, y0, T, N, order=2, starting=True, tol=1e-13, maxit=50):
    """Implicit Lubich convolution quadrature for D^alpha (y - T) = f(t, y):
        sum_k w_{n-k} u_k + sum_j W_{n,j} u_j = h^alpha f(t_n, T_n + u_n),   u = y - T,
    with starting weights W_{n,j}, j = 1..s, exact for t^q, q in {k + l alpha < order} (s <= 4).
    The first s steps are solved as one coupled system, then step by step (Newton)."""
    h = T / N
    t = h * np.arange(N + 1)
    Tp = taylor_part(t, y0, alpha)
    d = Tp.shape[1]
    F = _wrap(f, d)
    w = lubich_weights(alpha, N, order)
    ha = h ** alpha
    qs = _starting_exponents(alpha, order) if starting else []
    s = min(len(qs), N)
    W = np.zeros((N + 1, max(s, 1)))
    if s:
        js = np.arange(1, s + 1, dtype=float)
        Vinv = np.linalg.inv(np.array([js ** q for q in qs[:s]]))
        kk = np.arange(N + 1, dtype=float)
        kq = [kk ** q for q in qs[:s]]
        for n in range(1, N + 1):
            wr = w[n - 1::-1][:n]
            R = np.array([fsum(np.concatenate(([gamma(q + 1) / gamma(q + 1 - alpha) * n ** (q - alpha)],
                                               -wr * kq[i][1: n + 1]))) for i, q in enumerate(qs[:s])])
            W[n, :s] = Vinv @ R
    u = np.zeros((N + 1, d))

    def res(n, un, uu):
        acc = w[0] * un + (w[n - 1:0:-1] @ uu[1:n] if n > 1 else 0.0)
        if s:
            acc = acc + W[n, :s] @ uu[1: s + 1]
        return acc - ha * F(t[n], Tp[n] + un)

    if s:
        def Gs(x):
            uu = u.copy(); uu[1: s + 1] = x.reshape(s, d)
            return np.concatenate([res(n, uu[n], uu) for n in range(1, s + 1)])
        u[1: s + 1] = _newton(Gs, np.zeros(s * d), tol, maxit).reshape(s, d)
    for n in range(s + 1, N + 1):
        u[n] = _newton(lambda x: res(n, x, u), u[n - 1], tol, maxit)
    y = u + Tp
    return t, (y[:, 0] if d == 1 else y)


# ----------------------------------------------------------------- variable order
def vo_pece(f, alpha_fun, y0, T, N, corrector_iterations=1):
    """Variable-order PECE (A. H. Salas): the step n -> n+1 uses the PECE weights with the
    order frozen at alpha_{n+1} = alpha(t_{n+1}), 0 < alpha(t) < 1."""
    h = T / N
    t = h * np.arange(N + 1)
    y0v = np.atleast_1d(np.asarray(y0, dtype=float)); d = y0v.size
    F = _wrap(f, d)
    y = np.empty((N + 1, d)); fv = np.empty((N + 1, d))
    y[0] = y0v; fv[0] = F(0.0, y[0])
    for n in range(N):
        al = float(alpha_fun(t[n + 1]))
        if not 0 < al < 1:
            raise ValueError("alpha(t) must stay in (0, 1)")
        j = np.arange(n + 1, dtype=float)
        b = (n + 1 - j) ** al - (n - j) ** al
        a = (n - j + 2) ** (al + 1) + (n - j) ** (al + 1) - 2 * (n - j + 1) ** (al + 1)
        a[0] = n ** (al + 1) - (n - al) * (n + 1) ** al
        yp = y0v + h ** al / gamma(al + 1) * (b @ fv[: n + 1])
        hist = y0v + h ** al / gamma(al + 2) * (a @ fv[: n + 1])
        for _ in range(corrector_iterations):
            yp = hist + h ** al / gamma(al + 2) * F(t[n + 1], yp)
        y[n + 1] = yp; fv[n + 1] = F(t[n + 1], yp)
    return t, (y[:, 0] if d == 1 else y)


# ----------------------------------------------------------------- Caputo operators
def caputo_l1(u, alpha, h):
    """L1 approximation of D_C^alpha u(t_n), n = 1..N, from nodal values (0 < alpha < 1)."""
    u = np.asarray(u, dtype=float); N = len(u) - 1
    k = np.arange(N + 1, dtype=float); b = (k + 1) ** (1 - alpha) - k ** (1 - alpha)
    du = np.diff(u)
    c = h ** (-alpha) / gamma(2 - alpha)
    return np.array([c * (b[:n][::-1] @ du[:n]) for n in range(1, N + 1)])


def caputo_strong_fdm(u, alpha, h, N, d0=None):
    """Strong-FDM approximation of D_C^alpha u(t_n), n = 1..N (A. H. Salas), alpha > 0 not integer.

    With p = ceil(alpha) + 1 (one integration by parts),
        D^a u(t) = t^(p-1-a)/G(p-a) u^(p-1)(0) + 1/G(p-a) int_0^t (t-s)^(p-1-a) u^(p)(s) ds;
    u^(p) is frozen on each step, replaced by the strong p-th difference on the p+1 equally
    spaced points (j-1+r/p) h, and the kernel is integrated exactly:
        D^a u(t_n) ~ t_n^(p-1-a)/G(p-a) u^(p-1)(0)
             + p^p/(h^p G(p+1-a)) sum_j [((n-j+1)h)^(p-a) - ((n-j)h)^(p-a)] Delta_p u_j,
        Delta_p u_j = sum_r (-1)^(p-r) C(p,r) u((j-1+r/p) h).
    p = 2 (0<a<1): midpoint second difference;  p = 3 (1<a<2): third difference at thirds.
    u: callable;  d0 = u^(p-1)(0) (one-sided difference estimate when omitted)."""
    p = int(ceil(alpha)) + 1
    if d0 is None:
        e = 1e-4
        d0 = (u(e) - u(0.0)) / e if p == 2 else (u(2 * e) - 2 * u(e) + u(0.0)) / e ** 2
    jj = np.arange(1, N + 1)
    D = sum((-1) ** (p - r) * comb(p, r) * np.array([u((j - 1 + r / p) * h) for j in jj], dtype=float)
            for r in range(p + 1))
    c = p ** p / (h ** p * gamma(p + 1 - alpha))
    out = np.empty(N)
    for n in range(1, N + 1):
        j = jj[:n]
        wgt = ((n - j + 1) * h) ** (p - alpha) - ((n - j) * h) ** (p - alpha)
        out[n - 1] = (n * h) ** (p - 1 - alpha) / gamma(p - alpha) * d0 + c * (wgt @ D[:n])
    return out


def residual(f, alpha, t, y, operator="strong_fdm", d0=None):
    """Residual r_n = D_C^alpha y(t_n) - f(t_n, y_n), n = 1..N, of a computed scalar solution.
    operator="strong_fdm" (quintic-spline interpolant for the interior points; second order for
    smooth solutions) or "l1" (nodal values; robust for weakly singular solutions, 0<alpha<1)."""
    from scipy.interpolate import make_interp_spline
    h = t[1] - t[0]; N = len(t) - 1
    D = caputo_l1(y, alpha, h) if operator == "l1" else \
        caputo_strong_fdm(make_interp_spline(t, y, k=5), alpha, h, N, d0=d0)
    return D - np.array([f(tn, yn) for tn, yn in zip(t[1:], y[1:])])


def mittag_leffler(alpha, z, beta=1.0):
    """E_{alpha,beta}(z) by its power series with compensated summation
    (reference values for |z| <= 6, alpha >= 0.3)."""
    z = float(z)
    if z == 0:
        return 1.0 / gamma(beta)
    terms = []
    for k in range(4000):
        lt = k * log(abs(z)) - lgamma(alpha * k + beta)
        terms.append(exp(lt) * (1 if z > 0 or k % 2 == 0 else -1))
        if k > 10 and lt < -60:
            break
    return fsum(terms)
