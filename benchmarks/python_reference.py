"""Python reference values for the Wolfram Language cross-check (mathematica/RunTests.wls).
Problem 'relax': D^a y = -y;  problem 'bench': Diethelm-Ford-Freed benchmark on [0,1]."""
from pathlib import Path
import csv
from math import gamma
from fractional_ode_lab import pece, trapezoidal, lubich


def bench(a):
    return lambda t, y: (40320 / gamma(9 - a) * t ** (8 - a) - 3 * gamma(5 + a / 2) / gamma(5 - a / 2) * t ** (4 - a / 2)
                         + 9 / 4 * gamma(a + 1) + (1.5 * t ** (a / 2) - t ** 4) ** 3 - abs(y) ** 1.5)


rows = []
for a in (0.3, 0.5, 0.8, 1.5, 1.8):
    y0r = [1.0, 0.0] if a > 1 else 1.0
    y0b = [0.0, 0.0] if a > 1 else 0.0
    for method, fun in (("PECE", lambda f, y0, T, N: pece(f, a, y0, T, N)),
                        ("PECE3", lambda f, y0, T, N: pece(f, a, y0, T, N, corrector_iterations=3)),
                        ("TRAP", lambda f, y0, T, N: trapezoidal(f, a, y0, T, N)),
                        ("LUB2", lambda f, y0, T, N: lubich(f, a, y0, T, N))):
        rows.append(["relax", method, a, 100, repr(float(fun(lambda t, y: -y, y0r, 2.0, 100)[1][-1]))])
        rows.append(["bench", method, a, 100, repr(float(fun(bench(a), y0b, 1.0, 100)[1][-1]))])
out = Path(__file__).resolve().parent / "python_reference.csv"
with out.open("w", newline="") as f:
    w = csv.writer(f); w.writerow(["problem", "method", "alpha", "N", "y_final"]); w.writerows(rows)
print("wrote", out, len(rows), "rows")
