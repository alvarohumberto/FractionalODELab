"""FractionalODELab: Caputo fractional ODE solvers (PECE, Lubich CQ, variable order)
and Caputo-derivative operators (L1, Strong-FDM)."""
from .core import (taylor_part, pece, trapezoidal, lubich_weights, lubich, vo_pece,
                   caputo_l1, caputo_strong_fdm, residual, mittag_leffler)
__version__ = "0.1.0"
