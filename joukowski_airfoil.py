
"""
Flow Over Joukowski Airfoil
by Conformal Mapping
"""
import numpy as np
import math
from matplotlib import pyplot as plt

# Joukowski transform parameters
c = 1.0
h, k = 0.1, 0.1    # center of circle in z plane
R = math.sqrt((c - h)**2 + k**2)         # circle radius

# ===== Joukoski transform definitions ======

# curve in z plane
n = 500
x = np.linspace(-R + h, R + h, n)
yu = np.sqrt(R**2 - (x - h)**2) + k    # upper curve
yl = -np.sqrt(R**2 - (x - h)**2) + k   # lower curve

# hack to fix NaNs in yu, yl if sqrt(very small number) occurs
yu[np.argwhere(np.isnan(yu))] = k
yl[np.argwhere(np.isnan(yl))] = k

zu = x + yu * 1j   # upper curve
zl = x + yl * 1j   # lower curve

# zeta plane curve
zeta_l = zl + c**2 / zl
zeta_u = zu + c**2 / zu


# ====== plot z plane and zeta plane curves ========

plt.figure('zeta plane')
plt.plot(zeta_l.real, zeta_l.imag)
plt.plot(zeta_u.real, zeta_u.imag)
plt.axis('equal'), plt.grid(True)

plt.figure('z plane')
plt.plot(x, yu), plt.plot(x, yl)
plt.axis('equal'), plt.grid(True)


## generate grid
#N = 50           # number of points in each direction
#
## axis limits of graphs
#x_start, x_end = -6.0, 6.0
#y_start, y_end = -6.0, 6.0
#x = np.linspace(x_start, x_end, N)
#y = np.linspace(y_start, y_end, N)
#X, Y = np.meshgrid(x, y)
#
#
## Joukowski transform definition
#c = 10.0
#z = X + Y * 1j
#xi = z + c**2 / 2