
"""
Flow Over Joukowski Airfoil
by Conformal Mapping
"""
import numpy as np
import math
from matplotlib import pyplot as plt

"""
For cambered airfoil use parameters
c = 1.0
h, k = 0.1, 0.1                          
R_0 =math.sqrt((c - h)**2 + k**2)

For symmetric use
c = 1.0                      
h, k = -0.15, 0
R_0 = 1.15     
"""
# Joukowski transform parameters
c = 1.0                                  # transform parameter
h, k = -0.15, 0                          # center of circle in z plane
R_0 = 1.15                               # circle radius

# ===== Joukoski transform definitions ======

# curve in z plane
n = 500
x = np.linspace(-R_0 + h, R_0 + h, n)
yu = np.sqrt(R_0**2 - (x - h)**2) + k    # upper curve
yl = -np.sqrt(R_0**2 - (x - h)**2) + k   # lower curve

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

#===== generating the grid ===========================

# grid in polar coordinates
Rlim = 5                          # domain limit in r 
Nr, Ntheta = 100, 145             # number of grid points in r and theta
r = np.linspace(R_0, Rlim, Nr)
theta = np.linspace(0, 2 * np.pi, Ntheta)   
R, T = np.meshgrid(r, theta)

# convert polar grid to cartesian 
X = R * np.cos(T) + h
Y = R * np.sin(T) + k
Z = X + Y*1j
plt.figure('grid')
plt.scatter(X, Y, s=1)
plt.axis('equal')

# Joukoski transform on grid
zeta_grid = Z + c**2 / Z
plt.figure('zeta plane grid')
plt.scatter(zeta_grid.real, zeta_grid.imag, s=1)
plt.axis('equal')


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