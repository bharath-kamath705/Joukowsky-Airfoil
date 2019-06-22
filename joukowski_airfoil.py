
"""
Flow Over Joukowski Airfoil
by Conformal Mapping
"""
import numpy as np
import math
from matplotlib import pyplot as plt

"""
For cambered airfoil use parameters
==================================================================
c = 1.0                             # transform parameter
h, k = 0.1, 0.1                     # center of circle in z plane                   
R_0 =math.sqrt((c - h)**2 + k**2)   # circle radius

For symmetric airfoil
---------------------
c = 1.0                             # transform parameter
h, k = -0.15, 0                     # center of circle in z plane
R_0 = 1.15                          # circle radius 
-------------------------------------------------------------------
"""



# Joukowski transform parameters
c = 1.0                                  # transform parameter
h, k = -0.15, 0                          # center of circle in z plane
R_0 = 1.15                               # circle radius

# ===== Joukoski transform curves ============================

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


# plot z plane and zeta plane curves
fig1, ax1 = plt.subplots(1,2)
fig1.suptitle('$z$ and $ \\xi $ plane curves', fontsize=22)

ax1[0].plot(x, yu), ax1[0].plot(x, yl)
ax1[0].axis('equal'), ax1[0].grid(True)
ax1[0].set_xlabel('$x$', fontsize=20)
ax1[0].set_ylabel('$iy$', fontsize=20)

ax1[1].plot(zeta_l.real, zeta_l.imag)
ax1[1].plot(zeta_u.real, zeta_u.imag)
ax1[1].axis('equal'), ax1[1].grid(True)
ax1[1].set_xlabel('$ Re(\\xi) $', fontsize=20)
ax1[1].set_ylabel('$Img (\\xi) $', fontsize=20)
#----------------------------------------------------------

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

# Joukoski transform on grid
Z = X + Y*1j
zeta_grid = Z + c**2 / Z

# plot z plane and zeta plane grids
fig2, ax2 = plt.subplots(1,2)
fig2.suptitle('$z$ and $ \\xi $ plane grid', fontsize=22)

ax2[0].scatter(X, Y, s=1)
ax2[0].axis('equal')
ax2[0].set_xlabel('$x$', fontsize=20)
ax2[0].set_ylabel('$iy$', fontsize=20)

ax2[1].scatter(zeta_grid.real, zeta_grid.imag, s=1)
ax2[1].axis('equal')
ax2[1].set_xlabel('$ Re(\\xi) $', fontsize=20)
ax2[1].set_ylabel('$Img (\\xi) $', fontsize=20)
#----------------------------------------------------------
