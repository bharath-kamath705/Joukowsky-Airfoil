
"""
Flow Over Joukowski Airfoil
by Conformal Mapping
"""
import numpy as np
import math
from matplotlib import pyplot as plt
import pflow

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
c = 1.0                             # transform parameter
h, k = -0.15, 0                     # center of circle in z plane
R_0 = 1.15                          # circle radius 

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

# z plane
ax1[0].plot(x, yu), ax1[0].plot(x, yl)
ax1[0].axis('equal'), ax1[0].grid(True)
ax1[0].set_xlabel('$x$', fontsize=20)
ax1[0].set_ylabel('$iy$', fontsize=20)

# zeta plane
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

# z plane
ax2[0].scatter(X, Y, s=1)
ax2[0].axis('equal')
ax2[0].set_xlabel('$x$', fontsize=20)
ax2[0].set_ylabel('$iy$', fontsize=20)

# zeta plane
ax2[1].scatter(zeta_grid.real, zeta_grid.imag, s=1)
ax2[1].axis('equal')
ax2[1].set_xlabel('$ Re(\\xi) $', fontsize=20)
ax2[1].set_ylabel('$Img (\\xi) $', fontsize=20)
#----------------------------------------------------------

# ============ Solving flow over the airfoil ===============
U = 1.0                           # Uniform flow velocity
aoa = 20.0 * math.pi / 180         # angle of attack
Dstr = R_0**2 * 2 * math.pi * U   # doublet strength


# grid in the zp (z prime) reference frame
Xp = (X - h) * np.cos(aoa) + (Y - k) * np.sin(aoa)
Yp = (Y - k) * np.cos(aoa) - (X - h) * np.sin(aoa)

# Kutta condition (stagnation point at trailing edge)
Vstr = -Yp[0, 0] * 4 * np.pi * U

# velocity field in zp plane
up = pflow.vortex([0], [0], [Vstr], Xp, Yp)[0] + pflow.doublet([0], [0],
                 [Dstr], Xp, Yp)[0] + pflow.freestream(U, 0, Xp, Yp)[0]

vp = pflow.vortex([0], [0], [Vstr], Xp, Yp)[1] + pflow.doublet([0], [0],
                 [Dstr], Xp, Yp)[1] + pflow.freestream(U, 0, Xp, Yp)[1]

# stream function 
psi = pflow.vortex([0], [0], [Vstr], Xp, Yp)[2] + pflow.doublet([0], [0], 
                  [Dstr], Xp, Yp)[2] + pflow.freestream(U, 0, Xp, Yp)[2]

# velocity field in the z plane
u = up * np.cos(-aoa) + vp * np.sin(-aoa)
v = vp * np.cos(-aoa) - up * np.sin(-aoa)

# velocity field in zeta plane
dzeta_dz = 1 - (c/Z)**2
V_zeta = (u - v * 1j) / dzeta_dz
u_zeta = V_zeta.real
v_zeta = -V_zeta.imag


# pressure coefficients
cp_z = pflow.cp_get(u, v, U)                 # z plane
cp_zeta = pflow.cp_get(u_zeta, v_zeta, U)    # zeta plane


# plot z plane and zeta plane streamlines
fig3, ax3 = plt.subplots(1,2)
fig3.suptitle('$z$ and $ \\xi $ plane streamlilnes', fontsize=22)

# z plane
ax3[0].contour(X, Y, psi, 30)            # streamlines
ax3[0].plot(x, yu), ax3[0].plot(x, yl)   # cylinder curves
ax3[0].set_xlabel('$x$', fontsize=20)
ax3[0].set_ylabel('$iy$', fontsize=20)
ax3[0].set_xlim([-Rlim, Rlim])
ax3[0].set_ylim([-Rlim, Rlim])


# zeta plane
ax3[1].contour(zeta_grid.real, zeta_grid.imag, psi, 30)  # streamlines
ax3[1].plot(zeta_l.real, zeta_l.imag)                    # airfoil curve
ax3[1].plot(zeta_u.real, zeta_u.imag)                    # airfoil curve
ax3[1].set_xlabel('$ Re(\\xi) $', fontsize=20)
ax3[1].set_ylabel('$Img (\\xi) $', fontsize=20)
ax3[1].set_xlim([-Rlim, Rlim])
ax3[1].set_ylim([-Rlim, Rlim])

# plot z plane and zeta plane velocity fields

fig4, ax4 = plt.subplots(1,2)
fig4.suptitle('$z$ and $ \\xi $ plane velocity field', fontsize=22)

# z plane
ax4[0].quiver(X, Y, u, v)            # streamlines
ax4[0].plot(x, yu), ax4[0].plot(x, yl)   # cylinder curves
ax4[0].set_xlabel('$x$', fontsize=20)
ax4[0].set_ylabel('$iy$', fontsize=20)
ax4[0].set_xlim([-Rlim, Rlim])
ax4[0].set_ylim([-Rlim, Rlim])


# zeta plane
ax4[1].quiver(zeta_grid.real, zeta_grid.imag, u_zeta, v_zeta)  # streamlines
ax4[1].plot(zeta_l.real, zeta_l.imag)                    # airfoil curve
ax4[1].plot(zeta_u.real, zeta_u.imag)                    # airfoil curve
ax4[1].set_xlabel('$ Re(\\xi) $', fontsize=20)
ax4[1].set_ylabel('$Img (\\xi) $', fontsize=20)
ax4[1].set_xlim([-Rlim, Rlim])
ax4[1].set_ylim([-Rlim, Rlim])

# plot pressure coefficient contours
fig5, ax5 = plt.subplots(1,2)
fig5.suptitle('Pressure coefficient $ C_p $ ', fontsize=22)
# z plane
contf1 = ax5[0].contourf(X, Y, cp_z, levels=np.linspace(-1, 1, 500), extend='both')
ax5[0].set_xlabel('$x$', fontsize=20)
ax5[0].set_ylabel('$iy$', fontsize=20)
ax5[0].set_xlim([-Rlim, Rlim])
ax5[0].set_ylim([-Rlim, Rlim])
cbar1 = plt.colorbar(contf1)
cbar1.set_label('$C_p$', fontsize=20)

# zeta plane
contf2 = ax5[1].contourf(zeta_grid.real, zeta_grid.imag, cp_zeta, levels=np.linspace(-1, 1, 100), extend='both')
ax5[1].set_xlabel('$ Re(\\xi) $', fontsize=20)
ax5[1].set_ylabel('$Img (\\xi) $', fontsize=20)
ax5[1].set_xlim([-Rlim, Rlim])
ax5[1].set_ylim([-Rlim, Rlim])