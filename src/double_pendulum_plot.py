########################### Libraries ###########################
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import os
########################### Stereographic to Cartesian conversion ###########################
def stereographic_to_cartesian(xi, eta):
    denom = 1 + xi**2 + eta**2
    x = 2 * xi / denom
    y = 2 * eta / denom
    z = ( -1 + xi**2 + eta**2 ) / denom
    return x, y, z

########################### Data for double pendulum and plotting / creation of .mp4 file ###########################
file2 = np.loadtxt('../data/double_pendulum_output.txt',delimiter=' ')
t2 = file2[:,0]
xi1_p1 = file2[:,1]
eta1_p1 = file2[:,2]
xi2_p1 = file2[:,3]
eta2_p1 = file2[:,4]
x1_p1, y1_p1, z1_p1 = stereographic_to_cartesian(xi1_p1, eta1_p1)
x2_p1, y2_p1, z2_p1 = stereographic_to_cartesian(xi2_p1, eta2_p1)
x2_p1 += x1_p1
y2_p1 += y1_p1
z2_p1 += z1_p1

## Trajectory plot for double pendulum
plt.figure(figsize=(6,6))
ax = plt.figure().add_subplot(111, projection='3d')
ax.set_xlim([-2, 2])
ax.set_ylim([-2, 2])
ax.set_zlim([-2, 2])
ax.plot([0], [0], [0], color='black', label='Pivot Point', marker='o')
ax.plot(x1_p1, y1_p1, z1_p1, color='blue', label='Mass 1 Trajectory')
ax.plot(x2_p1, y2_p1, z2_p1, color='orange', label='Mass 2 Trajectory')
ax.set_title('3D Trajectory of the Double Pendulum')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.legend()
plt.tight_layout()
plt.savefig('../image/trajectory_2pend.png', dpi=300)

########################### Comparison between 2 double pendulums ###########################
if os.path.exists('../data/double_pendulum_output_bis.txt'):
    file3 = np.loadtxt('../data/double_pendulum_output_bis.txt',delimiter=' ')
    t3 = file3[:,0]
    xi1_p2 = file3[:,1]
    eta1_p2 = file3[:,2]
    xi2_p2 = file3[:,3]
    eta2_p2 = file3[:,4]
    x1_p2, y1_p2, z1_p2 = stereographic_to_cartesian(xi1_p2, eta1_p2)
    x2_p2, y2_p2, z2_p2 = stereographic_to_cartesian(xi2_p2, eta2_p2)
    x2_p2 += x1_p2
    y2_p2 += y1_p2
    z2_p2 += z1_p2

    err_1 = np.sqrt((x1_p1 - x1_p2)**2 + (y1_p1 - y1_p2)**2 + (z1_p1 - z1_p2)**2)
    err_2 = np.sqrt((x2_p1 - x2_p2)**2 + (y2_p1 - y2_p2)**2 + (z2_p1 - z2_p2)**2)

    plt.figure(figsize=(10,5))
    plt.plot(t2, err_1, label='Error Mass 1', color='blue')
    plt.plot(t2, err_2, label='Error Mass 2', color='orange')
    plt.yscale('log')
    plt.title('Error between Two Double Pendulum Simulations')
    plt.xlabel('Time (s)')
    plt.ylabel('Error (m)')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('../image/error_double_pendulum.png', dpi=300)