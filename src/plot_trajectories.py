########################### Libraries ###########################
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

########################### Data for one pendulum and plotting / creation of .mp4 file ###########################
file = np.loadtxt('res_1pend.txt',delimiter=' ')
t = file[:,0]
Theta = file[:,1]
dTheta = file[:,2]
Phi = file[:,3]
dPhi = file[:,4]
rho = 1.0  # Length of the pendulum

## Plotting phase space
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
ax0 = plt.subplot(gs[0])
ax0.plot(Theta, dTheta, color='blue')
ax0.set_title('Phase Space: Θ vs dΘ')
ax0.set_xlabel('Θ (rad)')
ax0.set_ylabel('dΘ (rad/s)')
ax0.grid()

ax1 = plt.subplot(gs[1])
ax1.plot(Phi % 2*np.pi, dPhi, color='red')
ax1.set_title('Phase Space: Φ vs dΦ')
ax1.set_xlabel('Φ (rad)')
ax1.set_ylabel('dΦ (rad/s)')
ax1.grid()
plt.tight_layout()
plt.savefig('image/phase_space_1pend.png', dpi=300)

## Trajectory plot
plt.figure(figsize=(6,6))
x = rho * np.sin(Theta) * np.cos(Phi)
y = rho * np.sin(Theta) * np.sin(Phi)
z = -rho * np.cos(Theta)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, color='green')
ax.set_title('3D Trajectory of the Pendulum')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.tight_layout()
plt.savefig('image/trajectory_1pend.png', dpi=300)

########################### Data for double pendulum and plotting / creation of .mp4 file ###########################
file2 = np.loadtxt('cart_2pend.txt',delimiter=' ')
t2 = file2[:,0]
x1_p1 = file2[:,1]
y1_p1 = file2[:,2]
z1_p1 = file2[:,3]
x2_p1 = file2[:,4]
y2_p1 = file2[:,5]
z2_p1 = file2[:,6]

## Trajectory plot for double pendulum
plt.figure(figsize=(6,6))
ax = plt.figure().add_subplot(111, projection='3d')
ax.plot(x1_p1, y1_p1, z1_p1, color='blue', label='Mass 1 Trajectory')
ax.plot(x2_p1, y2_p1, z2_p1, color='orange', label='Mass 2 Trajectory')
ax.set_title('3D Trajectory of the Double Pendulum')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.legend()
plt.tight_layout()
plt.savefig('image/trajectory_2pend.png', dpi=300)

########################### Comparison between 2 double pendulums ###########################
file3 = np.loadtxt('cart_2pend_bis.txt',delimiter=' ')
t3 = file3[:,0]
x1_p2 = file3[:,1]
y1_p2 = file3[:,2]
z1_p2 = file3[:,3]
x2_p2 = file3[:,4]
y2_p2 = file3[:,5]
z2_p2 = file3[:,6]

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
plt.savefig('image/error_double_pendulum.png', dpi=300)