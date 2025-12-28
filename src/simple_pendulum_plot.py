########################### Libraries ###########################
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

########################### Data for one pendulum and plotting / creation of .mp4 file ###########################
file = np.loadtxt('../data/simple_pendulum_output.txt',delimiter=' ')
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
plt.savefig('../image/phase_space_1pend.png', dpi=300)

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
plt.savefig('../image/trajectory_1pend.png', dpi=300)