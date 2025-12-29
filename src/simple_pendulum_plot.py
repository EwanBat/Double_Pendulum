########################### Libraries ###########################
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation

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

########################### Data for multiple pendulums and plotting, animation ###########################
# nb_simulations = 5  # Number of simulations
# multi_output_file = '../data/simple_pendulum_output'
# Coord_files = np.array([np.loadtxt(f'{multi_output_file}_sim_{i}.txt', delimiter=' ') for i in range(nb_simulations)])
# t_multi = Coord_files[0][:,0]
# Thetas_multi = Coord_files[:, :, 1]
# dThetas_multi = Coord_files[:, :, 2]
# Phis_multi = Coord_files[:, :, 3]
# dPhis_multi = Coord_files[:, :, 4]

# xs_multi = rho * np.sin(Thetas_multi) * np.cos(Phis_multi)
# ys_multi = rho * np.sin(Thetas_multi) * np.sin(Phis_multi)
# zs_multi = -rho * np.cos(Thetas_multi)

# ## Plotting in phase space for multiple pendulums and the map of the difference for close pendulums
# fig = plt.figure(figsize=(12,6))
# gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
# colors = cm.get_cmap('viridis')(np.linspace(0, 1, nb_simulations))

# ax_ps = plt.subplot(gs[0])
# ax_ps.set_title('Phase Space: Theta vs dPhi for Multiple Pendulums')
# ax_ps.set_xlabel('Theta (rad)')
# ax_ps.set_ylabel('dPhi (rad/s)')
# ax_ps.grid()
# ax_ps.set_xlim(Thetas_multi.min() - 0.1, Thetas_multi.max() + 0.1)
# ax_ps.set_ylim(dPhis_multi.min() - 0.1, dPhis_multi.max() + 0.1)

# ax_map = plt.subplot(gs[1])
# ax_map.set_title('Difference Map between Close Pendulums')
# ax_map.set_xlabel('Theta (rad)')
# ax_map.set_ylabel('dPhi (rad/s)')
# ax_map.grid()

# # Initialisation des artistes pour l'animation
# points_ps = []
# points_map = []

# for i in range(nb_simulations):
#     point, = ax_ps.plot([], [], 'o', color=colors[i], markersize=8, label=f'Pendulum {i+1}')
#     points_ps.append(point)

# if nb_simulations >= 2:
#     for i in range(nb_simulations - 1):
#         point, = ax_map.plot([], [], 'o', color=colors[i], markersize=6, alpha=0.6)
#         points_map.append(point)

# ax_ps.legend(loc='upper right', fontsize=8)
# time_text = ax_ps.text(0.02, 0.98, '', transform=ax_ps.transAxes, 
#                        verticalalignment='top', fontsize=10, 
#                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# def update(frame):
#     # Mise à jour du texte temporel
#     time_text.set_text(f'Time: {t_multi[frame]:.2f} s')
    
#     # Mise à jour de l'espace des phases
#     for i in range(nb_simulations):
#         points_ps[i].set_data([Thetas_multi[i, frame]], [dPhis_multi[i, frame]])
    
#     # Mise à jour de la carte des différences
#     if nb_simulations >= 2:
#         for i in range(nb_simulations - 1):
#             delta_theta = np.abs(Thetas_multi[i, frame] - Thetas_multi[i+1, frame])
#             delta_dphi = np.abs(dPhis_multi[i, frame] - dPhis_multi[i+1, frame])
#             points_map[i].set_data([delta_theta], [delta_dphi])
    
#     artists = points_ps + points_map + [time_text]
#     return artists

# ani = FuncAnimation(fig, update, frames=len(t_multi), blit=True, interval=20, repeat=True)
# ani.save('../image/phase_space_multiple_pendulums.gif', writer='pillow', dpi=150, fps=30)
# plt.close()

# print("Animation saved successfully!")
