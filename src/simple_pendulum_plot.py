########################### Libraries ###########################
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation

########################### Data loading and processing ###########################
def load_simple_pendulum_data(filepath, rho=1.0):
    """Load simple pendulum data from file and convert to Cartesian coordinates."""
    file = np.loadtxt(filepath, delimiter=' ')
    t = file[:, 0]
    theta = file[:, 1]
    dtheta = file[:, 2]
    phi = file[:, 3]
    dphi = file[:, 4]
    
    # Convert to Cartesian coordinates
    x = rho * np.sin(theta) * np.cos(phi)
    y = rho * np.sin(theta) * np.sin(phi)
    z = -rho * np.cos(theta)
    
    return {
        't': t, 'theta': theta, 'dtheta': dtheta, 'phi': phi, 'dphi': dphi,
        'x': x, 'y': y, 'z': z, 'rho': rho
    }

########################### Plotting functions ###########################
def plot_phase_space(data, output_path='../image/phase_space_1pend.png', dpi=300):
    """Plot phase space for the simple pendulum."""
    fig = plt.figure(figsize=(10, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
    
    ax0 = plt.subplot(gs[0])
    ax0.plot(data['theta'], data['phi'], color='purple')
    ax0.set_title('Phase Space: Theta vs Phi')
    ax0.set_xlabel('Theta (rad)')
    ax0.set_ylabel('Phi (rad)')
    ax0.grid()
    
    ax1 = plt.subplot(gs[1])
    ax1.plot(data['dtheta'], data['dphi'], color='brown')
    ax1.set_title('Phase Space: dTheta vs dPhi')
    ax1.set_xlabel('dTheta (rad/s)')
    ax1.set_ylabel('dPhi (rad/s)')
    ax1.grid()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()

def plot_3d_trajectory(data, output_path='../image/trajectory_1pend.png', dpi=300):
    """Plot 3D trajectory of the simple pendulum."""
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(data['x'], data['y'], data['z'], color='green')
    ax.set_title('3D Trajectory of the Pendulum')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()

########################### Main function ###########################
def main():
    """Main function to generate all plots."""
    # Load data
    data = load_simple_pendulum_data('../data/simple_pendulum_output.txt', rho=1.0)
    
    # Generate plots
    plot_phase_space(data)
    plot_3d_trajectory(data)
    
    print("Tous les graphiques ont été générés avec succès!")

if __name__ == '__main__':
    main()
