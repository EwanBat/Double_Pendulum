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

########################### Data loading and processing ###########################
def load_pendulum_data(filepath):
    """Load pendulum data from file and convert to Cartesian coordinates."""
    file = np.loadtxt(filepath, delimiter=' ')
    t = file[:, 0]
    xi1 = file[:, 1]
    eta1 = file[:, 2]
    xi2 = file[:, 3]
    eta2 = file[:, 4]
    
    x1, y1, z1 = stereographic_to_cartesian(xi1, eta1)
    x2, y2, z2 = stereographic_to_cartesian(xi2, eta2)
    
    # Adjust second mass position relative to first mass
    x2 += x1
    y2 += y1
    z2 += z1
    
    return {
        't': t, 'xi1': xi1, 'eta1': eta1, 'xi2': xi2, 'eta2': eta2,
        'x1': x1, 'y1': y1, 'z1': z1, 'x2': x2, 'y2': y2, 'z2': z2
    }

########################### Plotting functions ###########################
def plot_3d_trajectory(data, output_path='../image/trajectory_2pend.png', dpi=300):
    """Plot 3D trajectory of the double pendulum."""
    plt.figure(figsize=(6, 6))
    ax = plt.figure().add_subplot(111, projection='3d')
    ax.set_xlim([-2, 2])
    ax.set_ylim([-2, 2])
    ax.set_zlim([-2, 2])
    ax.plot([0], [0], [0], color='black', label='Pivot Point', marker='o')
    ax.plot(data['x1'], data['y1'], data['z1'], color='blue', label='Mass 1 Trajectory')
    ax.plot(data['x2'], data['y2'], data['z2'], color='orange', label='Mass 2 Trajectory')
    ax.set_title('3D Trajectory of the Double Pendulum')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()

def plot_error_between_pendulums(data1, data2, output_path='../image/error_2doublependulum.png', dpi=300):
    """Plot error between two double pendulums over time."""
    error_mass1 = np.sqrt((data1['x1'] - data2['x1'])**2 + 
                          (data1['y1'] - data2['y1'])**2 + 
                          (data1['z1'] - data2['z1'])**2)
    error_mass2 = np.sqrt((data1['x2'] - data2['x2'])**2 + 
                          (data1['y2'] - data2['y2'])**2 + 
                          (data1['z2'] - data2['z2'])**2)
    
    plt.figure(figsize=(8, 5))
    plt.plot(data1['t'], error_mass1, label='Error Mass 1', color='blue')
    plt.plot(data1['t'], error_mass2, label='Error Mass 2', color='orange')
    plt.yscale('log')
    plt.title('Error Between Two Double Pendulums Over Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Error (log scale)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()

def plot_phase_space(data, output_path='../image/phase_space_2pend.png', dpi=300):
    """Plot phase space for both masses."""
    gs = gridspec.GridSpec(1, 2)
    fig = plt.figure(figsize=(12, 5))
    
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(data['xi1'], data['eta1'], color='blue')
    ax1.set_title('Phase Space of Mass 1')
    ax1.set_xlabel('ξ₁')
    ax1.set_ylabel('η₁')
    ax1.grid()
    
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(data['xi2'], data['eta2'], color='orange')
    ax2.set_title('Phase Space of Mass 2')
    ax2.set_xlabel('ξ₂')
    ax2.set_ylabel('η₂')
    ax2.grid()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()

########################### Main function ###########################
def main():
    """Main function to generate all plots."""
    # Load data for both pendulums
    data1 = load_pendulum_data('../data/double_pendulum_output.txt')
    data2 = load_pendulum_data('../data/double_pendulum_output_2.txt')
    
    # Generate plots
    plot_3d_trajectory(data1)
    plot_error_between_pendulums(data1, data2)
    plot_phase_space(data1)
    
    print("Tous les graphiques ont été générés avec succès!")

if __name__ == '__main__':
    main()