"""
Thermodynamics and Phase Separation in Liquid Mixtures - Visualization
Phase Diagram and 3D Free Energy Surface
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D
import warnings
warnings.filterwarnings('ignore')

# Global Style Settings - Dark Theme
plt.style.use('dark_background')
plt.rcParams.update({
    'figure.facecolor': '#000000',
    'axes.facecolor': '#000000',
    'axes.edgecolor': '#3a3a5a',
    'axes.labelcolor': '#e0e0ff',
    'text.color': '#e0e0ff',
    'xtick.color': '#8080a0',
    'ytick.color': '#8080a0',
    'grid.color': '#2a2a4a',
    'font.size': 12,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'font.family': 'sans-serif',
})

# Custom Color Scheme
COLORS = {
    'primary': '#00d4ff',
    'secondary': '#ff6b6b',
    'accent': '#a855f7',
    'success': '#22c55e',
    'warning': '#f59e0b',
}

def plot_phase_diagram_and_free_energy():
    """Plot Phase Diagram and 3D Free Energy Surface"""
    fig = plt.figure(figsize=(18, 8), facecolor='#000000')
    gs = fig.add_gridspec(1, 2, wspace=0.3)
    
    # Left: Phase Diagram (T vs c)
    ax1 = fig.add_subplot(gs[0])
    T = np.linspace(0.1, 1.2, 100)
    c_spinodal_low = 0.5 - 0.3 * np.sqrt(1 - T/1.2)
    c_spinodal_high = 0.5 + 0.3 * np.sqrt(1 - T/1.2)
    c_binodal_low = 0.5 - 0.45 * np.sqrt(1 - T/1.2)
    c_binodal_high = 0.5 + 0.45 * np.sqrt(1 - T/1.2)
    
    # Fill regions
    ax1.fill_betweenx(T, c_binodal_low, c_spinodal_low, alpha=0.5, 
                      color=COLORS['warning'], label='Metastable Region')
    ax1.fill_betweenx(T, c_spinodal_high, c_binodal_high, alpha=0.5, 
                      color=COLORS['warning'])
    ax1.fill_betweenx(T, c_spinodal_low, c_spinodal_high, alpha=0.6, 
                      color=COLORS['secondary'], label='Spinodal Region')
    
    # Plot lines
    ax1.plot(c_binodal_low, T, '-', color=COLORS['primary'], linewidth=3, label='Binodal Line')
    ax1.plot(c_binodal_high, T, '-', color=COLORS['primary'], linewidth=3)
    ax1.plot(c_spinodal_low, T, '--', color=COLORS['secondary'], linewidth=2.5, label='Spinodal Line')
    ax1.plot(c_spinodal_high, T, '--', color=COLORS['secondary'], linewidth=2.5)
    
    # Critical Point
    ax1.scatter([0.5], [1.2], s=300, marker='*', c=[COLORS['accent']], 
                edgecolors='white', linewidths=2.5, zorder=5, label='Critical Point')
    ax1.annotate('Critical Point\n(Tc, cc)', (0.5, 1.2), textcoords="offset points", 
                 xytext=(60, -10), ha='left',
                 fontsize=13, fontweight='bold', color=COLORS['accent'],
                 bbox=dict(boxstyle='round,pad=0.5', facecolor='#1a1a2e', edgecolor=COLORS['accent'], linewidth=2))
    
    ax1.set_xlabel('Concentration c', fontweight='bold')
    ax1.set_ylabel('Temperature T / Tc', fontweight='bold')
    ax1.set_title('Binary Phase Diagram', fontweight='bold', fontsize=18, pad=15)
    ax1.legend(loc='upper right', framealpha=0.9, fontsize=11, fancybox=True, shadow=True)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1.35)
    ax1.grid(True, alpha=0.3, linestyle='--')
    
    # Right: 3D Free Energy Surface
    ax2 = fig.add_subplot(gs[1], projection='3d', facecolor='#000000')
    
    # Construct free energy surface f(c, T)
    c_grid = np.linspace(0, 1, 120)
    T_grid = np.linspace(0.3, 1.2, 100)
    C, TT = np.meshgrid(c_grid, T_grid)
    
    # Temperature dependent free energy
    Tc = 1.2
    with np.errstate(divide='ignore', invalid='ignore'):
        entropy = TT * (C * np.log(C + 1e-10) + (1 - C) * np.log(1 - C + 1e-10))
    F = (1 - TT/Tc) * C**2 * (1 - C)**2 + 0.1 * entropy
    F = np.nan_to_num(F, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Create gradient colormap
    surf_cmap = LinearSegmentedColormap.from_list(
        'surface', [
            '#0a0a1a',  # Very dark blue
            '#1e40af',  # Dark blue
            '#3b82f6',  # Blue
            '#22d3ee',  # Cyan
            '#a5f3fc',  # Light cyan
            '#fbbf24',  # Yellow
        ], N=256
    )
    
    # Plot surface
    surf = ax2.plot_surface(C, TT, F, cmap=surf_cmap, alpha=0.9,
                            linewidth=0, antialiased=True,
                            rstride=1, cstride=1, shade=True)
    
    # Add contour projection at the bottom
    contour_offset = F.min() - 0.03
    ax2.contour(C, TT, F, zdir='z', offset=contour_offset, 
                levels=20, cmap='cool', alpha=0.7, linewidths=1.5)
    
    # Labels and title
    ax2.set_xlabel('Concentration c', labelpad=12, fontweight='bold')
    ax2.set_ylabel('Temperature T/Tc', labelpad=12, fontweight='bold')
    ax2.set_zlabel('Free Energy F', labelpad=12, fontweight='bold')
    ax2.set_title('Free Energy Surface F(c, T)', fontweight='bold', fontsize=18, pad=20)
    
    # Adjust viewing angle
    ax2.view_init(elev=25, azim=-65)
    
    # Customize 3D pane appearance
    ax2.xaxis.pane.fill = False
    ax2.yaxis.pane.fill = False
    ax2.zaxis.pane.fill = False
    ax2.xaxis.pane.set_edgecolor('#3a3a5a')
    ax2.yaxis.pane.set_edgecolor('#3a3a5a')
    ax2.zaxis.pane.set_edgecolor('#3a3a5a')
    ax2.grid(True, alpha=0.2)
    
    # Add colorbar
    cbar = fig.colorbar(surf, ax=ax2, shrink=0.6, aspect=15, pad=0.1)
    cbar.set_label('Free Energy', rotation=270, labelpad=20, fontweight='bold')
    
    # Main title
    fig.suptitle('Cahn-Hilliard Model', 
                 fontsize=20, fontweight='bold', color='white', y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return fig

# Main Entry Point
if __name__ == "__main__":
    fig = plot_phase_diagram_and_free_energy()
    plt.savefig('phase_diagram_and_free_energy.png', dpi=300, bbox_inches='tight')
