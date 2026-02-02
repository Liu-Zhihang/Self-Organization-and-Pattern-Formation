import numpy as np
import matplotlib.pyplot as plt

# Set black theme style
plt.style.use('dark_background')

def plot_all_in_one(save_fig=False):
    """
    Plot all charts in one figure with black theme background
    
    Args:
        save_fig (bool): Whether to save the figure
    """
    # Create a figure with three subplots
    fig = plt.figure(figsize=(18, 6))
    
    # Subplot 1: Potential functions with different parameters
    ax1 = fig.add_subplot(131)
    
    # Plot multiple potential curves
    u = np.linspace(-2.5, 2.5, 500)
    
    # Case r < 0 (monostable)
    r1, h1 = -1, 0.2
    F1 = -0.5 * r1 * u**2 + 0.25 * u**4 - h1 * u
    ax1.plot(u, F1, lw=2, label=f'$r={r1:.1f}, h={h1:.1f}$')
    
    # Case r > 0, h = 0 (symmetric bistable)
    r2, h2 = 1, 0
    F2 = -0.5 * r2 * u**2 + 0.25 * u**4 - h2 * u
    ax1.plot(u, F2, lw=2, label=f'$r={r2:.1f}, h={h2:.1f}$')
    
    # Case r > 0, h != 0 (asymmetric bistable)
    r3, h3 = 1, 0.1
    F3 = -0.5 * r3 * u**2 + 0.25 * u**4 - h3 * u
    ax1.plot(u, F3, lw=2, label=f'$r={r3:.1f}, h={h3:.1f}$')
    
    ax1.set_title('Potential Landscapes $F(u)$')
    ax1.set_xlabel('$u$')
    ax1.set_ylabel('$F(u)$')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Subplot 2: Hysteresis loop
    ax2 = fig.add_subplot(132)
    r = 1  # For hysteresis phenomenon, r must be positive
    h_vals = np.linspace(-0.5, 0.5, 400)
    u_stars = []
    
    for h in h_vals:
        # Find roots of f(u) = r*u - u^3 + h = 0
        coeffs = [-1, 0, r, h]
        roots = np.roots(coeffs)
        # Filter for real roots
        real_roots = roots[np.isreal(roots)].real
        for root in sorted(real_roots):
            u_stars.append((h, root))
            
    h_coords, u_coords = zip(*u_stars)
    
    # Stability analysis: f'(u) = r - 3u^2
    stability = r - 3 * np.array(u_coords)**2
    stable_h = [h for h, s in zip(h_coords, stability) if s < 0]
    stable_u = [u for u, s in zip(u_coords, stability) if s < 0]
    unstable_h = [h for h, s in zip(h_coords, stability) if s >= 0]
    unstable_u = [u for u, s in zip(u_coords, stability) if s >= 0]
    
    ax2.plot(stable_h, stable_u, '.', color='cyan', markersize=4, label='Stable Fixed Points')
    ax2.plot(unstable_h, unstable_u, '--', color='orange', dashes=(5, 5), lw=2, label='Unstable Fixed Points')
    
    # Add jump arrows
    h_c = (2 / (3 * np.sqrt(3))) * r**(3/2)
    u_sn_pos = np.sqrt(r/3)
    u_sn_neg = -np.sqrt(r/3)
    
    # Jump up
    ax2.arrow(h_c, u_sn_pos, 0, -2*u_sn_pos, head_width=0.02, head_length=0.1, 
              fc='white', ec='white', lw=1.5)
    # Jump down
    ax2.arrow(-h_c, u_sn_neg, 0, 2*u_sn_neg, head_width=0.02, head_length=0.1, 
              fc='white', ec='white', lw=1.5)
    
    ax2.set_title(f'Hysteresis Loop for $r={r:.2f}$')
    ax2.set_xlabel('$h$')
    ax2.set_ylabel('Fixed Points $u^*$')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Subplot 3: Cusp bifurcation region diagram
    ax3 = fig.add_subplot(133)
    r_vals = np.linspace(0, 1.5, 200)
    h_c = (2 / (3 * np.sqrt(3))) * r_vals**(3/2)
    
    ax3.plot(r_vals, h_c, 'cyan', lw=2, label='Saddle-Node Bifurcation Line')
    ax3.plot(r_vals, -h_c, 'cyan', lw=2)
    
    ax3.fill_between(r_vals, h_c, -h_c, color='lightblue', alpha=0.2)
    
    ax3.text(1.0, 0, 'Bistable Region', horizontalalignment='center', color='white')
    ax3.text(0.5, 0.3, 'Monostable Region', horizontalalignment='center', color='white')
    ax3.text(0.5, -0.3, 'Monostable Region', horizontalalignment='center', color='white')
    
    ax3.set_title('Cusp Bifurcation Diagram')
    ax3.set_xlabel('$r$')
    ax3.set_ylabel('$h$')
    ax3.axhline(0, color='grey', linestyle='--', alpha=0.5)
    ax3.axvline(0, color='grey', linestyle='--', alpha=0.5)
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure (if needed)
    if save_fig:
        plt.savefig('cusp_bifurcation_combined.png', dpi=300, bbox_inches='tight', facecolor='black')
        print("Figure saved as cusp_bifurcation_combined.png")
    
    plt.show()

def plot_potential(r, h):
    """
    Plots the potential F(u) for the cusp bifurcation.
    
    Args:
        r (float): The control parameter r.
        h (float): The control parameter h.
    """
    u = np.linspace(-2.5, 2.5, 500)
    F = -0.5 * r * u**2 + 0.25 * u**4 - h * u
    
    plt.figure(figsize=(8, 6))
    plt.plot(u, F, lw=2)
    plt.title(f'Potential Landscape $F(u)$ for $r={r:.2f}$, $h={h:.2f}$')
    plt.xlabel('$u$')
    plt.ylabel('$F(u)$')
    plt.grid(True)
    plt.ylim(min(F) - 1, max(F) + 1 if r < 0 else 1)
    plt.show()

def plot_hysteresis_loop(r):
    """
    Plots the fixed points u* as a function of h for a given r > 0,
    illustrating the hysteresis loop.
    
    Args:
        r (float): The control parameter r. Must be positive.
    """
    if r <= 0:
        print("Hysteresis is only observed for r > 0.")
        return
        
    h_vals = np.linspace(-0.5, 0.5, 400)
    u_stars = []
    
    for h in h_vals:
        # Find roots of f(u) = r*u - u^3 + h = 0
        coeffs = [-1, 0, r, h]
        roots = np.roots(coeffs)
        # Filter for real roots
        real_roots = roots[np.isreal(roots)].real
        for root in sorted(real_roots):
            u_stars.append((h, root))
            
    h_coords, u_coords = zip(*u_stars)
    
    # Stability analysis: f'(u) = r - 3u^2
    stability = r - 3 * np.array(u_coords)**2
    stable_h = [h for h, s in zip(h_coords, stability) if s < 0]
    stable_u = [u for u, s in zip(u_coords, stability) if s < 0]
    unstable_h = [h for h, s in zip(h_coords, stability) if s >= 0]
    unstable_u = [u for u, s in zip(u_coords, stability) if s >= 0]
    
    plt.figure(figsize=(8, 6))
    plt.plot(stable_h, stable_u, 'b.', markersize=4, label='Stable Fixed Points')
    plt.plot(unstable_h, unstable_u, 'r--', dashes=(5, 5), lw=2, label='Unstable Fixed Points')
    
    # Add arrows to show the jumps
    h_c = (2 / (3 * np.sqrt(3))) * r**(3/2)
    u_sn_pos = np.sqrt(r/3)
    u_sn_neg = -np.sqrt(r/3)
    
    # Jump up
    plt.arrow(h_c, u_sn_pos, 0, -2*u_sn_pos, head_width=0.02, head_length=0.1, fc='k', ec='k', lw=1.5)
    # Jump down
    plt.arrow(-h_c, u_sn_neg, 0, 2*u_sn_neg, head_width=0.02, head_length=0.1, fc='k', ec='k', lw=1.5)
    
    plt.title(f'Hysteresis Loop for $r={r:.2f}$')
    plt.xlabel('$h$')
    plt.ylabel('Fixed Points $u^*$')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_cusp_region():
    """
    Plots the cusp bifurcation diagram in the (r, h) parameter space.
    """
    r_vals = np.linspace(0, 1.5, 200)
    h_c = (2 / (3 * np.sqrt(3))) * r_vals**(3/2)
    
    plt.figure(figsize=(8, 6))
    plt.plot(r_vals, h_c, 'k-', lw=2, label='Saddle-Node Bifurcation Line')
    plt.plot(r_vals, -h_c, 'k-', lw=2)
    
    plt.fill_between(r_vals, h_c, -h_c, color='lightblue', alpha=0.5)
    
    plt.text(1.0, 0, 'Bistable Region', horizontalalignment='center')
    plt.text(0.5, 0.3, 'Monostable Region', horizontalalignment='center')
    plt.text(0.5, -0.3, 'Monostable Region', horizontalalignment='center')
    
    plt.title('Cusp Bifurcation Diagram')
    plt.xlabel('$r$')
    plt.ylabel('$h$')
    plt.axhline(0, color='grey', linestyle='--')
    plt.axvline(0, color='grey', linestyle='--')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    # Generate combined figure and save
    plot_all_in_one(save_fig=True)