import numpy as np
import matplotlib.pyplot as plt

def plot_saddle_node_bifurcation():
    """
    Generate and display the bifurcation diagram and potential landscape for a saddle-node bifurcation.
    Uses 'dark_background' style.
    """
    # --- Set dark background theme ---
    plt.style.use('dark_background')

    # --- Create figure and subplots ---
    fig = plt.figure(figsize=(14, 6))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    
    # --- Figure 1: Bifurcation Diagram (u* vs mu) ---
    # du/dt = -mu + u^2  =>  u* = +/- sqrt(mu)
    mu_vals = np.linspace(0, 2, 400)
    u_stable = -np.sqrt(mu_vals)
    u_unstable = np.sqrt(mu_vals)
    
    # Plot stable branch (cyan)
    ax1.plot(mu_vals, u_stable, 'c-', linewidth=2, label='Stable Fixed Point')
    # Plot unstable branch (magenta)
    ax1.plot(mu_vals, u_unstable, 'm--', linewidth=2, label='Unstable Fixed Point')
    # Plot line for no fixed points when mu < 0 (gray)
    ax1.plot(np.linspace(-2, 0, 100), np.zeros(100), 'gray', linewidth=0.5, alpha=0.7)
    
    ax1.set_xlabel('Control Parameter $\\mu$')
    ax1.set_ylabel('Fixed Points $u^*$')
    ax1.set_title('Saddle-Node Bifurcation Diagram')
    ax1.axhline(0, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    ax1.axvline(0, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    
    # Set legend text color to white
    legend1 = ax1.legend()
    for text in legend1.get_texts():
        text.set_color('white')
        
    ax1.grid(True, linestyle='--', alpha=0.3)
    ax1.set_xlim(-2, 2)
    ax1.set_ylim(-2, 2)

    # --- Figure 2: Potential Landscape V(u) = mu*u - (1/3)*u^3 ---
    u_range = np.linspace(-2.5, 2.5, 500)
    mu_cases = [-1, 0, 1]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Blue, Orange, Green
    labels = [f'$\\mu = -1$ (No Fixed Points)', f'$\\mu = 0$ (Bifurcation Point)', f'$\\mu = 1$ (Two Fixed Points)']
    
    for mu, color, label in zip(mu_cases, colors, labels):
        V = mu * u_range - (1/3) * u_range**3
        ax2.plot(u_range, V, color=color, linewidth=2.5, label=label)
        
        # Mark fixed points
        if mu > 0:
            u_s = -np.sqrt(mu)  # Stable point
            u_u = np.sqrt(mu)   # Unstable point
            V_s = mu * u_s - (1/3) * u_s**3
            V_u = mu * u_u - (1/3) * u_u**3
            # Stable point (valley) - solid cyan
            ax2.plot(u_s, V_s, 'o', color='cyan', markersize=8)
            # Unstable point (peak) - hollow magenta
            ax2.plot(u_u, V_u, 'o', markerfacecolor='none', markeredgecolor='magenta', markersize=8)
        elif mu == 0:
            V_0 = 0
            # Inflection point - yellow
            ax2.plot(0, V_0, 'o', color='yellow', markersize=8)

    ax2.set_xlabel('State Variable $u$')
    ax2.set_ylabel('Potential Energy $V(u)$')
    ax2.set_title('Evolution of Potential Landscape with $\\mu$')
    
    # Set legend text color to white
    legend2 = ax2.legend()
    for text in legend2.get_texts():
        text.set_color('white')
        
    ax2.grid(True, linestyle='--', alpha=0.3)
    ax2.set_ylim(-1.5, 1.5)
    ax2.axhline(0, color='gray', linestyle=':', linewidth=1, alpha=0.7)

    plt.tight_layout()
    plt.savefig('saddle_node_bifurcation.png', dpi=300, bbox_inches='tight')
    plt.show()

# --- Run Simulation ---
if __name__ == "__main__":
    plot_saddle_node_bifurcation()
