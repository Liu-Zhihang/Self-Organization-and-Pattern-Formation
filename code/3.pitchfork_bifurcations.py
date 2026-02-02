import numpy as np
import matplotlib.pyplot as plt

def plot_pitchfork_bifurcations():
    """
    Generates and displays plots for supercritical and subcritical 
    pitchfork bifurcations using a dark background.
    """
    
    # --- Use dark background style ---
    plt.style.use('dark_background')

    # --- Create figure and axes ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # --- 1. Supercritical Pitchfork Bifurcation ---
    ax1, ax2 = axes[0, 0], axes[0, 1]
    r_vals_super = np.linspace(-2, 2, 500)
    
    # --- ax1: Supercritical Bifurcation Diagram ---
    # u* = 0 (Trivial branch)
    # Stable for r < 0 (cyan)
    ax1.plot(r_vals_super[r_vals_super < 0], np.zeros_like(r_vals_super[r_vals_super < 0]), 'c-', linewidth=2.5, label='Stable')
    # Unstable for r > 0 (magenta)
    ax1.plot(r_vals_super[r_vals_super >= 0], np.zeros_like(r_vals_super[r_vals_super >= 0]), 'm--', linewidth=2.5, label='Unstable')
    
    # u* = +/- sqrt(r) (Non-trivial branches)
    r_positive = r_vals_super[r_vals_super >= 0]
    u_plus = np.sqrt(r_positive)
    u_minus = -np.sqrt(r_positive)
    ax1.plot(r_positive, u_plus, 'c-', linewidth=2.5)
    ax1.plot(r_positive, u_minus, 'c-', linewidth=2.5)
    
    ax1.set_title('Supercritical Pitchfork Diagram')
    ax1.set_xlabel('Control Parameter $r$')
    ax1.set_ylabel('Fixed Point $u^*$')
    ax1.grid(True, linestyle='--', alpha=0.3)
    legend1 = ax1.legend()
    for text in legend1.get_texts(): text.set_color('white')
    ax1.axhline(0, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    ax1.axvline(0, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    ax1.set_ylim(-2.5, 2.5)
    ax1.set_xlim(-2, 2)

    # --- ax2: Supercritical Potential V(u) = -r/2 * u^2 + 1/4 * u^4 ---
    u_range = np.linspace(-2.5, 2.5, 500)
    r_cases = [(-1, '#1f77b4', 'r = -1 (Single Well)'), 
               (0, '#ff7f0e', 'r = 0 (Flat Bottom)'), 
               (2, '#2ca02c', 'r = 2 (Double Well)')]
    
    for r_val, color, label in r_cases:
        V = -0.5 * r_val * u_range**2 + 0.25 * u_range**4
        ax2.plot(u_range, V, label=label, color=color, linewidth=2.5)
    ax2.set_title('Supercritical Potential Landscapes')
    ax2.set_xlabel('State Variable $u$')
    ax2.set_ylabel('Potential $V(u)$')
    ax2.grid(True, linestyle='--', alpha=0.3)
    legend2 = ax2.legend()
    for text in legend2.get_texts(): text.set_color('white')
    ax2.set_ylim(-1.5, 2)
    ax2.axhline(0, color='gray', linestyle=':', linewidth=1, alpha=0.5)

    # --- 2. Subcritical Pitchfork Bifurcation (stabilized) ---
    ax3, ax4 = axes[1, 0], axes[1, 1]
    
    # --- ax3: Subcritical Bifurcation Diagram ---
    # f(u,r) = ru + u^3 - u^5 = 0
    
    # Trivial branch: u = 0
    # Stable for r < 0 (cyan), Unstable for r > 0 (magenta)
    ax3.plot(r_vals_super[r_vals_super < 0], np.zeros_like(r_vals_super[r_vals_super < 0]), 'c-', linewidth=2.5, label='Stable')
    ax3.plot(r_vals_super[r_vals_super >= 0], np.zeros_like(r_vals_super[r_vals_super >= 0]), 'm--', linewidth=2.5, label='Unstable')

    # Non-trivial branches: u^4 - u^2 - r = 0 => u^2 = (1 +/- sqrt(1+4r))/2
    # Outer stable branch: u^2 = (1 + sqrt(1+4r))/2 (exists for r > -0.25)
    r_outer = np.linspace(-0.25, 1, 400)
    u_outer = np.sqrt((1 + np.sqrt(1 + 4 * r_outer)) / 2)
    ax3.plot(r_outer, u_outer, 'c-', linewidth=2.5)
    ax3.plot(r_outer, -u_outer, 'c-', linewidth=2.5)
    
    # Inner unstable branch: u^2 = (1 - sqrt(1+4r))/2 (exists for -0.25 < r < 0)
    r_inner = np.linspace(-0.25, 0, 100, endpoint=False)
    u_inner = np.sqrt((1 - np.sqrt(1 + 4 * r_inner)) / 2)
    ax3.plot(r_inner, u_inner, 'm--', linewidth=2.5)
    ax3.plot(r_inner, -u_inner, 'm--', linewidth=2.5)
    
    ax3.set_title('Subcritical Pitchfork Diagram (Stabilized)')
    ax3.set_xlabel('Control Parameter $r$')
    ax3.set_ylabel('Fixed Point $u^*$')
    ax3.grid(True, linestyle='--', alpha=0.3)
    legend3 = ax3.legend()
    for text in legend3.get_texts(): text.set_color('white')
    ax3.axhline(0, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    ax3.axvline(0, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    ax3.axvline(-0.25, color='yellow', linestyle=':', linewidth=1, alpha=0.5, label='Saddle-Node $r_s$')
    ax3.set_ylim(-2.5, 2.5)
    ax3.set_xlim(-1, 1)

    # --- ax4: Subcritical Potential V(u) = -r/2*u^2 - 1/4*u^4 + 1/6*u^6 ---
    u_range_sub = np.linspace(-2, 2, 500)
    # r_s = -0.25 is where the non-trivial branch is born
    r_cases_sub = [(-0.5, '#1f77b4', 'r = -0.5 (Bistable)'),
                   (-0.25, 'yellow', 'r = -0.25 (Saddle-Node)'),
                   (0.1, '#2ca02c', 'r = 0.1 (u=0 Unstable)')]
    
    for r_val, color, label in r_cases_sub:
        V_sub = -0.5 * r_val * u_range_sub**2 - 0.25 * u_range_sub**4 + (1/6) * u_range_sub**6
        ax4.plot(u_range_sub, V_sub, label=label, color=color, linewidth=2.5)
        
    ax4.set_title('Subcritical Potential Landscapes (Stabilized)')
    ax4.set_xlabel('State Variable $u$')
    ax4.set_ylabel('Potential $V(u)$')
    ax4.grid(True, linestyle='--', alpha=0.3)
    legend4 = ax4.legend()
    for text in legend4.get_texts(): text.set_color('white')
    ax4.set_ylim(-0.1, 0.1)
    ax4.axhline(0, color='gray', linestyle=':', linewidth=1, alpha=0.5)

    plt.tight_layout()
    plt.savefig('pitchfork_bifurcations.png', dpi=300, bbox_inches='tight')  # 保存图片
    plt.show()

# --- Run the simulation ---
if __name__ == "__main__":
    plot_pitchfork_bifurcations()
