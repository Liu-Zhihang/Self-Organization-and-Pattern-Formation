import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Set dark background theme
plt.style.use('dark_background')

# ================= PART 1: Replicator Dynamics (ALVE) =================
def run_replicator_simulation():
    # Antisymmetric interaction matrix (from lecture Slide 5)
    A = np.array([[ 0,  2,  1,  0, -1],
                  [-2,  0,  2,  1, -3],
                  [-1, -2,  0, -1,  4],
                  [ 0, -1,  1,  0, -4],
                  [ 1,  3, -4,  4,  0]])
    S = A.shape[0]

    def alve_dynamics(t, x):
        # ALVE equation: dx_i/dt = x_i * (Ax)_i
        # Note: In numerical computation, x may slightly deviate from the simplex due to errors, but this does not affect short-term demonstration
        fitness = A @ x
        return x * fitness

    # Initial condition: random distribution
    np.random.seed(42)
    x0 = np.random.rand(S)
    x0 /= np.sum(x0)

    # Solve ODE
    t_span = (0, 15)
    sol = solve_ivp(alve_dynamics, t_span, x0, t_eval=np.linspace(0, 15, 500))

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Add a small epsilon to avoid log(0) error
    epsilon = 1e-10
    
    # Color cycle to distinguish different strategies
    colors = plt.cm.tab10(np.linspace(0, 1, S))

    for i in range(S):
        # Differentiate line styles based on final survival: solid for survivors, dashed for extinct
        is_survivor = sol.y[i, -1] > 1e-3
        style = '-' if is_survivor else '--'
        status = 'Condensate' if is_survivor else 'Depleted'
        label = f'Strategy {i+1} ({status})'
        
        ax.plot(sol.t, np.log(sol.y[i] + epsilon), 
                linestyle=style, label=label, linewidth=2, color=colors[i])
    
    ax.set_title('Replicator Dynamics: Condensation & Depletion', fontsize=16, color='white')
    ax.set_xlabel('Time', fontsize=14, color='white')
    ax.set_ylabel(r'$\ln(x_i)$', fontsize=14, color='white')
    ax.grid(True, alpha=0.2, linestyle='--')
    
    # Set legend text color
    legend = ax.legend(fontsize=10, loc='upper right')
    plt.setp(legend.get_texts(), color='white')
    
    plt.tight_layout()

    plt.savefig('replicator_dynamics.png') 
    plt.close() 

# ================= PART 2: 1D Diffusion Equation =================
def run_diffusion_simulation():
    L = 100.0
    Nx = 100
    D = 5.0
    T = 20.0
    
    dx = L / (Nx - 1)
    # Explicit difference stability condition: dt <= dx^2 / (2D)
    dt = 0.2 * dx**2 / D 
    
    x = np.linspace(0, L, Nx)
    
    # Initial condition: narrow Gaussian distribution at center
    u = np.exp(-((x - L/2)**2) / 20.0)
    
    # Store snapshots for plotting
    u_history = [u.copy()]
    t_snapshots = [0.0]
    
    steps = int(T / dt)
    plot_interval = steps // 10  # Record 10 time point snapshots
    
    for n in range(steps):
        u_new = u.copy()
        # Explicit finite difference scheme: u_new = u + alpha * (u_left - 2u + u_right)
        alpha = D * dt / dx**2
        u_new[1:-1] = u[1:-1] + alpha * (u[2:] - 2*u[1:-1] + u[:-2])
        
        # No flux boundary conditions (Neumann BCs): u[-1] = u[-2], u[0] = u[1]
        # This simulates particles bouncing off walls
        u_new[0] = u_new[1]
        u_new[-1] = u_new[-2]
        
        u = u_new
        
        if (n + 1) % plot_interval == 0:
            u_history.append(u.copy())
            t_snapshots.append((n + 1) * dt)

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Use gradient colors to show temporal evolution: purple (early) to yellow (late)
    colors = plt.cm.plasma(np.linspace(0, 1, len(u_history)))
    
    for i, u_snapshot in enumerate(u_history):
        label = f't={t_snapshots[i]:.1f}' if i == 0 or i == len(u_history)-1 else None
        ax.plot(x, u_snapshot, color=colors[i], lw=2, label=label)
    
    ax.set_title(f'1D Diffusion (No Flux BC, D={D})', fontsize=16, color='white')
    ax.set_xlabel('Position x', fontsize=14, color='white')
    ax.set_ylabel('Concentration C(x)', fontsize=14, color='white')
    
    # Add colorbar to indicate time progression
    sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma, norm=plt.Normalize(vmin=0, vmax=T))
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Time', color='white', fontsize=12)
    cbar.ax.yaxis.set_tick_params(color='white')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
    
    ax.grid(True, alpha=0.2, linestyle='--')
    legend = ax.legend(loc='upper right', fontsize=10)
    plt.setp(legend.get_texts(), color='white')
    
    plt.tight_layout()
    plt.savefig('diffusion_equation.png') 

if __name__ == "__main__":
    run_replicator_simulation()
    run_diffusion_simulation()
