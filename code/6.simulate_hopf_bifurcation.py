import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --- Set dark background style ---
plt.style.use('dark_background')

# Define the Hopf bifurcation normal form equations 
def hopf_system(t, Z, mu, rho, omega, lambda_):
    """Hopf bifurcation normal form in Cartesian coordinates."""
    u, v = Z
    c_sq = u**2 + v**2
    u_dot = mu * u - omega * v + (rho * u - lambda_ * v) * c_sq
    v_dot = mu * v + omega * u + (rho * v + lambda_ * u) * c_sq
    return [u_dot, v_dot]

# --- Parameters for Simulation ---
omega = 1.0     # Linear frequency
rho = -1.0      # Coefficient for supercritical (rho < 0)
lambda_ = 0.0   # Non-linear frequency correction (set to 0 for simplicity)
t_span = [0, 40]  # Simulation time (increased for better viz)
y0_in = [0.1, 0]   # Initial condition inside the limit cycle
y0_out = [2.0, 0]  # Initial condition outside the limit cycle

# --- Case 1: mu = -0.1 (Stable Spiral) ---
mu_neg = -0.1
sol_neg = solve_ivp(hopf_system, t_span, y0_out, 
                    args=(mu_neg, rho, omega, lambda_), 
                    dense_output=True, method='RK45')

# --- Case 2: mu = +0.1 (Limit Cycle) ---
mu_pos = 0.1
# Analytical limit cycle radius: c = sqrt(-mu/rho)
c_star = np.sqrt(-mu_pos / rho)
# Simulate trajectory starting inside
sol_pos_in = solve_ivp(hopf_system, t_span, y0_in, 
                       args=(mu_pos, rho, omega, lambda_), 
                       dense_output=True, method='RK45')
# Simulate trajectory starting outside
sol_pos_out = solve_ivp(hopf_system, t_span, y0_out, 
                        args=(mu_pos, rho, omega, lambda_), 
                        dense_output=True, method='RK45')

# --- Plotting ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7), subplot_kw={'aspect': 'equal'})
fig.suptitle(r'Supercritical Hopf Bifurcation ($\rho = -1$)', color='white')

# Plot 1: mu = -0.1 (Stable Spiral)
ax1.set_title(rf'Case 1: $\mu = {mu_neg}$ (Stable Spiral)', color='white')
ax1.plot(sol_neg.y[0], sol_neg.y[1], 'cyan')
ax1.plot(sol_neg.y[0][0], sol_neg.y[1][0], 'co', label='Start (y0=[2,0])') # Start point
ax1.plot(0, 0, 'ro', label='Stable Fixed Point') # End point
ax1.set_xlabel('u')
ax1.set_ylabel('v')
leg1 = ax1.legend()
for text in leg1.get_texts(): text.set_color('white')
ax1.grid(True, color='gray', linestyle='--', alpha=0.5)
ax1.tick_params(colors='white')
for spine in ax1.spines.values(): spine.set_edgecolor('white')

# Plot 2: mu = +0.1 (Stable Limit Cycle)
ax2.set_title(rf'Case 2: $\mu = {mu_pos}$ (Stable Limit Cycle)', color='white')
# Plot the trajectories
ax2.plot(sol_pos_in.y[0], sol_pos_in.y[1], 'cyan', label='Start inside ($y_0=[0.1, 0]$)')
ax2.plot(sol_pos_out.y[0], sol_pos_out.y[1], 'lime', label='Start outside ($y_0=[2.0, 0]$)')
# Plot the analytical limit cycle
theta = np.linspace(0, 2*np.pi, 100)
ax2.plot(c_star * np.cos(theta), c_star * np.sin(theta), 'r--', linewidth=2, label=rf'Limit Cycle (c={c_star:.3f})')
ax2.plot(0, 0, 'rx', label='Unstable Fixed Point')
ax2.set_xlabel('u')
ax2.set_ylabel('v')
leg2 = ax2.legend()
for text in leg2.get_texts(): text.set_color('white')
ax2.grid(True, color='gray', linestyle='--', alpha=0.5)
ax2.tick_params(colors='white')
for spine in ax2.spines.values(): spine.set_edgecolor('white')

plt.savefig('6.simulate_hopf_bifurcation.png', bbox_inches='tight')
