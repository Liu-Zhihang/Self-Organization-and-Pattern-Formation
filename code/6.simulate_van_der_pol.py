import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --- Set dark background style ---
plt.style.use('dark_background')

# Van der Pol equation as a 1st order system
# Let u = x, v = x_dot
# Then u_dot = v
# And v_dot = -u - gamma * (u^2 - 1) * v = gamma * (1 - u^2) * v - u
def van_der_pol(t, Z, gamma):
    """Van der Pol oscillator 1st order system."""
    u, v = Z
    u_dot = v
    v_dot = gamma * (1 - u**2) * v - u
    return [u_dot, v_dot]

# Time span
t_span = [0, 100]
# Initial condition
y0 = [2.0, 0.0] # Start at x=2, x_dot=0

# --- Case 1: gamma = 0.5 (Weakly non-linear, near-harmonic) ---
gamma_weak = 0.5
sol_weak = solve_ivp(van_der_pol, t_span, y0, args=(gamma_weak,), 
                     dense_output=True, method='RK45')

# --- Case 2: gamma = 10 (Strongly non-linear, relaxation oscillation) ---
gamma_strong = 10.0
sol_strong = solve_ivp(van_der_pol, t_span, y0, args=(gamma_strong,), 
                       dense_output=True, method='RK45')

# --- Plotting Phase Portraits ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
fig.suptitle(r'Van der Pol Oscillator Phase Portraits ($x$ vs $\dot{x}$)', color='white')

# Plot 1: gamma = 0.5
ax1.set_title(rf'$\gamma = {gamma_weak}$ (Near-harmonic oscillation)', color='white')
ax1.plot(sol_weak.y[0], sol_weak.y[1], 'cyan')
ax1.set_xlabel('x (Position)', color='white')
ax1.set_ylabel(r'$\dot{x}$ (Velocity)', color='white')
ax1.axis('equal')
ax1.grid(True, color='gray', linestyle='--', alpha=0.5)
ax1.tick_params(colors='white')
for spine in ax1.spines.values(): spine.set_edgecolor('white')

# Plot 2: gamma = 10.0
ax2.set_title(rf'$\gamma = {gamma_strong}$ (Relaxation Oscillation)', color='white')
ax2.plot(sol_strong.y[0], sol_strong.y[1], 'lime')
ax2.set_xlabel('x (Position)', color='white')
ax2.set_ylabel(r'$\dot{x}$ (Velocity)', color='white')
ax2.grid(True, color='gray', linestyle='--', alpha=0.5)
ax2.tick_params(colors='white')
for spine in ax2.spines.values(): spine.set_edgecolor('white')

fig.savefig('van_der_pol_phase_portraits.png', dpi=300, bbox_inches='tight', facecolor='black')
plt.close(fig)  # Close the figure to free memory

# --- Plotting Time Series for gamma = 10 ---
fig_ts, ax_ts = plt.subplots(figsize=(12, 5))
ax_ts.set_title(rf'$\gamma = {gamma_strong}$ Time Series (Relaxation)', color='white')
# Use the solution from t=50 onwards to show the stable cycle
t_eval = np.linspace(50, 100, 500)
sol_strong_cycle = sol_strong.sol(t_eval)
ax_ts.plot(t_eval, sol_strong_cycle[0], 'lime', label='x(t) (Position)')
ax_ts.set_xlabel('Time', color='white')
ax_ts.set_ylabel('x', color='white')
leg = ax_ts.legend()
for text in leg.get_texts(): text.set_color('white')
ax_ts.grid(True, color='gray', linestyle='--', alpha=0.5)
ax_ts.tick_params(colors='white')
for spine in ax_ts.spines.values(): spine.set_edgecolor('white')


fig_ts.savefig('van_der_pol_time_series.png', dpi=300, bbox_inches='tight', facecolor='black')
