import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def replicator_equation(t, x, b, c):
    """
    Replicator equation for the Prisoner's Dilemma.
    
    Parameters:
    t (float): Time (required by scipy.integrate).
    x (list): A list containing the fraction of cooperators, e.g., [0.5].
    b (float): Benefit from cooperation.
    c (float): Cost of cooperation.
    
    Returns:
    list: A list containing the rate of change of the cooperator fraction, [dx/dt].
    """
    x_val = x[0]
    if not (0 <= x_val <= 1):
        return [0]
    
    # Calculate fitnesses as derived in the notes
    f_coop = (b - c) * x_val - c * (1 - x_val)
    f_def = b * x_val
    
    # Replicator equation in the form: dx/dt = x * (1 - x) * (f_coop - f_def)
    dxdt = x_val * (1 - x_val) * (f_coop - f_def)
    
    return [dxdt]

# --- Set dark background style ---
plt.style.use('dark_background')

# Set game parameters (Prisoner's Dilemma condition: b > c > 0)
b = 3.0  # Benefit
c = 1.0  # Cost

# Set simulation parameters
t_span = [0, 10]
initial_conditions = [0.1, 0.3, 0.5, 0.7, 0.9]

# Create the plot
plt.figure(figsize=(10, 6))

# Solve and plot for each initial condition
for x0 in initial_conditions:
    sol = solve_ivp(
        replicator_equation, 
        t_span, 
        [x0], 
        args=(b, c), 
        dense_output=True,
        t_eval=np.linspace(t_span[0], t_span[1], 200)
    )
    plt.plot(sol.t, sol.y[0], lw=2.5, label=f'$x(0) = {x0}$')

# --- Style the plot for dark theme ---
plt.title("Evolution of Cooperation in Prisoner's Dilemma", fontsize=16)
plt.xlabel('Time $t$', fontsize=14)
plt.ylabel('Fraction of Cooperators $x(t)$', fontsize=14)

# Highlight the fixed points
plt.axhline(0, color='magenta', linestyle='--', lw=2, label='Stable Fixed Point (All Defectors)')
plt.axhline(1, color='gray', linestyle=':', lw=2, label='Unstable Fixed Point (All Cooperators)')

plt.grid(True, linestyle=':', alpha=0.5, color='gray')
plt.legend(fontsize=12)
plt.ylim(-0.05, 1.05)
plt.show()
