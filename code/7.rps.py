import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D

# Set dark background style
plt.style.use('dark_background')

def rps_system_deriv(t, y, a_RS, a_SP, a_PR):
    """
    Defines the 3D ODE system for the Rock-Paper-Scissors replicator dynamics.
    y = [R, P, S]
    """
    R, P, S = y
    
    # Growth driven by winning - Decay driven by losing
    # d[R]/dt = a_RS * R * S - a_PR * R * P
    dR_dt = a_RS * R * S - a_PR * R * P
    
    # d[P]/dt = a_PR * P * R - a_SP * P * S
    dP_dt = a_PR * P * R - a_SP * P * S
    
    # d[S]/dt = a_SP * S * P - a_RS * S * R
    dS_dt = a_SP * S * P - a_RS * S * R
    
    return [dR_dt, dP_dt, dS_dt]

# --- Simulation Parameters ---
# Interaction rates (Set to 1.0 for a symmetric cycle)
a_RS = 1.0
a_SP = 1.0
a_PR = 1.0
params = (a_RS, a_SP, a_PR)

# Initial conditions (Must sum to 1.0 for population fractions)
# Starting slightly off-center to induce oscillations
y0 = [0.2, 0.3, 0.5] 
t_span = [0, 50]

# Solve the ODE
sol = solve_ivp(rps_system_deriv, t_span, y0, args=params, 
                dense_output=True, method='RK45', rtol=1e-9)

R_traj = sol.y[0]
P_traj = sol.y[1]
S_traj = sol.y[2]

# --- Plotting ---
fig = plt.figure(figsize=(14, 6))

# Plot 1: Time Series
ax1 = fig.add_subplot(1, 2, 1)
ax1.plot(sol.t, R_traj, 'r-', linewidth=2, label='Rock [R]')
ax1.plot(sol.t, P_traj, 'g-', linewidth=2, label='Paper [P]')
ax1.plot(sol.t, S_traj, 'b-', linewidth=2, label='Scissors [S]')
ax1.set_xlabel('Time', color='white')
ax1.set_ylabel('Population Fraction', color='white')
ax1.set_title('RPS Dynamics: Time Series', color='white')
ax1.legend()
ax1.grid(True, linestyle=':', alpha=0.4)

# Plot 2: Phase Space Trajectory (3D)
ax2 = fig.add_subplot(1, 2, 2, projection='3d')
# Plot the trajectory
ax2.plot(R_traj, P_traj, S_traj, 'w-', linewidth=1.5)
# Mark Start and End
ax2.scatter(R_traj[0], P_traj[0], S_traj[0], color='white', s=50, label='Start')
ax2.scatter(R_traj[-1], P_traj[-1], S_traj[-1], color='red', s=50, label='End')

# Labels and Style
ax2.set_xlabel('Rock', color='white')
ax2.set_ylabel('Paper', color='white')
ax2.set_zlabel('Scissors', color='white')
ax2.set_title('RPS Dynamics: Phase Portrait (Simplex)', color='white')

# Adjust 3D axis pane colors for dark theme visibility
ax2.w_xaxis.set_pane_color((0.2, 0.2, 0.2, 1.0))
ax2.w_yaxis.set_pane_color((0.2, 0.2, 0.2, 1.0))
ax2.w_zaxis.set_pane_color((0.2, 0.2, 0.2, 1.0))
ax2.grid(False) # Clean look
ax2.legend()

plt.suptitle(f'Evolutionary Game Theory: RPS Cycle ($a_{{ij}}=1.0$)', color='white', fontsize=14)
plt.tight_layout()
plt.savefig('rps_dynamics.png', bbox_inches='tight')
