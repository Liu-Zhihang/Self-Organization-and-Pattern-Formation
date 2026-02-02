import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --- Set dark background style ---
plt.style.use('dark_background')

# Define the system of ODEs from the lecture 
# Using the equation u_dot = u - u^2*v (consistent with board drawing)
def system(t, Z):
    """Dynamical system for the nullcline example."""
    u, v = Z
    # Ensure u is not exactly zero for v=1/u calculation
    eps = 1e-9
    if abs(u) < eps:
        u_dot = 0.0 # On the u=0 nullcline
    else:
        u_dot = u - (u**2) * v
    v_dot = u - v
    return [u_dot, v_dot]

# Create a grid for the phase portrait
u_range = np.linspace(-2, 4, 20)
v_range = np.linspace(-1, 3, 20)
U, V = np.meshgrid(u_range, v_range)

# Calculate the vector field (dU, dV) at each grid point
dU, dV = np.zeros_like(U), np.zeros_like(V)
ni, nj = U.shape
for i in range(ni):
    for j in range(nj):
        dU[i,j], dV[i,j] = system(0, [U[i,j], V[i,j]])

# --- Create the plot ---
fig, ax = plt.subplots(figsize=(10, 8))
ax.set_title('Phase Portrait with Nullclines (Lecture Example)', color='white')
ax.set_xlabel('u', color='white')
ax.set_ylabel('v', color='white')

# 1. Plot the vector field using streamplot
# Streamplot shows the flow direction
ax.streamplot(U, V, dU, dV, color='gray', linewidth=0.7, density=1.5, broken_streamlines=False)

# 2. Plot the Nullclines
# v-NC (v_dot = 0): v = u
v_nc1 = u_range
ax.plot(u_range, v_nc1, 'cyan', linestyle='-', linewidth=2, label=r'v-NC: $\dot{v}=0$ ($v=u$)')

# u-NCs (u_dot = 0): u = 0 and v = 1/u
# u=0 (the y-axis)
ax.axvline(0, color='magenta', linestyle='-', linewidth=2, label=r'u-NC: $\dot{u}=0$ ($u=0$)')

# v = 1/u (plot in two parts)
u_nc2_pos = np.linspace(0.1, 4, 100) # Avoid u=0
v_nc2_pos = 1 / u_nc2_pos
ax.plot(u_nc2_pos, v_nc2_pos, 'magenta', linestyle='--', linewidth=2, label=r'u-NC: $\dot{u}=0$ ($v=1/u$)')

u_nc2_neg = np.linspace(-2, -0.1, 100) # Avoid u=0
v_nc2_neg = 1 / u_nc2_neg
ax.plot(u_nc2_neg, v_nc2_neg, 'magenta', linestyle='--', linewidth=2)

# 3. Find and plot fixed points (intersections of NCs)
# (0,0) is an intersection
ax.plot(0, 0, 'wo', markersize=8, markeredgecolor='black', label='Fixed Point (0, 0)')
# (1,1) is an intersection (v=u and v=1/u)
ax.plot(1, 1, 'wo', markersize=8, markeredgecolor='black', label='Fixed Point (1, 1)')
# (-1,-1) is an intersection
ax.plot(-1, -1, 'wo', markersize=8, markeredgecolor='black', label='Fixed Point (-1, -1)')

# 4. Plot example trajectories using solve_ivp
t_span = [0, 10]
sol1 = solve_ivp(system, t_span, [0.8, 0.5], dense_output=True, method='RK45')
ax.plot(sol1.y[0], sol1.y[1], 'lime', linewidth=2.5, label='Trajectory 1')

sol2 = solve_ivp(system, t_span, [1.5, 2.5], dense_output=True, method='RK45')
ax.plot(sol2.y[0], sol2.y[1], 'yellow', linewidth=2.5, label='Trajectory 2')

leg = ax.legend()
for text in leg.get_texts():
    text.set_color('white')
    
ax.set_ylim(-1, 3)
ax.set_xlim(-2, 4)
ax.grid(True, color='gray', linestyle='--', alpha=0.5)
ax.tick_params(colors='white')
for spine in ax.spines.values():
    spine.set_edgecolor('white')
plt.savefig('nullclines_and_phase_portrait.png')