import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ==========================================
# 1. Physics & Simulation Parameters
# ==========================================
N = 200           # Grid size
dx = 1.0          # Spatial resolution
kappa = 2.0       # Stiffness (Increased slightly to make interface smoother)
dt = 0.05         # Time step (Satisfies CFL: 0.05 < 1^2 / (4*2) = 0.125)

# Potential parameters
r = 1.0
u = 1.0
# Theoretical Interface Width: xi = sqrt(2*kappa/r) = sqrt(4) = 2.0
# Theoretical Critical Radius: R should be >> xi (e.g., R=40)

total_steps = 4000
steps_per_frame = 50

# ==========================================
# 2. Initialization: Spherical Droplet
# ==========================================
# We explicitly set up the "Droplet Scenario" from Section 2
# Background is phi = -1, Droplet is phi = +1
phi = -1.0 * np.ones((N, N))

# Create coordinate grid
Y, X = np.ogrid[:N, :N]
center = N // 2
R_initial = 60.0

# Create a sharp circular interface
mask = (X - center)**2 + (Y - center)**2 < R_initial**2
phi[mask] = 1.0

# Note: The simulation will naturally "smooth" this sharp edge 
# into a tanh profile within the first few steps.

# ==========================================
# 3. Solver Functions
# ==========================================
def laplacian(field, dx):
    """Discrete Laplacian with Periodic Boundary Conditions."""
    left = np.roll(field, -1, axis=0)
    right = np.roll(field, 1, axis=0)
    up = np.roll(field, -1, axis=1)
    down = np.roll(field, 1, axis=1)
    return (left + right + up + down - 4 * field) / dx**2

def free_energy_derivative(phi, r, u):
    """Derivative of the double-well potential."""
    return -r * phi + u * phi**3

# ==========================================
# 4. Visualization Setup
# ==========================================
fig, ax = plt.subplots(figsize=(6, 6))
# Use a diverging colormap to clearly see the interface (white/green region)
im = ax.imshow(phi, cmap='RdBu_r', vmin=-1.1, vmax=1.1, animated=True)

# Add a contour line to visualize the exact interface position (phi=0)
contour = ax.contour(phi, levels=[0], colors='yellow', linewidths=1)

time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, color='black', fontsize=12)
ax.set_title("Curvature-Driven Flow: Shrinking Droplet")
ax.axis('off')

def update(frame):
    global phi
    # Physics update loop
    for _ in range(steps_per_frame):
        lap = laplacian(phi, dx)
        reaction = -free_energy_derivative(phi, r, u)
        # Allen-Cahn Equation
        phi += dt * (kappa * lap + reaction)
    
    # Visualization update
    im.set_array(phi)
    
    # Update contour lines (remove old ones and draw new)
    for c in ax.collections:
        if c != im: c.remove()
    ax.contour(phi, levels=[0], colors='yellow', linewidths=1.5)
    
    current_time = (frame + 1) * steps_per_frame * dt
    time_text.set_text(f"Time: {current_time:.1f}")
    
    return im, time_text

# Run animation and save as video
ani = animation.FuncAnimation(fig, update, frames=200, interval=20, blit=False)
ani.save('droplet_shrinkage.mp4', writer='ffmpeg', fps=30)