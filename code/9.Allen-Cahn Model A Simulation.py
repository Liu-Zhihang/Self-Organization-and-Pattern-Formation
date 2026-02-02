import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap

# Use dark background style to match lecture slides
plt.style.use('dark_background')

# ==========================================
# 1. Physics & Simulation Parameters
# ==========================================
# Spatial parameters
N = 128                 # Grid size (N x N)
L = 100.0               # Physical size of the system
dx = L / N              # Spatial resolution
x = np.linspace(0, L, N)
y = np.linspace(0, L, N)

# Time parameters
dt = 0.01               # Time step
total_steps = 4000      # Total simulation steps
plot_interval = 50      # Frame update interval

# Ginzburg-Landau Model Parameters
kappa = 1.0             # Stiffness / Diffusion coefficient
r = 1.0                 # Control parameter (r > 0 for double well potential)
u = 1.0                 # Nonlinear coefficient (Quartic term)

# Stability Check (CFL Condition)
stability_limit = dx**2 / (4 * kappa)
print(f"Stability Check: dt = {dt}, Max Stable dt = {stability_limit:.4f}")
if dt > stability_limit:
    print("Warning: dt exceeds stability limit! Simulation may diverge.")

# ==========================================
# 2. Initial Condition (Quench)
# ==========================================
# Simulate a quench from high T (phi ~ 0) to low T.
# Initial state is small random noise around 0.
np.random.seed(42)  # Fixed seed for reproducibility
noise_strength = 0.1
phi = noise_strength * (2 * np.random.rand(N, N) - 1)  # Uniform noise in [-0.1, 0.1]

# ==========================================
# 3. Numerical Solver Kernel
# ==========================================
def laplacian_pbc(f, dx):
    """
    Compute 2D discrete Laplacian with Periodic Boundary Conditions (PBC).
    Using np.roll for efficient neighbor indexing.
    """
    left  = np.roll(f, -1, axis=0)
    right = np.roll(f,  1, axis=0)
    up    = np.roll(f, -1, axis=1)
    down  = np.roll(f,  1, axis=1)
    
    return (left + right + up + down - 4 * f) / (dx**2)

def update_step(phi, dt, dx, kappa, r, u):
    """
    Execute one Forward Euler time step.
    Equation: dphi/dt = kappa * grad^2 phi + (r*phi - u*phi^3)
    """
    # Diffusion term (Elastic relaxation)
    lap = laplacian_pbc(phi, dx)
    
    # Reaction term (Local potential drive)
    # Derivative f'(phi) = -r*phi + u*phi^3 => Reaction = -f' = r*phi - u*phi^3
    reaction = r * phi - u * phi**3
    
    # Update field
    phi_new = phi + dt * (kappa * lap + reaction)
    return phi_new

# ==========================================
# 4. Visualization Setup
# ==========================================
# Custom Colormap: Replicating the lecture's visual style
# -1 (Spin Down) -> Blue
#  0 (Interface) -> Green/Yellowish (Transition region)
# +1 (Spin Up)   -> Red
colors = [(0, 0, 1), (0.0, 0.8, 0.0), (1, 0, 0)]  # Blue -> Green -> Red
n_bins = 100
cmap_name = 'spin_decomposition_map'
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

fig, ax = plt.subplots(figsize=(6, 6))
# Set vmin/vmax slightly > 1 to avoid color saturation at stable points
im = ax.imshow(phi, cmap=cm, vmin=-1.2, vmax=1.2, interpolation='bilinear', origin='lower')
ax.set_title(f"Model A Dynamics (Step 0)")
ax.axis('off')

# Text for time step display
time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, color='white', fontsize=12)

def animate(frame_idx):
    global phi
    # Perform multiple physics steps per animation frame for speed
    steps_per_frame = 20
    for _ in range(steps_per_frame):
        phi = update_step(phi, dt, dx, kappa, r, u)
    
    im.set_array(phi)
    current_step = frame_idx * steps_per_frame
    time_text.set_text(f'Step: {current_step}')
    ax.set_title(f"Model A Dynamics (Step {current_step})")
    return [im, time_text]

# Create animation object
ani = animation.FuncAnimation(fig, animate, frames=total_steps//20, interval=30, blit=False)

# Save the animation as a video file
print("Saving simulation as video...")
ani.save('model_a_dynamics.mp4', writer='ffmpeg', fps=30, dpi=150)
print("Video saved as 'model_a_dynamics.mp4'.")

print("Starting simulation visualization...")
plt.show()
