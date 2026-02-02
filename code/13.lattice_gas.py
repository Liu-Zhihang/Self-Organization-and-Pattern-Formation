# =========================
# Model: lattice gas (Ising-like) with Kawasaki dynamics (particle-conserving)
# sigma=1: solute (P), sigma=0: solvent (S)
# Energy: H = -J * sum_<ij> s_i s_j with s_i = 2*sigma_i-1 in {+1,-1}
# Mapping to chi: chi ~ 2 z J / (k_B T)  (up to conventions); here we use a monotone proxy J(T)=J0*chi/chi_c
# =========================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches 
plt.style.use("dark_background")

def flory_huggins_f(phi, chi):
    eps = 1e-12
    phi = np.clip(phi, eps, 1 - eps)
    return phi*np.log(phi) + (1-phi)*np.log(1-phi) + chi*phi*(1-phi)

def to_spin(sigma):
    return 2*sigma - 1

def neighbor_sum(s, i, j):
    # periodic boundary conditions
    n = s.shape[0]
    return s[(i-1)%n, j] + s[(i+1)%n, j] + s[i, (j-1)%n] + s[i, (j+1)%n]

def kawasaki_sweep(sigma, betaJ, rng):
    # One Monte Carlo sweep: attempt N^2 swaps of unlike neighbors
    n = sigma.shape[0]
    for _ in range(n*n):
        i = rng.integers(0, n)
        j = rng.integers(0, n)
        # pick a random neighbor
        di, dj = [(1,0),(-1,0),(0,1),(0,-1)][rng.integers(0,4)]
        i2, j2 = (i+di)%n, (j+dj)%n
        if sigma[i,j] == sigma[i2,j2]:
            continue
        # compute ΔE for swapping sigma(i,j) and sigma(i2,j2) using spin form
        s = to_spin(sigma)
        si, sj = s[i,j], s[i2,j2]
        # local fields excluding the pair bond double-count issues handled by exact local diff
        hi = neighbor_sum(s, i, j) - sj
        hj = neighbor_sum(s, i2, j2) - si
        # Energy contribution for site i (excluding bond to j2): -J si hi ; after swap si' = sj -> -J sj hi
        # similarly for site j2: -J sj hj ; after swap sj' = si -> -J si hj
        # Bond between i and j2: -J si sj ; after swap -J sj si (same) -> cancels
        dE = - (sj*hi + si*hj) + (si*hi + sj*hj)  # multiplied by J later via betaJ
        # Metropolis accept
        if dE <= 0 or rng.random() < np.exp(-betaJ * dE):
            sigma[i,j], sigma[i2,j2] = sigma[i2,j2], sigma[i,j]
    return sigma

N = 64                 # lattice size N x N
phi0 = 0.45             # overall solute fraction (conserved)
seed = 42
rng = np.random.default_rng(seed)

# initialize sigma with fixed composition
phi_flat = np.zeros(N*N, dtype=np.int8)
phi_flat[:int(phi0*N*N)] = 1
rng.shuffle(phi_flat)
sigma = phi_flat.reshape(N, N)

# chi schedule for animation (controls free-energy landscape + interaction strength)
chi_c = 2.0
chi_values = np.linspace(0.8, 3.2, 120)

# link chi -> betaJ (monotone) : choose a scale so that around chi_c we see coarsening
J0 = 0.35
betaJ_values = J0 * (chi_values / chi_c)

# per-frame MC sweeps (increase for smoother coarsening, but slower)
mc_sweeps_per_frame = 4

fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 5), facecolor="black")
axL.set_title("Lattice Gas: Kawasaki Dynamics (Conserved)")
axL.axis("off")

solvent_color = '#ADD8E6'  # Light blue
solute_color = '#90EE90'   # Light green

# Function to draw lattice points (replaces draw_lattice with circles)
def draw_lattice(sigma):
    grid_size = sigma.shape[0]
    spacing = 0.6
    start_x = 1
    start_y = 1.5
    
    # Clear previous artists
    axL.clear()
    axL.set_title("Lattice Gas: Kawasaki Dynamics (Conserved)")
    axL.axis("off")
    
    # Plot points for each site
    x_coords = []
    y_coords = []
    colors = []
    
    for i in range(grid_size):
        for j in range(grid_size):
            x = start_x + i * spacing
            y = start_y + j * spacing
            state = sigma[i, j]
            color = solute_color if state == 1 else solvent_color
            x_coords.append(x)
            y_coords.append(y)
            colors.append(color)
    
    # Use scatter plot for efficiency and simplicity
    axL.scatter(x_coords, y_coords, c=colors, s=40, edgecolors='white', linewidth=0.5)

draw_lattice(sigma)

phi_grid = np.linspace(0.001, 0.999, 700)
lineR, = axR.plot([], [], lw=2)
dotR, = axR.plot([], [], marker="o", markersize=6, linestyle="")
axR.set_xlim(0, 1)
axR.set_ylim(-0.8, 0.2)
axR.set_xlabel("Volume Fraction $\\phi$")
axR.set_ylabel("Dimensionless Free Energy $f(\\phi)/k_BT$")
titleR = axR.set_title("")
axR.grid(alpha=0.15)

def update(frame):
    global sigma
    chi = float(chi_values[frame])
    betaJ = float(betaJ_values[frame])
    # evolve lattice configuration at current interaction strength
    for _ in range(mc_sweeps_per_frame):
        sigma = kawasaki_sweep(sigma, betaJ, rng)
    
    axL.clear()
    axL.set_title("Lattice Gas: Kawasaki Dynamics (Conserved)")
    axL.axis("off")
    draw_lattice(sigma)
    
    # free energy curve (shifted) at current chi
    f = flory_huggins_f(phi_grid, chi)
    lineR.set_data(phi_grid, f)
    # mark the conserved average composition on the free-energy curve
    f0 = flory_huggins_f(np.array([phi0]), chi)[0] 
    dotR.set_data([phi0], [f0])
    titleR.set_text(f"Free Energy Landscape: $\\chi = {chi:.2f}$ (right dot is system average $\\phi$)")

    return lineR, dotR, titleR  

ani = animation.FuncAnimation(fig, update, frames=len(chi_values), interval=80, blit=False)
ani.save("lattice_gas_vs_free_energy.gif", writer="pillow", fps=12)
plt.close(fig)
print("Saved: lattice_gas_vs_free_energy.gif")

