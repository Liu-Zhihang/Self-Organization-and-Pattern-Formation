import numpy as np
import matplotlib.pyplot as plt

# Set black background style
plt.style.use('dark_background')
plt.rcParams['figure.facecolor'] = 'black'
plt.rcParams['axes.facecolor'] = 'black'
plt.rcParams['savefig.facecolor'] = 'black'
plt.rcParams['savefig.edgecolor'] = 'none'

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 12

def plot_nucleation_barrier():
    """
    Function 1: Visualization of Classical Nucleation Theory (CNT)
    Plotting the Gibbs Free Energy change (Delta F) vs Droplet Radius (R).
    Demonstrates the competition between Surface Cost (R^2) and Volume Gain (R^3).
    """
    # 1. Define Physical Parameters (Arbitrary Simulation Units)
    gamma = 1.0       # Surface Tension (Cost per unit area)
    delta_f = 0.5     # Bulk Free Energy Difference per unit volume (|f_new - f_old|)
    
    # 2. Define Radius Domain
    R = np.linspace(0, 10, 300)
    
    # 3. Calculate Energy Components
    # Surface Term (Positive Cost): proportional to Area ~ R^2
    F_surf = 4 * np.pi * R**2 * gamma
    
    # Volume Term (Negative Gain): proportional to Volume ~ R^3
    # Note: The gain reduces the free energy, hence the negative sign.
    F_vol = - (4/3) * np.pi * R**3 * delta_f
    
    # Total Free Energy Change
    F_total = F_surf + F_vol
    
    # 4. Determine Critical Radius (Rc) Analytical Solution
    # dF/dR = 8*pi*R*gamma - 4*pi*R^2*delta_f = 0  => Rc = 2*gamma/delta_f
    Rc = 2 * gamma / delta_f
    # Calculate the Barrier Height (Activation Energy)
    F_max = 4 * np.pi * Rc**2 * gamma - (4/3) * np.pi * Rc**3 * delta_f
    
    # 5. Plotting
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Plot components
    ax.plot(R, F_surf, '--', color='#e74c3c', alpha=0.6, linewidth=2, label='Surface Cost ($+4\pi R^2 \gamma$)')
    ax.plot(R, F_vol, '--', color='#3498db', alpha=0.6, linewidth=2, label='Volume Gain ($-4/3\pi R^3 |\Delta f|$)')
    
    # Plot Total Energy
    ax.plot(R, F_total, '-', color='#2c3e50', linewidth=3, label='Total Free Energy $\Delta F(R)$')
    
    # Annotate Critical Radius (Nucleation Barrier)
    ax.axvline(Rc, color='#27ae60', linestyle=':', alpha=0.8, linewidth=2)
    ax.scatter([Rc], [F_max], color='#27ae60', s=150, zorder=5, edgecolors='white')
    
    # Add descriptive text and arrows
    ax.annotate(f'Critical Radius $R_c={Rc:.1f}$', xy=(Rc, F_max), xytext=(Rc+1.5, F_max+20),
                arrowprops=dict(facecolor='#27ae60', shrink=0.05), fontsize=12, color='#27ae60')
    
    # Calculate appropriate y-values from the F_total array for annotation
    # Find the index in the R array closest to the desired radius
    stable_radius_idx = min(np.argmin(np.abs(R - (Rc+2))), len(F_total)-1)
    unstable_radius_idx = min(np.argmin(np.abs(R - (Rc/2))), len(F_total)-1)
    
    ax.annotate('Stable Growth Region\n(Spontaneous)', xy=(Rc+2, F_total[stable_radius_idx]), 
                xytext=(Rc+3, 50), arrowprops=dict(facecolor='black', arrowstyle='->'), fontsize=10)
    
    ax.annotate('Unstable Region\n(Dissolution)', xy=(Rc/2, F_total[unstable_radius_idx]), 
                xytext=(0.5, 50), arrowprops=dict(facecolor='black', arrowstyle='->'), fontsize=10)

    # Axis labels and Title
    ax.set_title('Classical Nucleation Theory: The Free Energy Barrier', fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Droplet Radius $R$', fontsize=14)
    ax.set_ylabel('Free Energy Change $\Delta F$', fontsize=14)
    ax.axhline(0, color='gray', linewidth=1)
    
    # Set limits for better visibility
    ax.set_ylim(-150, F_max * 1.5)
    ax.set_xlim(0, 10)
    ax.legend(fontsize=12, loc='lower right', frameon=True)
    
    plt.tight_layout()
    plt.savefig('nucleation_barrier.png', dpi=300, bbox_inches='tight')
    plt.close()

def plot_concentration_profile():
    """
    Function 2: Concentration Profile around a Single Droplet
    Visualizes the Gibbs-Thomson effect and the quasi-static diffusion field (1/r decay).
    Corresponds to the solution of Laplace equation: c(r) = c_inf + (c_R - c_inf) * (R/r)
    """
    # 1. Define Parameters
    R_droplet = 2.0         # Radius of the droplet
    c_inf = 0.4             # Far-field concentration (Supersaturated background)
    
    # Gibbs-Thomson Shift:
    # Small droplets have higher surface concentration. 
    # Let's assume the shift raises c_surface ABOVE c_inf for a dissolving droplet,
    # or BELOW c_inf for a growing droplet.
    # Scenario: Growing Droplet (Flux goes inward) -> c_inf > c_surface > c_eq
    # Scenario: Gibbs-Thomson Effect simply sets c_surface. 
    # Let's visualize a case where c_surface is elevated due to curvature.
    c_eq_flat = 0.3         # Equilibrium concentration for flat interface
    capillary_length = 0.5  # l_gamma
    
    # c_surface (Gibbs-Thomson) = c_eq_flat * (1 + l_gamma / R)
    c_surface = c_eq_flat * (1 + capillary_length / R_droplet)
    c_in = 0.9              # Concentration inside the droplet (high density phase)

    # 2. Define Spatial Domain
    # Inside the droplet: 0 to R
    r_in = np.linspace(0, R_droplet, 100)
    # Outside the droplet: R to Far field
    r_out = np.linspace(R_droplet, 10.0, 400)
    
    # 3. Calculate Profiles
    # Inside: Constant profile (solution to Laplace eq regular at origin)
    c_profile_in = np.full_like(r_in, c_in)
    
    # Outside: 1/r decay behavior
    # c(r) = c_inf + (c_surface - c_inf) * (R / r)
    c_profile_out = c_inf + (c_surface - c_inf) * (R_droplet / r_out)
    
    # 4. Plotting
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Plot Inside
    ax.plot(r_in, c_profile_in, '-', color='#f39c12', linewidth=3, label='Droplet Interior (High Density)')
    
    # Plot Outside
    ax.plot(r_out, c_profile_out, '-', color='#8e44ad', linewidth=3, label='Matrix (Diffusion Field ~ $1/r$)')
    
    # Plot Interface Boundary
    ax.axvline(R_droplet, color='black', linestyle='--', linewidth=1.5, label='Interface Position $R$')
    
    # Annotations for key concentration levels
    # c_in
    ax.text(0.5, c_in + 0.02, '$c_{in}$ (Droplet Phase)', color='#f39c12', fontweight='bold')
    
    # c_surface (Gibbs-Thomson)
    ax.plot([R_droplet], [c_surface], 'o', color='red', zorder=10)
    ax.annotate(f'$c(R) = c_{{out}} + \delta c$\n(Gibbs-Thomson)', 
                xy=(R_droplet, c_surface), xytext=(R_droplet+1.5, c_surface+0.1),
                arrowprops=dict(facecolor='red', arrowstyle='->'), color='red')
    
    # c_inf
    ax.axhline(c_inf, color='gray', linestyle=':', linewidth=2, label='Far Field $c_{\infty}$')
    ax.text(9, c_inf - 0.04, '$c_{\infty}$', color='gray', fontsize=12)
    
    # Visualizing the Gradient (Flux)
    ax.annotate('Flux $J$ (Growth)', xy=(R_droplet+1, (c_surface+c_inf)/2), xytext=(R_droplet+3, (c_surface+c_inf)/2 + 0.05),
                arrowprops=dict(facecolor='blue', arrowstyle='->', lw=2), fontsize=12, color='blue')

    # Titles and Labels
    ax.set_title('Concentration Profile: Gibbs-Thomson Effect & Diffusion Field', fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Radial Distance $r$', fontsize=14)
    ax.set_ylabel('Concentration $c(r)$', fontsize=14)
    ax.legend(fontsize=12, loc='center right')
    ax.grid(True, linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig('concentration_profile.png', dpi=300, bbox_inches='tight')
    plt.close()

from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

def ch3d_simulate(
    N=64,
    L=64.0,
    M=1.0,
    r=1.0,
    u=1.0,
    kappa=1.0,
    dt=0.2,
    n_steps=1800,
    sample_steps=(0, 300, 600, 900, 1200, 1500, 1800),
    seed=42,
    mean_c=0.12,
    noise_amp=0.02,
):
    """
    3D Cahn–Hilliard (semi-implicit spectral), periodic boundary conditions

      ∂t c = M ∇² ( -r c + u c^3 - κ ∇² c )

    mean_c != 0 makes domains droplet-like (off-critical quench)
    """
    rng = np.random.default_rng(seed)
    dx = L / N

    c = mean_c + noise_amp * rng.standard_normal((N, N, N))

    k = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    kx, ky, kz = np.meshgrid(k, k, k, indexing="ij")
    k2 = kx**2 + ky**2 + kz**2
    k4 = k2**2

    denom = 1.0 - dt * (M * r * k2 - M * kappa * k4)
    denom[0, 0, 0] = 1.0

    snapshots = {}
    snapshots[0] = c.copy()

    for step in range(1, n_steps + 1):
        c_hat = np.fft.fftn(c)
        c3_hat = np.fft.fftn(c**3)

        numer = c_hat - dt * (M * u * k2) * c3_hat
        c_hat_new = numer / denom
        c_hat_new[0, 0, 0] = c_hat[0, 0, 0]  # exact mass conservation

        c = np.real(np.fft.ifftn(c_hat_new))

        if step in sample_steps:
            snapshots[step] = c.copy()

    return snapshots, dx

def save_ch3d_slices_gif(snapshots, dx, out_gif="ch3d_slices.gif"):
    """
    Make a GIF of central z-slices over time (pastel colormap on black)
    """
    steps = sorted(snapshots.keys())
    fields = [snapshots[s] for s in steps]
    N = fields[0].shape[0]
    z0 = N // 2

    slices = [f[:, :, z0].T for f in fields]  # transpose for nicer orientation

    vlim = np.percentile(np.abs(np.stack(slices)), 99)

    fig, ax = plt.subplots(figsize=(6, 6))
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")

    # pastel, publication-like on dark background
    cmap = plt.get_cmap("coolwarm")

    im = ax.imshow(slices[0], origin="lower", vmin=-vlim, vmax=vlim, cmap=cmap)
    ax.set_title(f"3D CH: central slice (step={steps[0]})", fontsize=12, color="white", pad=12)
    ax.set_axis_off()
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.ax.yaxis.set_tick_params(color="white")
    plt.setp(plt.getp(cb.ax.axes, "yticklabels"), color="white")

    def update(i):
        im.set_data(slices[i])
        ax.set_title(f"3D CH: central slice (step={steps[i]})", fontsize=12, color="white", pad=12)
        return (im,)

    anim = FuncAnimation(fig, update, frames=len(steps), interval=250, blit=False)
    anim.save(out_gif, writer=PillowWriter(fps=4))
    plt.close(fig)

def save_ch3d_droplet_rotate_gif(snapshots, dx, out_gif="ch3d_droplet_rotate.gif"):
    """
    Make a 3D rotating scatter GIF of the final field (droplet-like domains)
    """
    last_step = max(snapshots.keys())
    c = snapshots[last_step]
    N = c.shape[0]
    L = N * dx

    # threshold for one phase; choose above mean to show droplets
    thr = np.percentile(c, 70)
    idx = np.argwhere(c > thr)

    rng = np.random.default_rng(0)
    max_pts = 70000
    if idx.shape[0] > max_pts:
        sel = rng.choice(idx.shape[0], size=max_pts, replace=False)
        idx = idx[sel]

    xs = (idx[:, 0] + 0.5) * dx
    ys = (idx[:, 1] + 0.5) * dx
    zs = (idx[:, 2] + 0.5) * dx
    vals = c[idx[:, 0], idx[:, 1], idx[:, 2]]

    # pastel colormap on black
    cmap = plt.get_cmap("viridis")

    fig = plt.figure(figsize=(7.2, 6.2))
    fig.patch.set_facecolor("black")
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("black")

    sc = ax.scatter(xs, ys, zs, c=vals, s=1, alpha=0.35, cmap=cmap)
    ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_zlim(0, L)
    ax.set_xlabel("x", color="white")
    ax.set_ylabel("y", color="white")
    ax.set_zlabel("z", color="white")
    ax.set_title(f"3D Cahn–Hilliard droplet-like domains (step={last_step})", color="white", pad=12)

    # tick colors
    ax.tick_params(colors="white")

    cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.08)
    cb.ax.yaxis.set_tick_params(color="white")
    plt.setp(plt.getp(cb.ax.axes, "yticklabels"), color="white")

    # rotate animation (top-journal style: clean, minimal)
    def update(frame):
        azim = frame
        elev = 18 + 10 * np.sin(np.deg2rad(frame))
        ax.view_init(elev=elev, azim=azim)
        return (sc,)

    frames = list(range(0, 360, 6))
    anim = FuncAnimation(fig, update, frames=frames, interval=60, blit=False)
    anim.save(out_gif, writer=PillowWriter(fps=18))
    plt.close(fig)

if __name__ == "__main__":
    print("Generating original two figures (unchanged)...")
    plot_nucleation_barrier()
    plot_concentration_profile()

    print("Running 3D Cahn–Hilliard simulation and exporting GIFs...")
    snapshots, dx = ch3d_simulate(
        N=64,
        L=64.0,
        M=1.0,
        r=1.0,
        u=1.0,
        kappa=1.0,
        dt=0.2,
        n_steps=1800,
        sample_steps=(0, 300, 600, 900, 1200, 1500, 1800),
        seed=42,
        mean_c=0.12,     # off-critical to encourage droplet morphology
        noise_amp=0.02,
    )

    save_ch3d_slices_gif(snapshots, dx, out_gif="ch3d_slices.gif")
    save_ch3d_droplet_rotate_gif(snapshots, dx, out_gif="ch3d_droplet_rotate.gif")

    print("Done. Saved:")
    print(" - nucleation_barrier.png")
    print(" - concentration_profile.png")
    print(" - ch3d_slices.gif")
    print(" - ch3d_droplet_rotate.gif")
