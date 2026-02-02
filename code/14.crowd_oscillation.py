
from __future__ import annotations

import os
import math
from dataclasses import dataclass
from typing import Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# 0) Utility functions

def laplacian_periodic(field: np.ndarray) -> np.ndarray:
    """2D Laplacian with periodic boundaries: ∇²f ≈ neighbors - 4f"""
    return (
        -4.0 * field
        + np.roll(field, 1, axis=0)
        + np.roll(field, -1, axis=0)
        + np.roll(field, 1, axis=1)
        + np.roll(field, -1, axis=1)
    )

def moving_average_1d(x: np.ndarray, w: int) -> np.ndarray:
    if w <= 1:
        return x
    k = np.ones(w, dtype=float) / float(w)
    return np.convolve(x, k, mode="same")

def fft_spectrum_sv(vx: np.ndarray, vy: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Velocity spectrum: S_v(ω) = |FFT(vx)|^2 + |FFT(vy)|^2
    Returns: omega_sorted (rad/s), Sv_sorted
    """
    n = len(vx)
    fx = np.fft.fft(vx - vx.mean())
    fy = np.fft.fft(vy - vy.mean())
    Sv = (np.abs(fx) ** 2 + np.abs(fy) ** 2) / float(n)  # Simple normalization (shape only)

    freq = np.fft.fftfreq(n, d=dt)  # Hz
    omega = 2.0 * np.pi * freq      # rad/s

    idx = np.argsort(omega)
    return omega[idx], Sv[idx]

# 1) Part A: Cool GIF (2D toy spatial model)

@dataclass
class GifParams:
    # Grid
    N: int = 64
    dt: float = 0.005
    relax_steps: int = 1200
    steps: int = 3600

    frame_stride: int = 40
    quiver_stride: int = 5
    fps: int = 20

    # Local dynamics parameters (aligned with author's magnitude):contentReference[oaicite:4]{index=4}
    k: float = 0.027
    gamma_p: float = 18.0
    beta: float = 1.10
    epsilon: float = 0.025     # eta/gamma_p (author script same name eps):contentReference[oaicite:5]{index=5}
    sigma_p: float = 0.50      # toy noise

    # Toy spatial coupling (visually form clusters/large-scale flows)
    Du: float = 0.30
    Dp: float = 0.30

    # Particles
    show_particles: bool = True
    n_particles: int = 2000
    particle_speed_scale: float = 1.0

    seed: int = 0

def simulate_lattice_for_gif(prm: GifParams) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns: vx_frames, vy_frames, speed_frames
    shape: (n_frames, N, N)
    """
    rng = np.random.default_rng(prm.seed)
    N = prm.N
    dt = prm.dt

    # Initial small perturbation
    ux = 1e-3 * rng.standard_normal((N, N))
    uy = 1e-3 * rng.standard_normal((N, N))
    px = 1e-3 * rng.standard_normal((N, N))
    py = 1e-3 * rng.standard_normal((N, N))

    n_frames = prm.steps // prm.frame_stride
    vx_frames = np.zeros((n_frames, N, N), dtype=np.float32)
    vy_frames = np.zeros((n_frames, N, N), dtype=np.float32)
    sp_frames = np.zeros((n_frames, N, N), dtype=np.float32)

    frame_idx = 0
    total_steps = prm.relax_steps + prm.steps

    for step in range(total_steps):
        lap_ux = laplacian_periodic(ux)
        lap_uy = laplacian_periodic(uy)

        # v = du/dt
        vx = -prm.k * ux + px + prm.Du * lap_ux
        vy = -prm.k * uy + py + prm.Du * lap_uy

        # Update u
        ux = ux + dt * vx
        uy = uy + dt * vy

        # p dynamics: toy version of "windsock + nonlinear saturation + odd-friction"
        # We use an equivalent structure (consistent with your previous script) to ensure the visual effect of "collective oscillation".
        p2 = px * px + py * py
        pv = px * vx + py * vy
        odd_x = -p2 * vx + pv * px
        odd_y = -p2 * vy + pv * py

        lap_px = laplacian_periodic(px)
        lap_py = laplacian_periodic(py)

        dpx = (
            -prm.gamma_p * px
            + prm.gamma_p * prm.beta * (1.0 - prm.epsilon * p2) * vx
            + odd_x
            + prm.Dp * lap_px
        )
        dpy = (
            -prm.gamma_p * py
            + prm.gamma_p * prm.beta * (1.0 - prm.epsilon * p2) * vy
            + odd_y
            + prm.Dp * lap_py
        )

        px = px + dt * dpx + prm.sigma_p * math.sqrt(dt) * rng.standard_normal((N, N))
        py = py + dt * dpy + prm.sigma_p * math.sqrt(dt) * rng.standard_normal((N, N))

        if step >= prm.relax_steps:
            if (step - prm.relax_steps) % prm.frame_stride == 0 and frame_idx < n_frames:
                speed = np.sqrt(vx * vx + vy * vy)
                vx_frames[frame_idx] = vx.astype(np.float32)
                vy_frames[frame_idx] = vy.astype(np.float32)
                sp_frames[frame_idx] = speed.astype(np.float32)
                frame_idx += 1

    return vx_frames, vy_frames, sp_frames

def render_gif(vx_frames: np.ndarray, vy_frames: np.ndarray, sp_frames: np.ndarray, prm: GifParams, out_gif: str):
    n_frames, N, _ = vx_frames.shape

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_title("Crowd oscillations (toy spatial model)")
    ax.set_xlim(0, N)
    ax.set_ylim(0, N)
    ax.set_aspect("equal")

    vmax = float(np.percentile(sp_frames, 99))
    im = ax.imshow(
        sp_frames[0],
        origin="lower",
        interpolation="nearest",
        extent=(0, N, 0, N),
        vmin=0.0,
        vmax=vmax,
        cmap="RdYlBu_r",  # Blue→Red, close to Fig3f visual perception
    )
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("|v| (a.u.)")

    # quiver
    s = prm.quiver_stride
    x = np.arange(0.5, N, s)
    y = np.arange(0.5, N, s)
    X, Y = np.meshgrid(x, y)

    U0 = vx_frames[0][::s, ::s]
    V0 = vy_frames[0][::s, ::s]
    quiv = ax.quiver(X, Y, U0, V0, angles="xy", scale_units="xy", scale=1.5, width=0.004, color="k")

    # Particles (pure visual enhancement)
    scat = None
    if prm.show_particles:
        rng = np.random.default_rng(prm.seed + 123)
        px = rng.uniform(0, N, size=prm.n_particles)
        py = rng.uniform(0, N, size=prm.n_particles)
        scat = ax.scatter(px, py, s=2, alpha=0.6, c="k")

    dt_frame = prm.dt * prm.frame_stride

    def update(i: int):
        im.set_data(sp_frames[i])

        Ui = vx_frames[i][::s, ::s]
        Vi = vy_frames[i][::s, ::s]
        quiv.set_UVC(Ui, Vi)

        if scat is not None:
            pos = scat.get_offsets()
            gi = np.clip(pos[:, 1].astype(int), 0, N - 1)
            gj = np.clip(pos[:, 0].astype(int), 0, N - 1)
            vel_x = vx_frames[i][gi, gj]
            vel_y = vy_frames[i][gi, gj]
            pos[:, 0] = (pos[:, 0] + prm.particle_speed_scale * dt_frame * vel_x) % N
            pos[:, 1] = (pos[:, 1] + prm.particle_speed_scale * dt_frame * vel_y) % N
            scat.set_offsets(pos)
            return im, quiv, scat

        return im, quiv

    anim = FuncAnimation(fig, update, frames=n_frames, interval=1000 / prm.fps, blit=False)
    anim.save(out_gif, writer=PillowWriter(fps=prm.fps))
    plt.close(fig)

# 2) Part B: Fig.3f-like (0D mean-field scan beta)

@dataclass
class Fig3FParams:
    # Scan range (fewer points = faster)
    beta_ratio_min: float = 0.5
    beta_ratio_max: float = 2.5
    n_beta: int = 35

    # Mean-field parameters (aligned with author script magnitude):contentReference[oaicite:6]{index=6}
    k: float = 0.027
    gamma_p: float = 18.0
    epsilon: float = 0.025   # eta/gamma_p
    sigma: float = 2.0       # Noise amplitude (author script sigma):contentReference[oaicite:7]{index=7}

    # Numerical integration (significantly shortened for speedup)
    dt: float = 0.002
    n_relax: int = 20000
    n_steps: int = 60000
    n_run: int = 10

    # Frequency domain settings
    omega_max: float = 1.2        # rad/s (matching paper Fig3f vertical axis magnitude)
    smooth_window: int = 7        # Light smoothing of spectrum
    peak_exclude_omega: float = 0.10  # Avoid treating ω≈0 energy as peaks

    seed: int = 0

def mf_dynamics(Y: np.ndarray, k: float, beta: float, gamma_p: float, eps: float) -> np.ndarray:
    """
    We adopt the same mean-field structure as the author script (Eq.(S17)-style),
    directly using the expression form from the author's public code (but our own implementation).:contentReference[oaicite:8]{index=8}
    Y = [ux, uy, px, py]
    """
    ux, uy, px, py = Y
    prod_p = px * py
    norm_p2 = px * px + py * py
    dux = -k * ux + px
    duy = -k * uy + py

    dpx = (
        -gamma_p * (1.0 - beta + beta * eps * norm_p2) * px
        -k * beta * gamma_p * (1.0 - eps * norm_p2) * ux
        -k * prod_p * uy
        +k * (py * py) * ux
    )
    dpy = (
        -gamma_p * (1.0 - beta + beta * eps * norm_p2) * py
        -k * beta * gamma_p * (1.0 - eps * norm_p2) * uy
        -k * prod_p * ux
        +k * (px * px) * uy
    )
    return np.array([dux, duy, dpx, dpy], dtype=float)

def mf_integrate_one_run(prm: Fig3FParams, beta: float, rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray]:
    """
    Integrate mean-field using Euler–Maruyama (fast enough).
    Noise is added only to the p equation, aligning with the author script approach.:contentReference[oaicite:9]{index=9}
    Returns: vx(t), vy(t) (steady-state segment)
    """
    Y = 1e-10 * rng.standard_normal(4)  # Small perturbation

    dt = prm.dt
    noise_amp = prm.sigma

    # Relaxation phase
    for _ in range(prm.n_relax):
        dW = np.array([0.0, 0.0, rng.standard_normal() * noise_amp * math.sqrt(dt), rng.standard_normal() * noise_amp * math.sqrt(dt)])
        Y = Y + mf_dynamics(Y, prm.k, beta, prm.gamma_p, prm.epsilon) * dt + dW

    # Recording phase
    vx = np.zeros(prm.n_steps, dtype=float)
    vy = np.zeros(prm.n_steps, dtype=float)

    for i in range(prm.n_steps):
        dW = np.array([0.0, 0.0, rng.standard_normal() * noise_amp * math.sqrt(dt), rng.standard_normal() * noise_amp * math.sqrt(dt)])
        Y = Y + mf_dynamics(Y, prm.k, beta, prm.gamma_p, prm.epsilon) * dt + dW
        ux, uy, px, py = Y
        vx[i] = -prm.k * ux + px
        vy[i] = -prm.k * uy + py

    return vx, vy

def compute_fig3f_like(prm: Fig3FParams, out_png: str):
    # Critical beta_c formula: beta_c = 1 + k/gamma_p (author script IC condition also uses this):contentReference[oaicite:10]{index=10}
    beta_c = 1.0 + prm.k / prm.gamma_p

    beta_ratio = np.linspace(prm.beta_ratio_min, prm.beta_ratio_max, prm.n_beta)
    betas = beta_ratio * beta_c

    # Unified frequency axis: build once with FFTfreq
    n = prm.n_steps
    dt = prm.dt
    freq = np.fft.fftfreq(n, d=dt)
    omega = 2.0 * np.pi * freq

    pos = (omega >= 0.0) & (omega <= prm.omega_max)
    omega_pos = omega[pos]

    heat = np.zeros((prm.n_beta, omega_pos.size), dtype=float)
    omega0 = np.zeros(prm.n_beta, dtype=float)

    rng0 = np.random.default_rng(prm.seed)

    for i, be in enumerate(betas):
        Sv_acc = np.zeros_like(omega, dtype=float)

        # Average over multiple runs (but small n_run for speedup)
        for r in range(prm.n_run):
            # Use different sub-seeds for each run
            rng = np.random.default_rng(rng0.integers(0, 2**32 - 1))
            vx, vy = mf_integrate_one_run(prm, be, rng)
            w_sorted, Sv_sorted = fft_spectrum_sv(vx, vy, prm.dt)

            # Re-calculate unsorted Sv instead of interpolating back sorted values
            fx = np.fft.fft(vx - vx.mean())
            fy = np.fft.fft(vy - vy.mean())
            Sv = (np.abs(fx) ** 2 + np.abs(fy) ** 2) / float(n)
            Sv_acc += Sv

        Sv_mean = Sv_acc / float(prm.n_run)

        # Take positive frequencies and apply light smoothing
        Sv_pos = Sv_mean[pos]
        Sv_pos = moving_average_1d(Sv_pos, prm.smooth_window)

        # Normalize (to make heatmap resemble paper's "bright band" contrast; does not affect ω0)
        Sv_pos_norm = Sv_pos / (Sv_pos.max() + 1e-12)

        heat[i, :] = np.log10(Sv_pos_norm + 1e-12)

        # Extract dominant peak ω0: exclude near-zero peaks
        if be < beta_c:
            omega0[i] = 0.0
        else:
            mask_peak = omega_pos >= prm.peak_exclude_omega
            if np.any(mask_peak):
                j = np.argmax(Sv_pos_norm[mask_peak])
                omega0[i] = float(omega_pos[mask_peak][j])
            else:
                omega0[i] = 0.0

        print(f"[fig3f] {i+1}/{prm.n_beta} beta/beta_c={beta_ratio[i]:.2f} omega0={omega0[i]:.3f}")

    # --------- Plot Fig3f-like ---------
    fig, ax = plt.subplots(figsize=(7.2, 3.4))

    im = ax.imshow(
        heat.T,
        origin="lower",
        aspect="auto",
        extent=[beta_ratio[0], beta_ratio[-1], omega_pos[0], omega_pos[-1]],
        cmap="RdYlBu_r",   # Close to paper color scheme style
        vmin=np.percentile(heat, 3),
        vmax=np.percentile(heat, 99),
    )

    ax.set_xlabel(r'$\beta/\beta_c$')
    ax.set_ylabel(r'$\omega\ (\mathrm{rad}\ s^{-1})$')
    ax.set_title("Numerics (own): Fig.3f-like spectrum heatmap")

    # Threshold dashed line (beta/beta_c=1)
    ax.axvline(1.0, linestyle="--", linewidth=2.0, color="k")

    # Overlay ω0 black dots
    ax.plot(beta_ratio, omega0, "ko", markersize=4)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label(r'$\log_{10}\,S_v(\omega)$ (norm.)')

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close(fig)

# 3) Main: One-click output GIF + Fig3f

def main():
    out_gif = "crowd_oscillation.gif"
    out_fig3f = "fig3f_like_numerics.png"

    # ---- A) GIF ----
    gp = GifParams(
        N=64,
        dt=0.005,
        relax_steps=1200,
        steps=3600,
        beta=1.10,
        seed=0,
    )
    vx_f, vy_f, sp_f = simulate_lattice_for_gif(gp)
    render_gif(vx_f, vy_f, sp_f, gp, out_gif=out_gif)
    print(f"[done] saved {out_gif}")

    # ---- B) Fig3f-like ----
    fp = Fig3FParams(
        beta_ratio_min=0.5,
        beta_ratio_max=2.5,
        n_beta=35,
        dt=0.002,
        n_relax=20000,
        n_steps=60000,
        n_run=10,
        omega_max=1.2,
        smooth_window=7,
        peak_exclude_omega=0.10,
        seed=0,
    )
    compute_fig3f_like(fp, out_png=out_fig3f)
    print(f"[done] saved {out_fig3f}")

if __name__ == "__main__":
    main()
