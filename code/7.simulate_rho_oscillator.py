import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Use a dark theme to mimic the blackboard-style lecture notes
plt.style.use('dark_background')


def rho_2d_ode(t, y, mu, n):
    """
    Right-hand side of the 2D reduced RhoGTPase model.

    Parameters
    ----------
    t : float
        Time (not used explicitly, system is autonomous).
    y : array_like, shape (2,)
        State vector y = [M_T, M_D].
        M_T : active Rho on the membrane.
        M_D : inactive Rho on the membrane.
    mu : float
        Relative inactivation rate (strength of negative feedback).
    n : float
        Total protein mass (conservation: M_T + M_D + C_D = n).

    Returns
    -------
    dydt : list of float
        Time derivatives [dM_T/dt, dM_D/dt].
    """
    M_T, M_D = y

    # 2D ODE system (after eliminating C_D = n - M_T - M_D)
    dM_T = -mu * M_T + (1.0 + M_T**2) * M_D
    dM_D = (n - M_T - M_D) - (1.0 + M_T**2) * M_D

    return [dM_T, dM_D]


def simulate_rho(mu, n, y0, t_max, n_points=5000):
    """
    Integrate the Rho 2D ODE system for given parameters.

    Parameters
    ----------
    mu : float
        Relative inactivation rate.
    n : float
        Total protein mass.
    y0 : array_like, shape (2,)
        Initial condition [M_T(0), M_D(0)].
    t_max : float
        Final integration time.
    n_points : int, optional
        Number of time points for sampling the solution.

    Returns
    -------
    t : 1D ndarray
        Time points.
    sol : 2D ndarray
        Solution array with shape (2, len(t)).
        sol[0, :] : M_T(t)
        sol[1, :] : M_D(t)
    """
    t_span = (0.0, t_max)
    t_eval = np.linspace(t_span[0], t_span[1], n_points)

    sol = solve_ivp(
        rho_2d_ode,
        t_span,
        y0,
        t_eval=t_eval,
        args=(mu, n),
        rtol=1e-8,
        atol=1e-10
    )

    return sol.t, sol.y


def set_phase_limits(ax, M_T, M_D):
    """
    Set aesthetically pleasing limits and aspect ratio
    for a phase portrait panel.
    """
    mt_min, mt_max = M_T.min(), M_T.max()
    md_min, md_max = M_D.min(), M_D.max()

    # Add a small margin around the data range
    pad_x = 0.05 * (mt_max - mt_min) if mt_max > mt_min else 0.1
    pad_y = 0.05 * (md_max - md_min) if md_max > md_min else 0.1

    ax.set_xlim(mt_min - pad_x, mt_max + pad_x)
    ax.set_ylim(md_min - pad_y, md_max + pad_y)
    ax.set_aspect('equal', adjustable='box')


# -----------------------------
# Parameter choices and simulation
# -----------------------------
mu = 10.5              # fixed relative inactivation rate
n_stable = 15.0        # parameter in the stable-focus regime (tau < 0)
n_oscillatory = 25.0   # parameter in the oscillatory regime (tau > 0)

# Same initial condition for both parameter sets (inside the physical triangle)
y0 = [0.1, 0.1]

# Integrate the system in both regimes
t_stable, y_stable = simulate_rho(mu, n_stable, y0, t_max=200.0)
t_osc,    y_osc    = simulate_rho(mu, n_oscillatory, y0, t_max=400.0)

M_T_stable, M_D_stable = y_stable
M_T_osc,    M_D_osc    = y_osc

# For the oscillatory case, remove the initial transient and keep only
# the late-time motion on the limit cycle for the phase portrait
idx_transient = len(t_osc) // 2
M_T_osc_cycle = M_T_osc[idx_transient:]
M_D_osc_cycle = M_D_osc[idx_transient:]

# For the time series in the oscillatory regime, show only the last window_T
window_T = 100.0
t_max_osc = t_osc.max()
mask_osc_window = t_osc > (t_max_osc - window_T)

# -----------------------------
# Create a figure with phase portraits and time traces
# -----------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 8))  # no shared x-axis

ax_phase_stable = axes[0, 0]
ax_phase_osc    = axes[0, 1]
ax_time_stable  = axes[1, 0]
ax_time_osc     = axes[1, 1]

# --- (A) Phase portrait: stable focus ---
ax_phase_stable.plot(
    M_T_stable, M_D_stable,
    linewidth=2.0,
    label='trajectory'
)
ax_phase_stable.set_xlabel(r'$M_T$ (active membrane)')
ax_phase_stable.set_ylabel(r'$M_D$ (inactive membrane)')
ax_phase_stable.set_title(
    r'Stable focus: $\mu = 10.5,\ n = 15$',
    fontsize=11
)
ax_phase_stable.grid(alpha=0.3, linestyle=':')
set_phase_limits(ax_phase_stable, M_T_stable, M_D_stable)

# --- (B) Phase portrait: limit cycle (only late-time motion) ---
ax_phase_osc.plot(
    M_T_osc_cycle, M_D_osc_cycle,
    linewidth=2.0,
    label='trajectory'
)
ax_phase_osc.set_xlabel(r'$M_T$ (active membrane)')
ax_phase_osc.set_ylabel(r'$M_D$ (inactive membrane)')
ax_phase_osc.set_title(
    r'Limit cycle: $\mu = 10.5,\ n = 25$',
    fontsize=11
)
ax_phase_osc.grid(alpha=0.3, linestyle=':')
set_phase_limits(ax_phase_osc, M_T_osc_cycle, M_D_osc_cycle)

# --- (C) Time traces in the stable regime ---
ax_time_stable.plot(
    t_stable, M_T_stable,
    linewidth=1.8,
    label=r'$M_T(t)$'
)
ax_time_stable.plot(
    t_stable, M_D_stable,
    linewidth=1.2,
    linestyle='--',
    label=r'$M_D(t)$'
)
ax_time_stable.set_xlabel('time')
ax_time_stable.set_ylabel('concentration')
ax_time_stable.set_title(
    r'Time series (stable focus)',
    fontsize=11
)
ax_time_stable.legend(loc='upper right', fontsize=9)
ax_time_stable.grid(alpha=0.3, linestyle=':')

# --- (D) Time traces in the oscillatory regime (last window_T only) ---
ax_time_osc.plot(
    t_osc[mask_osc_window], M_T_osc[mask_osc_window],
    linewidth=1.8,
    label=r'$M_T(t)$'
)
ax_time_osc.plot(
    t_osc[mask_osc_window], M_D_osc[mask_osc_window],
    linewidth=1.2,
    linestyle='--',
    label=r'$M_D(t)$'
)
ax_time_osc.set_xlabel('time')
ax_time_osc.set_ylabel('concentration')
ax_time_osc.set_title(
    rf'Time series (limit cycle, last {window_T:.0f} time units)',
    fontsize=11
)
ax_time_osc.legend(loc='upper right', fontsize=9)
ax_time_osc.grid(alpha=0.3, linestyle=':')

# Global figure title (English caption)
fig.suptitle(
    "Figure 4.5  Numerical simulation of the 2D Rho oscillator\n"
    "Top: phase portraits. Bottom: time traces in stable vs. oscillatory regimes.",
    fontsize=12
)

# Make sure layout works well on a dark background
fig.tight_layout(rect=[0, 0.03, 1, 0.92])

# Save the figure with a dark-friendly background
plt.savefig(
    "rho_oscillator_phase_time.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
