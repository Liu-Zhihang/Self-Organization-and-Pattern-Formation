import numpy as np
import matplotlib.pyplot as plt

# --- Set dark background style ---
plt.style.use('dark_background')

# Define parameters
r_range = np.linspace(-2, 2, 400)
u = 1.0

# --- Calculate fixed points ---

# 1. For r < 0 (Ordered Phase)
r_neg = r_range[r_range < 0]
# The two new stable branches
phi_stable_pos = np.sqrt(-r_neg / u)
phi_stable_neg = -np.sqrt(-r_neg / u)
# The phi=0 branch becomes unstable
phi_unstable_zero = np.zeros_like(r_neg)

# 2. For r >= 0 (Disordered Phase)
r_pos = r_range[r_range >= 0]
# The phi=0 branch is stable
phi_stable_zero = np.zeros_like(r_pos)


# --- Create the plot ---
plt.figure(figsize=(10, 6))

# --- Plot the branches with visible colors ---

# Plot stable fixed points (solid bright lines)
plt.plot(r_neg, phi_stable_pos, 'cyan', lw=2.5, label='Stable Fixed Points ($\phi^* = \pm\sqrt{-r/u}$)')
plt.plot(r_neg, phi_stable_neg, 'cyan', lw=2.5)
plt.plot(r_pos, phi_stable_zero, 'cyan', lw=2.5)

# Plot unstable fixed point (dashed bright line)
plt.plot(r_neg, phi_unstable_zero, 'magenta', linestyle='--', lw=2, label='Unstable Fixed Point ($\phi^* = 0$)')

# --- Style the plot ---
plt.title('Supercritical Pitchfork Bifurcation in Ginzburg-Landau Model', fontsize=16)
plt.xlabel('Control Parameter $r \propto (T - T_c)$', fontsize=14)
plt.ylabel('Fixed Points $\phi^*$', fontsize=14)

# Add reference lines
plt.axhline(0, color='gray', linestyle='-', lw=1)
plt.axvline(0, color='gray', linestyle=':', lw=1, label='Bifurcation Point (r=0)')

plt.grid(True, linestyle=':', alpha=0.5, color='gray')
plt.legend(fontsize=12)
plt.ylim(-2, 2)
plt.xlim(-2, 2)
plt.show()