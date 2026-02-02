import numpy as np
import matplotlib.pyplot as plt

def ginzburg_landau_potential(phi, r, u):
    """
    Calculate the Ginzburg-Landau free energy potential (v=0).
    """
    return 0.5 * r * phi**2 + 0.25 * u * phi**4

# --- Use dark background as a base ---
plt.style.use('dark_background')

# --- Create a single figure and axes with default dark background ---
fig, ax = plt.subplots(figsize=(10, 8))

# --- Parameters ---
u = 1.0
phi_range = np.linspace(-2.5, 2.5, 400) # Focus on the central region

# --- Plot the T > Tc (r > 0) case (White Curve) ---
r_pos = 1.0
f_pos = ginzburg_landau_potential(phi_range, r_pos, u)
# Plot in white, mimicking the professor's drawing
ax.plot(phi_range, f_pos, lw=2, color='white', linestyle='--')
# Add text annotation
ax.text(0.8, 2.0, '$T > T_c$ ($r > 0$)', color='white', fontsize=16)


# --- Plot the T < Tc (r < 0) case (Yellow Curve) ---
r_neg = -1.0
f_neg = ginzburg_landau_potential(phi_range, r_neg, u)
# Plot the potential curve in yellow
ax.plot(phi_range, f_neg, lw=3, color='yellow')
# Add text annotation
ax.text(1.5, -0.5, '$T < T_c$ ($r < 0$)', color='yellow', fontsize=16)

# --- Calculate minima for T < Tc ---
minima_phi = np.sqrt(-r_neg / u)
minima_f = ginzburg_landau_potential(minima_phi, r_neg, u)

# --- Add text annotations like the blackboard ---
# Add "f" (for the y-axis)
ax.plot([0, 0], [minima_f - 0.5, 8.0], color='white', lw=1)
ax.arrow(0, 8.0, 0, 0.5, color='white', head_width=0.08, head_length=0.3)
ax.text(-0.25, 7.5, '$f$', color='white', fontsize=18)

# Add "phi" (for the x-axis)
ax.plot([-2.5, 2.5], [0, 0], color='white', lw=1)
ax.arrow(2.5, 0, 0.2, 0, color='white', head_width=0.1, head_length=0.1)
ax.text(2.6, -0.1, '$\phi$', color='white', fontsize=18)

# Add "phi_+" and "phi_-" for the T < Tc minima
ax.text(minima_phi + 0.1, minima_f + 0.3, '$\phi_+$', color='yellow', fontsize=18)
ax.text(-minima_phi - 0.5, minima_f + 0.3, '$\phi_-$', color='yellow', fontsize=18)
# Add small vertical ticks for the minima
ax.plot([minima_phi, minima_phi], [-0.1, 0.1], color='white', lw=1)
ax.plot([-minima_phi, -minima_phi], [-0.1, 0.1], color='white', lw=1)

# --- Add arrows (annotations) to show dynamics for T < Tc ---
arrow_props = dict(facecolor='white', edgecolor='none', arrowstyle='->', lw=2)
# Arrow from (phi > phi_+) downhill
ax.annotate('', xy=(minima_phi + 0.3, minima_f + 0.1), 
            xytext=(1.8, 1.0), 
            arrowprops=arrow_props)
# Arrow from (phi < -phi_+) downhill
ax.annotate('', xy=(-minima_phi - 0.3, minima_f + 0.1), 
            xytext=(-1.8, 1.0), 
            arrowprops=arrow_props)
# Arrow from (0 < phi < phi_+) downhill
ax.annotate('', xy=(0.3, minima_f), 
            xytext=(0.6, -0.1), 
            arrowprops=arrow_props)
# Arrow from (-phi_+ < phi < 0) downhill
ax.annotate('', xy=(-0.3, minima_f), 
            xytext=(-0.6, -0.1), 
            arrowprops=arrow_props)


# --- Clean up the plot ---
ax.grid(False) # Turn off grid
# 显示坐标轴刻度标签
# ax.set_xticklabels([]) 
# ax.set_yticklabels([]) 

# 显示坐标轴边框
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.spines['bottom'].set_visible(False)
# ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

# 设置坐标轴标签
ax.set_xlabel('$\phi$', color='white', fontsize=14)
ax.set_ylabel('$f$', color='white', fontsize=14)

# Set plot limits
ax.set_xlim(-2.8, 2.8)
ax.set_ylim(min(f_neg) - 0.5, 9.0) # Adjust y-limit to show both curves

plt.title('Ginzburg-Landau Potential: $T > T_c$ vs $T < T_c$', color='yellow', fontsize=18)
plt.savefig('Ginzburg-Landauyon.png')
plt.show()