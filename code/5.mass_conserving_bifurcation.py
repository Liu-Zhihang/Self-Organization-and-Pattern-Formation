import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')
kappa = 0.02
k_fb  = 1.6
def a(m):
    return kappa + k_fb * m**2 / (1.0 + m**2)
def ap(m):
    return k_fb * (2.0*m) / (1.0 + m**2)**2
def f(m, c):
    return a(m) * c - m
def c_nullcline(m):
    return m / a(m)
def g(m, n):
    return a(m) * (n - m) - m
def gm(m, n):
    return ap(m) * (n - m) - a(m) - 1.0
def sn_m_values(m_min=1e-6, m_max=2.5, num=5000):
    m = np.linspace(m_min, m_max, num)
    H = a(m) * (a(m) + 1.0) / ap(m) - m
    idx = np.where(np.sign(H[:-1]) * np.sign(H[1:]) < 0)[0]
    roots = []
    for i in idx:
        lo, hi = m[i], m[i+1]
        f_lo = a(lo) * (a(lo) + 1.0) / ap(lo) - lo
        f_hi = a(hi) * (a(hi) + 1.0) / ap(hi) - hi
        for _ in range(20):
            mid = 0.5 * (lo + hi)
            f_mid = a(mid) * (a(mid) + 1.0) / ap(mid) - mid
            if f_lo * f_mid <= 0:
                hi, f_hi = mid, f_mid
            else:
                lo, f_lo = mid, f_mid
        roots.append(0.5 * (lo + hi))
    return np.array(roots)
sn_ms = sn_m_values()
# Get the two saddle-node m values
sn_m1, sn_m2 = sn_ms[0], sn_ms[1]
sn_ns = sn_ms + (a(sn_ms) + 1.0) / ap(sn_ms)
sn_cs = sn_ns - sn_ms
n_low  = 0.7 * float(np.min(sn_ns))
n_mid  = np.mean(sn_ns) + 0.1
n_high = 1.2 * float(np.max(sn_ns))
n_list = [n_low, n_mid, n_high]
n_labels = ["$n_1$", "$n_2$", "$n_3$"]

def equilibria_on_n(n, m_grid=np.linspace(1e-6, 2.5, 2000)):
    y = g(m_grid, n)
    idx = np.where(np.sign(y[:-1]) * np.sign(y[1:]) <= 0)[0]
    roots = []
    for i in idx:
        lo, hi = m_grid[i], m_grid[i+1]
        f_lo, f_hi = y[i], y[i+1]
        m_star = 0.5 * (lo + hi)
        for _ in range(20):
            mid = 0.5 * (lo + hi)
            f_mid = g(mid, n)
            if f_lo * f_mid <= 0:
                hi, f_hi = mid, f_mid
            else:
                lo, f_lo = mid, f_mid
        m_star = 0.5 * (lo + hi)
        c_star = n - m_star
        stab = gm(m_star, n) < 0.0
        roots.append((m_star, c_star, stab))
    return roots
def plot_flow_arrows(n, eqs, ax, m_max_plot):
    m_stars_sorted = sorted([eq[0] for eq in eqs])
    arrow_props = dict(arrowstyle="->", color="red", lw=1.5, ls='-')
    if len(m_stars_sorted) == 1:
        m_stable = m_stars_sorted[0]
        m_start_1 = m_stable - 0.4
        m_end_1   = m_stable - 0.1
        if m_start_1 > 0.05:
            ax.annotate("", xy=(m_end_1, n - m_end_1), xytext=(m_start_1, n - m_start_1),
                         arrowprops=arrow_props)
        m_start_2 = m_stable + 0.4
        m_end_2   = m_stable + 0.1
        if m_start_2 < m_max_plot:
            ax.annotate("", xy=(m_end_2, n - m_end_2), xytext=(m_start_2, n - m_start_2),
                         arrowprops=arrow_props)
    elif len(m_stars_sorted) == 3:
        m_low_stable, m_unstable, m_high_stable = m_stars_sorted
        m_start_1 = m_unstable - 0.05
        m_end_1   = m_low_stable + 0.05
        ax.annotate("", xy=(m_end_1, n - m_end_1), xytext=(m_start_1, n - m_start_1),
                     arrowprops=arrow_props)
        m_start_2 = m_unstable + 0.05
        m_end_2   = m_high_stable - 0.05
        ax.annotate("", xy=(m_end_2, n - m_end_2), xytext=(m_start_2, n - m_start_2),
                     arrowprops=arrow_props)
        m_start_3 = m_low_stable - 0.1
        if m_start_3 > 0.05:
            ax.annotate("", xy=(m_low_stable-0.02, n - (m_low_stable-0.02)), xytext=(m_start_3, n - m_start_3),
                         arrowprops=arrow_props)
        m_start_4 = m_high_stable + 0.2
        if m_start_4 < m_max_plot:
            ax.annotate("", xy=(m_high_stable+0.02, n - (m_high_stable+0.02)), xytext=(m_start_4, n - m_start_4),
                         arrowprops=arrow_props)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
ax1.set_facecolor('black')
ax2.set_facecolor('black')
fig.patch.set_facecolor('black')

# Phase space with segmented f=0 nullcline
m_max = 2.0
m = np.linspace(1e-6, m_max, 3000)
c_nc = c_nullcline(m)

# Draw f=0 as solid for m<=sn_m1 and m>=sn_m2, dashed for sn_m1<m<sn_m2
mask_left  = m <= sn_m1
mask_mid   = (m > sn_m1) & (m < sn_m2)
mask_right = m >= sn_m2

ax1.plot(m[mask_left],  c_nc[mask_left],  linewidth=2, color='white', label="$f(m,c)=0$")
ax1.plot(m[mask_mid],   c_nc[mask_mid],   linestyle="--", linewidth=2, color='white')
ax1.plot(m[mask_right], c_nc[mask_right], linewidth=2, color='white')

# Plot n lines and equilibria
for i, n in enumerate(n_list):
    m_line = np.linspace(0.0, min(m_max, n), 600)
    c_line = n - m_line
    ax1.plot(m_line, c_line, linewidth=1.2, color='gray')
    # equilibria markers
    eqs = equilibria_on_n(n, np.linspace(1e-6, m_max, 2000))
    for (m_star, c_star, stable) in eqs:
        if stable:
            ax1.plot(m_star, c_star, marker="o", markersize=8, color='white', zorder=5)
        else:
            ax1.plot(m_star, c_star, marker="o", markersize=9, markerfacecolor="none",
                     markeredgewidth=2.0, color='white', zorder=6)
    # arrows along flow (skip vicinity of roots)
    pts = np.linspace(0.12, 0.88, 6) * (min(m_max, n))
    for mi in pts:
        ci = n - mi
        if abs(g(mi, n)) < 2e-3:  # avoid near equilibria
            continue
        s = 0.06 * min(m_max, n)  # step
        dx, dy = (s, -s) if g(mi, n) > 0 else (-s, s)
        ax1.annotate("", xy=(mi+dx, ci+dy), xytext=(mi, ci),
                     arrowprops=dict(arrowstyle="->", color="red", lw=1.2))

ax1.text(0.1, n_list[0] - 0.1, "$n_1$", fontsize=12, color='gray')
ax1.text(0.1, n_list[1] - 0.1, "$n_2$", fontsize=12, color='gray')
ax1.text(0.1, n_list[2] - 0.1, "$n_3$", fontsize=12, color='gray')
ax1.text(1.2, 2.0, "$f(m,c) = 0$", fontsize=12, color='white')
ax1.text(0.6, 1.6, "reactive flow", fontsize=12, rotation=-40, color='red')
ax1.set_xlim(0, m_max)
ax1.set_ylim(0, float(np.nanmax(c_nc))*1.05)
ax1.set_xlabel("m", fontsize=14, color='white')
ax1.set_ylabel("c", fontsize=14, color='white')
ax1.set_title("(m,c)-phase space (only unstable segment dashed)", fontsize=14, color='white')
ax1.grid(True, linestyle="--", alpha=0.3, color='gray')
ax1.legend(loc="best")
ax1.tick_params(colors='white')

# Bifurcation diagram
m_branch = np.linspace(1e-6, m_max, 2000)
c_branch = c_nullcline(m_branch)
n_branch = m_branch + c_branch
stab_branch = gm(m_branch, n_branch) < 0.0
mask = stab_branch.astype(int)
change_idx = np.where(mask[1:] != mask[:-1])[0] + 1
segments = np.split(np.arange(len(m_branch)), change_idx)

for seg in segments:
    if stab_branch[seg[0]]:
        style = "-"
    else:
        style = (0, (5, 3))
        
    ax2.plot(n_branch[seg], c_branch[seg], linestyle=style, linewidth=2, color='white')
ylim_top = 2.5
ax2.set_ylim(-0.1, ylim_top) 

for i, n in enumerate(n_list):
    ax2.text(n, -0.05, n_labels[i], ha='center', va='top', fontsize=12, color='white')
    
    eqs = equilibria_on_n(n, np.linspace(1e-6, m_max, 2000))
    eqs_sorted = sorted(eqs, key=lambda x: x[1])
    if len(eqs_sorted) == 1:
        c_star = eqs_sorted[0][1]
        ax2.annotate("", xy=(n, c_star - 0.1), xytext=(n, -0.1),
                    arrowprops=dict(arrowstyle="->", color="red", lw=1.5))
        ax2.plot(n, c_star, marker="o", markersize=9, color='white')

    elif len(eqs_sorted) == 3:
        c_low_stable = eqs_sorted[0][1]
        c_unstable = eqs_sorted[1][1]
        c_high_stable = eqs_sorted[2][1]

        ax2.annotate("", xy=(n, c_low_stable - 0.1), xytext=(n, -0.1),
                    arrowprops=dict(arrowstyle="->", color="red", lw=1.5))

        ax2.annotate("", xy=(n, c_low_stable + 0.1), xytext=(n, c_unstable - 0.1),
                    arrowprops=dict(arrowstyle="->", color="red", lw=1.5))

        if n != n_list[1]: 
            ax2.annotate("", xy=(n, ylim_top - 0.1), xytext=(n, c_high_stable + 0.1),
                        arrowprops=dict(arrowstyle="->", color="red", lw=1.5))
            ax2.plot(n, c_high_stable, marker="o", markersize=9, color='white')
        ax2.plot(n, c_unstable, marker="o", markersize=9, markerfacecolor="black", markeredgecolor='white', mew=1.5)
        ax2.plot(n, c_low_stable, marker="o", markersize=9, color='white')

if n_list[0] in [n_low, n_mid, n_high]:
    eqs_n1 = equilibria_on_n(n_list[0], np.linspace(1e-6, m_max, 2000))
    eqs_n1_sorted = sorted(eqs_n1, key=lambda x: x[1])
    if len(eqs_n1_sorted) == 1:
        c_star_n1 = eqs_n1_sorted[0][1]
        ax2.annotate("", xy=(n_list[0], c_star_n1 + 0.1), xytext=(n_list[0], ylim_top - 0.1),
                    arrowprops=dict(arrowstyle="->", color="red", lw=1.5))
if n_list[2] in [n_low, n_mid, n_high]:
    eqs_n3 = equilibria_on_n(n_list[2], np.linspace(1e-6, m_max, 2000))
    eqs_n3_sorted = sorted(eqs_n3, key=lambda x: x[1])
    if len(eqs_n3_sorted) == 1:
        c_star_n3 = eqs_n3_sorted[0][1]
        ax2.annotate("", xy=(n_list[2], c_star_n3 + 0.1), xytext=(n_list[2], ylim_top - 0.1),
                    arrowprops=dict(arrowstyle="->", color="red", lw=1.5))
ax2.text(0.2, 1.0, "$c^*(n)$", fontsize=14, rotation=50, color='white')
ax2.set_xlabel("n", fontsize=14, color='white')
ax2.set_ylabel("$c^*$", fontsize=14, color='white')
ax2.set_title("Bifurcation diagram", fontsize=14, color='white')
ax2.grid(True, linestyle="--", alpha=0.3, color='gray')
ax2.tick_params(colors='white')
plt.tight_layout()
plt.savefig("combined_diagram.png", dpi=200, bbox_inches='tight')
plt.show()