import numpy as np
import matplotlib.pyplot as plt

# Set dark background style to match lecture board
plt.style.use('dark_background')

def plot_rho_nullclines(mu, n, mt_max=6):
    """
    Plots the nullclines for the 2D RhoGTPase model.
    
    Args:
        mu (float): Control parameter (relative inactivation rate).
        n (float): Control parameter (total protein mass).
        mt_max (float): Maximum M_T value for plotting range.
    """
    
    # Create an array for M_T values
    M_T = np.linspace(0, mt_max, 500)
    
    # 1. Calculate M_D for the M_T-nullcline (dM_T/dt = 0)
    # Equation: -mu*M_T + (1+M_T^2)*M_D = 0  =>  M_D = mu*M_T / (1+M_T^2)
    M_D_nullcline_T = (mu * M_T) / (1 + M_T**2)
    
    # 2. Calculate M_D for the M_D-nullcline (dM_D/dt = 0)
    # Equation: (n - M_T - M_D) - (1+M_T^2)*M_D = 0 
    # Derived: M_D = (n - M_T) / (2 + M_T^2)
    M_D_nullcline_D = (n - M_T) / (2 + M_T**2)
    
    # Filter out non-physical negative values (concentration cannot be negative)
    M_D_nullcline_D[M_D_nullcline_D < 0] = np.nan
    
    # Plotting
    plt.figure(figsize=(10, 7))
    
    # Plot M_T-nullcline (Red solid line)
    plt.plot(M_T, M_D_nullcline_T, 'r-', linewidth=2.5, label=r'$\dot{M}_T = 0$ nullcline')
    
    # Plot M_D-nullcline (Cyan dashed line)
    plt.plot(M_T, M_D_nullcline_D, 'c--', linewidth=2.5, label=r'$\dot{M}_D = 0$ nullcline')
    
    # Labeling
    plt.xlabel(r'Active (Membrane) $M_T$', color='white', fontsize=12)
    plt.ylabel(r'Inactive (Membrane) $M_D$', color='white', fontsize=12)
    plt.title(f'RhoGTPase Nullclines ($\mu={mu}, n={n}$)', color='white', fontsize=14)
    
    # Styling
    legend = plt.legend(fontsize=12)
    plt.setp(legend.get_texts(), color='white')
    plt.grid(True, linestyle=':', alpha=0.4, color='gray')
    
    # Set axes limits to focus on the relevant interaction area
    plt.ylim(0, n / 2 + 0.5) 
    plt.xlim(0, mt_max)
    
    # Ensure spines are visible in dark mode
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_color('white')
    ax.tick_params(colors='white')

    plt.savefig('rho_nullclines.png', bbox_inches='tight')

# --- Example Usage ---
# Parameters chosen to replicate the crossing shown in the lecture board.
# mu=4.0 creates a pronounced "bell" shape.
# n=10.0 provides enough mass for the intersection to occur on the falling slope.
plot_rho_nullclines(mu=4.0, n=10.0, mt_max=6.0)