import numpy as np
import matplotlib.pyplot as plt

def plot_transcritical_bifurcation(lam, delta):
    """
    Plots the bifurcation diagram for the transcritical bifurcation
    in the infection model.
    
    Args:
        lam (float): Infection rate lambda.
        delta (float): Recovery rate delta.
    """
    n_c = delta / lam
    n_vals = np.linspace(0, 2 * n_c, 400)
    
    # Fixed point 1: m1* = 0
    m1_star = np.zeros_like(n_vals)
    
    # Fixed point 2: m2* = n - delta/lambda
    m2_star = n_vals - n_c
    
    # Stability analysis
    # m1* is stable for n < n_c, unstable for n > n_c
    n_stable1 = n_vals[n_vals < n_c]
    m_stable1 = np.zeros_like(n_stable1)
    
    n_unstable1 = n_vals[n_vals >= n_c]
    m_unstable1 = np.zeros_like(n_unstable1)
    
    # m2* is stable for n > n_c (and physically relevant)
    n_stable2 = n_vals[n_vals >= n_c]
    m_stable2 = n_stable2 - n_c
    
    plt.figure(figsize=(8, 6))
    
    # Plot stable branches
    plt.plot(n_stable1, m_stable1, 'b-', lw=2.5, label='Stable Fixed Points')
    plt.plot(n_stable2, m_stable2, 'b-', lw=2.5)
    
    # Plot unstable branches
    plt.plot(n_unstable1, m_unstable1, 'r--', dashes=(5, 5), lw=2, label='Unstable Fixed Points')
    # The other branch of m2* is unphysical (m<0) and unstable, so we don't plot it.
    
    plt.axvline(n_c, color='k', linestyle=':', label=f'Bifurcation Point $n_c = \\delta/\\lambda = {n_c:.2f}$')
    
    plt.title('Transcritical Bifurcation in the Infection Model')
    plt.xlabel('Total Population Density $n$')
    plt.ylabel('Infected Population $m^*$')
    plt.legend()
    plt.grid(True)
    plt.ylim(-0.1 * n_c, 1.1 * n_c)
    plt.xlim(0, 2 * n_c)
    plt.show()

if __name__ == '__main__':
    # Define model parameters
    lambda_rate = 0.1  # Infection rate
    delta_rate = 0.5  # Recovery rate
    
    print("Plotting the transcritical bifurcation diagram...")
    plot_transcritical_bifurcation(lam=lambda_rate, delta=delta_rate)