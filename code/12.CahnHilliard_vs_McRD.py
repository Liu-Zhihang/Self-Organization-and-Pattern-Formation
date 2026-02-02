"""
CahnHilliard vs McRD: Passive vs Active Phase Separation
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.colors import LinearSegmentedColormap
import os

# Parameters
N, L, dx = 128, 50.0, 50.0/128
dt_ch, dt_rd = 0.2, 0.001

# Cahn-Hilliard parameters
M, kappa, phi0 = 1.0, 0.5, 0.0
kx = 2 * np.pi * np.fft.fftfreq(N, d=dx)
ky = 2 * np.pi * np.fft.fftfreq(N, d=dx)
KX, KY = np.meshgrid(kx, ky, indexing='ij')
K2, K4 = KX**2 + KY**2, (KX**2 + KY**2)**2

# McRD parameters
Dm, Dc = 0.05, 5.0
k_on, k_off, k_ant = 0.5, 0.1, 8.0
rho_A, rho_B = 1.5, 1.5

def ch_step(phi, dt):
    nonlinear_hat = np.fft.fft2(phi**3 - phi)
    phi_hat = np.fft.fft2(phi)
    return np.real(np.fft.ifft2((phi_hat - M*dt*K2*nonlinear_hat) / (1 + M*kappa*dt*K4)))

def laplacian(f):
    return (np.roll(f,1,0) + np.roll(f,-1,0) + np.roll(f,1,1) + np.roll(f,-1,1) - 4*f) / dx**2

def rd_step(mA, mB, cA, cB, dt):
    fA = k_on*cA - k_off*mA - k_ant*mA*mB**2
    fB = k_on*cB - k_off*mB - k_ant*mB*mA**2
    mA = np.maximum(mA + dt*(Dm*laplacian(mA) + fA), 0)
    mB = np.maximum(mB + dt*(Dm*laplacian(mB) + fB), 0)
    cA = np.maximum(cA + dt*(Dc*laplacian(cA) - fA), 0)
    cB = np.maximum(cB + dt*(Dc*laplacian(cB) - fB), 0)
    return mA, mB, cA, cB

def run_simulation():
    np.random.seed(42)
    noise1, noise2 = 0.05*(2*np.random.rand(N,N)-1), 0.05*(2*np.random.rand(N,N)-1)
    phi = phi0 + noise1
    mA, mB = 0.5*rho_A + noise1, 0.5*rho_B + noise2
    cA, cB = 0.5*rho_A - noise1, 0.5*rho_B - noise2
    
    frames_ch, frames_rd, times = [phi.copy()], [(mA-mB).copy()], [0.0]
    
    for i in range(100):
        for _ in range(5): phi = ch_step(phi, dt_ch)
        for _ in range(500): mA, mB, cA, cB = rd_step(mA, mB, cA, cB, dt_rd)
        frames_ch.append(phi.copy())
        frames_rd.append((mA-mB).copy())
        times.append((i+1)*5*dt_ch)
        print(f"\rProgress: {i+1}%", end="")
    print("\nDone!")
    return frames_ch, frames_rd, times

def create_animation(frames_ch, frames_rd, times):
    cmap = LinearSegmentedColormap.from_list('ph', ['#2166ac','#f7f7f7','#b2182b'], N=256)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=100)
    vmax_rd = max(abs(np.min(frames_rd)), abs(np.max(frames_rd))) * 0.8
    
    im_ch = axes[0].imshow(frames_ch[0], cmap=cmap, vmin=-1, vmax=1, extent=[0,L,0,L], origin='lower')
    axes[0].set_title('Cahn-Hilliard (Model B)\nPassive Phase Separation', fontweight='bold')
    plt.colorbar(im_ch, ax=axes[0], shrink=0.8, label='$\\phi$')
    
    im_rd = axes[1].imshow(frames_rd[0], cmap=cmap, vmin=-vmax_rd, vmax=vmax_rd, extent=[0,L,0,L], origin='lower')
    axes[1].set_title('McRD (Active Matter)\nEffective Interfacial Tension', fontweight='bold')
    plt.colorbar(im_rd, ax=axes[1], shrink=0.8, label='$m_A - m_B$')
    
    title = fig.suptitle(f't = {times[0]:.1f}', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    def update(i):
        im_ch.set_array(frames_ch[i])
        im_rd.set_array(frames_rd[i])
        title.set_text(f't = {times[i]:.1f}')
        return [im_ch, im_rd, title]
    
    anim = FuncAnimation(fig, update, frames=len(times), interval=100, blit=True)
    output = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CahnHilliard_vs_McRD.mp4')
    anim.save(output, writer=FFMpegWriter(fps=10, bitrate=3000))
    plt.close()
    print(f"Saved: {output}")
    return output

if __name__ == '__main__':
    frames_ch, frames_rd, times = run_simulation()
    create_animation(frames_ch, frames_rd, times)
