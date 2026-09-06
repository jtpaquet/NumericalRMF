import numpy as np
import matplotlib.pyplot as plt
import pickle
from utils import plot_options



if __name__ == "__main__":
    np.set_printoptions(precision=2)
    plot_options["font.size"] = 24

    # Apply plot formatting
    import matplotlib as mpl
    mpl.rcParams.update(plot_options)


# ================================================================
# Run all three gamma values
# ================================================================
T_period = 2 * np.pi
dt = 0.001
steps_per_period = int(T_period / dt)
N_periods = 101
Nr = 32

gamma_c = 1.12 * 11.07
gam_list = [14.9, 16.6, 18.2]
gam_list = [16.6]
colors_gam = {14.9: 'blue', 16.6: 'green', 18.2: 'red'}

plot_times = [10, 40, 50, 60]
plot_times = [10, 30, 45, 60]
print(plot_times)

with open('all_results_Nr64_6283steps_fine.pkl', 'rb') as f:
    all_results = pickle.load(f)


# ================================================================
# ALPHA VS TIME
# ================================================================
print('\n--- Alpha vs Time ---')
fig_alpha, ax_alpha = plt.subplots(figsize=(14, 8))

for gam in gam_list:
    res = all_results[gam]
    ax_alpha.plot(res['times'], res['alphas'], color=colors_gam[gam], lw=2,
                  label=fr'$\gamma$={gam} ($\gamma/\gamma_c$={gam/gamma_c:.2f})')
    print(fr'$\gamma$={gam}: $\alpha$(50T)={res["alphas"][-1]:.4f}')

ax_alpha.axhline(1.0, color='k', ls='--', alpha=0.5)
ax_alpha.scatter([50], [0.48], color='blue', s=200, marker='*', zorder=5, label='Milroy targets')
ax_alpha.scatter([40], [0.98], color='green', s=200, marker='*', zorder=5)
ax_alpha.scatter([27], [0.98], color='red', s=200, marker='*', zorder=5)

ax_alpha.set_xlabel('Time (RMF periods)')
ax_alpha.set_ylabel(r'Penetration Factor $\alpha$')
ax_alpha.set_title(r'$\alpha$ vs Time - $\lambda$=11.07, $N_r$=32')
ax_alpha.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
ax_alpha.grid(True, alpha=0.3)
# ax_alpha.set_xlim(0, 50)
ax_alpha.set_ylim(0, 1.2)
plt.tight_layout()
plt.savefig('alpha_vs_time.pdf')
print('Saved alpha_vs_time.png')

# ================================================================
# FIELD LINE CONTOUR PLOTS — one figure per gamma
# ================================================================
Nr_plot = 100
Ntheta = 100
r_plot = np.linspace(0.02, 1.0, Nr_plot)
theta_plot = np.linspace(0, 2*np.pi, Ntheta)
R_mesh, Theta_mesh = np.meshgrid(r_plot, theta_plot)
X_mesh = R_mesh * np.cos(Theta_mesh)
Y_mesh = R_mesh * np.sin(Theta_mesh)
theta_circle = np.linspace(0, 2*np.pi, 200)

for gam in gam_list:
    res = all_results[gam]
    snapshots = res['snapshots']
    
    fig, axes = plt.subplots(1,4, figsize=(18, 8))
    axes = axes.flatten()
    
    for idx, t in enumerate(plot_times):
        ax = axes[idx]
        A_snap, B_snap, r_snap, tau_snap = snapshots[t]
        
        # Interpolate A to finer grid
        A_interp = np.interp(r_plot, r_snap, np.real(A_snap)) + \
                    1j * np.interp(r_plot, r_snap, np.imag(A_snap))
        
        # Flux function: ψ = Re[A * e^{i(θ-τ)}]
        psi = np.zeros((Ntheta, Nr_plot))
        for i in range(Nr_plot):
            psi[:, i] = np.real(A_interp[i] * np.exp(1j*(theta_plot - tau_snap)))
        
        # alpha from B
        alpha_t = np.abs(B_snap[-1] - B_snap[0])
        
        if t == 0:
            ax.text(0, 0, r'$t$ = 0\nNo field', ha='center', va='center', fontsize=14)
        else:
            psi_max = np.max(np.abs(psi))
            if psi_max > 0.01:
                levels = np.linspace(-psi_max*0.9, psi_max*0.9, 15)
                ax.contour(X_mesh, Y_mesh, psi, levels=30, colors='blue', linewidths=0.8)
        
        ax.plot(np.cos(theta_circle), np.sin(theta_circle), 'k-', lw=2)
        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-1.15, 1.15)
        ax.set_aspect('equal')
        ax.set_title(fr'$t$={t}T, $\alpha$={alpha_t:.3f}', fontsize=24)
        ax.set_xlabel('$x/R$')
        ax.set_ylabel('$y/R$')
    
    fig.suptitle(fr'RMF Penetration - $\gamma$={gam}, $\lambda$=11.07, $N_r$=64', fontsize=36)
    plt.tight_layout()
    fname = f'field_lines_gam{gam}.png'
    plt.savefig(fname, dpi=200)
    print(f'Saved {fname}')
    # plt.show()

plt.close('all')

# ================================================================
# 4. RADIAL PROFILES — one figure per gamma
# ================================================================
for gam in gam_list:
    res = all_results[gam]
    snapshots = res['snapshots']
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 8))
    axes = axes.flatten()
    
    for idx, t in enumerate(plot_times):
        ax = axes[idx]
        A_snap, B_snap, r_snap, tau_snap = snapshots[t]
        
        r_pos = r_snap[1:]  # skip r=0
        dr = r_snap[1] - r_snap[0]
        
        # B_r = |A|/r (transverse radial component amplitude)
        Br_profile = np.abs(A_snap[1:]) / r_pos
        
        # B_theta = |dA/dr| (transverse azimuthal component amplitude)
        dA_dr = np.zeros(len(r_snap), dtype=complex)
        for i in range(1, len(r_snap) - 1):
            dA_dr[i] = (A_snap[i+1] - A_snap[i-1]) / (2*dr)
        dA_dr[0] = (A_snap[1] - A_snap[0]) / dr
        dA_dr[-1] = (A_snap[-1] - A_snap[-2]) / dr
        Btheta_profile = np.abs(dA_dr[1:])
        
        # B_z = B (axial field from Hall current)
        Bz_profile = B_snap[1:]
        
        # Normalize by gamma (= B_omega in normalized units)
        Br_norm = Br_profile / gam
        Btheta_norm = Btheta_profile / gam
        Bz_norm = Bz_profile / abs(Bz_profile[0])
        
        ax.plot(r_pos, Br_norm, 'b-', lw=2, label=r'$|B_r|/B_\omega$')
        ax.plot(r_pos, Btheta_norm, 'g--', lw=2, label=r'$|B_\theta|/B_\omega$')
        ax.plot(r_pos, Bz_norm, 'r:', lw=2, label=r'$B_z/B_{z0}$')
        
        # alpha from B
        alpha_t = np.abs(B_snap[-1] - B_snap[0])
        
        ax.set_xlabel(r'$r/R$')
        # ax.set_ylabel('Normalized Field', fontsize=16)
        ax.set_title(fr'$t$={t}T, $\alpha$={alpha_t:.3f}')
        # ax.legend(loc='best', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        
        if t == 0:
            ax.set_ylim(-0.1, 0.5)
        else:
            ymax = max(np.max(Br_norm), np.max(Btheta_norm),
                       np.max(np.abs(Bz_norm))) * 1.2
            ax.set_ylim(min(-0.1, np.min(Bz_norm)*1.2), max(1.0, ymax))
        ax.set_ylim(-2.5, 2.5)

    for ax in axes:
        ax.tick_params(axis='both', labelsize=20)
    handles, labels = axes[0].get_legend_handles_labels()

    fig.legend(
        handles,
        labels,
        loc='center left',
        bbox_to_anchor=(1.02, 0.5),
        fontsize=24
    )

    fig.supylabel('Normalized Field', fontsize=24)
    fig.suptitle(
        fr'Radial Profiles - $\gamma$={gam} ($\gamma/\gamma_c$={gam/gamma_c:.2f}), $\lambda$=11.07'
    )

    plt.tight_layout()

    fname = f'radial_profiles_gam{gam}.pdf'
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    print(f'Saved {fname}')

plt.close('all')
print('\nDone.')