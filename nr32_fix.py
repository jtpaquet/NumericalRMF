import numpy as np
import matplotlib.pyplot as plt
import pickle

class HugrassRMF:
    def __init__(self, Nr=32, lam=11.07, gam=16.6):
        self.Nr = Nr
        self.lam = lam
        self.gam = gam
        self.lam2 = lam**2
        self.lam4 = lam**4
        
        self.dr = 1.0 / (Nr - 1)
        self.r = np.linspace(0, 1, Nr)
        
        self.A = np.zeros(Nr, dtype=complex)
        self.B = np.zeros(Nr)
        
        self.tau = 0.0
        self.T_period = 2 * np.pi
        self.rise_time = 3.0 * self.T_period
        
    def gamma_t(self):
        return self.gam * (1 - np.exp(-self.tau / self.rise_time))
    
    def A_external(self):
        return self.gamma_t() * np.exp(-1j * self.tau)
    
    def compute_rhs(self, A, B):
        r, dr, Nr = self.r, self.dr, self.Nr
        rhs_A = np.zeros(Nr, dtype=complex)
        rhs_B = np.zeros(Nr)
        
        # Compute lapA for all interior points
        lapA = np.zeros(Nr, dtype=complex)
        for i in range(1, Nr - 1):
            d2A = (A[i+1] - 2*A[i] + A[i-1]) / dr**2
            dA = (A[i+1] - A[i-1]) / (2*dr)
            lapA[i] = d2A + dA/r[i] - A[i]/(r[i]**2)
        
        # Compute product P = Im(lapA * conj(A)) for Hall term in B
        product = np.imag(lapA * np.conj(A))
        
        for i in range(1, Nr - 1):
            ri = r[i]
            
            # --- A equation ---
            diff_A = lapA[i] / (2*self.lam2)
            dB = (B[i+1] - B[i-1]) / (2*dr)
            hall_A = -1j * A[i] * dB / (2*ri)
            rhs_A[i] = diff_A + hall_A
            
            # --- B equation ---
            d2B = (B[i+1] - 2*B[i] + B[i-1]) / dr**2
            dB_dr = (B[i+1] - B[i-1]) / (2*dr)
            diff_B = (d2B + dB_dr/ri) / (2*self.lam2)
            
            # Hall term: (1/(lam4 * r)) * dP/dr
            # product = Im(L1 * conj(A1)) = P
            d_product = (product[i+1] - product[i-1]) / (2*dr)
            hall_B = d_product / (self.lam4 * ri)
            hall_B = 0.5 * d_product / (2 * self.lam4 * ri)
            rhs_B[i] = diff_B + hall_B
        
        # Axis BCs for RHS
        rhs_A[0] = 0
        rhs_B[0] = 2*(B[1] - B[0]) / dr**2 / (2*self.lam2)
        return rhs_A, rhs_B
    
    def apply_bc(self, A, B):
        A[0] = 0
        A_ext = self.A_external()
        # Robin BC: A + dA/dr = 2*gamma*exp(-i*tau) at r=1
        A[-1] = (2*A_ext + A[-2]/self.dr) / (1 + 1/self.dr)
        # Dirichlet for B at boundary
        B[-1] = 0
        return A, B
    
    def step(self, dt):
        # Predictor
        rhs_A, rhs_B = self.compute_rhs(self.A, self.B)
        A_star = self.A.copy()
        B_star = self.B.copy()
        A_star[1:-1] += dt * rhs_A[1:-1]
        B_star[0] += dt * rhs_B[0]
        B_star[1:-1] += dt * rhs_B[1:-1]
        
        self.tau += dt
        A_star, B_star = self.apply_bc(A_star, B_star)
        
        # Corrector
        rhs_A2, rhs_B2 = self.compute_rhs(A_star, B_star)
        self.A[1:-1] += 0.5*dt * (rhs_A[1:-1] + rhs_A2[1:-1])
        self.B[0] += 0.5*dt * (rhs_B[0] + rhs_B2[0])
        self.B[1:-1] += 0.5*dt * (rhs_B[1:-1] + rhs_B2[1:-1])
        self.A, self.B = self.apply_bc(self.A, self.B)
    
    def compute_alpha(self):
        # alpha = |B(r=1) - B(r=0)| = |0 - B[0]| = |B[0]|
        # B[0] goes negative for penetration, B[-1] = 0
        return np.abs(self.B[-1] - self.B[0])


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
colors_gam = {14.9: 'blue', 16.6: 'green', 18.2: 'red'}

plot_times = [10, 40, 50, 60]
print(plot_times)

all_results = {}

for gam in gam_list:
    print(f'\nRunning γ={gam} (γ/γc={gam/gamma_c:.2f})...')
    sim = HugrassRMF(Nr=Nr, lam=11.07, gam=gam)
    
    t_list, a_list = [], []
    snapshots = {}
    snapshots[0] = (sim.A.copy(), sim.B.copy(), sim.r.copy(), sim.tau)
    
    for period in range(N_periods):
        for _ in range(steps_per_period):
            sim.step(dt)
            
        
        alpha = sim.compute_alpha()
        t_list.append(period + 1)
        a_list.append(alpha)
        
        if (period + 1) in plot_times:
            snapshots[period + 1] = (sim.A.copy(), sim.B.copy(), sim.r.copy(), sim.tau)
        
        if (period + 1) % 10 == 0:
            print(f'  t={period+1}T: α={alpha:.4f}, B[0]={sim.B[0]:.4f}')
    
    all_results[gam] = {
        'times': np.array(t_list),
        'alphas': np.array(a_list),
        'snapshots': snapshots
    }


print(all_results)
with open(f"all_results_Nr{Nr}_{steps_per_period}steps.pkl", "wb") as f:
    pickle.dump(all_results, f)


# ================================================================
# ALPHA VS TIME
# ================================================================
print('\n--- Alpha vs Time ---')
fig_alpha, ax_alpha = plt.subplots(figsize=(10, 7))

for gam in gam_list:
    res = all_results[gam]
    ax_alpha.plot(res['times'], res['alphas'], color=colors_gam[gam], lw=2,
                  label=f'γ={gam} (γ/γc={gam/gamma_c:.2f})')
    print(f'γ={gam}: α(50T)={res["alphas"][-1]:.4f}')

ax_alpha.axhline(1.0, color='k', ls='--', alpha=0.5)
ax_alpha.scatter([50], [0.48], color='blue', s=100, marker='*', zorder=5, label='Milroy targets')
ax_alpha.scatter([40], [0.98], color='green', s=100, marker='*', zorder=5)
ax_alpha.scatter([27], [0.98], color='red', s=100, marker='*', zorder=5)

ax_alpha.set_xlabel('Time (RMF periods)', fontsize=13)
ax_alpha.set_ylabel('Penetration Factor α', fontsize=13)
ax_alpha.set_title('α vs Time — λ=11.07, Nr=32', fontsize=14)
ax_alpha.legend(loc='upper left', fontsize=11)
ax_alpha.grid(True, alpha=0.3)
# ax_alpha.set_xlim(0, 50)
ax_alpha.set_ylim(0, 1.2)
plt.tight_layout()
plt.savefig('alpha_vs_time.png', dpi=150)
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
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
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
            ax.text(0, 0, 't = 0\nNo field', ha='center', va='center', fontsize=14)
        else:
            psi_max = np.max(np.abs(psi))
            if psi_max > 0.01:
                levels = np.linspace(-psi_max*0.9, psi_max*0.9, 15)
                ax.contour(X_mesh, Y_mesh, psi, levels=30, colors='blue', linewidths=0.8)
        
        ax.plot(np.cos(theta_circle), np.sin(theta_circle), 'k-', lw=2)
        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-1.15, 1.15)
        ax.set_aspect('equal')
        ax.set_title(f't={t}T, α={alpha_t:.3f}', fontsize=12)
        ax.set_xlabel('x/R')
        ax.set_ylabel('y/R')
    
    fig.suptitle(f'RMF Field Lines — γ={gam}, λ=11.07, Nr=32', fontsize=14)
    plt.tight_layout()
    fname = f'field_lines_gam{gam}.png'
    plt.savefig(fname, dpi=150)
    print(f'Saved {fname}')

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
        Bz_norm = Bz_profile / gam
        
        ax.plot(r_pos, Br_norm, 'b-', lw=2, label=r'$|B_r|/B_\omega$')
        ax.plot(r_pos, Btheta_norm, 'g--', lw=2, label=r'$|B_\theta|/B_\omega$')
        ax.plot(r_pos, Bz_norm, 'r:', lw=2, label=r'$B_z/B_\omega$')
        
        # alpha from B
        alpha_t = np.abs(B_snap[-1] - B_snap[0])
        
        ax.set_xlabel('r/R', fontsize=11)
        ax.set_ylabel('Normalized Field', fontsize=11)
        ax.set_title(f't={t}T, α={alpha_t:.3f}', fontsize=12)
        ax.legend(loc='best', fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        
        if t == 0:
            ax.set_ylim(-0.1, 0.5)
        else:
            ymax = max(np.max(Br_norm), np.max(Btheta_norm),
                       np.max(np.abs(Bz_norm))) * 1.2
            ax.set_ylim(min(-0.1, np.min(Bz_norm)*1.2), max(1.0, ymax))
    
    fig.suptitle(f'Radial Profiles — γ={gam} (γ/γc={gam/gamma_c:.2f}), λ=11.07', fontsize=14)
    plt.tight_layout()
    fname = f'radial_profiles_gam{gam}.png'
    plt.savefig(fname, dpi=150)
    print(f'Saved {fname}')

plt.close('all')
print('\nDone.')