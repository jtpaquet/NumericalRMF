import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import paths
import pickle
import numpy as np
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import subprocess

def make_mp4():
    # ------------------------------------------------
    # Load data
    # ------------------------------------------------
    with open(paths.results("all_results_Nr64_6283steps_fine.pkl"), "rb") as f:
        all_results = pickle.load(f)

    gam = 16.6
    lambda_ = 11.07
    snapshots = all_results[gam]['snapshots']

    print(f"Available times: {list(snapshots.keys())[:10]} ...\n")

    # ------------------------------------------------
    # Select times
    # ------------------------------------------------
    dt = 0.5
    times = sorted(
        t for t in snapshots
        if 0.0 <= t <= 100.0
        and abs((t / dt) - round(t / dt)) < 1e-6
    )
    print(f"Making MP4 from {len(times)} snapshots")

    # ------------------------------------------------
    # Grids for contour plot
    # ------------------------------------------------
    Nr_plot = 100
    Ntheta = 100

    r_plot = np.linspace(0.02, 1.0, Nr_plot)
    theta_plot = np.linspace(0, 2*np.pi, Ntheta)

    R_mesh, Theta_mesh = np.meshgrid(r_plot, theta_plot)
    X_mesh = R_mesh * np.cos(Theta_mesh)
    Y_mesh = R_mesh * np.sin(Theta_mesh)

    theta_circle = np.linspace(0, 2*np.pi, 200)

    # ------------------------------------------------
    # Setup ffmpeg pipe
    # ------------------------------------------------
    width, height = 900, 900   # 2x2 layout
    fps = 30
    output_file = paths.animation("field_and_profiles_gam16p6.mp4")

    cmd = [
        "ffmpeg", "-y",
        "-f", "rawvideo",
        "-vcodec", "rawvideo",
        "-pix_fmt", "rgb24",
        "-s", f"{width}x{height}",
        "-r", str(fps),
        "-i", "-",
        "-an",
        "-c:v", "libx264",
        "-crf", "23",
        "-preset", "medium",
        "-pix_fmt", "yuv420p",
        output_file
    ]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)

    # ------------------------------------------------
    # Loop over frames
    # ------------------------------------------------
    for i, t in enumerate(times):
        A_snap, B_snap, r_snap, tau_snap = snapshots[t]

        fig, axes = plt.subplots(2, 2, figsize=(9, 9), dpi=100)
        ax_c, ax_br, ax_bt, ax_bz = axes.flatten()

        # =================================================
        # 1. CONTOUR PLOT (ψ)
        # =================================================
        A_interp = (
            np.interp(r_plot, r_snap, np.real(A_snap)) +
            1j * np.interp(r_plot, r_snap, np.imag(A_snap))
        )

        psi = np.zeros((Ntheta, Nr_plot))
        for j in range(Nr_plot):
            psi[:, j] = np.real(A_interp[j] * np.exp(1j * (theta_plot - tau_snap)))

        psi_max = np.max(np.abs(psi))
        if psi_max > 0.01:
            levels = np.linspace(-0.9 * psi_max, 0.9 * psi_max, 20)
            ax_c.contour(
                X_mesh, Y_mesh, psi,
                levels=levels,
                colors="black",
                linewidths=0.7,
                linestyles="solid"
            )

        ax_c.plot(np.cos(theta_circle), np.sin(theta_circle), "k-", lw=1.5)
        ax_c.set_aspect("equal")
        ax_c.set_xlim(-1.15, 1.15)
        ax_c.set_ylim(-1.15, 1.15)
        ax_c.set_title("Field lines", fontsize=11)
        ax_c.set_xlabel("x/R")
        ax_c.set_ylabel("y/R")

        # =================================================
        # RADIAL PROFILES
        # =================================================
        r_pos = r_snap[1:]
        dr = r_snap[1] - r_snap[0]

        # B_r = |A| / r
        Br = np.abs(A_snap[1:]) / r_pos

        # B_theta = |dA/dr|
        dA_dr = np.zeros_like(A_snap, dtype=complex)
        dA_dr[1:-1] = (A_snap[2:] - A_snap[:-2]) / (2 * dr)
        dA_dr[0] = (A_snap[1] - A_snap[0]) / dr
        dA_dr[-1] = (A_snap[-1] - A_snap[-2]) / dr
        Btheta = np.abs(dA_dr[1:])

        # B_z
        Bz = B_snap[1:]

        # Normalize by gamma
        Br /= gam
        Btheta /= gam
        # Bz /= lambda_

        # =================================================
        # 2. |B_r|
        # =================================================
        ax_br.plot(r_pos, Br, "b-", lw=2)
        ax_br.set_title(r"$|B_r|/B_\omega$", fontsize=11)
        ax_br.set_xlim(0, 1)
        ax_br.set_ylim(-0.1, 1.1)
        ax_br.grid(alpha=0.3)

        # =================================================
        # 3. |B_theta|
        # =================================================
        ax_bt.plot(r_pos, Btheta, "g-", lw=2)
        ax_bt.set_title(r"$|B_\theta|/B_\omega$", fontsize=11)
        ax_bt.set_xlim(0, 1)
        ax_bt.set_ylim(-0.1, 2.1)
        ax_bt.grid(alpha=0.3)

        # =================================================
        # 4. B_z
        # =================================================
        ax_bz.plot(r_pos, Bz, "r-", lw=2)
        ax_bz.set_title(r"$B_z/B_\omega$", fontsize=11)
        ax_bz.set_xlim(0, 1)
        ax_bz.set_ylim(-1, 1)
        ax_bz.grid(alpha=0.3)

        # =================================================
        # Titles
        # =================================================
        alpha_t = np.abs(B_snap[-1] - B_snap[0])
        fig.suptitle(
            f"γ = {gam},  t = {t:.2f},  α = {alpha_t:.3f}",
            fontsize=14
        )

        plt.tight_layout(rect=[0, 0, 1, 0.95])

        # ---- Send frame to ffmpeg ----
        fig.canvas.draw()
        frame = np.asarray(fig.canvas.renderer.buffer_rgba())[:, :, :3]
        proc.stdin.write(frame.tobytes())
        plt.close(fig)

        if i % 10 == 0:
            print(f"Frame {i}/{len(times)}")

    # ------------------------------------------------
    # Finish
    # ------------------------------------------------
    proc.stdin.close()
    proc.wait()
    print(f"\nSaved MP4 → {output_file}")


if __name__ == "__main__":
    make_mp4()
