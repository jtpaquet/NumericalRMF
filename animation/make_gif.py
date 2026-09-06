import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Use the non-interactive backend
import matplotlib.pyplot as plt
import subprocess

def make_mp4():
    # ------------------------------------------------
    # Load data
    # ------------------------------------------------
    with open("all_results_Nr64_6283steps_fine.pkl", "rb") as f:
        all_results = pickle.load(f)
    
    snapshots = all_results[16.6]['snapshots']
    print(f"Available times: {list(snapshots.keys())[:10]} ...")

    # ------------------------------------------------
    # Select times
    # ------------------------------------------------
    times = sorted(t for t in snapshots if 0.0 <= t <= 60.0)
    print(f"Making MP4 from {len(times)} snapshots")

    # ------------------------------------------------
    # Plot grid
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
    width, height = 500, 500  # matches figsize=(5,5) and dpi=100
    fps = 30
    output_file = "field_lines_gam16p6.mp4"

    cmd = [
        "ffmpeg",
        "-y",
        "-f", "rawvideo",
        "-vcodec", "rawvideo",
        "-pix_fmt", "rgb24",
        "-s", f"{width}x{height}",
        "-r", str(fps),
        "-i", "-",  # read from stdin
        "-an",
        "-c:v", "libx264",
        "-crf", "23",
        "-preset", "medium",
        "-pix_fmt", "yuv420p",
        output_file
    ]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)

    # ------------------------------------------------
    # Generate frames and pipe to ffmpeg
    # ------------------------------------------------
    for i, t in enumerate(times):
        A_snap, B_snap, r_snap, tau_snap = snapshots[t]

        fig, ax = plt.subplots(figsize=(5, 5), dpi=100)

        # ---- Interpolate A ----
        A_interp = (
            np.interp(r_plot, r_snap, np.real(A_snap)) +
            1j * np.interp(r_plot, r_snap, np.imag(A_snap))
        )

        # ---- Flux function ψ ----
        psi = np.zeros((Ntheta, Nr_plot))
        for j in range(Nr_plot):
            psi[:, j] = np.real(A_interp[j] * np.exp(1j * (theta_plot - tau_snap)))

        # ---- alpha from B ----
        alpha_t = np.abs(B_snap[-1] - B_snap[0])

        # ---- Plot ----
        psi_max = np.max(np.abs(psi))
        if psi_max > 0.01:
            levels = np.linspace(-0.9 * psi_max, 0.9 * psi_max, 20)
            ax.contour(X_mesh, Y_mesh, psi, levels=levels, colors="black", linewidths=0.8)

        ax.plot(np.cos(theta_circle), np.sin(theta_circle), "k-", lw=2)
        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-1.15, 1.15)
        ax.set_aspect("equal")
        ax.set_title(f"γ = 16.6,  t = {t:.2f},  α = {alpha_t:.3f}", fontsize=12)
        ax.set_xlabel("x/R")
        ax.set_ylabel("y/R")
        plt.tight_layout()

        # ---- Convert figure to RGB array ----
        fig.canvas.draw()
        frame = np.asarray(fig.canvas.renderer.buffer_rgba())
        frame = frame[:, :, :3]  # drop alpha channel
        proc.stdin.write(frame.tobytes())
        plt.close(fig)

        if i % 20 == 0:
            print(f"Frame {i}/{len(times)}")

    # ------------------------------------------------
    # Finish writing video
    # ------------------------------------------------
    proc.stdin.close()
    proc.wait()
    print(f"Saved MP4 → {output_file}")


if __name__ == "__main__":
    make_mp4()
