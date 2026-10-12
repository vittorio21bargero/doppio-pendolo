"""Make the figures in figures/ from the text files written by ./double_pendulum in data/.

Run from the repository root:  python scripts/plot_results.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

METHODS = ["Euler", "RK2", "RK4", "Verlet"]
COLORS = {"Euler": "#d62728", "RK2": "#ff7f0e", "RK4": "#1f77b4", "Verlet": "#2ca02c"}
H_COLORS = {2: "#d62728", 3: "#ff7f0e", 4: "#e6b800", 5: "#2ca02c", 6: "#1f77b4"}  # h = 10^-i

plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})


def load_grid(path):
    """Grid files: first row = (N, theta2 values), then rows (theta1, values...)."""
    raw = np.loadtxt(path)
    th2 = raw[0, 1:]
    th1 = raw[1:, 0]
    return th1, th2, raw[1:, 1:]


def energy():
    """Cumulative relative energy error xi_cum(t) for each method and time step h."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for ax, m in zip(axes.flat, METHODS):
        for i, c in H_COLORS.items():
            t, *_, err_cum = np.loadtxt(DATA / "energy" / f"{m}_dt_1e-{i}.txt", unpack=True)
            ax.semilogy(t[1:], err_cum[1:], ".", ms=3, color=c, label=rf"$h=10^{{-{i}}}$")
        ax.set_title(m if m != "Verlet" else "Position Verlet")
        ax.set_ylabel(r"$\xi_{cum}(t)$")
    for ax in axes[1]:
        ax.set_xlabel("t [s]")
    axes[0, 0].legend(frameon=False, fontsize=9, markerscale=3)
    fig.suptitle(r"Cumulative energy error, $\theta_1=\theta_2=\pi/2$ at rest")
    fig.tight_layout()
    fig.savefig(FIG / "energy_error.png", dpi=150)
    plt.close(fig)


def convergence():
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for m in METHODS:
        h, _, _, err = np.loadtxt(DATA / "convergence" / f"{m}.txt", unpack=True)
        ok = err > 1e-14  # below this, round-off dominates
        ax.loglog(h[ok], err[ok], "o-", ms=4, color=COLORS[m], label=m)

    # reference slopes h^p, drawn just below the method of that order
    for m, p, h_lo in [("Euler", 1, 1e-6), ("RK2", 2, 1e-6), ("RK4", 4, 3e-4)]:
        h, _, _, err = np.loadtxt(DATA / "convergence" / f"{m}.txt", unpack=True)
        h0 = 2.0 ** -6
        y0 = 0.15 * err[np.argmin(abs(h - h0))]
        h_ref = np.array([h_lo, h0])
        ax.loglog(h_ref, y0 * (h_ref / h0) ** p, "k--", lw=0.8)
        ax.text(h0 * 1.4, y0, rf"$\propto h^{p}$", va="top", fontsize=9)

    ax.set(xlabel="step size h", ylabel=r"$|\theta_1^{num}(10)-\theta_1^{exact}(10)|$",
           title="Convergence on the linearized system", ylim=(1e-16, 1e1))
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "convergence.png", dpi=150)
    plt.close(fig)


def lyapunov():
    """Largest Lyapunov exponent and first flip time along the symmetric initial condition."""
    theta, _, lam, tflip = np.loadtxt(DATA / "lyapunov" / "lyapunov_vs_theta.txt", unpack=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    ax1.plot(theta, lam, "o-", ms=3, lw=1, color="#1f77b4")
    ax1.axhline(0, color="k", lw=0.5)
    ax1.set(ylabel=r"$\lambda_{max}$ [s$^{-1}$]", title="Largest Lyapunov exponent")
    flipped = tflip > 0  # -1 = no flip within 10 s
    ax2.plot(theta[flipped], tflip[flipped], "o", ms=3.5, color="#d62728")
    ax2.set(ylabel="time of first flip [s]", title="Flip time", ylim=(0, 10))
    for ax in (ax1, ax2):
        ax.set_xlabel(r"$\theta_0$  (initial condition $\theta_1=\theta_2=\theta_0$, at rest)")
        ax.set_xlim(0, np.pi * 1.02)
    fig.tight_layout()
    fig.savefig(FIG / "lyapunov_vs_theta.png", dpi=150)
    plt.close(fig)


def maps():
    """Largest Lyapunov exponent and number of flips on a grid of initial conditions."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    ticks = [0, np.pi / 2, np.pi]
    labels = ["0", r"$\pi/2$", r"$\pi$"]

    th1, th2, lam = load_grid(DATA / "grids" / "lyapunov.txt")
    im = axes[0].pcolormesh(th2, th1, lam, cmap="viridis", shading="auto")
    fig.colorbar(im, ax=axes[0], label=r"$\lambda_{max}$ [s$^{-1}$]")
    axes[0].set_title("Largest Lyapunov exponent")

    th1, th2, flips = load_grid(DATA / "grids" / "flip_number.txt")
    im = axes[1].pcolormesh(th2, th1, flips, cmap="magma", shading="auto")
    fig.colorbar(im, ax=axes[1], label="number of flips")
    axes[1].set_title(r"Number of flips, $t_{max}=50$ s")

    for ax in axes:
        ax.set(xlabel=r"$\theta_{2,0}$", ylabel=r"$\theta_{1,0}$", aspect="equal")
        ax.set_xticks(ticks, labels)
        ax.set_yticks(ticks, labels)
    fig.tight_layout()
    fig.savefig(FIG / "maps.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    energy()
    convergence()
    if (DATA / "lyapunov" / "lyapunov_vs_theta.txt").exists():
        lyapunov()
    if (DATA / "grids" / "lyapunov.txt").exists():
        maps()
    print(f"figures written to {FIG}")
