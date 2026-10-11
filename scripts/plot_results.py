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

plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})


def load_grid(path):
    """Grid files: first row = (N, theta2 values), then rows (theta1, values...)."""
    raw = np.loadtxt(path)
    th2 = raw[0, 1:]
    th1 = raw[1:, 0]
    return th1, th2, raw[1:, 1:]


def energy():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

    for m in METHODS:
        t, _, _, _, dE, _ = np.loadtxt(DATA / "energy" / f"{m}_dt_1e-3.txt", unpack=True)
        ax1.semilogy(t[1:], dE[1:], color=COLORS[m], label=m)
    ax1.set(xlabel="t [s]", ylabel=r"$|E(t)-E_0|/E_0$",
            title=r"Energy error, $\Delta t=10^{-3}$, $\theta_1=\theta_2=\pi/2$")
    ax1.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.18))

    dts = [10.0**-i for i in range(2, 7)]
    for m in METHODS:
        max_err = [np.loadtxt(DATA / "energy" / f"{m}_dt_1e-{i}.txt")[:, 4].max() for i in range(2, 7)]
        ax2.loglog(dts, max_err, "o-", color=COLORS[m], label=m)
    ax2.set(xlabel=r"$\Delta t$", ylabel=r"max relative energy error over $t\in[0,10]$",
            title="Energy error vs time step")
    ax2.legend(frameon=False)

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
    theta, ftle, lam, tflip = np.loadtxt(DATA / "lyapunov" / "lyapunov_vs_theta.txt", unpack=True)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(theta, lam, "o-", ms=3, color="#1f77b4", label=r"$\lambda_{max}$ (Benettin)")
    ax.plot(theta, ftle, "s", ms=3, color="#7f7f7f", alpha=0.7, label="FTLE (single run, T=10)")
    ax.set(xlabel=r"initial angle $\theta_1=\theta_2=\theta$ [rad]", ylabel=r"Lyapunov exponent [s$^{-1}$]",
           title="Onset of chaos along the diagonal of initial conditions")
    ax.set_xticks([0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi],
                  ["0", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])
    ax.axhline(0, color="k", lw=0.5)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "lyapunov_vs_theta.png", dpi=150)
    plt.close(fig)


def maps():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    ticks = [0, np.pi / 2, np.pi]
    labels = ["0", r"$\pi/2$", r"$\pi$"]

    th1, th2, flips = load_grid(DATA / "grids" / "flip_number.txt")
    im = axes[0].pcolormesh(th2, th1, flips, cmap="magma", shading="auto")
    fig.colorbar(im, ax=axes[0], label="number of flips of the 2nd arm in 50 s")
    # energy bound: starting at rest, the 2nd arm can reach theta2 = pi only if 2cos(th1) + cos(th2) <= 1
    T2, T1 = np.meshgrid(np.linspace(0, np.pi, 400), np.linspace(0, np.pi, 400))
    axes[0].contour(T2, T1, 2 * np.cos(T1) + np.cos(T2), levels=[1], colors="w", linestyles="--", linewidths=1)
    axes[0].set_title("Flip number (dashed: energy bound)")

    th1, th2, lam = load_grid(DATA / "grids" / "lyapunov.txt")
    im = axes[1].pcolormesh(th2, th1, lam, cmap="viridis", shading="auto")
    fig.colorbar(im, ax=axes[1], label=r"$\lambda_{max}$ [s$^{-1}$]")
    axes[1].set_title("Largest Lyapunov exponent")

    for ax in axes:
        ax.set(xlabel=r"$\theta_2(0)$", ylabel=r"$\theta_1(0)$", aspect="equal")
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
