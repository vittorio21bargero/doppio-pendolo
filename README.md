# Chaos in the double pendulum

A double pendulum simulated from scratch in C++: four integrators written by hand and validated
against conservation laws and an exact solution, then used to map where in the space of initial
conditions the motion becomes chaotic.

*Course project for Numerical Algorithms, MSc in Physics of Complex Systems, University of Turin.*

![Largest Lyapunov exponent and number of flips over the plane of initial angles](figures/maps.png)

*Both arms start at rest from angles (θ₁,₀, θ₂,₀). **Left:** largest Lyapunov exponent. **Right:**
number of flips of the second arm in 50 s. Chaos does not depend only on the total energy but on
the specific combination of initial angles: an extended region of regular motion reaches deep into
the chaotic zone for θ₁,₀, θ₂,₀ ≳ π/2.*

## What is in the project

| | |
|---|---|
| **Model** | Point masses m₁ = m₂ = 1, rods L₁ = L₂ = 1, g = 9.81; equations of motion from the Lagrangian and their small-angle linearization ([`src/system.cpp`](src/system.cpp)) |
| **Integrators** | Euler, RK2 (midpoint), RK4, position Verlet ([`src/integrators.cpp`](src/integrators.cpp)) |
| **Validation** | cumulative energy error vs time step; global error against the exact normal-mode solution of the linearized system |
| **Chaos** | largest Lyapunov exponent with Benettin's algorithm, time of the first flip and number of flips ([`src/lyapunov.cpp`](src/lyapunov.cpp), [`src/flip_time.cpp`](src/flip_time.cpp)) |

## Results

**1. Energy conservation.** For an isolated pendulum the energy must stay constant, so its drift
measures the accuracy of each scheme. Euler and RK2 behave as expected for first- and second-order
methods. Position Verlet scales even more slowly than Euler: its good properties hold only when
the accelerations depend on positions alone, while here they also depend on the velocities.
RK4 is the most accurate, but its error is not monotonic in h: below h = 10⁻⁴ round-off errors
accumulate faster than truncation errors decrease, so h = 10⁻⁴ is the optimal step.

![Cumulative energy error](figures/energy_error.png)

**2. Convergence.** In the small-oscillation regime θ₁(t) is known exactly, and the error at
t = 10 s confirms the theoretical orders: 1 for Euler, 2 for RK2 and Verlet (which recovers its
order once the velocity-dependent terms become negligible), 4 for RK4, down to a floor of about
10⁻¹⁶ set by double precision.

![Convergence of the four integrators](figures/convergence.png)

**3. Transition to chaos.** With RK4 and h = 10⁻⁴, along the symmetric initial condition
θ₁ = θ₂ = θ₀ two indicators built on different principles, the divergence of nearby trajectories
(λ_max) and a macroscopic event (the first flip), identify the same threshold θ₀ ≈ 1.6, just above
π/2: below it the motion is regular and the pendulum never flips, above it λ_max becomes positive
and the flip time becomes finite but highly irregular.

![Lyapunov exponent and flip time along the symmetric initial condition](figures/lyapunov_vs_theta.png)

**4. Beyond the symmetric case.** On the full grid of initial angles (figure at the top) the
regular and chaotic regions have a non-trivial shape that the one-dimensional scan cannot capture.

## Repository structure

```
├── include/              headers
├── src/
│   ├── main.cpp          experiment selection (energy, convergence, lyapunov, grids)
│   ├── system.cpp        equations of motion, exact solution of the linearized system
│   ├── integrators.cpp   Euler, RK2, RK4, position Verlet
│   ├── simulate.cpp      integration loop with output callbacks, grids of initial conditions
│   ├── energy.cpp        mechanical energy and relative error
│   ├── lyapunov.cpp      Benettin's algorithm and FTLE
│   └── flip_time.cpp     flip time and number of flips
├── scripts/
│   └── plot_results.py   figures from the output tables
├── figures/              figures used in this README
└── Makefile
```

## Reproducibility

- **Requirements:** a C++17 compiler with OpenMP (tested with g++ 13); Python 3 with NumPy and Matplotlib for the figures.
- **Run:** `make`, then `./double_pendulum all` (or one of `energy`, `convergence`, `lyapunov`, `grids`)
  and `python scripts/plot_results.py`. Everything runs in a few minutes on a laptop.
- **Data:** none needed; the simulations write their output tables to `data/`.

## Author

Vittorio Bargero, MSc in Physics of Complex Systems, University of Turin.
