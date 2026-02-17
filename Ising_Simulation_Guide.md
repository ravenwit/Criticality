# Simulating the Ising Model

Simulating the Ising model is a cornerstone of statistical mechanics, especially for understanding phase transitions and critical phenomena. Since you are working at a Master’s level, the "names" or algorithms you choose should be evaluated based on their ability to handle **critical slowing down** and their adherence to **detailed balance**.

Here are the primary Monte Carlo (MC) algorithms and theoretical frameworks used for simulating the Ising model on various lattices (square, triangular, cubic, etc.).

---

## 1. Local Update Algorithms (Metropolis-Hastings)

The most fundamental approach involves flipping a single spin at a time.

* **The Metropolis Algorithm:** The standard Markov Chain Monte Carlo (MCMC) method. A spin is flipped, and the change in energy $\Delta E$ is calculated. The flip is accepted with a probability $p = \min(1, e^{-\beta \Delta E})$, where $\beta = 1/k_B T$.
* **Glauber Dynamics:** Similar to Metropolis but uses a slightly different transition probability: $p = \frac{1}{1 + e^{\beta \Delta E}}$. It is often used to study the time-dependent evolution of the system.

> **Note on Critical Slowing Down:** At the critical temperature $T_c$, the correlation length $\xi$ diverges. Local updates become highly inefficient because the system takes a long time to explore the phase space (the "autocorrelation time" $\tau$ scales as $\tau \propto \xi^z$, where $z \approx 2$).

---

## 2. Cluster Update Algorithms

To bypass critical slowing down, cluster algorithms flip large groups of correlated spins simultaneously.

* **Wolff Algorithm:** A single-cluster update. You start with a seed spin and build a cluster by adding neighbors with a probability $P_{add} = 1 - e^{-2\beta J}$ if they have the same orientation. The entire cluster is then flipped.
* **Swendsen-Wang Algorithm:** A multi-cluster variant that partitions the entire lattice into clusters based on bond probabilities and then flips each cluster randomly with a probability of 0.5.

---

## 3. Specialized Algorithms for Different Lattices

The geometry of the lattice (coordination number $z$) significantly changes the physics, such as the location of $T_c$ or the presence of **frustration**.

* **Wang-Landau Sampling:** Instead of sampling the Boltzmann distribution, this algorithm estimates the **Density of States** $g(E)$. It is particularly powerful for complex lattices or frustrated systems (like the Antiferromagnetic Ising model on a Triangular Lattice) where standard MC gets stuck in local minima.
* **Checkerboard (Red-Black) Decomposition:** Specifically for parallelizing the Metropolis algorithm on square or cubic lattices. You color the lattice like a chessboard; since "red" spins only depend on "black" neighbors, you can update all red spins simultaneously on a GPU.

---

## 4. Tensor Network Methods (Non-MC)

If you want to move beyond stochastic simulation toward more "rigid" numerical solutions:

* **TRG (Tensor Renormalization Group):** A deterministic method used to compute the partition function by iteratively coarse-graining the lattice. It is highly accurate for 2D lattices.
* **VUMPS (Variational Uniform Matrix Product States):** Used for studying the Ising model in the thermodynamic limit (infinite lattice) by representing the state as a tensor network.

---

### Comparison of Critical Temperatures ($T_c$)

When simulating, you will likely need to verify your code against known analytical solutions for different lattices. In the units where $J=1$ and $k_B=1$:

| Lattice Type | Analytical $T_c$ Formula | Approx $T_c$ Value |
| --- | --- | --- |
| **2D Square** | $2 / \ln(1+\sqrt{2})$ | ~2.269 |
| **2D Triangular** | $4 / \ln(3)$ | ~3.641 |
| **2D Honeycomb** | $2 / \ln(2+\sqrt{3})$ | ~1.519 |
| **3D Cubic** | No exact closed-form | ~4.511 |

### Summary Checklist for your Simulation:

1. **Define the Hamiltonian:** $H = -J \sum_{\langle i,j \rangle} s_i s_j - h \sum_i s_i$.
2. **Boundary Conditions:** Usually **Periodic Boundary Conditions (PBC)** are used to minimize edge effects.
3. **Observables:** Track Magnetization $M$, Magnetic Susceptibility $\chi$, Energy $E$, and Specific Heat $C_v$.
