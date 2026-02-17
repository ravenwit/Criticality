module IsingPhaseSim

using Statistics
using Random
using LinearAlgebra
using WGLMakie, Observables, Markdown

# ── Core simulation modules ──────────────────────────────────────────────────
include("lattices.jl")
include("ising.jl")
include("montecarlo.jl")
include("observables.jl")
include("analysis.jl")
include("gui.jl")

# ── Public API ───────────────────────────────────────────────────────────────
export
    # Lattice types & constructors
    AbstractLattice, BoundaryCondition, Periodic, Open,
    SquareLattice, TriangularLattice, HoneycombLattice, KagomeLattice, CubicLattice,
    create_lattice, LATTICE_NAMES,
    # Lattice interface
    nsites, neighbors, coordination, randomize!, all_up!, all_down!, boundary, spins,
    spin_matrix, lattice_name,
    # Ising model
    IsingModel, energy, energy_per_spin, delta_energy,
    magnetization, magnetization_per_spin, abs_magnetization_per_spin,
    # Monte Carlo
    metropolis_sweep!, metropolis_step!, wolff_step!, wolff_sweep!, mc_sweep!,
    run_temperature_sweep,
    # Observables
    ObservableSet, record!, nsamples,
    mean_energy, mean_abs_magnetization, mean_magnetization,
    specific_heat, susceptibility, binder_cumulant,
    SweepResult, compute_sweep_results,
    # Analysis
    autocorrelation, integrated_autocorrelation_time,
    spatial_correlation, correlation_length,
    finite_size_scaling,
    # GUI
    launch_gui

end # module
