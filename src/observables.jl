# ============================================================================
# Observables — Physical measurements and derived quantities
# ============================================================================

"""
Container for time-series observable data from a simulation.
All values are per-spin (intensive quantities).
"""
mutable struct ObservableSet
    energy::Vector{Float64}
    magnetization::Vector{Float64}
    abs_magnetization::Vector{Float64}
    energy_sq::Vector{Float64}
    mag_sq::Vector{Float64}
    mag_fourth::Vector{Float64}
    temperature::Float64
    nsites::Int

    function ObservableSet(T::Float64, N::Int; capacity::Int=1000)
        new(
            sizehint!(Float64[], capacity),
            sizehint!(Float64[], capacity),
            sizehint!(Float64[], capacity),
            sizehint!(Float64[], capacity),
            sizehint!(Float64[], capacity),
            sizehint!(Float64[], capacity),
            T, N
        )
    end
end

"""Record a measurement from the current lattice state."""
function record!(obs::ObservableSet, model::IsingModel, lat::AbstractLattice)
    N = obs.nsites
    e = energy(model, lat) / N
    m = magnetization(lat) / N
    push!(obs.energy, e)
    push!(obs.magnetization, m)
    push!(obs.abs_magnetization, abs(m))
    push!(obs.energy_sq, e^2)
    push!(obs.mag_sq, m^2)
    push!(obs.mag_fourth, m^4)
    nothing
end

"""Number of measurements taken."""
nsamples(obs::ObservableSet) = length(obs.energy)

# ── Derived thermodynamic quantities ────────────────────────────────────────

"""Mean energy per spin ⟨e⟩."""
mean_energy(obs::ObservableSet) = mean(obs.energy)

"""Mean absolute magnetization per spin ⟨|m|⟩."""
mean_abs_magnetization(obs::ObservableSet) = mean(obs.abs_magnetization)

"""Mean magnetization per spin ⟨m⟩."""
mean_magnetization(obs::ObservableSet) = mean(obs.magnetization)

"""
Specific heat per spin: Cᵥ = β²N(⟨e²⟩ - ⟨e⟩²)
                        = N/T² × Var(e)
"""
function specific_heat(obs::ObservableSet)
    β = 1.0 / obs.temperature
    N = obs.nsites
    e_mean = mean(obs.energy)
    e2_mean = mean(obs.energy_sq)
    return β^2 * N * (e2_mean - e_mean^2)
end

"""
Magnetic susceptibility per spin: χ = βN(⟨m²⟩ - ⟨|m|⟩²)
"""
function susceptibility(obs::ObservableSet)
    β = 1.0 / obs.temperature
    N = obs.nsites
    m2_mean = mean(obs.mag_sq)
    absm_mean = mean(obs.abs_magnetization)
    return β * N * (m2_mean - absm_mean^2)
end

"""
Binder cumulant: U₄ = 1 - ⟨m⁴⟩ / (3⟨m²⟩²)
Useful for locating Tc (crossing point for different L).
"""
function binder_cumulant(obs::ObservableSet)
    m2 = mean(obs.mag_sq)
    m4 = mean(obs.mag_fourth)
    return m2 ≈ 0.0 ? 0.0 : 1.0 - m4 / (3.0 * m2^2)
end

# ── Summary struct for temperature sweep results ───────────────────────────

struct SweepResult
    T::Vector{Float64}
    E::Vector{Float64}       # ⟨e⟩ per spin
    absM::Vector{Float64}    # ⟨|m|⟩ per spin
    Cv::Vector{Float64}      # specific heat
    chi::Vector{Float64}     # susceptibility
    U4::Vector{Float64}      # Binder cumulant
    L::Int
    lattice_name::String
end

"""Compute derived quantities from raw temperature sweep data."""
function compute_sweep_results(raw, L::Int, lattice_name::String)
    n = length(raw.T)
    Cv = zeros(n)
    chi = zeros(n)
    U4 = zeros(n)

    for i in 1:n
        T = raw.T[i]
        β = 1.0 / T

        # Determine N based on lattice type
        if lattice_name == "Honeycomb"
            N = 2 * L^2
        elseif lattice_name == "Kagome"
            N = 3 * L^2
        elseif lattice_name == "Cubic"
            N = L^3
        else
            N = L^2
        end

        # Specific heat: Cᵥ/N = β²N(⟨e²⟩ - ⟨e⟩²) where e is per spin
        Cv[i] = β^2 * N * (raw.E2[i] - raw.E[i]^2)

        # Susceptibility: χ/N = βN(⟨m²⟩ - ⟨|m|⟩²)
        chi[i] = β * N * (raw.M2[i] - raw.absM[i]^2)

        # Binder cumulant
        m2 = raw.M2[i]
        m4 = raw.M4[i]
        U4[i] = m2 ≈ 0.0 ? 0.0 : 1.0 - m4 / (3.0 * m2^2)
    end

    SweepResult(raw.T, raw.E, raw.absM, Cv, chi, U4, L, lattice_name)
end
