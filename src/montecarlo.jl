# ============================================================================
# Monte Carlo Algorithms — Metropolis & Wolff cluster
# ============================================================================

# ── Metropolis single-spin-flip ──────────────────────────────────────────────

"""
Perform one Metropolis sweep (N random single-spin-flip attempts).
Returns nothing, mutates lattice in-place.
"""
function metropolis_sweep!(lat::AbstractLattice, model::IsingModel, β::Float64)
    N = nsites(lat)
    s = spins(lat)

    # Precompute Boltzmann factors for possible ΔE values
    # For standard Ising on coordination-z lattice, ΔE ∈ {-2zJ, ..., +2zJ} in steps of 4J
    # But we compute on the fly for generality (non-uniform z under OBC)
    @inbounds for _ in 1:N
        i = rand(1:N)
        dE = delta_energy(model, lat, i)
        if dE ≤ 0.0 || rand() < exp(-β * dE)
            s[i] = -s[i]
        end
    end
    nothing
end

"""
Perform a single Metropolis step (one spin flip attempt).
"""
function metropolis_step!(lat::AbstractLattice, model::IsingModel, β::Float64)
    N = nsites(lat)
    s = spins(lat)
    i = rand(1:N)
    dE = delta_energy(model, lat, i)
    @inbounds if dE ≤ 0.0 || rand() < exp(-β * dE)
        s[i] = -s[i]
    end
    nothing
end

# ── Wolff cluster algorithm ─────────────────────────────────────────────────

"""
Perform one Wolff cluster update.
Returns the cluster size (number of flipped spins).
"""
function wolff_step!(lat::AbstractLattice, model::IsingModel, β::Float64)
    N = nsites(lat)
    s = spins(lat)
    P_add = 1.0 - exp(-2.0 * β * model.J)   # bond activation probability

    # Pick random seed
    seed = rand(1:N)
    @inbounds cluster_spin = s[seed]

    # BFS cluster growth
    cluster = falses(N)
    stack = Int[seed]
    cluster[seed] = true
    cluster_size = 0

    while !isempty(stack)
        site = pop!(stack)
        cluster_size += 1
        @inbounds for j in neighbors(lat, site)
            @inbounds if !cluster[j] && s[j] == cluster_spin && rand() < P_add
                cluster[j] = true
                push!(stack, j)
            end
        end
    end

    # Flip all spins in the cluster
    @inbounds for i in 1:N
        if cluster[i]
            s[i] = -s[i]
        end
    end

    return cluster_size
end

"""
Perform one Wolff sweep — enough cluster steps to flip ~N spins total.
"""
function wolff_sweep!(lat::AbstractLattice, model::IsingModel, β::Float64)
    N = nsites(lat)
    flipped = 0
    while flipped < N
        flipped += wolff_step!(lat, model, β)
    end
    nothing
end

# ── Unified sweep dispatcher ────────────────────────────────────────────────

"""
Perform one Monte Carlo sweep using the specified algorithm.
`algo` should be :metropolis or :wolff.
"""
function mc_sweep!(lat::AbstractLattice, model::IsingModel, β::Float64, algo::Symbol)
    if algo == :metropolis
        metropolis_sweep!(lat, model, β)
    elseif algo == :wolff
        wolff_sweep!(lat, model, β)
    else
        error("Unknown algorithm: $algo. Use :metropolis or :wolff")
    end
    nothing
end

# ── Simulation runner ────────────────────────────────────────────────────────

"""
Run a full simulation loop collecting observables.

Returns a NamedTuple with vectors of measured quantities per temperature.
"""
function run_temperature_sweep(
    lattice_name::String, L::Int, model::IsingModel;
    bc::BoundaryCondition = Periodic,
    algo::Symbol = :metropolis,
    T_range = range(1.0, 4.0, length=50),
    n_therm::Int = 5000,
    n_measure::Int = 5000,
    n_skip::Int = 5
)
    n_temps = length(T_range)
    results = (
        T = collect(T_range),
        E = zeros(n_temps),
        E2 = zeros(n_temps),
        M = zeros(n_temps),
        M2 = zeros(n_temps),
        M4 = zeros(n_temps),
        absM = zeros(n_temps),
    )

    for (ti, T) in enumerate(T_range)
        β = 1.0 / T
        lat = create_lattice(lattice_name, L; bc=bc)

        # Thermalization
        for _ in 1:n_therm
            mc_sweep!(lat, model, β, algo)
        end

        # Measurement
        N = nsites(lat)
        for k in 1:n_measure
            for _ in 1:n_skip
                mc_sweep!(lat, model, β, algo)
            end
            e = energy(model, lat) / N
            m = magnetization(lat) / N
            results.E[ti]    += e
            results.E2[ti]   += e^2
            results.M[ti]    += m
            results.M2[ti]   += m^2
            results.M4[ti]   += m^4
            results.absM[ti] += abs(m)
        end

        # Average
        results.E[ti]    /= n_measure
        results.E2[ti]   /= n_measure
        results.M[ti]    /= n_measure
        results.M2[ti]   /= n_measure
        results.M4[ti]   /= n_measure
        results.absM[ti] /= n_measure
    end

    return results
end
