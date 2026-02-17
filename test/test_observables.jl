using Test

# Load module if not already loaded
if !isdefined(@__MODULE__, :IsingPhaseSim)
    push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
    include(joinpath(@__DIR__, "..", "src", "IsingPhaseSim.jl"))
    using .IsingPhaseSim
end

@testset "Observable & Energy Tests" begin

    # ── Energy of known configurations ───────────────────────────────────
    @testset "All-up energy (Square PBC)" begin
        L = 8
        lat = SquareLattice(L; bc=Periodic)
        all_up!(lat)
        model = IsingModel(J=1.0, h=0.0)

        # For all-up: each bond contributes -J × (+1)(+1) = -J
        # Number of bonds = z×N/2 = 4×64/2 = 128
        # E = -J × 128 = -128
        E = energy(model, lat)
        @test E ≈ -1.0 * (coordination(lat) * nsites(lat) / 2)
    end

    @testset "All-up energy (Triangular PBC)" begin
        L = 6
        lat = TriangularLattice(L; bc=Periodic)
        all_up!(lat)
        model = IsingModel(J=1.0, h=0.0)

        # Bonds = z×N/2 = 6×36/2 = 108
        E = energy(model, lat)
        @test E ≈ -1.0 * (coordination(lat) * nsites(lat) / 2)
    end

    @testset "All-up energy (Honeycomb PBC)" begin
        L = 6
        lat = HoneycombLattice(L; bc=Periodic)
        all_up!(lat)
        model = IsingModel(J=1.0, h=0.0)

        E = energy(model, lat)
        @test E ≈ -1.0 * (coordination(lat) * nsites(lat) / 2)
    end

    @testset "All-up energy (Cubic PBC)" begin
        L = 4
        lat = CubicLattice(L; bc=Periodic)
        all_up!(lat)
        model = IsingModel(J=1.0, h=0.0)

        E = energy(model, lat)
        @test E ≈ -1.0 * (coordination(lat) * nsites(lat) / 2)
    end

    @testset "All-up magnetization" begin
        lat = SquareLattice(8; bc=Periodic)
        all_up!(lat)
        @test magnetization(lat) == nsites(lat)
        @test magnetization_per_spin(lat) ≈ 1.0
    end

    # ── ΔE consistency ───────────────────────────────────────────────────
    @testset "delta_energy consistency" begin
        for name in LATTICE_NAMES
            lat = create_lattice(name, 6; bc=Periodic)
            model = IsingModel(J=1.0, h=0.0)

            E_before = energy(model, lat)
            i = rand(1:nsites(lat))
            dE = delta_energy(model, lat, i)

            # Flip the spin manually and recompute
            spins(lat)[i] *= -1
            E_after = energy(model, lat)
            spins(lat)[i] *= -1  # flip back

            @test dE ≈ (E_after - E_before) atol=1e-10
        end
    end

    # ── External field ───────────────────────────────────────────────────
    @testset "External field energy" begin
        L = 4
        lat = SquareLattice(L; bc=Periodic)
        all_up!(lat)
        model = IsingModel(J=0.0, h=1.0)  # only field, no coupling

        # E = -h × Σ sᵢ = -1.0 × 16 = -16
        @test energy(model, lat) ≈ -1.0 * nsites(lat)
    end

    # ── Metropolis step maintains valid spins ────────────────────────────
    @testset "Metropolis step validity" begin
        lat = SquareLattice(8; bc=Periodic)
        model = IsingModel(J=1.0, h=0.0)
        β = 1.0

        for _ in 1:100
            metropolis_sweep!(lat, model, β)
        end
        @test all(s -> s ∈ (-1, 1), spins(lat))
    end

    # ── Wolff step maintains valid spins ─────────────────────────────────
    @testset "Wolff step validity" begin
        lat = SquareLattice(8; bc=Periodic)
        model = IsingModel(J=1.0, h=0.0)
        β = 1.0

        for _ in 1:50
            wolff_step!(lat, model, β)
        end
        @test all(s -> s ∈ (-1, 1), spins(lat))
    end

    # ── ObservableSet ────────────────────────────────────────────────────
    @testset "ObservableSet" begin
        lat = SquareLattice(8; bc=Periodic)
        model = IsingModel(J=1.0, h=0.0)
        obs = ObservableSet(2.0, nsites(lat))

        for _ in 1:10
            metropolis_sweep!(lat, model, 0.5)
            record!(obs, model, lat)
        end

        @test nsamples(obs) == 10
        @test !isnan(mean_energy(obs))
        @test !isnan(specific_heat(obs))
        @test !isnan(susceptibility(obs))
        @test !isnan(binder_cumulant(obs))
    end
end
