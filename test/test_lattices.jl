using Test

# We need to load the module manually for testing
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include(joinpath(@__DIR__, "..", "src", "IsingPhaseSim.jl"))
using .IsingPhaseSim

@testset "Lattice Tests" begin

    # ── Square Lattice ───────────────────────────────────────────────────
    @testset "Square Lattice PBC" begin
        L = 8
        lat = SquareLattice(L; bc=Periodic)
        @test nsites(lat) == L^2
        @test coordination(lat) == 4
        @test boundary(lat) == Periodic
        @test lattice_name(lat) == "Square"

        # Every site should have exactly 4 neighbors
        for i in 1:nsites(lat)
            @test length(neighbors(lat, i)) == 4
        end

        # Neighbor symmetry: if j ∈ neighbors(i) then i ∈ neighbors(j)
        for i in 1:nsites(lat)
            for j in neighbors(lat, i)
                @test i ∈ neighbors(lat, j)
            end
        end

        # No self-loops
        for i in 1:nsites(lat)
            @test i ∉ neighbors(lat, i)
        end

        # spin_matrix shape
        mat = spin_matrix(lat)
        @test size(mat) == (L, L)
    end

    @testset "Square Lattice OBC" begin
        L = 8
        lat = SquareLattice(L; bc=Open)
        @test nsites(lat) == L^2
        @test boundary(lat) == Open

        # Corner sites should have 2 neighbors
        corner_idx = 1  # site (1,1)
        @test length(neighbors(lat, corner_idx)) == 2

        # Edge sites (not corners) should have 3 neighbors
        edge_idx = 2  # site (1,2) — top edge
        @test length(neighbors(lat, edge_idx)) == 3

        # Interior sites should have 4 neighbors
        interior_idx = (L ÷ 2 - 1) * L + L ÷ 2  # somewhere in the middle
        @test length(neighbors(lat, interior_idx)) == 4

        # Symmetry still holds
        for i in 1:nsites(lat)
            for j in neighbors(lat, i)
                @test i ∈ neighbors(lat, j)
            end
        end
    end

    # ── Triangular Lattice ───────────────────────────────────────────────
    @testset "Triangular Lattice PBC" begin
        L = 6
        lat = TriangularLattice(L; bc=Periodic)
        @test nsites(lat) == L^2
        @test coordination(lat) == 6

        for i in 1:nsites(lat)
            @test length(neighbors(lat, i)) == 6
        end

        # Symmetry
        for i in 1:nsites(lat)
            for j in neighbors(lat, i)
                @test i ∈ neighbors(lat, j)
            end
        end
    end

    @testset "Triangular Lattice OBC" begin
        L = 6
        lat = TriangularLattice(L; bc=Open)
        @test nsites(lat) == L^2

        # Interior sites should have 6 neighbors
        mid = (L ÷ 2 - 1) * L + L ÷ 2
        @test length(neighbors(lat, mid)) == 6

        # Symmetry
        for i in 1:nsites(lat)
            for j in neighbors(lat, i)
                @test i ∈ neighbors(lat, j)
            end
        end
    end

    # ── Honeycomb Lattice ────────────────────────────────────────────────
    @testset "Honeycomb Lattice PBC" begin
        L = 6
        lat = HoneycombLattice(L; bc=Periodic)
        @test nsites(lat) == 2 * L^2
        @test coordination(lat) == 3

        for i in 1:nsites(lat)
            @test length(neighbors(lat, i)) == 3
        end

        # Symmetry
        for i in 1:nsites(lat)
            for j in neighbors(lat, i)
                @test i ∈ neighbors(lat, j)
            end
        end
    end

    # ── Kagome Lattice ───────────────────────────────────────────────────
    @testset "Kagome Lattice PBC" begin
        L = 6
        lat = KagomeLattice(L; bc=Periodic)
        @test nsites(lat) == 3 * L^2
        @test coordination(lat) == 4

        for i in 1:nsites(lat)
            @test length(neighbors(lat, i)) == 4
        end

        # Symmetry
        for i in 1:nsites(lat)
            for j in neighbors(lat, i)
                @test i ∈ neighbors(lat, j)
            end
        end
    end

    # ── Cubic Lattice ────────────────────────────────────────────────────
    @testset "Cubic Lattice PBC" begin
        L = 4
        lat = CubicLattice(L; bc=Periodic)
        @test nsites(lat) == L^3
        @test coordination(lat) == 6

        for i in 1:nsites(lat)
            @test length(neighbors(lat, i)) == 6
        end

        # Symmetry
        for i in 1:nsites(lat)
            for j in neighbors(lat, i)
                @test i ∈ neighbors(lat, j)
            end
        end
    end

    @testset "Cubic Lattice OBC" begin
        L = 4
        lat = CubicLattice(L; bc=Open)
        @test nsites(lat) == L^3

        # Corner should have 3 neighbors
        @test length(neighbors(lat, 1)) == 3

        # Symmetry
        for i in 1:nsites(lat)
            for j in neighbors(lat, i)
                @test i ∈ neighbors(lat, j)
            end
        end
    end

    # ── Spin operations ──────────────────────────────────────────────────
    @testset "Spin operations" begin
        lat = SquareLattice(4; bc=Periodic)

        all_up!(lat)
        @test all(spins(lat) .== 1)

        all_down!(lat)
        @test all(spins(lat) .== -1)

        randomize!(lat)
        @test all(s -> s ∈ (-1, 1), spins(lat))
    end

    # ── Factory function ─────────────────────────────────────────────────
    @testset "create_lattice factory" begin
        for name in LATTICE_NAMES
            lat = create_lattice(name, 4; bc=Periodic)
            @test lat isa AbstractLattice
            @test lattice_name(lat) == name
        end
    end
end
