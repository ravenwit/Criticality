# ============================================================================
# Lattice Module — Defines lattice geometries and neighbor topology
# ============================================================================

# ── Boundary condition enum ──────────────────────────────────────────────────
@enum BoundaryCondition Periodic Open

# ── Abstract type ────────────────────────────────────────────────────────────
abstract type AbstractLattice end

# Common interface (documented here, implemented per concrete type):
#   nsites(lat)        → Int           total number of sites
#   neighbors(lat, i)  → Vector{Int}   neighbor indices of site i
#   coordination(lat)  → Int           max coordination number z
#   randomize!(lat)    → nothing       fill spins ±1 randomly
#   all_up!(lat)       → nothing       set all spins to +1
#   boundary(lat)      → BoundaryCondition

nsites(lat::AbstractLattice)      = length(lat.spins)
boundary(lat::AbstractLattice)    = lat.bc
spins(lat::AbstractLattice)       = lat.spins

function randomize!(lat::AbstractLattice)
    lat.spins .= rand((-Int8(1), Int8(1)), nsites(lat))
    nothing
end

function all_up!(lat::AbstractLattice)
    lat.spins .= Int8(1)
    nothing
end

function all_down!(lat::AbstractLattice)
    lat.spins .= Int8(-1)
    nothing
end

# ── Helper: PBC wrapping ────────────────────────────────────────────────────
@inline wrap(i, L) = mod1(i, L)

# ── Helper: check if index is in bounds (for OBC) ──────────────────────────
@inline inbounds(i, j, L) = (1 ≤ i ≤ L) && (1 ≤ j ≤ L)
@inline inbounds3d(i, j, k, L) = (1 ≤ i ≤ L) && (1 ≤ j ≤ L) && (1 ≤ k ≤ L)

# ────────────────────────────────────────────────────────────────────────────
#  1. SQUARE LATTICE  (2D, z = 4)
# ────────────────────────────────────────────────────────────────────────────
struct SquareLattice <: AbstractLattice
    L::Int
    spins::Vector{Int8}
    neighbors_list::Vector{Vector{Int}}
    bc::BoundaryCondition
end

coordination(::SquareLattice) = 4
lattice_name(::SquareLattice) = "Square"

function SquareLattice(L::Int; bc::BoundaryCondition=Periodic)
    N = L * L
    sp = rand((-Int8(1), Int8(1)), N)
    nb = Vector{Vector{Int}}(undef, N)

    @inline idx(i, j) = (i - 1) * L + j   # row-major linear index

    for i in 1:L, j in 1:L
        site = idx(i, j)
        nbs = Int[]
        if bc == Periodic
            push!(nbs, idx(wrap(i-1,L), j))   # up
            push!(nbs, idx(wrap(i+1,L), j))   # down
            push!(nbs, idx(i, wrap(j-1,L)))   # left
            push!(nbs, idx(i, wrap(j+1,L)))   # right
        else  # Open
            (i > 1) && push!(nbs, idx(i-1, j))
            (i < L) && push!(nbs, idx(i+1, j))
            (j > 1) && push!(nbs, idx(i, j-1))
            (j < L) && push!(nbs, idx(i, j+1))
        end
        nb[site] = nbs
    end
    SquareLattice(L, sp, nb, bc)
end

"""Return 2D coordinates (row, col) for a site index on a square lattice."""
function site_coords(lat::SquareLattice, idx::Int)
    i = div(idx - 1, lat.L) + 1
    j = mod(idx - 1, lat.L) + 1
    return (i, j)
end

"""Return a 2D matrix view of spins for visualization."""
function spin_matrix(lat::SquareLattice)
    reshape(lat.spins, lat.L, lat.L)
end

# ────────────────────────────────────────────────────────────────────────────
#  2. TRIANGULAR LATTICE  (2D, z = 6)
# ────────────────────────────────────────────────────────────────────────────
struct TriangularLattice <: AbstractLattice
    L::Int
    spins::Vector{Int8}
    neighbors_list::Vector{Vector{Int}}
    bc::BoundaryCondition
end

coordination(::TriangularLattice) = 6
lattice_name(::TriangularLattice) = "Triangular"

function TriangularLattice(L::Int; bc::BoundaryCondition=Periodic)
    N = L * L
    sp = rand((-Int8(1), Int8(1)), N)
    nb = Vector{Vector{Int}}(undef, N)

    @inline idx(i, j) = (i - 1) * L + j

    # Triangular lattice on a skewed square grid:
    # Neighbors: ±(1,0), ±(0,1), +(1,1), -(1,1)  i.e. (i+1,j+1) and (i-1,j-1)
    offsets = [(1,0), (-1,0), (0,1), (0,-1), (1,1), (-1,-1)]

    for i in 1:L, j in 1:L
        site = idx(i, j)
        nbs = Int[]
        for (di, dj) in offsets
            ni, nj = i + di, j + dj
            if bc == Periodic
                push!(nbs, idx(wrap(ni, L), wrap(nj, L)))
            else
                inbounds(ni, nj, L) && push!(nbs, idx(ni, nj))
            end
        end
        nb[site] = nbs
    end
    TriangularLattice(L, sp, nb, bc)
end

function spin_matrix(lat::TriangularLattice)
    reshape(lat.spins, lat.L, lat.L)
end

# ────────────────────────────────────────────────────────────────────────────
#  3. HONEYCOMB LATTICE  (2D, z = 3)
# ────────────────────────────────────────────────────────────────────────────
struct HoneycombLattice <: AbstractLattice
    L::Int          # unit cells per direction → N = 2L²
    spins::Vector{Int8}
    neighbors_list::Vector{Vector{Int}}
    bc::BoundaryCondition
end

coordination(::HoneycombLattice) = 3
lattice_name(::HoneycombLattice) = "Honeycomb"

function HoneycombLattice(L::Int; bc::BoundaryCondition=Periodic)
    N = 2 * L * L    # 2 sublattices: A (even) and B (odd)
    sp = rand((-Int8(1), Int8(1)), N)
    nb = Vector{Vector{Int}}(undef, N)

    # Site indexing: for unit cell (i,j), sublattice s ∈ {0,1}
    # linear index = 2*((i-1)*L + (j-1)) + s + 1
    @inline idx(i, j, s) = 2 * ((i - 1) * L + (j - 1)) + s + 1

    for i in 1:L, j in 1:L
        # Sublattice A (s=0): connects to three B neighbors
        siteA = idx(i, j, 0)
        nbA = Int[]
        # Neighbor B in same unit cell
        push!(nbA, idx(i, j, 1))
        # Neighbor B in (i-1, j) and (i, j-1) cells
        if bc == Periodic
            push!(nbA, idx(wrap(i-1, L), j, 1))
            push!(nbA, idx(i, wrap(j-1, L), 1))
        else
            (i > 1) && push!(nbA, idx(i-1, j, 1))
            (j > 1) && push!(nbA, idx(i, j-1, 1))
        end
        nb[siteA] = nbA

        # Sublattice B (s=1): connects to three A neighbors
        siteB = idx(i, j, 1)
        nbB = Int[]
        # Neighbor A in same unit cell
        push!(nbB, idx(i, j, 0))
        # Neighbor A in (i+1, j) and (i, j+1) cells
        if bc == Periodic
            push!(nbB, idx(wrap(i+1, L), j, 0))
            push!(nbB, idx(i, wrap(j+1, L), 0))
        else
            (i < L) && push!(nbB, idx(i+1, j, 0))
            (j < L) && push!(nbB, idx(i, j+1, 0))
        end
        nb[siteB] = nbB
    end
    HoneycombLattice(L, sp, nb, bc)
end

function spin_matrix(lat::HoneycombLattice)
    # For visualization: arrange into a 2L × L matrix (A and B interleaved)
    L = lat.L
    mat = zeros(Int8, 2L, L)
    for i in 1:L, j in 1:L
        sA = 2 * ((i-1)*L + (j-1)) + 1
        sB = sA + 1
        mat[2i-1, j] = lat.spins[sA]
        mat[2i,   j] = lat.spins[sB]
    end
    mat
end

# ────────────────────────────────────────────────────────────────────────────
#  4. KAGOME LATTICE  (2D, z = 4)
# ────────────────────────────────────────────────────────────────────────────
struct KagomeLattice <: AbstractLattice
    L::Int          # unit cells per direction → N = 3L²
    spins::Vector{Int8}
    neighbors_list::Vector{Vector{Int}}
    bc::BoundaryCondition
end

coordination(::KagomeLattice) = 4
lattice_name(::KagomeLattice) = "Kagome"

function KagomeLattice(L::Int; bc::BoundaryCondition=Periodic)
    N = 3 * L * L
    sp = rand((-Int8(1), Int8(1)), N)
    nb = Vector{Vector{Int}}(undef, N)

    # 3-site unit cell on triangular grid. Sites 0,1,2 per cell (i,j).
    # 0 sits at the vertex, 1 on horizontal edge, 2 on diagonal edge.
    # Kagome = corner-sharing triangles.
    @inline idx(i, j, s) = 3 * ((i - 1) * L + (j - 1)) + s + 1

    for i in 1:L, j in 1:L
        # Site 0: neighbors are 1,2 in same cell; 2 in (i,j-1); 1 in (i-1,j)
        s0 = idx(i, j, 0)
        nb0 = Int[]
        push!(nb0, idx(i, j, 1))
        push!(nb0, idx(i, j, 2))
        if bc == Periodic
            push!(nb0, idx(i, wrap(j-1,L), 2))
            push!(nb0, idx(wrap(i-1,L), j, 1))
        else
            (j > 1) && push!(nb0, idx(i, j-1, 2))
            (i > 1) && push!(nb0, idx(i-1, j, 1))
        end
        nb[s0] = nb0

        # Site 1: neighbors are 0 in same cell; 2 in same cell; 0 in (i+1,j); 2 in (i+1,j-1)
        # Actually for Kagome: site 1 connects to 0 (same), 2 (same), 0 (i+1,j), 2 in (i, j-1) — adjusted
        s1 = idx(i, j, 1)
        nb1 = Int[]
        push!(nb1, idx(i, j, 0))
        push!(nb1, idx(i, j, 2))
        if bc == Periodic
            push!(nb1, idx(wrap(i+1,L), j, 0))
            push!(nb1, idx(wrap(i+1,L), wrap(j-1,L), 2))
        else
            (i < L) && push!(nb1, idx(i+1, j, 0))
            (i < L && j > 1) && push!(nb1, idx(i+1, j-1, 2))
        end
        nb[s1] = nb1

        # Site 2: neighbors are 0 (same), 1 (same), 0 in (i,j+1), 1 in (i-1,j+1)
        s2 = idx(i, j, 2)
        nb2 = Int[]
        push!(nb2, idx(i, j, 0))
        push!(nb2, idx(i, j, 1))
        if bc == Periodic
            push!(nb2, idx(i, wrap(j+1,L), 0))
            push!(nb2, idx(wrap(i-1,L), wrap(j+1,L), 1))
        else
            (j < L) && push!(nb2, idx(i, j+1, 0))
            (i > 1 && j < L) && push!(nb2, idx(i-1, j+1, 1))
        end
        nb[s2] = nb2
    end
    KagomeLattice(L, sp, nb, bc)
end

function spin_matrix(lat::KagomeLattice)
    # For visualization: arrange into a 3L × L matrix
    L = lat.L
    mat = zeros(Int8, 3L, L)
    for i in 1:L, j in 1:L
        base = 3 * ((i-1)*L + (j-1))
        mat[3i-2, j] = lat.spins[base+1]
        mat[3i-1, j] = lat.spins[base+2]
        mat[3i,   j] = lat.spins[base+3]
    end
    mat
end

# ────────────────────────────────────────────────────────────────────────────
#  5. CUBIC LATTICE  (3D, z = 6)
# ────────────────────────────────────────────────────────────────────────────
struct CubicLattice <: AbstractLattice
    L::Int
    spins::Vector{Int8}
    neighbors_list::Vector{Vector{Int}}
    bc::BoundaryCondition
end

coordination(::CubicLattice) = 6
lattice_name(::CubicLattice) = "Cubic"

function CubicLattice(L::Int; bc::BoundaryCondition=Periodic)
    N = L * L * L
    sp = rand((-Int8(1), Int8(1)), N)
    nb = Vector{Vector{Int}}(undef, N)

    @inline idx(i, j, k) = (i - 1) * L * L + (j - 1) * L + k

    offsets = [(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]

    for i in 1:L, j in 1:L, k in 1:L
        site = idx(i, j, k)
        nbs = Int[]
        for (di, dj, dk) in offsets
            ni, nj, nk = i + di, j + dj, k + dk
            if bc == Periodic
                push!(nbs, idx(wrap(ni,L), wrap(nj,L), wrap(nk,L)))
            else
                inbounds3d(ni, nj, nk, L) && push!(nbs, idx(ni, nj, nk))
            end
        end
        nb[site] = nbs
    end
    CubicLattice(L, sp, nb, bc)
end

function spin_matrix(lat::CubicLattice)
    # For 3D visualization: show a 2D slice at k = L÷2
    L = lat.L
    k_slice = max(1, L ÷ 2)
    mat = zeros(Int8, L, L)
    for i in 1:L, j in 1:L
        idx_val = (i - 1) * L * L + (j - 1) * L + k_slice
        mat[i, j] = lat.spins[idx_val]
    end
    mat
end

# ── Unified neighbors accessor ──────────────────────────────────────────────
@inline neighbors(lat::AbstractLattice, i::Int) = lat.neighbors_list[i]

# ── Factory function for GUI ────────────────────────────────────────────────
function create_lattice(name::String, L::Int; bc::BoundaryCondition=Periodic)
    if name == "Square"
        SquareLattice(L; bc=bc)
    elseif name == "Triangular"
        TriangularLattice(L; bc=bc)
    elseif name == "Honeycomb"
        HoneycombLattice(L; bc=bc)
    elseif name == "Kagome"
        KagomeLattice(L; bc=bc)
    elseif name == "Cubic"
        CubicLattice(L; bc=bc)
    else
        error("Unknown lattice type: $name")
    end
end

const LATTICE_NAMES = ["Square", "Triangular", "Honeycomb", "Kagome", "Cubic"]
