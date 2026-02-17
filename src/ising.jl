# ============================================================================
# Ising Model — Hamiltonian, local energy, magnetization
# ============================================================================

struct IsingModel
    J::Float64    # coupling constant (positive = ferromagnetic)
    h::Float64    # external magnetic field
end

IsingModel(; J=1.0, h=0.0) = IsingModel(J, h)

"""
Total energy: H = -J Σ_{⟨ij⟩} sᵢsⱼ - h Σᵢ sᵢ
Each bond counted once.
"""
function energy(model::IsingModel, lat::AbstractLattice)
    E_bond = 0.0
    E_field = 0.0
    N = nsites(lat)
    s = spins(lat)
    @inbounds for i in 1:N
        si = s[i]
        E_field += si
        for j in neighbors(lat, i)
            if j > i   # count each bond once
                E_bond += si * s[j]
            end
        end
    end
    return -model.J * E_bond - model.h * E_field
end

"""Energy per spin."""
energy_per_spin(model::IsingModel, lat::AbstractLattice) = energy(model, lat) / nsites(lat)

"""
Local energy change ΔE for flipping spin at site i:
  ΔE = 2 sᵢ (J Σⱼ sⱼ + h)
"""
function delta_energy(model::IsingModel, lat::AbstractLattice, i::Int)
    s = spins(lat)
    @inbounds si = s[i]
    nb_sum = 0
    @inbounds for j in neighbors(lat, i)
        nb_sum += s[j]
    end
    return 2 * si * (model.J * nb_sum + model.h)
end

"""Total magnetization M = Σᵢ sᵢ."""
function magnetization(lat::AbstractLattice)
    return sum(spins(lat))
end

"""Magnetization per spin."""
magnetization_per_spin(lat::AbstractLattice) = magnetization(lat) / nsites(lat)

"""Absolute magnetization per spin."""
abs_magnetization_per_spin(lat::AbstractLattice) = abs(magnetization(lat)) / nsites(lat)
