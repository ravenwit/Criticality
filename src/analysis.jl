# ============================================================================
# Critical Analysis — Autocorrelation, Correlation Length, FSS
# ============================================================================

using FFTW

# ── Autocorrelation ──────────────────────────────────────────────────────────

"""
Compute the normalized autocorrelation function C(t) of a time series x
using the FFT method (O(N log N) instead of O(N²)).

C(t) = (⟨x_k x_{k+t}⟩ - ⟨x⟩²) / (⟨x²⟩ - ⟨x⟩²)
"""
function autocorrelation(x::Vector{Float64}; max_lag::Int=0)
    N = length(x)
    max_lag = max_lag > 0 ? min(max_lag, N - 1) : N ÷ 4

    # Mean-center the data
    x_centered = x .- mean(x)
    var_x = sum(x_centered .^ 2) / N

    if var_x ≈ 0.0
        return zeros(max_lag + 1)
    end

    # FFT-based autocorrelation
    # Pad to avoid circular correlation artifacts
    n_padded = nextpow(2, 2N)
    x_padded = zeros(n_padded)
    x_padded[1:N] .= x_centered

    fx = fft(x_padded)
    acf_raw = real.(ifft(fx .* conj.(fx)))

    # Normalize
    C = zeros(max_lag + 1)
    for t in 0:max_lag
        C[t+1] = acf_raw[t+1] / (var_x * N)
    end

    return C
end

"""
Compute the integrated autocorrelation time:
  τ_int = 1/2 + Σ_{t=1}^{t_cut} C(t)

The sum is cut off when C(t) < 0 or |C(t)| < threshold to avoid noise.
"""
function integrated_autocorrelation_time(x::Vector{Float64}; max_lag::Int=0, threshold::Float64=0.05)
    C = autocorrelation(x; max_lag=max_lag)
    τ = 0.5
    for t in 2:length(C)
        if C[t] < threshold
            break
        end
        τ += C[t]
    end
    return τ
end

# ── Spatial Correlation Function ────────────────────────────────────────────

"""
Compute the spatial spin-spin correlation function G(r) for a 2D lattice.
Uses the row direction for distance r (averaged over columns and rows).

G(r) = ⟨sᵢ s_{i+r}⟩ - ⟨s⟩²

Works for Square and Triangular lattices. For Honeycomb/Kagome, uses
sublattice-A sites only.
"""
function spatial_correlation(lat::AbstractLattice)
    s = spins(lat)
    L = lat.L

    if lat isa CubicLattice
        return _spatial_corr_cubic(lat)
    end

    # Determine which sites to use and the effective matrix
    if lat isa HoneycombLattice
        # Use only sublattice A sites (every other site)
        N_eff = lat.L^2
        mat = zeros(Int8, L, L)
        for i in 1:L, j in 1:L
            idx = 2 * ((i-1)*L + (j-1)) + 1  # sublattice A
            mat[i, j] = s[idx]
        end
    elseif lat isa KagomeLattice
        # Use only sublattice 0 sites
        N_eff = lat.L^2
        mat = zeros(Int8, L, L)
        for i in 1:L, j in 1:L
            idx = 3 * ((i-1)*L + (j-1)) + 1  # sublattice 0
            mat[i, j] = s[idx]
        end
    else
        mat = spin_matrix(lat)
    end

    r_max = L ÷ 2
    G = zeros(r_max + 1)
    m_avg = mean(Float64.(mat))

    for r in 0:r_max
        corr = 0.0
        count = 0
        for i in 1:L, j in 1:L
            # Correlation along row direction with PBC
            j2 = mod1(j + r, L)
            corr += mat[i, j] * mat[i, j2]
            count += 1
        end
        G[r+1] = corr / count - m_avg^2
    end

    return G
end

function _spatial_corr_cubic(lat::CubicLattice)
    L = lat.L
    s = spins(lat)
    r_max = L ÷ 2
    G = zeros(r_max + 1)
    N = L^3
    m_avg = mean(Float64.(s))

    @inline idx(i, j, k) = (i-1)*L*L + (j-1)*L + k

    for r in 0:r_max
        corr = 0.0
        count = 0
        for i in 1:L, j in 1:L, k in 1:L
            k2 = mod1(k + r, L)
            corr += s[idx(i,j,k)] * s[idx(i,j,k2)]
            count += 1
        end
        G[r+1] = corr / count - m_avg^2
    end
    return G
end

"""
Estimate the correlation length ξ from the spatial correlation function
using the second-moment definition:

  ξ = (1/2sin(π/L)) × √(G̃(0)/G̃(q_min) - 1)

where G̃(q) is the Fourier transform of G(r) and q_min = 2π/L.
Falls back to exponential fit if second-moment method fails.
"""
function correlation_length(lat::AbstractLattice)
    G = spatial_correlation(lat)
    L = lat.L
    r_max = length(G) - 1

    # Second-moment method
    G_tilde_0 = sum(G)
    q_min = 2π / L
    G_tilde_q = 0.0
    for r in 0:r_max
        G_tilde_q += G[r+1] * cos(q_min * r)
    end

    if G_tilde_q > 0 && G_tilde_0 / G_tilde_q > 1
        ratio = G_tilde_0 / G_tilde_q - 1.0
        ξ = sqrt(ratio) / (2.0 * sin(π / L))
        return ξ
    end

    # Fallback: simple exponential fit from first few points
    # Find where G(r) first drops below G(0)/e
    if G[1] > 0
        threshold = G[1] / exp(1)
        for r in 1:r_max
            if G[r+1] < threshold
                # Linear interpolation
                ξ = r - 1 + (G[r] - threshold) / (G[r] - G[r+1] + 1e-20)
                return max(ξ, 0.1)
            end
        end
    end

    return Float64(r_max)
end

# ── Finite-Size Scaling ─────────────────────────────────────────────────────

"""
Finite-size scaling analysis.
Given sweep results for multiple system sizes, extract critical exponents.

Returns a NamedTuple with:
- Tc_estimate: estimated critical temperature from Binder crossing
- beta_nu:  β/ν from magnetization scaling
- gamma_nu: γ/ν from susceptibility scaling
"""
function finite_size_scaling(results::Vector{SweepResult})
    if length(results) < 2
        error("Need at least 2 system sizes for FSS")
    end

    # Sort by L
    sorted = sort(results, by=r -> r.L)

    # ── Estimate Tc from Binder cumulant crossings ──
    Tc_estimates = Float64[]
    for i in 1:length(sorted)-1
        r1, r2 = sorted[i], sorted[i+1]
        # Find crossing point of U4 curves
        tc = _find_binder_crossing(r1, r2)
        !isnan(tc) && push!(Tc_estimates, tc)
    end
    Tc = isempty(Tc_estimates) ? NaN : mean(Tc_estimates)

    # ── Extract β/ν from |m(Tc)| ~ L^{-β/ν} ──
    beta_nu = NaN
    if !isnan(Tc)
        logL = Float64[]
        logM = Float64[]
        for r in sorted
            idx = argmin(abs.(r.T .- Tc))
            m_tc = r.absM[idx]
            if m_tc > 0
                push!(logL, log(r.L))
                push!(logM, log(m_tc))
            end
        end
        if length(logL) ≥ 2
            # Linear regression: log|m| = -β/ν × logL + const
            beta_nu = -_linear_slope(logL, logM)
        end
    end

    # ── Extract γ/ν from χ(Tc) ~ L^{γ/ν} ──
    gamma_nu = NaN
    if !isnan(Tc)
        logL = Float64[]
        logChi = Float64[]
        for r in sorted
            idx = argmin(abs.(r.T .- Tc))
            chi_tc = r.chi[idx]
            if chi_tc > 0
                push!(logL, log(r.L))
                push!(logChi, log(chi_tc))
            end
        end
        if length(logL) ≥ 2
            gamma_nu = _linear_slope(logL, logChi)
        end
    end

    return (Tc=Tc, beta_nu=beta_nu, gamma_nu=gamma_nu)
end

"""Find Binder cumulant crossing between two system sizes."""
function _find_binder_crossing(r1::SweepResult, r2::SweepResult)
    # Interpolate to common T grid
    T_min = max(minimum(r1.T), minimum(r2.T))
    T_max = min(maximum(r1.T), maximum(r2.T))

    if T_min ≥ T_max
        return NaN
    end

    # Find where U4_1(T) - U4_2(T) changes sign
    for i in 1:length(r1.T)-1
        T = r1.T[i]
        if T < T_min || T > T_max
            continue
        end
        # Find closest index in r2
        j = argmin(abs.(r2.T .- T))
        T_next = r1.T[i+1]
        j_next = argmin(abs.(r2.T .- T_next))

        diff_now = r1.U4[i] - r2.U4[j]
        diff_next = r1.U4[i+1] - r2.U4[j_next]

        if diff_now * diff_next < 0  # sign change
            # Linear interpolation for crossing
            frac = diff_now / (diff_now - diff_next)
            return T + frac * (T_next - T)
        end
    end
    return NaN
end

"""Simple linear regression slope."""
function _linear_slope(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    x_mean = mean(x)
    y_mean = mean(y)
    num = sum((x .- x_mean) .* (y .- y_mean))
    den = sum((x .- x_mean) .^ 2)
    return den ≈ 0.0 ? NaN : num / den
end
