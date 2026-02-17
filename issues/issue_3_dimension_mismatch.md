# Issue: Makie DimensionMismatch (991 vs 0)

## Description
The GUI throws a `DimensionMismatch` error when plotting:
`DimensionMismatch: arrays could not be broadcast to a common size; got a dimension with lengths 991 and 0`.

## Stack Trace
```julia
DimensionMismatch: arrays could not be broadcast to a common size; got a dimension with lengths 991 and 0
Stacktrace:
  ...
  [7] convert_arguments(::MakieCore.PointBased, x::Vector{Int64}, y::Vector{Float64})
  ...
  [14] lines!(::Makie.Axis, ::Vararg{Any}; kw::@Kwargs{color::Symbol})
    @ MakieCore ...
  [15] (::IsingPhaseSim.var"#22#38"{…})(session::Bonito.Session{…})
```

## Cause
The error occurs because `lines!` is called with separate observables for x (`step_indices`) and y (`energy_history`). If these observables update asynchronously or one is cleared while the other retains data (due to a race condition or error during reset), Makie sees mismatched array lengths during the plot update. In this case, one array had 991 elements and the other 0.

## Fix Applied
Refactored `src/gui.jl` to use `Vector{Point2f}` for storing time-series data (`energy_history`, `mag_history`). This bundles (step, value) into a single atomic structure, ensuring that the x and y coordinates used by `lines!` are always consistent and of the same length. `step_indices` is removed.
