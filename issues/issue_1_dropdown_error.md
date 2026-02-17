# Issue: Bonito.Dropdown Boolean Attribute Error

## Description
When launching the GUI, `Bonito.Dropdown` throws an error regarding the `default` attribute.

## Stack Trace
```julia
ArgumentError: invalid index: nothing of type Nothing
Boolean attribute default expects a boolean! Found: String
Stacktrace:
  [1] error(s::String)
    @ Base ./error.jl:35
  ...
  [4] jsrender(session::Bonito.Session{Bonito.WebSocketConnection}, dropdown::Bonito.Dropdown)
    @ Bonito ~/.julia/packages/Bonito/OJTyF/src/widgets.jl:256
```

## Cause
The `default` keyword argument in `Bonito.Dropdown` seems to be interpreted as an HTML boolean attribute (like `<option default>`) rather than the selected value. The error "Boolean attribute default expects a boolean! Found: String" confirms this mismatch.

## Proposed Fix
Investigate `Bonito.Dropdown` API. Likely replace `default="Value"` with `index=Int` or using the observable value directly.
