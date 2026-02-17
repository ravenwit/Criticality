# Issue: Bonito Button Styling MethodError

## Description
When launching the GUI, `Bonito.Button` throws a `MethodError: Cannot convert an object of type String to an object of type Bonito.CSS`.

## Stack Trace
```julia
MethodError: Cannot `convert` an object of type String to an object of type Bonito.CSS
Stacktrace:
  [1] convert(::Type{OrderedCollections.OrderedDict{String, Bonito.CSS}}, d::Dict{String, String})
    @ OrderedCollections ...
  [2] Styles ...
  [3] Bonito.Styles(first::Bonito.Styles, rest::Dict{String, String}) ...
  [4] jsrender(session::Bonito.Session, button::Bonito.Button)
    @ Bonito ~/.julia/packages/Bonito/OJTyF/src/widgets.jl:64
```

## Cause
The `style` argument in `Bonito.Button` was passed as a `Dict{String, String}`, but `Bonito` expects a `Bonito.Styles` object or compatible types (like `Bonito.CSS`). The default `Styles` constructor does not handle a raw `Dict{String, String}` correctly when merging with defaults.

## Fix Applied
Updated `src/gui.jl` to construct styles using `Bonito.Styles` and `Bonito.CSS` correctly.

```julia
style=Bonito.Styles(
    Bonito.CSS("background-color" => "#4CAF50"),
    Bonito.CSS("color" => "white")
)
```
