"""
    SaddlePaths

A Julia module for defining and solving saddle path problems in dynamic optimization.

The module provides a domain-specific language (DSL) for specifying systems of 
differential equations with automatic classification of variables into states, 
costates, and parameters.

# Exports
- `@model` - Macro for defining saddle path models
- `SaddlePath` - Structure holding model specification
- `solve_saddle_path` - Solver for saddle path problems

# Example
```julia
using SaddlePaths

@model ramsey begin
    𝑑c = c(k) * (α * k^(α - 1) - δ - ρ)
    𝑑k = k^α - δ * k - c(k)
end
```
"""
module SaddlePaths
using SpectralKit

include("parser.jl")
include("computation.jl")

export @model, SaddlePath, solve_saddle_path

"""
    SaddlePath

Structure representing a saddle path problem with state variables, costate variables,
parameters, and ordinary differential equations (ODEs).

# Fields
- `state::Vector{Symbol}` - State variables (e.g., capital stock, assets)
- `costate::Vector{Symbol}` - Costate/control variables (e.g., consumption, shadow prices)
- `params::Vector{Symbol}` - Model parameters (Greek letters by convention)
- `ODEs::Vector` - Right-hand sides of differential equations as Julia expressions

# Example
```julia
julia> model = SaddlePath([:k], [:c], [:α, :δ, :ρ], 
                          [:(k^α - δ*k - c), :(c * (α*k^(α-1) - δ - ρ))])
```
"""
struct SaddlePath
    state::Vector{Symbol}
    costate::Vector{Symbol}
    params::Vector{Symbol}
    # FIXME: we may need to store state and co-state ODEs separately
    # FIXME: OR, we can use a different storage backend like a Dict(:k => :(Expr))
    ODEs::Vector
end


"""
    @model name begin
        𝑑var₁ = expr₁
        𝑑var₂ = expr₂
        ...
    end

Define a saddle path model using a system of differential equations.

The macro automatically classifies variables based on the following rules:
1. **States vs Costates**: Variables appearing as function calls (e.g., `c(k)`) are costates
2. **Parameters**: Symbols starting with Greek letters (α, β, γ, etc.) are parameters
3. **Other symbols**: Non-Greek symbols are treated as global functions or variables

# Arguments
- `name::Symbol` - Name for the model (will create a variable with this name)
- `block::Expr` - Block containing differential equations with `𝑑var = expr` syntax

# Returns
Creates a `SaddlePath` struct assigned to the given name.

# Examples

## Ramsey Model with Functional Notation
```julia
@model ramsey begin
    𝑑c = c(k) * (α * k^(α - 1) - δ - ρ)  # c(k) marks c as costate
    𝑑k = k^α - δ * k - c(k)
end
# Results in:
#   state = [:k]
#   costate = [:c]
#   params = [:α, :δ, :ρ]
```

## Solow Growth Model
```julia
@model solow begin
    𝑑k = s * k^α - δ * k  # s is not Greek, so not a parameter
end
# Results in:
#   state = [:k]
#   costate = []
#   params = [:α, :δ]
#   s is treated as a global variable
```

## Optimal Control Problem
```julia
@model control begin
    𝑑x = u(x) + α * x      # u(x) marks u as control/costate
    𝑑u = β * u - γ * x
end
# Results in:
#   state = [:x]
#   costate = [:u]
#   params = [:α, :β, :γ]
```

## Multiple State Variables
```julia
@model multisector begin
    𝑑k₁ = α * k₁^β - δ₁ * k₁ - c₁
    𝑑k₂ = α * k₂^β - δ₂ * k₂ - c₂
end
# Results in:
#   state = [:k₁, :k₂]  # All are states (no function call notation)
#   costate = []
#   params = [:α, :β, :δ₁, :δ₂]
```

## Using Global Functions
```julia
@model production begin
    𝑑k = f(k) - δ * k      # f is a global function
    𝑑c = c * g(k, α) - ρ * c  # g is another global function
end
# Results in:
#   state = [:k, :c]
#   params = [:α, :δ, :ρ]
#   f and g are treated as external functions
```

# Notes
- Use Unicode `𝑑` (\\itd + TAB in Julia REPL) for differentials
- Greek letters can be typed using LaTeX shortcuts (e.g., \\alpha + TAB → α)
- The order of ODEs in the output follows: states first, then costates
"""
macro model(name::Symbol, block::Expr)
    state, costate, params, odes = parse_differential_equations(block)
    
    # Create the SaddlePath struct
    quote
        $(esc(name)) = SaddlePath(
            $state,
            $costate,
            $params,
            $odes
        )
    end
end

"""
    solve_saddle_path(model::SaddlePath; kwargs...) -> Solution

Solve a saddle path problem for the given model.

This is a placeholder function that would implement numerical methods for finding
saddle-stable solutions to systems of differential equations. Common methods include:
- Shooting algorithms
- Time elimination (reverse shooting)
- Projection methods
- Relaxation algorithms

# Arguments
- `model::SaddlePath` - Model specification with states, costates, parameters, and ODEs
- `kwargs...` - Solver options (method, tolerance, initial conditions, etc.)

# Returns
- Solution object with trajectories and convergence information (placeholder)

# Example
```julia
@model ramsey begin
    𝑑c = c(k) * (α * k^(α - 1) - δ - ρ)
    𝑑k = k^α - δ * k - c(k)
end

# Would solve for optimal trajectories
solution = solve_saddle_path(ramsey, 
    initial_state = [k₀ = 1.0],
    params = Dict(:α => 0.3, :δ => 0.1, :ρ => 0.05),
    method = :shooting
)
```

# Notes
This is a placeholder implementation. A full solver would require:
- Steady state computation
- Linearization around steady state
- Eigenvalue/eigenvector analysis for stable manifold
- Numerical integration along saddle path
"""
function solve_saddle_path(model::SaddlePath; kwargs...)
    # Implementation would go here
    @info "Solving saddle path model" model.state model.costate model.params
    return nothing
end

"""
Store symbolic solutions for steady-state relationship between state variables and parameters.

# Examples
## State variable as a function of parameters only
ramsey_ss = SymbolicSolution(
    [:k],
    [:(((δ + ρ) / α)^(1/(α - 1)))]
) 

## Costate as a function of state and parameters
ramsey_ss = SymbolicSolution(
    [:c],
    [:(k^α - δ * k)]
) 

"""
struct SymbolicSolution
    variables::Vector{Symbol}
    solutions::Vector{Expr}
end

# FIXME: these quoted expressions are hard to read, lint them better
# FIXME: find more expressive name than :COSTATE for Case 3
"""
Find symbolic steady-state solutions

# Cases
## Case 1: all variables as function of parameters

solve_steady_state(ramsey_model, :ALL) = 
    SymbolicSolution(
        [:k, :c],
        [:(((δ + ρ) / α)^(1/(α - 1))),
        :((((δ + ρ) / α)^(1/(α - 1)))^α - δ * k)
        ]
    )

## Case 2: state variables as function of parameters

solve_steady_state(ramsey_model, :STATE) = 
    SymbolicSolution(
        [:k],
        [:(((δ + ρ) / α)^(1/(α - 1)))]
    )

## Case 3: co-state variables as function of states and parameters

solve_steady_state(ramsey_model, :COSTATE) = 
    SymbolicSolution(
        [:c],
        [:(k^α - δ * k)]
    ) 
"""
function solve_steady_state(model::SaddlePath)::SymbolicSolution
end


end # module