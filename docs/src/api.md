# API Reference

```@meta
CurrentModule = SaddlePaths
```

## Public API

### Macro

```@docs
@model
```

### Types

```@docs
SaddlePath
```

### Functions

```@docs
solve_saddle_path
```

## Internal Functions

These functions are used internally by the `@model` macro but can be useful for understanding or extending the package.

### Parsing Functions

```@docs
SaddlePaths.parse_differential_equations
SaddlePaths.find_function_calls!
SaddlePaths.collect_all_symbols!
SaddlePaths.is_greek_letter
```

## Working with Model Output

### Accessing Model Components

Once a model is created with `@model`, you can access its components:

```julia
model.state     # Vector{Symbol} - state variables
model.costate   # Vector{Symbol} - costate/control variables  
model.params    # Vector{Symbol} - parameters
model.ODEs      # Vector - right-hand side expressions
```

### Example: Converting to Functions

```julia
@model example begin
    𝑑k = α * k - δ * k - c
    𝑑c = c * (α - δ - ρ)
end

# Extract symbols
k, c = example.state[1], example.costate[1]
α, δ, ρ = example.params

# Create callable functions
f_k = @eval ($k, $c, $α, $δ) -> $(example.ODEs[1])
f_c = @eval ($k, $c, $α, $δ, $ρ) -> $(example.ODEs[2])

# Evaluate at specific point
k₀, c₀ = 1.0, 0.5
α₀, δ₀, ρ₀ = 0.3, 0.1, 0.05

dk_dt = f_k(k₀, c₀, α₀, δ₀)
dc_dt = f_c(k₀, c₀, α₀, δ₀, ρ₀)
```

### Example: Symbolic Manipulation

With appropriate symbolic packages, you can manipulate the expressions:

```julia
using Symbolics

@model symbolic_model begin
    𝑑x = α * x^2 - β * x
    𝑑y = γ * x * y
end

# The ODEs are stored as Julia expressions
expr = symbolic_model.ODEs[1]  # :(α * x^2 - β * x)

# Convert to symbolic form (with Symbolics.jl)
@variables x α β
symbolic_expr = eval(expr)

# Now you can differentiate, simplify, etc.
dx = Symbolics.derivative(symbolic_expr, x)  # 2α*x - β
```

### Example: Building ODE Systems

For numerical integration with DifferentialEquations.jl:

```julia
using DifferentialEquations

@model ode_model begin
    𝑑x = -α * x + β * y
    𝑑y = γ * x - δ * y
end

# Create ODE function
function model_ode!(du, u, p, t)
    x, y = u
    α, β, γ, δ = p
    
    du[1] = -α * x + β * y  # From ode_model.ODEs[1]
    du[2] = γ * x - δ * y   # From ode_model.ODEs[2]
end

# Set up and solve
u0 = [1.0, 0.5]  # Initial conditions
p = [0.5, 0.2, 0.3, 0.4]  # Parameters
tspan = (0.0, 10.0)

prob = ODEProblem(model_ode!, u0, p, tspan)
sol = solve(prob)
```

## Extending the Package

### Custom Variable Classification

To implement custom classification logic:

```julia
function my_custom_classifier(block::Expr)
    # Parse the block
    state, costate, params, odes = SaddlePaths.parse_differential_equations(block)
    
    # Apply custom logic
    # For example, treat all variables ending in "_s" as states
    for i in 1:length(costate)
        if endswith(string(costate[i]), "_s")
            push!(state, costate[i])
            deleteat!(costate, i)
        end
    end
    
    return state, costate, params, odes
end
```

### Custom Solver Integration

To integrate with a specific solver:

```julia
function solve_with_my_solver(model::SaddlePath; kwargs...)
    # Extract model components
    n_state = length(model.state)
    n_costate = length(model.costate)
    
    # Build Jacobian symbolically or numerically
    # ...
    
    # Find steady state
    # ...
    
    # Compute stable manifold
    # ...
    
    # Integrate along saddle path
    # ...
    
    return solution
end
```

## Type Stability and Performance

The `SaddlePath` struct stores expressions as `Vector{Any}` for flexibility. For performance-critical applications, consider:

1. **Compile expressions to functions**:
```julia
function compile_model(model::SaddlePath)
    # Generate specialized functions
    state_dims = length(model.state)
    costate_dims = length(model.costate)
    
    # Return compiled function
    return (du, u, p, t) -> begin
        # Specialized implementation
    end
end
```

2. **Use generated functions** for type stability:
```julia
@generated function evaluate_ode(::Val{model}, u, p) where model
    # Generate specialized code based on model
end
```

## Error Handling

The package validates input and provides informative errors:

```julia
# This will error - no differential operator
@model bad1 begin
    k = α * k - δ * k
end
# ERROR: No differential equations found

# This will error - invalid block structure  
@model bad2 α * k - δ * k
# ERROR: Expected a block expression
```