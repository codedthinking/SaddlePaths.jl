# DSL Reference

```@meta
CurrentModule = SaddlePaths
```

## The @model Macro

## Syntax Rules

### Differential Equations

Each equation must follow the pattern:
```julia
𝑑variable = expression
```

Where:
- `𝑑` is the Unicode differential operator (type `\itd` + TAB)
- `variable` is a single symbol (the differentiated variable)
- `expression` is any valid Julia expression

### Variable Classification Algorithm

The classification follows this precedence:

1. **Extract LHS variables**: All variables appearing after `𝑑` on the left-hand side
2. **Identify costates**: LHS variables that appear as function calls in any RHS
3. **Identify parameters**: Symbols starting with Greek letters
4. **Remaining symbols**: Treated as global functions or variables

### Function Call Notation for Costates

To mark a variable as a costate/control, use function call syntax:

```julia
# c is marked as costate
𝑑c = c(k) * (r - ρ)  

# Without function notation, c would be a state
𝑑c = c * (r - ρ)
```

This notation is inspired by the mathematical convention where control variables are often written as functions of state variables: $c(k(t))$.

### Greek Letter Convention

Parameters are identified by their first character:

```julia
# Recognized as parameters
α, β, γ, δ, ε, ζ, η, θ, κ, λ, μ, ν, ρ, σ, τ, φ, χ, ψ, ω

# Also with subscripts
α₁, β₂, δ_k

# Not parameters (Latin letters)
a, b, c, x, y, z
```

## Advanced Features

### Multiple States and Costates

```julia
@model multi_agent begin
    𝑑k₁ = f(k₁) - c₁(k₁, k₂)
    𝑑k₂ = g(k₂) - c₂(k₁, k₂)  
    𝑑c₁ = c₁(k₁, k₂) * r₁
    𝑑c₂ = c₂(k₁, k₂) * r₂
end
# States: [:k₁, :k₂]
# Costates: [:c₁, :c₂]
```

### External Functions

Functions not defined in the model are preserved as external references:

```julia
@model production begin
    𝑑k = CES(k, l, α, σ) - δ * k
end
# CES is treated as an external function
# Only δ, α, σ are classified as parameters
```

### Complex Expressions

The DSL preserves the full Julia AST:

```julia
@model complex begin
    𝑑x = exp(-ρ * t) * (α * x^β + γ * log(x))
    𝑑y = sqrt(x) / (1 + x^2) - δ * y
end
```

## Limitations and Considerations

### Current Limitations

1. **No implicit equations**: All equations must be explicit ODEs
2. **No algebraic constraints**: Only differential equations are supported
3. **Single differential operator**: Higher-order derivatives not supported
4. **No time dependence**: Autonomous systems only (though time can be added as state)

### Best Practices

1. **Use descriptive names**: `capital` instead of `k`, `consumption` instead of `c`
2. **Document parameters**: Add comments explaining parameter meanings
3. **Consistent notation**: Use either functional or standard notation throughout
4. **Group related equations**: Keep state and costate equations together

### Example with Best Practices

```julia
@model ramsey_cass_koopmans begin
    # State equation: capital accumulation
    𝑑capital = production(capital) - depreciation * capital - consumption(capital)
    
    # Costate equation: consumption Euler equation  
    𝑑consumption = consumption(capital) * (
        marginal_product(capital) - depreciation - time_preference
    ) / intertemporal_elasticity
end

# Document parameters
# depreciation: δ ∈ (0, 1)
# time_preference: ρ > 0  
# intertemporal_elasticity: σ > 0
```

## Output Structure

The `@model` macro creates a `SaddlePath` struct:



### Accessing Components

```julia
model.state     # Vector{Symbol} of state variables
model.costate   # Vector{Symbol} of costate variables
model.params    # Vector{Symbol} of parameters
model.ODEs      # Vector of Julia expressions (AST)
```

### Expression Manipulation

Since ODEs are stored as Julia expressions, you can manipulate them:

```julia
# Evaluate with specific values
eval(model.ODEs[1])

# Convert to functions
f = eval(:(($k, $c, $α, $δ, $ρ) -> $(model.ODEs[1])))

# Symbolic differentiation (with appropriate packages)
# diff(model.ODEs[1], :k)
```