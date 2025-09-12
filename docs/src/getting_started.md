# Getting Started

```@meta
CurrentModule = SaddlePaths
```

## Basic Usage

The SaddlePaths package provides a macro-based DSL for defining dynamic optimization models. Here's how to get started:

### Step 1: Import the Package

```julia
using SaddlePaths
```

### Step 2: Define Your Model

Use the `@model` macro with differential equations using Unicode notation:

```julia
@model mymodel begin
    𝑑state = f(state, control, params...)
    𝑑control = g(state, control, params...)
end
```

### Step 3: Access Model Components

The macro creates a `SaddlePath` struct with classified variables:

```julia
mymodel.state     # Vector of state variables
mymodel.costate   # Vector of costate/control variables  
mymodel.params    # Vector of parameters
mymodel.ODEs      # Vector of right-hand side expressions
```

## Unicode Input

The package uses Unicode characters for mathematical notation. In the Julia REPL or compatible editors:

- Type `\itd` + TAB to get `𝑑` (differential operator)
- Type `\alpha` + TAB to get `α` (Greek letters for parameters)
- Type `\_1` + TAB to get `₁` (subscripts for indexing)

Common Greek letters used in economics:
- `α` (\alpha) - productivity, share parameters
- `β` (\beta) - discount factor  
- `γ` (\gamma) - risk aversion, growth rate
- `δ` (\delta) - depreciation rate
- `ρ` (\rho) - time preference rate
- `σ` (\sigma) - elasticity of substitution
- `θ` (\theta) - various structural parameters

## Variable Classification Rules

The DSL automatically classifies variables based on syntax patterns:

### 1. State vs Costate Detection

Variables are classified as **costates** if they appear as function calls:

```julia
@model example1 begin
    𝑑c = c(k) * r     # c(k) syntax marks c as costate
    𝑑k = f(k) - c(k)  # k is state (appears on LHS but not as function)
end
```

Without function notation, all LHS variables are states:

```julia
@model example2 begin
    𝑑c = c * r        # c is state (no function notation)
    𝑑k = f(k) - c     # k is state
end
```

### 2. Parameter Detection

Symbols starting with Greek letters are automatically classified as parameters:

```julia
@model example3 begin
    𝑑k = α * k^β - δ * k  # α, β, δ are parameters
end
```

### 3. Global Functions and Variables

Non-Greek symbols that aren't states/costates are treated as external references:

```julia
@model example4 begin
    𝑑k = s * f(k) - δ * k  # s and f are external, δ is parameter
end
```

## Common Patterns

### Growth Models

```julia
@model solow begin
    𝑑k = s * k^α - δ * k
end
```

### Optimal Control

```julia
@model optimal_control begin
    𝑑x = u(x)        # u is control (costate)
    𝑑u = -H_x        # Hamiltonian derivative
end
```

### Multi-Sector Models

```julia
@model two_sector begin
    𝑑k₁ = f₁(k₁) - c₁ - δ₁ * k₁
    𝑑k₂ = f₂(k₂) - c₂ - δ₂ * k₂
end
```

## Next Steps

- See [DSL Reference](@ref) for complete syntax documentation
- Check [Examples](@ref) for more complex models
- Read [Theory](@ref) for mathematical background on saddle path problems