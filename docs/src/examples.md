# Examples

```@meta
CurrentModule = SaddlePaths
```

This page provides comprehensive examples of using SaddlePaths.jl for various economic models.

## Classic Economic Models

### Ramsey-Cass-Koopmans Model

The canonical optimal growth model:

```julia
@model ramsey begin
    𝑑c = c(k) * (α * k^(α - 1) - δ - ρ) / σ
    𝑑k = k^α - δ * k - c(k)
end

# Model components
# States: [:k] - capital stock
# Costates: [:c] - consumption (control variable)
# Parameters: [:α, :δ, :ρ, :σ]
```

Where:
- `k`: capital per capita
- `c`: consumption per capita  
- `α`: capital share in production
- `δ`: depreciation rate
- `ρ`: rate of time preference
- `σ`: inverse of intertemporal elasticity of substitution

### Solow-Swan Growth Model

The basic exogenous growth model:

```julia
@model solow begin
    𝑑k = s * k^α * (A * L)^(1 - α) / L - (n + g + δ) * k
end

# Simplified version with constant savings rate
@model solow_simple begin
    𝑑k = s * k^α - (n + δ) * k
end
```

### AK Model

Endogenous growth with constant returns:

```julia
@model ak_model begin
    𝑑c = c(k) * (A - δ - ρ) / σ
    𝑑k = A * k - δ * k - c(k)
end
```

## Multi-Sector Models

### Two-Sector Growth Model

```julia
@model two_sector begin
    # Capital in consumption goods sector
    𝑑k_c = α_c * k_c^(β_c) - δ_c * k_c - c(k_c, k_i)
    
    # Capital in investment goods sector  
    𝑑k_i = α_i * k_i^(β_i) - δ_i * k_i - i(k_c, k_i)
    
    # Consumption Euler equation
    𝑑c = c(k_c, k_i) * (r_c - ρ) / σ
    
    # Investment dynamics
    𝑑i = i(k_c, k_i) * (r_i - r_c)
end
```

### Uzawa-Lucas Model

Human capital accumulation:

```julia
@model uzawa_lucas begin
    # Physical capital
    𝑑k = (1 - u) * h * k^α - δ_k * k - c(k, h)
    
    # Human capital
    𝑑h = B * u * h - δ_h * h
    
    # Consumption
    𝑑c = c(k, h) * ((1 - u) * α * k^(α-1) * h - δ_k - ρ) / σ
    
    # Time allocation (control)
    𝑑u = u(k, h) * (B - (1 - α) * (1 - u) * k^α / h - δ_h + δ_k)
end
```

## Environmental and Resource Models

### Model with Pollution

```julia
@model pollution_model begin
    # Capital accumulation
    𝑑k = k^α * (1 - τ) - δ * k - c(k, p)
    
    # Pollution stock
    𝑑p = ε * k^α - η * p
    
    # Consumption with environmental quality
    𝑑c = c(k, p) * (α * k^(α-1) * (1 - τ) - δ - ρ + ψ * p) / σ
end
```

### Renewable Resource Model

```julia
@model renewable_resource begin
    # Resource stock
    𝑑R = r * R * (1 - R / K) - h(R)
    
    # Harvesting control
    𝑑h = h(R) * (p - c_h / R - ρ)
end
```

## Financial Models

### Portfolio Choice Model

```julia
@model portfolio begin
    # Wealth dynamics
    𝑑W = (r + ω(W) * (μ - r)) * W - c(W)
    
    # Consumption
    𝑑c = c(W) * ((r + ω(W) * (μ - r) - 0.5 * γ * ω(W)^2 * σ_m^2) - ρ) / γ
    
    # Portfolio weight
    𝑑ω = ω(W) * (μ - r) / (γ * σ_m^2)
end
```

### Stochastic Growth (Certainty Equivalent)

```julia
@model stochastic_growth begin
    # Expected capital dynamics
    𝑑k = k^α - δ * k - c(k) - 0.5 * σ^2 * k
    
    # Risk-adjusted consumption
    𝑑c = c(k) * (α * k^(α-1) - δ - ρ - γ * σ^2) / θ
end
```

## Game Theory and Strategic Interaction

### Differential Game (Open-Loop Nash)

```julia
@model differential_game begin
    # State: common resource
    𝑑x = f(x) - u₁(x) - u₂(x)
    
    # Player 1's control
    𝑑u₁ = u₁(x) * (p₁(x) - c₁ - λ₁)
    
    # Player 2's control  
    𝑑u₂ = u₂(x) * (p₂(x) - c₂ - λ₂)
    
    # Shadow prices
    𝑑λ₁ = λ₁ * (ρ₁ - f_x(x))
    𝑑λ₂ = λ₂ * (ρ₂ - f_x(x))
end
```

## Vintage Capital Models

### Putty-Clay Model

```julia
@model putty_clay begin
    # Aggregate capital (integral over vintages)
    𝑑K = I - δ * K
    
    # Investment in new vintage
    𝑑I = I(K, v) * (q - 1)
    
    # Shadow value of capital
    𝑑q = q * (r + δ) - F_K
    
    # Vintage-specific productivity  
    𝑑v = v(K) * technical_change
end
```

## Overlapping Generations (Continuous Time)

### Blanchard-Yaari Model

```julia
@model blanchard_yaari begin
    # Aggregate capital
    𝑑k = k^α - c_y - c_o - (n + δ) * k
    
    # Young generation consumption
    𝑑c_y = c_y(k) * (r - ρ - p) / σ
    
    # Old generation consumption  
    𝑑c_o = c_o(k) * (r + p - ρ) / σ
end
```

## Complex Dynamics Examples

### Model with Threshold Effects

```julia
@model threshold begin
    # Different dynamics based on state
    𝑑k = ifelse(k < k_bar,
                 A_L * k^α - δ_L * k - c(k),
                 A_H * k^α - δ_H * k - c(k))
    
    𝑑c = c(k) * (ifelse(k < k_bar,
                        α * A_L * k^(α-1) - δ_L,
                        α * A_H * k^(α-1) - δ_H) - ρ) / σ
end
```

### Model with Time-Varying Parameters

```julia
@model time_varying begin
    # Add time as state variable
    𝑑t = 1
    
    # Capital with time-dependent productivity
    𝑑k = A(t) * k^α - δ * k - c(k, t)
    
    # Consumption with time preference drift
    𝑑c = c(k, t) * (α * A(t) * k^(α-1) - δ - ρ(t)) / σ
end
```

## Using the Models

After defining a model, you can:

1. **Inspect the structure**:
```julia
println("State variables: ", model.state)
println("Control variables: ", model.costate)
println("Parameters: ", model.params)
```

2. **Extract equations for numerical solving**:
```julia
# Convert to function
state_eq = eval(:(($k, $c, $α, $δ) -> $(model.ODEs[1])))
```

3. **Pass to a solver** (when implemented):
```julia
solution = solve_saddle_path(model,
    params = Dict(:α => 0.3, :δ => 0.1, :ρ => 0.05, :σ => 2.0),
    initial_state = [k₀ = 1.0],
    horizon = 100
)
```

## Tips for Model Specification

1. **Use meaningful names**: Instead of single letters, consider `capital`, `consumption`
2. **Document assumptions**: Add comments about parameter ranges and model assumptions
3. **Check dimensions**: Ensure dimensional consistency in your equations
4. **Test limiting cases**: Verify your model reduces to known special cases
5. **Use consistent notation**: Stick to either functional `c(k)` or plain `c` notation