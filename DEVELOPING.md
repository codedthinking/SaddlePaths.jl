# SaddlePaths.jl Design Document

## 1. Goals & Scope

* **Ergonomic DSL**: `@model` macro with Unicode notation (𝒹 for time derivatives, Greek letters as parameters)
* **Symbolic core**: Build symbolic vector field $\dot{x}=f(\theta,x)$ with $x=(k,c)$ where $k$ are states and $c$ are co-states. Generate Jacobians/Hessians and JIT numerical functions using Symbolics/ModelingToolkit
* **Automatic steady state solving**: When equations use `0 =` instead of `𝒹x =`, solve for steady state symbolically/numerically
* **Stable manifold solver**: Approximate co-state policy $c(k)$ solving the HJB condition:
  $$\nabla c(k) \cdot \dot{k}(\theta,k,c(k)) = \dot{c}(\theta,k,c(k)), \quad c(k_{ss})=c_{ss}$$
  using Chebyshev bases and Smolyak sparse grids
* **Simulation**: Integrate $k'(t)=\dot{k}(\theta,k,c(k))$ and recover $c(t)=c(k(t))$ using DifferentialEquations.jl

## 2. Front-End DSL

### 2.1 Basic API

```julia
@model begin
  # States (automatically detected by 𝒹 prefix)
  𝒹k = α*k - c
  
  # Co-states must be explicitly flagged
  @costate 𝒹c = β*(α*k - c) - δ*c
end

# Compile model - system will detect states and co-states
M = compile_model()

# Analyze at parameters
A = analyze(M; θ=(α=0.3, β=0.99, δ=0.1))

# Solve for policy function
π = solve_policy(M; θ=(α=0.3, β=0.99, δ=0.1),
                 domain=:auto, order=7)

# Simulate trajectory  
traj = simulate(M, π; θ=(α=0.3, β=0.99, δ=0.1), k0=0.5, T=100.0)
```

### 2.2 Steady State Solving

The key innovation: use `0 =` to specify steady state conditions that the system should solve:

```julia
@model begin
  # Dynamic equations
  𝒹k = z*k^α - c - δ*k
  @costate 𝒹c = c*(ρ - z*α*k^(α-1) + δ)
  
  # Steady state conditions (0 = ... means solve for SS)
  0 = z*k_ss^α - c_ss - δ*k_ss
  0 = ρ - z*α*k_ss^(α-1) + δ
end

# System automatically solves for k_ss, c_ss given parameters
M = compile_model()
ss = solve_steady_state(M; θ=(α=0.36, δ=0.1, ρ=0.04, z=1.0))
# Returns: (k_ss=3.6, c_ss=0.26)
```

### 2.3 Variable Classification Rules

* **Parameters**: First character is Greek (α, β, δ, ρ, σ, etc.)
* **States**: Variables with 𝒹 prefix on LHS (𝒹k, 𝒹a, etc.)  
* **Co-states**: Variables with 𝒹 prefix preceded by `@costate` macro
* **Steady state variables**: Suffixed with `_ss` in `0 =` equations
* **Auxiliary variables**: Other Latin variables in the equations

### 2.4 Extended Example with Multiple States

```julia
@model begin
  # States
  𝒹k = z*k^α*h^(1-α) - c - δ_k*k
  𝒹h = ξ*h - δ_h*h
  
  # Co-states (Lagrange multipliers)
  @costate 𝒹λ_k = λ_k*(ρ - z*α*k^(α-1)*h^(1-α) + δ_k)
  @costate 𝒹λ_h = λ_h*(ρ - z*(1-α)*k^α*h^(-α) + δ_h - ξ)
  
  # Optimality condition
  c = λ_k^(-1/σ)
  
  # Steady state system (solved automatically)
  0 = z*k_ss^α*h_ss^(1-α) - c_ss - δ_k*k_ss
  0 = (ξ - δ_h)*h_ss
  0 = ρ - z*α*k_ss^(α-1)*h_ss^(1-α) + δ_k
  0 = c_ss - λ_k_ss^(-1/σ)
end
```

## 3. Macro Expansion & IR

The `@model` macro processes the DSL block:

1. **Parse equations**: Identify equation types
   - Dynamic: `𝒹x = ...` (state) or `@costate 𝒹x = ...` (co-state)
   - Steady state: `0 = ...` 
   - Algebraic: `x = ...`

2. **Classify variables**:
   - States: LHS of `𝒹x = ...` without `@costate`
   - Co-states: LHS of `@costate 𝒹x = ...`
   - Parameters: Greek first character
   - SS variables: Variables with `_ss` suffix in `0 =` equations

3. **Build symbolic system**:
   ```julia
   struct ModelIR
     params::Vector{Symbol}
     states::Vector{Symbol}
     costates::Vector{Symbol}
     dynamics::Dict{Symbol,Expr}      # x => ẋ expression
     steady_state_eqs::Vector{Expr}   # 0 = f(..._ss)
     algebraic_eqs::Vector{Expr}      # x = g(...)
   end
   ```

## 4. Steady State Solver

When `0 =` equations are present:

```julia
function solve_steady_state(M::Model; θ, initial_guess=nothing)
  # Extract steady state equations
  ss_vars = [k_ss, c_ss, ...]  # Variables with _ss suffix
  ss_eqs = M.steady_state_eqs   # 0 = f(ss_vars, θ)
  
  # Build nonlinear system
  F(x) = evaluate_ss_equations(ss_eqs, x, θ)
  
  # Solve using NLsolve/ModelingToolkit
  if initial_guess === nothing
    x0 = default_initial_guess(M, θ)
  end
  
  sol = nlsolve(F, x0)
  return (k_ss=sol[1], c_ss=sol[2], ...)
end
```

## 5. Stable Manifold Method

### 5.1 Problem Formulation

For states $k \in \mathbb{R}^{n_k}$ and co-states $c \in \mathbb{R}^{n_c}$:

$$\nabla c(k) \cdot \dot{k}(\theta,k,c(k)) = \dot{c}(\theta,k,c(k)), \quad c(k_{ss})=c_{ss}$$

### 5.2 Solution Strategy

1. **Linearization**: Get initial approximation from linearized dynamics
2. **Collocation**: Approximate $c(k) \approx \sum_j a_j \phi_j(k)$ using Chebyshev basis
3. **Residual minimization**: At collocation nodes $\{k_i\}$:
   $$R_i(a) = \nabla c_a(k_i) \cdot \dot{k}(\theta,k_i,c_a(k_i)) - \dot{c}(\theta,k_i,c_a(k_i))$$

## 6. Implementation Modules

```
SaddlePaths.jl/
├── src/
│   ├── SaddlePaths.jl           # Main module & exports
│   ├── dsl.jl                   # @model, @costate macros
│   ├── parser.jl                # Equation parsing & classification
│   ├── symbolic_backend.jl      # Symbolics/MTK integration
│   ├── steady_state.jl          # SS solver (0 = ... equations)
│   ├── stability.jl             # Eigenvalue analysis
│   ├── bases/
│   │   ├── chebyshev.jl        # Chebyshev polynomials
│   │   └── smolyak.jl          # Sparse grids
│   ├── policy.jl                # Stable manifold solver
│   └── simulate.jl              # ODE integration
```

## 7. Key Algorithms

### 7.1 DSL Parser with Co-state Detection

```julia
function parse_model(block)
  states = Symbol[]
  costates = Symbol[]
  dynamics = Dict{Symbol,Any}()
  ss_equations = []
  
  for expr in block.args
    if is_macro_call(expr, :costate)
      # Extract 𝒹x = rhs from @costate 𝒹x = rhs
      eq = expr.args[2]
      lhs, rhs = eq.args
      var = extract_var_from_deriv(lhs)  # 𝒹c → c
      push!(costates, var)
      dynamics[var] = rhs
      
    elseif is_derivative_eq(expr)
      # Regular state: 𝒹k = rhs
      lhs, rhs = expr.args
      var = extract_var_from_deriv(lhs)
      push!(states, var)
      dynamics[var] = rhs
      
    elseif is_steady_state_eq(expr)
      # Steady state: 0 = f(..._ss)
      push!(ss_equations, expr.args[2])
    end
  end
  
  return ModelIR(states, costates, dynamics, ss_equations)
end
```

### 7.2 Automatic Steady State Solution

```julia
function compile_steady_state_solver(ss_equations, params)
  # Convert symbolic equations to numerical function
  ss_vars = extract_ss_variables(ss_equations)
  
  # Build residual function
  function F!(resid, x, p)
    # Substitute x into ss_vars, p into params
    # Evaluate each equation, store in resid
  end
  
  # Return solver
  return (θ) -> begin
    prob = NonlinearProblem(F!, initial_guess, θ)
    sol = solve(prob)
    return NamedTuple(ss_vars .=> sol.u)
  end
end
```

## 8. Usage Examples

### 8.1 Neoclassical Growth Model

```julia
@model begin
  # Dynamics
  𝒹k = k^α - c - δ*k
  @costate 𝒹c = c*(ρ + δ - α*k^(α-1))
  
  # Steady state (solved automatically)
  0 = k_ss^α - c_ss - δ*k_ss
  0 = ρ + δ - α*k_ss^(α-1)
end

M = compile_model()
ss = solve_steady_state(M; θ=(α=0.36, δ=0.1, ρ=0.04))
π = solve_policy(M; θ=(α=0.36, δ=0.1, ρ=0.04))
```

### 8.2 Two-Sector Model

```julia
@model begin
  # States
  𝒹k_1 = i_1 - δ*k_1
  𝒹k_2 = i_2 - δ*k_2
  
  # Co-states (shadow prices)
  @costate 𝒹q_1 = q_1*(r + δ) - α*A_1*k_1^(α-1)
  @costate 𝒹q_2 = q_2*(r + δ) - α*A_2*k_2^(α-1)
  
  # Optimality
  q_1 = q_2  # Investment arbitrage
  c + i_1 + i_2 = A_1*k_1^α + A_2*k_2^α  # Resource constraint
  
  # Steady state system
  0 = δ*k_1_ss - i_1_ss
  0 = δ*k_2_ss - i_2_ss
  0 = r + δ - α*A_1*k_1_ss^(α-1)
  0 = r + δ - α*A_2*k_2_ss^(α-1)
end
```

## 9. Dependencies

- **Symbolics.jl / ModelingToolkit.jl**: Symbolic math & code generation
- **NLsolve.jl / NonlinearSolve.jl**: For steady state solving
- **DifferentialEquations.jl**: ODE integration
- **FastChebInterp.jl / ApproxFun.jl**: Chebyshev approximation
- **SparseGrids.jl**: Smolyak sparse grids (optional)