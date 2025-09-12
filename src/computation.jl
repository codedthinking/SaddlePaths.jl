# Import namespace management functions
using ..SaddlePaths: create_namespace, delete_namespace, set_variable!, 
                     get_variable, update_namespace!, eval_in_namespace,
                     vectorize_function_calls

"""
    NumericalSolution

Stores all components needed for iterative solution of saddle path problems.

Each instance maintains its own isolated namespace for variable evaluation,
ensuring multiple solutions can coexist without interference.

# Fields
- `model::SaddlePath` - The model specification
- `namespace::Symbol` - Unique identifier for this solution's variable namespace
- `params::Vector{<:Number}` - Parameter values
- `steady_state::Vector{<:Number}` - Steady state values
- `coefs::Vector{<:Number}` - Coefficients for approximation
- `error::AbstractFloat` - Current approximation error
- `iterations::Int` - Number of iterations performed
- `converged::Bool` - Convergence status

# Lifecycle
1. Create with `NumericalSolution(model, params)`
2. Initialize namespace automatically
3. Iterate to solve
4. Clean up with `cleanup!(solution)`

# Example
```julia
solution = NumericalSolution(model, [α=0.3, δ=0.1, ρ=0.05])
# ... perform iterations ...
cleanup!(solution)  # Releases namespace
```
"""
mutable struct NumericalSolution
    model::SaddlePath
    namespace::Symbol  # Unique namespace identifier
    params::Vector{T} where T <: Number
    steady_state::Vector{T} where T <: Number
    # grid::SpectralKit.AbstractGrid  # TODO: Add when SpectralKit is available
    # basis::SpectralKit.FunctionBasis  # TODO: Add when SpectralKit is available
    coefs::Vector{T} where T <: Number
    error::AbstractFloat
    iterations::Int
    converged::Bool
    
    """
        NumericalSolution(model::SaddlePath, param_values::Vector)
    
    Create a new numerical solution with its own isolated namespace.
    
    # Arguments
    - `model::SaddlePath` - Model specification
    - `param_values::Vector` - Parameter values in the order of model.params
    
    # Example
    ```julia
    model = @model ramsey begin
        𝑑c = c(k) * (α * k^(α - 1) - δ - ρ)
        𝑑k = k^α - δ * k - c(k)
    end
    
    solution = NumericalSolution(model, [0.3, 0.1, 0.05])  # α, δ, ρ
    ```
    """
    function NumericalSolution(model::SaddlePath, param_values::Vector{T}) where T <: Number
        # Create unique namespace for this solution
        ns_id = create_namespace()
        
        # Initialize parameter values in namespace
        for (name, val) in zip(model.params, param_values)
            set_variable!(ns_id, name, val)
        end
        
        # Create instance with default values
        new(model, ns_id, param_values, 
            zeros(T, length(model.state) + length(model.costate)),  # steady_state
            zeros(T, 0),  # coefs (empty initially)
            Inf,  # error
            0,    # iterations
            false # converged
        )
    end
end

"""
    cleanup!(solution::NumericalSolution)

Clean up resources associated with a numerical solution.

Removes the solution's namespace to free memory. Should be called
when the solution is no longer needed.

# Example
```julia
solution = NumericalSolution(model, params)
# ... use solution ...
cleanup!(solution)
```
"""
function cleanup!(solution::NumericalSolution)
    delete_namespace(solution.namespace)
end

"""
    evaluate_ode(solution::NumericalSolution, var_idx::Int, state_vals, costate_vals)

Evaluate an ODE at given state and costate values using the solution's namespace.

# Arguments
- `solution::NumericalSolution` - The numerical solution instance
- `var_idx::Int` - Index of the ODE to evaluate (1 to n_states + n_costates)
- `state_vals` - Current state variable values
- `costate_vals` - Current costate variable values

# Returns
- Result of the ODE evaluation

# Example
```julia
# Evaluate 𝑑k at current point
dk = evaluate_ode(solution, 1, [1.5], [0.8])
```
"""
function evaluate_ode(solution::NumericalSolution, var_idx::Int, state_vals, costate_vals)
    # Update namespace with current values
    update_namespace!(solution.namespace, solution.model, 
                      state_vals, costate_vals, solution.params)
    
    # Get the ODE expression
    ode_expr = solution.model.ODEs[var_idx]
    
    # Vectorize and evaluate in namespace
    vectorized_expr = vectorize_function_calls(ode_expr)
    return eval_in_namespace(solution.namespace, vectorized_expr)
end

"""
    evaluate_all_odes(solution::NumericalSolution, state_vals, costate_vals)

Evaluate all ODEs at given state and costate values.

# Returns
- Vector of ODE values [state_dots..., costate_dots...]

# Example
```julia
dots = evaluate_all_odes(solution, [1.5], [0.8])
# Returns [dk, dc] for Ramsey model
```
"""
function evaluate_all_odes(solution::NumericalSolution, state_vals, costate_vals)
    n_odes = length(solution.model.ODEs)
    results = Vector{eltype(state_vals)}(undef, n_odes)
    
    for i in 1:n_odes
        results[i] = evaluate_ode(solution, i, state_vals, costate_vals)
    end
    
    return results
end

"""
Functions to get the Chebyshev approximation and its derivative.
"""
# TODO: Implement when SpectralKit is available
# costate(s::NumericalSolution) = linear_combination(s.basis, s.coefs)
# 𝑑costate(s::NumericalSolution) = x -> linear_combination(s.basis, s.coefs, 𝑑(x))[1]

"""
Evaluate the c(k) expression using the approximation.
"""
function approx_eval(expr::Expr, var::Symbol, sol::NumericalSolution)
end

# TODO: Implement these functions
# state_dot(s::NumericalSolution)
# costate_dot(s::NumericalSolution)
