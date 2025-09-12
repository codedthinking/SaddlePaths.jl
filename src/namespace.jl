"""
    namespace.jl

Provides isolated variable namespaces for numerical solutions.

Each NumericalSolution instance gets its own submodule within Variables,
allowing multiple solutions to coexist without variable name conflicts.
"""

"""
    Variables

Parent module containing isolated namespaces for each NumericalSolution instance.
Each namespace is a dynamically generated submodule with a unique name.

# Structure
```
Variables
├── ns_abc123  # Namespace for solution 1
│   ├── k      # State variable
│   ├── c      # Costate variable
│   └── α      # Parameter
├── ns_def456  # Namespace for solution 2
│   ├── k      # Independent from solution 1's k
│   └── ...
```
"""
module Variables end

"""
    create_namespace() -> Symbol

Create a new isolated namespace within the Variables module.

Returns a unique symbol that identifies the namespace. This symbol
should be stored in the NumericalSolution struct for future reference.

# Returns
- `Symbol` - Unique identifier for the created namespace

# Example
```julia
ns_id = create_namespace()
# Returns something like :ns_7f8a9b2c

# The namespace is now available at Variables.ns_7f8a9b2c
```
"""
function create_namespace()::Symbol
    # Generate unique namespace name
    ns_name = Symbol("ns_", string(gensym())[3:end])  # Remove the ## prefix from gensym

    # Create the submodule within Variables
    Core.eval(Variables, :(module $ns_name end))

    return ns_name
end

"""
    delete_namespace(ns_id::Symbol)

Remove a namespace from the Variables module.

Should be called when a NumericalSolution is no longer needed
to free up memory and avoid namespace pollution.

# Arguments
- `ns_id::Symbol` - The namespace identifier returned by create_namespace()

# Example
```julia
delete_namespace(:ns_7f8a9b2c)
```
"""
function delete_namespace(ns_id::Symbol)
    # Julia doesn't allow true module deletion, but we can clear its contents
    # and remove references to make it eligible for garbage collection
    try
        ns = getfield(Variables, ns_id)
        for name in names(ns; all = true)
            if name != ns_id && name != :eval && name != Symbol("#eval")
                try
                    # Try to delete the binding if it exists
                    if isdefined(ns, name)
                        Base.delete_binding!(ns, name)
                    end
                catch
                    # Some names might be protected, set to nothing as fallback
                    try
                        Core.eval(ns, :($name = nothing))
                    catch
                        # Ignore errors
                    end
                end
            end
        end
    catch e
        @warn "Could not delete namespace $ns_id: $e"
    end
end

"""
    set_variable!(ns_id::Symbol, var_name::Symbol, value)

Set a variable value in the specified namespace.

# Arguments
- `ns_id::Symbol` - The namespace identifier
- `var_name::Symbol` - Name of the variable (e.g., :k, :c, :α)
- `value` - The value to assign

# Example
```julia
set_variable!(:ns_7f8a9b2c, :k, 1.5)
set_variable!(:ns_7f8a9b2c, :α, 0.3)
```
"""
function set_variable!(ns_id::Symbol, var_name::Symbol, value)
    ns = getfield(Variables, ns_id)
    Core.eval(ns, :($var_name = $value))
end

"""
    set_variable!(ns::Module, var_name::Symbol, value)

Set a variable value in a namespace module (used with with_namespace).

# Arguments
- `ns::Module` - The namespace module
- `var_name::Symbol` - Name of the variable (e.g., :k, :c, :α)
- `value` - The value to assign

# Example
```julia
with_namespace(ns_id) do ns
    set_variable!(ns, :k, 1.5)
end
```
"""
function set_variable!(ns::Module, var_name::Symbol, value)
    Core.eval(ns, :($var_name = $value))
end

"""
    get_variable(ns_id::Symbol, var_name::Symbol)

Retrieve a variable value from the specified namespace.

# Arguments
- `ns_id::Symbol` - The namespace identifier
- `var_name::Symbol` - Name of the variable to retrieve

# Returns
- The value of the variable

# Example
```julia
k_value = get_variable(:ns_7f8a9b2c, :k)
```
"""
function get_variable(ns_id::Symbol, var_name::Symbol)
    ns = getfield(Variables, ns_id)
    return getfield(ns, var_name)
end

"""
    get_variable(ns::Module, var_name::Symbol)

Retrieve a variable value from a namespace module (used with with_namespace).

# Arguments
- `ns::Module` - The namespace module
- `var_name::Symbol` - Name of the variable to retrieve

# Returns
- The value of the variable

# Example
```julia
with_namespace(ns_id) do ns
    k = get_variable(ns, :k)
end
```
"""
function get_variable(ns::Module, var_name::Symbol)
    return getfield(ns, var_name)
end

"""
    update_namespace!(ns_id::Symbol, model::SaddlePath, state_vals, costate_vals, param_vals)

Update all variables in a namespace with current values.

# Arguments
- `ns_id::Symbol` - The namespace identifier
- `model::SaddlePath` - Model specification with variable names
- `state_vals` - Vector of state variable values
- `costate_vals` - Vector of costate variable values  
- `param_vals` - Vector of parameter values

# Example
```julia
update_namespace!(ns_id, model, [1.5], [0.8], [0.3, 0.1, 0.05])
```
"""
function update_namespace!(ns_id::Symbol, model, state_vals, costate_vals, param_vals)
    # Set state variables
    for (name, val) in zip(model.state, state_vals)
        set_variable!(ns_id, name, val)
    end

    # Set costate variables
    for (name, val) in zip(model.costate, costate_vals)
        set_variable!(ns_id, name, val)
    end

    # Set parameters
    for (name, val) in zip(model.params, param_vals)
        set_variable!(ns_id, name, val)
    end
end

"""
    eval_in_namespace(ns_id::Symbol, expr::Expr)

Evaluate an expression in the specified namespace.

This is the key function that connects the parser output to the numerical backend.
It evaluates user-defined expressions (like ODEs) in the isolated namespace.

# Arguments
- `ns_id::Symbol` - The namespace identifier
- `expr::Expr` - Expression to evaluate (typically vectorized)

# Returns
- Result of the expression evaluation

# Example
```julia
# User expression: k^α - δ * k
expr = :(k^α - δ * k)
vectorized = vectorize_function_calls(expr)
result = eval_in_namespace(:ns_7f8a9b2c, vectorized)
```
"""
function eval_in_namespace(ns_id::Symbol, expr)
    ns = getfield(Variables, ns_id)
    return Core.eval(ns, expr)
end

"""
    with_namespace(f::Function, ns_id::Symbol)

Execute a function with a specific namespace context.

Useful for batch operations that need to access multiple variables.

# Arguments
- `f::Function` - Function to execute
- `ns_id::Symbol` - The namespace identifier

# Example
```julia
with_namespace(:ns_7f8a9b2c) do ns
    k = get_variable(ns, :k)
    c = get_variable(ns, :c)
    return k^2 + c
end
```
"""
function with_namespace(f::Function, ns_id::Symbol)
    ns = getfield(Variables, ns_id)
    return f(ns)
end

# Export the public interface
export Variables,
    create_namespace,
    delete_namespace,
    set_variable!,
    get_variable,
    update_namespace!,
    eval_in_namespace,
    with_namespace
