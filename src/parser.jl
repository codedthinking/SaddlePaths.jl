"""
    parse_differential_equations(block::Expr) -> (states, costates, params, odes)

Parse a block of differential equations to extract and classify variables.

# Classification Rules
1. **LHS Variables**: All variables on left-hand side (after 𝑑) are potential states/costates
2. **Costate Detection**: LHS variables that appear as function calls in RHS become costates
3. **Parameter Detection**: Symbols starting with Greek letters are parameters
4. **Other Symbols**: Non-Greek, non-state/costate symbols are global references

# Arguments
- `block::Expr` - A block expression containing differential equations

# Returns
- `state_vars::Vector{Symbol}` - Identified state variables
- `costate_vars::Vector{Symbol}` - Identified costate/control variables  
- `params::Vector{Symbol}` - Identified parameters (Greek letters)
- `odes::Vector` - Right-hand side expressions in order (states, then costates)

# Example
```julia
block = quote
    𝑑c = c(k) * (α * k^(α - 1) - δ - ρ)
    𝑑k = k^α - δ * k - c(k)
end
states, costates, params, odes = parse_differential_equations(block)
# states = [:k], costates = [:c], params = [:α, :δ, :ρ]
```
"""
function parse_differential_equations(block::Expr)
    @assert block.head == :block "Expected a block expression"
    
    # First pass: collect all differential variables and their RHS
    diff_eqs = Pair{Symbol,Any}[]
    lhs_vars = Symbol[]  # Variables on the LHS (states/costates)
    
    for expr in block.args
        # Skip line number nodes
        expr isa LineNumberNode && continue
        
        # Parse assignment expressions like 𝑑x = ...
        if expr isa Expr && expr.head == :(=)
            lhs = expr.args[1]
            rhs = expr.args[2]
            
            # Check if LHS is a differential (starts with 𝑑)
            if lhs isa Symbol
                lhs_str = string(lhs)
                if startswith(lhs_str, "𝑑")
                    # Extract variable name (remove 𝑑 prefix)
                    var_name = Symbol(lhs_str[nextind(lhs_str, 1):end])  # Skip the 𝑑 character properly
                    push!(diff_eqs, var_name => rhs)
                    push!(lhs_vars, var_name)
                end
            end
        end
    end
    
    # Second pass: analyze RHS to determine states, costates, and params
    state_vars = Symbol[]
    costate_vars = Symbol[]
    params = Symbol[]
    other_symbols = Symbol[]  # Functions or global variables
    
    # Find costate variables (those that appear as function calls like c(k))
    # FIXME: use better name, these are all function calls, not just costate
    costate_as_functions = Set{Symbol}()
    # FIXME: why is this a separate step? we only use this variable once in a for loop
    for (_, rhs) in diff_eqs
        find_function_calls!(costate_as_functions, rhs)
    end
    
    # Classify LHS variables
    # FIXME: this may even be more readable with list comprehension or filtering
    for var in lhs_vars
        if var in costate_as_functions
            push!(costate_vars, var)
        else
            push!(state_vars, var)
        end
    end
    
    # Collect all symbols from RHS and classify them
    all_rhs_symbols = Set{Symbol}()
    # FIXME: why is this a separate step? we only use this variable once in a for loop
    for (_, rhs) in diff_eqs
        collect_all_symbols!(all_rhs_symbols, rhs)
    end
    
    # Classify RHS symbols
    # FIXME: use more Julianic syntax? comprehension or guard clause or short-circuit logic
    for sym in all_rhs_symbols
        if sym in lhs_vars
            # Already classified as state or costate
            continue
        elseif is_greek_letter(sym)
            # Greek letters are parameters
            push!(params, sym)
        else
            # Everything else is a regular reference (function or global variable)
            push!(other_symbols, sym)
        end
    end

    # FIXME: this doesn't look necessary
    sort!(params)
    sort!(other_symbols)
    
    # Build ODEs array in correct order: state first, then costate
    odes = Any[]
    # FIXME: this looks way too complicated. we are extracting RHS expresssion from a list of Named Pairs
    for var in state_vars
        for eq in diff_eqs
            if eq.first == var
                push!(odes, eq.second)
                break
            end
        end
    end
    for var in costate_vars
        for eq in diff_eqs
            if eq.first == var
                push!(odes, eq.second)
                break
            end
        end
    end
    
    return state_vars, costate_vars, params, odes
end

"""
    find_function_calls!(funcs::Set{Symbol}, expr)

Recursively find all symbols that appear as function names in an expression tree.

This is used to identify costate variables, which by convention appear as functions
of state variables (e.g., `c(k)` indicates `c` is a costate).

# Arguments
- `funcs::Set{Symbol}` - Set to accumulate function names (modified in place)
- `expr` - Expression tree to search

# Example
```julia
expr = :(c(k) * α + f(x, y))
funcs = Set{Symbol}()
find_function_calls!(funcs, expr)
# funcs = Set([:c, :f])
```
"""
function find_function_calls!(funcs::Set{Symbol}, expr)
    # FIXME: i am worried, this looks too easy. 
    # - There may be different modules, like Math.exp(x)
    # - there are inline operators, like x + y
    # - (out of scope) there are pipes, like x|> log
    # - (out of scope?) there is auto-vectorization, like log.(x) 
    # BUT: these may not be relevant for costate definitions
    # we do want to be able to parse :(c(Math.log(x) + y))
    if expr isa Expr && expr.head == :call
        # This is a function call
        func_name = expr.args[1]
        if func_name isa Symbol
            # Check if this symbol also appears on LHS (is a variable, not a regular function)
            push!(funcs, func_name)
        end
        # Recursively check arguments
        for i in 2:length(expr.args)
            find_function_calls!(funcs, expr.args[i])
        end
    elseif expr isa Expr
        # Recursively check all subexpressions
        # FIXME: again, this looks suspiciously simple. check how we did this in Kezdi.jl
        for arg in expr.args
            find_function_calls!(funcs, arg)
        end
    end
end

"""
    collect_all_symbols!(symbols::Set{Symbol}, expr)

Collect all symbols from an expression tree, excluding operators and built-in functions.

# Arguments
- `symbols::Set{Symbol}` - Set to accumulate symbols (modified in place)
- `expr` - Expression tree to traverse

# Excluded Symbols
- Arithmetic operators: `+`, `-`, `*`, `/`, `^`
- Built-in functions: `exp`, `log`, `sqrt`
- Special symbols: `𝑑`

# Example
```julia
expr = :(α * k^β - δ * k + exp(γ))
symbols = Set{Symbol}()
collect_all_symbols!(symbols, expr)
# symbols = Set([:α, :k, :β, :δ, :γ])
```
"""
function collect_all_symbols!(symbols::Set{Symbol}, expr)
    if expr isa Symbol
        # Skip operators
        # FIXME: nonono! there are way more operators, check implementation in Kezdi.jl
        # FIXME: do not special case :exp, :log, :sqrt, these are function calls anyway
        if !in(expr, [:+, :-, :*, :/, :^, :𝑑, :exp, :log, :sqrt])
            push!(symbols, expr)
        end
    elseif expr isa Expr
        if expr.head == :call
            # For function calls, collect the function name and arguments
            func_name = expr.args[1]
            if func_name isa Symbol && !in(func_name, [:+, :-, :*, :/, :^, :exp, :log, :sqrt])
                push!(symbols, func_name)
            end
            for i in 2:length(expr.args)
                collect_all_symbols!(symbols, expr.args[i])
            end
        else
            # For other expressions, recursively process all arguments
            for arg in expr.args
                collect_all_symbols!(symbols, arg)
            end
        end
    end
end

"""
    is_greek_letter(sym::Symbol) -> Bool

Check if a symbol starts with a Greek letter (indicating it's a parameter).

In economics and optimization, Greek letters conventionally denote parameters
(α for productivity, δ for depreciation, ρ for discount rate, etc.).

# Arguments
- `sym::Symbol` - Symbol to check

# Returns
- `true` if the symbol starts with a Greek letter, `false` otherwise

# Examples
```julia
is_greek_letter(:α)      # true
is_greek_letter(:δ₁)     # true (starts with δ)
is_greek_letter(:beta)   # false (Latin letters)
is_greek_letter(:k)      # false
```

# Supported Greek Letters
Lowercase: α β γ δ ε ζ η θ ι κ λ μ ν ξ ο π ρ σ τ υ φ χ ψ ω
Uppercase: Α Β Γ Δ Ε Ζ Η Θ Ι Κ Λ Μ Ν Ξ Ο Π Ρ Σ Τ Υ Φ Χ Ψ Ω
"""
function is_greek_letter(sym::Symbol)
    # FIXME: it is more Julianic to use multiple dispatch as guard clause
    str = string(sym)
    isempty(str) && return false
    
    # Get first character as a string
    first_char = first(str)
    
    # Check if it's a Greek letter (lowercase or uppercase)
    # Common Greek letters used in economics
    greek_chars = Set(['α', 'β', 'γ', 'δ', 'ε', 'ζ', 'η', 'θ', 'ι', 'κ', 'λ', 'μ', 
                       'ν', 'ξ', 'ο', 'π', 'ρ', 'σ', 'τ', 'υ', 'φ', 'χ', 'ψ', 'ω',
                       'Α', 'Β', 'Γ', 'Δ', 'Ε', 'Ζ', 'Η', 'Θ', 'Ι', 'Κ', 'Λ', 'Μ',
                       'Ν', 'Ξ', 'Ο', 'Π', 'Ρ', 'Σ', 'Τ', 'Υ', 'Φ', 'Χ', 'Ψ', 'Ω'])
    
    return first_char in greek_chars
end

"""
Rewrite expression to use vector for states and costates.

# Examples

k^a - log(c) - d * k

should become

states[1]^params[1] - log(costates[1]) - params[2] * states[1]

# Notes
This has been done in Kezdi.jl, see solutions here https://github.com/codedthinking/Kezdi.jl/blob/0e37057e2b47ead1d3e249d69fc9bb2d76f8b4ca/src/codegen.jl#L148

"""
function rewrite_refs(e::Expr, m::SaddlePath)
    
end

"""
Vectorize all function and operator calls in the expression.

# Examples

k^a - log(c) - d * k

should become

states[1].^params[1] .- log.(costates[1]) .- params[2] .* states[1]

# Notes
This has been done in Kezdi.jl, see solutions here https://github.com/codedthinking/Kezdi.jl/blob/0e37057e2b47ead1d3e249d69fc9bb2d76f8b4ca/src/codegen.jl#L183

Given a user-given expression e, the actual expression to evaluate is

    e |> rewrite_refs(m) |> vectorize_function_calls(m)
"""
function vectorize_function_calls(e::Expr, m::SaddlePath)
    
end
