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

    # Find all function calls in RHS expressions
    function_calls = Set{Symbol}()
    foreach(eq -> find_function_calls!(function_calls, eq.second), diff_eqs)

    # Classify LHS variables based on whether they appear as functions
    costate_vars = [var for var in lhs_vars if var in function_calls]
    state_vars = [var for var in lhs_vars if var ∉ function_calls]

    # Collect all symbols from RHS
    all_rhs_symbols = Set{Symbol}()
    foreach(eq -> collect_all_symbols!(all_rhs_symbols, eq.second), diff_eqs)

    # Classify RHS symbols that aren't LHS variables
    remaining_symbols = setdiff(all_rhs_symbols, lhs_vars)
    params = sort([sym for sym in remaining_symbols if is_greek_letter(sym)])
    other_symbols = sort([sym for sym in remaining_symbols if !is_greek_letter(sym)])

    # Build ODEs array in correct order: state first, then costate
    # Create a lookup dictionary for efficiency
    eq_dict = Dict(diff_eqs)
    odes =
        vcat([eq_dict[var] for var in state_vars], [eq_dict[var] for var in costate_vars])

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
    if isfunctioncall(expr)
        if expr.head == :call
            # Regular function call f(x)
            func_name = expr.args[1]
            if func_name isa Symbol && !is_operator(func_name)
                push!(funcs, func_name)
            end
            # Recursively check arguments
            for i = 2:length(expr.args)
                find_function_calls!(funcs, expr.args[i])
            end
        elseif expr.head == Symbol(".") &&
               length(expr.args) >= 2 &&
               expr.args[2] isa Expr &&
               expr.args[2].head == :tuple
            # Vectorized function call f.(x)
            func_name = expr.args[1]
            if func_name isa Symbol && !is_operator(func_name)
                push!(funcs, func_name)
            end
            # Recursively check arguments in the tuple
            for arg in expr.args[2].args
                find_function_calls!(funcs, arg)
            end
        elseif expr.head in SYNTACTIC_OPERATORS
            # Syntactic operators like && ||
            for arg in expr.args
                find_function_calls!(funcs, arg)
            end
        end
    elseif expr isa Expr
        # Recursively check all subexpressions
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
- Arithmetic operators (using Base.isoperator)
- Special symbols: `𝑑`

# Example
```julia
expr = :(α * k^β - δ * k + exp(γ))
symbols = Set{Symbol}()
collect_all_symbols!(symbols, expr)
# symbols = Set([:α, :k, :β, :δ, :γ, :exp])
```
"""
function collect_all_symbols!(symbols::Set{Symbol}, expr)
    if expr isa Symbol
        # Skip operators and special symbols
        if !is_operator(expr) && expr != :𝑑
            push!(symbols, expr)
        end
    elseif expr isa Expr
        if expr.head == :call
            # For function calls, collect the function name and arguments
            func_name = expr.args[1]
            if func_name isa Symbol && !is_operator(func_name)
                push!(symbols, func_name)
            end
            for i = 2:length(expr.args)
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

# Greek letters used in economics
const GREEK_CHARS = Set([
    'α',
    'β',
    'γ',
    'δ',
    'ε',
    'ζ',
    'η',
    'θ',
    'ι',
    'κ',
    'λ',
    'μ',
    'ν',
    'ξ',
    'ο',
    'π',
    'ρ',
    'σ',
    'τ',
    'υ',
    'φ',
    'χ',
    'ψ',
    'ω',
    'Α',
    'Β',
    'Γ',
    'Δ',
    'Ε',
    'Ζ',
    'Η',
    'Θ',
    'Ι',
    'Κ',
    'Λ',
    'Μ',
    'Ν',
    'Ξ',
    'Ο',
    'Π',
    'Ρ',
    'Σ',
    'Τ',
    'Υ',
    'Φ',
    'Χ',
    'Ψ',
    'Ω',
])

# Syntactic operators that require special handling
const SYNTACTIC_OPERATORS = (
    :&&,
    :||,
    :+=,
    :-=,
    :*=,
    :/=,
    ://=,
    :\=,
    :^=,
    :÷=,
    :%=,
    :<<=,
    :>>=,
    :>>>=,
    :|=,
    :&=,
    :⊻=,
)

# Helper functions for expression analysis
"""
    is_operator(x) -> Bool

Check if a symbol is an operator (excluding dot reference).
"""
is_operator(::Any) = false
is_operator(x::Symbol) = Base.isoperator(x) && x != :.

"""
    is_dotted_operator(x) -> Bool

Check if a symbol is a broadcasting/vectorized operator (starts with .).
"""
is_dotted_operator(::Any) = false
is_dotted_operator(x::Symbol) =
    is_operator(x) && length(string(x)) > 0 && first(string(x)) == '.'

"""
    isfunctioncall(expr) -> Bool

Check if an expression is a function call, including operators and broadcasting.
"""
isfunctioncall(::Any) = false
isfunctioncall(x::Symbol) = is_operator(x)
isfunctioncall(x::Expr) =
    x.head == :call ||
    (
        x.head == Symbol(".") &&
        x.args[1] isa Symbol &&
        length(x.args) >= 2 &&
        x.args[2] isa Expr &&
        x.args[2].head == :tuple
    ) ||
    x.head in SYNTACTIC_OPERATORS

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
is_greek_letter(::Any) = false
is_greek_letter(sym::Symbol) = begin
    str = string(sym)
    !isempty(str) && first(str) in GREEK_CHARS
end

"""
    vectorize_function_calls(expr) -> Expr

Vectorize all function and operator calls in the expression for element-wise operations.

Transforms scalar operations to broadcasting operations by adding dots (.) to operators 
and function calls, enabling element-wise evaluation on arrays/grids.

# Arguments
- `expr` - Expression to vectorize

# Returns
- Vectorized expression with broadcasting operators

# Examples
```julia
vectorize_function_calls(:(k^α - log(c)))
# Returns: :(k.^α .- log.(c))

vectorize_function_calls(:(x + y * z))
# Returns: :(x .+ y .* z)
```

# Implementation Notes
- Operators (except syntactic ones like &&, ||) are dotted
- Function calls are dotted for broadcasting
- Already vectorized operations (e.g., f.(x)) are preserved
- Syntactic operators (&&, ||) are converted to their dotted equivalents (.&&, .||)

With the Variables module approach, the expression pipeline is simplified to:
    user_expr |> vectorize_function_calls |> (e -> Core.eval(Variables, e))
"""
function vectorize_function_calls(expr)
    # Base case: non-expression types pass through unchanged
    expr isa Expr || return expr

    # Check if this is a function call that needs vectorization
    if isfunctioncall(expr)
        # Already vectorized (broadcasting) - recurse on arguments but preserve the dot
        if expr.head == Symbol(".")
            return Expr(
                expr.head,
                vectorize_function_calls(expr.args[1]),
                Expr(:tuple, vectorize_function_calls.(expr.args[2].args)...),
            )
        end

        # Handle syntactic operators (&&, ||, etc.) - they become .&& .||
        if expr.head in SYNTACTIC_OPERATORS
            dotted_op = Symbol("." * String(expr.head))
            return Expr(dotted_op, vectorize_function_calls.(expr.args)...)
        end

        # Handle regular operators in call position: +(x, y) or x + y
        if expr.head == :call
            fname = expr.args[1]

            # Check if the function name is an operator that should be dotted
            if is_operator(fname) && !is_dotted_operator(fname)
                # Convert to dotted operator
                dotted_fname = Symbol("." * String(fname))
                return Expr(
                    :call,
                    dotted_fname,
                    vectorize_function_calls.(expr.args[2:end])...,
                )
            end

            # Regular function call f(x) becomes f.(x)
            if fname isa Symbol && !is_operator(fname)
                return Expr(
                    Symbol("."),
                    fname,
                    Expr(:tuple, vectorize_function_calls.(expr.args[2:end])...),
                )
            end

            # Otherwise recurse without changing structure
            return Expr(expr.head, fname, vectorize_function_calls.(expr.args[2:end])...)
        end
    end

    # For all other expression types, just recurse on arguments
    return Expr(expr.head, vectorize_function_calls.(expr.args)...)
end
