# SaddlePaths.jl Development Guide

## Project Overview

SaddlePaths.jl is a domain-specific language (DSL) for defining and solving saddle path problems in dynamic optimization, commonly used in economic modeling.

## Package Structure

```
SaddlePaths.jl/
├── src/
│   ├── SaddlePaths.jl         # Main module (exports public API)
│   ├── parser.jl              # DSL Parser implementation
│   └── computation.jl         # Numerical solver framework
├── test/
│   └── runtests.jl           # Test suite
├── docs/
│   ├── make.jl               # Documentation builder
│   └── src/
│       └── index.md          # Documentation homepage
├── Project.toml              # Package metadata and dependencies
├── Makefile                  # Build automation
└── README.md                 # User-facing documentation
```

## Module Structure

### SaddlePaths.jl (Main Module)
The entry point that exports the public API:
- `@model` macro for defining differential equation systems
- `SaddlePath` struct for model representation
- `solve_saddle_path` function for numerical solutions
- Integrates with SpectralKit.jl for numerical methods

### parser.jl (DSL Parser)
Implements the parsing logic for the `@model` macro:

**Core Functions:**
- `parse_differential_equations()` - Main parser that classifies variables from differential equations
- `find_function_calls!()` - Identifies costate variables via functional notation (e.g., `c(k)`)
- `collect_all_symbols!()` - Extracts all symbols from expressions
- `is_greek_letter()` - Detects Greek symbols as parameters

**Variable Classification Rules:**
1. Variables in function calls like `c(k)` become costates (control variables)
2. Greek letters (α, β, δ, etc.) are automatically parameters
3. Other left-hand side variables are state variables
4. Non-Greek, non-differential symbols are external references

**Planned Features:**
- `rewrite_refs()` - Transform symbolic expressions to use vector indexing
- `vectorize_function_calls()` - Add broadcasting for vectorized operations

### computation.jl (Numerical Solver)
Defines the numerical solution framework:

**Data Structures:**
- `NumericalSolution` - Mutable struct storing model, parameters, steady state, spectral grids, and convergence info
- `SymbolicSolution` - Stores analytical steady-state relationships

**Key Functions:**
- `costate()` - Returns Chebyshev approximation function
- `𝑑costate()` - Computes derivative of approximation
- `solve_steady_state()` - Finds symbolic steady-state solutions

**Solution Methods (Planned):**
- Spectral methods via SpectralKit integration
- Chebyshev polynomial approximation on sparse grids
- Iterative convergence tracking

## Usage Example

```julia
using SaddlePaths

# Define a Ramsey growth model
@model ramsey begin
    𝑑c = c(k) * (α * k^(α - 1) - δ - ρ)  # c(k) marks c as costate
    𝑑k = k^α - δ * k - c(k)              # k is state variable
end

# Model structure created:
# - state = [:k]
# - costate = [:c] 
# - params = [:α, :δ, :ρ]
```

## Design Philosophy

The module provides an elegant DSL that:
- Automatically classifies economic variables by convention
- Separates model specification from numerical implementation
- Integrates with high-performance spectral methods
- Maintains clean separation between parsing, symbolic manipulation, and computation

## Technical Notes

**Code Quality Issues to Address:**
- Parser uses overly simple recursion for expression traversal
- Missing support for module-qualified functions (e.g., `Math.exp`)
- Operator detection is incomplete (only basic arithmetic)
- Complex list comprehensions could improve readability
- Need better handling of vectorized operations

**Integration:**
- Designed to work with SpectralKit.jl for spectral methods
- Placeholder solver implementation awaiting full numerical backend
- Supports both symbolic and numerical solution approaches

## Testing

Run tests with:
```bash
make test
```

## Documentation

Build documentation with:
```bash
make docs
```

## Dependencies

- SpectralKit.jl (for spectral methods)
- (Additional dependencies to be determined)