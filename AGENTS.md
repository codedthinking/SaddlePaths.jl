# Julia Source Code Summary

This directory contains the **SaddlePaths.jl** module, a domain-specific language (DSL) for defining and solving saddle path problems in dynamic optimization, commonly used in economic modeling.

## Package Structure

```
SaddlePaths.jl/
├── src/
│   ├── SaddlePaths.jl         # Main module (exports public API)
│   ├── parser.jl              # DSL Parser implementation
│   └── computation.jl         # Numerical solver framework
├── test/
│   ├── runtests.jl           # Test suite
│   ├── test_saddlepaths.jl  # Model parsing tests
│   └── show_docs.jl          # Documentation display utility
├── docs/
│   ├── make.jl               # Documentation builder
│   └── src/
│       ├── index.md          # Documentation homepage
│       ├── getting_started.md # Quick start guide
│       ├── dsl.md            # DSL reference
│       ├── examples.md       # Usage examples
│       ├── api.md            # API reference
│       └── theory.md         # Mathematical theory
├── scrap/                     # Temporary test files and code snippets (not persisted)
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
- Includes parser.jl and computation.jl modules

### parser.jl (DSL Parser)
Implements the parsing logic for the `@model` macro (fully implemented):

**Core Functions:**
- `parse_differential_equations()` - Main parser that classifies variables from differential equations
- `find_function_calls!()` - Identifies costate variables via functional notation (e.g., `c(k)`)
- `collect_all_symbols!()` - Extracts all symbols from expressions
- `is_greek_letter()` - Detects Greek symbols as parameters
- `extract_lhs_variable()` - Extracts variable name from left-hand side of equation
- `classify_symbols()` - Final classification of symbols into categories

**Variable Classification Rules:**
1. Variables in function calls like `c(k)` become costates (control variables)
2. Greek letters (α, β, δ, etc.) are automatically parameters
3. Other left-hand side variables are state variables
4. Non-Greek, non-differential symbols are external references

**Planned Features (FIXMEs):**
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

**Implementation Status:**
- Parser is fully implemented with recursive expression traversal
- Supports automatic variable classification based on conventions
- Test suite includes multiple model examples (Ramsey, Solow, etc.)
- Documentation includes theory, examples, and API reference

**Path Configuration:**
- All scripts should be run from the project root directory with `julia --project=.`
- The package uses standard Julia package structure with Project.toml
- No manual LOAD_PATH modifications needed when using --project flag

**Development Notes:**
- The `scrap/` folder is for temporary test files and code snippets
- Files in `scrap/` are not tracked in version control and should not be persisted
- Use `scrap/` for experimental code, copy-pasted documentation, or test scripts
- May contain relevant code examples or documentation for reference

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