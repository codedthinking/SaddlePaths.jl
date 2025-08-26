# SaddlePaths.jl Development Guide

## Project Overview

SaddlePaths.jl is a Julia package for solving saddle-path problems in continuous-time macroeconomic models using stable manifold methods. The package provides:

1. An ergonomic DSL for model specification using Unicode notation
2. Automatic steady state solving from symbolic equations
3. Explicit co-state variable support via `@costate` macro
4. Symbolic differentiation and code generation
5. Policy function approximation via Chebyshev/Smolyak methods
6. Simulation capabilities

## Package Structure

```
SaddlePaths.jl/
├── src/
│   └── SaddlePaths.jl          # Main module file (exports, API)
├── test/
│   └── runtests.jl            # Test suite
├── docs/
│   ├── make.jl                # Documentation builder
│   └── src/
│       └── index.md           # Documentation homepage
├── Project.toml               # Package metadata and dependencies
├── Makefile                   # Build automation
└── README.md                  # User-facing documentation
```

## Implementation Plan

The package will be implemented incrementally with the following modules:

1. **dsl.jl**: `@model` and `@costate` macros for DSL parsing
2. **parser.jl**: Unicode parsing, symbol classification, equation type detection
3. **symbolic_backend.jl**: Symbolics.jl/ModelingToolkit.jl integration
4. **steady_state.jl**: Automatic steady state solver for `0 = ...` equations
5. **stability.jl**: Eigenvalue analysis and stability checking
6. **bases/** : Chebyshev and Smolyak basis functions
7. **policy.jl**: Stable manifold solver via collocation
8. **simulate.jl**: ODE integration using computed policies

## DSL Design

- State variables: Declared with `𝒹` prefix (e.g., `𝒹k = ...`)
- Co-state variables: Declared with `@costate 𝒹c = ...` 
- Parameters: First character is Greek (e.g., `α`, `β`, `δ`)
- Steady state equations: Use `0 = ...` to have system solve for steady state
- Auxiliary variables: Latin identifiers in algebraic equations

## Key Algorithms

1. **Stable Manifold Method**: Solve for policy function c(k) satisfying:
   ```
   ∇c(k) · ḱ(θ,k,c(k)) = ċ(θ,k,c(k))
   c(k_ss) = c_ss
   ```

2. **Chebyshev Collocation**: Approximate c(k) as basis expansion
3. **Smolyak Sparse Grids**: Handle high-dimensional state spaces

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

## Dependencies (to be added)

- Symbolics.jl / ModelingToolkit.jl
- DifferentialEquations.jl
- (Spectral approximation packages TBD)