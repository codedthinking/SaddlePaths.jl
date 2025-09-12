# SaddlePaths.jl

[![Build Status](https://github.com/koren/SaddlePaths.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/koren/SaddlePaths.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/koren/SaddlePaths.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/koren/SaddlePaths.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://koren.github.io/SaddlePaths.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://koren.github.io/SaddlePaths.jl/dev/)

A Julia package for solving saddle-path problems in continuous-time macroeconomic models using stable manifold methods.

## Features

- ✅ **Ergonomic DSL** for declaring continuous-time macro models with Unicode notation
- ✅ **Automatic Variable Classification** - States, costates, and parameters detected automatically
- ✅ **Isolated Namespaces** - Each solution instance has its own variable namespace
- ✅ **Modular Backend System** - Extensible architecture for different numerical methods
- ✅ **Comprehensive Testing** - 242 parser tests passing, TDD approach throughout
- 🔄 **Backend Abstraction** - Ready for spectral methods, finite differences, etc.
- 🔄 **Test-Driven Development** - Extensive test suite with 60+ computation tests written

## Installation

```julia
using Pkg
Pkg.add("SaddlePaths")
```

## Quick Example

```julia
using SaddlePaths

# Define a Ramsey growth model with automatic variable classification
@model ramsey begin
    # c(k) notation automatically marks c as costate variable
    # Greek letters (α, δ, ρ) automatically detected as parameters
    # k is automatically classified as state variable
    𝑑c = c(k) * (α * k^(α - 1) - δ - ρ)
    𝑑k = k^α - δ * k - c(k)
end

# Model structure is automatically created:
# - state = [:k]
# - costate = [:c]
# - params = [:α, :δ, :ρ]

# Create numerical solution with isolated namespace
solution = NumericalSolution(ramsey, [0.3, 0.1, 0.05])  # α, δ, ρ

# Evaluate ODEs at specific points
state_vals = [1.5]    # k = 1.5
costate_vals = [0.8]  # c = 0.8
dk, dc = evaluate_all_odes(solution, state_vals, costate_vals)

# Clean up namespace when done
cleanup!(solution)
```

## Current Implementation Status

- ✅ **Parser**: Fully implemented with 242 tests passing
- ✅ **Namespace Management**: Complete isolation system for multiple solutions
- ✅ **Backend Abstraction**: Modular system for numerical methods
- 🔄 **Numerical Solution**: Partial implementation with TDD tests
- 🔄 **Spectral Methods**: Backend ready, awaiting library integration

## Architecture

SaddlePaths.jl follows a clean, modular architecture:

- **Parser Module**: DSL for model specification with automatic variable classification
- **Namespace Module**: Isolated variable scopes for multiple solution instances
- **Backend Module**: Abstract interface for different numerical methods
- **Computation Module**: Numerical solution framework with ODE evaluation

## Testing

The project uses comprehensive test-driven development:

```bash
# Run all tests
make test

# Run specific test suites
julia --project=. test/test_parser.jl     # 242 tests passing
julia --project=. test/test_namespace.jl  # Namespace isolation
julia --project=. test/test_computation.jl # Numerical solution (TDD)
```

## Development Status

- **Lines of Code**: ~2,900 total
- **Test Coverage**: Parser fully tested, computation tests written
- **Dependencies**: None currently (ready for spectral method libraries)
- **Julia Version**: 1.10+

## Contributing

This project uses a test-driven development approach. New features should:
1. Add comprehensive tests first
2. Implement the functionality
3. Ensure all tests pass

## Documentation

For detailed documentation, see [the docs](https://koren.github.io/SaddlePaths.jl/dev/).

## License

MIT