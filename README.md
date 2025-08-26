# SaddlePaths.jl

[![Build Status](https://github.com/koren/SaddlePaths.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/koren/SaddlePaths.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/koren/SaddlePaths.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/koren/SaddlePaths.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://koren.github.io/SaddlePaths.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://koren.github.io/SaddlePaths.jl/dev/)

A Julia package for solving saddle-path problems in continuous-time macroeconomic models using stable manifold methods.

## Features

- Ergonomic DSL for declaring continuous-time macro models with Unicode notation
- Symbolic and numeric derivatives via Symbolics.jl/ModelingToolkit.jl
- Steady state analysis and stability checking
- Stable manifold policy function approximation using Chebyshev polynomials and Smolyak sparse grids
- Simulation of deterministic paths

## Installation

```julia
using Pkg
Pkg.add("SaddlePaths")
```

## Quick Example

```julia
using SaddlePaths

# Define a simple model
@model begin
  𝒹k = α*k - c
  𝒹c = β*(α*k - c) - δ*c
end

# Specify steady states
k_ss(α,β,δ) = α/δ
c_ss(α,β,δ) = α*k_ss(α,β,δ)

# Compile and solve
M = compile_model(@locals)
A = analyze(M; θ=(α=0.3, β=0.99, δ=0.1))
π = solve_policy(M; θ=(α=0.3, β=0.99, δ=0.1))

# Simulate
traj = simulate(M, π; θ=(α=0.3, β=0.99, δ=0.1), k0=0.5, T=100.0)
```

## Documentation

For detailed documentation, see [the docs](https://koren.github.io/SaddlePaths.jl/dev/).

## License

MIT