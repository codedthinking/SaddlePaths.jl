```@meta
CurrentModule = SaddlePaths
```

# SaddlePaths.jl

A Julia package for solving saddle-path problems in continuous-time macroeconomic models using stable manifold methods.

## Overview

SaddlePaths.jl provides an ergonomic domain-specific language (DSL) for declaring continuous-time macroeconomic models with Unicode notation, along with efficient numerical methods for computing stable manifold policy functions via Chebyshev polynomials and Smolyak sparse grids.

## Features

- **Ergonomic DSL**: Write models naturally using Unicode notation (e.g., `𝒹k` for time derivatives, Greek letters for parameters)
- **Symbolic Core**: Automatic generation of Jacobians and Hessians using Symbolics.jl/ModelingToolkit.jl
- **Steady State Analysis**: Built-in tools for analyzing steady states and checking stability
- **Stable Manifold Solver**: Approximation of policy functions `c(k)` using spectral methods
- **Simulation**: Integration of deterministic paths using DifferentialEquations.jl

## Installation

```julia
using Pkg
Pkg.add("SaddlePaths")
```

## Quick Start

```julia
using SaddlePaths

# Define a model using the DSL
@model begin
  𝒹k = α*k - c
  𝒹c = β*(α*k - c) - δ*c  
end

# Specify steady state functions
k_ss(α,β,δ) = α/δ
c_ss(α,β,δ) = α*k_ss(α,β,δ)

# Compile and analyze
M = compile_model(@locals)
A = analyze(M; θ=(α=0.3, β=0.99, δ=0.1))

# Solve for policy function
π = solve_policy(M; θ=(α=0.3, β=0.99, δ=0.1), 
                 domain=:auto, order=7, smolyak_level=3)

# Simulate trajectory
traj = simulate(M, π; θ=(α=0.3, β=0.99, δ=0.1), 
                k0=0.5, T=100.0)
```

## DSL Syntax

The package uses Unicode characters to distinguish between different types of variables:

- **State variables**: Declared using `𝒹` prefix (e.g., `𝒹k = ...`)
- **Parameters**: First character is Greek (e.g., `α`, `β`, `δ`)
- **Control variables**: Latin identifiers not declared as states

## API Reference

```@index
```

```@autodocs
Modules = [SaddlePaths]
```