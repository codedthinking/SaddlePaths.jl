# SaddlePaths.jl Documentation

```@meta
CurrentModule = SaddlePaths
```

## Overview

SaddlePaths.jl is a Julia package for defining and solving saddle path problems in dynamic optimization. It provides a domain-specific language (DSL) for specifying systems of differential equations with automatic classification of variables into states, costates, and parameters.

## Features

- **Intuitive DSL**: Define models using natural mathematical notation with Unicode characters
- **Automatic Classification**: Intelligently identifies state variables, costate/control variables, and parameters
- **Flexible Notation**: Support for both standard and functional notation for costates
- **Greek Letter Convention**: Automatic parameter detection based on Greek letter naming
- **Expression Preservation**: Maintains Julia AST for further symbolic or numeric processing

## Quick Example

```julia
using SaddlePaths

# Define a Ramsey growth model
@model ramsey begin
    𝑑c = c(k) * (α * k^(α - 1) - δ - ρ)
    𝑑k = k^α - δ * k - c(k)
end

# Access model components
println("States: ", ramsey.state)      # [:k]
println("Costates: ", ramsey.costate)  # [:c]
println("Parameters: ", ramsey.params) # [:α, :δ, :ρ]
```

## Installation

Since this is part of the babyboom project, add the source to your load path:

```julia
push!(LOAD_PATH, "path/to/babyboom/src")
using SaddlePaths
```

For standalone use, you would typically install via:

```julia
using Pkg
Pkg.add("SaddlePaths")  # Not yet registered
```

## Contents

```@contents
Pages = ["getting_started.md", "dsl.md", "examples.md", "api.md", "theory.md"]
Depth = 2
```

## Index

```@index
```