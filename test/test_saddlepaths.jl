#!/usr/bin/env julia

using SaddlePaths

# Test 1: Ramsey model with costate as function
@model ramsey_functional begin
    𝑑c = c(k) * (α * k^(α - 1) - δ - ρ)
    𝑑k = k^α - δ * k - c(k)
end

println("=== Ramsey Model (Functional Form) ===")
println("  state = ", ramsey_functional.state)
println("  costate = ", ramsey_functional.costate, " (c appears as c(k))")
println("  params = ", ramsey_functional.params, " (Greek letters)")
println("  ODEs = ", ramsey_functional.ODEs)

# Test 1b: Original Ramsey model (c as plain variable)
@model ramsey begin
    𝑑c = c * (α * k^(α - 1) - δ - ρ)
    𝑑k = k^α - δ * k - c
end

println("\n=== Ramsey Model (Standard Form) ===")
println("  state = ", ramsey.state, " (both are states since no function calls)")
println("  costate = ", ramsey.costate)
println("  params = ", ramsey.params, " (Greek letters)")
println("  ODEs = ", ramsey.ODEs)

# Test 2: Solow model (no costates)
@model solow begin
    𝑑k = s * k^α - δ * k
end

println("\n=== Solow Model ===")
println("  state = ", solow.state)
println("  costate = ", solow.costate)
println("  params = ", solow.params, " (s is not Greek, α and δ are)")
println("  ODEs = ", solow.ODEs)

# Test 3: Two-state model
@model twostates begin
    𝑑k₁ = α * k₁ - δ * k₁
    𝑑k₂ = β * k₂ - γ * k₂
end

println("\n=== Two-State Model ===")
println("  state = ", twostates.state)
println("  costate = ", twostates.costate)
println("  params = ", twostates.params, " (all Greek)")
println("  ODEs = ", twostates.ODEs)

# Test 4: Model with costate as function call
@model optimal_control begin
    𝑑x = u(x) + α * x
    𝑑u = β * u - γ * x
end

println("\n=== Optimal Control Model ===")
println("  state = ", optimal_control.state, " (x is state)")
println("  costate = ", optimal_control.costate, " (u appears as u(x), so it's costate)")
println("  params = ", optimal_control.params, " (Greek letters)")
println("  ODEs = ", optimal_control.ODEs)

# Test 5: Model with non-Greek parameter (should be treated as global/function)
@model mixed begin
    𝑑k = f(k) - δ * k
end

println("\n=== Mixed Model ===")
println("  state = ", mixed.state)
println("  costate = ", mixed.costate)
println("  params = ", mixed.params, " (only δ is Greek)")
println("  ODEs = ", mixed.ODEs)
println("  Note: 'f' is treated as a global function, not a parameter")