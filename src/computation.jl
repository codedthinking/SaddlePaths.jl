"""
NumericalSolution stores all the necessary components to iterate towards a solution and evaluate it.

The user would typically interact with a NumericalSolution through special getter and setter functions.
"""
mutable struct NumericalSolution
    model::SaddlePath
    params::Vector{T} where T <: Number
    steady_state::Vector{T} where T <: Number
    grid::SpectralKit.AbstractGrid
    basis::SpectralKit.FunctionBasis
    coefs::Vector{T} where T <: Number
    error::AbstractFloat
    iterations::Int
    converged::Bool
end

"""
Functions to get the Chebyshev approximation and its derivative.
"""
costate(s::NumericalSolution) = linear_combination(s.basis, s.coefs)
𝑑costate(s::NumericalSolution) = x -> linear_combination(s.basis, s.coefs, 𝑑(x))[1]

"""
Evaluate the c(k) expression using the approximation.
"""
function approx_eval(expr::Expr, var::Symbol, sol::NumericalSolution)
end

state_dot(s::NumericalSolution)
costate_dot(s::NumericalSolution)
