# Abstract backend type
abstract type AbstractBackend end

# Spectral backend implementation
struct SpectralBackend <: AbstractBackend
    basis::Any  # Placeholder for spectral basis
    grid::AbstractVector  # Evaluation grid
    coeffs::Vector{Float64}  # Coefficients for approximation
end

# Finite difference backend implementation
struct FiniteDifferenceBackend <: AbstractBackend
    grid::AbstractVector  # Grid points
    coeffs::Vector{Float64}  # Coefficients for approximation
end

# Projection backend implementation
struct ProjectionBackend <: AbstractBackend
    basis::Any  # Basis functions
    grid::AbstractVector  # Evaluation grid
    coeffs::Vector{Float64}  # Coefficients for approximation
end

# Core backend API functions
grid(backend::AbstractBackend) = error(
    "grid() not implemented for $(typeof(backend)). Implement this method in your concrete backend type.",
)

collocation_matrix(backend::AbstractBackend) = error(
    "collocation_matrix() not implemented for $(typeof(backend)). Implement this method in your concrete backend type.",
)

costate(backend::AbstractBackend, coeffs::AbstractVector) = error(
    "costate() not implemented for $(typeof(backend)). Implement this method in your concrete backend type.",
)

∇costate(backend::AbstractBackend, coeffs::AbstractVector) = error(
    "∇costate() not implemented for $(typeof(backend)). Implement this method in your concrete backend type.",
)
