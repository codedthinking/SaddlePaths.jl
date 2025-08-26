using SaddlePaths
using Test

@testset "SaddlePaths.jl" begin
    
    @testset "SaddlePath struct" begin
        @testset "Basic construction" begin
            # Create a simple saddle path problem
            sp = SaddlePath(
                states = [:k, :h],
                costates = [:λ_k, :λ_h],
                parameters = [:α, :β, :δ, :ρ]
            )
            
            @test sp isa SaddlePath
            @test length(states(sp)) == 2
            @test length(costates(sp)) == 2
            @test length(parameters(sp)) == 4
        end
        
        @testset "State accessors" begin
            sp = SaddlePath(
                states = [:k, :h],
                costates = [:λ_k, :λ_h],
                parameters = [:α, :β, :δ]
            )
            
            @test states(sp) == [:k, :h]
            @test :k in states(sp)
            @test :h in states(sp)
            @test :c ∉ states(sp)
        end
        
        @testset "Costate accessors" begin
            sp = SaddlePath(
                states = [:k],
                costates = [:λ_k, :λ_h],
                parameters = [:α]
            )
            
            @test costates(sp) == [:λ_k, :λ_h]
            @test :λ_k in costates(sp)
            @test :λ_h in costates(sp)
            @test :λ_c ∉ costates(sp)
        end
        
        @testset "Parameter accessors" begin
            sp = SaddlePath(
                states = [:k],
                costates = [:λ],
                parameters = [:α, :β, :δ, :ρ, :σ]
            )
            
            @test parameters(sp) == [:α, :β, :δ, :ρ, :σ]
            @test :α in parameters(sp)
            @test :σ in parameters(sp)
            @test :γ ∉ parameters(sp)
        end
        
        @testset "Dimension queries" begin
            sp = SaddlePath(
                states = [:k, :h, :a],
                costates = [:λ_k, :λ_h],
                parameters = [:α, :β]
            )
            
            @test dim_states(sp) == 3
            @test dim_costates(sp) == 2
            @test dim_parameters(sp) == 2
            @test dim_total(sp) == 5  # states + costates
        end
    end
    
    @testset "SaddlePath with dynamics" begin
        @testset "Construction with dynamics functions" begin
            # Define dynamics
            f_k = (k, c, θ) -> θ.α * k - c
            f_c = (k, c, θ) -> c * (θ.ρ - θ.α)
            
            sp = SaddlePath(
                states = [:k],
                costates = [:c],
                parameters = [:α, :ρ],
                dynamics = Dict(:k => f_k, :c => f_c)
            )
            
            @test haskey(dynamics(sp), :k)
            @test haskey(dynamics(sp), :c)
            @test length(dynamics(sp)) == 2
        end
        
        @testset "Evaluate dynamics" begin
            sp = SaddlePath(
                states = [:k],
                costates = [:c],
                parameters = [:α, :δ],
                dynamics = Dict(
                    :k => (k, c, θ) -> k^θ.α - c - θ.δ*k,
                    :c => (k, c, θ) -> c * (θ.α * k^(θ.α-1) - θ.δ)
                )
            )
            
            # Test evaluation at a point
            θ = (α=0.3, δ=0.1)
            k_val = 1.0
            c_val = 0.5
            
            @test evaluate_dynamics(sp, :k, k_val, c_val, θ) ≈ k_val^0.3 - c_val - 0.1*k_val
            @test evaluate_dynamics(sp, :c, k_val, c_val, θ) ≈ c_val * (0.3 * k_val^(-0.7) - 0.1)
        end
    end
    
    @testset "SaddlePath with steady state" begin
        @testset "Steady state equations" begin
            sp = SaddlePath(
                states = [:k],
                costates = [:c],
                parameters = [:α, :δ, :ρ],
                steady_state_equations = [
                    :(k_ss^α - c_ss - δ*k_ss),
                    :(ρ + δ - α*k_ss^(α-1))
                ]
            )
            
            @test length(steady_state_equations(sp)) == 2
            @test steady_state_variables(sp) == [:k_ss, :c_ss]
        end
        
        @testset "Steady state values" begin
            sp = SaddlePath(
                states = [:k],
                costates = [:c],
                parameters = [:α, :δ, :ρ],
                steady_state = Dict(:k => 2.0, :c => 0.5)
            )
            
            @test steady_state_value(sp, :k) == 2.0
            @test steady_state_value(sp, :c) == 0.5
            @test_throws KeyError steady_state_value(sp, :h)
        end
        
        @testset "Set steady state" begin
            sp = SaddlePath(
                states = [:k],
                costates = [:c],
                parameters = [:α]
            )
            
            set_steady_state!(sp, :k, 3.0)
            set_steady_state!(sp, :c, 1.0)
            
            @test steady_state_value(sp, :k) == 3.0
            @test steady_state_value(sp, :c) == 1.0
        end
    end
    
    @testset "Variable validation" begin
        @testset "Invalid variable names" begin
            # State names can't start with Greek letters
            @test_throws ArgumentError SaddlePath(
                states = [:α_k],
                costates = [:c],
                parameters = [:α]
            )
            
            # Parameter names must start with Greek letters
            @test_throws ArgumentError SaddlePath(
                states = [:k],
                costates = [:c],
                parameters = [:alpha]
            )
            
            # No duplicate names across categories
            @test_throws ArgumentError SaddlePath(
                states = [:k],
                costates = [:k],  # Duplicate!
                parameters = [:α]
            )
        end
        
        @testset "Valid variable names" begin
            # These should all work
            sp = SaddlePath(
                states = [:k, :k_1, :capital],
                costates = [:λ, :λ_1, :mu],
                parameters = [:α, :β_1, :δ_k]
            )
            
            @test states(sp) == [:k, :k_1, :capital]
            @test costates(sp) == [:λ, :λ_1, :mu]
            @test parameters(sp) == [:α, :β_1, :δ_k]
        end
    end
    
    @testset "Jacobian structure" begin
        @testset "Jacobian at steady state" begin
            sp = SaddlePath(
                states = [:k, :h],
                costates = [:λ_k, :λ_h],
                parameters = [:α, :δ],
                steady_state = Dict(:k => 1.0, :h => 2.0, :λ_k => 0.5, :λ_h => 0.3)
            )
            
            # Mock Jacobian for testing structure
            J = zeros(4, 4)
            J[1,1] = -0.1  # ∂ḱ/∂k
            J[1,3] = -1.0  # ∂ḱ/∂λ_k
            J[3,1] = 0.2   # ∂λ̇_k/∂k
            
            set_jacobian!(sp, J)
            
            @test size(jacobian(sp)) == (4, 4)
            @test jacobian(sp)[1,1] == -0.1
            
            # Test eigenvalue computation
            eigs = eigenvalues(sp)
            @test length(eigs) == 4
            
            # Test stability classification
            @test n_stable(sp) + n_unstable(sp) == 4
        end
    end
    
    @testset "Policy functions" begin
        @testset "Policy function storage" begin
            sp = SaddlePath(
                states = [:k],
                costates = [:c],
                parameters = [:α]
            )
            
            # Mock policy function
            policy = k -> 0.5 * k^0.6
            set_policy!(sp, :c, policy)
            
            @test has_policy(sp, :c) == true
            @test has_policy(sp, :k) == false
            
            # Evaluate policy
            @test policy_value(sp, :c, 2.0) ≈ 0.5 * 2.0^0.6
        end
        
        @testset "Policy gradient" begin
            sp = SaddlePath(
                states = [:k],
                costates = [:c],
                parameters = [:α]
            )
            
            # Policy with known gradient
            policy = k -> k^2
            policy_grad = k -> 2*k
            
            set_policy!(sp, :c, policy, gradient=policy_grad)
            
            @test policy_gradient(sp, :c, 3.0) ≈ 6.0
        end
    end
    
    @testset "Model compilation state" begin
        @testset "Compilation flags" begin
            sp = SaddlePath(
                states = [:k],
                costates = [:c],
                parameters = [:α]
            )
            
            @test is_compiled(sp) == false
            @test has_steady_state(sp) == false
            @test has_dynamics(sp) == false
            
            # After adding dynamics
            set_dynamics!(sp, :k, (k,c,θ) -> θ.α*k - c)
            set_dynamics!(sp, :c, (k,c,θ) -> c*θ.α)
            
            @test has_dynamics(sp) == true
            
            # After setting steady state
            set_steady_state!(sp, :k, 1.0)
            set_steady_state!(sp, :c, 0.3)
            
            @test has_steady_state(sp) == true
        end
    end
    
    @testset "Pretty printing" begin
        sp = SaddlePath(
            states = [:k, :h],
            costates = [:λ_k, :λ_h],
            parameters = [:α, :β, :δ]
        )
        
        # Test string representation
        str = string(sp)
        @test occursin("SaddlePath", str)
        @test occursin("states: 2", str)
        @test occursin("costates: 2", str)
        @test occursin("parameters: 3", str)
        
        # Test detailed summary
        summary_str = summary(sp)
        @test occursin("k", summary_str)
        @test occursin("λ_k", summary_str)
        @test occursin("α", summary_str)
    end
end