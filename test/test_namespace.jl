using Test
using SaddlePaths

# Import namespace functions for testing
import SaddlePaths: create_namespace, delete_namespace, set_variable!, get_variable,
                    update_namespace!, eval_in_namespace, with_namespace, Variables

@testset "Namespace Management" begin

    @testset "create_namespace" begin
        # Test namespace creation
        ns_id = create_namespace()
        @test ns_id isa Symbol
        @test startswith(string(ns_id), "ns_")
        @test length(string(ns_id)) > 3  # Should be longer than just "ns_"

        # Test that namespace exists
        @test hasfield(Variables, ns_id)
        ns = getfield(Variables, ns_id)
        @test ns isa Module

        # Clean up
        delete_namespace(ns_id)
    end

    @testset "set_variable! and get_variable" begin
        ns_id = create_namespace()

        # Test setting and getting various types of values
        set_variable!(ns_id, :k, 1.5)
        @test get_variable(ns_id, :k) == 1.5

        set_variable!(ns_id, :α, 0.3)
        @test get_variable(ns_id, :α) == 0.3

        set_variable!(ns_id, :name, "test")
        @test get_variable(ns_id, :name) == "test"

        set_variable!(ns_id, :arr, [1, 2, 3])
        @test get_variable(ns_id, :arr) == [1, 2, 3]

        # Test that variables are isolated between namespaces
        ns_id2 = create_namespace()
        set_variable!(ns_id2, :k, 2.0)
        @test get_variable(ns_id, :k) == 1.5  # Should still be 1.5
        @test get_variable(ns_id2, :k) == 2.0  # Should be 2.0

        # Clean up
        delete_namespace(ns_id)
        delete_namespace(ns_id2)
    end

    @testset "delete_namespace" begin
        ns_id = create_namespace()
        set_variable!(ns_id, :test_var, 42)

        # Verify namespace exists and has variable
        @test hasfield(Variables, ns_id)
        @test get_variable(ns_id, :test_var) == 42

        # Delete namespace
        delete_namespace(ns_id)

        # Namespace should still exist but variables should be cleared
        @test hasfield(Variables, ns_id)
        ns = getfield(Variables, ns_id)
        @test ns isa Module

        # Variable should be gone (or set to nothing)
        @test_throws UndefVarError get_variable(ns_id, :test_var)
    end

    @testset "eval_in_namespace" begin
        ns_id = create_namespace()

        # Set up variables
        set_variable!(ns_id, :x, 2.0)
        set_variable!(ns_id, :y, 3.0)

        # Test simple expression evaluation
        result = eval_in_namespace(ns_id, :(x + y))
        @test result == 5.0

        result = eval_in_namespace(ns_id, :(x * y))
        @test result == 6.0

        result = eval_in_namespace(ns_id, :(x^2 + y^2))
        @test result == 4.0 + 9.0

        # Test with function calls
        set_variable!(ns_id, :z, 1.0)
        result = eval_in_namespace(ns_id, :(sin(z) + cos(z)))
        @test result ≈ sin(1.0) + cos(1.0)

        # Test that expressions can modify variables
        eval_in_namespace(ns_id, :(x = x + 1))
        @test get_variable(ns_id, :x) == 3.0

        # Clean up
        delete_namespace(ns_id)
    end

    @testset "update_namespace!" begin
        # Create a mock model structure
        struct MockModel
            state::Vector{Symbol}
            costate::Vector{Symbol}
            params::Vector{Symbol}
        end

        model = MockModel([:k, :h], [:λ_k, :λ_h], [:α, :β, :δ])

        ns_id = create_namespace()

        # Update namespace with values
        state_vals = [1.5, 2.0]
        costate_vals = [0.8, 0.6]
        param_vals = [0.3, 0.2, 0.1]

        update_namespace!(ns_id, model, state_vals, costate_vals, param_vals)

        # Verify all variables were set correctly
        @test get_variable(ns_id, :k) == 1.5
        @test get_variable(ns_id, :h) == 2.0
        @test get_variable(ns_id, :λ_k) == 0.8
        @test get_variable(ns_id, :λ_h) == 0.6
        @test get_variable(ns_id, :α) == 0.3
        @test get_variable(ns_id, :β) == 0.2
        @test get_variable(ns_id, :δ) == 0.1

        # Clean up
        delete_namespace(ns_id)
    end

    @testset "with_namespace" begin
        ns_id = create_namespace()

        # Set up variables
        set_variable!(ns_id, :a, 10)
        set_variable!(ns_id, :b, 20)

        # Test function execution in namespace context
        result = with_namespace(ns_id) do ns
            a_val = get_variable(ns, :a)
            b_val = get_variable(ns, :b)
            return a_val + b_val
        end
        @test result == 30

        # Test more complex operations
        result = with_namespace(ns_id) do ns
            set_variable!(ns, :c, 5)
            return get_variable(ns, :a) * get_variable(ns, :c)
        end
        @test result == 50
        @test get_variable(ns_id, :c) == 5  # Verify variable was set

        # Clean up
        delete_namespace(ns_id)
    end

    @testset "Namespace isolation" begin
        # Test that multiple namespaces are truly isolated
        ns1 = create_namespace()
        ns2 = create_namespace()

        # Set same variable names with different values
        set_variable!(ns1, :x, 100)
        set_variable!(ns2, :x, 200)

        set_variable!(ns1, :y, "hello")
        set_variable!(ns2, :y, "world")

        # Verify isolation
        @test get_variable(ns1, :x) == 100
        @test get_variable(ns2, :x) == 200
        @test get_variable(ns1, :y) == "hello"
        @test get_variable(ns2, :y) == "world"

        # Test expression evaluation in each namespace
        result1 = eval_in_namespace(ns1, :(x + 10))
        result2 = eval_in_namespace(ns2, :(x + 10))
        @test result1 == 110
        @test result2 == 210

        # Clean up
        delete_namespace(ns1)
        delete_namespace(ns2)
    end

    @testset "Error handling" begin
        # Test accessing non-existent namespace
        @test_throws UndefVarError getfield(Variables, :nonexistent_ns)

        # Test accessing non-existent variable
        ns_id = create_namespace()
        @test_throws UndefVarError get_variable(ns_id, :nonexistent_var)

        # Test deleting non-existent namespace (should not error)
        delete_namespace(:nonexistent_ns)  # Should handle gracefully

        # Clean up
        delete_namespace(ns_id)
    end

    @testset "Complex expressions" begin
        ns_id = create_namespace()

        # Set up variables for complex expression
        set_variable!(ns_id, :k, 2.0)
        set_variable!(ns_id, :α, 0.5)
        set_variable!(ns_id, :δ, 0.1)
        set_variable!(ns_id, :c, 1.0)

        # Test complex economic expression (similar to Ramsey model)
        expr = :(k^α - δ * k - c)
        result = eval_in_namespace(ns_id, expr)
        expected = 2.0^0.5 - 0.1 * 2.0 - 1.0
        @test result ≈ expected

        # Test with vectorized operations
        set_variable!(ns_id, :k_vec, [1.0, 2.0, 3.0])
        result = eval_in_namespace(ns_id, :(k_vec .^ α))
        expected = [1.0^0.5, 2.0^0.5, 3.0^0.5]
        @test result ≈ expected

        # Clean up
        delete_namespace(ns_id)
    end

    @testset "Memory management" begin
        # Create multiple namespaces to test memory handling
        namespaces = Symbol[]
        for i in 1:10
            ns_id = create_namespace()
            push!(namespaces, ns_id)
            set_variable!(ns_id, :test_var, i)
        end

        # Verify all namespaces work independently
        for (i, ns_id) in enumerate(namespaces)
            @test get_variable(ns_id, :test_var) == i
        end

        # Clean up all namespaces
        for ns_id in namespaces
            delete_namespace(ns_id)
        end

        # Verify cleanup worked
        for ns_id in namespaces
            @test hasfield(Variables, ns_id)
            # Variables should be cleared
            @test_throws UndefVarError get_variable(ns_id, :test_var)
        end
    end

    @testset "Integration with SaddlePath" begin
        # Test with actual SaddlePath model
        model = @model ramsey begin
            𝑑c = c(k) * (α * k^(α - 1) - δ - ρ)
            𝑑k = k^α - δ * k - c(k)
        end

        ns_id = create_namespace()

        # Update namespace with model parameters
        param_vals = [0.3, 0.1, 0.05]  # α, δ, ρ
        for (name, val) in zip(model.params, param_vals)
            set_variable!(ns_id, name, val)
        end

        # Test parameter access
        @test get_variable(ns_id, :α) == 0.3
        @test get_variable(ns_id, :δ) == 0.1
        @test get_variable(ns_id, :ρ) == 0.05

        # Test expression evaluation with parameters
        result = eval_in_namespace(ns_id, :(α + δ + ρ))
        @test result == 0.3 + 0.1 + 0.05

        # Clean up
        delete_namespace(ns_id)
    end
end