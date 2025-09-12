using Test
using SaddlePaths

# Import internal functions for testing
import SaddlePaths: is_operator, is_dotted_operator, isfunctioncall, is_greek_letter,
                    collect_all_symbols!, find_function_calls!, parse_differential_equations,
                    vectorize_function_calls

@testset "Parser Functions" begin
    
    @testset "is_operator" begin
        # Test with operators
        @test is_operator(:+) == true
        @test is_operator(:-) == true
        @test is_operator(:*) == true
        @test is_operator(:/) == true
        @test is_operator(:^) == true
        @test is_operator(:%) == true
        @test is_operator(:÷) == true
        @test is_operator(:(==)) == true
        @test is_operator(:(!=)) == true
        @test is_operator(:(<=)) == true
        @test is_operator(:(>=)) == true
        @test is_operator(:(<)) == true
        @test is_operator(:(>)) == true
        @test is_operator(:(&&)) == true
        @test is_operator(:(||)) == true
        
        # Test dot operator exclusion
        @test is_operator(:.) == false
        
        # Test non-operators
        @test is_operator(:foo) == false
        @test is_operator(:α) == false
        @test is_operator(:log) == false
        
        # Test non-symbol inputs
        @test is_operator("not a symbol") == false
        @test is_operator(123) == false
        @test is_operator(nothing) == false
        @test is_operator(:(x + y)) == false
    end
    
    @testset "is_dotted_operator" begin
        # Test dotted operators
        @test is_dotted_operator(:.+) == true
        @test is_dotted_operator(:.-) == true
        @test is_dotted_operator(:.*) == true
        @test is_dotted_operator(:./) == true
        @test is_dotted_operator(:.^) == true
        @test is_dotted_operator(:(.==)) == true
        @test is_dotted_operator(:(.!=)) == true
        @test is_dotted_operator(:(.<=)) == true
        @test is_dotted_operator(:(.>=)) == true
        @test is_dotted_operator(:(.<)) == true
        @test is_dotted_operator(:(.>)) == true
        
        # Test non-dotted operators
        @test is_dotted_operator(:+) == false
        @test is_dotted_operator(:-) == false
        @test is_dotted_operator(:*) == false
        
        # Test non-operators
        @test is_dotted_operator(:.) == false
        @test is_dotted_operator(Symbol(".foo")) == false
        @test is_dotted_operator(:foo) == false
        
        # Test non-symbol inputs
        @test is_dotted_operator("not a symbol") == false
        @test is_dotted_operator(123) == false
        @test is_dotted_operator(nothing) == false
    end
    
    @testset "isfunctioncall" begin
        # Test regular function calls
        @test isfunctioncall(:(f(x))) == true
        @test isfunctioncall(:(log(x))) == true
        @test isfunctioncall(:(c(k))) == true
        @test isfunctioncall(:(foo(x, y, z))) == true
        
        # Test operator calls
        @test isfunctioncall(:(x + y)) == true
        @test isfunctioncall(:(x - y)) == true
        @test isfunctioncall(:(x * y)) == true
        @test isfunctioncall(:(x / y)) == true
        @test isfunctioncall(:(x^2)) == true
        
        # Test broadcasting calls
        @test isfunctioncall(:(f.(x))) == true
        @test isfunctioncall(:(log.(x))) == true
        
        # Test syntactic operators
        @test isfunctioncall(:(x && y)) == true
        @test isfunctioncall(:(x || y)) == true
        
        # Test operator symbols
        @test isfunctioncall(:+) == true
        @test isfunctioncall(:-) == true
        @test isfunctioncall(:*) == true
        @test isfunctioncall(:/) == true
        
        # Test non-function calls
        @test isfunctioncall(:x) == false
        @test isfunctioncall(:α) == false
        @test isfunctioncall(:(x)) == false  # Just a symbol in parens
        @test isfunctioncall(123) == false
        @test isfunctioncall("string") == false
        @test isfunctioncall(nothing) == false
    end
    
    @testset "is_greek_letter" begin
        # Test lowercase Greek letters
        @test is_greek_letter(:α) == true
        @test is_greek_letter(:β) == true
        @test is_greek_letter(:γ) == true
        @test is_greek_letter(:δ) == true
        @test is_greek_letter(:ε) == true
        @test is_greek_letter(:ζ) == true
        @test is_greek_letter(:η) == true
        @test is_greek_letter(:θ) == true
        @test is_greek_letter(:ι) == true
        @test is_greek_letter(:κ) == true
        @test is_greek_letter(:λ) == true
        @test is_greek_letter(:μ) == true
        @test is_greek_letter(:ν) == true
        @test is_greek_letter(:ξ) == true
        @test is_greek_letter(:ο) == true
        @test is_greek_letter(:π) == true
        @test is_greek_letter(:ρ) == true
        @test is_greek_letter(:σ) == true
        @test is_greek_letter(:τ) == true
        @test is_greek_letter(:υ) == true
        @test is_greek_letter(:φ) == true
        @test is_greek_letter(:χ) == true
        @test is_greek_letter(:ψ) == true
        @test is_greek_letter(:ω) == true
        
        # Test uppercase Greek letters
        @test is_greek_letter(:Α) == true
        @test is_greek_letter(:Β) == true
        @test is_greek_letter(:Γ) == true
        @test is_greek_letter(:Δ) == true
        @test is_greek_letter(:Π) == true
        @test is_greek_letter(:Σ) == true
        @test is_greek_letter(:Ω) == true
        
        # Test Greek letters with subscripts/additional characters
        @test is_greek_letter(:α₁) == true
        @test is_greek_letter(:β_bar) == true
        @test is_greek_letter(:δ_0) == true
        
        # Test non-Greek letters
        @test is_greek_letter(:a) == false
        @test is_greek_letter(:b) == false
        @test is_greek_letter(:x) == false
        @test is_greek_letter(:y) == false
        @test is_greek_letter(:k) == false
        @test is_greek_letter(:c) == false
        @test is_greek_letter(:foo) == false
        @test is_greek_letter(:bar123) == false
        
        # Test non-symbol inputs
        @test is_greek_letter("α") == false
        @test is_greek_letter(123) == false
        @test is_greek_letter(nothing) == false
        @test is_greek_letter(:(x + y)) == false
    end
    
    @testset "collect_all_symbols!" begin
        # Test simple symbol collection
        symbols = Set{Symbol}()
        collect_all_symbols!(symbols, :x)
        @test symbols == Set([:x])
        
        # Test operator exclusion
        symbols = Set{Symbol}()
        collect_all_symbols!(symbols, :+)
        @test isempty(symbols)
        
        # Test special symbol exclusion
        symbols = Set{Symbol}()
        collect_all_symbols!(symbols, :𝑑)
        @test isempty(symbols)
        
        # Test expression with multiple symbols
        symbols = Set{Symbol}()
        collect_all_symbols!(symbols, :(α * k^β - δ * k + exp(γ)))
        @test symbols == Set([:α, :k, :β, :δ, :exp, :γ])
        
        # Test function calls
        symbols = Set{Symbol}()
        collect_all_symbols!(symbols, :(f(x, y)))
        @test symbols == Set([:f, :x, :y])
        
        # Test nested function calls
        symbols = Set{Symbol}()
        collect_all_symbols!(symbols, :(log(exp(x))))
        @test symbols == Set([:log, :exp, :x])
        
        # Test complex expression
        symbols = Set{Symbol}()
        collect_all_symbols!(symbols, :(c(k) * (α * k^(α - 1) - δ - ρ)))
        @test symbols == Set([:c, :k, :α, :δ, :ρ])
        
        # Test expression with operators only
        symbols = Set{Symbol}()
        collect_all_symbols!(symbols, :(2 + 3 * 4))
        @test isempty(symbols)
        
        # Test empty expression handling
        symbols = Set{Symbol}()
        collect_all_symbols!(symbols, nothing)
        @test isempty(symbols)
    end
    
    @testset "find_function_calls!" begin
        # Test simple function call
        funcs = Set{Symbol}()
        find_function_calls!(funcs, :(f(x)))
        @test funcs == Set([:f])
        
        # Test multiple function calls
        funcs = Set{Symbol}()
        find_function_calls!(funcs, :(f(x) + g(y)))
        @test funcs == Set([:f, :g])
        
        # Test nested function calls
        funcs = Set{Symbol}()
        find_function_calls!(funcs, :(f(g(h(x)))))
        @test funcs == Set([:f, :g, :h])
        
        # Test vectorized function calls
        funcs = Set{Symbol}()
        find_function_calls!(funcs, :(f.(x)))
        @test funcs == Set([:f])
        
        # Test operators are not included
        funcs = Set{Symbol}()
        find_function_calls!(funcs, :(x + y))
        @test isempty(funcs)
        
        # Test complex expression from economics
        funcs = Set{Symbol}()
        find_function_calls!(funcs, :(c(k) * (α * k^(α - 1) - δ - ρ)))
        @test funcs == Set([:c])
        
        # Test expression with multiple function types
        funcs = Set{Symbol}()
        find_function_calls!(funcs, :(f(x) + log(y) - exp(z)))
        @test funcs == Set([:f, :log, :exp])
        
        # Test syntactic operators
        funcs = Set{Symbol}()
        find_function_calls!(funcs, :(x && y || z))
        @test isempty(funcs)
        
        # Test no function calls
        funcs = Set{Symbol}()
        find_function_calls!(funcs, :(x + y * z))
        @test isempty(funcs)
        
        # Test with symbols only
        funcs = Set{Symbol}()
        find_function_calls!(funcs, :x)
        @test isempty(funcs)
        
        # Test with nothing
        funcs = Set{Symbol}()
        find_function_calls!(funcs, nothing)
        @test isempty(funcs)
    end
    
    @testset "parse_differential_equations" begin
        # Test Ramsey model
        block = quote
            𝑑c = c(k) * (α * k^(α - 1) - δ - ρ)
            𝑑k = k^α - δ * k - c(k)
        end
        states, costates, params, odes = parse_differential_equations(block)
        @test states == [:k]
        @test costates == [:c]
        @test params == [:α, :δ, :ρ]
        @test length(odes) == 2
        @test odes[1] == :(k^α - δ * k - c(k))  # State equation first
        @test odes[2] == :(c(k) * (α * k^(α - 1) - δ - ρ))  # Costate equation second
        
        # Test Solow model (no costates)
        block = quote
            𝑑k = s * k^α - δ * k
        end
        states, costates, params, odes = parse_differential_equations(block)
        @test states == [:k]
        @test isempty(costates)
        @test params == [:α, :δ]
        @test length(odes) == 1
        @test odes[1] == :(s * k^α - δ * k)
        
        # Test multiple state variables
        block = quote
            𝑑x = α * x - β * y
            𝑑y = γ * x + δ * y
        end
        states, costates, params, odes = parse_differential_equations(block)
        @test states == [:x, :y]
        @test isempty(costates)
        @test params == [:α, :β, :γ, :δ]
        @test length(odes) == 2
        
        # Test with control variables
        block = quote
            𝑑x = u(x) + α * x
            𝑑u = β * u - γ * x
        end
        states, costates, params, odes = parse_differential_equations(block)
        @test states == [:x]
        @test costates == [:u]
        @test params == [:α, :β, :γ]
        @test length(odes) == 2
        @test odes[1] == :(u(x) + α * x)  # State equation first
        @test odes[2] == :(β * u - γ * x)  # Costate equation second
        
        # Test with global functions (not costates)
        block = quote
            𝑑k = f(k) - δ * k
            𝑑c = c * g(k) - ρ * c
        end
        states, costates, params, odes = parse_differential_equations(block)
        @test states == [:k, :c]
        @test isempty(costates)
        @test params == [:δ, :ρ]
        @test length(odes) == 2
        
        # Test empty block
        block = quote
        end
        states, costates, params, odes = parse_differential_equations(block)
        @test isempty(states)
        @test isempty(costates)
        @test isempty(params)
        @test isempty(odes)
        
        # Test with line numbers (should be ignored)
        block = quote
            𝑑x = α * x
        end
        # The quote block includes LineNumberNodes
        states, costates, params, odes = parse_differential_equations(block)
        @test states == [:x]
        @test params == [:α]
        
        # Test error on non-block expression
        @test_throws AssertionError parse_differential_equations(:(𝑑x = α * x))
        
        # Test with uppercase Greek parameters
        block = quote
            𝑑x = Γ * x - Δ * x^2
        end
        states, costates, params, odes = parse_differential_equations(block)
        @test states == [:x]
        @test params == [:Γ, :Δ]
        
        # Test preservation of equation order for states and costates
        block = quote
            𝑑a = a + b(a)
            𝑑c = c + d(c)
            𝑑b = b * α
            𝑑d = d * β
        end
        states, costates, params, odes = parse_differential_equations(block)
        @test states == [:a, :c]
        @test costates == [:b, :d]
        @test params == [:α, :β]
        # ODEs should be ordered: states first, then costates
        @test odes[1] == :(a + b(a))  # First state
        @test odes[2] == :(c + d(c))  # Second state
        @test odes[3] == :(b * α)     # First costate
        @test odes[4] == :(d * β)     # Second costate
    end
    
    @testset "vectorize_function_calls" begin
        # Test basic operators
        @test vectorize_function_calls(:(x + y)) == :(x .+ y)
        @test vectorize_function_calls(:(x - y)) == :(x .- y)
        @test vectorize_function_calls(:(x * y)) == :(x .* y)
        @test vectorize_function_calls(:(x / y)) == :(x ./ y)
        @test vectorize_function_calls(:(x^2)) == :(x.^2)
        
        # Test function calls
        @test vectorize_function_calls(:(log(x))) == :(log.(x))
        @test vectorize_function_calls(:(exp(x))) == :(exp.(x))
        @test vectorize_function_calls(:(sin(x))) == :(sin.(x))
        
        # Test compound expressions
        @test vectorize_function_calls(:(log(x) + log(y))) == :(log.(x) .+ log.(y))
        @test vectorize_function_calls(:(x^2 + y^2)) == :(x.^2 .+ y.^2)
        @test vectorize_function_calls(:(k^α - δ * k)) == :(k.^α .- δ .* k)
        
        # Test nested function calls
        @test vectorize_function_calls(:(log(exp(x)))) == :(log.(exp.(x)))
        @test vectorize_function_calls(:(sin(cos(x)))) == :(sin.(cos.(x)))
        
        # Test already vectorized operations (should be preserved)
        @test vectorize_function_calls(:(log.(x))) == :(log.(x))
        @test vectorize_function_calls(:(x .+ y)) == :(x .+ y)
        @test vectorize_function_calls(:(x.^2)) == :(x.^2)
        
        # Test syntactic operators
        @test vectorize_function_calls(:(x && y)) == :(x .&& y)
        @test vectorize_function_calls(:(x || y)) == :(x .|| y)
        
        # Test comparison operators
        @test vectorize_function_calls(:(x < y)) == :(x .< y)
        @test vectorize_function_calls(:(x > y)) == :(x .> y)
        @test vectorize_function_calls(:(x == y)) == :(x .== y)
        @test vectorize_function_calls(:(x != y)) == :(x .!= y)
        @test vectorize_function_calls(:(x <= y)) == :(x .<= y)
        @test vectorize_function_calls(:(x >= y)) == :(x .>= y)
        
        # Test with multiple arguments
        @test vectorize_function_calls(:(f(x, y))) == :(f.(x, y))
        @test vectorize_function_calls(:(max(x, y))) == :(max.(x, y))
        
        # Test complex economics expression
        @test vectorize_function_calls(:(c(k) * (α * k^(α - 1) - δ - ρ))) == 
              :(c.(k) .* (α .* k.^(α .- 1) .- δ .- ρ))
        
        # Test that symbols alone are not modified
        @test vectorize_function_calls(:x) == :x
        @test vectorize_function_calls(:α) == :α
        
        # Test with numbers
        @test vectorize_function_calls(:(1 + x)) == :(1 .+ x)
        @test vectorize_function_calls(:(x + 2.5)) == :(x .+ 2.5)
        
        # Test parentheses preservation
        @test vectorize_function_calls(:((x + y) * z)) == :((x .+ y) .* z)
        @test vectorize_function_calls(:(x * (y + z))) == :(x .* (y .+ z))
    end
end