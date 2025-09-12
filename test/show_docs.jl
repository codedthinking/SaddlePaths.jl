#!/usr/bin/env julia

push!(LOAD_PATH, "../src")
using SaddlePaths

# Show module documentation
println("=" ^ 60)
println("MODULE DOCUMENTATION")
println("=" ^ 60)
@doc SaddlePaths

println("\n" * "=" ^ 60)
println("SADDLEPATH STRUCT DOCUMENTATION")
println("=" ^ 60)
@doc SaddlePath

println("\n" * "=" ^ 60)
println("@MODEL MACRO DOCUMENTATION (excerpt)")
println("=" ^ 60)
# Get the docstring for @model
docs = Base.Docs.doc(SaddlePaths.var"@model")
docs_str = string(docs)
# Print first 50 lines
for (i, line) in enumerate(split(docs_str, '\n')[1:min(50, end)])
    println(line)
end