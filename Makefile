.PHONY: test docs install clean

# Run tests
test:
	julia --project=. -e "using Pkg; Pkg.test()"

# Build documentation
docs:
	julia --project=docs/ docs/make.jl

# Install dependencies
install:
	julia --project=. -e "using Pkg; Pkg.instantiate()"

# Clean generated files
clean:
	rm -rf docs/build/
	find . -name "*.cov" -delete
	find . -name "*.mem" -delete

# Format code
format:
	julia --project=. -e "using JuliaFormatter; format(\".\")"

# Run all checks (tests + docs)
check: test docs