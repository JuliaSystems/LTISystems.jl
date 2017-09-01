using Polynomials
using RationalFunctions
using PolynomialMatrices

using LTISystems

using Base.Test
using Base.Test.@testset

# conversions
include("conversions/tf2ss.jl")

# methods
include("methods/simulation.jl")

# types
include("types/system/ltisystem.jl")
include("types/system/transferfunction.jl")
include("types/system/statespace.jl")
include("types/system/mfd.jl")
