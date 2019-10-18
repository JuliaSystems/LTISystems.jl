using Polynomials
using RationalFunctions
using PolynomialMatrices
using LinearAlgebra

using LTISystems

using Test
#using Test.@testset

# conversions
include("conversions/tf2ss.jl")

# methods
include("methods/simulation.jl")

# types
include("types/system/ltisystem.jl")
include("types/system/transferfunction.jl")
include("types/system/statespace.jl")
include("types/system/mfd.jl")
