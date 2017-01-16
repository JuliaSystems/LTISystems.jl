using SystemsBase
using Base.Test
using PolynomialMatrices
using Polynomials

# continuous state space type tests
#include("continuousss.jl")
#include("dsisorational.jl")
#include("csisorational.jl")
#include("dsisozpk.jl")
#include("csisozpk.jl")
#include("dmimo.jl")
#include("cmimo.jl")

## MFDs
s = Poly([0, 1],:s)
N = PolyMatrix([s; -s])
D = PolyMatrix([0 -(s^3+4*s^2+5*s+2); (s+2)^2 s+2])
M = lfd(N,D)
@test islfd(M)
@test !isrfd(M)
