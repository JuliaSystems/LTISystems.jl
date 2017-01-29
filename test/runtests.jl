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

# Defining left MFDs
s = Poly([0, 1],:s)
N = PolyMatrix([s; -s])
D = PolyMatrix([0 -(s^3+4*s^2+5*s+2); (s+2)^2 s+2])
M = lfd(N,D)
@test islfd(M)
@test !isrfd(M)

# INCOMPLETE
# Conversion from SS to MFD
A = [-4 -4 0 -1 -2; 1 0 0 0 0; 0 0 -4 -5 -2; 0 0 1 0 0; 0 0 0 1 0]'
B = [1 0 0 0 0; -1 0 1 0 0]'
C = [0 1; 0 0; -1 0; 0 0; 0 0]'
D = zeros(2,2)
sys = ss(A,B,C,D)
m = rfd(sys)
m.N
m.D
cond([B A*B A^2*B A^3*B A^4*B])
