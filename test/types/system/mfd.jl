## MatrixFractionDescriptions
using Polynomials
using PolynomialMatrices

# Defining left MatrixFractionDescriptions
s = Poly([0, 1],:s)
N = PolyMatrix([s; -s])
D = PolyMatrix([0 -(s^3+4*s^2+5*s+2); (s+2)^2 s+2])
M = lfd(N,D)
@test islfd(M)
@test !isrfd(M)

# isproper
s = Poly([0, 1],:s)
N = PolyMatrix([2s^2+1 2])
D = PolyMatrix([s^3+s s; s^2+s+1 1])
M = rfd(N,D)

@test !is_col_proper(M.D)
#@test isproper(M)

# Defining right MatrixFractionDescriptions
N = PolyMatrix([s -s])
D = PolyMatrix([zero(s) -(s^3+4*s^2+5*s+2); (s+2)^2 s+2])
M = rfd(N,D)
@test !islfd(M)
@test isrfd(M)

# Defining lfd2rfd and rfd2lfd
s = variable("s")
Nₗ = PolyMatrix([s^2 zero(s); -4s s])
Dₗ = PolyMatrix([s^3+2s^2-1 s+1; -5s^2-13s-8 (s+1)*(s+4)])
Nᵣ = PolyMatrix([-s^2 -s; zero(s) -s])
Dᵣ = PolyMatrix([-s^3-2s^2+1 -(s+1)^2; (s+2)^2*(s+1) zero(s)])

rmfd = rfd(Nᵣ, Dᵣ)
lmfd = lfd(Nₗ, Dₗ)
#@test isapprox(rfd(lmfd), rmfd)
#@test isapprox(rfd(lfd(rmfd)), rmfd)

# INCOMPLETE
# Conversion from SS to MatrixFractionDescription
A = [-4 -4 0 -1 -2; 1 0 0 0 0; 0 0 -4 -5 -2; 0 0 1 0 0; 0 0 0 1 0]'
B = [1 0 0 0 0; -1 0 1 0 0]'
C = [0 1; 0 0; -1 0; 0 0; 0 0]'
D = zeros(Int,2,2)
sys = ss(A,B,C,D)
# mr = rfd(sys)
# mr.N
# mr.D
# cond([B A*B A^2*B A^3*B A^4*B])
#
# ml = lfd(sys)
# ml.N
# ml.D
#
# @test isapprox(mr, ml; rtol=10)
