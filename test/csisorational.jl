println("Starting continuous SISO rational transfer function type tests...")

# Construction
s1 = tf([1, 2], [1, 2, 1])
s2 = tf([1, 2], [1., 2., 1])
s3 = tf(Poly([Float64(2), 1.0]), Poly([1.0, Float64(2), 1]))
s4 = tf(5.0)
s5 = tf([1//2, 3//5], [1, 1, 2])

# print functions
show(s1)
showall(s1)
showcompact(s1)

@test typeof(s1)     == ControlCore.CSisoRational{Int,Int}
@test typeof(s2)     == ControlCore.CSisoRational{Int,Float64}
@test typeof(s3)     == ControlCore.CSisoRational{Float64,Float64}
@test typeof(s4)     == ControlCore.CSisoRational{Float64,Float64}
@test typeof(s5)     == ControlCore.CSisoRational{Rational{Int},Int}

# Conversions
@test promote_type(typeof(s1), Float64) == ControlCore.CSisoRational
@test promote_type(ControlCore.CSisoRational{Int}, Float64) == ControlCore.CSisoRational

@test convert(ControlCore.CSisoRational, one(Float64)) ≈ tf(one(Float64))

for s in [s1,s2,s3]
# I/O mapping
  @test numstates(s)    == 3
  @test numinputs(s)    == 1
  @test numoutputs(s)   == 1

# Dimension information
  @test ndims(s)        == 1
  @test size(s)         == 1
  @test size(s, 2)      == 1
  @test size(s, 1, 2)   == (1,1)

# Iteration interface
  @test s[1]            == s
  @test s[:]            == s
  @test_throws BoundsError s[2]

# poles and zeros
  @test zeros(s)         ≈ Array{Complex{eltype(s.num)},1}([-2+0im])
  @test poles(s)         ≈ Array{Complex{eltype(s.den)},1}([-1+0im, -1+0im])

# return vectors
  @test numvec(s)       == coeffs(s.num)[end:-1:1]
  @test denvec(s)       == coeffs(s.den)[end:-1:1]
  @test numpoly(s)      == s.num
  @test denpoly(s)      == s.den
  @test zpkdata(s)      == (zeros(s), poles(s), s.num[end]/s.den[end])
end

# addition

@test s1+s1 == tf([2, 4], [1, 2, 1])
@test s1+s2 == tf([2, 4], [1, 2, 1])
@test s1+s3 == tf([2, 4], [1, 2, 1])
@test s1+s4 == tf([5,11,7],[1,2,1])
@test s1+s5 ≈  tf([3//2, 23//5, 57//10, 23//5], [1, 3, 5, 5, 2])

@test s1+1 == tf([1, 3, 3], [1, 2, 1])
@test s1+0 == s1

@test s1+s4 ≈  s4+s1
@test 1+s1  ≈  s1+1

# subtraction

@test -s1   == tf([-1, -2], [1, 2, 1])
@test s1-s1 ≈  zero(s1)
@test s1-s2 ≈  zero(s1)
@test s1-s3 ≈  zero(s1)
@test s1-s4 == tf([-5,-9,-3.],[1.,2,1])
@test s1-s5 ≈  tf([1//2, 7//5, 23//10, 17//5], [1//1, 3//1, 5//1, 5//1, 2//1])

@test s1-1 == tf([-1, -1, 1], [1, 2, 1])
@test s1-0 == s1

@test s1-s4 ≈ -(s4-s1)
@test 1-s1  ≈ -(s1-1)

# multiplication

@test s1*s1 == tf([1., 4., 4.], [1., 4., 6., 4., 1.])
@test s1*s2 == tf([1., 4., 4.], [1., 4., 6., 4., 1.])
@test s1*s3 == tf([1., 4., 4.], [1., 4., 6., 4., 1.])
@test s1*s4 == tf([5., 10.],[1., 2., 1.])
@test s1*s5 ≈  tf([1//2, 8//5, 6//5], [1//1, 3//1, 5//1, 5//1, 2//1])

@test s1*1 == s1
@test s1*0 ≈  zero(s1)

@test s1*s4 ≈  s4*s1
@test 1*s1  == s1*1

# division

@test s1/s1 ≈  one(s1)
@test s1/s2 ≈  one(s1)
@test s1/s3 ≈  one(s1)
@test s1/s4 == tf([1., 2.], [5., 10., 5.])
@test s1/s5 ≈  tf([1, 3, 4, 4], [1//2, 8//5, 17//10, 3//5])

@test s1/1 == s1
@test s1/0 == tf([Inf, Inf], [1, 2, 1])

@test s1/s4 ≈  1/(s4/s1)
@test 1/s1  ≈  1/(s1/1)

# equality

@test s1 == s1
@test s1 == s2
@test s1 == s3
@test s1 != s4
@test s1 != s5

# isapprox

# TODO
