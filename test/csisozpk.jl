println("Starting continuous zpk transfer function type tests...")

# Construction
s1 = zpk([-2], [1, 1], 1)
s2 = zpk([-2], [1., 1.], 1)
s3 = zpk(5.0)
s4 = zpk([1//2], [1, 1], 1)
s5 = zpk([1], [1], 3)

# print functions
show(s1)
showall(s1)
showcompact(s1)

@test typeof(s1)     == ControlCore.CSisoZpk{Int,Int,Int}
@test typeof(s2)     == ControlCore.CSisoZpk{Int,Float64,Int}
@test typeof(s3)     == ControlCore.CSisoZpk{Int8,Int8,Float64}
@test typeof(s4)     == ControlCore.CSisoZpk{Rational{Int},Int,Int}

# Conversions
@test promote_type(typeof(s1), Float64)           == ControlCore.CSisoZpk
@test promote_type(ControlCore.CSisoZpk, Float64) == ControlCore.CSisoZpk

@test convert(ControlCore.CSisoZpk, one(Float64)) ≈  zpk(one(Float64))

for s in [s1,s2]
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
  @test zeros(s)        == [-2]
  @test poles(s)        == [1, 1]

# return vectors
  @test numvec(s)       == coeffs(poly(s.z))[end:-1:1]
  @test denvec(s)       == coeffs(poly(s.p))[end:-1:1]
  @test numpoly(s)      == poly(s.z)
  @test denpoly(s)      == poly(s.p)
  @test zpkdata(s)      == (s.z, s.p, s.k)
  @test samplingtime(s) == 0
end

# addition

@test s1+s1 == zpk([-2], [1, 1], 2)
@test s1+s2 == zpk([-2], [1, 1], 2)
@test s1+s3 ≈  zpk([(9+sqrt(59)*im)/10, (9-sqrt(59)*im)/10], [1,1], 5)
@test s1+s4 ≈  zpk([-3//4], [1, 1], 2)

@test s1+1 ≈  zpk([(1+sqrt(11)*im)/2, (1-sqrt(11)*im)/2], [1, 1], 1)
@test s1+0 == s1

@test s1+s4 == s4+s1
@test 1+s1  == s1+1

# subtraction

@test -s1   == zpk([-2], [1, 1], -1)
@test s1-s1 ≈  zero(s1)
@test s1-s2 ≈  zero(s1)
@test s1-s3 ≈  zpk([1/10*(11-sqrt(61)), 1/10*(11+sqrt(61))], [1, 1], -5)
@test s1-s4 ≈  zpk(Float64[], [1, 1], 2.5)

@test s1-1 ≈  zpk([(3+sqrt(13))/2, (3-sqrt(13))/2], [1, 1], -1)
@test s1+0 == s1

@test s1-s4 == -(s4-s1)
@test 1-s1  == -(s1-1)

# multiplication

@test s1*s1 == zpk([-2, -2], ones(4), 1)
@test s1*s2 == zpk([-2, -2], ones(4), 1)
@test s1*s3 == zpk([-2.], ones(2), 5.0)
@test s1*s4 == zpk([-2, 0.5], ones(4), 1)

@test s1*1 ≈  s1
@test s1*0 ≈  zero(s1)

@test s1*s4 ≈  s4*s1
@test 1*s1  == s1*1

# division

@test s1/s1 ≈  one(s1)
@test s1/s2 ≈  one(s1)
@test s1/s3 == zpk([-2.], ones(2), 0.2)
@test s1/s4 ≈  zpk([-2.], [.5], 1)

@test s1/1 ≈  s1
@test isinf((s1/0).k)

@test s1/s4 ≈  1/(s4/s1)
@test 1/s1  ≈  1/(s1/1)

# equality

@test s1 == s1
@test s1 == s2
@test s1 != s3
@test s1 != s4
@test s1 != s5

# isapprox

# TODO
