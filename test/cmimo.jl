println("Starting continuous MIMO transfer function type tests...")

# Construction
nsisosystem = 2
s1m = Array{ControlCore.SisoTf}(2,1)
s2m = Array{ControlCore.SisoTf}(2,1)
s3m = Array{ControlCore.SisoTf}(2,1)
s4m = Array{ControlCore.SisoTf}(2,1)

# zpk construction
s1m[1] = zpk([-2], [1, 1], 1)
s2m[1] = zpk([-2], [1., 1.], 1)
s3m[1] = zpk(5.0)
s4m[1] = zpk([1//2], [1, 1], 1)

zm = Array{Array{Int64},2}(2,1)
pm = Array{Array{Float64},2}(2,1)
km = Array{Float64}(2,1)
zm[1], pm[1], km[1] = zpkdata(s1m[1])
zm[2], pm[2], km[2] = zpkdata(s2m[1])

@test zpk(zm, pm, km) ≈ mimo(reshape([s1m[1] s2m[1]], 2, 1))
@test zpk(5*eye(2)) ≈ mimo( [[s3m[1] zero(s3m[1])]; [zero(s3m[1]) s3m[1]]])
@test_throws DomainError zpk(zm, pm.', km)

# tf construction
s1m[2] = tf([1, 2], [1, 2, 1])
s2m[2] = tf([1, 2], [1., 2., 1])
s3m[2] = tf(5.0)
s4m[2] = tf([1//2, 3//5], [1, 1, 2])

numm = Array{Array{Int64},2}(2,1)
denm = Array{Array{Float64},2}(2,1)
numm[1] = numvec(s1m[2])
numm[2] = numvec(s2m[2])
denm[1] = denvec(s1m[2])
denm[2] = denvec(s2m[2])

@test tf(numm, denm) ≈ mimo(reshape([s1m[2] s2m[2]],2,1))
@test tf(5*eye(2)) ≈ mimo( [[s3m[2] zero(s3m[2])]; [zero(s3m[2]) s3m[2]]])
@test_throws DomainError tf(numm, denm.')

# print functions
show(s1m[2])
showall(s1m[2])
showcompact(s1m[2])

for kk in 1:nsisosystem
  s1 = s1m[kk]
  s2 = s2m[kk]
  s3 = s3m[kk]
  s4 = s4m[kk]

  m1 = mimo([[s1 s1]; [s1 s1]])
  m2 = mimo([[s1 s2]; [s3 s4]])

  # Conversions
  @test promote_type(typeof(m1), typeof(m2)) == ControlCore.CMimo
  @test promote_type(typeof(m1), Array{Float64,2}) == ControlCore.CMimo
  @test promote_type(typeof(m1), Float64) == ControlCore.CMimo

  @test convert(ControlCore.CMimo, one(Float64))  ≈ one(m1)
  @test convert(ControlCore.CMimo, zero(Float64)) ≈ zero(m1)
  @test convert(ControlCore.CMimo, eye(2)) ≈ zpk(eye(2))

  for s in [m1,m2]
    # I/O mapping
    @test numstates(s)    == [[numstates(s[1]) numstates(s[3])] ;
      [numstates(s[2]) numstates(s[4])]]
    @test numinputs(s)    == 2
    @test numoutputs(s)   == 2

    # Dimension information
    @test ndims(s)        == 2
    @test size(s)         == (2,2)
    @test size(s, 2)      == 2
    @test size(s, 1, 2)   == (2,2)

    # Iteration interface
    @test s[1]            ≈  s[1]
    @test s[:]            ≈  mimo(s.m[:])
    @test_throws BoundsError s[5]
  end

  # addition

  @test m1+m2 ≈ mimo([[s1+s1 s1+s2] ; [s1+s3 s1+s4]])
  @test m1+m2 ≈ m2+m1

  @test m1+eye(2)     ≈ mimo([[s1+1 s1+0] ; [s1+0 s1+1]])
  @test m1+eye(2)     ≈ eye(2)+m1
  @test m1+zeros(2,2) ≈ m1
  @test m1+zeros(2,2) ≈ zeros(2,2)+m1

  # subtraction

  @test -m1   ≈ mimo(-getmatrix(m1))

  @test m1-m2 ≈ mimo([[s1-s1 s1-s2] ; [s1-s3 s1-s4]])
  @test m1-m2 ≈ -(m2-m1)

  @test m1-eye(2)     ≈ mimo([[s1-1 s1-0] ; [s1-0 s1-1]])
  @test m1-eye(2)     ≈ eye(2)-m1
  @test m1-zeros(2,2) ≈ m1
  @test m1-zeros(2,2) ≈ zeros(2,2)-m1

  # multiplication

  @test m1*m2 ≈ mimo([[s1*s1+s1*s3   s1*s2+s1*s4] ;
                 [s1*s1+s1*s3   s1*s2+s1*s4]])
  @test m1*m2 ≈ (m2.'*m1.').'

  @test m1*eye(2) ≈  m1
  @test m1*eye(2) ≈ eye(2)*m1
  @test m1*zeros(2,2) ≈ mimo(getmatrix(m1)*0)
  @test m1*zeros(2,2) ≈ zeros(2,2)*m1

  # division

  @test_throws DomainError  m1/m2

  @test m1/eye(2) ≈  m1

  # equality

  @test m1 == m1
  @test m1 != m2

  # isapprox

  # TODO

end
