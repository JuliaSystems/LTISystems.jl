# definition of discrete rational transfer function

immutable DSisoRational{T1<:Real, T2<:Real} <: DSisoTf
  num::Poly{T1}
  den::Poly{T2}
  Ts::Float64

  function call{T1<:Real, T2<:Real}(::Type{DSisoRational}, num::Poly{T1},
      den::Poly{T2}, Ts::Float64)

    Ts_ = Ts > zero(Float64) ? Ts : NaN
    new{T1,T2}(num, den, Ts_)
  end
end

# creation of discrete rational transfer functions

function tf{V1<:AbstractVector, V2<:AbstractVector, T3<:Real}(num::V1, den::V2,
    Ts::T3; var::Symbol=:z)
  @assert eltype(num) <: Number string("num must be vector of T<:Number elements")
  @assert eltype(den) <: Number string("den must be vector of T<:Number elements")
  @assert var == :q || var == :z || var == :zinv || var == :qinv
    string("var need to be either :q, :z, :qinv or :zinv")

  if var == :zinv || var == :qinv
    numlast = findlast(num)
    denlast = findlast(den)
    order = max(numlast,denlast)
    num_ = zeros(eltype(num),order)
    num_[1:numlast] = num[1:numlast]
    den_ = zeros(eltype(den),order)
    den_[1:denlast] = den[1:denlast]
    num = num_
    den = den_
  end

  pnum = Poly(num[end:-1:1])
  pden = Poly(den[end:-1:1])
  DSisoRational(pnum, pden, Float64(Ts))
end

tf{T1<:Real, T2<:Real, T3<:Real}(num::Poly{T1}, den::Poly{T2}, Ts::T3) = DSisoRational(num, den, Float64(Ts))
tf{T1<:Real, T2<:Real}(gain::T1, Ts::T2)                               = DSisoRational(Poly([gain]), Poly([one(T1)]), Float64(Ts))


# conversion and promotion

promote_rule{T11,T12,T21,T22}(::Type{DSisoRational{T11,T12}},
  ::Type{DSisoRational{T21,T22}})                              = DSisoRational
promote_rule{T<:Real,T11,T12}(::Type{DSisoRational{T11,T12}},
::Type{T})                                                     = DSisoRational
promote_rule{T1,T<:Real}(::Type{DSisoRational{T1}}, ::Type{T}) = DSisoRational

convert{T11,T12,T21,T22}(::Type{DSisoRational{T11,T12}},
  s::DSisoRational{T21,T22})                                   = s
convert{T<:Real}(::Type{DSisoRational}, x::T)                  = tf([x], [one(T)], NaN64)

# overloading identities

zero(::Type{DSisoRational})       = tf([zero(Int8)], [one(Int8)], NaN64)
zero(s::DSisoRational)            = tf([zero(Int8)], [one(Int8)], NaN64)
one(::Type{DSisoRational})        = tf([one(Int8)], [one(Int8)], NaN64)
one(s::DSisoRational)             = tf([one(Int8)], [one(Int8)], NaN64)

# interface implementation

zeros(s::DSisoRational)           = convert(Vector{Complex{Float64}}, roots(s.num))
poles(s::DSisoRational)           = convert(Vector{Complex{Float64}}, roots(s.den))
numvec(s::DSisoRational)          = reverse(coeffs(s.num))
denvec(s::DSisoRational)          = reverse(coeffs(s.den))
numpoly(s::DSisoRational)         = copy(s.num)
denpoly(s::DSisoRational)         = copy(s.den)
zpkdata(s::DSisoRational)         = (zeros(s), poles(s), s.num[end]/s.den[end])
samplingtime(s::DSisoRational)    = s.Ts
isdiscrete(s::DSisoRational)      = true
isdiscrete(::Type{DSisoRational}) = true

# overload printing functions

function show(io::IO, s::DSisoRational)
  println(io, "Discrete time rational transfer function model")
  println(io, "\ty = Gu")
  if s.Ts > 0
    println(io, "with Ts=", s.Ts, ".")
  elseif s.Ts == 0
    println(io, "with Ts=unspecified.")
  end
end

function showall(io::IO, s::DSisoRational)
  show(io, s)
  println(io, "")
  printtransferfunction(io::IO, s)
end

# overload mathematical operations

function +(s1::DSisoRational, s2::DSisoRational)
  Ts::Float64
  if s1.Ts == s2.Ts || isnan(s2.Ts)
    Ts = s1.Ts
  elseif isnan(s1.Ts)
    Ts = s2.Ts
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
  den1,den2,dengcd   = rmgcd(s1.den, s2.den)
  tf(s1.num*den2 + s2.num*den1, den1*den2*dengcd, Ts)
end
+{T<:Real}(s::DSisoRational, n::T)       = tf(s.num + n*s.den, copy(s.den), s.Ts)
+{T<:Real}(n::T, s::DSisoRational)       = s + n

.+{T<:Real}(s::DSisoRational, n::T)      = s + n
.+{T<:Real}(n::T, s::DSisoRational)      = s + n
.+(s1::DSisoRational, s2::DSisoRational) = +(s1,-s2)

-{T}(s::DSisoRational{T})                = tf(-s.num, copy(s.den), s.Ts)

-(s1::DSisoRational, s2::DSisoRational)  = +(s1,-s2)
-{T<:Real}(n::T, s::DSisoRational)       = +(n, -s)
-{T<:Real}(s::DSisoRational, n::T)       = +(s, -n)

.-{T<:Real}(s::DSisoRational, n::T)      = -(s, n)
.-{T<:Real}(n::T, s::DSisoRational)      = -(n, s)
.-(s1::DSisoRational, s2::DSisoRational) = -(s1, s2)

function *(s1::DSisoRational, s2::DSisoRational)
  Ts::Float64
  if s1.Ts == s2.Ts || isnan(s2.Ts)
    Ts = s1.Ts
  elseif isnan(s1.Ts)
    Ts = s2.Ts
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
  num1,den2,gcd1 = rmgcd(s1.num, s2.den)
  den1,num2,gcd2 = rmgcd(s1.den, s2.num)
  tf(num1*num2, den1*den2, Ts)
end
*{T<:Real}(s::DSisoRational, n::T)       = tf(s.num*n, copy(s.den), s.Ts)
*{T<:Real}(n::T, s::DSisoRational)       = *(s, n)

.*{T<:Real}(s::DSisoRational, n::T)      = *(s, n)
.*{T<:Real}(n::T, s::DSisoRational)      = *(n, s)
.*(s1::DSisoRational, s2::DSisoRational) = *(s1, s2)

/(s1::DSisoRational, s2::DSisoRational)  = s1*(1/s2)
/{T<:Real}(n::T, s::DSisoRational)       = tf(n*s.den, copy(s.num), s.Ts)
/{T<:Real}(s::DSisoRational, n::T)       = s*(1/n)

./{T<:Real}(n::T, s::DSisoRational)      = /(n, s)
./{T<:Real}(s::DSisoRational, n::T)      = /(s, n)
./(s1::DSisoRational, s2::DSisoRational) = /(s1, s2)

function ==(s1::DSisoRational, s2::DSisoRational)
  s1.Ts == s2.Ts && s1.num == s2.num &&
    (s1.den == s2.den || s1.num == zero(s1.num))
    # TODO scaling of num and den
end

!=(s1::DSisoRational, s2::DSisoRational) = !(s1 == s2)

function isapprox(s1::DSisoRational, s2::DSisoRational,
    rtol::Real=sqrt(eps()), atol::Real=0)
  sdiff = s2-s1
  return norm(sdiff.num) < rtol
end
