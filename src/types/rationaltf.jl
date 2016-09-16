immutable RationalTF{T,S,M1,M2} <: LtiSystem{T,S}
  num::M1
  den::M2
  Ts::Float64

  # Continuous-time, single-input-single-output rational transfer function model
  @compat function (::Type{RationalTF}){M1<:Poly,M2<:Poly}(num::M1, den::M2)
    @assert eltype(num) <: Real       "RationalTF: num must be a polynomial with real coefficients"
    @assert eltype(den) <: Real       "RationalTF: den must be a polynomial with real coefficients"
    @assert den != zero(den)          "RationalTF: den cannot be zero"
    @assert degree(num) ≤ degree(den) "RationalTF: system is not proper"

    new{Siso{true},Continuous{true},M1,M2}(num, den, zero(Float64))
  end

  # Discrete-time, single-input-single-output rational transfer function model
  @compat function (::Type{RationalTF}){M1<:Poly,M2<:Poly,M3<:Real}(num::M1,
    den::M2, Ts::M3)
    @assert eltype(num) <: Real         "RationalTF: num must be a polynomial with real coefficients"
    @assert eltype(den) <: Real         "RationalTF: den must be a polynomial with real coefficients"
    @assert den != zero(den)            "RationalTF: den cannot be zero"
    @assert degree(num) ≤ degree(den)   "RationalTF: system is not proper"
    @assert Ts ≥ zero(Ts) && !isinf(Ts) "StateSpace: Ts must be non-negative number"

    new{Siso{true},Continuous{false},M1,M2}(num, den, convert(Float64, Ts))
  end

  # Continuous-time, multi-input-multi-output rational transfer function model
  @compat function (::Type{RationalTF}){M1<:AbstractMatrix,
    M2<:AbstractMatrix}(num::M1, den::M2)
    @assert size(num) == size(den)    "RationalTF: num and den must have the same size"
    @assert eltype(num) <: Poly &&
      eltype(eltype(num)) <: Real     "RationalTF: num must be a matrix of polynomials with real coefficients"
    @assert eltype(den) <: Poly &&
      eltype(eltype(den)) <: Real     "RationalTF: den must be a matrix of polynomials with real coefficients"
    for idx in eachindex(num)
      @assert degree(num[idx]) ≤ degree(den[idx]) "RationalTF: system is not proper"
      @assert den[idx] != zero(den[idx])          "RationalTF: den cannot be zero"
    end

    new{Siso{false},Continuous{true},M1,M2}(num, den, zero(Float64))
  end

  # Discrete-time, multi-input-multi-output rational transfer function model
  @compat function (::Type{RationalTF}){M1<:AbstractMatrix,
    M2<:AbstractMatrix,M3<:Real}(num::M1, den::M2, Ts::M3)
    @assert size(num) == size(den)    "RationalTF: num and den must have the same size"
    @assert eltype(num) <: Poly &&
      eltype(eltype(num)) <: Real     "RationalTF: num must be a matrix of polynomials with real coefficients"
    @assert eltype(den) <: Poly &&
      eltype(eltype(den)) <: Real     "RationalTF: den must be a matrix of polynomials with real coefficients"
    for idx in eachindex(num)
      @assert degree(num[idx]) ≤ degree(den[idx]) "RationalTF: system is not proper"
      @assert den[idx] != zero(den[idx])          "RationalTF: den cannot be zero"
    end

    new{Siso{false},Continuous{false},M1,M2}(num, den, convert(Float64, Ts))
  end
end

# Legacy SisoTf
# # Printing functions
# summary(s::SisoTf) = string("tf(nx=", numstates(s), (isdiscrete(s) ?
#   string(",Ts=", samplingtime(s)) : ""), ")")
#
# showcompact(io::IO, s::SisoTf)  = print(io, summary(s))
# show(io::IO, s::SisoTf)         = print(io, summary(s))
# showall(io::IO, s::SisoTf)      = print(io, summary(s))
#
# function rmgcd{T1<:Number, T2<:Number}(p1::Poly{T1}, p2::Poly{T2})
#   R = promote_type(T1,T2)
#   gcdp1p2::Poly{R} = gcd(p1,p2)
#   gcdp1p2_::Poly{R} = gcdp1p2/gcdp1p2[end]
#   p1_::Poly = div(p1,gcdp1p2_)
#   p2_::Poly = div(p2,gcdp1p2_)
#   return (p1_,p2_,gcdp1p2_)
# end
#
# function rmcommon{T1<:AbstractVector, T2<:AbstractVector}(a::T1, b::T2)
#   T     = promote_type(eltype(T1),eltype(T2))
#   list  = Dict{T,Vector{Int}}()
#
#   l1    = length(a)
#   l2    = length(b)
#   s1    = l1 > l2 ? b : a
#   s2    = l1 > l2 ? a : b
#
#   for elem in s1
#     list[elem] = get(list, elem, Int[0, 0]) + [1, 0]
#   end
#
#   for elem in s2
#     try
#       list[elem] += [0, 1]
#     end
#   end
#
#   temp::Vector{T} = []
#   for key in keys(list)
#     append!(temp, fill(key, min(list[key]...)))
#   end
#   ar::T1 = copy(a)
#   br::T2 = copy(b)
#   for elem in temp
#     splice!(ar, findfirst(ar,elem))
#     splice!(br, findfirst(br,elem))
#   end
#   return (ar,br,temp)
# end

# Legacy CSisoRational
# immutable CSisoRational{T1<:Real, T2<:Real} <: CSisoTf
#   num::Poly{T1}
#   den::Poly{T2}
#
#   function call{T1<:Real, T2<:Real}(::Type{CSisoRational}, num::Poly{T1},
#       den::Poly{T2})
#
#     new{T1,T2}(num, den)
#   end
# end
#
# # creation of continuous rational transfer functions
#
# function tf{V1<:AbstractVector, V2<:AbstractVector}(num::V1, den::V2)
#   @assert eltype(num) <: Number string("num must be vector of T<:Number elements")
#   @assert eltype(den) <: Number string("den must be vector of T<:Number elements")
#   pnum = Poly(num[end:-1:1])
#   pden = Poly(den[end:-1:1])
#   CSisoRational(pnum, pden)
# end
#
# tf{T1<:Real, T2<:Real}(num::Poly{T1}, den::Poly{T2}) = CSisoRational(num, den)
# tf{T1<:Real}(gain::T1) = CSisoRational(Poly([gain]), Poly([one(T1)]))
#
# # conversion and promotion
#
# promote_rule{T11,T12,T21,T22}(::Type{CSisoRational{T11,T12}},
#   ::Type{CSisoRational{T21,T22}})                              = CSisoRational
# promote_rule{T<:Real,T11,T12}(::Type{CSisoRational{T11,T12}},
# ::Type{T})                                                     = CSisoRational
# promote_rule{T1,T<:Real}(::Type{CSisoRational{T1}}, ::Type{T}) = CSisoRational
#
# convert{T11,T12,T21,T22}(::Type{CSisoRational{T11,T12}},
#   s::CSisoRational{T21,T22})                                   = s
# convert{T<:Real}(::Type{CSisoRational}, x::T)                  = tf([x], [one(T)])
#
# # overloading identities
#
# zero(::Type{CSisoRational})       = tf([zero(Int8)], [one(Int8)])
# zero(s::CSisoRational)            = tf([zero(Int8)], [one(Int8)])
# one(::Type{CSisoRational})        = tf([one(Int8)], [one(Int8)])
# one(s::CSisoRational)             = tf([one(Int8)], [one(Int8)])
#
# # interface implementation
#
# zeros(s::CSisoRational)           = convert(Vector{Complex{Float64}}, roots(s.num))
# poles(s::CSisoRational)           = convert(Vector{Complex{Float64}}, roots(s.den))
# numvec(s::CSisoRational)          = reverse(coeffs(s.num))
# denvec(s::CSisoRational)          = reverse(coeffs(s.den))
# numpoly(s::CSisoRational)         = copy(s.num)
# denpoly(s::CSisoRational)         = copy(s.den)
# zpkdata(s::CSisoRational)         = (zeros(s), poles(s), s.num[end]/s.den[end])
# samplingtime(s::CSisoRational)    = zero(Float64)
# isdiscrete(s::CSisoRational)      = false
# isdiscrete(::Type{CSisoRational}) = false
#
# # overload printing functions
#
# function show(io::IO, s::CSisoRational)
#   println(io, "Continuous time rational transfer function model")
#   println(io, "\ty = Gu")
# end
#
# function showall(io::IO, s::CSisoRational)
#   show(io, s)
#   println(io, "")
#   printtransferfunction(io::IO, s)
# end
#
# # overload mathematical operations
#
# function +(s1::CSisoRational, s2::CSisoRational)
#   den1,den2,dengcd   = rmgcd(s1.den, s2.den)
#   tf(s1.num*den2 + s2.num*den1, den1*den2*dengcd)
# end
# +{T<:Real}(s::CSisoRational, n::T)       = tf(s.num + n*s.den, copy(s.den))
# +{T<:Real}(n::T, s::CSisoRational)       = s + n
#
# .+{T<:Real}(s::CSisoRational, n::T)      = s + n
# .+{T<:Real}(n::T, s::CSisoRational)      = s + n
# .+(s1::CSisoRational, s2::CSisoRational) = +(s1, s2)
#
# -(s::CSisoRational)                      = tf(-s.num, copy(s.den))
#
# -(s1::CSisoRational, s2::CSisoRational)  = +(s1,-s2)
# -{T<:Real}(n::T, s::CSisoRational)       = +(n, -s)
# -{T<:Real}(s::CSisoRational, n::T)       = +(s, -n)
#
# .-{T<:Real}(s::CSisoRational, n::T)      = -(s, n)
# .-{T<:Real}(n::T, s::CSisoRational)      = -(n, s)
# .-(s1::CSisoRational, s2::CSisoRational) = +(s1, -s2)
#
# function *(s1::CSisoRational, s2::CSisoRational)
#   num1,den2,gcd1   = rmgcd(s1.num, s2.den)
#   den1,num2,gcd2   = rmgcd(s1.den, s2.num)
#   tf(num1*num2, den1*den2)
# end
#
# *{T<:Real}(s::CSisoRational, n::T)       = tf(s.num*n, copy(s.den))
# *{T<:Real}(n::T, s::CSisoRational)       = *(s, n)
#
# .*{T<:Real}(s::CSisoRational, n::T)      = *(s, n)
# .*{T<:Real}(n::T, s::CSisoRational)      = *(n, s)
# .*(s1::CSisoRational, s2::CSisoRational) = *(s1, s2)
#
# /(s1::CSisoRational, s2::CSisoRational)  = s1*(1/s2)
# /{T<:Real}(n::T, s::CSisoRational)       = tf(n*s.den, copy(s.num))
# /{T<:Real}(s::CSisoRational, n::T)       = s*(1/n)
#
# ./{T<:Real}(n::T, s::CSisoRational)      = /(n, s)
# ./{T<:Real}(s::CSisoRational, n::T)      = /(s, n)
# ./(s1::CSisoRational, s2::CSisoRational) = /(s1, s2)
#
# function ==(s1::CSisoRational, s2::CSisoRational)
#   s1.num == s2.num && (s1.den == s2.den || s1.num == zero(s1.num))
# end
#
# !=(s1::CSisoRational, s2::CSisoRational) = !(s1 == s2)
#
# function isapprox(s1::CSisoRational, s2::CSisoRational,
#   rtol::Real=sqrt(eps()), atol::Real=0)
#   sdiff = s2-s1
#   norm(sdiff.num) < rtol
# end

# Legacy DSisoRational
# immutable DSisoRational{T1<:Real, T2<:Real} <: DSisoTf
#   num::Poly{T1}
#   den::Poly{T2}
#   Ts::Float64
#
#   function call{T1<:Real, T2<:Real}(::Type{DSisoRational}, num::Poly{T1},
#       den::Poly{T2}, Ts::Float64)
#
#     Ts_ = Ts > zero(Float64) ? Ts : NaN
#     new{T1,T2}(num, den, Ts_)
#   end
# end
#
# # creation of discrete rational transfer functions
#
# function tf{V1<:AbstractVector, V2<:AbstractVector, T3<:Real}(num::V1, den::V2,
#     Ts::T3; var::Symbol=:z)
#   @assert eltype(num) <: Number string("num must be vector of T<:Number elements")
#   @assert eltype(den) <: Number string("den must be vector of T<:Number elements")
#   @assert var == :q || var == :z || var == :zinv || var == :qinv
#     string("var need to be either :q, :z, :qinv or :zinv")
#
#   if var == :zinv || var == :qinv
#     numlast = findlast(num)
#     denlast = findlast(den)
#     order = max(numlast,denlast)
#     num_ = zeros(eltype(num),order)
#     num_[1:numlast] = num[1:numlast]
#     den_ = zeros(eltype(den),order)
#     den_[1:denlast] = den[1:denlast]
#     num = num_
#     den = den_
#   end
#
#   pnum = Poly(num[end:-1:1])
#   pden = Poly(den[end:-1:1])
#   DSisoRational(pnum, pden, Float64(Ts))
# end
#
# tf{T1<:Real, T2<:Real, T3<:Real}(num::Poly{T1}, den::Poly{T2}, Ts::T3) = DSisoRational(num, den, Float64(Ts))
# tf{T1<:Real, T2<:Real}(gain::T1, Ts::T2)                               = DSisoRational(Poly([gain]), Poly([one(T1)]), Float64(Ts))
#
#
# # conversion and promotion
#
# promote_rule{T11,T12,T21,T22}(::Type{DSisoRational{T11,T12}},
#   ::Type{DSisoRational{T21,T22}})                              = DSisoRational
# promote_rule{T<:Real,T11,T12}(::Type{DSisoRational{T11,T12}},
# ::Type{T})                                                     = DSisoRational
# promote_rule{T1,T<:Real}(::Type{DSisoRational{T1}}, ::Type{T}) = DSisoRational
#
# convert{T11,T12,T21,T22}(::Type{DSisoRational{T11,T12}},
#   s::DSisoRational{T21,T22})                                   = s
# convert{T<:Real}(::Type{DSisoRational}, x::T)                  = tf([x], [one(T)], NaN64)
#
# # overloading identities
#
# zero(::Type{DSisoRational})       = tf([zero(Int8)], [one(Int8)], NaN64)
# zero(s::DSisoRational)            = tf([zero(Int8)], [one(Int8)], NaN64)
# one(::Type{DSisoRational})        = tf([one(Int8)], [one(Int8)], NaN64)
# one(s::DSisoRational)             = tf([one(Int8)], [one(Int8)], NaN64)
#
# # interface implementation
#
# zeros(s::DSisoRational)           = convert(Vector{Complex{Float64}}, roots(s.num))
# poles(s::DSisoRational)           = convert(Vector{Complex{Float64}}, roots(s.den))
# numvec(s::DSisoRational)          = reverse(coeffs(s.num))
# denvec(s::DSisoRational)          = reverse(coeffs(s.den))
# numpoly(s::DSisoRational)         = copy(s.num)
# denpoly(s::DSisoRational)         = copy(s.den)
# zpkdata(s::DSisoRational)         = (zeros(s), poles(s), s.num[end]/s.den[end])
# samplingtime(s::DSisoRational)    = s.Ts
# isdiscrete(s::DSisoRational)      = true
# isdiscrete(::Type{DSisoRational}) = true
#
# # overload printing functions
#
# function show(io::IO, s::DSisoRational)
#   println(io, "Discrete time rational transfer function model")
#   println(io, "\ty = Gu")
#   if s.Ts > 0
#     println(io, "with Ts=", s.Ts, ".")
#   elseif s.Ts == 0
#     println(io, "with Ts=unspecified.")
#   end
# end
#
# function showall(io::IO, s::DSisoRational)
#   show(io, s)
#   println(io, "")
#   printtransferfunction(io::IO, s)
# end
#
# # overload mathematical operations
#
# function +(s1::DSisoRational, s2::DSisoRational)
#   Ts::Float64
#   if s1.Ts == s2.Ts || isnan(s2.Ts)
#     Ts = s1.Ts
#   elseif isnan(s1.Ts)
#     Ts = s2.Ts
#   else
#     warn("Sampling time mismatch")
#     throw(DomainError())
#   end
#   den1,den2,dengcd   = rmgcd(s1.den, s2.den)
#   tf(s1.num*den2 + s2.num*den1, den1*den2*dengcd, Ts)
# end
# +{T<:Real}(s::DSisoRational, n::T)       = tf(s.num + n*s.den, copy(s.den), s.Ts)
# +{T<:Real}(n::T, s::DSisoRational)       = s + n
#
# .+{T<:Real}(s::DSisoRational, n::T)      = s + n
# .+{T<:Real}(n::T, s::DSisoRational)      = s + n
# .+(s1::DSisoRational, s2::DSisoRational) = +(s1,-s2)
#
# -{T}(s::DSisoRational{T})                = tf(-s.num, copy(s.den), s.Ts)
#
# -(s1::DSisoRational, s2::DSisoRational)  = +(s1,-s2)
# -{T<:Real}(n::T, s::DSisoRational)       = +(n, -s)
# -{T<:Real}(s::DSisoRational, n::T)       = +(s, -n)
#
# .-{T<:Real}(s::DSisoRational, n::T)      = -(s, n)
# .-{T<:Real}(n::T, s::DSisoRational)      = -(n, s)
# .-(s1::DSisoRational, s2::DSisoRational) = -(s1, s2)
#
# function *(s1::DSisoRational, s2::DSisoRational)
#   Ts::Float64
#   if s1.Ts == s2.Ts || isnan(s2.Ts)
#     Ts = s1.Ts
#   elseif isnan(s1.Ts)
#     Ts = s2.Ts
#   else
#     warn("Sampling time mismatch")
#     throw(DomainError())
#   end
#   num1,den2,gcd1 = rmgcd(s1.num, s2.den)
#   den1,num2,gcd2 = rmgcd(s1.den, s2.num)
#   tf(num1*num2, den1*den2, Ts)
# end
# *{T<:Real}(s::DSisoRational, n::T)       = tf(s.num*n, copy(s.den), s.Ts)
# *{T<:Real}(n::T, s::DSisoRational)       = *(s, n)
#
# .*{T<:Real}(s::DSisoRational, n::T)      = *(s, n)
# .*{T<:Real}(n::T, s::DSisoRational)      = *(n, s)
# .*(s1::DSisoRational, s2::DSisoRational) = *(s1, s2)
#
# /(s1::DSisoRational, s2::DSisoRational)  = s1*(1/s2)
# /{T<:Real}(n::T, s::DSisoRational)       = tf(n*s.den, copy(s.num), s.Ts)
# /{T<:Real}(s::DSisoRational, n::T)       = s*(1/n)
#
# ./{T<:Real}(n::T, s::DSisoRational)      = /(n, s)
# ./{T<:Real}(s::DSisoRational, n::T)      = /(s, n)
# ./(s1::DSisoRational, s2::DSisoRational) = /(s1, s2)
#
# function ==(s1::DSisoRational, s2::DSisoRational)
#   s1.Ts == s2.Ts && s1.num == s2.num &&
#     (s1.den == s2.den || s1.num == zero(s1.num))
#     # TODO scaling of num and den
# end
#
# !=(s1::DSisoRational, s2::DSisoRational) = !(s1 == s2)
#
# function isapprox(s1::DSisoRational, s2::DSisoRational,
#     rtol::Real=sqrt(eps()), atol::Real=0)
#   sdiff = s2-s1
#   return norm(sdiff.num) < rtol
# end
