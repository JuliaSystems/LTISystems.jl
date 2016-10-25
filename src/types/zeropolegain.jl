immutable ZeroPoleGain{T,S,M1,M2,M3} <: LtiSystem{T,S}
  Z::M1
  P::M2
  K::M3
  Ts::Float64

  # Continuous-time, single-input-single-output zero-pole-gain model
  @compat function (::Type{ZeroPoleGain}){M1<:AbstractVector,M2<:AbstractVector,
    M3<:Real}(z::M1, p::M2, k::M3)
    @assert eltype(z) <: Number   "ZeroPoleGain: z must be a vector of numbers"
    @assert eltype(p) <: Number   "ZeroPoleGain: p must be a vector of numbers"
    @assert length(z) ≤ length(p) "ZeroPoleGain: system is not proper"

    new{Siso{true},Continuous{true},M1,M2,M3}(z, p, k, zero(Float64))
  end

  # Discrete-time, single-input-single-output zero-pole-gain model
  @compat function (::Type{ZeroPoleGain}){M1<:AbstractVector,M2<:AbstractVector,
    M3<:Real, M4<:Real}(z::M1, p::M2, k::M3, Ts::M4)
    @assert eltype(z) <: Number         "ZeroPoleGain: z must be a vector of numbers"
    @assert eltype(p) <: Number         "ZeroPoleGain: p must be a vector of numbers"
    @assert length(z) ≤ length(p)       "ZeroPoleGain: system is not proper"
    @assert Ts ≥ zero(Ts) && !isinf(Ts) "ZeroPoleGain: Ts must be a non-negative number"

    new{Siso{true},Continuous{false},M1,M2,M3}(z, p, k, convert(Float64, Ts))
  end

  # Continuous-time, multi-input-multi-output zero-pole-gain model
  @compat function (::Type{ZeroPoleGain}){M1<:AbstractMatrix,M2<:AbstractMatrix,
    M3<:AbstractMatrix}(z::M1, p::M2, k::M3)
    @assert size(z) == size(p) == size(k) "ZeroPoleGain: dimensions of z, p and k must match"
    @assert eltype(z) <: AbstractVector &&
      eltype(eltype(z)) <: Number "ZeroPoleGain: z must be a matrix of vectors of numbers"
    @assert eltype(p) <: AbstractVector &&
      eltype(eltype(p)) <: Number "ZeroPoleGain: p must be a matrix of vectors of numbers"
    for idx in eachindex(z)
      @assert length(z[idx]) ≤ length(p[idx]) "ZeroPoleGain: system is not proper"
    end
    @assert eltype(k) <: Real "ZeroPoleGain: k must be a matrix of reals"

    new{Siso{false},Continuous{true},M1,M2,M3}(z, p, k, zero(Float64))
  end

  # Discrete-time, multi-input-multi-output zero-pole-gain model
  @compat function (::Type{ZeroPoleGain}){M1<:AbstractMatrix,M2<:AbstractMatrix,
    M3<:AbstractMatrix}(z::M1, p::M2, k::M3, Ts::Real)
    @assert size(z) == size(p) == size(k) "ZeroPoleGain: dimensions of z, p and k must match"
    @assert eltype(z) <: AbstractVector &&
      eltype(eltype(z)) <: Number "ZeroPoleGain: z must be a matrix of vectors of numbers"
    @assert eltype(p) <: AbstractVector &&
      eltype(eltype(p)) <: Number "ZeroPoleGain: p must be a matrix of vectors of numbers"
    for idx in eachindex(z)
      @assert length(z[idx]) ≤ length(p[idx]) "ZeroPoleGain: system is not proper"
    end
    @assert eltype(k) <: Real           "ZeroPoleGain: k must be a matrix of reals"
    @assert Ts ≥ zero(Ts) && !isinf(Ts) "ZeroPoleGain: Ts must be a non-negative number"

    new{Siso{false},Continuous{false},M1,M2,M3}(z, p, k, convert(Float64, Ts))
  end
end

# Legacy CSisoZpk
# immutable CSisoZpk{T1<:Real, T2<:Real, T3<:Real} <: CSisoTf
#   z::Vector{Complex{T1}}
#   p::Vector{Complex{T2}}
#   k::T3
#
#   function call{T1<:Number, T2<:Number, T3<:Real}(::Type{CSisoZpk},
#       z::Vector{T1}, p::Vector{T2}, k::T3)
#
#     new{real(T1),real(T2),T3}(z, p, k)
#   end
# end
#
# # creation of continuous zpk transfer functions
#
# """
#     zpk([z, p,] k[, Ts])
#
# Constructs a continuous transfer function with zeros `z`,
# poles `p` of type `Vector` and gain `k.
#
# A gain transfer function without zeros and poles is constructed by only
#   supplying `k`.
#
# A continuous MIMO transfer function can be constructed from a `Matrix` of `Vector`s.
#
# Discrete transfer functions are constructed using additional argument sampling
# time `Ts`.
#
# # Examples
# ```julia
# julia> zpk([1,0,3],[1,1,2],2)
# ZPK:
#     (s - 1.0)s(s - 3.0)
# 2.0 --------------------
#     (s - 1.0)^2(s - 2.0)
#
# Continuous-time zero-pole-gain model
# ```
#
# ```julia
# julia> zpk(2)
# ZPK:
#     1.0
# 2.0 ---
#     1.0
# ```
#
# ```julia
# julia> zpk([1,0,3],[1,1,2],2,1)
# ZPK:
#     (z - 1.0)z(z - 3.0)
# 2.0 --------------------
#     (z - 1.0)^2(z - 2.0)
#
#
# Sample Time: 1 (seconds)
# Discrete-time zero-pole-gain model
# ```
# """
# zpk{T1<:Number, T2<:Number, T3<:Real}(z::Vector{T1}, p::Vector{T2}, k::T3) = CSisoZpk(z, p, k)
# zpk{T1<:Real}(k::T1) = CSisoZpk(Vector{Int8}(), Vector{Int8}(), k)
#
# # conversion and promotion
#
# promote_rule{T11,T12,T13,T21,T22,T23}(::Type{CSisoZpk{T11,T12,T13}},
#   ::Type{CSisoZpk{T21,T22,T23}})                                 = CSisoZpk
# promote_rule{T11,T12,T13,T5<:Real}(::Type{CSisoZpk{T11,T12,T13}}, ::Type{T5}) =
#   CSisoZpk
# promote_rule{T<:Real}(::Type{CSisoZpk}, ::Type{T})               = CSisoZpk
#
# convert{T11,T12,T13}(::Type{CSisoZpk}, x::CSisoZpk{T11,T12,T13}) = x
# convert{T<:Real}(::Type{CSisoZpk}, x::T)                         = zpk(x)
#
# # overloading identities
#
# zero(::Type{CSisoZpk})       = zpk(zero(Int8))
# zero(s::CSisoZpk)            = zpk(zero(Int8))
# one(::Type{CSisoZpk})        = zpk(one(Int8))
# one(s::CSisoZpk)             = zpk(one(Int8))
#
# # interface implementation
#
# zeros(s::CSisoZpk)           = copy(s.z)
# poles(s::CSisoZpk)           = copy(s.p)
# numvec(s::CSisoZpk)          = reverse(coeffs(numpoly(s)))
# denvec(s::CSisoZpk)          = reverse(coeffs(denpoly(s)))
# numpoly(s::CSisoZpk)         = convert(Poly{real(eltype(s.z))}, poly(s.z))
# denpoly(s::CSisoZpk)         = convert(Poly{real(eltype(s.p))}, poly(s.p))
# zpkdata(s::CSisoZpk)         = (s.z, s.p, s.k)
# samplingtime(s::CSisoZpk)    = zero(Float64)
# isdiscrete(s::CSisoZpk)      = false
# isdiscrete(::Type{CSisoZpk}) = false
#
# # overload printing functions
#
# function show(io::IO, s::CSisoZpk)
#   println(io, "Continuous time zpk transfer function model")
#   println(io, "\ty = Gu")
# end
#
# function showall(io::IO, s::CSisoZpk)
#   show(io, s)
#   println(io, "")
#   printtransferfunction(io::IO, s)
# end
#
# # overload mathematical operations
#
# function +(s1::CSisoZpk, s2::CSisoZpk)
#   p1,p2,pcommon = rmcommon(copy(s1.p), copy(s2.p))
#   z1,z2,zcommon = rmcommon(copy(s1.z), copy(s2.z))
#   Z = s1.k*poly(z1)*poly(p2) + s2.k*poly(z2)*poly(p1)
#   z = vcat(convert(Vector{Complex{Float64}},roots(Z)), zcommon)
#   p = vcat(p1, p2, pcommon)
#   k = real(Z[end]) # Poly is now reverse order
#   zpk(z, p, k)
# end
#
# function +{T<:Real}(s::CSisoZpk, n::T)
#   Z = s.k*poly(s.z) + n*poly(s.p)
#   z::Array{Complex{Float64}} = roots(Z)
#   p = s.p
#   k = real(Z[end]) # Poly is now reverse order
#   CSisoZpk(z, p, k)
# end
# +{T<:Real}(n::T, s::CSisoZpk)  = +(s, n)
#
# .+{T<:Real}(s::CSisoZpk, n::T) = +(s, n)
# .+{T<:Real}(n::T, s::CSisoZpk) = +(n, s)
# .+(s1::CSisoZpk, s2::CSisoZpk) = +(s1, s2)
#
# -(s::CSisoZpk)                 = zpk(copy(s.z), copy(s.p), -s.k)
#
# -(s1::CSisoZpk, s2::CSisoZpk)  = +(s1, -s2)
# -{T<:Real}(s::CSisoZpk, n::T)  = +(s, -n)
# -{T<:Real}(n::T, s::CSisoZpk)  = +(n, -s)
#
# .-{T<:Real}(s::CSisoZpk, n::T) = +(s, -n)
# .-{T<:Real}(n::T, s::CSisoZpk) = +(n, -s)
# .-(s1::CSisoZpk, s2::CSisoZpk) = +(s1, -s2)
#
# function *(s1::CSisoZpk, s2::CSisoZpk)
#   z1,p2,pcommon = rmcommon(s1.z, s2.p)
#   p1,z2,zcommon = rmcommon(s1.p, s2.z)
#   z = vcat(z1, z2)
#   p = vcat(p1, p2)
#   k = s1.k*s2.k
#   zpk(z, p, k)
# end
# *{T<:Real}(s::CSisoZpk, n::T)  = zpk(copy(s.z), copy(s.p), n*s.k)
# *{T<:Real}(n::T, s::CSisoZpk)  = *(s, n)
#
# .*{T<:Real}(s::CSisoZpk, n::T) = *(n, s)
# .*{T<:Real}(n::T, s::CSisoZpk) = *(s, n)
# .*(s1::CSisoZpk, s2::CSisoZpk) = *(s1,s2)
#
# /{T<:Real}(n::T, s::CSisoZpk)  = zpk(copy(s.p), copy(s.z), n./s.k)
# /{T<:Real}(s::CSisoZpk, n::T)  = s*(1/n)
# /(s1::CSisoZpk, s2::CSisoZpk)  = s1*(1/s2)
#
# ./{T<:Real}(n::T, s::CSisoZpk) = n*(1/s)
# ./{T<:Real}(s::CSisoZpk, n::T) = s*(1/n)
# ./(s1::CSisoZpk, s2::CSisoZpk) = s1*(1/s2)
#
# function ==(s1::CSisoZpk, s2::CSisoZpk)
#   fields = [:z, :p, :k]
#   for field in fields
#     if getfield(s1, field) != getfield(s2, field)
#       return false
#     end
#   end
#   true
# end
#
# !=(s1::CSisoZpk, s2::CSisoZpk) = !(s1==s2)
#
# function isapprox(s1::CSisoZpk, s2::CSisoZpk,
#     rtol::Real=sqrt(eps()), atol::Real=0)
#   sdiff = s2-s1
#   norm(sdiff.k) < rtol
# end

# Legacy DSisoZpk
# immutable DSisoZpk{T1<:Real, T2<:Real, T3<:Real} <: DSisoTf
#   z::Vector{Complex{T1}}
#   p::Vector{Complex{T2}}
#   k::T3
#   Ts::Float64
#
#   function call{T1<:Number, T2<:Number, T3<:Real}(::Type{DSisoZpk},
#     z::Vector{T1}, p::Vector{T2}, k::T3, Ts::Float64)
#
#     Ts_ = Ts > zero(Float64) ? Ts : NaN
#     new{real(T1),real(T2),T3}(z, p, k, Ts_)
#   end
# end
#
# # creation of discrete zpk transfer functions
#
# zpk{T1<:Number, T2<:Number, T3<:Real, T4<:Real}(z::Vector{T1}, p::Vector{T2},
#   k::T3, Ts::T4)                       = DSisoZpk(z, p, k, Float64(Ts))
# zpk{T1<:Real, T2<:Real}(k::T1, Ts::T2) = DSisoZpk(Vector{Int8}(), Vector{Int8}(), k, Float64(Ts))
#
# # conversion and promotion
#
# promote_rule{T11,T12,T13,T5<:Real}(::Type{DSisoZpk{T11,T12,T13}},
#   ::Type{T5})                                                    = DSisoZpk
# promote_rule{T11,T12,T13,T21,T22,T23}(::Type{DSisoZpk{T11,T12,T13}},
#   ::Type{DSisoZpk{T21,T22,T23}})                                 = DSisoZpk
# promote_rule{T<:Real}(::Type{DSisoZpk}, ::Type{T})               = DSisoZpk
#
# convert{T11,T12,T13}(::Type{DSisoZpk}, x::DSisoZpk{T11,T12,T13}) = x
# convert{T<:Real}(::Type{DSisoZpk}, x::T)                         = zpk(x, NaN64)
#
# # overloading identities
#
#
# zero(::Type{DSisoZpk})       = zpk(zero(Int8), NaN64)
# zero(s::DSisoZpk)            = zpk(zero(Int8), NaN64)
# one(::Type{DSisoZpk})        = zpk(one(Int8), NaN64)
# one(s::DSisoZpk)             = zpk(one(Int8), NaN64)
#
# # interface implementation
#
# zeros(s::DSisoZpk)           = copy(s.z)
# poles(s::DSisoZpk)           = copy(s.p)
# numvec(s::DSisoZpk)          = reverse(coeffs(numpoly(s)))
# denvec(s::DSisoZpk)          = reverse(coeffs(denpoly(s)))
# numpoly(s::DSisoZpk)         = convert(Poly{real(eltype(s.z))}, poly(s.z))
# denpoly(s::DSisoZpk)         = convert(Poly{real(eltype(s.p))}, poly(s.p))
# zpkdata(s::DSisoZpk)         = (s.z, s.p, s.k)
# samplingtime(s::DSisoZpk)    = copy(s.Ts)
# isdiscrete(s::DSisoZpk)      = true
# isdiscrete(::Type{DSisoZpk}) = true
#
# # overload printing functions
#
# function show(io::IO, s::DSisoZpk)
#   println(io, "Discrete time zpk transfer function model")
#   println(io, "\ty = Gu")
#   if s.Ts > 0
#     println(io, "with Ts=", s.Ts, ".")
#   elseif s.Ts == 0
#     println(io, "with Ts=unspecified.")
#   end
# end
#
# function showall(io::IO, s::DSisoZpk)
#   show(io, s)
#   println(io, "")
#   printtransferfunction(io::IO, s)
# end
#
# # overload mathematical operations
#
# function +{T11,T12,T13,T21,T22,T23}(s1::DSisoZpk{T11,T12,T13}, s2::DSisoZpk{T21,T22,T23})
#   Ts::Float64
#   if s1.Ts == s2.Ts || isnan(s2.Ts)
#     Ts = s1.Ts
#   elseif isnan(s1.Ts)
#     Ts = s2.Ts
#   else
#     warn("Sampling time mismatch")
#     throw(DomainError())
#   end
#   p1,p2,pcommon = rmcommon(copy(s1.p), copy(s2.p))
#   z1,z2,zcommon = rmcommon(copy(s1.z), copy(s2.z))
#   Z = s1.k*poly(z1)*poly(p2) + s2.k*poly(z2)*poly(p1)
#   z = vcat(convert(Vector{Complex{Float64}},roots(Z)), zcommon)
#   p = vcat(p1, p2, pcommon)
#   k = real(Z[end]) # Poly is now reverse order
#   zpk(z, p, k, Ts)
# end
#
# function +{T11,T12,T13,T<:Real}(s::DSisoZpk{T11,T12,T13}, n::T)
#   Tk = promote_type(T11, T)
#   Tz = float(promote_type(T11, T))
#   Z = s.k*poly(s.z) + n*poly(s.p)
#   z = roots(Z)
#   p = copy(s.p)
#   k = real(Z[end])
#   zpk(z, p, k, s.Ts)::DSisoZpk{Tz,T13,Tk}
# end
# +{T<:Real}(n::T, s::DSisoZpk)  = +(s, n)
#
# .+{T<:Real}(s::DSisoZpk, n::T) = +(s, n)
# .+{T<:Real}(n::T, s::DSisoZpk) = +(n, s)
# .+(s1::DSisoZpk, s2::DSisoZpk) = +(s1, s2)
#
# -(s::DSisoZpk)                 = zpk(copy(s.z), copy(s.p), -s.k, s.Ts)
#
# -(s1::DSisoZpk, s2::DSisoZpk)  = +(s1, -s2)
# -{T<:Real}(s::DSisoZpk, n::T)  = +(s, -n)
# -{T<:Real}(n::T, s::DSisoZpk)  = +(n, -s)
#
# .-{T<:Real}(s::DSisoZpk, n::T) = +(s, -n)
# .-{T<:Real}(n::T, s::DSisoZpk) = +(n, -s)
# .-(s1::DSisoZpk, s2::DSisoZpk) = +(s1, -s2)
#
# function *(s1::DSisoZpk, s2::DSisoZpk)
#   Ts::Float64
#   if s1.Ts == s2.Ts || isnan(s2.Ts)
#     Ts = s1.Ts
#   elseif isnan(s1.Ts)
#     s2.Ts
#   else
#     warn("Sampling time mismatch")
#     throw(DomainError())
#   end
#   z1,p2,pcommon = rmcommon(s1.z, s2.p)
#   p1,z2,zcommon = rmcommon(s1.p, s2.z)
#   z = vcat(z1, z2)
#   p = vcat(p1, p2)
#   k = s1.k*s2.k
#   zpk(z, p, k, Ts)
# end
# *{T<:Real}(s::DSisoZpk, n::T)  = zpk(copy(s.z), copy(s.p), n*s.k, s.Ts)
# *{T<:Real}(n::T, s::DSisoZpk)  = *(s, n)
#
# .*{T<:Real}(s::DSisoZpk, n::T) = *(n, s)
# .*{T<:Real}(n::T, s::DSisoZpk) = *(s, n)
# .*(s1::DSisoZpk, s2::DSisoZpk) = *(s1, s2)
#
# /{T<:Real}(n::T, s::DSisoZpk)  = zpk(copy(s.p), copy(s.z), n./s.k, s.Ts)
# /{T<:Real}(s::DSisoZpk, n::T)  = s*(1/n)
# /(s1::DSisoZpk, s2::DSisoZpk)  = s1*(1/s2)
#
# ./{T<:Real}(n::T, s::DSisoZpk) = n*(1/s)
# ./{T<:Real}(s::DSisoZpk, n::T) = s*(1/n)
# ./(s1::DSisoZpk, s2::DSisoZpk) = s1*(1/s2)
#
# function ==(s1::DSisoZpk, s2::DSisoZpk)
#   s1.Ts == s2.Ts && s1.z == s2.z &&
#     s1.p == s2.p && s1.k == s2.k
# end
#
# !=(s1::DSisoZpk, s2::DSisoZpk) = !(s1==s2)
#
# function isapprox(s1::DSisoZpk, s2::DSisoZpk,
#     rtol::Real=sqrt(eps()), atol::Real=0)
#   sdiff = s2-s1
#   norm(sdiff.k) < rtol
# end
