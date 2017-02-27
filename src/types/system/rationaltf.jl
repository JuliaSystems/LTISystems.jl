immutable RationalTF{T,S,M} <: LtiSystem{T,S}
  mat::M
  nu::Int
  ny::Int
  Ts::Float64

  # Continuous-time, single-input-single-output rational transfer function model
  @compat function (::Type{RationalTF}){S,U<:Real,V<:Real}(r::RationalFunction{
    Var{:s},Conj{S},U,V})
    mat     = fill(r, 1, 1)
    ny, nu  = tfcheck(mat)
    new{Siso{true},Continuous{true},typeof(mat)}(mat, nu, ny, zero(Float64))
  end

  # Discrete-time, single-input-single-output rational transfer function model
  @compat function (::Type{RationalTF}){S,U<:Real,V<:Real}(r::RationalFunction{
    Var{:z},Conj{S},U,V}, Ts::Real)
    mat     = fill(r, 1, 1)
    ny, nu  = tfcheck(mat, Ts)
    new{Siso{true},Continuous{false},typeof(mat)}(mat, nu, ny, convert(Float64, Ts))
  end

  # Continuous-time, multi-input-multi-output rational transfer function model
  @compat function (::Type{RationalTF}){S,U<:Real,V<:Real}(
    mat::AbstractMatrix{RationalFunction{Var{:s},Conj{S},U,V}})
    ny, nu  = tfcheck(mat)
    new{Siso{false},Continuous{true},typeof(mat)}(mat, nu, ny, zero(Float64))
  end

  # Discrete-time, multi-input-multi-output rational transfer function model
  @compat function (::Type{RationalTF}){S,U<:Real,V<:Real}(
    mat::AbstractMatrix{RationalFunction{Var{:z},Conj{S},U,V}}, Ts::Real)
    ny, nu  = tfcheck(mat, Ts)
    new{Siso{false},Continuous{false},typeof(mat)}(mat, nu, ny, convert(Float64, Ts))
  end
end

@compat function (::Type{RationalTF})(r::RationalFunction)
  warn("RationalTF(r): r can only be a real-coefficient `RationalFunction` of variable `:s`")
  throw(DomainError())
end

@compat function (::Type{RationalTF})(r::RationalFunction, Ts::Real)
  warn("RationalTF(r, Ts): r can only be a real-coefficient `RationalFunction` of variable `:z`")
  throw(DomainError())
end

@compat function (::Type{RationalTF})(m::AbstractMatrix)
  warn("RationalTF(m): m can only be an `AbstractMatrix` of real-coefficient `RationalFunction` objects of variable `:s`")
  throw(DomainError())
end

@compat function (::Type{RationalTF})(m::AbstractMatrix, Ts::Real)
  warn("RationalTF(m, Ts): m can only be an `AbstractMatrix` of real-coefficient `RationalFunction` objects of variable `:z`")
  throw(DomainError())
end

# Enforce rational transfer function type invariance
function tfcheck{T,S,U<:Real,V<:Real}(mat::AbstractMatrix{RationalFunction{Var{T},
  Conj{S},U,V}}, Ts::Real = zero(Float64))
  # Check sampling time
  if Ts < zero(Ts) || isinf(Ts)
    warn("RationalTF: Ts must be a non-negative real number")
    throw(DomainError())
  end

  # Check input-output dimensions
  ny, nu = size(mat)

  if ny == 0 || nu == 0
    warn("RationalTF: `min(nu, ny) = 0`")
    throw(DomainError())
  end

  for idx in eachindex(mat)
    numpoly   = num(mat[idx])
    denpoly   = den(mat[idx])
    col, row  = divrem(idx-1, ny)

    if degree(numpoly) > degree(denpoly)
      warn("RationalTF: mat[$(row+1),$(col+1)] is not proper")
      throw(DomainError())
    end
    if denpoly == zero(denpoly)
      warn("RationalTF: mat[$(row+1),$(col+1)] has a zero denominator")
      throw(DomainError())
    end
  end

  return ny, nu
end

# Helper function
_reverse(x)         = reverse(x)
_reverse(x::Number) = [x]
typealias NV Union{Number,Vector}

# Outer constructors
## Continuous-time, single-input-single-output rational transfer function model
tf(r::RationalFunction) = RationalTF(r)
tf(num::NV, den::NV)    = RationalTF(RationalFunction(_reverse(num), _reverse(den), :s))
tf(num::NV, den::Poly)  = RationalTF(RationalFunction(_reverse(num), den))
tf(num::Poly, den::NV)  = RationalTF(RationalFunction(num, _reverse(den)))
tf(num::Poly, den::Poly)= RationalTF(RationalFunction(num, den))

## Discrete-time, single-input-single-output rational transfer function model
tf(r::RationalFunction, Ts::Real) = RationalTF(r, Ts)
tf(num::NV, den::NV, Ts::Real)    = RationalTF(RationalFunction(_reverse(num), _reverse(den), :z), Ts)
tf(num::NV, den::Poly, Ts::Real)  = RationalTF(RationalFunction(_reverse(num), den), Ts)
tf(num::Poly, den::NV, Ts::Real)  = RationalTF(RationalFunction(num, _reverse(den)), Ts)
tf(num::Poly, den::Poly, Ts::Real)= RationalTF(RationalFunction(num, den), Ts)

function tf(num::AbstractVector, den::AbstractVector, Ts::Real, var::Symbol)
  vars    = [:z̄,:q̄,:qinv,:zinv]
  if var ∉ vars
    warn("tf: $(var) ∉ ", vars)
    throw(DomainError())
  end

  numlast         = findlast(num)
  denlast         = findlast(den)
  order           = max(numlast, denlast)
  num_            = zeros(eltype(num), order)
  num_[1:numlast] = num[1:numlast]
  den_            = zeros(eltype(den), order)
  den_[1:denlast] = den[1:denlast]

  numpoly = Poly(reverse(num_), :z)
  denpoly = Poly(reverse(den_), :z)
  RationalTF(RationalFunction(numpoly, denpoly), Ts)
end

## Continuous-time, multi-input-multi-output rational transfer function model
tf(mat::AbstractMatrix)            = RationalTF(mat)

## Discrete-time, multi-input-multi-output rational transfer function model
tf(mat::AbstractMatrix, Ts::Real)  = RationalTF(mat, Ts)

# # Interfaces
# samplingtime(s::RationalTF) = s.Ts
#
# # Think carefully about how to implement numstates
# numstates(s::RationalTF)    = numstates(ss(s))
# # Currently, we only allow for proper systems
# numstates(s::RationalTF{Siso{true}}) = degree(den[1])
#
# numinputs(s::RationalTF)    = s.nu
# numoutputs(s::RationalTF)   = s.ny
#
# # Dimension information
# ndims(s::RationalTF{Siso{true}})  = 1
# ndims(s::RationalTF{Siso{false}}) = 2
# size(s::RationalTF)               = size(s.num)
# size(s::RationalTF, dim::Int)     = size(s.num, dim)
# size(s::RationalTF, dims::Int...) = size(s.num, dims)
#
# # Iteration interface (meaningful only for MIMO systems)
# # TODO
#
# # Slicing (`getindex`) of MIMO systems
# # TODO
#
# # Printing functions
# summary(s::RationalTF{Siso{true},Continuous{true}})   =
#   string("tf(nx=", numstates(s), ")")
# summary(s::RationalTF{Siso{true},Continuous{false}})  =
#   string("tf(nx=", numstates(s), ",Ts=", s.Ts, ")")
# summary(s::RationalTF{Siso{false},Continuous{true}})  =
#   string("tf(nx=", numstates(s), ",nu=", s.nu, ",ny=", s.ny, ")")
# summary(s::RationalTF{Siso{false},Continuous{false}}) =
#   string("tf(nx=", numstates(s), ",nu=", s.nu, ",ny=", s.ny, ",Ts=", s.Ts, ")")
#
# showcompact(io::IO, s::RationalTF) = print(io, summary(s))
#
# function show{T}(io::IO, s::RationalTF{T,Continuous{true}})
#   # TODO
# end
#
# function show{T}(io::IO, s::RationalTF{T,Continuous{false}})
#   # TODO
# end
#
# function showall(io::IO, s::RationalTF)
#   # TODO
# end
#
# # Conversion and promotion
# promote_rule{T<:Real,S}(::Type{T}, ::Type{RationalTF{Siso{true},S}}) =
#   RationalTF{Siso{true},S}
# promote_rule{T<:AbstractMatrix,S}(::Type{T}, ::Type{RationalTF{Siso{false},S}}) =
#   RationalTF{Siso{false},S}
#
# convert(::Type{RationalTF{Siso{true},Continuous{true}}}, g::Real)             =
#   tf(g)
# convert(::Type{RationalTF{Siso{true},Continuous{false}}}, g::Real)            =
#   tf(g, zero(Float64))
# convert(::Type{RationalTF{Siso{false},Continuous{true}}}, g::AbstractMatrix)  =
#   tf(g)
# convert(::Type{RationalTF{Siso{false},Continuous{false}}}, g::AbstractMatrix) =
#   tf(g, zero(Float64))
#
# # Multiplicative and additive identities (meaningful only for SISO)
# one(::Type{RationalTF{Siso{true},Continuous{true}}})    =
#   tf(one(Int8))
# one(::Type{RationalTF{Siso{true},Continuous{false}}})   =
#   tf(one(Int8), zero(Float64))
# zero(::Type{RationalTF{Siso{true},Continuous{true}}})   =
#   tf(zero(Int8))
# zero(::Type{RationalTF{Siso{true},Continuous{false}}})  =
#   tf(zero(Int8), zero(Float64))
#
# one(s::RationalTF{Siso{true},Continuous{true}})   =
#   tf(one(eltype(s.num)), one(eltype(s.den)))
# one(s::RationalTF{Siso{true},Continuous{false}})  =
#   tf(one(eltype(s.num)), one(eltype(s.den)), zero(Float64))
# zero(s::RationalTF{Siso{true},Continuous{true}})  =
#   tf(zero(eltype(s.num)), one(eltype(s.den)))
# zero(s::RationalTF{Siso{true},Continuous{false}}) =
#   tf(zero(eltype(s.num)), one(eltype(s.den)), zero(Float64))
#
# # Inverse of a rational transfer function model
# function tfinv(s::RationalTF{Siso{true}})
#   @assert degree(num[1]) == degree(den[1]) "inv(sys): inverse is not a proper system"
#   return den, num
# end
#
# function tfinv(s::RationalTF{Siso{false}})
#   @assert s.ny == s.nu "inv(sys): s.ny ≠ s.nu"
#   # TODO
#   # return num, den
# end
#
# function inv(s::RationalTF{Siso{true},Continuous{true}})
#   num, den = tfinv(s)
#   tf(num[1], den[1])
# end
#
# function inv(s::RationalTF{Siso{true},Continuous{false}})
#   num, den = tfinv(s)
#   tf(num[1], den[1], s.Ts)
# end
#
# function inv(s::RationalTF{Siso{false},Continuous{true}})
#   num, den = tfinv(s)
#   tf(num, den)
# end
#
# function inv(s::RationalTF{Siso{false},Continuous{false}})
#   num, den = tfinv(s)
#   tf(num, den, s.Ts)
# end
#
# # Invariant zeros of a rational transfer function model
# zeros(s::RationalTF{Siso{true}}) = roots(s.num[1])
#
# function zeros(s::RationalTF)
#   # TODO
# end
#
# # Transmission zeros of a rational transfer function model
# tzeros(s::RationalTF) = zeros(minreal(s))
#
# # Poles of a rational transfer function model
# poles(s::RationalTF{Siso{true}}) = roots(s.den[1])
#
# function poles(s::RationalTF)
#   # TODO
# end
#
# # Negative of a rational transfer function model
# -(s::RationalTF{Siso{true},Continuous{true}})   =
#   RationalTF(-s.num[1], s.den[1])
# -(s::RationalTF{Siso{true},Continuous{false}})  =
#   RationalTF(-s.num[1], s.den[1], s.Ts)
# -(s::RationalTF{Siso{false},Continuous{true}})  =
#   RationalTF(-s.num, s.den)
# -(s::RationalTF{Siso{false},Continuous{false}}) =
#   RationalTF(-s.num, s.den, s.Ts)
#
# # Addition
# function tfparallel{T1,T2,S}(s1::RationalTF{T1,S}, s2::RationalTF{T2,S})
#   @assert s1.Ts ≈ s2.Ts || s1.Ts == zero(Float64) ||
#     s2.Ts == zero(Float64) "parallel(s1,s2): Sampling time mismatch"
#   @assert size(s1) == size(s2) "parallel(s1,s2): size(s1) ≠ size(s2)"
#
#   # return num, den
# end
#
# function +(s1::RationalTF{Siso{true},Continuous{true}},
#   s2::RationalTF{Siso{true},Continuous{true}})
#   num, den = tfparallel(s1, s2)
#   RationalTF(num[1], den[1])
# end
#
# function +(s1::RationalTF{Siso{true},Continuous{false}},
#   s2::RationalTF{Siso{true},Continuous{false}})
#   num, den = tfparallel(s1, s2)
#   RationalTF(num[1], den[1], s1.Ts)
# end
#
# function +{T1,T2}(s1::RationalTF{T1,Continuous{true}},
#   s2::RationalTF{T2,Continuous{true}})
#   num, den = tfparallel(s1, s2)
#   RationalTF(num, den)
# end
#
# function +{T1,T2}(s1::RationalTF{T1,Continuous{false}},
#   s2::RationalTF{T2,Continuous{false}})
#   num, den = tfparallel(s1, s2)
#   RationalTF(num, den, s1.Ts)
# end
#
# .+(s1::RationalTF{Siso{true}}, s2::RationalTF{Siso{true}}) = +(s1, s2)
#
# +{T}(s::RationalTF{T,Continuous{true}}, g)  = +(s, tf(g))
# +{T}(s::RationalTF{T,Continuous{false}}, g) = +(s, tf(g, zero(Float64)))
# +{T}(g, s::RationalTF{T,Continuous{true}})  = +(tf(g), s)
# +{T}(g, s::RationalTF{T,Continuous{false}}) = +(tf(g, zero(Float64)), s)
#
# .+(s::RationalTF{Siso{true}}, g::Real)  = +(s, g)
# .+(g::Real, s::RationalTF{Siso{true}})  = +(g, s)
#
# # Subtraction
# -(s1::RationalTF, s2::RationalTF) = +(s1, -s2)
#
# .-(s1::RationalTF{Siso{true}}, s2::RationalTF{Siso{true}}) = -(s1, s2)
#
# -{T}(s::RationalTF{T,Continuous{true}}, g)  = -(s, tf(g))
# -{T}(s::RationalTF{T,Continuous{false}}, g) = -(s, tf(g, zero(Float64)))
# -{T}(g, s::RationalTF{T,Continuous{true}})  = -(tf(g), s)
# -{T}(g, s::RationalTF{T,Continuous{false}}) = -(tf(g, zero(Float64)), s)
#
# .-(s::RationalTF{Siso{true}}, g::Real)  = -(s, g)
# .-(g::Real, s::RationalTF{Siso{true}})  = -(g, s)
#
# # Multiplication
# function tfseries{T1,T2,S}(s1::RationalTF{T1,S}, s2::RationalTF{T2,S})
#   # Remark: s1*s2 implies u -> s2 -> s1 -> y
#   @assert s1.Ts ≈ s2.Ts || s1.Ts == zero(Float64) ||
#     s2.Ts == zero(Float64) "series(s1,s2): Sampling time mismatch"
#   @assert s1.nu == s2.ny "series(s1,s2): s1.nu ≠ s2.ny"
#
#   return num, den
# end
#
# function *(s1::RationalTF{Siso{true},Continuous{true}},
#   s2::RationalTF{Siso{true},Continuous{true}})
#   num, den = tfseries(s1, s2)
#   RationalTF(num[1], den[1])
# end
#
# function *(s1::RationalTF{Siso{true},Continuous{false}},
#   s2::RationalTF{Siso{true},Continuous{false}})
#   num, den = tfseries(s1, s2)
#   RationalTF(num[1], den[1], s1.Ts)
# end
#
# function *{T1,T2}(s1::RationalTF{T1,Continuous{true}},
#   s2::RationalTF{T2,Continuous{true}})
#   num, den = tfseries(s1, s2)
#   RationalTF(num, den)
# end
#
# function *{T1,T2}(s1::RationalTF{T1,Continuous{false}},
#   s2::RationalTF{T2,Continuous{false}})
#   num, den = tfseries(s1, s2)
#   RationalTF(num, den, s1.Ts)
# end
#
# .*(s1::RationalTF{Siso{true}}, s2::RationalTF{Siso{true}}) = *(s1, s2)
#
# *{T}(s::RationalTF{T,Continuous{true}}, g)  = *(s, tf(g))
# *{T}(s::RationalTF{T,Continuous{false}}, g) = *(s, tf(g, zero(Float64)))
# *{T}(g, s::RationalTF{T,Continuous{true}})  = *(tf(g), s)
# *{T}(g, s::RationalTF{T,Continuous{false}}) = *(tf(g, zero(Float64)), s)
#
# .*(s::RationalTF{Siso{true}}, g::Real)  = *(s, g)
# .*(g::Real, s::RationalTF{Siso{true}})  = *(g, s)
#
# # Division
# /(s1::RationalTF, s2::RationalTF) = *(s1, inv(s2))
#
# ./(s1::RationalTF{Siso{true}}, s2::RationalTF{Siso{true}}) = /(s1, s2)
#
# /{T}(s::RationalTF{T,Continuous{true}}, g)  = /(s, tf(g))
# /{T}(s::RationalTF{T,Continuous{false}}, g) = /(s, tf(g, zero(Float64)))
# /{T}(g, s::RationalTF{T,Continuous{true}})  = /(tf(g), s)
# /{T}(g, s::RationalTF{T,Continuous{false}}) = /(tf(g, zero(Float64)), s)
#
# ./(s::RationalTF{Siso{true}}, g::Real)  = /(s, g)
# ./(g::Real, s::RationalTF{Siso{true}})  = /(g, s)
