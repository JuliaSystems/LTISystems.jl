# Parameters:
#   T:  Val{:siso} or Val{:mimo}
#   S:  Val{:cont} or Val{:disc}
#   L:  Val{:lfd} or Val{:mfd}
#   M1: Type of numerator polynomial matrix
#   M2: Type of denominator polynomial matrix
#   N:  Numerator polynomial matrix
#   D:  Denominator polynomial matrix
#   nu: Number of inputs
#   ny: Number of outputs
#   Ts: Sampling time (= zero(Float64) for continuous-time systems)
#
# NOTE: Should SISO MFDs be based on 1x1 polynomial matrices (similarly to `RationalTF`)?
# NOTE: Should the constructors verify that D is nonsingular and that the MFD is proper?
immutable MFD{T,S,L,M1,M2}  <: LtiSystem{T,S}
  N::M1
  D::M2
  nu::Int
  ny::Int
  Ts::Float64

  # Continuous-time, single-input-single-output MFD model
  @compat function (::Type{MFD}){L,M1<:Polynomials.Poly,M2<:Polynomials.Poly}(N::M1, D::M2, ::Type{Val{L}})
    mfdcheck(N,D)
    new{Val{:siso},Val{:cont},Val{L},M1,M2}(N, D, 1, 1, zero(Float64))
  end

  # Discrete-time, single-input-single-output MFD model
  @compat function (::Type{MFD}){L,M1<:Polynomials.Poly,M2<:Polynomials.Poly}(N::M1, D::M2, Ts::Real, ::Type{Val{L}})
    mfdcheck(N,D,Ts)
    new{Val{:siso},Val{:disc},Val{L},M1,M2}(N, D, 1, 1, convert(Float64, Ts))
  end

  # Continuous-time, multi-input-multi-output MFD model
  @compat function (::Type{MFD}){L,M1<:PolynomialMatrices.PolyMatrix,
    M2<:PolynomialMatrices.PolyMatrix}(N::M1, D::M2, ::Type{Val{L}})
    ny, nu = mfdcheck(N, D, Val{L})
    new{Val{:mimo},Val{:cont},Val{L},M1,M2}(N, D, nu, ny, zero(Float64))
  end

  # Discrete-time, multi-input-multi-output MFD model
  @compat function (::Type{MFD}){L,M1<:PolynomialMatrices.PolyMatrix,
    M2<:PolynomialMatrices.PolyMatrix}(N::M1, D::M2, Ts::Real, ::Type{Val{L}})
    ny, nu = mfdcheck(N, D, Val{L}, Ts)
    new{Val{:mimo},Val{:disc},Val{L},M1,M2}(N, D, nu, ny, convert(Float64, Ts))
  end

end

function mfdcheck{T<:Real,S<:Real}(N::Poly{T}, D::Poly{S}, Ts::Real = zero(Float64))
  @assert Ts ≥ zero(Ts) && !isinf(Ts) "MFD: Ts must be non-negative real number"
end

# Enforce rational transfer function type invariance
function mfdcheck{M1<:PolynomialMatrices.PolyMatrix,M2<:PolynomialMatrices.PolyMatrix}(
  N::M1, D::M2, ::Type{Val{:lfd}})
  @assert size(N,1) == size(D,1) "MFD: size(N,1) ≠ size(D,1)"
  @assert size(D,1) == size(D,2) "MFD: size(D,1) ≠ size(D,2)"

  return size(N,1), size(N,2)
end
function mfdcheck{M1<:PolynomialMatrices.PolyMatrix,M2<:PolynomialMatrices.PolyMatrix}(
  N::M1, D::M2, ::Type{Val{:rfd}})
  @assert size(N,2) == size(D,2) "MFD: size(N,2) ≠ size(D,2)"
  @assert size(D,1) == size(D,2) "MFD: size(D,1) ≠ size(D,2)"

  return size(N,1), size(N,2)
end
function mfdcheck{T,M1<:PolynomialMatrices.PolyMatrix,M2<:PolynomialMatrices.PolyMatrix}(
  N::M1, D::M2, ::Type{Val{T}}, Ts::Real)
  @assert Ts ≥ zero(Ts) && !isinf(Ts) "MFD: Ts must be non-negative real number"
  return mfdcheck(N, D, Val{T})
end

# Outer constructors
lfd(N::Poly, D::Poly)           = MFD(N, D, Val{:lfd})
lfd(N::Poly, D::Poly, Ts::Real) = MFD(N, D, convert(Float64, Ts), Val{:lfd})
lfd{M1<:PolyMatrix,M2<:PolyMatrix}(N::M1, D::M2) =
  MFD(N, D, Val{:lfd})
lfd{M1<:PolyMatrix,M2<:PolyMatrix}(N::M1, D::M2, Ts::Real) =
  MFD(N, D, convert(Float64, Ts), Val{:lfd})

rfd(N::Poly, D::Poly) = MFD(N, D, Val{:rfd})
rfd(N::Poly, D::Poly, Ts::Real) = MFD(N, D, convert(Float64, Ts), Val{:rfd})
rfd(N::PolyMatrix, D::PolyMatrix) = MFD(N, D, Val{:rfd})
rfd(N::PolyMatrix, D::PolyMatrix, Ts::Real) = MFD(N, D, convert(Float64, Ts), Val{:rfd})

# Vector constructors
lfd{T1<:Real, T2<:Real}(N::AbstractVector{T1}, D::AbstractVector{T2}) =
  lfd(Poly(reverse(N), :s), Poly(reverse(D), :s))
lfd{T1<:Real, T2<:Real}(N::AbstractVector{T1}, D::AbstractVector{T2}, Ts::Real) =
  lfd(Poly(reverse(N), :s), Poly(reverse(D), :s), Ts)
rfd{T1<:Real, T2<:Real}(N::AbstractVector{T1}, D::AbstractVector{T2}) =
  rfd(Poly(reverse(N), :s), Poly(reverse(D), :s))
rfd{T1<:Real, T2<:Real}(N::AbstractVector{T1}, D::AbstractVector{T2}, Ts::Real) =
  rfd(Poly(reverse(N), :s), Poly(reverse(D), :s), Ts)

function lfd{T1<:Real, T2<:Real}(N::AbstractVector{T1}, D::AbstractVector{T2},
  Ts::Real, var::Symbol)
  vars    = [:z̄,:q̄,:qinv,:zinv]
  @assert var ∈ vars string("tf: var ∉ ", vars)

  nlast       = findlast(N)
  dlast       = findlast(D)
  order       = max(nlast, dlast)
  N_          = zeros(T1, order)
  N_[1:nlast] = N[1:nlast]
  D_          = zeros(T2, order)
  D_[1:dlast] = D[1:dlast]

  lfd(Poly(reverse(N_), :z), Poly(reverse(D_), :z), Ts)
end

function rfd{T1<:Real, T2<:Real}(N::AbstractVector{T1}, D::AbstractVector{T2},
  Ts::Real, var::Symbol)
  vars    = [:z̄,:q̄,:qinv,:zinv]
  @assert var ∈ vars string("tf: var ∉ ", vars)

  nlast       = findlast(N)
  dlast       = findlast(D)
  order       = max(nlast, dlast)
  N_          = zeros(T1, order)
  N_[1:nlast] = N[1:nlast]
  D_          = zeros(T2, order)
  D_[1:dlast] = D[1:dlast]

  rfd(Poly(reverse(N_), :z), Poly(reverse(D_), :z), Ts)
end

# Interfaces
samplingtime(s::MFD) = s.Ts
islfd{T,S,L}(s::MFD{T,S,Val{L}})  = false
islfd{T,S}(s::MFD{T,S,Val{:lfd}}) = true
isrfd{T,S,L}(s::MFD{T,S,Val{L}})  = !islfd(s)
num(s::MFD) = s.N
den(s::MFD) = s.D

# Think carefully about how to implement numstates
numstates(s::MFD)    = 1#numstates(ss(s))
# Currently, we only allow for proper systems
numstates(s::MFD{Val{:siso}}) = degree(s.D)

numinputs(s::MFD)    = s.nu
numoutputs(s::MFD)   = s.ny

# Dimension information
ndims(s::MFD{Val{:siso}})  = 1
ndims(s::MFD{Val{:mimo}}) = 2
size(s::MFD)               = size(s.N)
size(s::MFD, dim::Int)     = size(s.N, dim)
size(s::MFD, dims::Int...) = size(s.N, dims)

# Iteration interface (meaningful only for MIMO systems)
# TODO

# Slicing (`getindex`) of MIMO systems
# TODO

# Printing functions
summary(s::MFD{Val{:siso},Val{:cont}})   =
  string("tf(nx=", numstates(s), ")")
summary(s::MFD{Val{:siso},Val{:disc}})  =
  string("tf(nx=", numstates(s), ",Ts=", s.Ts, ")")
summary(s::MFD{Val{:mimo},Val{:cont}})  =
  string("tf(nx=", numstates(s), ",nu=", s.nu, ",ny=", s.ny, ")")
summary(s::MFD{Val{:mimo},Val{:disc}}) =
  string("tf(nx=", numstates(s), ",nu=", s.nu, ",ny=", s.ny, ",Ts=", s.Ts, ")")

showcompact(io::IO, s::MFD) = print(io, summary(s))

function show{T}(io::IO, s::MFD{T,Val{:cont}})
  # TODO
end

function show{T}(io::IO, s::MFD{T,Val{:disc}})
  # TODO
end

function showall(io::IO, s::MFD)
  # TODO
end

# # Conversion and promotion
# promote_rule{T<:Real,S}(::Type{T}, ::Type{MFD{Val{:siso},S}}) =
#   MFD{Val{:siso},S}
# promote_rule{T<:AbstractMatrix,S}(::Type{T}, ::Type{MFD{Val{:mimo},S}}) =
#   MFD{Val{:mimo},S}
#
# convert(::Type{RationalTF{Val{:siso},Val{:cont}}},
#   s::MFD{Val{:siso},Val{:cont}}) = tf(s.N, s.D)
#
# convert(::Type{MFD{Val{:siso},Val{:cont}}}, g::Real)             =
#   tf(g)
# convert(::Type{MFD{Val{:siso},Val{:disc}}}, g::Real)            =
#   tf(g, zero(Float64))
# convert(::Type{MFD{Val{:mimo},Val{:cont}}}, g::AbstractMatrix)  =
#   tf(g)
# convert(::Type{MFD{Val{:mimo},Val{:disc}}}, g::AbstractMatrix) =
#   tf(g, zero(Float64))
#
# # Multiplicative and additive identities (meaningful only for SISO)
# one(::Type{MFD{Val{:siso},Val{:cont}}})    =
#   tf(one(Int8))
# one(::Type{MFD{Val{:siso},Val{:disc}}})   =
#   tf(one(Int8), zero(Float64))
# zero(::Type{MFD{Val{:siso},Val{:cont}}})   =
#   tf(zero(Int8))
# zero(::Type{MFD{Val{:siso},Val{:disc}}})  =
#   tf(zero(Int8), zero(Float64))
#
# one(s::MFD{Val{:siso},Val{:cont}})   =
#   tf(one(eltype(s.num)), one(eltype(s.den)))
# one(s::MFD{Val{:siso},Val{:disc}})  =
#   tf(one(eltype(s.num)), one(eltype(s.den)), zero(Float64))
# zero(s::MFD{Val{:siso},Val{:cont}})  =
#   tf(zero(eltype(s.num)), one(eltype(s.den)))
# zero(s::MFD{Val{:siso},Val{:disc}}) =
#   tf(zero(eltype(s.num)), one(eltype(s.den)), zero(Float64))

lfd(s::MFD{Val{:siso},Val{:cont},Val{:rfd}}) = lfd(s.N,s.D)
lfd(s::MFD{Val{:siso},Val{:disc},Val{:rfd}}) = lfd(s.N,s.D,s.Ts)
lfd(s::MFD{Val{:mimo},Val{:cont},Val{:rfd}}) = lfd(_rfd2lfd(s)...)
lfd(s::MFD{Val{:mimo},Val{:disc},Val{:rfd}}) = lfd(_rfd2lfd(s)...,s.Ts)
lfd{T,S}(s::MFD{T,S,Val{:lfd}})              = s

rfd(s::MFD{Val{:siso},Val{:cont},Val{:lfd}}) = rfd(s.N,s.D)
rfd(s::MFD{Val{:siso},Val{:disc},Val{:lfd}}) = rfd(s.N,s.D,s.Ts)
rfd(s::MFD{Val{:mimo},Val{:cont},Val{:lfd}}) = rfd(_lfd2rfd(s)...)
rfd(s::MFD{Val{:mimo},Val{:disc},Val{:lfd}}) = rfd(_lfd2rfd(s)...,s.Ts)
rfd{T,S}(s::MFD{T,S,Val{:rfd}})              = s

function _lfd2rfd{T,S}(s::MFD{T,S,Val{:lfd}})
  n,m = size(s)
  p = PolyMatrix(hcat(-s.N, s.D))
  L,U = ltriang(p)

  D = U[1:m,n+1:n+m]
  N = U[m+1:n+m,n+1:n+m]
  return N,D
end

function _rfd2lfd{T,S}(s::MFD{T,S,Val{:rfd}})
  n,m = size(s)
  p = PolyMatrix(vcat(-s.D, s.N))
  _,U = rtriang(p)

  N = U[m+1:n+m,1:m]
  D = U[m+1:n+m,m+1:n+m]
  return N,D
end

## Comparison
=={T,S,L}(s1::MFD{T,S,L}, s2::MFD{T,S,L}) =
  (s1.N == s2.N) && (s1.D == s2.D) && (s1.Ts == s2.Ts)
=={T1,S1,L1,T2,S2,L2}(s1::MFD{T1,S1,L1}, s2::MFD{T2,S2,L2}) = false

hash(s::MFD, h::UInt)     = hash(s.D, hash(S.N, hash(S.Ts, h)))
isequal(s1::MFD, s2::MFD) = (hash(s1) == hash(s2))

function isapprox{T,S,L1,L2,M1,M2,M3,M4}(s1::MFD{T,S,L1,M1,M2}, s2::MFD{T,S,L2,M3,M4};
  rtol::Real=Base.rtoldefault(promote_type(eltype(eltype(s1.N)),eltype(eltype(s1.D)),eltype(eltype(s2.N)),eltype(eltype(s2.D)))),
  atol::Real=0, norm::Function=vecnorm)
  isapprox(s1.Ts, s2.Ts) || return false # quick exit
  lfd1 = lfd(s1)
  lfd2 = lfd(s2)

  D1, U = hermite(lfd1.D)
  N1    = lfd1.N*U
  D2, U = hermite(lfd2.D)
  N2    = lfd2.N*U
  return isapprox(D1, D2; rtol=rtol, atol=atol, norm=norm) && isapprox(N1,N2; rtol=rtol, atol=atol, norm=norm)
end
