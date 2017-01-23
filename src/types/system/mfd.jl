# Value type to differentiate between MFD of left and right type
immutable Lfd{T}
end

# Parameters:
#   T:  Siso{true} or Siso{false}
#   S:  Continuous{true} or Continuous{false}
#   L:  Lfd{true} or Lfd{false}
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
  @compat function (::Type{MFD}){L,M1<:Polynomials.Poly,M2<:Polynomials.Poly}(N::M1, D::M2, ::Type{Lfd{L}})
    mfdcheck(N,D)
    new{Siso{true},Continuous{true},Lfd{L},M1,M2}(N, D, 1, 1, zero(Float64))
  end

  # Discrete-time, single-input-single-output MFD model
  @compat function (::Type{MFD}){L,M1<:Polynomials.Poly,M2<:Polynomials.Poly}(N::M1, D::M2, Ts::Real, ::Type{Lfd{L}})
    mfdcheck(N,D,Ts)
    new{Siso{true},Continuous{false},Lfd{L},M1,M2}(N, D, 1, 1, convert(Float64, Ts))
  end

  # Continuous-time, multi-input-multi-output MFD model
  @compat function (::Type{MFD}){L,M1<:PolynomialMatrices.PolyMatrix,
    M2<:PolynomialMatrices.PolyMatrix}(N::M1, D::M2, ::Type{Lfd{L}})
    ny, nu = mfdcheck(N, D, Lfd{L})
    new{Siso{false},Continuous{true},Lfd{L},M1,M2}(N, D, nu, ny, zero(Float64))
  end

  # Discrete-time, multi-input-multi-output MFD model
  @compat function (::Type{MFD}){L,M1<:PolynomialMatrices.PolyMatrix,
    M2<:PolynomialMatrices.PolyMatrix}(N::M1, D::M2, Ts::Real, ::Type{Lfd{L}})
    ny, nu = mfdcheck(N, D, Lfd{L}, Ts)
    new{Siso{false},Continuous{false},Lfd{L},M1,M2}(N, D, nu, ny, convert(Float64, Ts))
  end

end

function mfdcheck{T<:Real,S<:Real}(N::Poly{T}, D::Poly{S}, Ts::Real = zero(Float64))
  @assert Ts ≥ zero(Ts) && !isinf(Ts) "MFD: Ts must be non-negative real number"
end

# Enforce rational transfer function type invariance
function mfdcheck{M1<:PolynomialMatrices.PolyMatrix,M2<:PolynomialMatrices.PolyMatrix}(
  N::M1, D::M2, ::Type{Lfd{true}})
  @assert size(N,1) == size(D,1) "MFD: size(N,1) ≠ size(D,1)"
  @assert size(D,1) == size(D,2) "MFD: size(D,1) ≠ size(D,2)"

  return size(N,1), size(N,2)
end
function mfdcheck{M1<:PolynomialMatrices.PolyMatrix,M2<:PolynomialMatrices.PolyMatrix}(
  N::M1, D::M2, ::Type{Lfd{false}})
  @assert size(N,2) == size(D,2) "MFD: size(N,2) ≠ size(D,2)"
  @assert size(D,1) == size(D,2) "MFD: size(D,1) ≠ size(D,2)"

  return size(N,1), size(N,2)
end
function mfdcheck{T,M1<:PolynomialMatrices.PolyMatrix,M2<:PolynomialMatrices.PolyMatrix}(
  N::M1, D::M2, ::Type{Lfd{T}}, Ts::Real)
  @assert Ts ≥ zero(Ts) && !isinf(Ts) "MFD: Ts must be non-negative real number"
  return mfdcheck(N, D, Lfd{T})
end

# Outer constructors
lfd(N::Poly, D::Poly)           = MFD(N, D, Lfd{true})
lfd(N::Poly, D::Poly, Ts::Real) = MFD(N, D, convert(Float64, Ts), Lfd{true})
lfd{M1<:PolyMatrix,M2<:PolyMatrix}(N::M1, D::M2) =
  MFD(N, D, Lfd{true})
lfd{M1<:PolyMatrix,M2<:PolyMatrix}(N::M1, D::M2, Ts::Real) =
  MFD(N, D, convert(Float64, Ts), Lfd{true})

rfd(N::Poly, D::Poly) = MFD(N, D, Lfd{false})
rfd(N::Poly, D::Poly, Ts::Real) = MFD(N, D, convert(Float64, Ts), Lfd{false})
rfd(N::PolyMatrix, D::PolyMatrix) = MFD(N, D, Lfd{false})
rfd(N::PolyMatrix, D::PolyMatrix, Ts::Real) = MFD(N, D, convert(Float64, Ts), Lfd{false})

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
islfd{T,S,L}(s::MFD{T,S,Lfd{L}}) = L::Bool
isrfd{T,S,L}(s::MFD{T,S,Lfd{L}}) = !L::Bool
num(s::MFD) = s.N
den(s::MFD) = s.D

# Think carefully about how to implement numstates
numstates(s::MFD)    = 1#numstates(ss(s))
# Currently, we only allow for proper systems
numstates(s::MFD{Siso{true}}) = degree(s.D)

numinputs(s::MFD)    = s.nu
numoutputs(s::MFD)   = s.ny

# Dimension information
ndims(s::MFD{Siso{true}})  = 1
ndims(s::MFD{Siso{false}}) = 2
size(s::MFD)               = size(s.N)
size(s::MFD, dim::Int)     = size(s.N, dim)
size(s::MFD, dims::Int...) = size(s.N, dims)

# Iteration interface (meaningful only for MIMO systems)
# TODO

# Slicing (`getindex`) of MIMO systems
# TODO

# Printing functions
summary(s::MFD{Siso{true},Continuous{true}})   =
  string("tf(nx=", numstates(s), ")")
summary(s::MFD{Siso{true},Continuous{false}})  =
  string("tf(nx=", numstates(s), ",Ts=", s.Ts, ")")
summary(s::MFD{Siso{false},Continuous{true}})  =
  string("tf(nx=", numstates(s), ",nu=", s.nu, ",ny=", s.ny, ")")
summary(s::MFD{Siso{false},Continuous{false}}) =
  string("tf(nx=", numstates(s), ",nu=", s.nu, ",ny=", s.ny, ",Ts=", s.Ts, ")")

showcompact(io::IO, s::MFD) = print(io, summary(s))

function show{T}(io::IO, s::MFD{T,Continuous{true}})
  # TODO
end

function show{T}(io::IO, s::MFD{T,Continuous{false}})
  # TODO
end

function showall(io::IO, s::MFD)
  # TODO
end

# Conversion and promotion
promote_rule{T<:Real,S}(::Type{T}, ::Type{MFD{Siso{true},S}}) =
  MFD{Siso{true},S}
promote_rule{T<:AbstractMatrix,S}(::Type{T}, ::Type{MFD{Siso{false},S}}) =
  MFD{Siso{false},S}

convert(::Type{RationalTF{Siso{true},Continuous{true}}},
  s::MFD{Siso{true},Continuous{true}}) = tf(s.N, s.D)

convert(::Type{MFD{Siso{true},Continuous{true}}}, g::Real)             =
  tf(g)
convert(::Type{MFD{Siso{true},Continuous{false}}}, g::Real)            =
  tf(g, zero(Float64))
convert(::Type{MFD{Siso{false},Continuous{true}}}, g::AbstractMatrix)  =
  tf(g)
convert(::Type{MFD{Siso{false},Continuous{false}}}, g::AbstractMatrix) =
  tf(g, zero(Float64))

# Multiplicative and additive identities (meaningful only for SISO)
one(::Type{MFD{Siso{true},Continuous{true}}})    =
  tf(one(Int8))
one(::Type{MFD{Siso{true},Continuous{false}}})   =
  tf(one(Int8), zero(Float64))
zero(::Type{MFD{Siso{true},Continuous{true}}})   =
  tf(zero(Int8))
zero(::Type{MFD{Siso{true},Continuous{false}}})  =
  tf(zero(Int8), zero(Float64))

one(s::MFD{Siso{true},Continuous{true}})   =
  tf(one(eltype(s.num)), one(eltype(s.den)))
one(s::MFD{Siso{true},Continuous{false}})  =
  tf(one(eltype(s.num)), one(eltype(s.den)), zero(Float64))
zero(s::MFD{Siso{true},Continuous{true}})  =
  tf(zero(eltype(s.num)), one(eltype(s.den)))
zero(s::MFD{Siso{true},Continuous{false}}) =
  tf(zero(eltype(s.num)), one(eltype(s.den)), zero(Float64))
