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
    new{Val{:mimo},Val{:cont},Val{L},M1,M2}(PolyMatrix(coeffs(N), size(N), Val{:s}),
      PolyMatrix(coeffs(D), size(D), Val{:s}), nu, ny, zero(Float64))
    # should we do better checks than just converting the variable to the correct one?
  end

  # Discrete-time, multi-input-multi-output MFD model
  @compat function (::Type{MFD}){L,M1<:PolynomialMatrices.PolyMatrix,
    M2<:PolynomialMatrices.PolyMatrix}(N::M1, D::M2, Ts::Real, ::Type{Val{L}})
    ny, nu = mfdcheck(N, D, Val{L}, Ts)
    new{Val{:mimo},Val{:disc},Val{L},M1,M2}(PolyMatrix(coeffs(N), size(N), Val{:z}),
      PolyMatrix(coeffs(D), size(D), Val{:z}), nu, ny, convert(Float64, Ts))
    # should we do better checks than just converting the variable to the correct one?
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
  lfd(Poly(reverse(N), :z), Poly(reverse(D), :z), Ts)
rfd{T1<:Real, T2<:Real}(N::AbstractVector{T1}, D::AbstractVector{T2}) =
  rfd(Poly(reverse(N), :s), Poly(reverse(D), :s))
rfd{T1<:Real, T2<:Real}(N::AbstractVector{T1}, D::AbstractVector{T2}, Ts::Real) =
  rfd(Poly(reverse(N), :z), Poly(reverse(D), :z), Ts)

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
numstates(s::MFD)               = 1 #numstates(ss(s))
# Currently, we only allow for proper systems
numstates(s::MFD{Val{:siso}})   = degree(s.D)

numinputs(s::MFD)               = s.nu
numoutputs(s::MFD)              = s.ny

# Dimension information
ndims(s::MFD{Val{:siso}})       = 1
ndims(s::MFD{Val{:mimo}})       = 2
size(s::MFD)                    = size(s.N)
size(s::MFD, d)                 = size(s.N, d)

## Iteration interface
start(s::MFD{Val{:mimo}})       = start(s.N)
next(s::MFD{Val{:mimo}}, state) = (s[state], state+1)
done(s::MFD{Val{:mimo}}, state) = done(s.N, state)

eltype{S,M1}(::Type{MFD{Val{:mimo},Val{S},M1}}) =
  MFD{Val{:siso},Val{S},M1}

length(s::MFD{Val{:mimo}})      = length(s.N)
size(s::MFD)                    = size(s.N)
size(s::MFD, d)                 = size(s.N, d)

# Indexing of MIMO systems
function getindex(s::MFD{Val{:mimo},Val{:cont},Val{:lfd}}, I...)
  @boundscheck checkbounds(s.N, I...)
  lfd(ss(s)[I...])
end

function getindex(s::MFD{Val{:mimo},Val{:cont},Val{:rfd}}, I...)
  @boundscheck checkbounds(s.N, I...)
  rfd(ss(s)[I...])
end

endof(s::MFD{Val{:mimo}})               = endof(s.N)

# Conversion and promotion
promote_rule{T<:Real,S,L}(::Type{T}, ::Type{MFD{Val{:siso},S,L}}) =
  MFD{Val{:siso},S,L}
promote_rule{T<:AbstractMatrix,S,L}(::Type{T}, ::Type{MFD{Val{:mimo},S,L}}) =
  MFD{Val{:mimo},S,L}

convert(::Type{MFD{Val{:siso},Val{:cont},Val{:lfd}}}, g::Real)            =
  lfd(Poly(g,:s), Poly(one(g),:s))
convert(::Type{MFD{Val{:siso},Val{:disc},Val{:lfd}}}, g::Real)            =
  lfd(Poly(g,:z), Poly(one(g),:z), zero(Float64))
convert(::Type{MFD{Val{:mimo},Val{:cont},Val{:lfd}}}, g::AbstractMatrix)  =
  lfd(PolyMatrix(g, Val{:s}), PolyMatrix(eye(eltype(g), size(g,1), size(g,1))))
convert(::Type{MFD{Val{:mimo},Val{:disc},Val{:lfd}}}, g::AbstractMatrix)  =
  lfd(PolyMatrix(g, Val{:z}), PolyMatrix(eye(eltype(g), size(g,1), size(g,1))), zero(Float64))

convert(::Type{MFD{Val{:siso},Val{:cont},Val{:rfd}}}, g::Real)            =
  rfd(Poly(g,:s), Poly(one(g),:s))
convert(::Type{MFD{Val{:siso},Val{:disc},Val{:rfd}}}, g::Real)            =
  rfd(Poly(g,:z), Poly(one(g),:z), zero(Float64))
convert(::Type{MFD{Val{:mimo},Val{:cont},Val{:rfd}}}, g::AbstractMatrix)  =
  rfd(PolyMatrix(g, Val{:s}), PolyMatrix(eye(eltype(g), size(g,2), size(g,2))))
convert(::Type{MFD{Val{:mimo},Val{:disc},Val{:rfd}}}, g::AbstractMatrix)  =
  rfd(PolyMatrix(g, Val{:z}), PolyMatrix(eye(eltype(g), size(g,2), size(g,2))), zero(Float64))

# conversions betwwen lfd and rfd
convert{S,T,L}(::Type{RationalTF{Val{S},Val{T},Val{:lfd}}},
  s::MFD{Val{S},Val{T},Val{L}}) = lfd(s)
convert{S,T,L}(::Type{RationalTF{Val{S},Val{T},Val{:rfd}}},
  s::MFD{Val{S},Val{T},Val{L}}) = rfd(s)

# Multiplicative and additive identities (meaningful only for SISO)
one{M1,M2}(::Type{MFD{Val{:siso},Val{:cont},Val{:lfd},M1,M2}})  =
  lfd(one(Int8))
one{M1,M2}(::Type{MFD{Val{:siso},Val{:disc},Val{:lfd},M1,M2}})  =
  lfd(one(Int8), zero(Float64))
zero{M1,M2}(::Type{MFD{Val{:siso},Val{:cont},Val{:lfd},M1,M2}}) =
  lfd(zero(Int8))
zero{M1,M2}(::Type{MFD{Val{:siso},Val{:disc},Val{:lfd},M1,M2}}) =
  lfd(zero(Int8), zero(Float64))

one{M1,M2}(::Type{MFD{Val{:siso},Val{:cont},Val{:rfd},M1,M2}})  =
  rfd(one(Int8))
one{M1,M2}(::Type{MFD{Val{:siso},Val{:disc},Val{:rfd},M1,M2}})  =
  rfd(one(Int8), zero(Float64))
zero{M1,M2}(::Type{MFD{Val{:siso},Val{:cont},Val{:rfd},M1,M2}}) =
  rfd(zero(Int8))
zero{M1,M2}(::Type{MFD{Val{:siso},Val{:disc},Val{:rfd},M1,M2}}) =
  tf(zero(Int8), zero(Float64))

one{T,S,L,M1,M2}(s::MFD{Val{T},Val{S},Val{L},M1,M2})   =
  one(typeof(s))
zero{T,S,L,M1,M2}(s::MFD{Val{T},Val{S},Val{L},M1,M2})  =
  zero(typeof(s))

# conversions between lfd and rfd
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

function inv{M<:MFD}(s::M)
  _mfdinvcheck(s)
  _inv(s)
end

_inv{T,L}(s::MFD{Val{T},Val{:cont},Val{L}}) = lfd(copy(s.D), copy(s.N))
_inv{T,L}(s::MFD{Val{T},Val{:disc},Val{L}}) = lfd(copy(s.D), copy(s.N), s.Ts)

function _mfdinvcheck(s::MFD)
  if s.ny ≠ s.nu
    warn("inv(sys): s.ny ≠ s.nu")
    throw(DomainError())
  end

  if norm(det(s.N)) ≈ zero(eltype(s.N))
    warn("inv(sys): sys is not invertible")
    throw(DomainError())
  end
  _inv(s)
end

# Negative of a transfer-function model
-{T}(s::MFD{Val{T},Val{:cont},Val{:lfd}}) = lfd(-s.N, copy(s.D))
-{T}(s::MFD{Val{T},Val{:disc},Val{:lfd}}) = lfd(-s.N, copy(s.D), s.Ts)
-{T}(s::MFD{Val{T},Val{:cont},Val{:rfd}}) = rfd(-s.N, copy(s.D))
-{T}(s::MFD{Val{T},Val{:disc},Val{:rfd}}) = rfd(-s.N, copy(s.D), s.Ts)

# Addition (parallel)
function _mfdparallelcheck{T1,T2,S,L}(s1::MFD{Val{T1},Val{S},Val{L}},
  s2::MFD{Val{T2},Val{S},Val{L}})
  if s1.Ts ≉ s2.Ts && s1.Ts ≠ zero(s1.Ts) && s2.Ts ≠ zero(s2.Ts)
    warn("parallel(s1,s2): Sampling time mismatch")
    throw(DomainError())
  end

  if size(s1,1) ≠ size(s2,1)
    warn("parallel(s1,s2): size(s1,1) ≠ size(s2,1)")
    throw(DomainError())
  end
end

function _mfdparallel{S,L}(s₁::MFD{Val{:siso},Val{S},Val{:L}},
  s₂::MFD{Val{:siso},Val{S},Val{L}})
  R   = gcrd(s₁.D, s₂.D)
  D₁  = div(s₁.D, R)
  D   = D₁*s₂.D        # only include common part R once
  D₂  = div(s₂.D, R)
  N   = s₁.N*D₂ + s₂.N*D₁
  N, D, max(s1.Ts, s2.Ts)
end

function _mfdparallel{S}(s₁::MFD{Val{:mimo},Val{S},Val{:lfd}},
  s₂::MFD{Val{:mimo},Val{S},Val{lfd}})
  R, V₁, V₂ = gcrd(s₁.D, s₂.D)
  detV₁, adjV₁ = inv(V₁)
  detV₂, adjV₂ = inv(V₂)
  N₁ = adjV₁*s₁.N/detV₁(0)
  N₂ = adjV₂*s₂.N/detV₂(0)
  N₁+N₂, R, max(s1.Ts, s2.Ts)
end

function _mfdparallel{S}(s₁::MFD{:mimo},Val{S},Val{:lfd}},
  s₂::MFD{Val{:mimo},Val{S},Val{lfd}})
  R, V₁, V₂ = gcrd(s₁, s₂)
  detV₁, adjV₁ = inv(V₁)
  detV₂, adjV₂ = inv(V₂)
  N₁ = adjV₁*s₁.N/detV₁(0)
  N₂ = adjV₂*s₂.N/detV₂(0)
  N₁+N₂, R, max(s1.Ts, s2.Ts)
end

function +{T1,T2,L}(s1::MFD{Val{T1},Val{:cont},Val{L}},
  s2::MFD{Val{T2},Val{:cont},Val{L}})
  _mfdparallelcheck(s1, s2)
  N, D, _ = _mfdparallel(s1, s2)
  MFD(N, D, Val{L})
end

function +{T1,T2,L}(s1::MFD{Val{T1},Val{:disc},Val{L}},
  s2::MFD{Val{T2},Val{:disc},Val{L}})
  _mfdparallelcheck(s1, s2)
  N, D, Ts = _mfdparallel(s1, s2)
  MFD(N, D, Ts, Val{L})
end

function +{T1,T2,L1,L2}(s1::MFD{Val{T1},Val{:cont},Val{L1}},
  s2::MFD{Val{T2},Val{:cont},Val{L2}})
  N, D, _ = _mfdparallel(lfd(s1), lfd(s2))
  lfd(N, D, Val{L})
end

function +{T1,T2,L1,L2}(s1::MFD{Val{T1},Val{:disc},Val{L1}},
  s2::MFD{Val{T2},Val{:disc},Val{L2}})
  N, D, Ts = _mfdparallel(lfd(s1), lfd(s2))
  lfd(N, D, Ts, Val{L})
end

.+(s1::MFD{Val{:siso}}, s2::MFD{Val{:siso}}) = +(s1, s2)

+{T,S}(s::MFD{Val{T},Val{S}}, g::Union{Real,AbstractMatrix}) =
  +(s, convert(typeof(s), g))
+{T,S}(g::Union{Real,AbstractMatrix}, s::MFD{Val{T},Val{S}}) =
  +(convert(typeof(s), g), s)

.+(s::MFD{Val{:siso}}, g::Real) = +(s, g)
.+(g::Real, s::MFD{Val{:siso}}) = +(g, s)

# Subtraction
-(s1::MFD, s2::MFD) = +(s1, -s2)

.-(s1::MFD{Val{:siso}}, s2::MFD{Val{:siso}}) = -(s1, s2)

-{T,S}(s::MFD{Val{T},Val{S}}, g::Union{Real,AbstractMatrix}) =
  -(s, convert(typeof(s), g))
-{T,S}(g::Union{Real,AbstractMatrix}, s::MFD{Val{T},Val{S}}) =
  -(convert(typeof(s), g), s)

.-(s::MFD{Val{:siso}}, g::Real)    = -(s, g)
.-(g::Real, s::MFD{Val{:siso}})    = -(g, s)

# Multiplication
function _mfdseries{T1,T2,S}(s1::MFD{Val{T1},Val{S}},
  s2::MFD{Val{T2},Val{S}})
  # Remark: s1*s2 implies u -> s2 -> s1 -> y

  if s1.Ts ≉ s2.Ts && s1.Ts ≠ zero(s1.Ts) && s2.Ts ≠ zero(s2.Ts)
    warn("series(s1,s2): Sampling time mismatch")
    throw(DomainError())
  end

  if size(s1,2) ≠ size(s2,1)
    warn("parallel(s1,s2): size(s1,2) ≠ size(s2,1)")
    throw(DomainError())
  end
end

function _mfdseries{T1,T2,S}(s1::MFD{Val{T1},Val{S}},
  s2::MFD{Val{T2},Val{S}})

  return N, D, max(s1.Ts, s2.Ts)
end

function *(s1::MFD{Val{:siso},Val{:cont}},
  s2::MFD{Val{:siso},Val{:cont}})
  a, b, c, d, _ = _ssseries(s1, s2)
  MFD(a, b, c, d[1])
end

function *(s1::MFD{Val{:siso},Val{:disc}},
  s2::MFD{Val{:siso},Val{:disc}})
  a, b, c, d, Ts = _ssseries(s1, s2)
  MFD(a, b, c, d[1], Ts)
end

function *{T1,T2}(s1::MFD{Val{T1},Val{:cont}},
  s2::MFD{Val{T2},Val{:cont}})
  a, b, c, d, _ = _ssseries(s1, s2)
  MFD(a, b, c, d)
end

function *{T1,T2}(s1::MFD{Val{T1},Val{:disc}},
  s2::MFD{Val{T2},Val{:disc}})
  a, b, c, d, Ts = _ssseries(s1, s2)
  MFD(a, b, c, d, Ts)
end

.*(s1::MFD{Val{:siso}}, s2::MFD{Val{:siso}}) = *(s1, s2)

*{T,S}(s::MFD{Val{T},Val{S}}, g::Union{Real,AbstractMatrix}) =
  *(s, convert(typeof(s), g))
*{T,S}(g::Union{Real,AbstractMatrix}, s::MFD{Val{T},Val{S}}) =
  *(convert(typeof(s), g), s)

.*(s::MFD{Val{:siso}}, g::Real)    = *(s, g)
.*(g::Real, s::MFD{Val{:siso}})    = *(g, s)

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
