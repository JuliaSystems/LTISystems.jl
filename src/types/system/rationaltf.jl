immutable RationalTF{T,S,U,V} <: LtiSystem{T,S}
  mat::V
  nu::Int
  ny::Int
  Ts::Float64

  # Continuous-time, single-input-single-output rational transfer function model
  function (::Type{RationalTF}){S,U<:Real,V<:Real}(r::RationalFunction{
    Val{:s},Val{S},U,V})
    mat     = fill(r, 1, 1)
    ny, nu  = _tfcheck(mat)
    new{Val{:siso},Val{:cont},Val{S},typeof(mat)}(mat, nu, ny, zero(Float64))
  end

  # Discrete-time, single-input-single-output rational transfer function model
  function (::Type{RationalTF}){S,U<:Real,V<:Real}(r::RationalFunction{
    Val{:z},Val{S},U,V}, Ts::Real)
    mat     = fill(r, 1, 1)
    ny, nu  = _tfcheck(mat, Ts)
    new{Val{:siso},Val{:disc},Val{S},typeof(mat)}(mat, nu, ny, convert(Float64, Ts))
  end

  # Continuous-time, multi-input-multi-output rational transfer function model
  function (::Type{RationalTF}){S,U<:Real,V<:Real}(
    mat::AbstractMatrix{RationalFunction{Val{:s},Val{S},U,V}})
    ny, nu  = _tfcheck(mat)
    new{Val{:mimo},Val{:cont},Val{S},typeof(mat)}(mat, nu, ny, zero(Float64))
  end
  function (::Type{RationalTF}){T<:Real,S}(mat::AbstractMatrix{T},
    t::Type{Val{S}} = Val{:notc})
    m = map(x->RationalFunction(x, :s, t), mat)
    ny, nu  = _tfcheck(m)
    new{Val{:mimo},Val{:cont},t,typeof(m)}(m, nu, ny, zero(Float64))
  end

  # Discrete-time, multi-input-multi-output rational transfer function model
  function (::Type{RationalTF}){S,U<:Real,V<:Real}(
    mat::AbstractMatrix{RationalFunction{Val{:z},Val{S},U,V}}, Ts::Real)
    ny, nu  = _tfcheck(mat, Ts)
    new{Val{:mimo},Val{:disc},Val{S},typeof(mat)}(mat, nu, ny, convert(Float64, Ts))
  end
  function (::Type{RationalTF}){T<:Real,S}(mat::AbstractMatrix{T},
    Ts::Real, t::Type{Val{S}} = Val{:notc})
    m = map(x->RationalFunction(x, :z, t), mat)
    ny, nu  = _tfcheck(m)
    new{Val{:mimo},Val{:disc},t,typeof(m)}(m, nu, ny, convert(Float64, Ts))
  end

  # Function calls
  # Evaluate system models at given complex numbers
  (sys::RationalTF{Val{:siso}})(x::Number)                                    =
    (f = sys.mat[1]; convert(Complex128, f(x)))
  (sys::RationalTF{Val{:siso}}){M<:Number}(X::AbstractArray{M})               =
    (f = sys.mat[1]; Complex128[f(x) for x in X])
  (sys::RationalTF{Val{:mimo}})(x::Number)                                    =
    Complex128[f(x) for f in sys.mat]
  (sys::RationalTF{Val{:mimo}}){M<:Number}(X::AbstractArray{M})               =
    Complex128[f(x) for f in sys.mat, x in X]
  function (sys::RationalTF){T<:Real}(; ω::Union{T, AbstractArray{T}} = Float64[])
    if isempty(ω)
      warn("sys(): Provide an argument for the function call. Refer to `?freqresp`.")
      throw(DomainError())
    end
    freqresp(sys, ω)
  end
end

# Warn the user in other type constructions
function (::Type{RationalTF})(r::RationalFunction)
  warn("tf(r): `r` can only be a real-coefficient `RationalFunction` of variable `:s`")
  throw(DomainError())
end

function (::Type{RationalTF})(r::RationalFunction, Ts::Real)
  warn("tf(r, Ts): `r` can only be a real-coefficient `RationalFunction` of variable `:z`")
  throw(DomainError())
end

function (::Type{RationalTF})(m::AbstractMatrix)
  warn("tf(m): `m` can only be an `AbstractMatrix` of real-coefficient `RationalFunction` objects of variable `:s`")
  throw(DomainError())
end

function (::Type{RationalTF})(m::AbstractMatrix, Ts::Real)
  warn("tf(m, Ts): `m` can only be an `AbstractMatrix` of real-coefficient `RationalFunction` objects of variable `:z`")
  throw(DomainError())
end

# Enforce rational transfer function type invariance
function _tfcheck{T,S,U<:Real,V<:Real}(mat::AbstractMatrix{RationalFunction{Val{T},
  Val{S},U,V}}, Ts::Real = zero(Float64))
  # Check sampling time
  if Ts < zero(Ts) || isinf(Ts)
    warn("tf(m, Ts): `Ts` must be a non-negative real number")
    throw(DomainError())
  end

  # Check input-output dimensions
  ny, nu = size(mat)

  if ny == 0 || nu == 0
    warn("tf(m[, Ts]): `min(nu, ny) = 0`")
    throw(DomainError())
  end

  # # Check properness, etc.
  # for idx in 1:length(mat)
  #   col, row        = divrem(idx-1, ny)
  #   numdeg, dendeg  = degree(mat[idx])
  #   denpoly         = den(mat[idx])
  #
  #   if numdeg > dendeg
  #     warn("tf(m[, Ts]): `m[$(row+1),$(col+1)]` is not proper")
  #     throw(DomainError())
  #   end
  #   if denpoly ≈ zero(denpoly)
  #     warn("tf(m[, Ts]): `m[$(row+1),$(col+1)]` has a zero denominator")
  #     throw(DomainError())
  #   end
  # end

  return ny, nu
end

# Helper function
_reverse(x)         = reverse(x)
_reverse(x::Number) = [x]
const NV            = Union{Number,Vector}

# Outer constructors
# Continuous-time, single-input-single-output rational transfer function model
tf(r::RationalFunction) = RationalTF(r)
tf(num::NV, den::NV)    = RationalTF(RationalFunction(_reverse(num), _reverse(den), :s))
tf(num::NV, den::Poly)  = RationalTF(RationalFunction(_reverse(num), den))
tf(num::Poly, den::NV)  = RationalTF(RationalFunction(num, _reverse(den)))
tf(num::Poly, den::Poly)= RationalTF(RationalFunction(num, den))

# Discrete-time, single-input-single-output rational transfer function model
tf(r::RationalFunction, Ts::Real) = RationalTF(r, Ts)
tf(num::NV, den::NV, Ts::Real)    = RationalTF(RationalFunction(_reverse(num), _reverse(den), :z), Ts)
tf(num::NV, den::Poly, Ts::Real)  = RationalTF(RationalFunction(_reverse(num), den), Ts)
tf(num::Poly, den::NV, Ts::Real)  = RationalTF(RationalFunction(num, _reverse(den)), Ts)
tf(num::Poly, den::Poly, Ts::Real)= RationalTF(RationalFunction(num, den), Ts)

function tf(num::AbstractVector, den::AbstractVector, Ts::Real, var::Symbol)
  vars    = [:z̄,:q̄,:qinv,:zinv]
  if var ∉ vars
    warn("tf(num,den,Ts,var): `var` ∉ ", vars)
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

# Continuous-time, multi-input-multi-output rational transfer function model
tf{S,U,V}(mat::AbstractMatrix{RationalFunction{Val{:s},Val{S},U,V}})            =
  RationalTF(mat)

# Discrete-time, multi-input-multi-output rational transfer function model
tf{S,U,V}(mat::AbstractMatrix{RationalFunction{Val{:z},Val{S},U,V}}, Ts::Real)  =
  RationalTF(mat, Ts)

# Continuous-time, multi-input-multi-output rational transfer function model
tf{T<:Real,S}(mat::AbstractMatrix{T}, t::Type{Val{S}} = Val{:notc})             =
  RationalTF(mat, t)

# Discrete-time, multi-input-multi-output rational transfer function model
tf{T<:Real,S}(mat::AbstractMatrix{T}, Ts::Real, t::Type{Val{S}} = Val{:notc})   =
  RationalTF(mat, Ts, t)

# Interfaces
samplingtime(s::RationalTF) = s.Ts

numinputs(s::RationalTF)    = s.nu
numoutputs(s::RationalTF)   = s.ny

# Iteration interface
start(s::RationalTF{Val{:mimo}})        = start(s.mat)
next(s::RationalTF{Val{:mimo}}, state)  = (s[state], state+1)
done(s::RationalTF{Val{:mimo}}, state)  = done(s.mat, state)

eltype{S,U,V}(::Type{RationalTF{Val{:mimo},Val{S},Val{U},V}}) =
  RationalTF{Val{:siso},Val{S},Val{U},Matrix{eltype(V)}}

length(s::RationalTF{Val{:mimo}}) = length(s.mat)
size(s::RationalTF)               = size(s.mat)
size(s::RationalTF, d)            = size(s.mat, d)

# Indexing of MIMO systems
function getindex(s::RationalTF{Val{:mimo},Val{:cont}}, row::Int, col::Int)
  (1 ≤ row ≤ s.ny && 1 ≤ col ≤ s.nu) || throw(BoundsError(s.mat, (row,col)))
  RationalTF(s.mat[row,col])
end

function getindex(s::RationalTF{Val{:mimo},Val{:disc}}, row::Int, col::Int)
  (1 ≤ row ≤ s.ny && 1 ≤ col ≤ s.nu) || throw(BoundsError(s.mat, (row,col)))
  RationalTF(s.mat[row,col], s.Ts)
end

function getindex(s::RationalTF{Val{:mimo}}, idx::Int)
  (1 ≤ idx ≤ length(s.mat)) || throw(BoundsError(s.mat, idx))
  col, row  = divrem(idx-1, s.ny)
  s[row+1, col+1]
end

function getindex(s::RationalTF{Val{:mimo},Val{:cont}}, rows::AbstractVector{Int},
  cols::AbstractVector{Int})
  1 ≤ minimum(rows) ≤ maximum(rows) ≤ s.ny || throw(BoundsError(s.mat, rows))
  1 ≤ minimum(cols) ≤ maximum(cols) ≤ s.nu || throw(BoundsError(s.mat, cols))

  RationalTF(view(s.mat, rows, cols))
end

function getindex(s::RationalTF{Val{:mimo},Val{:disc}}, rows::AbstractVector{Int},
  cols::AbstractVector{Int})
  1 ≤ minimum(rows) ≤ maximum(rows) ≤ s.ny || throw(BoundsError(s.mat, rows))
  1 ≤ minimum(cols) ≤ maximum(cols) ≤ s.nu || throw(BoundsError(s.mat, cols))

  RationalTF(view(s.mat, rows, cols), s.Ts)
end

function getindex(s::RationalTF{Val{:mimo}}, indices::AbstractVector{Int})
  1 ≤ minimum(indices) ≤ maximum(indices) ≤ length(s.mat) || throw(BoundsError(s.mat, indices))

  temp  = map(x->divrem(x-1, s.ny), indices)
  cols  = map(x->x[1]+1, temp)
  rows  = map(x->x[2]+1, temp)

  s[rows, cols]
end

getindex(s::RationalTF{Val{:mimo}}, rows, ::Colon)    = s[rows, 1:s.nu]
getindex(s::RationalTF{Val{:mimo}}, ::Colon, cols)    = s[1:s.ny, cols]
getindex(s::RationalTF{Val{:mimo}}, ::Colon)          = s[1:end]
getindex(s::RationalTF{Val{:mimo}}, ::Colon, ::Colon) = s[1:s.ny,1:s.nu]
endof(s::RationalTF{Val{:mimo}})                      = endof(s.mat)

# Conversion and promotion
promote_rule{T1<:Real,T2,S,U,V}(::Type{T1}, ::Type{RationalTF{Val{T2},Val{S},Val{U},V}}) =
  RationalTF{Val{T2},Val{S},Val{U}}
promote_rule{T<:AbstractMatrix,S,U,V}(::Type{T}, ::Type{RationalTF{Val{:mimo},Val{S},Val{U},V}}) =
  RationalTF{Val{:mimo},Val{S},Val{U}}

convert{U,V}(::Type{RationalTF{Val{:siso},Val{:cont},Val{U},V}}, g::Real) =
  tf(RationalFunction(g, :s, Val{U}))
convert{U,V}(::Type{RationalTF{Val{:siso},Val{:disc},Val{U},V}}, g::Real) =
  tf(RationalFunction(g, :z, Val{U}), zero(Float64))
convert{U,V}(::Type{RationalTF{Val{:mimo},Val{:cont},Val{U},V}}, g::Real) =
  tf(fill(g,1,1), Val{U})
convert{U,V}(::Type{RationalTF{Val{:mimo},Val{:disc},Val{U},V}}, g::Real) =
  tf(fill(g,1,1), zero(Float64), Val{U})
convert{U,V}(::Type{RationalTF{Val{:mimo},Val{:cont},Val{U},V}}, g::AbstractMatrix) =
  tf(g)
convert{U,V}(::Type{RationalTF{Val{:mimo},Val{:disc},Val{U},V}}, g::AbstractMatrix) =
  tf(g, zero(Float64))

# Multiplicative and additive identities (meaningful only for SISO)
one{U,V}(::Type{RationalTF{Val{:siso},Val{:cont},Val{U},V}})  =
  tf(one(eltype(V)))
one{U,V}(::Type{RationalTF{Val{:siso},Val{:disc},Val{U},V}})  =
  tf(one(eltype(V)), zero(Float64))
zero{U,V}(::Type{RationalTF{Val{:siso},Val{:cont},Val{U},V}}) =
  tf(zero(eltype(V)))
zero{U,V}(::Type{RationalTF{Val{:siso},Val{:disc},Val{U},V}}) =
  tf(zero(eltype(V)), zero(Float64))

one(s::RationalTF{Val{:siso},Val{:cont}})   = one(typeof(s))
one(s::RationalTF{Val{:siso},Val{:disc}})   = one(typeof(s))
zero(s::RationalTF{Val{:siso},Val{:cont}})  = zero(typeof(s))
zero(s::RationalTF{Val{:siso},Val{:disc}})  = zero(typeof(s))

# Inverse of a transfer-function model
function _tfinv(s::RationalTF)
  if s.ny ≠ s.nu
    warn("inv(sys): s.ny ≠ s.nu")
    throw(DomainError())
  end

  try
    return inv(s.mat)
  catch err
    warn("inv(sys): sys is not invertible")
    throw(DomainError())
  end
end

function inv(s::RationalTF{Val{:siso},Val{:cont}})
  mat = _tfinv(s)
  RationalTF(mat[1])
end

function inv(s::RationalTF{Val{:siso},Val{:disc}})
  mat = _tfinv(s)
  RationalTF(mat[1], s.Ts)
end

function inv(s::RationalTF{Val{:mimo},Val{:cont}})
  mat = _tfinv(s)
  RationalTF(mat)
end

function inv(s::RationalTF{Val{:mimo},Val{:disc}})
  mat = _tfinv(s)
  RationalTF(mat, s.Ts)
end

# Invariant zeros of a transfer-function model
zeros(s::RationalTF{Val{:siso}}) = roots(Base.num(s.mat[1]))
function zeros(s::RationalTF{Val{:mimo}})
  # TODO: Implement an efficient version. For now, the below seems to be working
  poles(1/det(s.mat))
end

# Transmission zeros of a transfer-function model
tzeros(s::RationalTF) = zeros(minreal(s))

# Poles of a transfer-function model
poles(s::RationalTF{Val{:siso}}) = roots(Base.den(s.mat[1]))
function poles(s::RationalTF{Val{:mimo}})
  poleset = Set{Complex128}()
  for r in s.mat
    push!(poleset, poles(r)...)
  end
  return [pole for pole in poleset]
end

# Negative of a transfer-function model
-(s::RationalTF{Val{:siso},Val{:cont}}) = RationalTF(-s.mat[1])
-(s::RationalTF{Val{:siso},Val{:disc}}) = RationalTF(-s.mat[1], s.Ts)
-(s::RationalTF{Val{:mimo},Val{:cont}}) = RationalTF(-s.mat)
-(s::RationalTF{Val{:mimo},Val{:disc}}) = RationalTF(-s.mat, s.Ts)

# Addition
function _tfparallel{T1,T2,S,U}(s1::RationalTF{Val{T1},Val{S},Val{U}},
  s2::RationalTF{Val{T2},Val{S},Val{U}})
  if s1.Ts ≉ s2.Ts && s1.Ts ≠ zero(s1.Ts) && s2.Ts ≠ zero(s2.Ts)
    warn("parallel(s1,s2): Sampling time mismatch")
    throw(DomainError())
  end

  if size(s1) ≠ size(s2)
    warn("parallel(s1,s2): size(s1) ≠ size(s2)")
    throw(DomainError())
  end

  return s1.mat + s2.mat, max(s1.Ts, s2.Ts)
end

function +{U}(s1::RationalTF{Val{:siso},Val{:cont},Val{U}},
  s2::RationalTF{Val{:siso},Val{:cont},Val{U}})
  mat, _ = _tfparallel(s1, s2)
  RationalTF(mat[1])
end

function +{U}(s1::RationalTF{Val{:siso},Val{:disc},Val{U}},
  s2::RationalTF{Val{:siso},Val{:disc},Val{U}})
  mat, Ts = _tfparallel(s1, s2)
  RationalTF(mat[1], Ts)
end

function +{T1,T2,U}(s1::RationalTF{Val{T1},Val{:cont},Val{U}},
  s2::RationalTF{Val{T2},Val{:cont},Val{U}})
  mat, _ = _tfparallel(s1, s2)
  RationalTF(mat)
end

function +{T1,T2,U}(s1::RationalTF{Val{T1},Val{:disc},Val{U}},
  s2::RationalTF{Val{T2},Val{:disc},Val{U}})
  mat, Ts = _tfparallel(s1, s2)
  RationalTF(mat, Ts)
end

.+(s1::RationalTF{Val{:siso}}, s2::RationalTF{Val{:siso}}) = +(s1, s2)

+{T,S,U}(s::RationalTF{Val{T},Val{S},Val{U}}, g::Union{Real,AbstractMatrix})  =
  +(s, convert(typeof(s), g))
+{T,S,U}(g::Union{Real,AbstractMatrix}, s::RationalTF{Val{T},Val{S},Val{U}})  =
  +(convert(typeof(s), g), s)

.+(s::RationalTF{Val{:siso}}, g::Real)    = +(s, g)
.+(g::Real, s::RationalTF{Val{:siso}})    = +(g, s)

# Subtraction
-(s1::RationalTF, s2::RationalTF) = +(s1, -s2)

.-(s1::RationalTF{Val{:siso}}, s2::RationalTF{Val{:siso}}) = -(s1, s2)

-{T,S,U}(s::RationalTF{Val{T},Val{S},Val{U}}, g::Union{Real,AbstractMatrix})  =
  -(s, convert(typeof(s), g))
-{T,S,U}(g::Union{Real,AbstractMatrix}, s::RationalTF{Val{T},Val{S},Val{U}})  =
  -(convert(typeof(s), g), s)

.-(s::RationalTF{Val{:siso}}, g::Real)    = -(s, g)
.-(g::Real, s::RationalTF{Val{:siso}})    = -(g, s)

# Multiplication
function _tfseries{T1,T2,S,U}(s1::RationalTF{Val{T1},Val{S},Val{U}},
  s2::RationalTF{Val{T2},Val{S},Val{U}})
  # Remark: s1*s2 implies u -> s2 -> s1 -> y

  if s1.Ts ≉ s2.Ts && s1.Ts ≠ zero(s1.Ts) && s2.Ts == zero(s2.Ts)
    warn("series(s1,s2): Sampling time mismatch")
    throw(DomainError())
  end

  if s1.nu ≠ s2.ny
    warn("series(s1,s2): s1.nu ≠ s2.ny")
    throw(DomainError())
  end

  # NOTE: Trick borrowed from Julia v0.6 implementation of `matmul`.
  el1 = one(eltype(s1.mat))
  el2 = one(eltype(s2.mat))
  el  = el1*el2 + el1*el2

  temp::AbstractMatrix{typeof(el)} = s1.mat * s2.mat

  return temp, max(s1.Ts, s2.Ts)
end

function *{U}(s1::RationalTF{Val{:siso},Val{:cont},Val{U}},
  s2::RationalTF{Val{:siso},Val{:cont},Val{U}})
  mat, _ = _tfseries(s1, s2)
  RationalTF(mat[1])
end

function *{U}(s1::RationalTF{Val{:siso},Val{:disc},Val{U}},
  s2::RationalTF{Val{:siso},Val{:disc},Val{U}})
  mat, Ts = _tfseries(s1, s2)
  RationalTF(mat[1], Ts)
end

function *{T1,T2,U}(s1::RationalTF{Val{T1},Val{:cont},Val{U}},
  s2::RationalTF{Val{T2},Val{:cont},Val{U}})
  mat, _ = _tfseries(s1, s2)
  RationalTF(mat)
end

function *{T1,T2,U}(s1::RationalTF{Val{T1},Val{:disc},Val{U}},
  s2::RationalTF{Val{T2},Val{:disc},Val{U}})
  mat, Ts = _tfseries(s1, s2)
  RationalTF(mat, Ts)
end

.*(s1::RationalTF{Val{:siso}}, s2::RationalTF{Val{:siso}}) = *(s1, s2)

*{T,S,U}(s::RationalTF{Val{T},Val{S},Val{U}}, g::Union{Real,AbstractMatrix})  =
  *(s, convert(typeof(s), g))
*{T,S,U}(g::Union{Real,AbstractMatrix}, s::RationalTF{Val{T},Val{S},Val{U}})  =
  *(convert(typeof(s), g), s)

.*(s::RationalTF{Val{:siso}}, g::Real)    = *(s, g)
.*(g::Real, s::RationalTF{Val{:siso}})    = *(g, s)

# Division
/(s1::RationalTF, s2::RationalTF)         = *(s1, inv(s2))

./(s1::RationalTF{Val{:siso}}, s2::RationalTF{Val{:siso}}) = /(s1, s2)

/{T,S,U}(s::RationalTF{Val{T},Val{S},Val{U}}, g::Union{Real,AbstractMatrix})  =
  /(s, convert(typeof(s), g))
/{T,S,U}(g::Union{Real,AbstractMatrix}, s::RationalTF{Val{T},Val{S},Val{U}})  =
  /(convert(typeof(s), g), s)

./(s::RationalTF{Val{:siso}}, g::Real)    = /(s, g)
./(g::Real, s::RationalTF{Val{:siso}})    = /(g, s)
