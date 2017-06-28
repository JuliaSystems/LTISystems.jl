immutable StateSpace{T,S,M1,M2,M3,M4} <: LtiSystem{T,S}
  A::M1
  B::M2
  C::M3
  D::M4
  nx::Int
  nu::Int
  ny::Int
  Ts::Float64

  # Constructors
  # Continuous-time, single-input-single-output state-space model
  function (::Type{StateSpace}){M1<:AbstractMatrix,M2<:AbstractMatrix,
    M3<:AbstractMatrix}(A::M1, B::M2, C::M3, D::Real)
    d = fill(D,1,1)
    nx, nu, ny = _sscheck(A, B, C, d)
    new{Val{:siso},Val{:cont},M1,M2,M3,typeof(d)}(A, B, C, d, nx, nu, ny, zero(Float64))
  end

  # Discrete-time, single-input-single-output state-space model
  function (::Type{StateSpace}){M1<:AbstractMatrix,M2<:AbstractMatrix,
    M3<:AbstractMatrix}(A::M1, B::M2, C::M3, D::Real, Ts::Real)
    d = fill(D,1,1)
    nx, nu, ny = _sscheck(A, B, C, d, Ts)
    new{Val{:siso},Val{:disc},M1,M2,M3,typeof(d)}(A, B, C, d, nx, nu, ny,
      convert(Float64, Ts))
  end

  # Continuous-time, multi-input-multi-output state-space model
  function (::Type{StateSpace}){M1<:AbstractMatrix,M2<:AbstractMatrix,
    M3<:AbstractMatrix,M4<:AbstractMatrix}(A::M1, B::M2, C::M3, D::M4)
    nx, nu, ny = _sscheck(A, B, C, D)
    new{Val{:mimo},Val{:cont},M1,M2,M3,M4}(A, B, C, D, nx, nu, ny, zero(Float64))
  end

  # Discrete-time, multi-input-multi-output state-space model
  function (::Type{StateSpace}){M1<:AbstractMatrix,M2<:AbstractMatrix,
    M3<:AbstractMatrix,M4<:AbstractMatrix}(A::M1, B::M2, C::M3, D::M4, Ts::Real)
    nx, nu, ny = _sscheck(A, B, C, D, Ts)
    new{Val{:mimo},Val{:disc},M1,M2,M3,M4}(A, B, C, D, nx, nu, ny,
      convert(Float64, Ts))
  end

  # Function calls
  ## Time response
  function (sys::StateSpace)(t::Real, x::AbstractVector, dx::AbstractVector, u)
    A, B, C, D = sys.A, sys.B, sys.C, sys.D
    ucalc = u(t, x)
    ucalc = isa(ucalc, Real) ? [ucalc] : ucalc
    if !isa(ucalc, AbstractVector) || length(ucalc) ≠ sys.nu || !(eltype(ucalc) <: Real)
      warn("sys(t,x,dx,u): u(t,x) has to be an `AbstractVector` of length $(sys.nu), containing `Real` values")
      throw(DomainError())
    end
    dx[:] = A*x + B*ucalc
  end
  function (sys::StateSpace)(t::Real, x::DiffEqBase.DEDataArray, dx::AbstractVector, u)
    A, B, C, D = sys.A, sys.B, sys.C, sys.D
    ucalc = u(t, x.x)
    ucalc = isa(ucalc, Real) ? [ucalc] : ucalc
    if !isa(ucalc, AbstractVector) || length(ucalc) ≠ sys.nu || !(eltype(ucalc) <: Real)
      warn("sys(t,x,dx,u): u(t,x) has to be an `AbstractVector` of length $(sys.nu), containing `Real` values")
      throw(DomainError())
    end
    x.u   = ucalc
    x.y   = C*x.x + D*x.u
    dx[:] = A*x.x + B*x.u
  end

  ## Frequency response
  (sys::StateSpace{Val{:siso}})(x::Number)                                    =
    _eval(sys, x)[1]
  (sys::StateSpace{Val{:siso}}){M<:Number}(X::AbstractArray{M})               =
    reshape(_eval(sys, X), size(X))
  (sys::StateSpace{Val{:mimo}})(x::Number)                                    =
    _eval(sys, x)
  (sys::StateSpace{Val{:mimo}}){M<:Number}(X::AbstractArray{M})               =
    _eval(sys, X)
  function (sys::StateSpace){T<:Real}(; ω::Union{T, AbstractArray{T}} = Float64[])
    if isempty(ω)
      warn("sys(): Provide an argument for the function call. Refer to `?freqresp`.")
      throw(DomainError())
    end
    freqresp(sys, ω)
  end
end

# Enforce state-space type invariance
function _sscheck(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix,
  D::AbstractMatrix, Ts::Real = zero(Float64))
  na, ma  = size(A)
  nb, mb  = size(B)
  nc, mc  = size(C)
  nd, md  = size(D)

  if Ts < zero(Ts) || isinf(Ts)
    warn("StateSpace: Ts must be non-negative real number")
    throw(DomainError())
  end

  if na ≠ ma || !(eltype(A) <: Real)
    warn("StateSpace: A must be a square matrix of real numbers")
    throw(DomainError())
  end

  if !(eltype(B) <: Real)
    warn("StateSpace: B must be a matrix of real numbers")
    throw(DomainError())
  end

  if !(eltype(C) <: Real)
    warn("StateSpace: C must be a matrix of real numbers")
    throw(DomainError())
  end

  if !(eltype(D) <: Real)
    warn("StateSpace: D must be a matrix of real numbers")
    throw(DomainError())
  end

  if na ≠ nb
    warn("StateSpace: A and B must have the same number of rows")
    throw(DomainError())
  end

  if ma ≠ mc
    warn("StateSpace: A and C must have the same number of columns")
    throw(DomainError())
  end

  if nc ≠ nd || nc < 1
    warn("StateSpace: C and D must have the same number (≥1) of rows")
    throw(DomainError())
  end

  if mb ≠ md || mb < 1
    warn("StateSpace: B and D must have the same number (≥1) of columns")
    throw(DomainError())
  end

  return na, mb, nc
end

# Outer constructors
# SISO
ss(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::Real = zero(Float64)) =
  StateSpace(A, B, C, D)

ss(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::Real, Ts::Real)    =
  StateSpace(A, B, C, D, Ts)

ss{T<:Real}(D::T)           = StateSpace(zeros(T,0,0), zeros(T,0,1), zeros(T,1,0), D)
ss{T<:Real}(D::T, Ts::Real) = StateSpace(zeros(T,0,0), zeros(T,0,1), zeros(T,1,0), D, Ts)

# MIMO
ss(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::AbstractMatrix)    =
  StateSpace(A, B, C, D)

function ss(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix,
  D::AbstractVector)
  @assert isempty(D) "ss(A,B,C,D): D can only be an empty vector"
  d = spzeros(Float64, size(C,1), size(B,2))
  StateSpace(A, B, C, d)
end

ss(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::AbstractMatrix,
  Ts::Real) = StateSpace(A, B, C, D, Ts)

function ss(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix,
  D::AbstractVector, Ts::Real)
  @assert isempty(D) "ss(A,B,C,D): D can only be an empty vector"
  d = spzeros(Float64, size(C,1), size(B,2))
  StateSpace(A, B, C, d, Ts)
end

ss{T<:Real}(D::AbstractMatrix{T})           = StateSpace(zeros(T,0,0),
  zeros(T,0,size(D,2)), zeros(T,size(D,1),0), D)
ss{T<:Real}(D::AbstractMatrix{T}, Ts::Real) = StateSpace(zeros(T,0,0),
  zeros(T,0,size(D,2)), zeros(T,size(D,1),0), D, Ts)

# # Catch-all for convenience when dealing with scalars, vectors, etc.
# function _reshape(A::Union{Real,VecOrMat}, B::Union{Real,VecOrMat},
#   C::Union{Real,VecOrMat})
#   a = isa(A, Real) ? fill(A,1,1) : reshape(A, size(A,1,2)...)
#   b = isa(B, Real) ? fill(B,1,1) : reshape(B, size(B,1,2)...)
#   c = isa(C, Real) ? fill(C,1,1) : reshape(C, size(C,1,2)...)
#   return a, b, c
# end
# ss(A::Union{Real,VecOrMat}, B::Union{Real,VecOrMat},
#   C::Union{Real,VecOrMat})              = ss(_reshape(A, B, C)...)
# ss(A::Union{Real,VecOrMat}, B::Union{Real,VecOrMat},
#   C::Union{Real,VecOrMat}, D)           = ss(_reshape(A, B, C)..., D)
# ss(A::Union{Real,VecOrMat}, B::Union{Real,VecOrMat},
#   C::Union{Real,VecOrMat}, D, Ts::Real) = ss(_reshape(A, B, C)..., D, Ts)

# Interfaces
samplingtime(s::StateSpace) = s.Ts

numstates(s::StateSpace)    = s.nx
numinputs(s::StateSpace)    = s.nu
numoutputs(s::StateSpace)   = s.ny

# Iteration interface
start(s::StateSpace{Val{:mimo}})        = start(s.D)
next(s::StateSpace{Val{:mimo}}, state)  = (s[state], state+1)
done(s::StateSpace{Val{:mimo}}, state)  = done(s.D, state)

eltype{S,M1}(::Type{StateSpace{Val{:mimo},Val{S},M1}}) =
  StateSpace{Val{:siso},Val{S},M1}

length(s::StateSpace{Val{:mimo}}) = length(s.D)
size(s::StateSpace)               = size(s.D)
size(s::StateSpace, d)            = size(s.D, d)

# Indexing of MIMO systems
function getindex(s::StateSpace{Val{:mimo},Val{:cont}}, row::Int, col::Int)
  (1 ≤ row ≤ s.ny && 1 ≤ col ≤ s.nu) || throw(BoundsError(s.D, (row,col)))
  StateSpace(s.A, view(s.B, :, col:col), view(s.C, row:row, :), s.D[row, col])
end

function getindex(s::StateSpace{Val{:mimo},Val{:disc}}, row::Int, col::Int)
  (1 ≤ row ≤ s.ny && 1 ≤ col ≤ s.nu) || throw(BoundsError(s.D, (row,col)))
  StateSpace(s.A, view(s.B, :, col:col), view(s.C, row:row, :), s.D[row, col], s.Ts)
end

function getindex(s::StateSpace{Val{:mimo}}, idx::Int)
  (1 ≤ idx ≤ length(s.D)) || throw(BoundsError(s.D, idx))
  col, row  = divrem(idx-1, s.ny)
  s[row+1, col+1]
end

function getindex(s::StateSpace{Val{:mimo},Val{:cont}}, rows::AbstractVector{Int},
  cols::AbstractVector{Int})
  1 ≤ minimum(rows) ≤ maximum(rows) ≤ s.ny || throw(BoundsError(s.D, rows))
  1 ≤ minimum(cols) ≤ maximum(cols) ≤ s.nu || throw(BoundsError(s.D, cols))

  StateSpace(s.A, view(s.B, :, cols), view(s.C, rows, :), view(s.D, rows, cols))
end

function getindex(s::StateSpace{Val{:mimo},Val{:disc}}, rows::AbstractVector{Int},
  cols::AbstractVector{Int})
  1 ≤ minimum(rows) ≤ maximum(rows) ≤ s.ny || throw(BoundsError(s.D, rows))
  1 ≤ minimum(cols) ≤ maximum(cols) ≤ s.nu || throw(BoundsError(s.D, cols))

  StateSpace(s.A, view(s.B, :, cols), view(s.C, rows, :), view(s.D, rows, cols), s.Ts)
end

function getindex(s::StateSpace{Val{:mimo}}, indices::AbstractVector{Int})
  1 ≤ minimum(indices) ≤ maximum(indices) ≤ length(s.D) || throw(BoundsError(s.D, indices))

  temp  = map(x->divrem(x-1, s.ny), indices)
  cols  = map(x->x[1]+1, temp)
  rows  = map(x->x[2]+1, temp)

  s[rows, cols]
end

getindex(s::StateSpace{Val{:mimo}}, rows, ::Colon)    = s[rows, 1:s.nu]
getindex(s::StateSpace{Val{:mimo}}, ::Colon, cols)    = s[1:s.ny, cols]
getindex(s::StateSpace{Val{:mimo}}, ::Colon)          = s[1:end]
getindex(s::StateSpace{Val{:mimo}}, ::Colon, ::Colon) = s[1:s.ny,1:s.nu]
endof(s::StateSpace{Val{:mimo}})                      = endof(s.D)

# Multiplicative and additive identities (meaningful only for SISO)
one{M1,M2,M3,M4}(::Type{StateSpace{Val{:siso},Val{:cont},M1,M2,M3,M4}})   =
  StateSpace(zeros(eltype(M1),0,0), zeros(eltype(M2),0,1), zeros(eltype(M3),1,0),
  one(eltype(M4)))
one{M1,M2,M3,M4}(::Type{StateSpace{Val{:siso},Val{:disc},M1,M2,M3,M4}})   =
  StateSpace(zeros(eltype(M1),0,0), zeros(eltype(M2),0,1), zeros(eltype(M3),1,0),
  one(eltype(M4)), zero(Float64))
zero{M1,M2,M3,M4}(::Type{StateSpace{Val{:siso},Val{:cont},M1,M2,M3,M4}})  =
  StateSpace(zeros(eltype(M1),0,0), zeros(eltype(M2),0,1), zeros(eltype(M3),1,0),
  zero(eltype(M4)))
zero{M1,M2,M3,M4}(::Type{StateSpace{Val{:siso},Val{:disc},M1,M2,M3,M4}})  =
  StateSpace(zeros(eltype(M1),0,0), zeros(eltype(M2),0,1), zeros(eltype(M3),1,0),
  zero(eltype(M4)), zero(Float64))

one(s::StateSpace{Val{:siso},Val{:cont}})   = one(typeof(s))
one(s::StateSpace{Val{:siso},Val{:disc}})   = one(typeof(s))
zero(s::StateSpace{Val{:siso},Val{:cont}})  = zero(typeof(s))
zero(s::StateSpace{Val{:siso},Val{:disc}})  = zero(typeof(s))

# Inverse of a state-space model
function _ssinv(s::StateSpace)
  if s.ny ≠ s.nu
    warn("inv(sys): s.ny ≠ s.nu")
    throw(DomainError())
  end

  try
    Dinv = inv(s.D);
    Ainv = s.A - s.B*Dinv*s.C;
    Binv = s.B*Dinv
    Cinv = -Dinv*s.C
    return Ainv, Binv, Cinv, Dinv
  catch err
    warn("inv(sys): sys is not invertible")
    throw(DomainError())
  end
end

function inv(s::StateSpace{Val{:siso},Val{:cont}})
  Ainv, Binv, Cinv, Dinv = _ssinv(s)
  StateSpace(Ainv, Binv, Cinv, Dinv[1])
end

function inv(s::StateSpace{Val{:siso},Val{:disc}})
  Ainv, Binv, Cinv, Dinv = _ssinv(s)
  StateSpace(Ainv, Binv, Cinv, Dinv[1], s.Ts)
end

function inv(s::StateSpace{Val{:mimo},Val{:cont}})
  Ainv, Binv, Cinv, Dinv = _ssinv(s)
  StateSpace(Ainv, Binv, Cinv, Dinv)
end

function inv(s::StateSpace{Val{:mimo},Val{:disc}})
  Ainv, Binv, Cinv, Dinv = _ssinv(s)
  StateSpace(Ainv, Binv, Cinv, Dinv, s.Ts)
end

# Invariant zeros of a state-space model
function zeros(s::StateSpace)
  Ar, Br, Cr, Dr, mr, nr, pr        = reduce(s.A, s.B, s.C, s.D)
  if nr == 0
    return Complex{Float64}[]
  end
  Arc, Brc, Crc, Drc, mrc, nrc, prc = reduce(Ar.', Cr.', Br.', Dr.')
  if nrc == 0
    return Complex{Float64}[]
  end

  svdobj  = svdfact([Crc Drc], thin = false)
  W       = flipdim(svdobj.Vt', 2)
  Af      = [Arc Brc]*W[:, 1:nrc]

  if mrc == 0
    zerovalues = eigfact(Af).values
    return zerovalues
  else
    Bf    = W[1:nrc,1:nrc]
    zerovalues = eigfact(Af, Bf).values
    return zerovalues
  end
end

# Transmission zeros of a state-space model
tzeros(s::StateSpace) = zeros(minreal(s))

# Poles of a state-space model
function poles(s::StateSpace)
  Aₘ, _, _, _ = minreal(s.A, s.B, s.C, s.D)
  return eigfact(Aₘ).values
end

# Negative of a state-space model
-(s::StateSpace{Val{:siso},Val{:cont}}) = StateSpace(s.A, s.B, -s.C, -s.D[1])
-(s::StateSpace{Val{:siso},Val{:disc}}) = StateSpace(s.A, s.B, -s.C, -s.D[1], s.Ts)
-(s::StateSpace{Val{:mimo},Val{:cont}}) = StateSpace(s.A, s.B, -s.C, -s.D)
-(s::StateSpace{Val{:mimo},Val{:disc}}) = StateSpace(s.A, s.B, -s.C, -s.D, s.Ts)

# Addition
function _ssparallel{T1,T2,S}(s1::StateSpace{Val{T1},Val{S}},
  s2::StateSpace{Val{T2},Val{S}})
  if s1.Ts ≉ s2.Ts && s1.Ts ≠ zero(s1.Ts) && s2.Ts ≠ zero(s2.Ts)
    warn("parallel(s1,s2): Sampling time mismatch")
    throw(DomainError())
  end

  if size(s1) ≠ size(s2)
    warn("parallel(s1,s2): size(s1) ≠ size(s2)")
    throw(DomainError())
  end

  T = promote_type(eltype(s1.A), eltype(s2.A))
  a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B, s2.B)
  c = hcat(s1.C, s2.C)
  d = s1.D + s2.D

  return a, b, c, d, max(s1.Ts, s2.Ts)
end

function +(s1::StateSpace{Val{:siso},Val{:cont}},
  s2::StateSpace{Val{:siso},Val{:cont}})
  a, b, c, d, _ = _ssparallel(s1, s2)
  StateSpace(a, b, c, d[1])
end

function +(s1::StateSpace{Val{:siso},Val{:disc}},
  s2::StateSpace{Val{:siso},Val{:disc}})
  a, b, c, d, Ts = _ssparallel(s1, s2)
  StateSpace(a, b, c, d[1], Ts)
end

function +{T1,T2}(s1::StateSpace{Val{T1},Val{:cont}},
  s2::StateSpace{Val{T2},Val{:cont}})
  a, b, c, d, _ = _ssparallel(s1, s2)
  StateSpace(a, b, c, d)
end

function +{T1,T2}(s1::StateSpace{Val{T1},Val{:disc}},
  s2::StateSpace{Val{T2},Val{:disc}})
  a, b, c, d, Ts = _ssparallel(s1, s2)
  StateSpace(a, b, c, d, Ts)
end

.+(s1::StateSpace{Val{:siso}}, s2::StateSpace{Val{:siso}}) = +(s1, s2)

+{T}(s::StateSpace{Val{T},Val{:disc}}, g::Union{Real,AbstractMatrix}) =
  +(s, ss(g, samplingtime(s)))
+{T}(g::Union{Real,AbstractMatrix}, s::StateSpace{Val{T},Val{:disc}}) =
  +(ss(g, samplingtime(s)), s)
+{T}(s::StateSpace{Val{T},Val{:cont}}, g::Union{Real,AbstractMatrix}) =
  +(s, ss(g))
+{T}(g::Union{Real,AbstractMatrix}, s::StateSpace{Val{T},Val{:cont}}) =
  +(ss(g), s)

.+(s::StateSpace{Val{:siso}}, g::Real)    = +(s, g)
.+(g::Real, s::StateSpace{Val{:siso}})    = +(g, s)

# Subtraction
-(s1::StateSpace, s2::StateSpace) = +(s1, -s2)

.-(s1::StateSpace{Val{:siso}}, s2::StateSpace{Val{:siso}}) = -(s1, s2)

-(s::StateSpace, g::Union{Real,AbstractMatrix}) = +(s, -g)
-(g::Union{Real,AbstractMatrix}, s::StateSpace) = +(g, -s)

.-(s::StateSpace{Val{:siso}}, g::Real)    = -(s, g)
.-(g::Real, s::StateSpace{Val{:siso}})    = -(g, s)

# Multiplication
function _ssseries{T1,T2,S}(s1::StateSpace{Val{T1},Val{S}},
  s2::StateSpace{Val{T2},Val{S}})
  # Remark: s1*s2 implies u -> s2 -> s1 -> y

  if s1.Ts ≉ s2.Ts && s1.Ts ≠ zero(s1.Ts) && s2.Ts ≠ zero(s2.Ts)
    warn("series(s1,s2): Sampling time mismatch")
    throw(DomainError())
  end

  if s1.nu ≠ s2.ny
    warn("series(s1,s2): s1.nu ≠ s2.ny")
    throw(DomainError())
  end

  T = promote_type(eltype(s1.A), eltype(s1.B), eltype(s2.A), eltype(s2.C))

  a = vcat(hcat(s1.A, s1.B*s2.C),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B*s2.D, s2.B)
  c = hcat(s1.C, s1.D*s2.C)
  d = s1.D * s2.D

  return a, b, c, d, max(s1.Ts, s2.Ts)
end

function *(s1::StateSpace{Val{:siso},Val{:cont}},
  s2::StateSpace{Val{:siso},Val{:cont}})
  a, b, c, d, _ = _ssseries(s1, s2)
  StateSpace(a, b, c, d[1])
end

function *(s1::StateSpace{Val{:siso},Val{:disc}},
  s2::StateSpace{Val{:siso},Val{:disc}})
  a, b, c, d, Ts = _ssseries(s1, s2)
  StateSpace(a, b, c, d[1], Ts)
end

function *{T1,T2}(s1::StateSpace{Val{T1},Val{:cont}},
  s2::StateSpace{Val{T2},Val{:cont}})
  a, b, c, d, _ = _ssseries(s1, s2)
  StateSpace(a, b, c, d)
end

function *{T1,T2}(s1::StateSpace{Val{T1},Val{:disc}},
  s2::StateSpace{Val{T2},Val{:disc}})
  a, b, c, d, Ts = _ssseries(s1, s2)
  StateSpace(a, b, c, d, Ts)
end

.*(s1::StateSpace{Val{:siso}}, s2::StateSpace{Val{:siso}}) = *(s1, s2)

*{T}(s::StateSpace{Val{T},Val{:disc}}, g::Union{Real,AbstractMatrix}) =
  *(s, ss(g, samplingtime(s)))
*{T}(g::Union{Real,AbstractMatrix}, s::StateSpace{Val{T},Val{:disc}}) =
  *(ss(g, samplingtime(s)), s)
*{T}(s::StateSpace{Val{T},Val{:cont}}, g::Union{Real,AbstractMatrix}) =
  *(s, ss(g))
*{T}(g::Union{Real,AbstractMatrix}, s::StateSpace{Val{T},Val{:cont}}) =
  *(ss(g), s)

.*(s::StateSpace{Val{:siso}}, g::Real)    = *(s, g)
.*(g::Real, s::StateSpace{Val{:siso}})    = *(g, s)

# Division
/(s1::StateSpace, s2::StateSpace)         = *(s1, inv(s2))

./(s1::StateSpace{Val{:siso}}, s2::StateSpace{Val{:siso}}) = /(s1, s2)

/(s::StateSpace, g::Union{Real,AbstractMatrix}) =
  *(s, inv(g))
/(g::Union{Real,AbstractMatrix}, s::StateSpace) =
  *(g, inv(s))

./(s::StateSpace{Val{:siso}}, g::Real)    = /(s, g)
./(g::Real, s::StateSpace{Val{:siso}})    = /(g, s)
