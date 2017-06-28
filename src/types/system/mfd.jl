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
# NOTE: Should SISO MatrixFractionDescriptions be based on 1x1 polynomial matrices (similarly to `TransferFunction`)?
# NOTE: Should the constructors verify that D is nonsingular and that the MatrixFractionDescription is proper?
immutable MatrixFractionDescription{T,S,L,M1,M2}  <: LtiSystem{T,S}
  N::M1
  D::M2
  nu::Int
  ny::Int
  Ts::Float64

  # Continuous-time, single-input-single-output MatrixFractionDescription model
  function (::Type{MatrixFractionDescription}){L,M1<:Polynomials.Poly,M2<:Polynomials.Poly}(N::M1, D::M2, ::Type{Val{L}})
    mfdcheck(N,D)
    nN,nD = length(N),length(D)
    _N = PolyMatrix(reshape(coeffs(N),1,1,nN), (1,1,nN), Val{:s})
    _D = PolyMatrix(reshape(coeffs(D),1,1,nD), (1,1,nD), Val{:s})
    new{Val{:siso},Val{:cont},Val{L},typeof(_N),typeof(_D)}(_N, _D, 1, 1, zero(Float64))
  end

  # Discrete-time, single-input-single-output MatrixFractionDescription model
  function (::Type{MatrixFractionDescription}){L,M1<:Polynomials.Poly,M2<:Polynomials.Poly}(N::M1, D::M2, Ts::Real, ::Type{Val{L}})
    mfdcheck(N,D,Ts)
    nN,nD = length(N),length(D)
    _N = PolyMatrix(reshape(coeffs(N),1,1,nN), (1,1,nN), Val{:s})
    _D = PolyMatrix(reshape(coeffs(D),1,1,nD), (1,1,nD), Val{:s})
    new{Val{:siso},Val{:disc},Val{L},typeof(_N),typeof(_D)}(_N, _D, 1, 1, convert(Float64, Ts))
  end

  # Continuous-time, multi-input-multi-output MatrixFractionDescription model
  function (::Type{MatrixFractionDescription}){L,M1<:PolynomialMatrices.PolyMatrix,
    M2<:PolynomialMatrices.PolyMatrix}(N::M1, D::M2, ::Type{Val{L}})
    ny, nu = mfdcheck(N, D, Val{L})
    _N = PolyMatrix(coeffs(N), size(N), Val{:s})
    _D = PolyMatrix(coeffs(D), size(D), Val{:s})
    new{Val{:mimo},Val{:cont},Val{L},typeof(_N),typeof(_D)}(_N, _D, nu, ny,
                                                            zero(Float64))
    # should we do better checks than just converting the variable to the correct one?
  end

  # Discrete-time, multi-input-multi-output MatrixFractionDescription model
  function (::Type{MatrixFractionDescription}){L,M1<:PolynomialMatrices.PolyMatrix,
    M2<:PolynomialMatrices.PolyMatrix}(N::M1, D::M2, Ts::Real, ::Type{Val{L}})
    ny, nu = mfdcheck(N, D, Val{L}, Ts)
    _N = PolyMatrix(coeffs(N), size(N), Val{:z})
    _D = PolyMatrix(coeffs(D), size(D), Val{:z})
    new{Val{:mimo},Val{:disc},Val{L},typeof(_N),typeof(_D)}(_N, _D, nu, ny,
                                                            convert(Float64, Ts))
    # should we do better checks than just converting the variable to the correct one?
  end

  # Function calls
  # TODO needs to be properly set up according to e.g. Kailath Linear systems 6.4
  # to properly handle all column/row reduced cases
  function (sys::MatrixFractionDescription{T,S,Val{:lfd}}){T,S}(t::Real, x::DiffEqBase.DEDataArray, dx::AbstractVector, u)
    ucalc = u(t, x.y)
    ucalc = isa(ucalc, Real) ? [ucalc] : ucalc
    if !isa(ucalc, AbstractVector) || length(ucalc) ≠ sys.nu || !(eltype(ucalc) <: Real)
      warn("sys(t,x,dx,u): u(t,y) has to be an `AbstractVector` of length $(sys.nu), containing `Real` values")
      throw(DomainError())
    end

    ny  = numoutputs(sys)
    N   = num(sys)
    D   = den(sys)

    degs, Dhc = high_row_deg_matrix(D)
    Dhci = inv(Dhc)
    D = Dhci*D

    xidx = (cumsum(vcat(0,degs))+1)[1:end-1]
    dxidx = cumsum(degs,1)[:]

    dx[:] = zeros(dx)
    # "A matrix"
    for (k,v) in coeffs(D)
      for row in find(degs .> k)
        dx[dxidx[row]-k] += -dot(v[row,:], x.x[xidx])
      end
    end

    for row in 1:ny
      nrowx = degs[row]-1
      dx[xidx[row]-1+(1:nrowx)] += x.x[xidx[row]+(1:nrowx)]
    end

    # "B matrix"
    for (k,v) in coeffs(N)
      for row in find(degs .> k)
        dx[dxidx[row]-k] += dot(v[row,:], ucalc)
      end
    end

    # take care of direct term
    if degree(N) == degree(D)
      tmp = coeffs(N)[degree(N)]*ucalc
      for (k,v) in coeffs(D)
        for row in find(degs .> k)
          dx[dxidx[row]-k] += -dot(v[row,:], tmp[xidx])
        end
      end
    end

    # "C matrix"
    x.y = degree(N) == degree(D) ? coeffs(N)[degree(N)]*ucalc : zeros(ny)
    x.y += Dhci*x[xidx]

    x.u[:]  = ucalc
  end

  function (sys::MatrixFractionDescription{T,S,Val{:rfd}}){T,S}(t::Real, x::DiffEqBase.DEDataArray, dx::AbstractVector, u)
    ucalc = u(t, x.y)
    ucalc = isa(ucalc, Real) ? [ucalc] : ucalc
    if !isa(ucalc, AbstractVector) || length(ucalc) ≠ sys.nu || !(eltype(ucalc) <: Real)
      warn("sys(t,x,dx,u): u(t,y) has to be an `AbstractVector` of length $(sys.nu), containing `Real` values")
      throw(DomainError())
    end

    ny  = numoutputs(sys)
    nu  = numinputs(sys)
    N   = num(sys)
    D   = den(sys)

    degs, Dhc = high_col_deg_matrix(D)
    Dhci = inv(Dhc)
    D = Dhci*D

    xidx = cumsum(degs,2)[:]
    dxidx = (cumsum(hcat(0,degs),2)+1)[1:end-1]

    x.y = zeros(x.y)
    dx[:] = zeros(dx)

    # "A matrix"
    for (k,v) in coeffs(D)
      for col in find(degs .> k)
        dx[dxidx] += -v[:,col]*x.x[xidx[col]-k]
      end
    end
    for col in 1:nu
      ncolx = degs[col]-1
      dx[dxidx[col]+(1:ncolx)] += x.x[dxidx[col]-1+(1:ncolx)]
    end

    # "B matrix"
    dx[dxidx] += Dhci*ucalc

    # "C matrix"
    for (k,v) in coeffs(N)
      for col in find(degs .> k)
        x.y[:] += v[:,col]*x.x[xidx[col]-k]
      end
    end
    x.u[:]  = ucalc

    # take care of direct term
    if degree(N) == degree(D)
      b0 = coeffs(N)[degree(N)]
      tmp = ucalc
      for (k,v) in coeffs(D)
        for col in find(degs .> k)
          tmp += -v[:,col]*x.x[xidx[col]-k]
        end
      end
      x.y += b0*tmp
    end
  end
end

function mfdcheck{T<:Real,S<:Real}(N::Poly{T}, D::Poly{S}, Ts::Real = zero(Float64))
  @assert Ts ≥ zero(Ts) && !isinf(Ts) "MatrixFractionDescription: Ts must be non-negative real number"
end

isproper(s::MatrixFractionDescription{Val{:siso}})  = degree(s.D) ≥ degree(s.N)

function isproper{S}(s::MatrixFractionDescription{Val{:mimo},Val{S},Val{:rfd}})
  if is_col_proper(s.D)
    return all(col_degree(s.D) .>= col_degree(s.N))
  else
    isproper(tf(s))
  end
end

function isproper{S}(s::MatrixFractionDescription{Val{:mimo},Val{S},Val{:lfd}})
  if is_row_proper(s.D)
    return all(row_degree(s.D) .>= row_degree(s.N))
  else
    isproper(tf(s))
  end
end

isstrictlyproper(s::MatrixFractionDescription{Val{:siso}}) = degree(s.D) > degree(s.N)

function isstrictlyproper{S}(s::MatrixFractionDescription{Val{:mimo},Val{S},Val{:rfd}})
  if is_col_proper(s.D)
    return all(col_degree(s.D) .>= col_degree(s.N))
  else
    isproper(tf(s))
  end
end

function isstrictlyproper{S}(s::MatrixFractionDescription{Val{:mimo},Val{S},Val{:lfd}})
  if is_row_proper(s.D)
    return all(row_degree(s.D) .>= row_degree(s.N))
  else
    isstrictlyproper(tf(s))
  end
end

# Enforce rational transfer function type invariance
function mfdcheck{M1<:PolynomialMatrices.PolyMatrix,M2<:PolynomialMatrices.PolyMatrix}(
  N::M1, D::M2, ::Type{Val{:lfd}})
  @assert size(N,1) == size(D,1) "MatrixFractionDescription: size(N,1) ≠ size(D,1)"
  @assert size(D,1) == size(D,2) "MatrixFractionDescription: size(D,1) ≠ size(D,2)"

  return size(N,1), size(N,2)
end
function mfdcheck{M1<:PolynomialMatrices.PolyMatrix,M2<:PolynomialMatrices.PolyMatrix}(
  N::M1, D::M2, ::Type{Val{:rfd}})
  @assert size(N,2) == size(D,2) "MatrixFractionDescription: size(N,2) ≠ size(D,2)"
  @assert size(D,1) == size(D,2) "MatrixFractionDescription: size(D,1) ≠ size(D,2)"

  return size(N,1), size(N,2)
end
function mfdcheck{T,M1<:PolynomialMatrices.PolyMatrix,M2<:PolynomialMatrices.PolyMatrix}(
  N::M1, D::M2, ::Type{Val{T}}, Ts::Real)
  @assert Ts ≥ zero(Ts) && !isinf(Ts) "MatrixFractionDescription: Ts must be non-negative real number"
  return mfdcheck(N, D, Val{T})
end

# Outer constructors
lfd(N::Poly, D::Poly)           = MatrixFractionDescription(N, D, Val{:lfd})
lfd(N::Poly, D::Poly, Ts::Real) = MatrixFractionDescription(N, D, convert(Float64, Ts), Val{:lfd})
lfd{M1<:PolyMatrix,M2<:PolyMatrix}(N::M1, D::M2) =
  MatrixFractionDescription(N, D, Val{:lfd})
lfd{M1<:PolyMatrix,M2<:PolyMatrix}(N::M1, D::M2, Ts::Real) =
  MatrixFractionDescription(N, D, convert(Float64, Ts), Val{:lfd})

rfd(N::Poly, D::Poly) = MatrixFractionDescription(N, D, Val{:rfd})
rfd(N::Poly, D::Poly, Ts::Real) = MatrixFractionDescription(N, D, convert(Float64, Ts), Val{:rfd})
rfd(N::PolyMatrix, D::PolyMatrix) = MatrixFractionDescription(N, D, Val{:rfd})
rfd(N::PolyMatrix, D::PolyMatrix, Ts::Real) = MatrixFractionDescription(N, D, convert(Float64, Ts), Val{:rfd})

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
samplingtime(s::MatrixFractionDescription)              = s.Ts
islfd{T,S,L}(s::MatrixFractionDescription{T,S,Val{L}})  = false
islfd{T,S}(s::MatrixFractionDescription{T,S,Val{:lfd}}) = true
isrfd{T,S,L}(s::MatrixFractionDescription{T,S,Val{L}})  = !islfd(s)
num(s::MatrixFractionDescription)                       = s.N
den(s::MatrixFractionDescription)                       = s.D

# Think carefully about how to implement numstates
numstates{T,S}(s::MatrixFractionDescription{T,S,Val{:lfd}}) = sum(row_degree(s.D))
numstates{T,S}(s::MatrixFractionDescription{T,S,Val{:rfd}}) = sum(col_degree(s.D))
# Currently, we only allow for proper systems
numinputs(s::MatrixFractionDescription)               = s.nu
numoutputs(s::MatrixFractionDescription)              = s.ny

# Dimension information
size(s::MatrixFractionDescription)                    = size(s.N)
size(s::MatrixFractionDescription, d)                 = size(s.N, d)
length(s::MatrixFractionDescription{Val{:mimo}})      = length(s.N)

# conversion between 1×1 mimo and siso
function siso{L}(s::MatrixFractionDescription{Val{:mimo},Val{:cont},Val{L}})
  if size(s) != (1,1)
    warn("siso(s): system is not 1×1")
    throw(DomainError())
  end
  MatrixFractionDescription(s.N[1], s.D[1], Val{L})
end

function siso{L}(s::MatrixFractionDescription{Val{:mimo},Val{:disc},Val{L}})
  if size(s) != (1,1)
    warn("siso(s): system is not 1×1")
    throw(DomainError())
  end
  MatrixFractionDescription(s.N[1], s.D[1], s.Ts, Val{L})
end

function mimo{L}(s::MatrixFractionDescription{Val{:siso},Val{:cont},Val{L}})
  MatrixFractionDescription(PolyMatrix(s.N, (1,1), Val{:s}) , PolyMatrix(s.D, (1,1), Val{:s}), Val{L})
end

function mimo{L}(s::MatrixFractionDescription{Val{:siso},Val{:disc},Val{L}})
  MatrixFractionDescription(PolyMatrix(s.N, (1,1), Val{:z}) , PolyMatrix(s.D, (1,1), Val{:z}), s.Ts, Val{L})
end

# ## Iteration interface
# start(s::MatrixFractionDescription{Val{:mimo}})       = start(s.N)
# next(s::MatrixFractionDescription{Val{:mimo}}, state) = (s[state], state+1)
# done(s::MatrixFractionDescription{Val{:mimo}}, state) = done(s.N, state)
#
# eltype{S,M1}(::Type{MatrixFractionDescription{Val{:mimo},Val{S},M1}}) =
#   MatrixFractionDescription{Val{:siso},Val{S},M1}
#
# # Indexing of MIMO systems
# function getindex(s::MatrixFractionDescription{Val{:mimo},Val{:cont},Val{:lfd}}, I...)
#   @boundscheck checkbounds(s.N, I...)
#   lfd(ss(s)[I...])
# end
#
# function getindex(s::MatrixFractionDescription{Val{:mimo},Val{:cont},Val{:rfd}}, I...)
#   @boundscheck checkbounds(s.N, I...)
#   rfd(ss(s)[I...])
# end
#
# endof(s::MatrixFractionDescription{Val{:mimo}}) = endof(s.N)
#
# Conversion and promotion
promote_rule{T<:Real,S,L}(::Type{T}, ::Type{MatrixFractionDescription{Val{:siso},S,L}}) =
  MatrixFractionDescription{Val{:siso},S,L}
promote_rule{T<:AbstractMatrix,S,L}(::Type{T}, ::Type{MatrixFractionDescription{Val{:mimo},S,L}}) =
  MatrixFractionDescription{Val{:mimo},S,L}

convert(::Type{MatrixFractionDescription{Val{:siso},Val{:cont},Val{:lfd}}}, g::Real)            =
  lfd(Poly(g,:s), Poly(one(g),:s))
convert(::Type{MatrixFractionDescription{Val{:siso},Val{:disc},Val{:lfd}}}, g::Real)            =
  lfd(Poly(g,:z), Poly(one(g),:z), zero(Float64))
convert(::Type{MatrixFractionDescription{Val{:mimo},Val{:cont},Val{:lfd}}}, g::AbstractMatrix)  =
  lfd(PolyMatrix(g, Val{:s}), PolyMatrix(eye(eltype(g), size(g,1), size(g,1))))
convert(::Type{MatrixFractionDescription{Val{:mimo},Val{:disc},Val{:lfd}}}, g::AbstractMatrix)  =
  lfd(PolyMatrix(g, Val{:z}), PolyMatrix(eye(eltype(g), size(g,1), size(g,1))), zero(Float64))

convert(::Type{MatrixFractionDescription{Val{:siso},Val{:cont},Val{:rfd}}}, g::Real)            =
  rfd(Poly(g,:s), Poly(one(g),:s))
convert(::Type{MatrixFractionDescription{Val{:siso},Val{:disc},Val{:rfd}}}, g::Real)            =
  rfd(Poly(g,:z), Poly(one(g),:z), zero(Float64))
convert(::Type{MatrixFractionDescription{Val{:mimo},Val{:cont},Val{:rfd}}}, g::AbstractMatrix)  =
  rfd(PolyMatrix(g, Val{:s}), PolyMatrix(eye(eltype(g), size(g,2), size(g,2))))
convert(::Type{MatrixFractionDescription{Val{:mimo},Val{:disc},Val{:rfd}}}, g::AbstractMatrix)  =
  rfd(PolyMatrix(g, Val{:z}), PolyMatrix(eye(eltype(g), size(g,2), size(g,2))), zero(Float64))

# conversions between lfd and rfd
convert{S,T,L}(::Type{MatrixFractionDescription{Val{S},Val{T},Val{:lfd}}},
  s::MatrixFractionDescription{Val{S},Val{T},Val{L}}) = lfd(s)
convert{S,T,L}(::Type{MatrixFractionDescription{Val{S},Val{T},Val{:rfd}}},
  s::MatrixFractionDescription{Val{S},Val{T},Val{L}}) = rfd(s)

# Multiplicative and additive identities (meaningful only for SISO)
one{M1,M2}(::Type{MatrixFractionDescription{Val{:siso},Val{:cont},Val{:lfd},M1,M2}})  =
  lfd(one(Int8))
one{M1,M2}(::Type{MatrixFractionDescription{Val{:siso},Val{:disc},Val{:lfd},M1,M2}})  =
  lfd(one(Int8), zero(Float64))
zero{M1,M2}(::Type{MatrixFractionDescription{Val{:siso},Val{:cont},Val{:lfd},M1,M2}}) =
  lfd(zero(Int8))
zero{M1,M2}(::Type{MatrixFractionDescription{Val{:siso},Val{:disc},Val{:lfd},M1,M2}}) =
  lfd(zero(Int8), zero(Float64))

one{M1,M2}(::Type{MatrixFractionDescription{Val{:siso},Val{:cont},Val{:rfd},M1,M2}})  =
  rfd(one(Int8))
one{M1,M2}(::Type{MatrixFractionDescription{Val{:siso},Val{:disc},Val{:rfd},M1,M2}})  =
  rfd(one(Int8), zero(Float64))
zero{M1,M2}(::Type{MatrixFractionDescription{Val{:siso},Val{:cont},Val{:rfd},M1,M2}}) =
  rfd(zero(Int8))
zero{M1,M2}(::Type{MatrixFractionDescription{Val{:siso},Val{:disc},Val{:rfd},M1,M2}}) =
  tf(zero(Int8), zero(Float64))

one{T,S,L,M1,M2}(s::MatrixFractionDescription{Val{T},Val{S},Val{L},M1,M2})   =
  one(typeof(s))
zero{T,S,L,M1,M2}(s::MatrixFractionDescription{Val{T},Val{S},Val{L},M1,M2})  =
  zero(typeof(s))

# conversions between lfd and rfd
lfd(s::MatrixFractionDescription{Val{:siso},Val{:cont},Val{:rfd}}) = lfd(s.N,s.D)
lfd(s::MatrixFractionDescription{Val{:siso},Val{:disc},Val{:rfd}}) = lfd(s.N,s.D,s.Ts)
lfd(s::MatrixFractionDescription{Val{:mimo},Val{:cont},Val{:rfd}}) = lfd(_rfd2lfd(s)...)
lfd(s::MatrixFractionDescription{Val{:mimo},Val{:disc},Val{:rfd}}) = lfd(_rfd2lfd(s)...,s.Ts)
lfd{T,S}(s::MatrixFractionDescription{T,S,Val{:lfd}})              = s

rfd(s::MatrixFractionDescription{Val{:siso},Val{:cont},Val{:lfd}}) = rfd(s.N,s.D)
rfd(s::MatrixFractionDescription{Val{:siso},Val{:disc},Val{:lfd}}) = rfd(s.N,s.D,s.Ts)
rfd(s::MatrixFractionDescription{Val{:mimo},Val{:cont},Val{:lfd}}) = rfd(_lfd2rfd(s)...)
rfd(s::MatrixFractionDescription{Val{:mimo},Val{:disc},Val{:lfd}}) = rfd(_lfd2rfd(s)...,s.Ts)
rfd{T,S}(s::MatrixFractionDescription{T,S,Val{:rfd}})              = s

function _lfd2rfd{T,S}(s::MatrixFractionDescription{T,S,Val{:lfd}})
  n,m = size(s)
  p   = hcat(-s.N, s.D)
  _,U = ltriang(p)

  D = U[1:m,n+1:n+m]
  N = U[m+1:n+m,n+1:n+m]
  return N,D
end

function _rfd2lfd{T,S}(s::MatrixFractionDescription{T,S,Val{:rfd}})
  n,m = size(s)
  p   = vcat(-s.D, s.N)
  _,U = rtriang(p)

  N = U[m+1:n+m,1:m]
  D = U[m+1:n+m,m+1:n+m]
  return N,D
end

function inv{M<:MatrixFractionDescription}(s::M)
  _mfdinvcheck(s)
  _inv(s)
end

_inv{T,L}(s::MatrixFractionDescription{Val{T},Val{:cont},Val{L}}) = MatrixFractionDescription(copy(s.D), copy(s.N), Val{L})
_inv{T,L}(s::MatrixFractionDescription{Val{T},Val{:disc},Val{L}}) = MatrixFractionDescription(copy(s.D), copy(s.N), s.Ts, Val{L})

function _mfdinvcheck(s::MatrixFractionDescription)
  if s.ny ≠ s.nu
    warn("inv(sys): s.ny ≠ s.nu")
    throw(DomainError())
  end

  if fastrank(s.N) ≠ s.nu
    warn("inv(sys): sys is not invertible")
    throw(DomainError())
  end
end

# Negative of a transfer-function model
-{T}(s::MatrixFractionDescription{Val{T},Val{:cont},Val{:lfd}}) = lfd(-s.N, copy(s.D))
-{T}(s::MatrixFractionDescription{Val{T},Val{:disc},Val{:lfd}}) = lfd(-s.N, copy(s.D), s.Ts)
-{T}(s::MatrixFractionDescription{Val{T},Val{:cont},Val{:rfd}}) = rfd(-s.N, copy(s.D))
-{T}(s::MatrixFractionDescription{Val{T},Val{:disc},Val{:rfd}}) = rfd(-s.N, copy(s.D), s.Ts)

# Addition (parallel)
function _mfdparallelcheck{T1,T2,S,L}(s₁::MatrixFractionDescription{Val{T1},Val{S},Val{L}},
  s₂::MatrixFractionDescription{Val{T2},Val{S},Val{L}})
  if s₁.Ts ≉ s₂.Ts && s₁.Ts ≠ zero(s₁.Ts) && s₂.Ts ≠ zero(s₂.Ts)
    warn("parallel(s₁,s₂): Sampling time mismatch")
    throw(DomainError())
  end

  if size(s₁,1) ≠ size(s₂,1)
    warn("parallel(s₁,s₂): size(s₁,1) ≠ size(s₂,1)")
    throw(DomainError())
  end
end

# siso version
function _mfdparallel{S,L}(s₁::MatrixFractionDescription{Val{:siso},Val{S},Val{:L}},
  s₂::MatrixFractionDescription{Val{:siso},Val{S},Val{L}})
  R   = gcd(s₁.D, s₂.D)
  D₁  = div(s₁.D, R)
  D   = D₁*s₂.D        # only include common part R once
  D₂  = div(s₂.D, R)
  N   = s₁.N*D₂ + s₂.N*D₁
  N, D, max(s₁.Ts, s₂.Ts)
end

# mimo lfd version
function _mfdparallel{S}(s₁::MatrixFractionDescription{Val{:mimo},Val{S},Val{:lfd}},
  s₂::MatrixFractionDescription{Val{:mimo},Val{S},Val{:lfd}})
  R, V₁, V₂ = gcrd(s₁.D, s₂.D)
  detV₁, adjV₁ = inv(V₁)
  detV₂, adjV₂ = inv(V₂)
  N₁ = adjV₁*s₁.N/detV₁(0)
  N₂ = adjV₂*s₂.N/detV₂(0)
  N₁+N₂, R, max(s₁.Ts, s₂.Ts)
end

# mimo rfd version
function _mfdparallel{S}(s₁::MatrixFractionDescription{Val{:mimo},Val{S},Val{:rfd}},
  s₂::MatrixFractionDescription{Val{:mimo},Val{S},Val{:rfd}})
  L, V₁, V₂ = gcld(s₁.D, s₂.D)
  detV₁, adjV₁ = inv(V₁)
  detV₂, adjV₂ = inv(V₂)
  N₁ = s₁.N*adjV₁/detV₁(0)
  N₂ = s₂.N*adjV₂/detV₂(0)
  N₁+N₂, L, max(s₁.Ts, s₂.Ts)
end

# mimo mixed lfd/rfd version
function _mfdparallel{S,L1,L2}(s₁::MatrixFractionDescription{Val{:mimo},Val{S},Val{L1}},
  s₂::MatrixFractionDescription{Val{:mimo},Val{S},Val{L2}})
  _mfdparallel(lfd(s₁), lfd(s₂))
end

# siso and mimo of dimensions 1×1
function _mfdparallel{T1,T2,S,L1,L2}(s₁::MatrixFractionDescription{Val{T1},Val{S},Val{L1}},
  s₂::MatrixFractionDescription{Val{T2},Val{S},Val{L2}})
  _mfdparallel(mimo(s₁), mimo(s₂))
end

function +{T1,T2,L1,L2}(s₁::MatrixFractionDescription{Val{T1},Val{:cont},Val{L1}},
  s₂::MatrixFractionDescription{Val{T2},Val{:cont},Val{L2}})
  _mfdparallelcheck(s₁, s₂)
  N, D, _ = _mfdparallel(s₁, s₂)
  MatrixFractionDescription(N, D, Val{L1})
end

function +{T1,T2,L1,L2}(s₁::MatrixFractionDescription{Val{T1},Val{:disc},Val{L1}},
  s₂::MatrixFractionDescription{Val{T2},Val{:disc},Val{L2}})
  _mfdparallelcheck(s₁, s₂)
  N, D, Ts = _mfdparallel(s₁, s₂)
  MatrixFractionDescription(N, D, Ts, Val{L1})
end

.+(s₁::MatrixFractionDescription{Val{:siso}}, s₂::MatrixFractionDescription{Val{:siso}}) = +(s₁, s₂)

+{T,S}(s::MatrixFractionDescription{Val{T},Val{S}}, g::Union{Real,AbstractMatrix}) =
  +(s, convert(typeof(s), g))
+{T,S}(g::Union{Real,AbstractMatrix}, s::MatrixFractionDescription{Val{T},Val{S}}) =
  +(convert(typeof(s), g), s)

.+(s::MatrixFractionDescription{Val{:siso}}, g::Real) = +(s, g)
.+(g::Real, s::MatrixFractionDescription{Val{:siso}}) = +(g, s)

# Subtraction
-(s₁::MatrixFractionDescription, s₂::MatrixFractionDescription) = +(s₁, -s₂)

.-(s₁::MatrixFractionDescription{Val{:siso}}, s₂::MatrixFractionDescription{Val{:siso}}) = -(s₁, s₂)

-{T,S}(s::MatrixFractionDescription{Val{T},Val{S}}, g::Union{Real,AbstractMatrix}) =
  -(s, convert(typeof(s), g))
-{T,S}(g::Union{Real,AbstractMatrix}, s::MatrixFractionDescription{Val{T},Val{S}}) =
  -(convert(typeof(s), g), s)

.-(s::MatrixFractionDescription{Val{:siso}}, g::Real)    = -(s, g)
.-(g::Real, s::MatrixFractionDescription{Val{:siso}})    = -(g, s)

# Multiplication
function _mfdseriescheck{T1,T2,S}(s₁::MatrixFractionDescription{Val{T1},Val{S}},
  s₂::MatrixFractionDescription{Val{T2},Val{S}})
  # Remark: s₁*s₂ implies u -> s₂ -> s₁ -> y

  if s₁.Ts ≉ s₂.Ts && s₁.Ts ≠ zero(s₁.Ts) && s₂.Ts ≠ zero(s₂.Ts)
    warn("series(s₁,s₂): Sampling time mismatch")
    throw(DomainError())
  end

  if size(s₁,2) ≠ size(s₂,1)
    warn("series(s₁,s₂): size(s₁,2) ≠ size(s₂,1)")
    throw(DomainError())
  end
end

# siso version
function _mfdseries{S,L1,L2}(s₁::MatrixFractionDescription{Val{:siso},Val{S},Val{L1}},
  s₂::MatrixFractionDescription{Val{:siso},Val{S},Val{L2}})
  R₁  = gcd(s₁.D, s₂.N)
  R₂  = gcd(s₂.D, s₁.N)
  D   = div(s₁.D, R₁)*div(s₂.D, R₂)
  N   = div(s₂.N, R₁)*div(s₁.N, R₂)
  N, D, max(s₁.Ts, s₂.Ts)
end

# mimo lfd version
function _mfdseries{S}(s₁::MatrixFractionDescription{Val{:mimo},Val{S},Val{:lfd}},
  s₂::MatrixFractionDescription{Val{:mimo},Val{S},Val{lfd}})
  # D₁^-1 N₁ D₂^-1 N₂
  sᵢ  = lfd(rfd(s₁.N, s₂.D))
  # D₁^-1 Dᵢ^-1 Nᵢ N₂
  D   = Dᵢ*s₁.D
  N   = Nᵢ*N₂
  # ensure coprimeness
  L, V₁, V₂ = gcld(N,D)
  V₁, V₂, max(s₁.Ts, s₂.Ts)
end

# mimo rfd version
function _mfdseries{S}(s₁::MatrixFractionDescription{Val{:mimo},Val{S},Val{:rfd}},
  s₂::MatrixFractionDescription{Val{:mimo},Val{S},Val{:rfd}})
  # N₁ D₁^-1 N₂ D₂^-1
  sᵢ  = rfd(lfd(s₂.N, s₁.D))
  # N₁ Nᵢ Dᵢ^-1 D₂^-1
  D   = s₂.D*Dᵢ
  N   = N₁*Nᵢ
  # ensure coprimeness
  R, V₁, V₂ = gcrd(N,D)
  V₁, V₂, max(s₁.Ts, s₂.Ts)
end

# mimo mixed lfd/rfd versions
function _mfdseries{S}(s₁::MatrixFractionDescription{Val{:mimo},Val{S},Val{:lfd}},
  s₂::MatrixFractionDescription{Val{:mimo},Val{S},Val{:rfd}})
  # N₁ D₁^-1 D₂^-1 N₂
  Dᵢ  = s₂.D*s₁.D
  sᵢ  = lfd(rfd(s₂.N, Dᵢ))
  #  Dᵢ^-1 Nᵢ N₂
  N   = sᵢ.N*N₂
  # ensure coprimeness
  L, V₁, V₂ = gcld(N, sᵢ.D)
  V₁, V₂, max(s₁.Ts, s₂.Ts)
end

function _mfdseries{S}(s₁::MatrixFractionDescription{Val{:mimo},Val{S},Val{:rfd}},
  s₂::MatrixFractionDescription{Val{:mimo},Val{S},Val{:lfd}})
  # D₁^-1 N₁ N₂ D₂^-1
  Nᵢ  = s₁.N*s₂.N
  sᵢ  = rfd(lfd(Nᵢ, s₁.D))
  #  Nᵢ Dᵢ^-1
  D   = D₂*sᵢ.D
  # ensure coprimeness
  R, V₁, V₂ = gcrd(sᵢ.N, D)
  V₁, V₂, max(s₁.Ts, s₂.Ts)
end

# siso and mimo of dimensions 1×1
function _mfdseries{T1,T2,S,L1,L2}(s₁::MatrixFractionDescription{Val{T1},Val{S},Val{L1}},
  s₂::MatrixFractionDescription{Val{T2},Val{S},Val{L2}})
  _mfdseries(mimo(s₁), mimo(s₂))
end

function *{T1,T2,L1,L2}(s₁::MatrixFractionDescription{Val{T1},Val{:cont},Val{L1}},
  s₂::MatrixFractionDescription{Val{T2},Val{:cont},Val{L2}})
  _mfdseriescheck(s₁, s₂)
  N, D, _ = _mfdseries(s₁, s₂)
  MatrixFractionDescription(N, D, Val{L1})
end

function *{T1,T2,L1,L2}(s₁::MatrixFractionDescription{Val{T1},Val{:disc},Val{L1}},
  s₂::MatrixFractionDescription{Val{T2},Val{:disc},Val{L2}})
  _mfdseriescheck(s₁, s₂)
  N, D, Ts = _mfdseries(s₁, s₂)
  MatrixFractionDescription(N, D, Ts, Val{L1})
end

.*(s₁::MatrixFractionDescription{Val{:siso}}, s₂::MatrixFractionDescription{Val{:siso}}) = *(s₁, s₂)

*{T,S}(s::MatrixFractionDescription{Val{T},Val{S}}, g::Union{Real,AbstractMatrix}) =
  *(s, convert(typeof(s), g))
*{T,S}(g::Union{Real,AbstractMatrix}, s::MatrixFractionDescription{Val{T},Val{S}}) =
  *(convert(typeof(s), g), s)

.*(s::MatrixFractionDescription{Val{:siso}}, g::Real)    = *(s, g)
.*(g::Real, s::MatrixFractionDescription{Val{:siso}})    = *(g, s)

## Comparison
=={T,S,L}(s₁::MatrixFractionDescription{T,S,L}, s₂::MatrixFractionDescription{T,S,L}) =
  (s₁.N == s₂.N) && (s₁.D == s₂.D) && (s₁.Ts == s₂.Ts)
=={T1,S1,L1,T2,S2,L2}(s₁::MatrixFractionDescription{T1,S1,L1}, s₂::MatrixFractionDescription{T2,S2,L2}) = false

hash(s::MatrixFractionDescription, h::UInt)     = hash(s.D, hash(S.N, hash(S.Ts, h)))
isequal(s₁::MatrixFractionDescription, s₂::MatrixFractionDescription) = (hash(s₁) == hash(s₂))

function isapprox{T,S,L1,L2,M1,M2,M3,M4}(s₁::MatrixFractionDescription{T,S,L1,M1,M2}, s₂::MatrixFractionDescription{T,S,L2,M3,M4};
  rtol::Real=Base.rtoldefault(promote_type(eltype(mattype(s₁.N)),eltype(mattype(s₁.D)),eltype(mattype(s₂.N)),eltype(mattype(s₂.D)))),
  atol::Real=0, norm::Function=vecnorm)
  isapprox(s₁.Ts, s₂.Ts) || return false # quick exit
  lfd1 = lfd(s₁)
  lfd2 = lfd(s₂)

  D1, U = hermite(lfd1.D)
  N1    = lfd1.N*U
  D2, U = hermite(lfd2.D)
  N2    = lfd2.N*U
  return isapprox(D1, D2; rtol=rtol, atol=atol, norm=norm) && isapprox(N1,N2; rtol=rtol, atol=atol, norm=norm)
end
