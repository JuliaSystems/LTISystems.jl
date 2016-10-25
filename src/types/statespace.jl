immutable StateSpace{T,S,M1,M2,M3,M4} <: LtiSystem{T,S}
  A::M1
  B::M2
  C::M3
  D::M4
  nx::Int
  nu::Int
  ny::Int
  Ts::Float64

  # Continuous-time, single-input-single-output state-space model
  @compat function (::Type{StateSpace}){M1<:AbstractMatrix,M2<:AbstractMatrix,
    M3<:AbstractMatrix, M4<:Real}(A::M1, B::M2, C::M3, D::M4)
    d = fill(D,1,1)
    nx, nu, ny = sscheck(A, B, C, d)
    new{Siso{true},Continuous{true},M1,M2,M3,Matrix{M4}}(A, B, C, d, nx, nu, ny,
      zero(Float64))
  end

  # Discrete-time, single-input-single-output state-space model
  @compat function (::Type{StateSpace}){M1<:AbstractMatrix,M2<:AbstractMatrix,
    M3<:AbstractMatrix, M4<:Real}(A::M1, B::M2, C::M3, D::M4, Ts::Real)
    d = fill(D,1,1)
    nx, nu, ny = sscheck(A, B, C, d, Ts)
    new{Siso{true},Continuous{false},M1,M2,M3,Matrix{M4}}(A, B, C, d, nx, nu, ny,
      convert(Float64, Ts))
  end

  # Continuous-time, multi-input-multi-output state-space model
  @compat function (::Type{StateSpace}){M1<:AbstractMatrix,M2<:AbstractMatrix,
    M3<:AbstractMatrix, M4<:AbstractMatrix}(A::M1, B::M2, C::M3, D::M4)
    nx, nu, ny = sscheck(A, B, C, D)
    new{Siso{false},Continuous{true},M1,M2,M3,M4}(A, B, C, D, nx, nu, ny,
      zero(Float64))
  end

  # Discrete-time, multi-input-multi-output state-space model
  @compat function (::Type{StateSpace}){M1<:AbstractMatrix,M2<:AbstractMatrix,
    M3<:AbstractMatrix, M4<:AbstractMatrix}(A::M1, B::M2, C::M3, D::M4, Ts::Real)
    nx, nu, ny = sscheck(A, B, C, D, Ts)
    new{Siso{false},Continuous{false},M1,M2,M3,M4}(A, B, C, D, nx, nu, ny,
      convert(Float64, Ts))
  end
end

# Enforce state-space type invariance
function sscheck(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix,
  D::AbstractMatrix, Ts::Real = zero(Float64))
  na, ma  = size(A)
  nb, mb  = size(B)
  nc, mc  = size(C)
  nd, md  = size(D)

  @assert na == ma                    "StateSpace: A must be square"
  @assert eltype(A) <: Real           "StateSpace: A must be a matrix of real numbers"
  @assert na == nb                    "StateSpace: A and B must have the same number of rows"
  @assert eltype(B) <: Real           "StateSpace: B must be a matrix of real numbers"
  @assert ma == mc                    "StateSpace: A and C must have the same number of columns"
  @assert eltype(C) <: Real           "StateSpace: C must be a matrix of real numbers"
  @assert nc == nd && nc ≥ 1          "StateSpace: C and D must have the same number (≥1) of rows"
  @assert mb == md && mb ≥ 1          "StateSpace: B and D must have the same number (≥1) of columns"
  @assert eltype(D) <: Real           "StateSpace: D must be a matrix of real numbers"
  @assert Ts ≥ zero(Ts) && !isinf(Ts) "StateSpace: Ts must be non-negative real number"

  nx = na
  nu = mb
  ny = nc

  return nx, nu, ny
end

# Outer constructors
ss(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::Real = zero(Int8)) =
  StateSpace(A, B, C, D)

ss(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::Real, Ts::Real)    =
  StateSpace(A, B, C, D, Ts)

ss(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::AbstractMatrix)    =
  StateSpace(A, B, C, D)

function ss(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix,
  D::AbstractVector)
  @assert isempty(D) "ss(A,B,C,D): D can only be an empty vector"
  d = spzeros(Int8, size(C,1), size(B,2))
  StateSpace(A, B, C, d)
end

ss(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::AbstractMatrix,
  Ts::Real) = StateSpace(A, B, C, D, Ts)

function ss(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix,
  D::AbstractVector, Ts::Real)
  @assert isempty(D) "ss(A,B,C,D): D can only be an empty vector"
  d = spzeros(Int8, size(C,1), size(B,2))
  StateSpace(A, B, C, d, Ts)
end

ss(D::Real)           = StateSpace(zeros(0,0), zeros(0,1), zeros(1,0), D)
ss(D::Real, Ts::Real) = StateSpace(zeros(0,0), zeros(0,1), zeros(1,0), D, Ts)

ss(D::AbstractMatrix)           = StateSpace(zeros(0,0), zeros(0, size(D,2)),
  zeros(size(D,1), 0), D)
ss(D::AbstractMatrix, Ts::Real) = StateSpace(zeros(0,0), zeros(0, size(D,2)),
  zeros(size(D,1), 0), D, Ts)

# Interfaces
isdiscrete{T,S}(s::StateSpace{T,Continuous{S}}) = !S

samplingtime(s::StateSpace) = s.Ts

numstates(s::StateSpace)    = s.nx
numinputs(s::StateSpace)    = s.nu
numoutputs(s::StateSpace)   = s.ny

# Dimension information
ndims(s::StateSpace{Siso{true}})  = 1
ndims(s::StateSpace{Siso{false}}) = 2
size(s::StateSpace)               = size(s.D)
size(s::StateSpace, dim::Int)     = size(s.D, dim)
size(s::StateSpace, dims::Int...) = size(s.D, dims)

# Iteration interface (meaningful only for MIMO systems)
done(s::StateSpace{Siso{false}}, state::Int)  = done(s.D, state)
next(s::StateSpace{Siso{false}}, state::Int)  = (s[state], state+1)
eltype{S,M1,M2,M3,M4}(::Type{StateSpace{Siso{false},S,M1,M2,M3,M4}}) =
  StateSpace{Siso{true},S,M1}
length(s::StateSpace{Siso{false}})      = length(s.D)
eachindex(s::StateSpace{Siso{false}})   = eachindex(s.D)
start(s::StateSpace{Siso{false}})       = start(s.D)
endof(s::StateSpace{Siso{false}})       = endof(s.D)

# Slicing (`getindex`) of MIMO systems
function getindex(s::StateSpace{Siso{false},Continuous{true}}, row::Int, col::Int)
  @assert 1 ≤ row ≤ s.ny "s[idx,]: idx out of bounds"
  @assert 1 ≤ col ≤ s.nu "s[,idx]: idx out of bounds"
  ss(s.A, view(s.B, :, col:col), view(s.C, row:row, :), s.D[row, col])
end

function getindex(s::StateSpace{Siso{false},Continuous{false}}, row::Int, col::Int)
  @assert 1 ≤ row ≤ s.ny "s[idx,]: idx out of bounds"
  @assert 1 ≤ col ≤ s.nu "s[,idx]: idx out of bounds"
  ss(s.A, view(s.B, :, col:col), view(s.C, row:row, :), s.D[row, col], s.Ts)
end

function getindex(s::StateSpace{Siso{false}}, idx::Int)
  @assert 1 ≤ idx ≤ length(s.D) "s[idx]: idx out of bounds"
  col, row  = divrem(idx-1, s.ny)
  col       += 1
  row       += 1
  s[row, col]
end

getindex(s::StateSpace{Siso{false}}, ::Colon)           = s
getindex(s::StateSpace{Siso{false}}, ::Colon, ::Colon)  = s

function getindex(s::StateSpace{Siso{false},Continuous{true}}, ::Colon, cols)
  @assert 1 ≤ min(cols...) ≤ max(cols...) ≤ s.nu "s[:,cols]: cols out of bounds"
  ss(s.A, view(s.B, :, cols), s.C, view(s.D, :, cols))
end

function getindex(s::StateSpace{Siso{false},Continuous{false}}, ::Colon, cols)
  @assert 1 ≤ min(cols...) ≤ max(cols...) ≤ s.nu "s[:,cols]: cols out of bounds"
  ss(s.A, view(s.B, :, cols), s.C, view(s.D, :, cols), s.Ts)
end

function getindex(s::StateSpace{Siso{false},Continuous{true}}, rows, ::Colon)
  @assert 1 ≤ min(rows...) ≤ max(rows...) ≤ s.ny "s[rows,:]: rows out of bounds"
  ss(s.A, s.B, view(s.C, rows, :), view(s.D, rows, :))
end

function getindex(s::StateSpace{Siso{false},Continuous{false}}, rows, ::Colon)
  @assert 1 ≤ min(rows...) ≤ max(rows...) ≤ s.ny "s[rows,:]: rows out of bounds"
  ss(s.A, s.B, view(s.C, rows, :), view(s.D, rows, :), s.Ts)
end

function getindex(s::StateSpace{Siso{false},Continuous{true}}, rows, cols)
  @assert 1 ≤ min(rows...) ≤ max(rows...) ≤ s.ny "s[rows,cols]: rows out of bounds"
  @assert 1 ≤ min(cols...) ≤ max(cols...) ≤ s.nu "s[rows,cols]: cols out of bounds"
  ss(s.A, view(s.B, :, cols), view(s.C, rows, :), view(s.D, rows, cols))
end

function getindex(s::StateSpace{Siso{false},Continuous{false}}, rows, cols)
  @assert 1 ≤ min(rows...) ≤ max(rows...) ≤ s.ny "s[rows,cols]: rows out of bounds"
  @assert 1 ≤ min(cols...) ≤ max(cols...) ≤ s.nu "s[rows,cols]: cols out of bounds"
  ss(s.A, view(s.B, :, cols), view(s.C, rows, :), view(s.D, rows, cols), s.Ts)
end

# Printing functions
summary(s::StateSpace{Siso{true},Continuous{true}})   =
  string("ss(nx=", s.nx, ")")
summary(s::StateSpace{Siso{true},Continuous{false}})  =
  string("ss(nx=", s.nx, ",Ts=", s.Ts, ")")
summary(s::StateSpace{Siso{false},Continuous{true}})  =
  string("ss(nx=", s.nx, ",nu=", s.nu, ",ny=", s.ny, ")")
summary(s::StateSpace{Siso{false},Continuous{false}}) =
  string("ss(nx=", s.nx, ",nu=", s.nu, ",ny=", s.ny, ",Ts=", s.Ts, ")")

showcompact(io::IO, s::StateSpace) = print(io, summary(s))

function show{T}(io::IO, s::StateSpace{T,Continuous{true}})
  println(io, "Continuous time state space model")
  println(io, "\tẋ = Ax + Bu")
  println(io, "\ty = Cx + Du")
  print(io, "with nx=", s.nx, ", nu=", s.nu, ", ny=", s.ny, ".")
end

function show{T}(io::IO, s::StateSpace{T,Continuous{false}})
  println(io, "Discrete time state space model")
  println(io, "\tx[k+1] = Ax[k] + Bu[k]")
  println(io, "\ty[k]   = Cx[k] + Du[k]")
  print(io, "with nx=", s.nx, ", nu=", s.nu, ", ny=", s.ny, ", Ts=", s.Ts, ".")
end

function showall(io::IO, s::StateSpace)
  show(io, s)
  println(io, "System matrix (A):")
  println(io, s.A)
  println(io, "Input matrix (B):")
  println(io, s.B)
  println(io, "Output matrix (C):")
  println(io, s.C)
  println(io, "Feedforward matrix (D):")
  print(io, s.D)
end

# Conversion and promotion
promote_rule{T<:Real,S}(::Type{T}, ::Type{StateSpace{Siso{true},S}}) =
  StateSpace{Siso{true},S}
promote_rule{T<:AbstractMatrix,S}(::Type{T}, ::Type{StateSpace{Siso{false},S}}) =
  StateSpace{Siso{false},S}

convert(::Type{StateSpace{Siso{true},Continuous{true}}}, g::Real)             =
  ss(g)
convert(::Type{StateSpace{Siso{true},Continuous{false}}}, g::Real)            =
  ss(g, zero(Float64))
convert(::Type{StateSpace{Siso{false},Continuous{true}}}, g::AbstractMatrix)  =
  ss(g)
convert(::Type{StateSpace{Siso{false},Continuous{false}}}, g::AbstractMatrix) =
  ss(g, zero(Float64))

# Multiplicative and additive identities (meaningful only for SISO)
one(::Type{StateSpace{Siso{true},Continuous{true}}})    =
  ss(one(Int8))
one(::Type{StateSpace{Siso{true},Continuous{false}}})   =
  ss(one(Int8), zero(Float64))
zero(::Type{StateSpace{Siso{true},Continuous{true}}})   =
  ss(zero(Int8))
zero(::Type{StateSpace{Siso{true},Continuous{false}}})  =
  ss(zero(Int8), zero(Float64))

one(s::StateSpace{Siso{true},Continuous{true}})   = StateSpace(similar(s.A,0,0),
  similar(s.B,0,1), similar(s.C,1,0), one(eltype(s.D)))
one(s::StateSpace{Siso{true},Continuous{false}})  = StateSpace(similar(s.A,0,0),
  similar(s.B,0,1), similar(s.C,1,0), one(eltype(s.D)), s.Ts)
zero(s::StateSpace{Siso{true},Continuous{true}})  = StateSpace(similar(s.A,0,0),
  similar(s.B,0,1), similar(s.C,1,0), zero(eltype(s.D)))
zero(s::StateSpace{Siso{true},Continuous{false}}) = StateSpace(similar(s.A,0,0),
  similar(s.B,0,1), similar(s.C,1,0), zero(eltype(s.D)), s.Ts)

# Inverse of a state-space model
function ssinv(s::StateSpace)
  @assert s.ny == s.nu "inv(sys): s.ny ≠ s.nu"
  try
    Dinv = inv(s.D);
    Ainv = s.A - s.B*Dinv*s.C;
    Binv = s.B*Dinv
    Cinv = -Dinv*s.C
    return Ainv, Binv, Cinv, Dinv
  catch
    warn("inv(sys): D is not invertible")
    throw(DomainError())
  end
end

function inv(s::StateSpace{Siso{true},Continuous{true}})
  Ainv, Binv, Cinv, Dinv = ssinv(s)
  StateSpace(Ainv, Binv, Cinv, Dinv[1])
end

function inv(s::StateSpace{Siso{true},Continuous{false}})
  Ainv, Binv, Cinv, Dinv = ssinv(s)
  StateSpace(Ainv, Binv, Cinv, Dinv[1], s.Ts)
end

function inv(s::StateSpace{Siso{false},Continuous{true}})
  Ainv, Binv, Cinv, Dinv = ssinv(s)
  StateSpace(Ainv, Binv, Cinv, Dinv)
end

function inv(s::StateSpace{Siso{false},Continuous{false}})
  Ainv, Binv, Cinv, Dinv = ssinv(s)
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
    return zerovalues::Vector{Complex{Float64}}
  else
    Bf    = W[1:nrc,1:nrc]
    zerovalues = eigfact(Af, Bf).values
    return zerovalues::Vector{Complex{Float64}}
  end
end

# Transmission zeros of a state-space model
tzeros(s::StateSpace) = zeros(minreal(s))

# Poles of a state-space model
function poles(s::StateSpace)
  Am, Bm, Cm, Dm, = minreal(s.A, s.B, s.C, s.D)
  return eigfact(Am).values::Vector{Complex{Float64}}
end

# Negative of a state-space model
-(s::StateSpace{Siso{true},Continuous{true}})   =
  StateSpace(s.A, s.B, -s.C, -s.D[1])
-(s::StateSpace{Siso{true},Continuous{false}})  =
  StateSpace(s.A, s.B, -s.C, -s.D[1], s.Ts)
-(s::StateSpace{Siso{false},Continuous{true}})  =
  StateSpace(s.A, s.B, -s.C, -s.D)
-(s::StateSpace{Siso{false},Continuous{false}}) =
  StateSpace(s.A, s.B, -s.C, -s.D, s.Ts)

# Addition
function ssparallel{T1,T2,S}(s1::StateSpace{T1,S}, s2::StateSpace{T2,S})
  @assert s1.Ts ≈ s2.Ts || s1.Ts == zero(Float64) ||
    s2.Ts == zero(Float64) "parallel(s1,s2): Sampling time mismatch"
  @assert size(s1) == size(s2) "parallel(s1,s2): size(s1) ≠ size(s2)"

  T = promote_type(eltype(s1.A), eltype(s2.A))
  a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B, s2.B)
  c = hcat(s1.C, s2.C)
  d = s1.D + s2.D
  return a, b, c, d
end

function +(s1::StateSpace{Siso{true},Continuous{true}},
  s2::StateSpace{Siso{true},Continuous{true}})
  a, b, c, d = ssparallel(s1, s2)
  StateSpace(a, b, c, d[1])
end

function +(s1::StateSpace{Siso{true},Continuous{false}},
  s2::StateSpace{Siso{true},Continuous{false}})
  a, b, c, d = ssparallel(s1, s2)
  StateSpace(a, b, c, d[1], s1.Ts)
end

function +{T1,T2}(s1::StateSpace{T1,Continuous{true}},
  s2::StateSpace{T2,Continuous{true}})
  a, b, c, d = ssparallel(s1, s2)
  StateSpace(a, b, c, d)
end

function +{T1,T2}(s1::StateSpace{T1,Continuous{false}},
  s2::StateSpace{T2,Continuous{false}})
  a, b, c, d = ssparallel(s1, s2)
  StateSpace(a, b, c, d, s1.Ts)
end

.+(s1::StateSpace{Siso{true}}, s2::StateSpace{Siso{true}}) = +(s1, s2)

+{T}(s::StateSpace{T,Continuous{true}}, g)  = +(s, ss(g))
+{T}(s::StateSpace{T,Continuous{false}}, g) = +(s, ss(g, zero(Float64)))
+{T}(g, s::StateSpace{T,Continuous{true}})  = +(ss(g), s)
+{T}(g, s::StateSpace{T,Continuous{false}}) = +(ss(g, zero(Float64)), s)

.+(s::StateSpace{Siso{true}}, g::Real)  = +(s, g)
.+(g::Real, s::StateSpace{Siso{true}})  = +(g, s)

# Subtraction
-(s1::StateSpace, s2::StateSpace) = +(s1, -s2)

.-(s1::StateSpace{Siso{true}}, s2::StateSpace{Siso{true}}) = -(s1, s2)

-{T}(s::StateSpace{T,Continuous{true}}, g)  = -(s, ss(g))
-{T}(s::StateSpace{T,Continuous{false}}, g) = -(s, ss(g, zero(Float64)))
-{T}(g, s::StateSpace{T,Continuous{true}})  = -(ss(g), s)
-{T}(g, s::StateSpace{T,Continuous{false}}) = -(ss(g, zero(Float64)), s)

.-(s::StateSpace{Siso{true}}, g::Real)  = -(s, g)
.-(g::Real, s::StateSpace{Siso{true}})  = -(g, s)

# Multiplication
function ssseries{T1,T2,S}(s1::StateSpace{T1,S}, s2::StateSpace{T2,S})
  # Remark: s1*s2 implies u -> s2 -> s1 -> y
  @assert s1.Ts ≈ s2.Ts || s1.Ts == zero(Float64) ||
    s2.Ts == zero(Float64) "series(s1,s2): Sampling time mismatch"
  @assert s1.nu == s2.ny "series(s1,s2): s1.nu ≠ s2.ny"

  T = promote_type(eltype(s1.A), eltype(s1.B), eltype(s2.A), eltype(s2.C))

  a = vcat(hcat(s1.A, s1.B*s2.C),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B*s2.D, s2.B)
  c = hcat(s1.C, s1.D*s2.C)
  d = s1.D * s2.D
  return a, b, c, d
end

function *(s1::StateSpace{Siso{true},Continuous{true}},
  s2::StateSpace{Siso{true},Continuous{true}})
  a, b, c, d = ssseries(s1, s2)
  StateSpace(a, b, c, d[1])
end

function *(s1::StateSpace{Siso{true},Continuous{false}},
  s2::StateSpace{Siso{true},Continuous{false}})
  a, b, c, d = ssseries(s1, s2)
  StateSpace(a, b, c, d[1], s1.Ts)
end

function *{T1,T2}(s1::StateSpace{T1,Continuous{true}},
  s2::StateSpace{T2,Continuous{true}})
  a, b, c, d = ssseries(s1, s2)
  StateSpace(a, b, c, d)
end

function *{T1,T2}(s1::StateSpace{T1,Continuous{false}},
  s2::StateSpace{T2,Continuous{false}})
  a, b, c, d = ssseries(s1, s2)
  StateSpace(a, b, c, d, s1.Ts)
end

.*(s1::StateSpace{Siso{true}}, s2::StateSpace{Siso{true}}) = *(s1, s2)

*{T}(s::StateSpace{T,Continuous{true}}, g)  = *(s, ss(g))
*{T}(s::StateSpace{T,Continuous{false}}, g) = *(s, ss(g, zero(Float64)))
*{T}(g, s::StateSpace{T,Continuous{true}})  = *(ss(g), s)
*{T}(g, s::StateSpace{T,Continuous{false}}) = *(ss(g, zero(Float64)), s)

.*(s::StateSpace{Siso{true}}, g::Real)  = *(s, g)
.*(g::Real, s::StateSpace{Siso{true}})  = *(g, s)

# Division
/(s1::StateSpace, s2::StateSpace) = *(s1, inv(s2))

./(s1::StateSpace{Siso{true}}, s2::StateSpace{Siso{true}}) = /(s1, s2)

/{T}(s::StateSpace{T,Continuous{true}}, g)  = /(s, ss(g))
/{T}(s::StateSpace{T,Continuous{false}}, g) = /(s, ss(g, zero(Float64)))
/{T}(g, s::StateSpace{T,Continuous{true}})  = /(ss(g), s)
/{T}(g, s::StateSpace{T,Continuous{false}}) = /(ss(g, zero(Float64)), s)

./(s::StateSpace{Siso{true}}, g::Real)  = /(s, g)
./(g::Real, s::StateSpace{Siso{true}})  = /(g, s)
