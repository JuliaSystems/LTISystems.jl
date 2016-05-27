# type definition
immutable CMimoSs{T<:Real,M1<:AbstractMatrix{T},M2<:AbstractMatrix{T},
  M3<:AbstractMatrix{T},M4<:AbstractMatrix{T}} <: MimoSystem
  A::M1
  B::M2
  C::M3
  D::M4
  nx::Int
  nu::Int
  ny::Int

  function call{M1<:AbstractMatrix,M2<:AbstractMatrix,M3<:AbstractMatrix,
    M4<:AbstractMatrix}(::Type{CMimoSs}, A::M1, B::M2, C::M3, D::M4)
    @assert eltype(A) <: Real string("A must be a matrix of T<:Real elements")
    @assert eltype(B) <: Real string("B must be a matrix of T<:Real elements")
    @assert eltype(C) <: Real string("C must be a matrix of T<:Real elements")
    @assert eltype(D) <: Real string("D must be a matrix of T<:Real elements")

    na, ma  = size(A,1,2)
    nb, mb  = size(B,1,2)
    nc, mc  = size(C,1,2)

    d       = isempty(D) ? sparse(Int[],Int[],Int8[],nc,mb) : D
    M5      = typeof(d)
    nd, md  = size(d,1,2)

    T       = promote_type(eltype(A), eltype(B), eltype(C), eltype(d))

    if na != ma
      warn("A must be square")
      throw(DomainError())
    elseif nb != na
      warn("B must have the same row size as that of A")
      throw(DomainError())
    elseif mc != ma
      warn("C must have the same column size as that of A")
      throw(DomainError())
    elseif md != mb
      warn("D must have the same column size as that of B")
      throw(DomainError())
    elseif nd != nc
      warn("D must have the same row size as that of C")
      throw(DomainError())
    end

    new{T,M1,M2,M3,M5}(A, B, C, d, na, mb, nc)
  end
end

# interface implementation
isdiscrete(s::CMimoSs)          = false
isdiscrete(::Type{CMimoSs})     = false
samplingtime(s::CMimoSs)        = NaN64

# I/O mapping
numstates(s::CMimoSs)           = s.nx
numinputs(s::CMimoSs)           = s.nu
numoutputs(s::CMimoSs)          = s.ny

# Dimension information
ndims(s::CMimoSs)               = 2
size(s::CMimoSs)                = size(s.D)
size(s::CMimoSs, dim::Int)      = size(s.D, dim)
size(s::CMimoSs, dims::Int...)  = size(s.D, dims)

# overload iteration interface
done(s::CMimoSs, state::Int)                          = done(s.D, state)
eltype{T,M1,M2,M3,M4}(::Type{CMimoSs{T,M1,M2,M3,M4}}) = CSisoSs{T}
length(s::CMimoSs)                                    = length(s.D)
eachindex(s::CMimoSs)                                 = eachindex(s.D)
endof(s::CMimoSs)                                     = endof(s.D)

# overload slicing functions
function getindex(s::CMimoSs, idx::Int)
  if idx < 1 || idx > length(s.D)
    warn("s[idx]: Trying to access idx < 1 or idx > length(s.D)")
    throw(BoundsError(s.D, idx))
  end

  col, row = divrem(idx-1, s.ny)
  col += 1
  row += 1

  CSisoSs(s.A, sub(s.B, :, col:col), sub(s.C, row:row, :),
    sub(s.D, row:row, col:col))
end

function getindex(s::CMimoSs, row::Int, col::Int)
  if row < 1 || row > s.ny
    warn("s[i,]: Trying to access non-existent outputs")
    throw(BoundsError(s.C, row))
  elseif col < 1 || col > s.nu
    warn("s[,j]: Trying to access non-existent inputs")
    throw(BoundsError(s.B, col))
  end

  CSisoSs(s.A, sub(s.B, :, col:col), sub(s.C, row:row, :),
    sub(s.D, row:row, col:col))
end

getindex(s::CMimoSs, ::Colon, ::Colon)  = s

getindex(s::CMimoSs, rows, cols)        = CMimoSs(s.A, sub(s.B, :, cols),
  sub(s.C, rows, :), sub(s.D, rows, cols))

getindex(s::CMimoSs, ::Colon, cols)     = CMimoSs(s.A, sub(s.B, :, cols), s.C,
  sub(s.D, :, cols))

getindex(s::CMimoSs, rows, ::Colon)     = CMimoSs(s.A, s.B, sub(s.C, rows, :),
  sub(s.D, rows, :))

# overload printing functions
summary(s::CMimoSs) = string("ss(nx=", s.nx, ",nu=", s.nu, ",ny=", s.ny, ")")

showcompact(io::IO, s::CMimoSs) = print(io, summary(s))

function show(io::IO, s::CMimoSs)
  println(io, "Continuous time state space model")
  println(io, "\tẋ = Ax + Bu")
  println(io, "\ty = Cx + Du")
  print(io, "with nx=", s.nx, ", nu=", s.nu, ", ny=", s.ny, ".")
end

function showall(io::IO, s::CMimoSs)
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

# creation of continuous state space types
ss{T1<:AbstractMatrix, T2<:AbstractMatrix, T3<:AbstractMatrix,
  T4<:AbstractMatrix}(A::T1, B::T2, C::T3, D::T4) = CMimoSs(A, B, C, D)

ss{T<:AbstractMatrix}(g::T) = CMimoSs(zeros(Int8,0,0), zeros(Int8,0,size(g,2)),
  zeros(Int8,size(g,1),0), g)

# conversion and promotion
promote_rule{T1,M1,M2,M3,M4,T2<:Real}(::Type{CMimoSs{T1,M1,M2,M3,M4}},
  ::Type{T2}) = CMimoSs
convert{T<:Real}(::Type{CMimoSs}, g::T) = ss(fill(g,1,1))

promote_rule{T1,M1,M2,M3,M4,T2<:AbstractMatrix}(::Type{CMimoSs{T1,M1,M2,M3,M4}},
  ::Type{T2}) = CMimoSs
convert{T<:AbstractMatrix}(::Type{CMimoSs}, g::T) = ss(g)

# overloading identities
one{T}(s::CMimoSs{T})                               = ss(ones(T,1,1))
one{T,M1,M2,M3,M4}(::Type{CMimoSs{T,M1,M2,M3,M4}})  = ss(ones(T,1,1))
zero{T}(s::CMimoSs{T})                              = ss(zeros(T,1,1))
zero{T,M1,M2,M3,M4}(::Type{CMimoSs{T,M1,M2,M3,M4}}) = ss(zeros(T,1,1))

# overload inv and zeros
function inv{T}(s::CMimoSs{T})
  ny, nu = size(s.D, 1, 2)
  @assert ny == nu string("inv(sys): D must be square")

  try
    Dinv = inv(s.D);
    Ainv = s.A - s.B*Dinv*s.C;
    Binv = s.B*Dinv
    Cinv = -Dinv*s.C

    CMimoSs(Ainv, Binv, Cinv, Dinv)
  catch
    warn("inv(sys): D is not invertible")
    throw(DomainError())
  end
end

function zeros{T}(s::CMimoSs{T})
  Ar, Br, Cr, Dr, mr, nr, pr        = reduce(s.A, s.B, s.C, s.D)
  if nr == 0
    return (Complex{Float64}[], mr::Int)
  end
  Arc, Brc, Crc, Drc, mrc, nrc, prc = reduce(Ar.', Cr.', Br.', Dr.')
  if nrc == 0
    return (Complex{Float64}[], mrc::Int)
  end

  svdobj  = svdfact([Crc Drc], thin = false)
  W       = flipdim(svdobj.Vt', 2)
  Af      = [Arc Brc]*W[:, 1:nrc]

  if mrc == 0
    zerovalues = eigfact(Af).values
    # return (zerovalues::Vector{Complex{Float64}}, mrc::Int)
    return zerovalues
  else
    Bf    = W[1:nrc,1:nrc]
    zerovalues = eigfact(Af, Bf).values
    # return (zerovalues::Vector{Complex{Float64}}, mrc::Int)
    return zerovalues
  end
end

function poles{T}(s::CMimoSs{T})
  Am, Bm, Cm, Dm, = minreal(s.A, s.B, s.C, s.D)
  return eigfact(Am).values
end

# overload mathematical operations
-(s::CMimoSs) = CMimoSs(s.A, s.B, -s.C, -s.D)

function +{T1, T2}(s1::CMimoSs{T1}, s2::CMimoSs{T2})
  # Ensure systems have same shapes
  if size(s1) != size(s2)
    warn("s1+s2: size(s1) != size(s2)")
    throw(DomainError())
  end

  T = promote_type(T1, T2)

  a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B, s2.B)
  c = hcat(s1.C, s2.C)
  d = s1.D + s2.D

  CMimoSs(a,b,c,d)
end

+{T<:Real}(s::CMimoSs, g::T)  = CMimoSs(copy(s.A), copy(s.B), copy(s.C), s.D + g)
+{T<:Real}(g::T, s::CMimoSs)  = +(s, g)

function +{T<:AbstractMatrix}(s::CMimoSs, g::T)
  if size(s.D) != size(g)
    warn("s+g: size(s.D) != size(g)")
    throw(DomainError())
  end

  CMimoSs(copy(s.A), copy(s.B), copy(s.C), s.D + g)
end
+{T<:AbstractMatrix}(g::T, s::CMimoSs) = +(s, g)

.+(s1::CMimoSs, s2::CMimoSs)  = +(s1, s2)
.+{T<:Real}(s::CMimoSs, g::T) = +(s, g)
.+{T<:Real}(g::T, s::CMimoSs) = +(s, g)

function -{T1, T2}(s1::CMimoSs{T1}, s2::CMimoSs{T2})
  # Ensure systems have same shapes
  if size(s1) != size(s2)
    warn("s1-s2: size(s1) != size(s2)")
    throw(DomainError())
  end

  T = promote_type(T1, T2)

  a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B, s2.B)
  c = hcat(s1.C, -s2.C)
  d = s1.D - s2.D

  CMimoSs(a,b,c,d)
end

-{T<:Real}(s::CMimoSs, g::T)  = CMimoSs(copy(s.A), copy(s.B), copy(s.C), s.D - g)
-{T<:Real}(g::T, s::CMimoSs)  = +(g, -s)

function -{T<:AbstractMatrix}(s::CMimoSs, g::T)
  if size(s.D) != size(g)
    warn("s-g: size(s.D) != size(g)")
    throw(DomainError())
  end

  CMimoSs(copy(s.A), copy(s.B), copy(s.C), s.D - g)
end
-{T<:AbstractMatrix}(g::T, s::CMimoSs) = +(g, -s)

.-(s1::CMimoSs, s2::CMimoSs)  = -(s1, s2)
.-{T<:Real}(s::CMimoSs, g::T) = -(s, g)
.-{T<:Real}(g::T, s::CMimoSs) = +(g, -s)

function *{T1, T2}(s1::CMimoSs{T1}, s2::CMimoSs{T2})
  # Remark: s1*s2 implies u -> s2 -> s1 -> y
  if s1.nu != s2.ny
    warn("s1*s2: s1.nu != s2.ny")
    throw(DomainError())
  end

  T = promote_type(T1, T2)

  a = vcat(hcat(s1.A, s1.B*s2.C),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B*s2.D, s2.B)
  c = hcat(s1.C, s1.D*s2.C)
  d = s1.D * s2.D

  CMimoSs(a,b,c,d)
end

*{T<:Real}(s::CMimoSs, g::T) = CMimoSs(copy(s.A), s.B*g, copy(s.C), s.D*g)
*{T<:Real}(g::T, s::CMimoSs) = CMimoSs(copy(s.A), copy(s.B), g*s.C, g*s.D)

function *{T<:AbstractMatrix}(s::CMimoSs, g::T)
  if s.nu != size(g, 1)
    warn("s*g: s.nu != size(g, 1)")
    throw(DomainError())
  end

  CMimoSs(copy(s.A), s.B*g, copy(s.C), s.D*g)
end

function *{T<:AbstractMatrix}(g::T, s::CMimoSs)
  if s.ny != size(g, 2)
    warn("g*s: s.ny != size(g, 2)")
    throw(DomainError())
  end

  CMimoSs(copy(s.A), copy(s.B), g*s.C, g*s.D)
end

.*(s1::CMimoSs, s2::CMimoSs)  = *(s1, s2)
.*{T<:Real}(s::CMimoSs, g::T) = *(s, g)
.*{T<:Real}(g::T, s::CMimoSs) = *(g, s)

/(s1::CMimoSs, s2::CMimoSs)   = *(s1, inv(s2))

/{T<:Real}(s::CMimoSs, g::T)  = CMimoSs(copy(s.A), s.B/g, copy(s.C), s.D/g)
/{T<:Real}(g::T, s::CMimoSs)  = *(g, inv(s))

function /{T<:AbstractMatrix}(s::CMimoSs, g::T)
  ginv = inv(g)

  CMimoSs(copy(s.A), s.B*ginv, copy(s.C), s.D*ginv)
end
/{T<:AbstractMatrix}(g::T, s::CMimoSs)  = *(g, inv(s))

./(s1::CMimoSs, s2::CMimoSs)  = /(s1, s2)
./{T<:Real}(s::CMimoSs, g::T) = /(s, g)
./{T<:Real}(g::T, s::CMimoSs) = /(g, s)

function ==(s1::CMimoSs, s2::CMimoSs)
  # TODO: Implement
  throw(ErrorException("==(s1,s2) for CMimoSs is not implemented"))
end

!=(s1::CMimoSs, s2::CMimoSs) = !(s1 == s2)

function isapprox(s1::CMimoSs, s2::CMimoSs)
  # TODO: Implement
  throw(ErrorException("isapprox(s1,s2) CMimoSs is not implemented"))
end
