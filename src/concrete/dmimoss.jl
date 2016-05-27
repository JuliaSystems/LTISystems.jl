# type definition
immutable DMimoSs{T<:Real,M1<:AbstractMatrix{T},M2<:AbstractMatrix{T},
  M3<:AbstractMatrix{T},M4<:AbstractMatrix{T}} <: MimoSystem
  A::M1
  B::M2
  C::M3
  D::M4
  nx::Int
  nu::Int
  ny::Int
  Ts::Float64

  function call{T1<:Real,M1<:AbstractMatrix,M2<:AbstractMatrix,M3<:AbstractMatrix,
    M4<:AbstractMatrix}(::Type{DMimoSs}, A::M1, B::M2, C::M3, D::M4, Ts::T1)
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
    Ts_     = convert(Float64, Ts > zero(Ts) ? Ts : NaN)

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

    new{T,M1,M2,M3,M5}(A, B, C, d, na, mb, nc, Ts_)
  end
end

# interface implementation
isdiscrete(s::DMimoSs)          = true
isdiscrete(::Type{DMimoSs})     = true
samplingtime(s::DMimoSs)        = s.Ts

# I/O mapping
numstates(s::DMimoSs)           = s.nx
numinputs(s::DMimoSs)           = s.nu
numoutputs(s::DMimoSs)          = s.ny

# Dimension information
ndims(s::DMimoSs)               = 2
size(s::DMimoSs)                = size(s.D)
size(s::DMimoSs, dim::Int)      = size(s.D, dim)
size(s::DMimoSs, dims::Int...)  = size(s.D, dims)

# overload iteration interface
done(s::DMimoSs, state::Int)                          = done(s.D, state)
eltype{T,M1,M2,M3,M4}(::Type{DMimoSs{T,M1,M2,M3,M4}}) = DSisoSs{T}
length(s::DMimoSs)                                    = length(s.D)
eachindex(s::DMimoSs)                                 = eachindex(s.D)
endof(s::DMimoSs)                                     = endof(s.D)

# overload slicing functions
function getindex(s::DMimoSs, idx::Int)
  if idx < 1 || idx > length(s.D)
    warn("s[idx]: Trying to access idx < 1 or idx > length(s.D)")
    throw(BoundsError(s.D, idx))
  end

  col, row = divrem(idx-1, s.ny)
  col += 1
  row += 1

  DSisoSs(s.A, sub(s.B, :, col:col), sub(s.C, row:row, :),
    sub(s.D, row:row, col:col), s.Ts)
end

function getindex(s::DMimoSs, row::Int, col::Int)
  if row < 1 || row > s.ny
    warn("s[i,]: Trying to access non-existent outputs")
    throw(BoundsError(s.C, row))
  elseif col < 1 || col > s.nu
    warn("s[,j]: Trying to access non-existent inputs")
    throw(BoundsError(s.B, col))
  end

  DSisoSs(s.A, sub(s.B, :, col:col), sub(s.C, row:row, :),
    sub(s.D, row:row, col:col), s.Ts)
end

getindex(s::DMimoSs, ::Colon, ::Colon)  = s

getindex(s::DMimoSs, rows, cols)        = DMimoSs(s.A, sub(s.B, :, cols),
  sub(s.C, rows, :), sub(s.D, rows, cols), s.Ts)

getindex(s::DMimoSs, ::Colon, cols)     = DMimoSs(s.A, sub(s.B, :, cols), s.C,
  sub(s.D, :, cols), s.Ts)

getindex(s::DMimoSs, rows, ::Colon)     = DMimoSs(s.A, s.B, sub(s.C, rows, :),
  sub(s.D, rows, :), s.Ts)

# overload printing functions
summary(s::DMimoSs) = string("ss(nx=", s.nx, ",nu=", s.nu, ",ny=", s.ny, ",Ts=",
  s.Ts, ")")

showcompact(io::IO, s::DMimoSs) = print(io, summary(s))

function show(io::IO, s::DMimoSs)
  println(io, "Discrete time state space model")
  println(io, "\tx[k+1] = Ax[k] + Bu[k]")
  println(io, "\ty[k]   = Cx[k] + Du[k]")
  print(io, "with nx=", s.nx, ", nu=", s.nu, ", ny=", s.ny, ", Ts=", s.Ts, ".")
end

function showall(io::IO, s::DMimoSs)
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
  T4<:AbstractMatrix, T5<:Real}(A::T1, B::T2, C::T3, D::T4, Ts::T5) =
  DMimoSs(A, B, C, D, Ts)

ss{T1<:AbstractMatrix, T2<:Real}(g::T1, Ts::T2) = DMimoSs(zeros(Int8,0,0),
  zeros(Int8,0,size(g,2)), zeros(Int8,size(g,1),0), g, Ts)

# conversion and promotion
promote_rule{T1,M1,M2,M3,M4,T2<:Real}(::Type{DMimoSs{T1,M1,M2,M3,M4}},
  ::Type{T2}) = DMimoSs
convert{T<:Real}(::Type{DMimoSs}, g::T) = ss(fill(g,1,1), 0)

promote_rule{T1,M1,M2,M3,M4,T2<:AbstractMatrix}(::Type{DMimoSs{T1,M1,M2,M3,M4}},
  ::Type{T2}) = DMimoSs
convert{T<:AbstractMatrix}(::Type{DMimoSs}, g::T) = ss(g, 0)

# overloading identities
one{T}(s::DMimoSs{T})                               = ss(ones(T,1,1), 0)
one{T,M1,M2,M3,M4}(::Type{DMimoSs{T,M1,M2,M3,M4}})  = ss(ones(T,1,1), 0)
zero{T}(s::DMimoSs{T})                              = ss(zeros(T,1,1), 0)
zero{T,M1,M2,M3,M4}(::Type{DMimoSs{T,M1,M2,M3,M4}}) = ss(zeros(T,1,1), 0)

# overload inv and zeros
function inv{T}(s::DMimoSs{T})
  ny, nu = size(s.D, 1, 2)
  @assert ny == nu string("inv(sys): D must be square")

  try
    Dinv = inv(s.D);
    Ainv = s.A - s.B*Dinv*s.C;
    Binv = s.B*Dinv
    Cinv = -Dinv*s.C

    DMimoSs(Ainv, Binv, Cinv, Dinv, s.Ts)
  catch
    warn("inv(sys): D is not invertible")
    throw(DomainError())
  end
end

function zeros{T}(s::DMimoSs{T})
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

function poles{T}(s::DMimoSs{T})
  Am, Bm, Cm, Dm, = minreal(s.A, s.B, s.C, s.D)
  return eigfact(Am).values
end

# overload mathematical operations
-(s::DMimoSs) = DMimoSs(s.A, s.B, -s.C, -s.D)

function +{T1, T2}(s1::DMimoSs{T1}, s2::DMimoSs{T2})
  # Ensure systems have same sampling times and shapes
  if !isnan(s1.Ts) && !isnan(s2.Ts) && s1.Ts ≉ s2.Ts
    warn("s1+s2: Sampling time mismatch")
    throw(DomainError())
  elseif size(s1) != size(s2)
    warn("s1+s2: size(s1) != size(s2)")
    throw(DomainError())
  end

  T = promote_type(T1, T2)

  a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B, s2.B)
  c = hcat(s1.C, s2.C)
  d = s1.D + s2.D

  DMimoSs(a,b,c,d,s1.Ts)
end

+{T<:Real}(s::DMimoSs, g::T)  = DMimoSs(copy(s.A), copy(s.B), copy(s.C),
  s.D + g, s.Ts)
+{T<:Real}(g::T, s::DMimoSs)  = +(s, g)

function +{T<:AbstractMatrix}(s::DMimoSs, g::T)
  # Ensure systems have same shapes
  if size(s.D) != size(g)
    warn("s+g: size(s.D) != size(g)")
    throw(DomainError())
  end

  DMimoSs(copy(s.A), copy(s.B), copy(s.C), s.D + g, s.Ts)
end
+{T<:AbstractMatrix}(g::T, s::DMimoSs) = +(s, g)

.+(s1::DMimoSs, s2::DMimoSs)  = +(s1, s2)
.+{T<:Real}(s::DMimoSs, g::T) = +(s, g)
.+{T<:Real}(g::T, s::DMimoSs) = +(s, g)

function -{T1, T2}(s1::DMimoSs{T1}, s2::DMimoSs{T2})
  # Ensure systems have same sampling times and shapes
  if !isnan(s1.Ts) && !isnan(s2.Ts) && s1.Ts ≉ s2.Ts
    warn("s1-s2: Sampling time mismatch")
    throw(DomainError())
  elseif size(s1) != size(s2)
    warn("s1-s2: size(s1) != size(s2)")
    throw(DomainError())
  end

  T = promote_type(T1, T2)

  a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B, s2.B)
  c = hcat(s1.C, -s2.C)
  d = s1.D - s2.D

  DMimoSs(a,b,c,d,s1.Ts)
end

-{T<:Real}(s::DMimoSs, g::T)  = DMimoSs(copy(s.A), copy(s.B), copy(s.C),
  s.D - g, s.Ts)
-{T<:Real}(g::T, s::DMimoSs)  = +(g, -s)

function -{T<:AbstractMatrix}(s::DMimoSs, g::T)
  if size(s.D) != size(g)
    warn("s-g: size(s.D) != size(g)")
    throw(DomainError())
  end

  DMimoSs(copy(s.A), copy(s.B), copy(s.C), s.D - g, s.Ts)
end
-{T<:AbstractMatrix}(g::T, s::DMimoSs) = +(g, -s)

.-(s1::DMimoSs, s2::DMimoSs)  = -(s1, s2)
.-{T<:Real}(s::DMimoSs, g::T) = -(s, g)
.-{T<:Real}(g::T, s::DMimoSs) = +(g, -s)

function *{T1, T2}(s1::DMimoSs{T1}, s2::DMimoSs{T2})
  # Remark: s1*s2 implies u -> s2 -> s1 -> y
  if !isnan(s1.Ts) && !isnan(s2.Ts) && s1.Ts ≉ s2.Ts
    warn("s1*s2: Sampling time mismatch")
    throw(DomainError())
  elseif s1.nu != s2.ny
    warn("s1*s2: s1.nu != s2.ny")
    throw(DomainError())
  end

  T = promote_type(T1, T2)

  a = vcat(hcat(s1.A, s1.B*s2.C),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B*s2.D, s2.B)
  c = hcat(s1.C, s1.D*s2.C)
  d = s1.D * s2.D

  DMimoSs(a,b,c,d,s1.Ts)
end

*{T<:Real}(s::DMimoSs, g::T) = DMimoSs(copy(s.A), s.B*g, copy(s.C), s.D*g, s.Ts)
*{T<:Real}(g::T, s::DMimoSs) = DMimoSs(copy(s.A), copy(s.B), g*s.C, g*s.D, s.Ts)

function *{T<:AbstractMatrix}(s::DMimoSs, g::T)
  if s.nu != size(g, 1)
    warn("s*g: s.nu != size(g, 1)")
    throw(DomainError())
  end

  DMimoSs(copy(s.A), s.B*g, copy(s.C), s.D*g, s.Ts)
end

function *{T<:AbstractMatrix}(g::T, s::DMimoSs)
  if s.ny != size(g, 2)
    warn("g*s: s.ny != size(g, 2)")
    throw(DomainError())
  end

  DMimoSs(copy(s.A), copy(s.B), g*s.C, g*s.D, s.Ts)
end

.*(s1::DMimoSs, s2::DMimoSs)  = *(s1, s2)
.*{T<:Real}(s::DMimoSs, g::T) = *(s, g)
.*{T<:Real}(g::T, s::DMimoSs) = *(g, s)

/(s1::DMimoSs, s2::DMimoSs)   = *(s1, inv(s2))

/{T<:Real}(s::DMimoSs, g::T)  = DMimoSs(copy(s.A), s.B/g, copy(s.C), s.D/g, s.Ts)
/{T<:Real}(g::T, s::DMimoSs)  = *(g, inv(s))

function /{T<:AbstractMatrix}(s::DMimoSs, g::T)
  ginv = inv(g)

  DMimoSs(copy(s.A), s.B*ginv, copy(s.C), s.D*ginv)
end
/{T<:AbstractMatrix}(g::T, s::DMimoSs)  = *(g, inv(s))

./(s1::DMimoSs, s2::DMimoSs)  = /(s1, s2)
./{T<:Real}(s::DMimoSs, g::T) = /(s, g)
./{T<:Real}(g::T, s::DMimoSs) = /(g, s)

function ==(s1::DMimoSs, s2::DMimoSs)
  # TODO: Implement
  !isnan(s1.Ts) && !isnan(s2.Ts) && !isapprox(s1.Ts, s2.Ts) && return false
  throw(ErrorException("==(s1,s2) for DMimoSs is not implemented"))
end

!=(s1::DMimoSs, s2::DMimoSs) = !(s1 == s2)

function isapprox(s1::DMimoSs, s2::DMimoSs)
  # TODO: Implement
  !isnan(s1.Ts) && !isnan(s2.Ts) && !isapprox(s1.Ts, s2.Ts) && return false
  throw(ErrorException("isapprox(s1,s2) DMimoSs is not implemented"))
end
