# type definition
immutable CSisoSs{T<:Real,M1<:AbstractMatrix{T},M2<:AbstractMatrix{T},
  M3<:AbstractMatrix{T},M4<:AbstractMatrix{T}} <: SisoSs{T}
  A::M1
  B::M2
  C::M3
  D::M4
  nx::Int

  function call{M1<:AbstractMatrix,M2<:AbstractMatrix,M3<:AbstractMatrix,
    M4<:AbstractMatrix}(::Type{CSisoSs}, A::M1, B::M2, C::M3, D::M4)
    @assert eltype(A) <: Real string("A must be a matrix of T<:Real elements")
    @assert eltype(B) <: Real string("B must be a matrix of T<:Real elements")
    @assert eltype(C) <: Real string("C must be a matrix of T<:Real elements")
    @assert eltype(D) <: Real string("D must be a matrix of T<:Real elements")

    T       = promote_type(eltype(A), eltype(B), eltype(C), eltype(D))
    na, ma  = size(A,1,2)
    nb, mb  = size(B,1,2)
    nc, mc  = size(C,1,2)
    nd, md  = size(D,1,2)

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

    new{T,M1,M2,M3,M4}(A, B, C, D, na)
  end
end

# interface implementation
isdiscrete(s::CSisoSs)          = false
isdiscrete(::Type{CSisoSs})     = false
samplingtime(s::CSisoSs)        = NaN64

# I/O mapping
numstates(s::CSisoSs)           = s.nx

# overload slicing functions
function getindex(s::CSisoSs, idx::Int)
  if idx != 1
    warn("s[idx]: Trying to access idx != 1")
    throw(BoundsError(s.D, idx))
  end

  return s
end

# overload printing functions
summary(s::CSisoSs)             = string("ss(nx=", s.nx, ")")
showcompact(io::IO, s::CSisoSs) = print(io, summary(s))

function show(io::IO, s::CSisoSs)
  println(io, "Continuous time state space model")
  println(io, "\tẋ = Ax + Bu")
  println(io, "\ty = Cx + Du")
  print(io, "with nx=", s.nx, ".")
end

function showall(io::IO, s::CSisoSs)
  show(io, s)
  println(io)
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
ss{T1<:AbstractMatrix, T2<:AbstractMatrix, T3<:AbstractMatrix, T4<:Real}(A::T1,
  B::T2, C::T3, D::T4 = zero(Int8)) = CSisoSs(A, B, C, fill(D,1,1))

ss{T<:Real}(g::T) = CSisoSs(zeros(Int8,0,0), zeros(Int8,0,1), zeros(Int8,1,0),
  fill(g,1,1))

# conversion and promotion
promote_rule{T1,M1,M2,M3,M4,T2<:Real}(::Type{CSisoSs{T1,M1,M2,M3,M4}},
  ::Type{T2}) = CSisoSs
convert{T<:Real}(::Type{CSisoSs}, g::T) = ss(g)

# overloading identities
one{T}(s::CSisoSs{T})                               = ss(one(T))
one{T,M1,M2,M3,M4}(::Type{CSisoSs{T,M1,M2,M3,M4}})  = ss(one(T))
one(::Type{CSisoSs})                                = ss(one(Int8))
zero{T}(s::CSisoSs{T})                              = ss(zero(T))
zero{T,M1,M2,M3,M4}(::Type{CSisoSs{T,M1,M2,M3,M4}}) = ss(zero(T))
zero(::Type{CSisoSs})                               = ss(zero(Int8))

# overload inv and zeros
function inv{T}(s::CSisoSs{T})
  if s.D[1] == zero(eltype(s.D))
    warn("inv(sys): D is not invertible")
    throw(DomainError())
  end

  Dinv = inv(s.D);
  Ainv = s.A - s.B*Dinv*s.C;
  Binv = s.B*Dinv
  Cinv = -Dinv*s.C

  CSisoSs(Ainv, Binv, Cinv, Dinv)
end

function zeros{T}(s::CSisoSs{T})
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

function poles{T}(s::CSisoSs{T})
  Am, Bm, Cm, Dm, = minreal(s.A, s.B, s.C, s.D)
  return eigfact(Am).values
end

# overload mathematical operations
-(s::CSisoSs) = CSisoSs(s.A, s.B, -s.C, -s.D)

function +{T1, T2}(s1::CSisoSs{T1}, s2::CSisoSs{T2})
  T = promote_type(T1, T2)

  a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B, s2.B)
  c = hcat(s1.C, s2.C)
  d = s1.D + s2.D

  CSisoSs(a,b,c,d)
end

+{T<:Real}(s::CSisoSs, g::T)  = CSisoSs(copy(s.A), copy(s.B), copy(s.C), s.D + g)
+{T<:Real}(g::T, s::CSisoSs)  = +(s, g)

.+(s1::CSisoSs, s2::CSisoSs)  = +(s1, s2)
.+{T<:Real}(s::CSisoSs, g::T) = +(s, g)
.+{T<:Real}(g::T, s::CSisoSs) = +(s, g)

function -{T1, T2}(s1::CSisoSs{T1}, s2::CSisoSs{T2})
  T = promote_type(T1, T2)

  a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B, s2.B)
  c = hcat(s1.C, -s2.C)
  d = s1.D - s2.D

  CSisoSs(a,b,c,d)
end

-{T<:Real}(s::CSisoSs, g::T)  = CSisoSs(copy(s.A), copy(s.B), copy(s.C), s.D - g)
-{T<:Real}(g::T, s::CSisoSs)  = +(g, -s)

.-(s1::CSisoSs, s2::CSisoSs)  = -(s1, s2)
.-{T<:Real}(s::CSisoSs, g::T) = -(s, g)
.-{T<:Real}(g::T, s::CSisoSs) = +(g, -s)

function *{T1, T2}(s1::CSisoSs{T1}, s2::CSisoSs{T2})
  # Remark: "y = (s1*s2) u" implies "u -> s2 -> s1 -> y"

  T = promote_type(T1, T2)

  a = vcat(hcat(s1.A, s1.B*s2.C'),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B*s2.D, s2.B)
  c = hcat(s1.C, s1.D*s2.C)
  d = s1.D * s2.D

  CSisoSs(a,b,c,d)
end

*{T<:Real}(s::CSisoSs, g::T)  = CSisoSs(copy(s.A), s.B*g, copy(s.C), s.D*g)
*{T<:Real}(g::T, s::CSisoSs)  = CSisoSs(copy(s.A), copy(s.B), g*s.C, g*s.D)

.*(s1::CSisoSs, s2::CSisoSs)  = *(s1, s2)
.*{T<:Real}(s::CSisoSs, g::T) = *(s, g)
.*{T<:Real}(g::T, s::CSisoSs) = *(g, s)

/(s1::CSisoSs, s2::CSisoSs)   = *(s1, inv(s2))

/{T<:Real}(s::CSisoSs, g::T)  = CSisoSs(copy(s.A), s.B/g, copy(s.C), s.D/g)
/{T<:Real}(g::T, s::CSisoSs)  = *(g, inv(s))

./(s1::CSisoSs, s2::CSisoSs)  = /(s1, s2)
./{T<:Real}(s::CSisoSs, g::T) = /(s, g)
./{T<:Real}(g::T, s::CSisoSs) = /(g, s)

function ==(s1::CSisoSs, s2::CSisoSs)
  # TODO: Implement
  throw(ErrorException("==(s1,s2) for CSisoSs is not implemented"))
end

!=(s1::CSisoSs, s2::CSisoSs) = !(s1 == s2)

function isapprox(s1::CSisoSs, s2::CSisoSs)
  # TODO: Implement
  throw(ErrorException("isapprox(s1,s2) for CSisoSs is not implemented"))
end
