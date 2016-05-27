# type definition
immutable DSisoSs{T1<:Real,M1<:AbstractMatrix{T1},M2<:AbstractMatrix{T1},
  M3<:AbstractMatrix{T1},M4<:AbstractMatrix{T1}} <: SisoSs{T1}
  A::M1
  B::M2
  C::M3
  D::M4
  nx::Int
  Ts::Float64

  function call{T1<:Real,M1<:AbstractMatrix,M2<:AbstractMatrix,M3<:AbstractMatrix,
    M4<:AbstractMatrix}(::Type{DSisoSs}, A::M1, B::M2, C::M3, D::M4, Ts::T1)
    @assert eltype(A) <: Real string("A must be a matrix of T<:Real elements")
    @assert eltype(B) <: Real string("B must be a matrix of T<:Real elements")
    @assert eltype(C) <: Real string("C must be a matrix of T<:Real elements")
    @assert eltype(D) <: Real string("D must be a matrix of T<:Real elements")

    T       = promote_type(eltype(A), eltype(B), eltype(C), eltype(D))
    na, ma  = size(A,1,2)
    nb, mb  = size(B,1,2)
    nc, mc  = size(C,1,2)
    nd, md  = size(D,1,2)
    Ts_     = convert(Float64, Ts > zero(Ts) ? Ts : NaN)

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

    new{T,M1,M2,M3,M4}(A, B, C, D, na, Ts_)
  end
end

# interface implementation
isdiscrete(s::DSisoSs)      = true
isdiscrete(::Type{DSisoSs}) = true
samplingtime(s::DSisoSs)    = s.Ts

# I/O mapping
numstates(s::DSisoSs)       = s.nx

# overload slicing functions
function getindex(s::DSisoSs, idx::Int)
  if idx != 1
    warn("s[idx]: Trying to access idx != 1")
    throw(BoundsError(s.D, idx))
  end

  return s
end

# overload printing functions
summary(s::DSisoSs)             = string("ss(nx=", s.nx, ",Ts=", s.Ts, ")")
showcompact(io::IO, s::DSisoSs) = print(io, summary(s))

function show(io::IO, s::DSisoSs)
  println(io, "Discrete time state space model")
  println(io, "\tx[k+1] = Ax[k] + Bu[k]")
  println(io, "\ty[k]   = Cx[k] + Du[k]")
  print(io, "with nx=", s.nx, ", Ts=", s.Ts, ".")
end

function showall(io::IO, s::DSisoSs)
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
ss{T1<:AbstractMatrix, T2<:AbstractMatrix, T3<:AbstractMatrix, T4<:Real,
  T5<:Real}(A::T1, B::T2, C::T3, D::T4, Ts::T5) = DSisoSs(A, B, C, fill(D,1,1), Ts)

ss{T1<:Real, T2<:Real}(g::T1, Ts::T2) = DSisoSs(zeros(Int8,0,0), zeros(Int8,0,1),
  zeros(Int8,1,0), fill(g,1,1), Ts)

# conversion and promotion
promote_rule{T1,M1,M2,M3,M4,T2<:Real}(::Type{DSisoSs{T1,M1,M2,M3,M4}},
  ::Type{T2}) = DSisoSs
convert{T<:Real}(::Type{DSisoSs}, g::T) = ss(g, 0)

# overloading identities
one{T}(s::DSisoSs{T})                               = ss(one(T), 0)
one{T,M1,M2,M3,M4}(::Type{DSisoSs{T,M1,M2,M3,M4}})  = ss(one(T), 0)
one(::Type{DSisoSs})                                = ss(one(Int8), 0)
zero{T}(s::DSisoSs{T})                              = ss(zero(T), 0)
zero{T,M1,M2,M3,M4}(::Type{DSisoSs{T,M1,M2,M3,M4}}) = ss(zero(T), 0)
zero(::Type{DSisoSs})                               = ss(zero(Int8), 0)

# overload inv and zeros
function inv{T}(s::DSisoSs{T})
  if s.D[1] == zero(eltype(s.D))
    warn("inv(sys): D is not invertible")
    throw(DomainError())
  end

  Dinv = inv(s.D);
  Ainv = s.A - s.B*Dinv*s.C;
  Binv = s.B*Dinv
  Cinv = -Dinv*s.C

  DSisoSs(Ainv, Binv, Cinv, Dinv, s.Ts)
end

function zeros{T}(s::DSisoSs{T})
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

function poles{T}(s::DSisoSs{T})
  Am, Bm, Cm, Dm, = minreal(s.A, s.B, s.C, s.D)
  return eigfact(Am).values
end

# overload mathematical operations
-(s::DSisoSs) = DSisoSs(s.A, s.B, -s.C, -s.D)

function +{T1, T2}(s1::DSisoSs{T1}, s2::DSisoSs{T2})
  if !isnan(s1.Ts) && !isnan(s2.Ts) && s1.Ts ≉ s2.Ts
    warn("s1+s2: Sampling time mismatch")
    throw(DomainError())
  end
  T = promote_type(T1, T2)

  a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B, s2.B)
  c = hcat(s1.C, s2.C)
  d = s1.D + s2.D

  DSisoSs(a,b,c,d,s1.Ts)
end

+{T<:Real}(s::DSisoSs, g::T)  = DSisoSs(copy(s.A), copy(s.B), copy(s.C),
  s.D + g, s.Ts)
+{T<:Real}(g::T, s::DSisoSs)  = +(s, g)

.+(s1::DSisoSs, s2::DSisoSs)  = +(s1, s2)
.+{T<:Real}(s::DSisoSs, g::T) = +(s, g)
.+{T<:Real}(g::T, s::DSisoSs) = +(s, g)

function -{T1, T2}(s1::DSisoSs{T1}, s2::DSisoSs{T2})
  if !isnan(s1.Ts) && !isnan(s2.Ts) && s1.Ts ≉ s2.Ts
    warn("s1-s2: Sampling time mismatch")
    throw(DomainError())
  end
  T = promote_type(T1, T2)

  a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B, s2.B)
  c = hcat(s1.C, -s2.C)
  d = s1.D - s2.D

  DSisoSs(a,b,c,d,s1.Ts)
end

-{T<:Real}(s::DSisoSs, g::T)  = DSisoSs(copy(s.A), copy(s.B), copy(s.C),
  s.D - g, s.Ts)
-{T<:Real}(g::T, s::DSisoSs)  = +(g, -s)

.-(s1::DSisoSs, s2::DSisoSs)  = -(s1, s2)
.-{T<:Real}(s::DSisoSs, g::T) = -(s, g)
.-{T<:Real}(g::T, s::DSisoSs) = +(g, -s)

function *{T1, T2}(s1::DSisoSs{T1}, s2::DSisoSs{T2})
  # Remark: "y = (s1*s2) u" implies "u -> s2 -> s1 -> y"
  if !isnan(s1.Ts) && !isnan(s2.Ts) && s1.Ts ≉ s2.Ts
    warn("s1*s2: Sampling time mismatch")
    throw(DomainError())
  end
  T = promote_type(T1, T2)

  a = vcat(hcat(s1.A, s1.B*s2.C'),
        hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B*s2.D, s2.B)
  c = hcat(s1.C, s1.D*s2.C)
  d = s1.D * s2.D

  DSisoSs(a,b,c,d,s1.Ts)
end

*{T<:Real}(s::DSisoSs, g::T)  = DSisoSs(copy(s.A), s.B*g, copy(s.C), s.D*g, s.Ts)
*{T<:Real}(g::T, s::DSisoSs)  = DSisoSs(copy(s.A), copy(s.B), g*s.C, g*s.D, s.Ts)

.*(s1::DSisoSs, s2::DSisoSs)  = *(s1, s2)
.*{T<:Real}(s::DSisoSs, g::T) = *(s, g)
.*{T<:Real}(g::T, s::DSisoSs) = *(g, s)

/(s1::DSisoSs, s2::DSisoSs)   = *(s1, inv(s2))

/{T<:Real}(s::DSisoSs, g::T)  = DSisoSs(copy(s.A), s.B/g, copy(s.C), s.D/g, s.Ts)
/{T<:Real}(g::T, s::DSisoSs)  = *(g, inv(s))

./(s1::DSisoSs, s2::DSisoSs)  = /(s1, s2)
./{T<:Real}(s::DSisoSs, g::T) = /(s, g)
./{T<:Real}(g::T, s::DSisoSs) = /(g, s)

function ==(s1::DSisoSs, s2::DSisoSs)
  # TODO: Implement
  !isnan(s1.Ts) && !isnan(s2.Ts) && !isapprox(s1.Ts, s2.Ts) && return false
  throw(ErrorException("==(s1,s2) for DSisoSs is not implemented"))
end

!=(s1::DSisoSs, s2::DSisoSs) = !(s1 == s2)

function isapprox(s1::DSisoSs, s2::DSisoSs)
  # TODO: Implement
  !isnan(s1.Ts) && !isnan(s2.Ts) && !isapprox(s1.Ts, s2.Ts) && return false
  throw(ErrorException("isapprox(s1,s2) for DSisoSs is not implemented"))
end
