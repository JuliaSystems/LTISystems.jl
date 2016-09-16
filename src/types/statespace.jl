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
    na, ma  = size(A,1,2)
    nb, mb  = size(B,1,2)
    nc, mc  = size(C,1,2)

    @assert na == ma          "StateSpace: A must be square"
    @assert eltype(A) <: Real "StateSpace: A must be a matrix of real numbers"
    @assert na == nb          "StateSpace: A and B must have the same number of rows"
    @assert eltype(B) <: Real "StateSpace: B must be a matrix of real numbers"
    @assert ma == mc          "StateSpace: A and C must have the same number of columns"
    @assert eltype(C) <: Real "StateSpace: C must be a matrix of real numbers"
    @assert mb == nc == 1     "StateSpace: B must have 1 column, C must have 1 row"

    new{Siso{true},Continuous{true},M1,M2,M3,Matrix{M4}}(A, B, C, fill(D,1,1),
      na, mb, nc, zero(Float64))
  end

  # Discrete-time, single-input-single-output state-space model
  @compat function (::Type{StateSpace}){M1<:AbstractMatrix,M2<:AbstractMatrix,
    M3<:AbstractMatrix, M4<:Real, M5<:Real}(A::M1, B::M2, C::M3, D::M4, Ts::M5)
    na, ma  = size(A,1,2)
    nb, mb  = size(B,1,2)
    nc, mc  = size(C,1,2)

    @assert na == ma                    "StateSpace: A must be square"
    @assert eltype(A) <: Real           "StateSpace: A must be a matrix of real numbers"
    @assert na == nb                    "StateSpace: A and B must have the same number of rows"
    @assert eltype(B) <: Real           "StateSpace: B must be a matrix of real numbers"
    @assert ma == mc                    "StateSpace: A and C must have the same number of columns"
    @assert eltype(C) <: Real           "StateSpace: C must be a matrix of real numbers"
    @assert mb == nc == 1               "StateSpace: B must have 1 column, C must have 1 row"
    @assert Ts ≥ zero(Ts) && !isinf(Ts) "StateSpace: Ts must be non-negative number"

    new{Siso{true},Continuous{false},M1,M2,M3,Matrix{M4}}(A, B, C, fill(D,1,1),
      na, mb, nc, convert(Float64, Ts))
  end

  # Continuous-time, multi-input-multi-output state-space model
  @compat function (::Type{StateSpace}){M1<:AbstractMatrix,M2<:AbstractMatrix,
    M3<:AbstractMatrix, M4<:AbstractMatrix}(A::M1, B::M2, C::M3, D::M4)
    na, ma  = size(A,1,2)
    nb, mb  = size(B,1,2)
    nc, mc  = size(C,1,2)
    nd, md  = size(D,1,2)

    @assert na == ma          "StateSpace: A must be square"
    @assert eltype(A) <: Real "StateSpace: A must be a matrix of real numbers"
    @assert na == nb          "StateSpace: A and B must have the same number of rows"
    @assert eltype(B) <: Real "StateSpace: B must be a matrix of real numbers"
    @assert ma == mc          "StateSpace: A and C must have the same number of columns"
    @assert eltype(C) <: Real "StateSpace: C must be a matrix of real numbers"
    @assert nc == nd          "StateSpace: C and D must have the same number of rows"
    @assert mb == md          "StateSpace: B and D must have the same number of columns"
    @assert eltype(D) <: Real "StateSpace: D must be a matrix of real numbers"

    new{Siso{false},Continuous{true},M1,M2,M3,M4}(A, B, C, D, na, mb, nc,
      zero(Float64))
  end

  # Discrete-time, multi-input-multi-output state-space model
  @compat function (::Type{StateSpace}){M1<:AbstractMatrix,M2<:AbstractMatrix,
    M3<:AbstractMatrix, M4<:AbstractMatrix, M5<:Real}(A::M1, B::M2, C::M3,
    D::M4, Ts::M5)
    na, ma  = size(A,1,2)
    nb, mb  = size(B,1,2)
    nc, mc  = size(C,1,2)
    nd, md  = size(D,1,2)

    @assert na == ma                    "StateSpace: A must be square"
    @assert eltype(A) <: Real           "StateSpace: A must be a matrix of real numbers"
    @assert na == nb                    "StateSpace: A and B must have the same number of rows"
    @assert eltype(B) <: Real           "StateSpace: B must be a matrix of real numbers"
    @assert ma == mc                    "StateSpace: A and C must have the same number of columns"
    @assert eltype(C) <: Real           "StateSpace: C must be a matrix of real numbers"
    @assert nc == nd                    "StateSpace: C and D must have the same number of rows"
    @assert mb == md                    "StateSpace: B and D must have the same number of columns"
    @assert eltype(D) <: Real           "StateSpace: D must be a matrix of real numbers"
    @assert Ts ≥ zero(Ts) && !isinf(Ts) "StateSpace: Ts must be non-negative number"

    new{Siso{false},Continuous{false},M1,M2,M3,M4}(A, B, C, D, na, mb, nc,
      convert(Float64, Ts))
  end
end

# Legacy CMimoSs
# immutable CMimoSs{T<:Real,M1<:AbstractMatrix{T},M2<:AbstractMatrix{T},
#   M3<:AbstractMatrix{T},M4<:AbstractMatrix{T}} <: MimoSystem
#   A::M1
#   B::M2
#   C::M3
#   D::M4
#   nx::Int
#   nu::Int
#   ny::Int
#
#   function call{M1<:AbstractMatrix,M2<:AbstractMatrix,M3<:AbstractMatrix,
#     M4<:AbstractMatrix}(::Type{CMimoSs}, A::M1, B::M2, C::M3, D::M4)
#     @assert eltype(A) <: Real string("A must be a matrix of T<:Real elements")
#     @assert eltype(B) <: Real string("B must be a matrix of T<:Real elements")
#     @assert eltype(C) <: Real string("C must be a matrix of T<:Real elements")
#     @assert eltype(D) <: Real string("D must be a matrix of T<:Real elements")
#
#     na, ma  = size(A,1,2)
#     nb, mb  = size(B,1,2)
#     nc, mc  = size(C,1,2)
#
#     d       = isempty(D) ? sparse(Int[],Int[],Int8[],nc,mb) : D
#     M5      = typeof(d)
#     nd, md  = size(d,1,2)
#
#     T       = promote_type(eltype(A), eltype(B), eltype(C), eltype(d))
#
#     if na != ma
#       warn("A must be square")
#       throw(DomainError())
#     elseif nb != na
#       warn("B must have the same row size as that of A")
#       throw(DomainError())
#     elseif mc != ma
#       warn("C must have the same column size as that of A")
#       throw(DomainError())
#     elseif md != mb
#       warn("D must have the same column size as that of B")
#       throw(DomainError())
#     elseif nd != nc
#       warn("D must have the same row size as that of C")
#       throw(DomainError())
#     end
#
#     new{T,M1,M2,M3,M5}(A, B, C, d, na, mb, nc)
#   end
# end
#
# # interface implementation
# isdiscrete(s::CMimoSs)          = false
# isdiscrete(::Type{CMimoSs})     = false
# samplingtime(s::CMimoSs)        = NaN64
#
# # I/O mapping
# numstates(s::CMimoSs)           = s.nx
# numinputs(s::CMimoSs)           = s.nu
# numoutputs(s::CMimoSs)          = s.ny
#
# # Dimension information
# ndims(s::CMimoSs)               = 2
# size(s::CMimoSs)                = size(s.D)
# size(s::CMimoSs, dim::Int)      = size(s.D, dim)
# size(s::CMimoSs, dims::Int...)  = size(s.D, dims)
#
# # overload iteration interface
# done(s::CMimoSs, state::Int)                          = done(s.D, state)
# eltype{T,M1,M2,M3,M4}(::Type{CMimoSs{T,M1,M2,M3,M4}}) = CSisoSs{T}
# length(s::CMimoSs)                                    = length(s.D)
# eachindex(s::CMimoSs)                                 = eachindex(s.D)
# endof(s::CMimoSs)                                     = endof(s.D)
#
# # overload slicing functions
# function getindex(s::CMimoSs, idx::Int)
#   if idx < 1 || idx > length(s.D)
#     warn("s[idx]: Trying to access idx < 1 or idx > length(s.D)")
#     throw(BoundsError(s.D, idx))
#   end
#
#   col, row = divrem(idx-1, s.ny)
#   col += 1
#   row += 1
#
#   CSisoSs(s.A, sub(s.B, :, col:col), sub(s.C, row:row, :),
#     sub(s.D, row:row, col:col))
# end
#
# function getindex(s::CMimoSs, row::Int, col::Int)
#   if row < 1 || row > s.ny
#     warn("s[i,]: Trying to access non-existent outputs")
#     throw(BoundsError(s.C, row))
#   elseif col < 1 || col > s.nu
#     warn("s[,j]: Trying to access non-existent inputs")
#     throw(BoundsError(s.B, col))
#   end
#
#   CSisoSs(s.A, sub(s.B, :, col:col), sub(s.C, row:row, :),
#     sub(s.D, row:row, col:col))
# end
#
# getindex(s::CMimoSs, ::Colon, ::Colon)  = s
#
# getindex(s::CMimoSs, rows, cols)        = CMimoSs(s.A, sub(s.B, :, cols),
#   sub(s.C, rows, :), sub(s.D, rows, cols))
#
# getindex(s::CMimoSs, ::Colon, cols)     = CMimoSs(s.A, sub(s.B, :, cols), s.C,
#   sub(s.D, :, cols))
#
# getindex(s::CMimoSs, rows, ::Colon)     = CMimoSs(s.A, s.B, sub(s.C, rows, :),
#   sub(s.D, rows, :))
#
# # overload printing functions
# summary(s::CMimoSs) = string("ss(nx=", s.nx, ",nu=", s.nu, ",ny=", s.ny, ")")
#
# showcompact(io::IO, s::CMimoSs) = print(io, summary(s))
#
# function show(io::IO, s::CMimoSs)
#   println(io, "Continuous time state space model")
#   println(io, "\tẋ = Ax + Bu")
#   println(io, "\ty = Cx + Du")
#   print(io, "with nx=", s.nx, ", nu=", s.nu, ", ny=", s.ny, ".")
# end
#
# function showall(io::IO, s::CMimoSs)
#   show(io, s)
#   println(io, "System matrix (A):")
#   println(io, s.A)
#   println(io, "Input matrix (B):")
#   println(io, s.B)
#   println(io, "Output matrix (C):")
#   println(io, s.C)
#   println(io, "Feedforward matrix (D):")
#   print(io, s.D)
# end
#
# # creation of continuous state space types
# ss{T1<:AbstractMatrix, T2<:AbstractMatrix, T3<:AbstractMatrix,
#   T4<:AbstractMatrix}(A::T1, B::T2, C::T3, D::T4) = CMimoSs(A, B, C, D)
#
# ss{T<:AbstractMatrix}(g::T) = CMimoSs(zeros(Int8,0,0), zeros(Int8,0,size(g,2)),
#   zeros(Int8,size(g,1),0), g)
#
# # conversion and promotion
# promote_rule{T1,M1,M2,M3,M4,T2<:Real}(::Type{CMimoSs{T1,M1,M2,M3,M4}},
#   ::Type{T2}) = CMimoSs
# convert{T<:Real}(::Type{CMimoSs}, g::T) = ss(fill(g,1,1))
#
# promote_rule{T1,M1,M2,M3,M4,T2<:AbstractMatrix}(::Type{CMimoSs{T1,M1,M2,M3,M4}},
#   ::Type{T2}) = CMimoSs
# convert{T<:AbstractMatrix}(::Type{CMimoSs}, g::T) = ss(g)
#
# # overloading identities
# one{T}(s::CMimoSs{T})                               = ss(ones(T,1,1))
# one{T,M1,M2,M3,M4}(::Type{CMimoSs{T,M1,M2,M3,M4}})  = ss(ones(T,1,1))
# zero{T}(s::CMimoSs{T})                              = ss(zeros(T,1,1))
# zero{T,M1,M2,M3,M4}(::Type{CMimoSs{T,M1,M2,M3,M4}}) = ss(zeros(T,1,1))
#
# # overload inv and zeros
# function inv{T}(s::CMimoSs{T})
#   ny, nu = size(s.D, 1, 2)
#   @assert ny == nu string("inv(sys): D must be square")
#
#   try
#     Dinv = inv(s.D);
#     Ainv = s.A - s.B*Dinv*s.C;
#     Binv = s.B*Dinv
#     Cinv = -Dinv*s.C
#
#     CMimoSs(Ainv, Binv, Cinv, Dinv)
#   catch
#     warn("inv(sys): D is not invertible")
#     throw(DomainError())
#   end
# end
#
# function zeros{T}(s::CMimoSs{T})
#   Ar, Br, Cr, Dr, mr, nr, pr        = reduce(s.A, s.B, s.C, s.D)
#   if nr == 0
#     return (Complex{Float64}[], mr::Int)
#   end
#   Arc, Brc, Crc, Drc, mrc, nrc, prc = reduce(Ar.', Cr.', Br.', Dr.')
#   if nrc == 0
#     return (Complex{Float64}[], mrc::Int)
#   end
#
#   svdobj  = svdfact([Crc Drc], thin = false)
#   W       = flipdim(svdobj.Vt', 2)
#   Af      = [Arc Brc]*W[:, 1:nrc]
#
#   if mrc == 0
#     zerovalues = eigfact(Af).values
#     # return (zerovalues::Vector{Complex{Float64}}, mrc::Int)
#     return zerovalues
#   else
#     Bf    = W[1:nrc,1:nrc]
#     zerovalues = eigfact(Af, Bf).values
#     # return (zerovalues::Vector{Complex{Float64}}, mrc::Int)
#     return zerovalues
#   end
# end
#
# function poles{T}(s::CMimoSs{T})
#   Am, Bm, Cm, Dm, = minreal(s.A, s.B, s.C, s.D)
#   return eigfact(Am).values
# end
#
# # overload mathematical operations
# -(s::CMimoSs) = CMimoSs(s.A, s.B, -s.C, -s.D)
#
# function +{T1, T2}(s1::CMimoSs{T1}, s2::CMimoSs{T2})
#   # Ensure systems have same shapes
#   if size(s1) != size(s2)
#     warn("s1+s2: size(s1) != size(s2)")
#     throw(DomainError())
#   end
#
#   T = promote_type(T1, T2)
#
#   a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
#         hcat(zeros(T, s2.nx, s1.nx), s2.A))
#   b = vcat(s1.B, s2.B)
#   c = hcat(s1.C, s2.C)
#   d = s1.D + s2.D
#
#   CMimoSs(a,b,c,d)
# end
#
# +{T<:Real}(s::CMimoSs, g::T)  = CMimoSs(copy(s.A), copy(s.B), copy(s.C), s.D + g)
# +{T<:Real}(g::T, s::CMimoSs)  = +(s, g)
#
# function +{T<:AbstractMatrix}(s::CMimoSs, g::T)
#   if size(s.D) != size(g)
#     warn("s+g: size(s.D) != size(g)")
#     throw(DomainError())
#   end
#
#   CMimoSs(copy(s.A), copy(s.B), copy(s.C), s.D + g)
# end
# +{T<:AbstractMatrix}(g::T, s::CMimoSs) = +(s, g)
#
# .+(s1::CMimoSs, s2::CMimoSs)  = +(s1, s2)
# .+{T<:Real}(s::CMimoSs, g::T) = +(s, g)
# .+{T<:Real}(g::T, s::CMimoSs) = +(s, g)
#
# function -{T1, T2}(s1::CMimoSs{T1}, s2::CMimoSs{T2})
#   # Ensure systems have same shapes
#   if size(s1) != size(s2)
#     warn("s1-s2: size(s1) != size(s2)")
#     throw(DomainError())
#   end
#
#   T = promote_type(T1, T2)
#
#   a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
#         hcat(zeros(T, s2.nx, s1.nx), s2.A))
#   b = vcat(s1.B, s2.B)
#   c = hcat(s1.C, -s2.C)
#   d = s1.D - s2.D
#
#   CMimoSs(a,b,c,d)
# end
#
# -{T<:Real}(s::CMimoSs, g::T)  = CMimoSs(copy(s.A), copy(s.B), copy(s.C), s.D - g)
# -{T<:Real}(g::T, s::CMimoSs)  = +(g, -s)
#
# function -{T<:AbstractMatrix}(s::CMimoSs, g::T)
#   if size(s.D) != size(g)
#     warn("s-g: size(s.D) != size(g)")
#     throw(DomainError())
#   end
#
#   CMimoSs(copy(s.A), copy(s.B), copy(s.C), s.D - g)
# end
# -{T<:AbstractMatrix}(g::T, s::CMimoSs) = +(g, -s)
#
# .-(s1::CMimoSs, s2::CMimoSs)  = -(s1, s2)
# .-{T<:Real}(s::CMimoSs, g::T) = -(s, g)
# .-{T<:Real}(g::T, s::CMimoSs) = +(g, -s)
#
# function *{T1, T2}(s1::CMimoSs{T1}, s2::CMimoSs{T2})
#   # Remark: s1*s2 implies u -> s2 -> s1 -> y
#   if s1.nu != s2.ny
#     warn("s1*s2: s1.nu != s2.ny")
#     throw(DomainError())
#   end
#
#   T = promote_type(T1, T2)
#
#   a = vcat(hcat(s1.A, s1.B*s2.C),
#         hcat(zeros(T, s2.nx, s1.nx), s2.A))
#   b = vcat(s1.B*s2.D, s2.B)
#   c = hcat(s1.C, s1.D*s2.C)
#   d = s1.D * s2.D
#
#   CMimoSs(a,b,c,d)
# end
#
# *{T<:Real}(s::CMimoSs, g::T) = CMimoSs(copy(s.A), s.B*g, copy(s.C), s.D*g)
# *{T<:Real}(g::T, s::CMimoSs) = CMimoSs(copy(s.A), copy(s.B), g*s.C, g*s.D)
#
# function *{T<:AbstractMatrix}(s::CMimoSs, g::T)
#   if s.nu != size(g, 1)
#     warn("s*g: s.nu != size(g, 1)")
#     throw(DomainError())
#   end
#
#   CMimoSs(copy(s.A), s.B*g, copy(s.C), s.D*g)
# end
#
# function *{T<:AbstractMatrix}(g::T, s::CMimoSs)
#   if s.ny != size(g, 2)
#     warn("g*s: s.ny != size(g, 2)")
#     throw(DomainError())
#   end
#
#   CMimoSs(copy(s.A), copy(s.B), g*s.C, g*s.D)
# end
#
# .*(s1::CMimoSs, s2::CMimoSs)  = *(s1, s2)
# .*{T<:Real}(s::CMimoSs, g::T) = *(s, g)
# .*{T<:Real}(g::T, s::CMimoSs) = *(g, s)
#
# /(s1::CMimoSs, s2::CMimoSs)   = *(s1, inv(s2))
#
# /{T<:Real}(s::CMimoSs, g::T)  = CMimoSs(copy(s.A), s.B/g, copy(s.C), s.D/g)
# /{T<:Real}(g::T, s::CMimoSs)  = *(g, inv(s))
#
# function /{T<:AbstractMatrix}(s::CMimoSs, g::T)
#   ginv = inv(g)
#
#   CMimoSs(copy(s.A), s.B*ginv, copy(s.C), s.D*ginv)
# end
# /{T<:AbstractMatrix}(g::T, s::CMimoSs)  = *(g, inv(s))
#
# ./(s1::CMimoSs, s2::CMimoSs)  = /(s1, s2)
# ./{T<:Real}(s::CMimoSs, g::T) = /(s, g)
# ./{T<:Real}(g::T, s::CMimoSs) = /(g, s)
#
# function ==(s1::CMimoSs, s2::CMimoSs)
#   # TODO: Implement
#   throw(ErrorException("==(s1,s2) for CMimoSs is not implemented"))
# end
#
# !=(s1::CMimoSs, s2::CMimoSs) = !(s1 == s2)
#
# function isapprox(s1::CMimoSs, s2::CMimoSs)
#   # TODO: Implement
#   throw(ErrorException("isapprox(s1,s2) CMimoSs is not implemented"))
# end

# Legacy CSisoSs
# immutable CSisoSs{T<:Real,M1<:AbstractMatrix{T},M2<:AbstractMatrix{T},
#   M3<:AbstractMatrix{T},M4<:AbstractMatrix{T}} <: SisoSs{T}
#   A::M1
#   B::M2
#   C::M3
#   D::M4
#   nx::Int
#
#   function call{M1<:AbstractMatrix,M2<:AbstractMatrix,M3<:AbstractMatrix,
#     M4<:AbstractMatrix}(::Type{CSisoSs}, A::M1, B::M2, C::M3, D::M4)
#     @assert eltype(A) <: Real string("A must be a matrix of T<:Real elements")
#     @assert eltype(B) <: Real string("B must be a matrix of T<:Real elements")
#     @assert eltype(C) <: Real string("C must be a matrix of T<:Real elements")
#     @assert eltype(D) <: Real string("D must be a matrix of T<:Real elements")
#
#     T       = promote_type(eltype(A), eltype(B), eltype(C), eltype(D))
#     na, ma  = size(A,1,2)
#     nb, mb  = size(B,1,2)
#     nc, mc  = size(C,1,2)
#     nd, md  = size(D,1,2)
#
#     if na != ma
#       warn("A must be square")
#       throw(DomainError())
#     elseif nb != na
#       warn("B must have the same row size as that of A")
#       throw(DomainError())
#     elseif mc != ma
#       warn("C must have the same column size as that of A")
#       throw(DomainError())
#     elseif md != mb
#       warn("D must have the same column size as that of B")
#       throw(DomainError())
#     elseif nd != nc
#       warn("D must have the same row size as that of C")
#       throw(DomainError())
#     end
#
#     new{T,M1,M2,M3,M4}(A, B, C, D, na)
#   end
# end
#
# # interface implementation
# isdiscrete(s::CSisoSs)          = false
# isdiscrete(::Type{CSisoSs})     = false
# samplingtime(s::CSisoSs)        = NaN64
#
# # I/O mapping
# numstates(s::CSisoSs)           = s.nx
#
# # overload slicing functions
# function getindex(s::CSisoSs, idx::Int)
#   if idx != 1
#     warn("s[idx]: Trying to access idx != 1")
#     throw(BoundsError(s.D, idx))
#   end
#
#   return s
# end
#
# # overload printing functions
# summary(s::CSisoSs)             = string("ss(nx=", s.nx, ")")
# showcompact(io::IO, s::CSisoSs) = print(io, summary(s))
#
# function show(io::IO, s::CSisoSs)
#   println(io, "Continuous time state space model")
#   println(io, "\tẋ = Ax + Bu")
#   println(io, "\ty = Cx + Du")
#   print(io, "with nx=", s.nx, ".")
# end
#
# function showall(io::IO, s::CSisoSs)
#   show(io, s)
#   println(io)
#   println(io, "System matrix (A):")
#   println(io, s.A)
#   println(io, "Input matrix (B):")
#   println(io, s.B)
#   println(io, "Output matrix (C):")
#   println(io, s.C)
#   println(io, "Feedforward matrix (D):")
#   print(io, s.D)
# end
#
# # creation of continuous state space types
# ss{T1<:AbstractMatrix, T2<:AbstractMatrix, T3<:AbstractMatrix, T4<:Real}(A::T1,
#   B::T2, C::T3, D::T4 = zero(Int8)) = CSisoSs(A, B, C, fill(D,1,1))
#
# ss{T<:Real}(g::T) = CSisoSs(zeros(Int8,0,0), zeros(Int8,0,1), zeros(Int8,1,0),
#   fill(g,1,1))
#
# # conversion and promotion
# promote_rule{T1,M1,M2,M3,M4,T2<:Real}(::Type{CSisoSs{T1,M1,M2,M3,M4}},
#   ::Type{T2}) = CSisoSs
# convert{T<:Real}(::Type{CSisoSs}, g::T) = ss(g)
#
# # overloading identities
# one{T}(s::CSisoSs{T})                               = ss(one(T))
# one{T,M1,M2,M3,M4}(::Type{CSisoSs{T,M1,M2,M3,M4}})  = ss(one(T))
# one(::Type{CSisoSs})                                = ss(one(Int8))
# zero{T}(s::CSisoSs{T})                              = ss(zero(T))
# zero{T,M1,M2,M3,M4}(::Type{CSisoSs{T,M1,M2,M3,M4}}) = ss(zero(T))
# zero(::Type{CSisoSs})                               = ss(zero(Int8))
#
# # overload inv and zeros
# function inv{T}(s::CSisoSs{T})
#   if s.D[1] == zero(eltype(s.D))
#     warn("inv(sys): D is not invertible")
#     throw(DomainError())
#   end
#
#   Dinv = inv(s.D);
#   Ainv = s.A - s.B*Dinv*s.C;
#   Binv = s.B*Dinv
#   Cinv = -Dinv*s.C
#
#   CSisoSs(Ainv, Binv, Cinv, Dinv)
# end
#
# function zeros{T}(s::CSisoSs{T})
#   Ar, Br, Cr, Dr, mr, nr, pr        = reduce(s.A, s.B, s.C, s.D)
#   if nr == 0
#     return (Complex{Float64}[], mr::Int)
#   end
#   Arc, Brc, Crc, Drc, mrc, nrc, prc = reduce(Ar.', Cr.', Br.', Dr.')
#   if nrc == 0
#     return (Complex{Float64}[], mrc::Int)
#   end
#
#   svdobj  = svdfact([Crc Drc], thin = false)
#   W       = flipdim(svdobj.Vt', 2)
#   Af      = [Arc Brc]*W[:, 1:nrc]
#
#   if mrc == 0
#     zerovalues = eigfact(Af).values
#     # return (zerovalues::Vector{Complex{Float64}}, mrc::Int)
#     return zerovalues
#   else
#     Bf    = W[1:nrc,1:nrc]
#     zerovalues = eigfact(Af, Bf).values
#     # return (zerovalues::Vector{Complex{Float64}}, mrc::Int)
#     return zerovalues
#   end
# end
#
# function poles{T}(s::CSisoSs{T})
#   Am, Bm, Cm, Dm, = minreal(s.A, s.B, s.C, s.D)
#   return eigfact(Am).values
# end
#
# # overload mathematical operations
# -(s::CSisoSs) = CSisoSs(s.A, s.B, -s.C, -s.D)
#
# function +{T1, T2}(s1::CSisoSs{T1}, s2::CSisoSs{T2})
#   T = promote_type(T1, T2)
#
#   a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
#         hcat(zeros(T, s2.nx, s1.nx), s2.A))
#   b = vcat(s1.B, s2.B)
#   c = hcat(s1.C, s2.C)
#   d = s1.D + s2.D
#
#   CSisoSs(a,b,c,d)
# end
#
# +{T<:Real}(s::CSisoSs, g::T)  = CSisoSs(copy(s.A), copy(s.B), copy(s.C), s.D + g)
# +{T<:Real}(g::T, s::CSisoSs)  = +(s, g)
#
# .+(s1::CSisoSs, s2::CSisoSs)  = +(s1, s2)
# .+{T<:Real}(s::CSisoSs, g::T) = +(s, g)
# .+{T<:Real}(g::T, s::CSisoSs) = +(s, g)
#
# function -{T1, T2}(s1::CSisoSs{T1}, s2::CSisoSs{T2})
#   T = promote_type(T1, T2)
#
#   a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
#         hcat(zeros(T, s2.nx, s1.nx), s2.A))
#   b = vcat(s1.B, s2.B)
#   c = hcat(s1.C, -s2.C)
#   d = s1.D - s2.D
#
#   CSisoSs(a,b,c,d)
# end
#
# -{T<:Real}(s::CSisoSs, g::T)  = CSisoSs(copy(s.A), copy(s.B), copy(s.C), s.D - g)
# -{T<:Real}(g::T, s::CSisoSs)  = +(g, -s)
#
# .-(s1::CSisoSs, s2::CSisoSs)  = -(s1, s2)
# .-{T<:Real}(s::CSisoSs, g::T) = -(s, g)
# .-{T<:Real}(g::T, s::CSisoSs) = +(g, -s)
#
# function *{T1, T2}(s1::CSisoSs{T1}, s2::CSisoSs{T2})
#   # Remark: "y = (s1*s2) u" implies "u -> s2 -> s1 -> y"
#
#   T = promote_type(T1, T2)
#
#   a = vcat(hcat(s1.A, s1.B*s2.C'),
#         hcat(zeros(T, s2.nx, s1.nx), s2.A))
#   b = vcat(s1.B*s2.D, s2.B)
#   c = hcat(s1.C, s1.D*s2.C)
#   d = s1.D * s2.D
#
#   CSisoSs(a,b,c,d)
# end
#
# *{T<:Real}(s::CSisoSs, g::T)  = CSisoSs(copy(s.A), s.B*g, copy(s.C), s.D*g)
# *{T<:Real}(g::T, s::CSisoSs)  = CSisoSs(copy(s.A), copy(s.B), g*s.C, g*s.D)
#
# .*(s1::CSisoSs, s2::CSisoSs)  = *(s1, s2)
# .*{T<:Real}(s::CSisoSs, g::T) = *(s, g)
# .*{T<:Real}(g::T, s::CSisoSs) = *(g, s)
#
# /(s1::CSisoSs, s2::CSisoSs)   = *(s1, inv(s2))
#
# /{T<:Real}(s::CSisoSs, g::T)  = CSisoSs(copy(s.A), s.B/g, copy(s.C), s.D/g)
# /{T<:Real}(g::T, s::CSisoSs)  = *(g, inv(s))
#
# ./(s1::CSisoSs, s2::CSisoSs)  = /(s1, s2)
# ./{T<:Real}(s::CSisoSs, g::T) = /(s, g)
# ./{T<:Real}(g::T, s::CSisoSs) = /(g, s)
#
# function ==(s1::CSisoSs, s2::CSisoSs)
#   # TODO: Implement
#   throw(ErrorException("==(s1,s2) for CSisoSs is not implemented"))
# end
#
# !=(s1::CSisoSs, s2::CSisoSs) = !(s1 == s2)
#
# function isapprox(s1::CSisoSs, s2::CSisoSs)
#   # TODO: Implement
#   throw(ErrorException("isapprox(s1,s2) for CSisoSs is not implemented"))
# end

# Legacy DMimoSs
# immutable DMimoSs{T<:Real,M1<:AbstractMatrix{T},M2<:AbstractMatrix{T},
#   M3<:AbstractMatrix{T},M4<:AbstractMatrix{T}} <: MimoSystem
#   A::M1
#   B::M2
#   C::M3
#   D::M4
#   nx::Int
#   nu::Int
#   ny::Int
#   Ts::Float64
#
#   function call{T1<:Real,M1<:AbstractMatrix,M2<:AbstractMatrix,M3<:AbstractMatrix,
#     M4<:AbstractMatrix}(::Type{DMimoSs}, A::M1, B::M2, C::M3, D::M4, Ts::T1)
#     @assert eltype(A) <: Real string("A must be a matrix of T<:Real elements")
#     @assert eltype(B) <: Real string("B must be a matrix of T<:Real elements")
#     @assert eltype(C) <: Real string("C must be a matrix of T<:Real elements")
#     @assert eltype(D) <: Real string("D must be a matrix of T<:Real elements")
#
#     na, ma  = size(A,1,2)
#     nb, mb  = size(B,1,2)
#     nc, mc  = size(C,1,2)
#
#     d       = isempty(D) ? sparse(Int[],Int[],Int8[],nc,mb) : D
#     M5      = typeof(d)
#     nd, md  = size(d,1,2)
#     Ts_     = convert(Float64, Ts > zero(Ts) ? Ts : NaN)
#
#     T       = promote_type(eltype(A), eltype(B), eltype(C), eltype(d))
#
#     if na != ma
#       warn("A must be square")
#       throw(DomainError())
#     elseif nb != na
#       warn("B must have the same row size as that of A")
#       throw(DomainError())
#     elseif mc != ma
#       warn("C must have the same column size as that of A")
#       throw(DomainError())
#     elseif md != mb
#       warn("D must have the same column size as that of B")
#       throw(DomainError())
#     elseif nd != nc
#       warn("D must have the same row size as that of C")
#       throw(DomainError())
#     end
#
#     new{T,M1,M2,M3,M5}(A, B, C, d, na, mb, nc, Ts_)
#   end
# end
#
# # interface implementation
# isdiscrete(s::DMimoSs)          = true
# isdiscrete(::Type{DMimoSs})     = true
# samplingtime(s::DMimoSs)        = s.Ts
#
# # I/O mapping
# numstates(s::DMimoSs)           = s.nx
# numinputs(s::DMimoSs)           = s.nu
# numoutputs(s::DMimoSs)          = s.ny
#
# # Dimension information
# ndims(s::DMimoSs)               = 2
# size(s::DMimoSs)                = size(s.D)
# size(s::DMimoSs, dim::Int)      = size(s.D, dim)
# size(s::DMimoSs, dims::Int...)  = size(s.D, dims)
#
# # overload iteration interface
# done(s::DMimoSs, state::Int)                          = done(s.D, state)
# eltype{T,M1,M2,M3,M4}(::Type{DMimoSs{T,M1,M2,M3,M4}}) = DSisoSs{T}
# length(s::DMimoSs)                                    = length(s.D)
# eachindex(s::DMimoSs)                                 = eachindex(s.D)
# endof(s::DMimoSs)                                     = endof(s.D)
#
# # overload slicing functions
# function getindex(s::DMimoSs, idx::Int)
#   if idx < 1 || idx > length(s.D)
#     warn("s[idx]: Trying to access idx < 1 or idx > length(s.D)")
#     throw(BoundsError(s.D, idx))
#   end
#
#   col, row = divrem(idx-1, s.ny)
#   col += 1
#   row += 1
#
#   DSisoSs(s.A, sub(s.B, :, col:col), sub(s.C, row:row, :),
#     sub(s.D, row:row, col:col), s.Ts)
# end
#
# function getindex(s::DMimoSs, row::Int, col::Int)
#   if row < 1 || row > s.ny
#     warn("s[i,]: Trying to access non-existent outputs")
#     throw(BoundsError(s.C, row))
#   elseif col < 1 || col > s.nu
#     warn("s[,j]: Trying to access non-existent inputs")
#     throw(BoundsError(s.B, col))
#   end
#
#   DSisoSs(s.A, sub(s.B, :, col:col), sub(s.C, row:row, :),
#     sub(s.D, row:row, col:col), s.Ts)
# end
#
# getindex(s::DMimoSs, ::Colon, ::Colon)  = s
#
# getindex(s::DMimoSs, rows, cols)        = DMimoSs(s.A, sub(s.B, :, cols),
#   sub(s.C, rows, :), sub(s.D, rows, cols), s.Ts)
#
# getindex(s::DMimoSs, ::Colon, cols)     = DMimoSs(s.A, sub(s.B, :, cols), s.C,
#   sub(s.D, :, cols), s.Ts)
#
# getindex(s::DMimoSs, rows, ::Colon)     = DMimoSs(s.A, s.B, sub(s.C, rows, :),
#   sub(s.D, rows, :), s.Ts)
#
# # overload printing functions
# summary(s::DMimoSs) = string("ss(nx=", s.nx, ",nu=", s.nu, ",ny=", s.ny, ",Ts=",
#   s.Ts, ")")
#
# showcompact(io::IO, s::DMimoSs) = print(io, summary(s))
#
# function show(io::IO, s::DMimoSs)
#   println(io, "Discrete time state space model")
#   println(io, "\tx[k+1] = Ax[k] + Bu[k]")
#   println(io, "\ty[k]   = Cx[k] + Du[k]")
#   print(io, "with nx=", s.nx, ", nu=", s.nu, ", ny=", s.ny, ", Ts=", s.Ts, ".")
# end
#
# function showall(io::IO, s::DMimoSs)
#   show(io, s)
#   println(io, "System matrix (A):")
#   println(io, s.A)
#   println(io, "Input matrix (B):")
#   println(io, s.B)
#   println(io, "Output matrix (C):")
#   println(io, s.C)
#   println(io, "Feedforward matrix (D):")
#   print(io, s.D)
# end
#
# # creation of continuous state space types
# ss{T1<:AbstractMatrix, T2<:AbstractMatrix, T3<:AbstractMatrix,
#   T4<:AbstractMatrix, T5<:Real}(A::T1, B::T2, C::T3, D::T4, Ts::T5) =
#   DMimoSs(A, B, C, D, Ts)
#
# ss{T1<:AbstractMatrix, T2<:Real}(g::T1, Ts::T2) = DMimoSs(zeros(Int8,0,0),
#   zeros(Int8,0,size(g,2)), zeros(Int8,size(g,1),0), g, Ts)
#
# # conversion and promotion
# promote_rule{T1,M1,M2,M3,M4,T2<:Real}(::Type{DMimoSs{T1,M1,M2,M3,M4}},
#   ::Type{T2}) = DMimoSs
# convert{T<:Real}(::Type{DMimoSs}, g::T) = ss(fill(g,1,1), 0)
#
# promote_rule{T1,M1,M2,M3,M4,T2<:AbstractMatrix}(::Type{DMimoSs{T1,M1,M2,M3,M4}},
#   ::Type{T2}) = DMimoSs
# convert{T<:AbstractMatrix}(::Type{DMimoSs}, g::T) = ss(g, 0)
#
# # overloading identities
# one{T}(s::DMimoSs{T})                               = ss(ones(T,1,1), 0)
# one{T,M1,M2,M3,M4}(::Type{DMimoSs{T,M1,M2,M3,M4}})  = ss(ones(T,1,1), 0)
# zero{T}(s::DMimoSs{T})                              = ss(zeros(T,1,1), 0)
# zero{T,M1,M2,M3,M4}(::Type{DMimoSs{T,M1,M2,M3,M4}}) = ss(zeros(T,1,1), 0)
#
# # overload inv and zeros
# function inv{T}(s::DMimoSs{T})
#   ny, nu = size(s.D, 1, 2)
#   @assert ny == nu string("inv(sys): D must be square")
#
#   try
#     Dinv = inv(s.D);
#     Ainv = s.A - s.B*Dinv*s.C;
#     Binv = s.B*Dinv
#     Cinv = -Dinv*s.C
#
#     DMimoSs(Ainv, Binv, Cinv, Dinv, s.Ts)
#   catch
#     warn("inv(sys): D is not invertible")
#     throw(DomainError())
#   end
# end
#
# function zeros{T}(s::DMimoSs{T})
#   Ar, Br, Cr, Dr, mr, nr, pr        = reduce(s.A, s.B, s.C, s.D)
#   if nr == 0
#     return (Complex{Float64}[], mr::Int)
#   end
#   Arc, Brc, Crc, Drc, mrc, nrc, prc = reduce(Ar.', Cr.', Br.', Dr.')
#   if nrc == 0
#     return (Complex{Float64}[], mrc::Int)
#   end
#
#   svdobj  = svdfact([Crc Drc], thin = false)
#   W       = flipdim(svdobj.Vt', 2)
#   Af      = [Arc Brc]*W[:, 1:nrc]
#
#   if mrc == 0
#     zerovalues = eigfact(Af).values
#     # return (zerovalues::Vector{Complex{Float64}}, mrc::Int)
#     return zerovalues
#   else
#     Bf    = W[1:nrc,1:nrc]
#     zerovalues = eigfact(Af, Bf).values
#     # return (zerovalues::Vector{Complex{Float64}}, mrc::Int)
#     return zerovalues
#   end
# end
#
# function poles{T}(s::DMimoSs{T})
#   Am, Bm, Cm, Dm, = minreal(s.A, s.B, s.C, s.D)
#   return eigfact(Am).values
# end
#
# # overload mathematical operations
# -(s::DMimoSs) = DMimoSs(s.A, s.B, -s.C, -s.D)
#
# function +{T1, T2}(s1::DMimoSs{T1}, s2::DMimoSs{T2})
#   # Ensure systems have same sampling times and shapes
#   if !isnan(s1.Ts) && !isnan(s2.Ts) && s1.Ts ≉ s2.Ts
#     warn("s1+s2: Sampling time mismatch")
#     throw(DomainError())
#   elseif size(s1) != size(s2)
#     warn("s1+s2: size(s1) != size(s2)")
#     throw(DomainError())
#   end
#
#   T = promote_type(T1, T2)
#
#   a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
#         hcat(zeros(T, s2.nx, s1.nx), s2.A))
#   b = vcat(s1.B, s2.B)
#   c = hcat(s1.C, s2.C)
#   d = s1.D + s2.D
#
#   DMimoSs(a,b,c,d,s1.Ts)
# end
#
# +{T<:Real}(s::DMimoSs, g::T)  = DMimoSs(copy(s.A), copy(s.B), copy(s.C),
#   s.D + g, s.Ts)
# +{T<:Real}(g::T, s::DMimoSs)  = +(s, g)
#
# function +{T<:AbstractMatrix}(s::DMimoSs, g::T)
#   # Ensure systems have same shapes
#   if size(s.D) != size(g)
#     warn("s+g: size(s.D) != size(g)")
#     throw(DomainError())
#   end
#
#   DMimoSs(copy(s.A), copy(s.B), copy(s.C), s.D + g, s.Ts)
# end
# +{T<:AbstractMatrix}(g::T, s::DMimoSs) = +(s, g)
#
# .+(s1::DMimoSs, s2::DMimoSs)  = +(s1, s2)
# .+{T<:Real}(s::DMimoSs, g::T) = +(s, g)
# .+{T<:Real}(g::T, s::DMimoSs) = +(s, g)
#
# function -{T1, T2}(s1::DMimoSs{T1}, s2::DMimoSs{T2})
#   # Ensure systems have same sampling times and shapes
#   if !isnan(s1.Ts) && !isnan(s2.Ts) && s1.Ts ≉ s2.Ts
#     warn("s1-s2: Sampling time mismatch")
#     throw(DomainError())
#   elseif size(s1) != size(s2)
#     warn("s1-s2: size(s1) != size(s2)")
#     throw(DomainError())
#   end
#
#   T = promote_type(T1, T2)
#
#   a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
#         hcat(zeros(T, s2.nx, s1.nx), s2.A))
#   b = vcat(s1.B, s2.B)
#   c = hcat(s1.C, -s2.C)
#   d = s1.D - s2.D
#
#   DMimoSs(a,b,c,d,s1.Ts)
# end
#
# -{T<:Real}(s::DMimoSs, g::T)  = DMimoSs(copy(s.A), copy(s.B), copy(s.C),
#   s.D - g, s.Ts)
# -{T<:Real}(g::T, s::DMimoSs)  = +(g, -s)
#
# function -{T<:AbstractMatrix}(s::DMimoSs, g::T)
#   if size(s.D) != size(g)
#     warn("s-g: size(s.D) != size(g)")
#     throw(DomainError())
#   end
#
#   DMimoSs(copy(s.A), copy(s.B), copy(s.C), s.D - g, s.Ts)
# end
# -{T<:AbstractMatrix}(g::T, s::DMimoSs) = +(g, -s)
#
# .-(s1::DMimoSs, s2::DMimoSs)  = -(s1, s2)
# .-{T<:Real}(s::DMimoSs, g::T) = -(s, g)
# .-{T<:Real}(g::T, s::DMimoSs) = +(g, -s)
#
# function *{T1, T2}(s1::DMimoSs{T1}, s2::DMimoSs{T2})
#   # Remark: s1*s2 implies u -> s2 -> s1 -> y
#   if !isnan(s1.Ts) && !isnan(s2.Ts) && s1.Ts ≉ s2.Ts
#     warn("s1*s2: Sampling time mismatch")
#     throw(DomainError())
#   elseif s1.nu != s2.ny
#     warn("s1*s2: s1.nu != s2.ny")
#     throw(DomainError())
#   end
#
#   T = promote_type(T1, T2)
#
#   a = vcat(hcat(s1.A, s1.B*s2.C),
#         hcat(zeros(T, s2.nx, s1.nx), s2.A))
#   b = vcat(s1.B*s2.D, s2.B)
#   c = hcat(s1.C, s1.D*s2.C)
#   d = s1.D * s2.D
#
#   DMimoSs(a,b,c,d,s1.Ts)
# end
#
# *{T<:Real}(s::DMimoSs, g::T) = DMimoSs(copy(s.A), s.B*g, copy(s.C), s.D*g, s.Ts)
# *{T<:Real}(g::T, s::DMimoSs) = DMimoSs(copy(s.A), copy(s.B), g*s.C, g*s.D, s.Ts)
#
# function *{T<:AbstractMatrix}(s::DMimoSs, g::T)
#   if s.nu != size(g, 1)
#     warn("s*g: s.nu != size(g, 1)")
#     throw(DomainError())
#   end
#
#   DMimoSs(copy(s.A), s.B*g, copy(s.C), s.D*g, s.Ts)
# end
#
# function *{T<:AbstractMatrix}(g::T, s::DMimoSs)
#   if s.ny != size(g, 2)
#     warn("g*s: s.ny != size(g, 2)")
#     throw(DomainError())
#   end
#
#   DMimoSs(copy(s.A), copy(s.B), g*s.C, g*s.D, s.Ts)
# end
#
# .*(s1::DMimoSs, s2::DMimoSs)  = *(s1, s2)
# .*{T<:Real}(s::DMimoSs, g::T) = *(s, g)
# .*{T<:Real}(g::T, s::DMimoSs) = *(g, s)
#
# /(s1::DMimoSs, s2::DMimoSs)   = *(s1, inv(s2))
#
# /{T<:Real}(s::DMimoSs, g::T)  = DMimoSs(copy(s.A), s.B/g, copy(s.C), s.D/g, s.Ts)
# /{T<:Real}(g::T, s::DMimoSs)  = *(g, inv(s))
#
# function /{T<:AbstractMatrix}(s::DMimoSs, g::T)
#   ginv = inv(g)
#
#   DMimoSs(copy(s.A), s.B*ginv, copy(s.C), s.D*ginv)
# end
# /{T<:AbstractMatrix}(g::T, s::DMimoSs)  = *(g, inv(s))
#
# ./(s1::DMimoSs, s2::DMimoSs)  = /(s1, s2)
# ./{T<:Real}(s::DMimoSs, g::T) = /(s, g)
# ./{T<:Real}(g::T, s::DMimoSs) = /(g, s)
#
# function ==(s1::DMimoSs, s2::DMimoSs)
#   # TODO: Implement
#   !isnan(s1.Ts) && !isnan(s2.Ts) && !isapprox(s1.Ts, s2.Ts) && return false
#   throw(ErrorException("==(s1,s2) for DMimoSs is not implemented"))
# end
#
# !=(s1::DMimoSs, s2::DMimoSs) = !(s1 == s2)
#
# function isapprox(s1::DMimoSs, s2::DMimoSs)
#   # TODO: Implement
#   !isnan(s1.Ts) && !isnan(s2.Ts) && !isapprox(s1.Ts, s2.Ts) && return false
#   throw(ErrorException("isapprox(s1,s2) DMimoSs is not implemented"))
# end

# Legacy DSisoSs
# immutable DSisoSs{T1<:Real,M1<:AbstractMatrix{T1},M2<:AbstractMatrix{T1},
#   M3<:AbstractMatrix{T1},M4<:AbstractMatrix{T1}} <: SisoSs{T1}
#   A::M1
#   B::M2
#   C::M3
#   D::M4
#   nx::Int
#   Ts::Float64
#
#   function call{T1<:Real,M1<:AbstractMatrix,M2<:AbstractMatrix,M3<:AbstractMatrix,
#     M4<:AbstractMatrix}(::Type{DSisoSs}, A::M1, B::M2, C::M3, D::M4, Ts::T1)
#     @assert eltype(A) <: Real string("A must be a matrix of T<:Real elements")
#     @assert eltype(B) <: Real string("B must be a matrix of T<:Real elements")
#     @assert eltype(C) <: Real string("C must be a matrix of T<:Real elements")
#     @assert eltype(D) <: Real string("D must be a matrix of T<:Real elements")
#
#     T       = promote_type(eltype(A), eltype(B), eltype(C), eltype(D))
#     na, ma  = size(A,1,2)
#     nb, mb  = size(B,1,2)
#     nc, mc  = size(C,1,2)
#     nd, md  = size(D,1,2)
#     Ts_     = convert(Float64, Ts > zero(Ts) ? Ts : NaN)
#
#     if na != ma
#       warn("A must be square")
#       throw(DomainError())
#     elseif nb != na
#       warn("B must have the same row size as that of A")
#       throw(DomainError())
#     elseif mc != ma
#       warn("C must have the same column size as that of A")
#       throw(DomainError())
#     elseif md != mb
#       warn("D must have the same column size as that of B")
#       throw(DomainError())
#     elseif nd != nc
#       warn("D must have the same row size as that of C")
#       throw(DomainError())
#     end
#
#     new{T,M1,M2,M3,M4}(A, B, C, D, na, Ts_)
#   end
# end
#
# # interface implementation
# isdiscrete(s::DSisoSs)      = true
# isdiscrete(::Type{DSisoSs}) = true
# samplingtime(s::DSisoSs)    = s.Ts
#
# # I/O mapping
# numstates(s::DSisoSs)       = s.nx
#
# # overload slicing functions
# function getindex(s::DSisoSs, idx::Int)
#   if idx != 1
#     warn("s[idx]: Trying to access idx != 1")
#     throw(BoundsError(s.D, idx))
#   end
#
#   return s
# end
#
# # overload printing functions
# summary(s::DSisoSs)             = string("ss(nx=", s.nx, ",Ts=", s.Ts, ")")
# showcompact(io::IO, s::DSisoSs) = print(io, summary(s))
#
# function show(io::IO, s::DSisoSs)
#   println(io, "Discrete time state space model")
#   println(io, "\tx[k+1] = Ax[k] + Bu[k]")
#   println(io, "\ty[k]   = Cx[k] + Du[k]")
#   print(io, "with nx=", s.nx, ", Ts=", s.Ts, ".")
# end
#
# function showall(io::IO, s::DSisoSs)
#   show(io, s)
#   println(io)
#   println(io, "System matrix (A):")
#   println(io, s.A)
#   println(io, "Input matrix (B):")
#   println(io, s.B)
#   println(io, "Output matrix (C):")
#   println(io, s.C)
#   println(io, "Feedforward matrix (D):")
#   print(io, s.D)
# end
#
# # creation of continuous state space types
# ss{T1<:AbstractMatrix, T2<:AbstractMatrix, T3<:AbstractMatrix, T4<:Real,
#   T5<:Real}(A::T1, B::T2, C::T3, D::T4, Ts::T5) = DSisoSs(A, B, C, fill(D,1,1), Ts)
#
# ss{T1<:Real, T2<:Real}(g::T1, Ts::T2) = DSisoSs(zeros(Int8,0,0), zeros(Int8,0,1),
#   zeros(Int8,1,0), fill(g,1,1), Ts)
#
# # conversion and promotion
# promote_rule{T1,M1,M2,M3,M4,T2<:Real}(::Type{DSisoSs{T1,M1,M2,M3,M4}},
#   ::Type{T2}) = DSisoSs
# convert{T<:Real}(::Type{DSisoSs}, g::T) = ss(g, 0)
#
# # overloading identities
# one{T}(s::DSisoSs{T})                               = ss(one(T), 0)
# one{T,M1,M2,M3,M4}(::Type{DSisoSs{T,M1,M2,M3,M4}})  = ss(one(T), 0)
# one(::Type{DSisoSs})                                = ss(one(Int8), 0)
# zero{T}(s::DSisoSs{T})                              = ss(zero(T), 0)
# zero{T,M1,M2,M3,M4}(::Type{DSisoSs{T,M1,M2,M3,M4}}) = ss(zero(T), 0)
# zero(::Type{DSisoSs})                               = ss(zero(Int8), 0)
#
# # overload inv and zeros
# function inv{T}(s::DSisoSs{T})
#   if s.D[1] == zero(eltype(s.D))
#     warn("inv(sys): D is not invertible")
#     throw(DomainError())
#   end
#
#   Dinv = inv(s.D);
#   Ainv = s.A - s.B*Dinv*s.C;
#   Binv = s.B*Dinv
#   Cinv = -Dinv*s.C
#
#   DSisoSs(Ainv, Binv, Cinv, Dinv, s.Ts)
# end
#
# function zeros{T}(s::DSisoSs{T})
#   Ar, Br, Cr, Dr, mr, nr, pr        = reduce(s.A, s.B, s.C, s.D)
#   if nr == 0
#     return (Complex{Float64}[], mr::Int)
#   end
#   Arc, Brc, Crc, Drc, mrc, nrc, prc = reduce(Ar.', Cr.', Br.', Dr.')
#   if nrc == 0
#     return (Complex{Float64}[], mrc::Int)
#   end
#
#   svdobj  = svdfact([Crc Drc], thin = false)
#   W       = flipdim(svdobj.Vt', 2)
#   Af      = [Arc Brc]*W[:, 1:nrc]
#
#   if mrc == 0
#     zerovalues = eigfact(Af).values
#     # return (zerovalues::Vector{Complex{Float64}}, mrc::Int)
#     return zerovalues
#   else
#     Bf    = W[1:nrc,1:nrc]
#     zerovalues = eigfact(Af, Bf).values
#     # return (zerovalues::Vector{Complex{Float64}}, mrc::Int)
#     return zerovalues
#   end
# end
#
# function poles{T}(s::DSisoSs{T})
#   Am, Bm, Cm, Dm, = minreal(s.A, s.B, s.C, s.D)
#   return eigfact(Am).values
# end
#
# # overload mathematical operations
# -(s::DSisoSs) = DSisoSs(s.A, s.B, -s.C, -s.D)
#
# function +{T1, T2}(s1::DSisoSs{T1}, s2::DSisoSs{T2})
#   if !isnan(s1.Ts) && !isnan(s2.Ts) && s1.Ts ≉ s2.Ts
#     warn("s1+s2: Sampling time mismatch")
#     throw(DomainError())
#   end
#   T = promote_type(T1, T2)
#
#   a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
#         hcat(zeros(T, s2.nx, s1.nx), s2.A))
#   b = vcat(s1.B, s2.B)
#   c = hcat(s1.C, s2.C)
#   d = s1.D + s2.D
#
#   DSisoSs(a,b,c,d,s1.Ts)
# end
#
# +{T<:Real}(s::DSisoSs, g::T)  = DSisoSs(copy(s.A), copy(s.B), copy(s.C),
#   s.D + g, s.Ts)
# +{T<:Real}(g::T, s::DSisoSs)  = +(s, g)
#
# .+(s1::DSisoSs, s2::DSisoSs)  = +(s1, s2)
# .+{T<:Real}(s::DSisoSs, g::T) = +(s, g)
# .+{T<:Real}(g::T, s::DSisoSs) = +(s, g)
#
# function -{T1, T2}(s1::DSisoSs{T1}, s2::DSisoSs{T2})
#   if !isnan(s1.Ts) && !isnan(s2.Ts) && s1.Ts ≉ s2.Ts
#     warn("s1-s2: Sampling time mismatch")
#     throw(DomainError())
#   end
#   T = promote_type(T1, T2)
#
#   a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
#         hcat(zeros(T, s2.nx, s1.nx), s2.A))
#   b = vcat(s1.B, s2.B)
#   c = hcat(s1.C, -s2.C)
#   d = s1.D - s2.D
#
#   DSisoSs(a,b,c,d,s1.Ts)
# end
#
# -{T<:Real}(s::DSisoSs, g::T)  = DSisoSs(copy(s.A), copy(s.B), copy(s.C),
#   s.D - g, s.Ts)
# -{T<:Real}(g::T, s::DSisoSs)  = +(g, -s)
#
# .-(s1::DSisoSs, s2::DSisoSs)  = -(s1, s2)
# .-{T<:Real}(s::DSisoSs, g::T) = -(s, g)
# .-{T<:Real}(g::T, s::DSisoSs) = +(g, -s)
#
# function *{T1, T2}(s1::DSisoSs{T1}, s2::DSisoSs{T2})
#   # Remark: "y = (s1*s2) u" implies "u -> s2 -> s1 -> y"
#   if !isnan(s1.Ts) && !isnan(s2.Ts) && s1.Ts ≉ s2.Ts
#     warn("s1*s2: Sampling time mismatch")
#     throw(DomainError())
#   end
#   T = promote_type(T1, T2)
#
#   a = vcat(hcat(s1.A, s1.B*s2.C'),
#         hcat(zeros(T, s2.nx, s1.nx), s2.A))
#   b = vcat(s1.B*s2.D, s2.B)
#   c = hcat(s1.C, s1.D*s2.C)
#   d = s1.D * s2.D
#
#   DSisoSs(a,b,c,d,s1.Ts)
# end
#
# *{T<:Real}(s::DSisoSs, g::T)  = DSisoSs(copy(s.A), s.B*g, copy(s.C), s.D*g, s.Ts)
# *{T<:Real}(g::T, s::DSisoSs)  = DSisoSs(copy(s.A), copy(s.B), g*s.C, g*s.D, s.Ts)
#
# .*(s1::DSisoSs, s2::DSisoSs)  = *(s1, s2)
# .*{T<:Real}(s::DSisoSs, g::T) = *(s, g)
# .*{T<:Real}(g::T, s::DSisoSs) = *(g, s)
#
# /(s1::DSisoSs, s2::DSisoSs)   = *(s1, inv(s2))
#
# /{T<:Real}(s::DSisoSs, g::T)  = DSisoSs(copy(s.A), s.B/g, copy(s.C), s.D/g, s.Ts)
# /{T<:Real}(g::T, s::DSisoSs)  = *(g, inv(s))
#
# ./(s1::DSisoSs, s2::DSisoSs)  = /(s1, s2)
# ./{T<:Real}(s::DSisoSs, g::T) = /(s, g)
# ./{T<:Real}(g::T, s::DSisoSs) = /(g, s)
#
# function ==(s1::DSisoSs, s2::DSisoSs)
#   # TODO: Implement
#   !isnan(s1.Ts) && !isnan(s2.Ts) && !isapprox(s1.Ts, s2.Ts) && return false
#   throw(ErrorException("==(s1,s2) for DSisoSs is not implemented"))
# end
#
# !=(s1::DSisoSs, s2::DSisoSs) = !(s1 == s2)
#
# function isapprox(s1::DSisoSs, s2::DSisoSs)
#   # TODO: Implement
#   !isnan(s1.Ts) && !isnan(s2.Ts) && !isapprox(s1.Ts, s2.Ts) && return false
#   throw(ErrorException("isapprox(s1,s2) for DSisoSs is not implemented"))
# end
