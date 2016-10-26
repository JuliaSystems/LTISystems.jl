immutable RationalTF{T,S,M1,M2} <: LtiSystem{T,S}
  num::M1
  den::M2
  nu::Int
  ny::Int
  Ts::Float64

  # Continuous-time, single-input-single-output rational transfer function model
  @compat function (::Type{RationalTF}){M1<:Poly,M2<:Poly}(num::M1, den::M2)
    n = fill(num,1,1)
    d = fill(den,1,1)
    ny, nu = tfcheck(n, d)
    new{Siso{true},Continuous{true},Matrix{M1},Matrix{M2}}(n, d, nu, ny, zero(Float64))
  end

  # Discrete-time, single-input-single-output rational transfer function model
  @compat function (::Type{RationalTF}){M1<:Poly,M2<:Poly}(num::M1, den::M2, Ts::Real)
    n = fill(num,1,1)
    d = fill(den,1,1)
    ny, nu = tfcheck(n, d, Ts)
    new{Siso{true},Continuous{false},Matrix{M1},Matrix{M2}}(n, d, nu, ny,
      convert(Float64, Ts))
  end

  # Continuous-time, multi-input-multi-output rational transfer function model
  @compat function (::Type{RationalTF}){M1<:AbstractMatrix,
    M2<:AbstractMatrix}(num::M1, den::M2)
    ny, nu = tfcheck(num, den)
    new{Siso{false},Continuous{true},M1,M2}(num, den, nu, ny, zero(Float64))
  end

  # Discrete-time, multi-input-multi-output rational transfer function model
  @compat function (::Type{RationalTF}){M1<:AbstractMatrix,
    M2<:AbstractMatrix}(num::M1, den::M2, Ts::Real)
    ny, nu = tfcheck(num, den, Ts)
    new{Siso{false},Continuous{false},M1,M2}(num, den, nu, ny, convert(Float64, Ts))
  end
end

# Enforce rational transfer function type invariance
function tfcheck(num::AbstractMatrix, den::AbstractMatrix, Ts::Real = zero(Float64))
  @assert size(num) == size(den)    "RationalTF: size(num) ≠ size(den)"
  @assert !isempty(num)             "RationalTF: min(nu, ny) = 0"
  @assert eltype(num) <: Poly &&
    eltype(num[1]) <: Real     "RationalTF: num polynomial(s) do not have real coefficients"
  @assert eltype(den) <: Poly &&
    eltype(den[1]) <: Real     "RationalTF: den polynomial(s) do not have real coefficients"
  for idx in eachindex(num)
    @assert degree(num[idx]) ≤ degree(den[idx]) "RationalTF: system is not proper"
    @assert den[idx] != zero(den[idx])          "RationalTF: den polynomial(s) cannot be zero"
  end
  @assert Ts ≥ zero(Ts) && !isinf(Ts) "RationalTF: Ts must be non-negative real number"

  return size(num)
end

# Outer constructors
tf(num::Poly, den::Poly)            = RationalTF(num, den)
tf(num::Poly, den::Poly, Ts::Real)  = RationalTF(num, den, Ts)

tf(num::AbstractMatrix, den::AbstractMatrix)            = RationalTF(num, den)
tf(num::AbstractMatrix, den::AbstractMatrix, Ts::Real)  = RationalTF(num, den, Ts)

function tf{T1<:Real, T2<:Real}(num::AbstractVector{T1}, den::AbstractVector{T2})
  numpoly = Poly(reverse(num), :s)
  denpoly = Poly(reverse(den), :s)
  RationalTF(numpoly, denpoly)
end

function tf{T1<:Real, T2<:Real}(num::AbstractVector{T1}, den::AbstractVector{T2},
  Ts::Real)
  numpoly = Poly(reverse(num), :z)
  denpoly = Poly(reverse(den), :z)
  RationalTF(numpoly, denpoly, Ts)
end

function tf{T1<:Real, T2<:Real}(num::AbstractVector{T1}, den::AbstractVector{T2},
  Ts::Real, var::Symbol)
  vars    = [:z̄,:q̄,:qinv,:zinv]
  @assert var ∈ vars string("tf: var ∉ ", vars)

  numlast         = findlast(num)
  denlast         = findlast(den)
  order           = max(numlast, denlast)
  num_            = zeros(eltype(num), order)
  num_[1:numlast] = num[1:numlast]
  den_            = zeros(eltype(den), order)
  den_[1:denlast] = den[1:denlast]

  numpoly = Poly(reverse(num_), :z)
  denpoly = Poly(reverse(den_), :z)
  RationalTF(numpoly, denpoly, Ts)
end

# Interfaces
isdiscrete{T,S}(s::RationalTF{T,Continuous{S}}) = !S

samplingtime(s::RationalTF) = s.Ts

# Think carefully about how to implement numstates
numstates(s::RationalTF)    = numstates(ss(s))
# Currently, we only allow for proper systems
numstates(s::RationalTF{Siso{true}}) = degree(den[1])

numinputs(s::RationalTF)    = s.nu
numoutputs(s::RationalTF)   = s.ny

# Dimension information
ndims(s::RationalTF{Siso{true}})  = 1
ndims(s::RationalTF{Siso{false}}) = 2
size(s::RationalTF)               = size(s.num)
size(s::RationalTF, dim::Int)     = size(s.num, dim)
size(s::RationalTF, dims::Int...) = size(s.num, dims)

# Iteration interface (meaningful only for MIMO systems)
# TODO

# Slicing (`getindex`) of MIMO systems
# TODO

# Printing functions
summary(s::RationalTF{Siso{true},Continuous{true}})   =
  string("tf(nx=", numstates(s), ")")
summary(s::RationalTF{Siso{true},Continuous{false}})  =
  string("tf(nx=", numstates(s), ",Ts=", s.Ts, ")")
summary(s::RationalTF{Siso{false},Continuous{true}})  =
  string("tf(nx=", numstates(s), ",nu=", s.nu, ",ny=", s.ny, ")")
summary(s::RationalTF{Siso{false},Continuous{false}}) =
  string("tf(nx=", numstates(s), ",nu=", s.nu, ",ny=", s.ny, ",Ts=", s.Ts, ")")

showcompact(io::IO, s::RationalTF) = print(io, summary(s))

function show{T}(io::IO, s::RationalTF{T,Continuous{true}})
  # TODO
end

function show{T}(io::IO, s::RationalTF{T,Continuous{false}})
  # TODO
end

function showall(io::IO, s::RationalTF)
  # TODO
end

# Conversion and promotion
promote_rule{T<:Real,S}(::Type{T}, ::Type{RationalTF{Siso{true},S}}) =
  RationalTF{Siso{true},S}
promote_rule{T<:AbstractMatrix,S}(::Type{T}, ::Type{RationalTF{Siso{false},S}}) =
  RationalTF{Siso{false},S}

convert(::Type{RationalTF{Siso{true},Continuous{true}}}, g::Real)             =
  tf(g)
convert(::Type{RationalTF{Siso{true},Continuous{false}}}, g::Real)            =
  tf(g, zero(Float64))
convert(::Type{RationalTF{Siso{false},Continuous{true}}}, g::AbstractMatrix)  =
  tf(g)
convert(::Type{RationalTF{Siso{false},Continuous{false}}}, g::AbstractMatrix) =
  tf(g, zero(Float64))

# Multiplicative and additive identities (meaningful only for SISO)
one(::Type{RationalTF{Siso{true},Continuous{true}}})    =
  tf(one(Int8))
one(::Type{RationalTF{Siso{true},Continuous{false}}})   =
  tf(one(Int8), zero(Float64))
zero(::Type{RationalTF{Siso{true},Continuous{true}}})   =
  tf(zero(Int8))
zero(::Type{RationalTF{Siso{true},Continuous{false}}})  =
  tf(zero(Int8), zero(Float64))

one(s::RationalTF{Siso{true},Continuous{true}})   =
  tf(one(eltype(s.num)), one(eltype(s.den)))
one(s::RationalTF{Siso{true},Continuous{false}})  =
  tf(one(eltype(s.num)), one(eltype(s.den)), zero(Float64))
zero(s::RationalTF{Siso{true},Continuous{true}})  =
  tf(zero(eltype(s.num)), one(eltype(s.den)))
zero(s::RationalTF{Siso{true},Continuous{false}}) =
  tf(zero(eltype(s.num)), one(eltype(s.den)), zero(Float64))

# Inverse of a rational transfer function model
function tfinv(s::RationalTF{Siso{true}})
  @assert degree(num[1]) == degree(den[1]) "inv(sys): inverse is not a proper system"
  return den, num
end

function tfinv(s::RationalTF{Siso{false}})
  @assert s.ny == s.nu "inv(sys): s.ny ≠ s.nu"
  # TODO
  # return num, den
end

function inv(s::RationalTF{Siso{true},Continuous{true}})
  num, den = tfinv(s)
  tf(num[1], den[1])
end

function inv(s::RationalTF{Siso{true},Continuous{false}})
  num, den = tfinv(s)
  tf(num[1], den[1], s.Ts)
end

function inv(s::RationalTF{Siso{false},Continuous{true}})
  num, den = tfinv(s)
  tf(num, den)
end

function inv(s::RationalTF{Siso{false},Continuous{false}})
  num, den = tfinv(s)
  tf(num, den, s.Ts)
end

# Invariant zeros of a rational transfer function model
zeros(s::RationalTF{Siso{true}}) = roots(s.num[1])

function zeros(s::RationalTF)
  # TODO
end

# Transmission zeros of a rational transfer function model
tzeros(s::RationalTF) = zeros(minreal(s))

# Poles of a rational transfer function model
poles(s::RationalTF{Siso{true}}) = roots(s.den[1])

function poles(s::RationalTF)
  # TODO
end

# Negative of a rational transfer function model
-(s::RationalTF{Siso{true},Continuous{true}})   =
  RationalTF(-s.num[1], s.den[1])
-(s::RationalTF{Siso{true},Continuous{false}})  =
  RationalTF(-s.num[1], s.den[1], s.Ts)
-(s::RationalTF{Siso{false},Continuous{true}})  =
  RationalTF(-s.num, s.den)
-(s::RationalTF{Siso{false},Continuous{false}}) =
  RationalTF(-s.num, s.den, s.Ts)

# Addition
function tfparallel{T1,T2,S}(s1::RationalTF{T1,S}, s2::RationalTF{T2,S})
  @assert s1.Ts ≈ s2.Ts || s1.Ts == zero(Float64) ||
    s2.Ts == zero(Float64) "parallel(s1,s2): Sampling time mismatch"
  @assert size(s1) == size(s2) "parallel(s1,s2): size(s1) ≠ size(s2)"

  # return num, den
end

function +(s1::RationalTF{Siso{true},Continuous{true}},
  s2::RationalTF{Siso{true},Continuous{true}})
  num, den = tfparallel(s1, s2)
  RationalTF(num[1], den[1])
end

function +(s1::RationalTF{Siso{true},Continuous{false}},
  s2::RationalTF{Siso{true},Continuous{false}})
  num, den = tfparallel(s1, s2)
  RationalTF(num[1], den[1], s1.Ts)
end

function +{T1,T2}(s1::RationalTF{T1,Continuous{true}},
  s2::RationalTF{T2,Continuous{true}})
  num, den = tfparallel(s1, s2)
  RationalTF(num, den)
end

function +{T1,T2}(s1::RationalTF{T1,Continuous{false}},
  s2::RationalTF{T2,Continuous{false}})
  num, den = tfparallel(s1, s2)
  RationalTF(num, den, s1.Ts)
end

.+(s1::RationalTF{Siso{true}}, s2::RationalTF{Siso{true}}) = +(s1, s2)

+{T}(s::RationalTF{T,Continuous{true}}, g)  = +(s, tf(g))
+{T}(s::RationalTF{T,Continuous{false}}, g) = +(s, tf(g, zero(Float64)))
+{T}(g, s::RationalTF{T,Continuous{true}})  = +(tf(g), s)
+{T}(g, s::RationalTF{T,Continuous{false}}) = +(tf(g, zero(Float64)), s)

.+(s::RationalTF{Siso{true}}, g::Real)  = +(s, g)
.+(g::Real, s::RationalTF{Siso{true}})  = +(g, s)

# Subtraction
-(s1::RationalTF, s2::RationalTF) = +(s1, -s2)

.-(s1::RationalTF{Siso{true}}, s2::RationalTF{Siso{true}}) = -(s1, s2)

-{T}(s::RationalTF{T,Continuous{true}}, g)  = -(s, tf(g))
-{T}(s::RationalTF{T,Continuous{false}}, g) = -(s, tf(g, zero(Float64)))
-{T}(g, s::RationalTF{T,Continuous{true}})  = -(tf(g), s)
-{T}(g, s::RationalTF{T,Continuous{false}}) = -(tf(g, zero(Float64)), s)

.-(s::RationalTF{Siso{true}}, g::Real)  = -(s, g)
.-(g::Real, s::RationalTF{Siso{true}})  = -(g, s)

# Multiplication
function tfseries{T1,T2,S}(s1::RationalTF{T1,S}, s2::RationalTF{T2,S})
  # Remark: s1*s2 implies u -> s2 -> s1 -> y
  @assert s1.Ts ≈ s2.Ts || s1.Ts == zero(Float64) ||
    s2.Ts == zero(Float64) "series(s1,s2): Sampling time mismatch"
  @assert s1.nu == s2.ny "series(s1,s2): s1.nu ≠ s2.ny"

  return num, den
end

function *(s1::RationalTF{Siso{true},Continuous{true}},
  s2::RationalTF{Siso{true},Continuous{true}})
  num, den = tfseries(s1, s2)
  RationalTF(num[1], den[1])
end

function *(s1::RationalTF{Siso{true},Continuous{false}},
  s2::RationalTF{Siso{true},Continuous{false}})
  num, den = tfseries(s1, s2)
  RationalTF(num[1], den[1], s1.Ts)
end

function *{T1,T2}(s1::RationalTF{T1,Continuous{true}},
  s2::RationalTF{T2,Continuous{true}})
  num, den = tfseries(s1, s2)
  RationalTF(num, den)
end

function *{T1,T2}(s1::RationalTF{T1,Continuous{false}},
  s2::RationalTF{T2,Continuous{false}})
  num, den = tfseries(s1, s2)
  RationalTF(num, den, s1.Ts)
end

.*(s1::RationalTF{Siso{true}}, s2::RationalTF{Siso{true}}) = *(s1, s2)

*{T}(s::RationalTF{T,Continuous{true}}, g)  = *(s, tf(g))
*{T}(s::RationalTF{T,Continuous{false}}, g) = *(s, tf(g, zero(Float64)))
*{T}(g, s::RationalTF{T,Continuous{true}})  = *(tf(g), s)
*{T}(g, s::RationalTF{T,Continuous{false}}) = *(tf(g, zero(Float64)), s)

.*(s::RationalTF{Siso{true}}, g::Real)  = *(s, g)
.*(g::Real, s::RationalTF{Siso{true}})  = *(g, s)

# Division
/(s1::RationalTF, s2::RationalTF) = *(s1, inv(s2))

./(s1::RationalTF{Siso{true}}, s2::RationalTF{Siso{true}}) = /(s1, s2)

/{T}(s::RationalTF{T,Continuous{true}}, g)  = /(s, tf(g))
/{T}(s::RationalTF{T,Continuous{false}}, g) = /(s, tf(g, zero(Float64)))
/{T}(g, s::RationalTF{T,Continuous{true}})  = /(tf(g), s)
/{T}(g, s::RationalTF{T,Continuous{false}}) = /(tf(g, zero(Float64)), s)

./(s::RationalTF{Siso{true}}, g::Real)  = /(s, g)
./(g::Real, s::RationalTF{Siso{true}})  = /(g, s)

# Legacy SisoTf
# # Printing functions
# summary(s::SisoTf) = string("tf(nx=", numstates(s), (isdiscrete(s) ?
#   string(",Ts=", samplingtime(s)) : ""), ")")
#
# showcompact(io::IO, s::SisoTf)  = print(io, summary(s))
# show(io::IO, s::SisoTf)         = print(io, summary(s))
# showall(io::IO, s::SisoTf)      = print(io, summary(s))
#
# function rmgcd{T1<:Number, T2<:Number}(p1::Poly{T1}, p2::Poly{T2})
#   R = promote_type(T1,T2)
#   gcdp1p2::Poly{R} = gcd(p1,p2)
#   gcdp1p2_::Poly{R} = gcdp1p2/gcdp1p2[end]
#   p1_::Poly = div(p1,gcdp1p2_)
#   p2_::Poly = div(p2,gcdp1p2_)
#   return (p1_,p2_,gcdp1p2_)
# end
#
# function rmcommon{T1<:AbstractVector, T2<:AbstractVector}(a::T1, b::T2)
#   T     = promote_type(eltype(T1),eltype(T2))
#   list  = Dict{T,Vector{Int}}()
#
#   l1    = length(a)
#   l2    = length(b)
#   s1    = l1 > l2 ? b : a
#   s2    = l1 > l2 ? a : b
#
#   for elem in s1
#     list[elem] = get(list, elem, Int[0, 0]) + [1, 0]
#   end
#
#   for elem in s2
#     try
#       list[elem] += [0, 1]
#     end
#   end
#
#   temp::Vector{T} = []
#   for key in keys(list)
#     append!(temp, fill(key, min(list[key]...)))
#   end
#   ar::T1 = copy(a)
#   br::T2 = copy(b)
#   for elem in temp
#     splice!(ar, findfirst(ar,elem))
#     splice!(br, findfirst(br,elem))
#   end
#   return (ar,br,temp)
# end

# Legacy CSisoRational
# immutable CSisoRational{T1<:Real, T2<:Real} <: CSisoTf
#   num::Poly{T1}
#   den::Poly{T2}
#
#   function call{T1<:Real, T2<:Real}(::Type{CSisoRational}, num::Poly{T1},
#       den::Poly{T2})
#
#     new{T1,T2}(num, den)
#   end
# end
#
# # creation of continuous rational transfer functions
#
# # conversion and promotion
#
# promote_rule{T11,T12,T21,T22}(::Type{CSisoRational{T11,T12}},
#   ::Type{CSisoRational{T21,T22}})                              = CSisoRational
# promote_rule{T<:Real,T11,T12}(::Type{CSisoRational{T11,T12}},
# ::Type{T})                                                     = CSisoRational
# promote_rule{T1,T<:Real}(::Type{CSisoRational{T1}}, ::Type{T}) = CSisoRational
#
# convert{T11,T12,T21,T22}(::Type{CSisoRational{T11,T12}},
#   s::CSisoRational{T21,T22})                                   = s
# convert{T<:Real}(::Type{CSisoRational}, x::T)                  = tf([x], [one(T)])
#
# # interface implementation
#
# zeros(s::CSisoRational)           = convert(Vector{Complex{Float64}}, roots(s.num))
# poles(s::CSisoRational)           = convert(Vector{Complex{Float64}}, roots(s.den))
# numvec(s::CSisoRational)          = reverse(coeffs(s.num))
# denvec(s::CSisoRational)          = reverse(coeffs(s.den))
# numpoly(s::CSisoRational)         = copy(s.num)
# denpoly(s::CSisoRational)         = copy(s.den)
# zpkdata(s::CSisoRational)         = (zeros(s), poles(s), s.num[end]/s.den[end])
#
# # overload printing functions
#
# function show(io::IO, s::CSisoRational)
#   println(io, "Continuous time rational transfer function model")
#   println(io, "\ty = Gu")
# end
#
# function showall(io::IO, s::CSisoRational)
#   show(io, s)
#   println(io, "")
#   printtransferfunction(io::IO, s)
# end
#
# # overload mathematical operations
#
# function +(s1::CSisoRational, s2::CSisoRational)
#   den1,den2,dengcd   = rmgcd(s1.den, s2.den)
#   tf(s1.num*den2 + s2.num*den1, den1*den2*dengcd)
# end
# +{T<:Real}(s::CSisoRational, n::T)       = tf(s.num + n*s.den, copy(s.den))
# +{T<:Real}(n::T, s::CSisoRational)       = s + n
#
# .+{T<:Real}(s::CSisoRational, n::T)      = s + n
# .+{T<:Real}(n::T, s::CSisoRational)      = s + n
# .+(s1::CSisoRational, s2::CSisoRational) = +(s1, s2)
#
# -(s::CSisoRational)                      = tf(-s.num, copy(s.den))
#
# -(s1::CSisoRational, s2::CSisoRational)  = +(s1,-s2)
# -{T<:Real}(n::T, s::CSisoRational)       = +(n, -s)
# -{T<:Real}(s::CSisoRational, n::T)       = +(s, -n)
#
# .-{T<:Real}(s::CSisoRational, n::T)      = -(s, n)
# .-{T<:Real}(n::T, s::CSisoRational)      = -(n, s)
# .-(s1::CSisoRational, s2::CSisoRational) = +(s1, -s2)
#
# function *(s1::CSisoRational, s2::CSisoRational)
#   num1,den2,gcd1   = rmgcd(s1.num, s2.den)
#   den1,num2,gcd2   = rmgcd(s1.den, s2.num)
#   tf(num1*num2, den1*den2)
# end
#
# *{T<:Real}(s::CSisoRational, n::T)       = tf(s.num*n, copy(s.den))
# *{T<:Real}(n::T, s::CSisoRational)       = *(s, n)
#
# .*{T<:Real}(s::CSisoRational, n::T)      = *(s, n)
# .*{T<:Real}(n::T, s::CSisoRational)      = *(n, s)
# .*(s1::CSisoRational, s2::CSisoRational) = *(s1, s2)
#
# /(s1::CSisoRational, s2::CSisoRational)  = s1*(1/s2)
# /{T<:Real}(n::T, s::CSisoRational)       = tf(n*s.den, copy(s.num))
# /{T<:Real}(s::CSisoRational, n::T)       = s*(1/n)
#
# ./{T<:Real}(n::T, s::CSisoRational)      = /(n, s)
# ./{T<:Real}(s::CSisoRational, n::T)      = /(s, n)
# ./(s1::CSisoRational, s2::CSisoRational) = /(s1, s2)
#
# function ==(s1::CSisoRational, s2::CSisoRational)
#   s1.num == s2.num && (s1.den == s2.den || s1.num == zero(s1.num))
# end
#
# !=(s1::CSisoRational, s2::CSisoRational) = !(s1 == s2)
#
# function isapprox(s1::CSisoRational, s2::CSisoRational,
#   rtol::Real=sqrt(eps()), atol::Real=0)
#   sdiff = s2-s1
#   norm(sdiff.num) < rtol
# end

# Legacy DSisoRational
# immutable DSisoRational{T1<:Real, T2<:Real} <: DSisoTf
#   num::Poly{T1}
#   den::Poly{T2}
#   Ts::Float64
#
#   function call{T1<:Real, T2<:Real}(::Type{DSisoRational}, num::Poly{T1},
#       den::Poly{T2}, Ts::Float64)
#
#     Ts_ = Ts > zero(Float64) ? Ts : NaN
#     new{T1,T2}(num, den, Ts_)
#   end
# end
#
# # creation of discrete rational transfer functions
#
# function tf{V1<:AbstractVector, V2<:AbstractVector, T3<:Real}(num::V1, den::V2,
#     Ts::T3; var::Symbol=:z)
#   @assert eltype(num) <: Number string("num must be vector of T<:Number elements")
#   @assert eltype(den) <: Number string("den must be vector of T<:Number elements")
#   @assert var == :q || var == :z || var == :zinv || var == :qinv
#     string("var need to be either :q, :z, :qinv or :zinv")
#
#   if var == :zinv || var == :qinv
#     numlast = findlast(num)
#     denlast = findlast(den)
#     order = max(numlast,denlast)
#     num_ = zeros(eltype(num),order)
#     num_[1:numlast] = num[1:numlast]
#     den_ = zeros(eltype(den),order)
#     den_[1:denlast] = den[1:denlast]
#     num = num_
#     den = den_
#   end
#
#   pnum = Poly(num[end:-1:1])
#   pden = Poly(den[end:-1:1])
#   DSisoRational(pnum, pden, Float64(Ts))
# end
#
# tf{T1<:Real, T2<:Real, T3<:Real}(num::Poly{T1}, den::Poly{T2}, Ts::T3) = DSisoRational(num, den, Float64(Ts))
# tf{T1<:Real, T2<:Real}(gain::T1, Ts::T2)                               = DSisoRational(Poly([gain]), Poly([one(T1)]), Float64(Ts))
#
#
# # conversion and promotion
#
# promote_rule{T11,T12,T21,T22}(::Type{DSisoRational{T11,T12}},
#   ::Type{DSisoRational{T21,T22}})                              = DSisoRational
# promote_rule{T<:Real,T11,T12}(::Type{DSisoRational{T11,T12}},
# ::Type{T})                                                     = DSisoRational
# promote_rule{T1,T<:Real}(::Type{DSisoRational{T1}}, ::Type{T}) = DSisoRational
#
# convert{T11,T12,T21,T22}(::Type{DSisoRational{T11,T12}},
#   s::DSisoRational{T21,T22})                                   = s
# convert{T<:Real}(::Type{DSisoRational}, x::T)                  = tf([x], [one(T)], NaN64)
#
# # interface implementation
#
# zeros(s::DSisoRational)           = convert(Vector{Complex{Float64}}, roots(s.num))
# poles(s::DSisoRational)           = convert(Vector{Complex{Float64}}, roots(s.den))
# numvec(s::DSisoRational)          = reverse(coeffs(s.num))
# denvec(s::DSisoRational)          = reverse(coeffs(s.den))
# numpoly(s::DSisoRational)         = copy(s.num)
# denpoly(s::DSisoRational)         = copy(s.den)
# zpkdata(s::DSisoRational)         = (zeros(s), poles(s), s.num[end]/s.den[end])
#
# # overload printing functions
#
# function show(io::IO, s::DSisoRational)
#   println(io, "Discrete time rational transfer function model")
#   println(io, "\ty = Gu")
#   if s.Ts > 0
#     println(io, "with Ts=", s.Ts, ".")
#   elseif s.Ts == 0
#     println(io, "with Ts=unspecified.")
#   end
# end
#
# function showall(io::IO, s::DSisoRational)
#   show(io, s)
#   println(io, "")
#   printtransferfunction(io::IO, s)
# end
#
# # overload mathematical operations
#
# function +(s1::DSisoRational, s2::DSisoRational)
#   Ts::Float64
#   if s1.Ts == s2.Ts || isnan(s2.Ts)
#     Ts = s1.Ts
#   elseif isnan(s1.Ts)
#     Ts = s2.Ts
#   else
#     warn("Sampling time mismatch")
#     throw(DomainError())
#   end
#   den1,den2,dengcd   = rmgcd(s1.den, s2.den)
#   tf(s1.num*den2 + s2.num*den1, den1*den2*dengcd, Ts)
# end
# +{T<:Real}(s::DSisoRational, n::T)       = tf(s.num + n*s.den, copy(s.den), s.Ts)
# +{T<:Real}(n::T, s::DSisoRational)       = s + n
#
# .+{T<:Real}(s::DSisoRational, n::T)      = s + n
# .+{T<:Real}(n::T, s::DSisoRational)      = s + n
# .+(s1::DSisoRational, s2::DSisoRational) = +(s1,-s2)
#
# -{T}(s::DSisoRational{T})                = tf(-s.num, copy(s.den), s.Ts)
#
# -(s1::DSisoRational, s2::DSisoRational)  = +(s1,-s2)
# -{T<:Real}(n::T, s::DSisoRational)       = +(n, -s)
# -{T<:Real}(s::DSisoRational, n::T)       = +(s, -n)
#
# .-{T<:Real}(s::DSisoRational, n::T)      = -(s, n)
# .-{T<:Real}(n::T, s::DSisoRational)      = -(n, s)
# .-(s1::DSisoRational, s2::DSisoRational) = -(s1, s2)
#
# function *(s1::DSisoRational, s2::DSisoRational)
#   Ts::Float64
#   if s1.Ts == s2.Ts || isnan(s2.Ts)
#     Ts = s1.Ts
#   elseif isnan(s1.Ts)
#     Ts = s2.Ts
#   else
#     warn("Sampling time mismatch")
#     throw(DomainError())
#   end
#   num1,den2,gcd1 = rmgcd(s1.num, s2.den)
#   den1,num2,gcd2 = rmgcd(s1.den, s2.num)
#   tf(num1*num2, den1*den2, Ts)
# end
# *{T<:Real}(s::DSisoRational, n::T)       = tf(s.num*n, copy(s.den), s.Ts)
# *{T<:Real}(n::T, s::DSisoRational)       = *(s, n)
#
# .*{T<:Real}(s::DSisoRational, n::T)      = *(s, n)
# .*{T<:Real}(n::T, s::DSisoRational)      = *(n, s)
# .*(s1::DSisoRational, s2::DSisoRational) = *(s1, s2)
#
# /(s1::DSisoRational, s2::DSisoRational)  = s1*(1/s2)
# /{T<:Real}(n::T, s::DSisoRational)       = tf(n*s.den, copy(s.num), s.Ts)
# /{T<:Real}(s::DSisoRational, n::T)       = s*(1/n)
#
# ./{T<:Real}(n::T, s::DSisoRational)      = /(n, s)
# ./{T<:Real}(s::DSisoRational, n::T)      = /(s, n)
# ./(s1::DSisoRational, s2::DSisoRational) = /(s1, s2)
#
# function ==(s1::DSisoRational, s2::DSisoRational)
#   s1.Ts == s2.Ts && s1.num == s2.num &&
#     (s1.den == s2.den || s1.num == zero(s1.num))
#     # TODO scaling of num and den
# end
#
# !=(s1::DSisoRational, s2::DSisoRational) = !(s1 == s2)
#
# function isapprox(s1::DSisoRational, s2::DSisoRational,
#     rtol::Real=sqrt(eps()), atol::Real=0)
#   sdiff = s2-s1
#   return norm(sdiff.num) < rtol
# end
