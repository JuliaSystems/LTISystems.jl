# definition of discrete mimo transfer function

DSiso = Union{DSisoTf}

immutable DMimo{T<:DSiso, M<:AbstractArray{T}} <: MimoSystem{T}
  m::M
  ny::Int
  nu::Int
  function call{M<:AbstractArray}(::Type{DMimo}, m::M)
    @assert 0 < ndims(m)
    @assert eltype(m) <: DSiso

    T = eltype(M)
    new{T,M}(m, size(m, 1), size(m, 2))
  end
end

# creation of discrete mimo transfer functions

function tf{M1<:AbstractArray, M2<:AbstractArray, T1<:Real}(num::M1,
  den::M2, Ts::T1)
  @assert eltype(eltype(num)) <: Real
  @assert eltype(eltype(den)) <: Real

  if size(num, 1, 2) != size(den, 1, 2)
    warn("num and den dimensions must match")
    throw(DomainError())
  end
  m = similar(num, DSisoRational)
  for idx in eachindex(m)
    m[idx] = tf(num[idx], den[idx], Ts)
  end
  DMimo(m)
end

function tf{M<:AbstractArray, T1<:Real}(gain::M, Ts::T1)
  @assert eltype(eltype(gain)) <: Real

  m = similar(gain, DSisoRational{eltype(gain)})
  for idx in eachindex(m)
    m[idx] = tf(gain[idx], Ts)
  end
  DMimo(m)
end

function zpk{M1<:AbstractArray, M2<:AbstractArray, M3<:AbstractArray, T1<:Real}(
  z::M1, p::M2, k::M3, Ts::T1)
  @assert eltype(eltype(z)) <: Number
  @assert eltype(eltype(p)) <: Number
  @assert eltype(k)         <: Real

  # Validate input and output dimensions match
  ny, nu = size(z, 1, 2)
  if (ny, nu) != size(p, 1, 2) || (ny, nu) != size(k, 1, 2)
    warn("z, p and k dimensions must match")
    throw(DomainError())
  end
  m = similar(z, DSisoZpk)
  for idx in eachindex(m)
    m[idx] = zpk(z[idx], p[idx], k[idx], Float64(Ts))
  end
  DMimo(m)
end

function zpk{M1<:AbstractArray, T1<:Real}(gain::M1, Ts::T1)
  @assert eltype(gain) <: Real

  m = similar(gain, DSisoZpk)
  for idx in eachindex(m)
    m[idx] = zpk(gain[idx], Ts)
  end
  DMimo(m)
end

# conversion and promotion

promote_rule{T<:Real,T1,T2}(::Type{DMimo{T1,T2}}, ::Type{T}) = DMimo
convert{T<:Real,T1,T2}(::Type{DMimo{T1,T2}}, x::T) =
  zpk(sparse([1],[1],[x]), zero(Float64))

promote_rule{T<:Real}(::Type{DMimo}, ::Type{T}) = DMimo
convert{T<:Real}(::Type{DMimo}, x::T) = zpk(sparse([1],[1],[x]), zero(Float64))

promote_rule{M1<:AbstractArray,T1,T2}(::Type{DMimo{T1,T2}}, ::Type{M1}) = DMimo
promote_rule{M1<:AbstractArray}(::Type{DMimo}, ::Type{M1}) = DMimo
convert{M1<:AbstractArray}(::Type{DMimo}, x::M1) = zpk(x, zero(Float64))

# overloading identities

one(s::DMimo)                     = zpk(sparse([1],[1],[Int8(1)]), zero(Float64))
one(::Type{DMimo})                = zpk(sparse([1],[1],[Int8(1)]), zero(Float64))
one{T}(::Type{DMimo{T}})          = zpk(sparse([1],[1],[Int8(1)]), zero(Float64))
one{T1,T2}(::Type{DMimo{T1,T2}})  = zpk(sparse([1],[1],[Int8(1)]), zero(Float64))
zero(s::DMimo)                    = zpk(sparse([1],[1],[Int8(0)]), zero(Float64))
zero(::Type{DMimo})               = zpk(sparse([1],[1],[Int8(0)]), zero(Float64))
zero{T}(::Type{DMimo{T}})         = zpk(sparse([1],[1],[Int8(0)]), zero(Float64))
zero{T1,T2}(::Type{DMimo{T1,T2}}) = zpk(sparse([1],[1],[Int8(0)]), zero(Float64))

# interface implementation

samplingtime(s::DMimo)    = map(samplingtime, s.m)
getmatrix(s::DMimo)       = s.m
isdiscrete(s::DMimo)      = true
isdiscrete(::Type{DMimo}) = true
transpose(s::DMimo)       = mimo(deepcopy(s.m).')

# overload printing functions

showcompact(io::IO, s::DMimo) = print(io, summary(s))

function show(io::IO, s::DMimo)
  println(io, "Discrete time transfer function model")
  println(io, "\ty = Gu")
end

function showall(io::IO, s::DMimo)
  show(io, s)
  println(io, "")
  for subs in s
    printtransferfunction(io::IO, subs)
    println(io, "")
  end
end

# overload mathematical operations

function +(s1::DMimo, s2::DMimo)
  if size(s1.m) != size(s2.m)
    warn("Systems have different shapes")
    throw(DomainError)
  end
  m = s1.m + s2.m
  mimo(m)
end

+{T<:Real}(s::DMimo, n::T)          = mimo(s.m+n)
+{T<:Real}(n::T, s::DMimo)          = +(s, n)

function +{T<:Real}(s::DMimo, n::Matrix{T})
  if size(s.m) != size(n)
    warn("Systems have different shapes")
    throw(DomainError)
  end
  mimo(s.m + n)
end
+{T<:Real}(n::Matrix{T}, s::DMimo)  = +(s, n)

.+{T<:Real}(n::T, s::DMimo)         = +(n, s)
.+{T<:Real}(s::DMimo, n::T)         = +(s, n)

.+{T<:Real}(n::Matrix{T}, s::DMimo) = +(n, s)
.+{T<:Real}(s::DMimo, n::Matrix{T}) = +(s, n)
.+(s1::DMimo, s2::DMimo)            = +(s1, s2)

-(s::DMimo)                         = mimo(-s.m)

-{T<:Real}(s::DMimo, n::Matrix{T})  = +(s, -n)
-{T<:Real}(n::Matrix{T}, s::DMimo)  = +(-n, s)
-(s1::DMimo,s2::DMimo)              = +(s1,-s2)
-{T<:Real}(s::DMimo, n::T)          = +(s,-n)
-{T<:Real}(n::T, s::DMimo)          = +(-n, s)

.-{T<:Real}(n::T, s::DMimo)         = +(-n, s)
.-{T<:Real}(s::DMimo, n::T)         = +(s,-n)

.-{T<:Real}(n::Matrix{T}, s::DMimo) = -(n, s)
.-{T<:Real}(s::DMimo, n::Matrix{T}) = -(s, n)

function *(s1::DMimo, s2::DMimo)
  if size(s1.m, 2) != size(s2.m, 1)
    warn("s1*s2: s1 must have same number of inputs as s2 has outputs")
    throw(DomainError())
  end
  mimo(deepcopy(s1.m)*deepcopy(s2.m))
end

*{T<:Real}(s::DMimo, n::T)          = mimo(s.m*n)
*{T<:Real}(n::T, s::DMimo)          = *(s, n)

*{T<:Real}(s::DMimo, n::Matrix{T})  = mimo(s.m*n)
*{T<:Real}(n::Matrix{T}, s::DMimo)  = mimo(n*s.m)

.*{T<:Real}(n::T, s::DMimo)         = *(n, s)
.*{T<:Real}(s::DMimo, n::T)         = *(s, n)
.*{T<:Real}(n::Matrix{T}, s::DMimo) = mimo(n.*s.m)
.*{T<:Real}(s::DMimo, n::Matrix{T}) = mimo(s.m.*n)
.*(s1::DMimo, s2::DMimo)            = mimo(s1.m.*s2.m)

function /{T<:Real}(n::T, s::DMimo)
  warn("MIMO transfer function inversion is not implemented yet")
  throw(DomainError())
end
/(s1::DMimo, s2::DMimo)             = s1*(1/s2)
/{T<:Real}(s::DMimo, n::T)          = s*(1/n)
function /{T<:Real}(n::Matrix{T}, s::DMimo)
  warn("MIMO transfer function inversion is not implemented yet")
  throw(DomainError())
end
/{T<:Real}(s::DMimo, n::Matrix{T})  = mimo(s.m/n)

./{T<:Real}(s::DMimo, n::T)         = s/n
./{T<:Real}(n::T, s::DMimo)         = n/s
./{T<:Real}(s::DMimo, n::Matrix{T}) = s/n
./{T<:Real}(n::Matrix{T}, s::DMimo) = n/s

==(s1::DMimo, s2::DMimo)            = s1.m == s2.m

!=(s1::DMimo, s2::DMimo)            = !(s1.m == s2.m)

function isapprox(s1::DMimo, s2::DMimo)
  if size(s1.m) != size(s2.m)
    return false
  end
  a = 0
  for i in eachindex(s1.m)
    a += s1.m[i] ≈ s2.m[i]
  end
  a < length(s1.m) ? false : true
end
