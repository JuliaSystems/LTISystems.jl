abstract MimoSystem{T<:SisoSystem} <: LtiSystem

# common constructor

function mimo{M<:AbstractArray}(m::M)
  @assert eltype(m) <: SisoSystem
  if length(m) > 0
    if all(map(isdiscrete,m))
      DMimo(m)
    elseif all(map(iscontinuous,m))
      CMimo(m)
    else
      warn("Do not support mixed continuous and discrete MIMO")
      throw(DomainError)
    end
  else
    warn("Do not support MIMO construction of empty collection")
    throw(DomainError)
  end
end

# I/O mapping
numstates(s::MimoSystem)            = map(numstates, getmatrix(s)::AbstractArray)
numinputs(s::MimoSystem)            = size(getmatrix(s)::AbstractArray, 2)
numoutputs(s::MimoSystem)           = size(getmatrix(s)::AbstractArray, 1)

# Dimension information
ndims(s::MimoSystem)                = ndims(getmatrix(s)::AbstractArray)
size(s::MimoSystem)                 = size(getmatrix(s)::AbstractArray)
size(s::MimoSystem, dim::Int)       = size(getmatrix(s)::AbstractArray, dim)
size(s::MimoSystem, dims::Int...)   = map(x->size(s, x), dims)

# Iteration interface
start(::MimoSystem)                 = 1
next(s::MimoSystem, state::Int)     = (s[state], state+1)
done(s::MimoSystem, state::Int)     = state > length(s)
eltype{T<:SisoSystem}(::Type{MimoSystem{T}})  = T
length(s::MimoSystem)               = length(getmatrix(s)::AbstractArray)
eachindex(s::MimoSystem)            = eachindex(getmatrix(s)::AbstractArray)
endof(s::MimoSystem)                = endof(getmatrix(s)::AbstractArray)

getindex{T<:SisoSystem}(s::MimoSystem{T}, idx::Int)           =
  getindex(getmatrix(s)::AbstractArray, idx)
getindex{T<:SisoSystem}(s::MimoSystem{T}, row::Int, col::Int) =
  getindex(getmatrix(s)::AbstractArray, row, col)
getindex{T<:SisoSystem}(s::MimoSystem{T}, ::Colon, ::Colon)   = s
getindex{T<:SisoSystem}(s::MimoSystem{T}, ::Colon, cols)      =
  mimo(getindex(getmatrix(s)::AbstractArray, :, cols))
getindex{T<:SisoSystem}(s::MimoSystem{T}, rows::Int, ::Colon) =
  mimo(getindex(getmatrix(s)::AbstractArray, rows, :))
getindex{T<:SisoSystem}(s::MimoSystem{T}, ::Colon)          = mimo(s.m[:])

# Common type interface

zeros(s::MimoSystem)   = map(poles, getmatrix(s)::AbstractArray)
poles(s::MimoSystem)   = map(poles, getmatrix(s)::AbstractArray)
numvec(s::MimoSystem)  = map(numvec, getmatrix(s)::AbstractArray)
denvec(s::MimoSystem)  = map(denvec, getmatrix(s)::AbstractArray)
numpoly(s::MimoSystem) = map(numpoly, getmatrix(s)::AbstractArray)
denpoly(s::MimoSystem) = map(denpoly, getmatrix(s)::AbstractArray)
zpkdata(s::MimoSystem) = map(zpkdata, getmatrix(s)::AbstractArray)

# Printing functions
summary(s::MimoSystem) = string("mimo(nx=", numstates(s), ",ny=", numoutputs(s),
  ",nu=", numinputs(s), (isdiscrete(s) ? string(",Ts=", samplingtime(s)) : ""),
  ")")

showcompact(io::IO, s::MimoSystem)  = print(io, summary(s))
show(io::IO, s::MimoSystem)         = print(io, summary(s))
showall(io::IO, s::MimoSystem)      = print(io, summary(s))
