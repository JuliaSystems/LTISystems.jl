# Value type to differentiate between SISO and MIMO systems
immutable Siso{T}
end

# Value type to differentiate between continuous- and discrete-time systems
immutable Continuous{T}
end

abstract LtiSystem{T,S}

# Printing functions
summary(s::LtiSystem{Siso{true},Continuous{true}})    =
  string("siso(nx=", numstates(s)::Int, ")")
summary(s::LtiSystem{Siso{true},Continuous{false}})   =
  string("siso(nx=", numstates(s)::Int, ",Ts=", samplingtime(s)::Float64, ")")
summary(s::LtiSystem{Siso{false},Continuous{true}})   =
  string("mimo(nx=", numstates(s)::Int, ",ny=", numoutputs(s)::Int, ",nu=",
    numinputs(s)::Int, ")")
summary(s::LtiSystem{Siso{false},Continuous{false}})  =
  string("mimo(nx=", numstates(s)::Int, ",ny=", numoutputs(s)::Int, ",nu=",
    numinputs(s)::Int, ",Ts=", samplingtime(s)::Float64, ")")

showcompact(io::IO, s::LtiSystem)     = print(io, summary(s))
show(io::IO, s::LtiSystem)            = print(io, summary(s))
showall(io::IO, s::LtiSystem)         = print(io, summary(s))

# Sampling time related functions
iscontinuous{T,S}(::LtiSystem{T,Continuous{S}})     = S::Bool
isdiscrete{T,S}(::LtiSystem{T,Continuous{S}})       = !S::Bool
samplingtime{T}(::LtiSystem{T,Continuous{true}})    = zero(Float64)
samplingtime{T}(s::LtiSystem{T,Continuous{false}})  = s.Ts::Float64

# Legacy SISO system
# abstract SisoSystem{T<:Real} <: LtiSystem
#
# # I/O mapping
# numstates(s::SisoSystem)              = degree(numpoly(s)) + degree(denpoly(s))
# numinputs(s::SisoSystem)              = 1
# numoutputs(s::SisoSystem)             = 1
#
# # Dimension information
# ndims(s::SisoSystem)                  = 1
# size(s::SisoSystem)                   = 1
# size(s::SisoSystem, dim::Int)         = 1
# size(s::SisoSystem, dims::Int...)     = map(x->size(s, x), dims)
#
# # Iteration interface
# start(::SisoSystem)                   = 1
# next(s::SisoSystem, state::Int)       = (s[state], state+1)
# done(s::SisoSystem, state::Int)       = state > 1
# eltype{T<:SisoSystem}(::Type{T})      = T
# length(s::SisoSystem)                 = 1
# eachindex(s::SisoSystem)              = 1
# endof(s::SisoSystem)                  = 1
#
# getindex(s::SisoSystem, idx::Int)     = (idx == 1) ? s :
#   (warn("s[idx]: Trying to access idx != 1"); throw(BoundsError()))
# getindex(s::SisoSystem, ::Colon)      = s
#
# # Printing functions
# summary(s::SisoSystem) = string("siso(nx=", numstates(s), (isdiscrete(s) ?
#   string(",Ts=", samplingtime(s)) : ""), ")")
#
# showcompact(io::IO, s::SisoSystem)    = print(io, summary(s))
# show(io::IO, s::SisoSystem)           = print(io, summary(s))
# showall(io::IO, s::SisoSystem)        = print(io, summary(s))

# Legacy MIMO system
# abstract MimoSystem{T<:SisoSystem} <: LtiSystem
#
# # common constructor
#
# function mimo{M<:AbstractArray}(m::M)
#   @assert eltype(m) <: SisoSystem
#   if length(m) > 0
#     if all(map(isdiscrete,m))
#       DMimo(m)
#     elseif all(map(iscontinuous,m))
#       CMimo(m)
#     else
#       warn("Do not support mixed continuous and discrete MIMO")
#       throw(DomainError)
#     end
#   else
#     warn("Do not support MIMO construction of empty collection")
#     throw(DomainError)
#   end
# end
#
# # I/O mapping
# numstates(s::MimoSystem)            = map(numstates, getmatrix(s)::AbstractArray)
# numinputs(s::MimoSystem)            = size(getmatrix(s)::AbstractArray, 2)
# numoutputs(s::MimoSystem)           = size(getmatrix(s)::AbstractArray, 1)
#
# # Dimension information
# ndims(s::MimoSystem)                = ndims(getmatrix(s)::AbstractArray)
# size(s::MimoSystem)                 = size(getmatrix(s)::AbstractArray)
# size(s::MimoSystem, dim::Int)       = size(getmatrix(s)::AbstractArray, dim)
# size(s::MimoSystem, dims::Int...)   = map(x->size(s, x), dims)
#
# # Iteration interface
# start(::MimoSystem)                 = 1
# next(s::MimoSystem, state::Int)     = (s[state], state+1)
# done(s::MimoSystem, state::Int)     = state > length(s)
# eltype{T<:SisoSystem}(::Type{MimoSystem{T}})  = T
# length(s::MimoSystem)               = length(getmatrix(s)::AbstractArray)
# eachindex(s::MimoSystem)            = eachindex(getmatrix(s)::AbstractArray)
# endof(s::MimoSystem)                = endof(getmatrix(s)::AbstractArray)
#
# getindex{T<:SisoSystem}(s::MimoSystem{T}, idx::Int)           =
#   getindex(getmatrix(s)::AbstractArray, idx)
# getindex{T<:SisoSystem}(s::MimoSystem{T}, row::Int, col::Int) =
#   getindex(getmatrix(s)::AbstractArray, row, col)
# getindex{T<:SisoSystem}(s::MimoSystem{T}, ::Colon, ::Colon)   = s
# getindex{T<:SisoSystem}(s::MimoSystem{T}, ::Colon, cols)      =
#   mimo(getindex(getmatrix(s)::AbstractArray, :, cols))
# getindex{T<:SisoSystem}(s::MimoSystem{T}, rows::Int, ::Colon) =
#   mimo(getindex(getmatrix(s)::AbstractArray, rows, :))
# getindex{T<:SisoSystem}(s::MimoSystem{T}, ::Colon)          = mimo(s.m[:])
#
# # Common type interface
#
# zeros(s::MimoSystem)   = map(poles, getmatrix(s)::AbstractArray)
# poles(s::MimoSystem)   = map(poles, getmatrix(s)::AbstractArray)
# numvec(s::MimoSystem)  = map(numvec, getmatrix(s)::AbstractArray)
# denvec(s::MimoSystem)  = map(denvec, getmatrix(s)::AbstractArray)
# numpoly(s::MimoSystem) = map(numpoly, getmatrix(s)::AbstractArray)
# denpoly(s::MimoSystem) = map(denpoly, getmatrix(s)::AbstractArray)
# zpkdata(s::MimoSystem) = map(zpkdata, getmatrix(s)::AbstractArray)
#
# # Printing functions
# summary(s::MimoSystem) = string("mimo(nx=", numstates(s), ",ny=", numoutputs(s),
#   ",nu=", numinputs(s), (isdiscrete(s) ? string(",Ts=", samplingtime(s)) : ""),
#   ")")
#
# showcompact(io::IO, s::MimoSystem)  = print(io, summary(s))
# show(io::IO, s::MimoSystem)         = print(io, summary(s))
# showall(io::IO, s::MimoSystem)      = print(io, summary(s))

# # Comparison of `LtiSystem`s
# function ==(s1::LtiSystem, s2::LtiSystem)
#   # TODO: Implement
# end
#
# !=(s1::LtiSystem, s2::LtiSystem) = !(s1 == s2)
#
# function isapprox(s1::LtiSystem, s2::LtiSystem)
#   # TODO: Implement
# end
