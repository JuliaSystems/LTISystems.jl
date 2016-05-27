abstract SisoSystem{T<:Real} <: LtiSystem

# I/O mapping
numstates(s::SisoSystem)              = degree(numpoly(s)) + degree(denpoly(s))
numinputs(s::SisoSystem)              = 1
numoutputs(s::SisoSystem)             = 1

# Dimension information
ndims(s::SisoSystem)                  = 1
size(s::SisoSystem)                   = 1
size(s::SisoSystem, dim::Int)         = 1
size(s::SisoSystem, dims::Int...)     = map(x->size(s, x), dims)

# Iteration interface
start(::SisoSystem)                   = 1
next(s::SisoSystem, state::Int)       = (s[state], state+1)
done(s::SisoSystem, state::Int)       = state > 1
eltype{T<:SisoSystem}(::Type{T})      = T
length(s::SisoSystem)                 = 1
eachindex(s::SisoSystem)              = 1
endof(s::SisoSystem)                  = 1

getindex(s::SisoSystem, idx::Int)     = (idx == 1) ? s :
  (warn("s[idx]: Trying to access idx != 1"); throw(BoundsError()))
getindex(s::SisoSystem, ::Colon)      = s

# Printing functions
summary(s::SisoSystem) = string("siso(nx=", numstates(s), (isdiscrete(s) ?
  string(",Ts=", samplingtime(s)) : ""), ")")

showcompact(io::IO, s::SisoSystem)    = print(io, summary(s))
show(io::IO, s::SisoSystem)           = print(io, summary(s))
showall(io::IO, s::SisoSystem)        = print(io, summary(s))
