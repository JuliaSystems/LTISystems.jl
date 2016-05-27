abstract SisoSs{T<:Real} <: SisoSystem{T}

# Printing functions
summary(s::SisoSs) = string("ss(nx=", numstates(s), (isdiscrete(s) ?
  string(",Ts=", samplingtime(s)) : ""), ")")

showcompact(io::IO, s::SisoSs)  = print(io, summary(s))
show(io::IO, s::SisoSs)         = print(io, summary(s))
showall(io::IO, s::SisoSs)      = print(io, summary(s))
