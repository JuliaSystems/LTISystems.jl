abstract CSisoTf <: SisoTf

# Sampling information
isdiscrete(s::CSisoTf)          = false
isdiscrete(::Type{CSisoTf})     = false
samplingtime(s::CSisoTf)        = zero(T)

# Printing functions
summary(s::CSisoTf)             = string("tf(nx=", numstates(s), ")")

showcompact(io::IO, s::CSisoTf) = print(io, summary(s))
show(io::IO, s::CSisoTf)        = print(io, summary(s))
showall(io::IO, s::CSisoTf)     = print(io, summary(s))
