abstract LtiSystem

# Printing functions
summary(s::LtiSystem) = string("lti(nx=", numstates(s), ",ny=", numoutputs(s),
  ",nu=", numinputs(s), (isdiscrete(s) ? string(",Ts=", samplingtime(s)) : ""),
  ")")

showcompact(io::IO, s::LtiSystem)     = print(io, summary(s))
show(io::IO, s::LtiSystem)            = print(io, summary(s))
showall(io::IO, s::LtiSystem)         = print(io, summary(s))

iscontinuous(s::LtiSystem)            = !isdiscrete(s)
iscontinuous{T<:LtiSystem}(::Type{T}) = !isdiscrete(T)

# LtiSystem math
+(s1::LtiSystem, s2::LtiSystem)  = +(tf(s1), tf(s2))
+{T<:Real}(s::LtiSystem, n::T)   = +(tf(s), n)
+{T<:Real}(n::T, s::LtiSystem)   = +(n, tf(s))
.+{T<:Real}(s::LtiSystem, n::T)  = +(tf(s), n)
.+{T<:Real}(n::T, s::LtiSystem)  = +(n, tf(s))
.+(s1::LtiSystem, s2::LtiSystem) = +(tf(s1), tf(s2))

-(s::LtiSystem)                  = -tf(s)
-(s1::LtiSystem, s2::LtiSystem)  = +(tf(s1), -tf(s2))
-{T<:Real}(s::LtiSystem, n::T)   = +(tf(s), -n)
-{T<:Real}(n::T, s::LtiSystem)   = +(n, -tf(s))
.-{T<:Real}(s::LtiSystem, n::T)  = -(tf(s), n)
.-{T<:Real}(n::T, s::LtiSystem)  = -(n, tf(s))
.-(s1::LtiSystem, s2::LtiSystem) = -(tf(s1), tf(s2))

*(s1::LtiSystem, s2::LtiSystem)  = *(tf(s1), tf(s2))
*{T<:Real}(s::LtiSystem, n::T)   = *(tf(s), n)
*{T<:Real}(n::T, s::LtiSystem)   = *(n, tf(s))
.*{T<:Real}(s::LtiSystem, n::T)  = *(tf(s), n)
.*{T<:Real}(n::T, s::LtiSystem)  = *(n, tf(s))
.*(s1::LtiSystem, s2::LtiSystem) = *(tf(s1), tf(s2))

/(s1::LtiSystem, s2::LtiSystem)  = *(tf(s1),1/tf(s2))
/{T<:Real}(n::T, s::LtiSystem)   = /(n, tf(s))
/{T<:Real}(s::LtiSystem, n::T)   = /(tf(s), n)
./{T<:Real}(n::T, s::LtiSystem)  = /(n, tf(s))
./{T<:Real}(s::LtiSystem, n::T)  = /(tf(s), n)
./(s1::LtiSystem, s2::LtiSystem) = /(tf(s1), tf(s2))
