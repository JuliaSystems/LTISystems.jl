# LTI system abstraction
@compat abstract type LtiSystem{T,S} end

# Interfaces
issiso(::LtiSystem{Val{:siso}}) = true
ismimo(::LtiSystem{Val{:mimo}}) = false

# conversion between 1×1 mimo and siso
function siso(s::LtiSystem{Val{:mimo}})
  warn("siso(s) not implemented for s::$(typeof(s))")
  throw(MethodError(siso, Tuple{typeof(s)}))
end
function mimo(s::LtiSystem{Val{:siso}})
  warn("mimo(s) not implemented for s::$(typeof(s))")
  throw(MethodError(mimo, Tuple{typeof(s)}))
end

# Sampling-time
iscontinuous(::LtiSystem{Val{T},Val{:cont}}) where {T} = true
iscontinuous(::LtiSystem{Val{T},Val{:disc}}) where {T} = false
isdiscrete(s::LtiSystem)                               = !iscontinuous(s)
samplingtime(::LtiSystem{Val{T},Val{:disc}}) where {T} = zero(Float64)
