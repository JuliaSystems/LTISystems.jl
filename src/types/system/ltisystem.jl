# LTI system abstraction
@compat abstract type LtiSystem{T,S} end

# Interfaces
issiso(::LtiSystem{Val{:siso}}) = true
ismimo(::LtiSystem{Val{:mimo}}) = false

# Sampling-time
iscontinuous{T}(::LtiSystem{Val{T},Val{:cont}}) = true
iscontinuous{T}(::LtiSystem{Val{T},Val{:disc}}) = false
isdiscrete(s::LtiSystem)                        = !iscontinuous(s)

function samplingtime(s::LtiSystem)
  warn("samplingtime(s) not implemented for s::$(typeof(s))")
  throw(MethodError(samplingtime, Tuple{typeof(s)}))
end

# State, input, output information
function numstates(s::LtiSystem)
  warn("numstates(s) not implemented for s::$(typeof(s))")
  throw(MethodError(numstates, Tuple{typeof(s)}))
end

function numinputs(s::LtiSystem)
  warn("numinputs(s) not implemented for s::$(typeof(s))")
  throw(MethodError(numinputs, Tuple{typeof(s)}))
end

function numoutputs(s::LtiSystem)
  warn("numoutputs(s) not implemented for s::$(typeof(s))")
  throw(MethodError(numoutputs, Tuple{typeof(s)}))
end

# Iteration interface (meaningful in MIMO)
function start(s::LtiSystem{Val{:mimo}})
  warn("start(s) not implemented for s::$(typeof(s))")
  throw(MethodError(start, Tuple{typeof(s)}))
end

function next(s::LtiSystem{Val{:mimo}}, state)
  warn("next(s, state) not implemented for s::$(typeof(s))")
  throw(MethodError(next, Tuple{typeof(s), typeof(state)}))
end

function done(s::LtiSystem{Val{:mimo}}, state)
  warn("done(s, state) not implemented for s::$(typeof(s))")
  throw(MethodError(done, Tuple{typeof(s), typeof(state)}))
end

function eltype(s::Type{LtiSystem{Val{:mimo}}})
  warn("eltype(s) not implemented for s::$(typeof(s))")
  throw(MethodError(eltype, Tuple{typeof(s)}))
end

function length(s::LtiSystem{Val{:mimo}})
  warn("length(s) not implemented for s::$(typeof(s))")
  throw(MethodError(length, Tuple{typeof(s)}))
end

function size(s::LtiSystem{Val{:mimo}}, d)
  warn("size(s, d) not implemented for s::$(typeof(s))")
  throw(MethodError(size, Tuple{typeof(s), typeof(d)}))
end

# Indexing
function getindex(s::LtiSystem{Val{:mimo}}, i)
  warn("getindex(s, i) not implemented for s::$(typeof(s))")
  throw(MethodError(getindex, Tuple{typeof(s), typeof(i)}))
end

function setindex!(s::LtiSystem{Val{:mimo}}, v, i)
  warn("setindex!(s, v, i) not implemented for s::$(typeof(s))")
  throw(MethodError(setindex!, Tuple{typeof(s), typeof(v), typeof(i)}))
end

function endof(s::LtiSystem{Val{:mimo}})
  warn("endof(s) not implemented for s::$(typeof(s))")
  throw(MethodError(endof, Tuple{typeof(s)}))
end
