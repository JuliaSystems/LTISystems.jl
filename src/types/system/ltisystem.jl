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
  warn("next(s, state) not implemented for (s::$(typeof(s)), state::$(typeof(state)))")
  throw(MethodError(next, Tuple{typeof(s),typeof(state)}))
end

function done(s::LtiSystem{Val{:mimo}}, state)
  warn("done(s, state) not implemented for (s::$(typeof(s)), state::$(typeof(state)))")
  throw(MethodError(done, Tuple{typeof(s),typeof(state)}))
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
  warn("size(s, d) not implemented for (s::$(typeof(s)), d::$(typeof(d)))")
  throw(MethodError(size, Tuple{typeof(s),typeof(d)}))
end

# Indexing
function getindex(s::LtiSystem{Val{:mimo}}, i)
  warn("getindex(s, i) not implemented for (s::$(typeof(s)), i::$(typeof(i)))")
  throw(MethodError(getindex, Tuple{typeof(s),typeof(i)}))
end

function setindex!(s::LtiSystem{Val{:mimo}}, v, i)
  warn("setindex!(s, v, i) not implemented for (s::$(typeof(s)), v::$(typeof(v)), i::$(typeof(i)))")
  throw(MethodError(setindex!, Tuple{typeof(s),typeof(v),typeof(i)}))
end

function endof(s::LtiSystem{Val{:mimo}})
  warn("endof(s) not implemented for s::$(typeof(s))")
  throw(MethodError(endof, Tuple{typeof(s)}))
end

# Simple analysis things
function poles(s::LtiSystem)
  warn("poles(s) not implemented for s::$(typeof(s))")
  throw(MethodError(poles, Tuple{typeof(s)}))
end

function zeros(s::LtiSystem)
  warn("zeros(s) not implemented for s::$(typeof(s))")
  throw(MethodError(zeros, Tuple{typeof(s)}))
end

function tzeros(s::LtiSystem)
  warn("tzeros(s) not implemented for s::$(typeof(s))")
  throw(MethodError(tzeros, Tuple{typeof(s)}))
end

function zpkdata(s::LtiSystem)
  warn("zpkdata(s) not implemented for s::$(typeof(s))")
  throw(MethodError(zpkdata, Tuple{typeof(s)}))
end

function numvec(s::LtiSystem)
  warn("numvec(s) not implemented for s::$(typeof(s))")
  throw(MethodError(numvec, Tuple{typeof(s)}))
end

function denvec(s::LtiSystem)
  warn("denvec(s) not implemented for s::$(typeof(s))")
  throw(MethodError(denvec, Tuple{typeof(s)}))
end

function numpoly(s::LtiSystem)
  warn("numpoly(s) not implemented for s::$(typeof(s))")
  throw(MethodError(numpoly, Tuple{typeof(s)}))
end

function denpoly(s::LtiSystem)
  warn("denpoly(s) not implemented for s::$(typeof(s))")
  throw(MethodError(denpoly, Tuple{typeof(s)}))
end

# Constructors among different types
function mfd(s::LtiSystem)
  warn("mfd(s) not implemented for s::$(typeof(s))")
  throw(MethodError(mfd, Tuple{typeof(s)}))
end

function ss(s::LtiSystem)
  warn("ss(s) not implemented for s::$(typeof(s))")
  throw(MethodError(ss, Tuple{typeof(s)}))
end

function tf(s::LtiSystem)
  warn("tf(s) not implemented for s::$(typeof(s))")
  throw(MethodError(tf, Tuple{typeof(s)}))
end

# Methods
function series{T1,T2,S}(s1::LtiSystem{Val{T1},Val{S}}, s2::LtiSystem{Val{T2},Val{S}})
  warn("series(s1, s2) not implemented for (s1::$(typeof(s1)), s2::$(typeof(s2)))")
  throw(MethodError(series, Tuple{typeof(s1),typeof(s2)}))
end

function series(s1::LtiSystem, s2::LtiSystem)
  warn("series(s1, s2) cannot be defined for (s1::$(typeof(s1)), s2::$(typeof(s2)))")
  throw(MethodError(series, Tuple{typeof(s1),typeof(s2)}))
end

function parallel{T1,T2,S}(s1::LtiSystem{Val{T1},Val{S}}, s2::LtiSystem{Val{T2},Val{S}})
  warn("parallel(s1, s2) not implemented for (s1::$(typeof(s1)), s2::$(typeof(s2)))")
  throw(MethodError(parallel, Tuple{typeof(s1),typeof(s2)}))
end

function parallel(s1::LtiSystem, s2::LtiSystem)
  warn("parallel(s1, s2) cannot be defined for (s1::$(typeof(s1)), s2::$(typeof(s2)))")
  throw(MethodError(parallel, Tuple{typeof(s1),typeof(s2)}))
end

function feedback{T1,T2,S}(s1::LtiSystem{Val{T1},Val{S}}, s2::LtiSystem{Val{T2},Val{S}},
  neg::Bool = true)
  warn("feedback(s1, s2, neg) not implemented for (s1::$(typeof(s1)), s2::$(typeof(s2)), neg::$(typeof(neg)))")
  throw(MethodError(feedback, Tuple{typeof(s1),typeof(s2),typeof(neg)}))
end

function feedback(s1::LtiSystem, s2::LtiSystem, neg::Bool = true)
  warn("feedback(s1, s2, neg) cannot be defined for (s1::$(typeof(s1)), s2::$(typeof(s2)), neg::$(typeof(neg)))")
  throw(MethodError(feedback, Tuple{typeof(s1),typeof(s2),typeof(neg)}))
end

function minreal(s::LtiSystem)
  warn("minreal(s) not implemented for s::$(typeof(s))")
  throw(MethodError(minreal, Tuple{typeof(s)}))
end

function freqresp(s::LtiSystem, ω)
  warn("freqresp(s, ω) not implemented for (s::$(typeof(s1)), ω::$(typeof(ω)))")
  throw(MethodError(freqresp, Tuple{typeof(s),typeof(ω)}))
end

function evalfr(s::LtiSystem, ω)
  warn("evalfr(s, ω) not implemented for (s::$(typeof(s1)), ω::$(typeof(ω)))")
  throw(MethodError(evalfr, Tuple{typeof(s),typeof(ω)}))
end
