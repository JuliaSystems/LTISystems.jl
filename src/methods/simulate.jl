# TODO: For now, keep SimType{T} in `dense` format. Later on, think of relaxing
#       this thing, *i.e.*, SimType{T<:Real,V1<:AbstractVector{T},etc.}.

type SimType{T<:Real} <: DEDataVector{T}
  x::Vector{T}
  y::Vector{T}
  u::Vector{T}
end

function (::Type{SimType})(x::Vector, y::Vector, u::Vector)
  T = promote_type(eltype(x), eltype(y), eltype(u))
  SimType(convert(Vector{T}, x), convert(Vector{T}, y), convert(Vector{T}, u))
end

function simulate{T}(sys::LtiSystem{Val{T},Val{:cont}}, tspan;
  input = (t,x)->zeros(numinputs(sys)), alg::AbstractODEAlgorithm = Tsit5(),
  initial::AbstractVector = zeros(numstates(sys)), tstops = Float64[], kwargs...)

  f     = (t,x,dx)->sys(t,x,dx,input)
  simvar= SimType(initial, zeros(numoutputs(sys)), zeros(numinputs(sys)))
  f(tspan[1], simvar, zeros(initial))

  tstops= [tstops..., discontinuities(input, tspan)...]
  prob  = ODEProblem(f, simvar, tspan)
  sln   = OrdinaryDiffEq.solve(prob, alg; tstops = tstops, kwargs...)

  tvals = sln.t
  xvals = Matrix{eltype(sln[1].x)}(length(tvals), numstates(sys))
  yvals = Matrix{eltype(sln[1].y)}(length(tvals), numoutputs(sys))
  uvals = Matrix{eltype(sln[1].u)}(length(tvals), numinputs(sys))
  for idx in 1:length(tvals)
    xvals[idx,:] = sln[idx]
    yvals[idx,:] = sln[idx].y
    uvals[idx,:] = sln[idx].u
  end

  TimeResponse(tvals,xvals,yvals,uvals)
end

function simulate{T}(sys::LtiSystem{Val{T},Val{:disc}}, tspan;
  input = (t,x)->zeros(numinputs(sys)),
  initial::AbstractVector = zeros(numstates(sys)), kwargs...)

  f      = (t,x,dx)->sys(t,x,dx,input)
  simvar = SimType(initial, zeros(numoutputs(sys)), zeros(numinputs(sys)))
  temp   = zeros(simvar.x)
  f(0.0, simvar, temp)  # ensure first step is correct
  simvar.x = temp       # alternative is adding extra time step, see below

  dt    = samplingtime(sys)
  prob  = DiscreteProblem(f, simvar, tspan)
  sln   = OrdinaryDiffEq.solve(prob, FunctionMap(scale_by_time=false);
    kwargs..., dt = dt)

  tvals = sln.t
  xvals = Matrix{eltype(sln[1].x)}(length(tvals), numstates(sys))
  yvals = Matrix{eltype(sln[1].y)}(length(tvals), numoutputs(sys))
  uvals = Matrix{eltype(sln[1].u)}(length(tvals), numinputs(sys))
  for idx in 1:length(tvals)
    xvals[idx,:] = sln[idx]
    yvals[idx,:] = sln[idx].y
    uvals[idx,:] = sln[idx].u
  end

  TimeResponse(tvals,xvals,yvals,uvals,Val{:disc})
end
# to ensure that first step is correct
# simvar is initialized above with correct first sample of
# input and dx.
# The alternative is to manually add extra time step and use
# (tspan[1]-dt, tspan[2]) in simulation.
# hacky fix of initial step with diffEq.

simulate(input::Function, sys::LtiSystem, tspan; kwargs...) =
  simulate(sys, tspan; input = input, kwargs...)

step(sys::LtiSystem, tspan; kwargs...) = simulate(sys, tspan; kwargs...,
  input = Step(zeros(numinputs(sys)), ones(numinputs(sys))))
ramp(sys::LtiSystem, tspan; kwargs...) = simulate(sys, tspan; kwargs...,
  input = Ramp(zeros(numinputs(sys)), ones(numinputs(sys))))
