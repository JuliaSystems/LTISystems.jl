# TODO: For now, keep SimType{T} in `dense` format. Later on, think of relaxing
#       this thing, *i.e.*, SimType{T<:Real,V1<:AbstractVector{T},etc.}.

type SimType{T<:Real} <: DEDataArray{T}
  x::Vector{T}
  y::Vector{T}
  u::Vector{T}
end

function (::Type{SimType})(x::Vector, y::Vector, u::Vector)
  T = promote_type(eltype(x), eltype(y), eltype(u))
  SimType(convert(Vector{T}, x), convert(Vector{T}, y), convert(Vector{T}, u))
end

function simulate{T}(sys::LtiSystem{Val{T},Val{:cont}}, tspan;
  input = (t,x)->zeros(numinputs(sys)), alg::OrdinaryDiffEqAlgorithm = Tsit5(),
  initial::AbstractVector = zeros(numstates(sys)), tstops = Float64[], kwargs...)

  f     = (t,x,dx)->sys(t,x,dx,input)
  simvar= SimType(initial, zeros(numoutputs(sys)), zeros(numinputs(sys)))

  tstops= [tstops..., discontinuities(input, tspan)...]
  prob  = ODEProblem(f, simvar, tspan)
  sln   = OrdinaryDiffEq.solve(prob, alg; tstops = tstops, kwargs...)

  tvals = sln.t
  xvals = Matrix{eltype(tvals)}(length(tvals), numstates(sys))
  yvals = Matrix{eltype(tvals)}(length(tvals), numoutputs(sys))
  uvals = Matrix{eltype(tvals)}(length(tvals), numinputs(sys))
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

  f     = (t,x,dx)->sys(t,x,dx,input)
  simvar= SimType(initial, zeros(numoutputs(sys)), zeros(numinputs(sys)))

  prob  = DiscreteProblem(f, simvar, tspan)
  dt    = samplingtime(sys)
  sln   = OrdinaryDiffEq.solve(prob, FunctionMap(scale_by_time=false);
    kwargs..., dt = samplingtime(sys))

  tvals = sln.t
  xvals = Matrix{eltype(tvals)}(length(tvals), numstates(sys))
  yvals = Matrix{eltype(tvals)}(length(tvals), numoutputs(sys))
  uvals = Matrix{eltype(tvals)}(length(tvals), numinputs(sys))
  for idx in 1:length(tvals)
    xvals[idx,:] = sln[idx]
    yvals[idx,:] = sln[idx].y
    uvals[idx,:] = sln[idx].u
  end

  TimeResponse(tvals,xvals,yvals,uvals,Val{:disc})
end

simulate(input::Function, sys::LtiSystem, tspan; kwargs...) =
  simulate(sys, tspan; input = input, kwargs...)

step(sys::LtiSystem, tspan; kwargs...) = simulate(sys, tspan, kwargs...,
  input = Step(zeros(numinputs(sys)), ones(numinputs(sys))))
ramp(sys::LtiSystem, tspan; kwargs...) = simulate(sys, tspan, kwargs...,
  input = Ramp(zeros(numinputs(sys)), ones(numinputs(sys))))
