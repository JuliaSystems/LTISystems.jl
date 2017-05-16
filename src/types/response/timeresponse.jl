type TimeResponse{T,S} <: SystemResponse
  t::Vector{T}
  x::Matrix{T}
  y::Matrix{T}
  u::Matrix{T}
  tsloc::Int

  function (::Type{TimeResponse}){S}(t::Vector, x::Matrix, y::Matrix, u::Matrix,
    ::Type{Val{S}} = Val{:cont})
    (length(t) == size(x, 1) == size(y, 1) == size(u, 1)) ||
      DimensionMismatch("TimeResponse(t,x,y,u): dimension mistmatch")

    T = promote_type(eltype(t), eltype(x), eltype(y), eltype(u))
    if !(T <: Real)
      warn("TimeResponse(t,x,y,u): variables must be `Real` values")
      DomainError()
    end

    new{T,Val{S}}(convert(Vector{T}, t), convert(Matrix{T}, x), convert(Matrix{T}, y),
      convert(Matrix{T}, u), 0)
  end
end

start(res::TimeResponse)        = (res.tsloc = 1)
next(res::TimeResponse, state)  = (state += 1; res.tsloc = state; (res, state))
done(res::TimeResponse, state)  = (state ≥ length(res.t))

@recipe function f{T,S}(res::TimeResponse{T,Val{S}}; outputs = 1:size(res.y, 2),
  inputs = Int[], states = Int[])

  numplots  = 0
  curplot   = 0
  endidx    = res.tsloc == 0 ? length(res.t) : res.tsloc
  outputset = IntSet()
  inputset  = IntSet()
  stateset  = IntSet()

  for idx in outputs
    if 1 ≤ idx ≤ size(res.y, 2)
      push!(outputset, idx)
    else
      warn("$(idx) in outputs is out of range. Discarding...")
    end
  end

  for idx in inputs
    if 1 ≤ idx ≤ size(res.u, 2)
      push!(inputset, idx)
    else
      warn("$(idx) in inputs is out of range. Discarding...")
    end
  end

  for idx in states
    if 1 ≤ idx ≤ size(res.x, 2)
      push!(stateset, idx)
    else
      warn("$(idx) in states is out of range. Discarding...")
    end
  end

  numplots += !isempty(outputset)
  numplots += !isempty(inputset)
  numplots += !isempty(stateset)

  @assert numplots > 0 "plot(res::TimeResponse): nothing to plot"

  outputlab = reshape(["\$y_{$idx}\$" for idx in outputset], 1, length(outputset))
  inputlab  = reshape(["\$u_{$idx}\$" for idx in inputset], 1, length(inputset))
  statelab  = reshape(["\$x_{$idx}\$" for idx in stateset], 1, length(stateset))

  ylabs = String[]

  if !isempty(outputset)
    push!(ylabs, "Outputs")
  end
  if !isempty(inputset)
    push!(ylabs, "Inputs")
  end
  if !isempty(stateset)
    push!(ylabs, "States")
  end

  t_plotted = res.t

  ylims = Tuple[]
  y_plotted = res.y[:, [outputset...]]
  u_plotted = res.u[:, [inputset...]]
  x_plotted = res.x[:, [stateset...]]

  if !isempty(y_plotted)
    push!(ylims, (1.1*minimum(y_plotted), 1.1*maximum(y_plotted)))
  end
  if !isempty(u_plotted)
    push!(ylims, (1.1*minimum(u_plotted), 1.1*maximum(u_plotted)))
  end
  if !isempty(x_plotted)
    push!(ylims, (1.1*minimum(x_plotted), 1.1*maximum(x_plotted)))
  end

  layout  --> (numplots,1)
  link    --> :x
  xlims   --> (minimum(t_plotted), maximum(t_plotted))
  ylims   --> reshape(ylims, 1, length(ylims))

  xguide  --> "t (sec)"
  yguide  --> reshape(ylabs, 1, length(ylabs))
  label   --> [outputlab inputlab statelab]

  if !isempty(outputset)
    curplot += 1
    @series begin
      subplot     :=  curplot
      seriestype  :=  ifelse(S == :cont, :path, :steppre)
      t_plotted[1:endidx], y_plotted[1:endidx,:]
    end
  end

  if !isempty(inputset)
    curplot += 1
    @series begin
      subplot     :=  curplot
      seriestype  :=  ifelse(S == :cont, :path, :steppre)
      t_plotted[1:endidx], u_plotted[1:endidx,:]
    end
  end

  if !isempty(stateset)
    curplot += 1
    @series begin
      subplot     :=  curplot
      seriestype  :=  ifelse(S == :cont, :path, :steppre)
      t_plotted[1:endidx], x_plotted[1:endidx,:]
    end
  end
end
