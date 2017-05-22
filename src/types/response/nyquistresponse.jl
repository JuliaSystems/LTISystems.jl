immutable NyquistResponse{T<:Real} <: SystemResponse
  freqs::Vector{T}  # rad/sec
  real::Array{T,3}   # -
  imag::Array{T,3} # -

  function (::Type{NyquistResponse}){T1<:Real,T2<:Real,T3<:Real}(
    freqs::AbstractVector{T1}, real::AbstractArray{T2,3}, imag::AbstractArray{T3,3})
    if size(real) ≠ size(imag)
      warn("NyquistResponse(freqs, real, imag): `real` and `imag` must have same dimensions")
      throw(DomainError())
    elseif length(freqs) ≠ size(real,3)
      warn("NyquistResponse(freqs, real, imag): `size(real,3) ≠ length(freqs)`")
      throw(DomainError())
    end

    T = promote_type(T1,T2,T3)

    new{T}(convert(Vector{T}, freqs), convert(Array{T,3}, real), convert(Array{T,3}, imag))
  end
end

_nyquist{T<:Real}(sys::LtiSystem{Val{:siso}}, ω::AbstractVector{T})  =
  reshape(freqresp(sys, ω), size(sys)..., length(ω))
_nyquist{T<:Real}(sys::LtiSystem{Val{:mimo}}, ω::AbstractVector{T})  =
  freqresp(sys, ω)

"""
    nyquist(sys, ω) -> nr

Return Nyquist response type `nr` for the system `sys` over the `Real` angular
frequency `Vector`, `ω`.

`nr` is a custom data type (`NyquistResponse`) containing

  * Frequency values (`nr.freqs`) in rad/sec, over which the frequency response
    is taken,
  * Real part (`nr.real`) of the frequency response, and,
  * Imaginary part (`nr.imag`) of the frequency response.

Plotting recipe is defined for `nr`, which allows for `plot(nr; <keyword arguments>)`
when `using Plots`.

# Arguments

  * `iopairs = Tuple{Int,Int}[]`: list of input-output pairs to draw the Nyquist
    plot of. When empty, all input-output pairs will be plotted. `iopairs` can
    be either a single `Tuple{Int,Int}` or a `Vector` of `Tuple{Int,Int}` values.
    Duplicate and out-of-range input-output pairs will be ignored, and,
  * `freqs::Symbol = :rads`: unit of frequencies, *i.e.*, `:rads` for rad/sec,
    `:Hz` for Herz.

**See also:** `freqresp`, `bode`.
"""
function nyquist{T}(sys::LtiSystem, ω::AbstractVector{T})
  # TODO: should we check for ω ≥ 0 ?
  fr = _nyquist(sys, ω)
  NyquistResponse(ω, real(fr), imag(fr))
end

nyquist{T}(sys::LtiSystem{Val{T},Val{:disc}}) =
  nyquist(sys, logspace(-6, log10(π/samplingtime(sys)), 1000))

@recipe function f(nr::NyquistResponse; iopairs = Tuple{Int,Int}[], freqs = :rads,
  mcircles = false)
  if !isa(iopairs, Vector{Tuple{Int,Int}}) && !isa(iopairs, Tuple{Int,Int})
    warn("plot(nr, <keyword arguments>): `iopairs` must be a `Tuple` or a `Vector{Tuple}`")
    throw(DomainError())
  elseif !isa(freqs, Symbol) && freqs ∉ [:rads, :Hz]
    warn("plot(nr, <keyword arguments>): `freqs` can be either `:rads` or `:Hz`")
    throw(DomainError())
  end

  iopairs = (iopairs == Tuple{Int,Int}) ? [iopairs] : iopairs

  iomap   = Dict{Int,Set{Int}}()
  ny, nu  = size(nr.real,1,2)

  for pair in iopairs
    outidx, inidx = pair
    if 1 ≤ outidx ≤ ny && 1 ≤ inidx ≤ nu
      iomap[outidx] = push!(get!(iomap, outidx, Set{Int}()), inidx)
    else
      warn("plot(nr, <keyword arguments>): $(pair) out of range. Discarding...")
    end
  end

  if isempty(iomap)
    for outidx = 1:ny
      iomap[outidx] = Set(1:nu)
    end
  end

  sortedmap = sort!([(key,[val for val in set]) for (key,set) in iomap],
                    lt = (l,r)->(l[1]≤r[1]))
  for pair in sortedmap
    sort!(pair[2])
  end

  # Define plotting rules
  layout      --> (1,1)
  link        --> :x
  xlabel      --> "Real Axis"
  ylabel      --> "Imaginary Axis"
  legend      --> :topright
  grid        --> true

  plottedfreqs = freqs == :rads ? nr.freqs : 0.5*nr.freqs/π

  for pair in sortedmap
    outidx = pair[1]
    for inidx in pair[2]
      @series begin
        primary := true
        label --> "\$G_{$(outidx),$(inidx)}\$"
        hover --> [plottedfreqs; reverse(plottedfreqs)]

        plottedreals  = nr.real[outidx, inidx, :]
        plottedimag   = nr.imag[outidx, inidx, :]

        [plottedreals; reverse(plottedreals)], [plottedimag;-reverse(plottedimag)]
      end
    end
  end

  if mcircles
    for pair in sortedmap
      outidx = pair[1]
      for inidx in pair[2]
        @series begin
          primary   := false
          label     --> ""
          linecolor --> :black
          linealpha --> 0.5
          linestyle --> :dot
          hover     --> false

          # iso modules
          modules = [-20., -10, -6, -4, -2, 2, 4, 6, 10, 20]
          moduler, modulec = _iso_modules(modules)

          # iso phases
          phases = [-90, -60, -45, -30, -15, 15, 30, 45, 60, 90]
          phaser, phasec = _iso_phases(phases)

          [moduler, phaser], [modulec, phasec]
        end
      end
    end
  end

end

# args specified in degrees
function _iso_phases{T1<:Real,T2<:Real}(args::AbstractVector{T1},
  w::AbstractVector{T2}=linspace(0, 2π, 200))

  N = tan(args*π/180)
  radius = sqrt(1+N.^2)./2N
  xc = -1/2
  yc = 1./2N

  c = cos(w)
  s = sin(w)

  reals = xc*ones(c*radius.') + c*radius.'
  imags = ones(s)*yc.' + s*radius.'
  reals, imags
end

# modules specified in db
function _iso_modules{T1<:Real,T2<:Real}(modules::AbstractVector{T1},
  w::AbstractVector{T2}=linspace(0, 2π, 200))

  M = exp10(modules/20)
  abs2(M)
  radius = M./(abs2(M)-1)

  xc = -M.*radius
  yc = 0

  c = cos(w)
  s = sin(w)

  reals = ones(c)*xc.' + c*radius.'
  imags = yc*ones(s*radius.') + s*radius.'
  reals, imags
end
