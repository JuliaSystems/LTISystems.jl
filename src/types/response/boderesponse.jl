immutable BodeResponse{T<:Real} <: SystemResponse
  freqs::Vector{T}  # rad/sec
  mag::Array{T,3}   # abs, i.e., |G|
  phase::Array{T,3} # radians

  function (::Type{BodeResponse}){T1<:Real,T2<:Real,T3<:Real}(
    freqs::AbstractVector{T1}, mag::AbstractArray{T2,3}, phase::AbstractArray{T3,3})
    if size(mag) ≠ size(phase)
      warn("BodeResponse(freqs, mag, phase): `mag` and `phase` must have same dimensions")
      throw(DomainError())
    elseif length(freqs) ≠ size(mag,3)
      warn("BodeResponse(freqs, mag, phase): `size(mag,3) ≠ length(freqs)`")
      throw(DomainError())
    elseif any(mag .≤ zero(T2))
      warn("BodeResponse(freqs, mag, phase): `mag` values cannot be negative")
      throw(DomainError())
    end

    T = promote_type(T1,T2,T3)

    new{T}(convert(Vector{T}, freqs), convert(Array{T,3}, mag), convert(Array{T,3}, phase))
  end
end

_bode{T<:Real}(sys::LtiSystem{Val{:siso}}, ω::AbstractVector{T})  =
  reshape(freqresp(sys, ω), size(sys)..., length(ω))
_bode{T<:Real}(sys::LtiSystem{Val{:mimo}}, ω::AbstractVector{T})  =
  freqresp(sys, ω)

"""
    bode(sys, ω) -> br

Return Bode response type `br` for the system `sys` over the `Real` angular frequency
`Vector`, `ω`.

`br` is a custom data type (`BodeResponse`) containing

  * Frequency values (`br.freqs`) in rad/sec, over which the frequency response
    is taken,
  * Magnitude (`br.mag`) in absolute values, and,
  * Phase shift (`br.phase`) in radians.

Plotting recipe is defined for `br`, which allows for `plot(br; <keyword arguments>)`
when `using Plots`.

# Arguments

  * `iopairs = Tuple{Int,Int}[]`: list of input-output pairs to draw the Bode
    plot of. When empty, all input-output pairs will be plotted. `iopairs` can
    be either a single `Tuple{Int,Int}` or a `Vector` of `Tuple{Int,Int}` values.
    Duplicate and out-of-range input-output pairs will be ignored,
  * `freqs::Symbol = :rads`: unit of frequencies, *i.e.*, `:rads` for rad/sec,
    `:Hz` for Herz,
  * `magnitude::Symbol = :dB`: unit of the magnitude plot, *i.e.*, `:dB` for
    decibels, `:abs` for absolute values, and,
  * `phase::Symbol = :deg`: unit of phase shift, *i.e.*, `:deg` for degrees,
    `:rad` for radians.

**See also:** `freqresp`, `nyquist`.
"""
function bode{T}(sys::LtiSystem, ω::AbstractVector{T})
  # TODO: should we check for ω ≥ 0 ?
  fr = _bode(sys, ω)
  BodeResponse(ω, abs(fr), unwrap!(angle(fr)))
end

bode{T}(sys::LtiSystem{Val{T},Val{:disc}}) =
  bode(sys, logspace(-6, log10(π/samplingtime(sys)), 1000))

@recipe function f(br::BodeResponse; iopairs = Tuple{Int,Int}[], freqs = :rads,
  magnitude = :dB, phase = :deg)
  if !isa(iopairs, Vector{Tuple{Int,Int}}) && !isa(iopairs, Tuple{Int,Int})
    warn("plot(br, <keyword arguments>): `iopairs` must be a `Tuple` or a `Vector{Tuple}`")
    throw(DomainError())
  elseif !isa(freqs, Symbol) && freqs ∉ [:rads, :Hz]
    warn("plot(br, <keyword arguments>): `freqs` can be either `:rads` or `:Hz`")
    throw(DomainError())
  elseif !isa(magnitude, Symbol) && magnitude ∉ [:dB, :abs]
    warn("plot(br, <keyword arguments>): `magnitude` can be either `:dB` or `:abs`")
    throw(DomainError())
  elseif !isa(phase, Symbol) && phase ∉ [:deg, :rad]
    warn("plot(br, <keyword arguments>): `phase` can be either `:deg` or `:rad`")
    throw(DomainError())
  end

  iopairs = (iopairs == Tuple{Int,Int}) ? [iopairs] : iopairs

  iomap   = Dict{Int,Set{Int}}()
  ny, nu  = size(br.mag,1,2)

  for pair in iopairs
    outidx, inidx = pair
    if 1 ≤ outidx ≤ ny && 1 ≤ inidx ≤ nu
      iomap[outidx] = push!(get!(iomap, outidx, Set{Int}()), inidx)
    else
      warn("plot(br, <keyword arguments>): $(pair) out of range. Discarding...")
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
  layout      --> (2,1)
  link        --> :x
  xscale      --> :log10
  xlabel      --> ifelse(freqs == :rads, "Frequency, \$\\omega\$ (rad/sec)",
                        "Frequency, \$f\$ (Hz)")
  legend      --> :topright
  grid        --> true

  plottedfreqs = freqs == :rads ? br.freqs : 0.5*br.freqs/π

  for pair in sortedmap
    outidx = pair[1]
    for inidx in pair[2]
      @series begin
        ylabel  --> ifelse(magnitude == :dB, "Magnitude (dB)", "Magnitude (\$\\vert G \\vert\$)")
        label   --> "\$G_{$(outidx),$(inidx)}\$"
        plottedmag = (magnitude == :dB) ? 20*log10(vec(br.mag[outidx, inidx, :])) :
                                          vec(br.mag[outidx, inidx, :])
        plottedfreqs, plottedmag
      end
      @series begin
        ylabel  --> ifelse(phase == :deg, "Phase (deg)", "Phase (rad)")
        label   --> "\$G_{$(outidx),$(inidx)}\$"
        plottedphase = (phase == :deg) ? 180/π*vec(br.phase[outidx, inidx, :]) :
                                         vec(br.phase[outidx, inidx, :])
        plottedfreqs, plottedphase
      end
    end
  end
end
