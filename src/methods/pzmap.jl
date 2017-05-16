immutable PZData{T}
  p::Vector{Complex128}
  z::Vector{Complex128}

  function (::Type{PZData{Val{U}}}){T<:Number,S<:Number,U}(p::AbstractVector{T},
    z::AbstractVector{S})
    new{Val{U}}(convert(Vector{Complex128}, p), convert(Vector{Complex128}, z))
  end
end

_pzmap(sys::LtiSystem{Val{:siso}}) = poles(sys), zeros(sys)
_pzmap(sys::LtiSystem{Val{:mimo}}) = poles(sys), tzeros(sys)

"""
    pzmap(sys) -> pzdata

Return pole-zero map `pzdata` for continuous- or discrete-time system model `sys`.

`pzdata` is a custom data type (`PZData`) containing

  * `poles(sys)` and `zeros(sys)`, when `sys` is a single-input, single-output
    system, and,
  * `poles(sys)` and `tzeros(sys)`, when `sys` is a multiple-input, multiple-output
    system.

Should the user want to plot `pzdata` herself, she can access the poles and zeros
of the system using the fields `pzdata.p` and `pzdata.z`, recpectively.

However, plotting recipes are defined for `pzdata`, which allows for
`plot(pzdata; <keyword arguments>)` when `using Plots`.

# Arguments

  * `stability::Bool = true`: when `true`, the stability border (imaginary axis
    for continuous-time and the unit circle for discrete-time systems) is drawn
    together with the poles and zeros.

**See also:** `poles`, `zeros`, `tzeros`.
"""
pzmap{T,S}(sys::LtiSystem{Val{T},Val{S}}) = PZData{Val{S}}(_pzmap(sys)...)

@recipe function f{T}(pzdata::PZData{Val{T}}; stability = true)
  @assert isa(stability, Bool) "plot(pzdata): `stability` must be a `Bool`"
  layout  := 1
  xscale  := :identity
  yscale  := :identity

  title         --> "Pole-Zero Map"
  xlabel        --> "Real Axis"
  ylabel        --> "Imaginary Axis"
  grid          --> true
  aspect_ratio  --> :equal
  legend        --> :topright

  realp   = real(pzdata.p)
  imagp   = imag(pzdata.p)
  realz   = real(pzdata.z)
  imagz   = imag(pzdata.z)

  maxy    = max(maxabs(imagp), maxabs(imagz))

  @series begin
    seriestype  :=  :scatter
    markershape --> :xcross
    label       --> "Pole"
    realp, imagp
  end

  @series begin
    seriestype  :=  :scatter
    markershape --> :circle
    label       --> "Zero"
    realz, imagz
  end

  if stability
    @series begin
      seriestype  :=  :path
      linestyle   --> :dash
      seriescolor --> :black
      label       --> ifelse(T == :disc, "Unit Circle", "Imaginary Axis")
      if T == :disc
        mycirc = [exp(1im*ω) for ω in linspace(0,2π,100)]
        real(mycirc), imag(mycirc)
      else
        zeros(100), linspace(-1.2*maxy,1.2*maxy,100)
      end
    end
  end
end
