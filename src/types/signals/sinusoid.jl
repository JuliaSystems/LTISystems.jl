immutable Sinusoid{T,N} <: AbstractSignal{T,N}
  magnitude::Vector{Vector{T}}
  dc::Vector{Vector{T}}
  shift::Vector{Vector{T}}
  frequency::Vector{Vector{T}}

  function (::Type{Sinusoid}){T1<:Real,T2<:Real,T3<:Real,T4<:Real}(mag::Vector{Vector{T1}},
    dc::Vector{Vector{T2}}    = [zeros(elem) for elem in mag],
    shift::Vector{Vector{T3}} = [zeros(elem) for elem in mag],
    freq::Vector{Vector{T4}}  = [ones(elem) for elem in mag])
    if !(length(mag) == length(dc) == length(shift) == length(freq))
      warn("Sinusoid(mag, dc, shift, freq): `mag`, `dc`, `shift` and `freq` must all have the same lengths")
      throw(DomainError())
    else
      for idx in 1:length(mag)
        if !(length(mag[idx]) == length(dc[idx]) == length(shift[idx]) == length(freq[idx]))
          warn("Sinusoid(mag, dc, shift, freq): `mag`, `dc`, `shift` and `freq`, at every index, must have vectors of the same length")
          throw(DomainError())
        end
      end
    end
    T = promote_type(T1, T2, T3, T4)
    N = length(mag)

    new{T,N}(convert(Vector{Vector{T}}, mag), convert(Vector{Vector{T}}, dc),
          convert(Vector{Vector{T}}, shift), convert(Vector{Vector{T}}, freq))
  end

  (::Type{Sinusoid}){T1<:Real,T2<:Real,T3<:Real,T4<:Real}(mag::Vector{T1},
    dc::Vector{T2} = zeros(mag), shift::Vector{T3} = zeros(mag),
    freq::Vector{T4} = ones(mag)) = Sinusoid([[elem] for elem in mag],
    [[elem] for elem in dc], [[elem] for elem in shift], [[elem] for elem in freq])

  (::Type{Sinusoid})(mag::Real, dc::Real = zero(mag), shift::Real = zero(mag),
    freq::Real = one(mag)) = Sinusoid([[mag]], [[dc]], [[shift]], [[freq]])
  (::Type{Sinusoid})() = Sinusoid(1.)

  function (s::Sinusoid)(t::Real, x = nothing)
    [sum(m .* sin(2π*f*t + s) + dc) for (m, f, s, dc) in zip(s.mag, s.freq, s.shift, s.dc)]
  end
end

function discontinuities{T1<:Real,T2<:Real}(s::Sinusoid, tspan::Tuple{T1,T2})
  resp = s(tspan[1])
  return (resp ≈ zeros(resp)) ? T1[] : [tspan[1]]
end

convert{T,N}(::Type{AbstractSignal{T,N}}, s::Sinusoid{T,N})       = s
convert{T1,T2,N}(::Type{AbstractSignal{T1,N}}, s::Sinusoid{T2,N}) =
  Sinusoid(convert(Vector{Vector{T1}}, s.magnitude),
           convert(Vector{Vector{T1}}, s.dc),
           convert(Vector{Vector{T1}}, s.shift),
           convert(Vector{Vector{T1}}, s.frequency))

-(s::Sinusoid) = Sinusoid(-s.magnitude, -s.dc, s.shift, s.frequency)

# Relation between `Real`s

+{T<:Real}(x::Union{T,AbstractVector{T}}, s::Sinusoid) =
  Sinusoid(s.magnitude, s.dc + x, s.shift, s.frequency)
+{T<:Real}(s::Sinusoid, x::Union{T,AbstractVector{T}}) = +(x,  s)
-{T<:Real}(x::Union{T,AbstractVector{T}}, s::Sinusoid) = +(x, -s)
-{T<:Real}(s::Sinusoid, x::Union{T,AbstractVector{T}}) = +(s, -x)

*(x::Real, s::Sinusoid) = Sinusoid(x*s.magnitude, x*s.dc, s.shift, s.frequency)
*(s::Sinusoid, x::Real) = *(x,   s)
/(s::Sinusoid, x::Real) = *(s, 1/x)
