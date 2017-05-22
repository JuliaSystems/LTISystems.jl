immutable Square{T,N} <: AbstractSignal{T,N}
  magnitude::Vector{T}
  width::Vector{T}
  dc::Vector{T}
  shift::Vector{T}
  period::Vector{T}

  function (::Type{Square}){T1,T2,T3,T4,T5}(mag::Vector{T1},
    width::Vector{T2} = ones(mag), dc::Vector{T3} = zeros(mag),
    shift::Vector{T4} = zeros(mag), per::Vector{T5} = width)
    if !(length(mag) == length(width) == length(dc) == length(shift) == length(per))
      warn("Square(mag, width, dc, shift, per): all inputs must have the same lengths")
      throw(DomainError())
    elseif any(width .< zero(T2))
      warn("Square(mag, width, dc, shift, per): `width` must have non-negative elements")
      throw(DomainError())
    elseif any(per .< width)
      warn("Square(mag, width, dc, shift, per): `width` cannot exceed `per`")
      throw(DomainError())
    end

    T, N = promote_type(T1, T2, T3, T4, T5), length(per)
    new{T,N}(convert(Vector{T}, mag), convert(Vector{T}, width),
            convert(Vector{T}, dc), convert(Vector{T}, shift),
            convert(Vector{T}, per))
  end

  (::Type{Square})(mag::Real, width::Real = one(mag), dc::Real = zero(mag),
    shift::Real = zero(mag), per::Real = width) = Square([mag], [width], [dc],
    [shift], [per])
  (::Type{Square})() = Square(1.)

  function (s::Square{T}){T}(t::Real, x = nothing)
    res   = divrem.(t - s.shift, s.period)
    tnorm = map(y->y[2], res)
    [(t < zero(t)) ? dc : (t < 0.5w) ? dc + m : (t < w) ? dc - m : dc for
      (t,dc,w,m) in zip(tnorm, s.dc, s.width, s.magnitude)]
  end
end

function discontinuities{T1<:Real,T2<:Real,T3<:Real}(s::Square{T1}, tspan::Tuple{T2,T3})
  tstops = Set{T1}()

  for i1 in 1:length(s.period)
    per   = s.period[i1]
    temp  = [s.shift[i1], s.shift[i1] + 0.5s.width[i1], s.shift[i1] + s.width[i1]]
    push!(tstops, temp...)

    d, r  = divrem(tspan[1], per)
    d     += sign(r)
    for i2 in 1:Int(abs(d))
      push!(tstops, (temp + i2*sign(d)*per)...)
    end

    d, r  = divrem(tspan[2], per)
    d     += sign(r)
    for i2 in 1:Int(abs(d))
      push!(tstops, (temp + i2*sign(d)*per)...)
    end
  end

  [tstop for tstop in tstops if tspan[1] ≤ tstop ≤ tspan[2]]
end

convert{T,N}(::Type{AbstractSignal{T,N}}, s::Square{T,N})       = s
convert{T1,T2,N}(::Type{AbstractSignal{T1,N}}, s::Square{T2,N}) =
  Square(convert(Vector{T1}, s.magnitude), convert(Vector{T1}, s.width),
    convert(Vector{T1}, s.dc), convert(Vector{T1}, s.shift),
    convert(Vector{T1}, s.period))

-(s::Square) = Square(-s.magnitude, s.width, -s.dc, s.shift, s.period)

# Relation between `Real`s

+{T<:Real}(x::Union{T,AbstractVector{T}}, s::Square)  =
  Square(s.magnitude, s.width, s.dc + x, s.shift, s.period)
+{T<:Real}(s::Square, x::Union{T,AbstractVector{T}})  = +(x,  r)
-{T<:Real}(x::Union{T,AbstractVector{T}}, s::Square)  = +(x, -r)
-{T<:Real}(s::Square, x::Union{T,AbstractVector{T}})  = +(r, -x)

*(x::Real, s::Square) = Square(x*s.magnitude, s.width, x*s.dc, s.shift, s.period)
*(s::Square, x::Real) = *(x,   r)
/(s::Square, x::Real) = *(r, 1/x)
