immutable Triangle{T,N} <: AbstractSignal{T,N}
  magnitude::Vector{T}
  width::Vector{T}
  dc::Vector{T}
  shift::Vector{T}
  period::Vector{T}

  function (::Type{Triangle}){T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real}(mag::Vector{T1},
    width::Vector{T2} = ones(mag), dc::Vector{T3} = zeros(mag),
    shift::Vector{T4} = zeros(mag), per::Vector{T5} = width)
    if !(length(mag) == length(width) == length(dc) == length(shift) == length(per))
      warn("Triangle(mag, width, dc, shift, per): all inputs must have the same lengths")
      throw(DomainError())
    elseif any(width .≤ zero(T2))
      warn("Triangle(mag, width, dc, shift, per): `width` must have positive elements")
      throw(DomainError())
    elseif any(per .< width)
      warn("Triangle(mag, width, dc, shift, per): `width` cannot exceed `per`")
      throw(DomainError())
    end

    T, N = promote_type(T1, T2, T3, T4, T5), length(per)
    new{T,N}(convert(Vector{T}, mag), convert(Vector{T}, width),
            convert(Vector{T}, dc), convert(Vector{T}, shift),
            convert(Vector{T}, per))
  end

  (::Type{Triangle})(mag::Real, width::Real = one(mag), dc::Real = zero(mag),
    shift::Real = zero(mag), per::Real = width) = Triangle([mag], [width], [dc],
    [shift], [per])
  (::Type{Triangle})() = Triangle(1.)

  function (tr::Triangle{T}){T}(t::Real, x = nothing)
    res   = divrem.(t - tr.shift, tr.period)
    tnorm = map(y->y[2], res)
    [(t < zero(t)) ? dc : (t < 0.25w) ? dc + 4m/w*t : (t < 0.75w) ? dc + m - 4m/w*(t-0.25w) :
      (t < w) ? dc - m + 4m/w*(t-0.75w) : dc for (t,dc,w,m) in zip(tnorm, tr.dc, tr.width,
      tr.magnitude)]
  end
end

function discontinuities{T1<:Real,T2<:Real,T3<:Real}(t::Triangle{T1}, tspan::Tuple{T2,T3})
  tstops = Set{T1}()

  for i1 in 1:length(t.period)
    per   = t.period[i1]
    temp  = [t.shift[i1], t.shift[i1] + 0.25t.width[i1], t.shift[i1] + 0.75t.width[i1],
            t.shift[i1] + t.width[i1]]
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

convert{T,N}(::Type{AbstractSignal{T,N}}, t::Triangle{T,N})       = t
convert{T1,T2,N}(::Type{AbstractSignal{T1,N}}, t::Triangle{T2,N}) =
  Triangle(convert(Vector{T1}, t.magnitude), convert(Vector{T1}, t.width),
    convert(Vector{T1}, t.dc), convert(Vector{T1}, t.shift),
    convert(Vector{T1}, t.period))

-(t::Triangle) = Triangle(-t.magnitude, t.width, -t.dc, t.shift, t.period)

# Relation between `Real`s

+{T<:Real}(x::Union{T,AbstractVector{T}}, r::Triangle)  =
  Triangle(t.magnitude, t.width, t.dc + x, t.shift, t.period)
+{T<:Real}(t::Triangle, x::Union{T,AbstractVector{T}})  = +(x,  t)
-{T<:Real}(x::Union{T,AbstractVector{T}}, t::Triangle)  = +(x, -t)
-{T<:Real}(t::Triangle, x::Union{T,AbstractVector{T}})  = +(t, -x)

*(x::Real, t::Triangle) = Triangle(x*t.magnitude, t.width, x*t.dc, t.shift, t.period)
*(t::Triangle, x::Real) = *(x,   t)
/(t::Triangle, x::Real) = *(t, 1/x)
