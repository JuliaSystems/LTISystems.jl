immutable Step{T,N} <: AbstractSignal{T,N}
  time::Vector{T}
  step::Vector{T}
  dc::Vector{T}

  function (::Type{Step}){T1<:Real,T2<:Real,T3<:Real}(st::Vector{T1}, sv::Vector{T2},
    dc::Vector{T3} = zeros(sv))
    if length(st) ≠ length(sv) || length(st) ≠ length(dc)
      warn("Step(st, sv, dc): `st`, `sv` and `dc` must all have the same lengths")
      throw(DomainError())
    end

    T, N = promote_type(T1,T2,T3), length(st)
    new{T,N}(convert(Vector{T}, st), convert(Vector{T}, sv), convert(Vector{T}, dc))
  end
  (::Type{Step})(st::Real, sv::Real, dc::Real = zero(sv)) = Step([st], [sv], [dc])
  (::Type{Step})()                                        = Step([0.], [1.], [0.])

  function (s::Step)(t::Real, x = nothing)
    [t < st ? dc : dc + sv for (st, sv, dc) in zip(s.time, s.step, s.dc)]
  end
end

discontinuities{T1<:Real,T2<:Real}(s::Step, tspan::Tuple{T1,T2}) =
  [st for st in s.time if tspan[1] ≤ st ≤ tspan[2]]

convert{T,N}(::Type{AbstractSignal{T,N}}, s::Step{T,N})       = s
convert{T1,T2,N}(::Type{AbstractSignal{T1,N}}, s::Step{T2,N}) =
  Step(convert(Vector{T1}, s.time), convert(Vector{T1}, s.step),
    convert(Vector{T1}, s.dc))

-(s::Step) = Step(s.time, -s.step, -s.dc)

# Relation between `Real`s

+{T<:Real}(x::Union{T,AbstractVector{T}}, s::Step) = Step(s.time, s.step, s.dc + x)
+{T<:Real}(s::Step, x::Union{T,AbstractVector{T}}) = +(x,  s)
-{T<:Real}(x::Union{T,AbstractVector{T}}, s::Step) = +(x, -s)
-{T<:Real}(s::Step, x::Union{T,AbstractVector{T}}) = +(s, -x)

*(x::Real, s::Step) = Step(s.time, x*s.step, x*s.dc)
*(s::Step, x::Real) = *(x,   s)
/(s::Step, x::Real) = *(s, 1/x)
