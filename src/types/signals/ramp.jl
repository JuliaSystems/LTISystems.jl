immutable Ramp{T,N} <: AbstractSignal{T,N}
  time::Vector{T}
  rate::Vector{T}
  dc::Vector{T}

  function (::Type{Ramp}){T1<:Real,T2<:Real,T3<:Real}(rt::Vector{T1}, rr::Vector{T2},
    dc::Vector{T3} = zeros(rr))
    if length(rt) ≠ length(rr) || length(rt) ≠ length(dc)
      warn("Ramp(rt, rr, dc): `rt`, `rr` and `dc` must all have the same lengths")
      throw(DomainError())
    end

    T, N = promote_type(T1,T2,T3), length(rt)
    new{T,N}(convert(Vector{T}, rt), convert(Vector{T}, rr), convert(Vector{T}, dc))
  end
  (::Type{Ramp})(rt::Real, rr::Real, dc::Real = zero(rr)) = Ramp([rt], [rr], [dc])
  (::Type{Ramp})()                                        = Ramp([0.], [1.], [0.])

  function (r::Ramp{T,N}){T,N}(t::Real, x = nothing)
    res = one(T) + one(t)*one(T)
    R   = typeof(res)

    R[t < rt ? dc : dc + (t-rt)*rr for (rt, rr, dc) in zip(r.time, r.rate, r.dc)]
  end
end

discontinuities{T1<:Real,T2<:Real}(r::Ramp, tspan::Tuple{T1,T2}) =
  [rt for rt in r.time if tspan[1] ≤ rt ≤ tspan[2]]

convert{T,N}(::Type{AbstractSignal{T,N}}, r::Ramp{T,N})       = r
convert{T1,T2,N}(::Type{AbstractSignal{T1,N}}, r::Ramp{T2,N}) =
  Ramp(convert(Vector{T1}, r.time), convert(Vector{T1}, r.rate),
    convert(Vector{T1}, r.dc))

-(r::Ramp) = Ramp(r.time, -r.rate, -r.dc)

# Relation between `Real`s

+{T<:Real}(x::Union{T,AbstractVector{T}}, r::Ramp) = Ramp(r.time, r.rate, r.dc + x)
+{T<:Real}(r::Ramp, x::Union{T,AbstractVector{T}}) = +(x,  r)
-{T<:Real}(x::Union{T,AbstractVector{T}}, r::Ramp) = +(x, -r)
-{T<:Real}(r::Ramp, x::Union{T,AbstractVector{T}}) = +(r, -x)

*(x::Real, r::Ramp) = Ramp(r.time, x*r.rate, x*r.dc)
*(r::Ramp, x::Real) = *(x,   r)
/(r::Ramp, x::Real) = *(r, 1/x)
