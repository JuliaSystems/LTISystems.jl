struct Ramp{T,N} <: AbstractSignal{T,N}
  time::Vector{T}
  rate::Vector{T}
  dc::Vector{T}

  function (::Type{Ramp})(rt::Vector{T1}, rr::Vector{T2},
    dc::Vector{T3} = zeros(rr)) where {T1<:Real,T2<:Real,T3<:Real}
    if length(rt) ≠ length(rr) || length(rt) ≠ length(dc)
      warn("Ramp(rt, rr, dc): `rt`, `rr` and `dc` must all have the same lengths")
      throw(DomainError())
    end

    T, N = promote_type(T1,T2,T3), length(rt)
    new{T,N}(convert(Vector{T}, rt), convert(Vector{T}, rr), convert(Vector{T}, dc))
  end
  (::Type{Ramp})(rt::Real, rr::Real, dc::Real = zero(rr)) = Ramp([rt], [rr], [dc])
  (::Type{Ramp})()                                        = Ramp([0.], [1.], [0.])

  function (r::Ramp{T,N})(t::Real, x = nothing) where {T,N}
    res = one(T) + one(t)*one(T)
    R   = typeof(res)

    R[t < rt ? dc : dc + (t-rt)*rr for (rt, rr, dc) in zip(r.time, r.rate, r.dc)]
  end
end

discontinuities(r::Ramp, tspan::Tuple{T1,T2}) where {T1<:Real,T2<:Real} =
  [rt for rt in r.time if tspan[1] ≤ rt ≤ tspan[2]]

convert(::Type{AbstractSignal{T,N}}, r::Ramp{T,N}) where {T,N} = r
convert(::Type{AbstractSignal{T1,N}}, r::Ramp{T2,N}) where {T1,T2,N} =
  Ramp(convert(Vector{T1}, r.time), convert(Vector{T1}, r.rate),
    convert(Vector{T1}, r.dc))

-(r::Ramp) = Ramp(r.time, -r.rate, -r.dc)

# Relation between `Real`s

+(x::Union{T,AbstractVector{T}}, r::Ramp) where {T<:Real} = Ramp(r.time, r.rate, r.dc + x)
+(r::Ramp, x::Union{T,AbstractVector{T}}) where {T<:Real} = +(x,  r)
-(x::Union{T,AbstractVector{T}}, r::Ramp) where {T<:Real} = +(x, -r)
-(r::Ramp, x::Union{T,AbstractVector{T}}) where {T<:Real} = +(r, -x)

*(x::Real, r::Ramp) = Ramp(r.time, x*r.rate, x*r.dc)
*(r::Ramp, x::Real) = *(x,   r)
/(r::Ramp, x::Real) = *(r, 1/x)
