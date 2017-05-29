immutable PRBS{T,N,S1,S2} <: AbstractSignal{T,N}
  mag::Vector{T}
  dc::Vector{T}
  internal::Vector{Tuple{S1,S1,S1,S1}}

  function (::Type{PRBS}){T1<:Real,T2<:Real,T3<:Integer,N}(mag::Vector{T1},
    dc::Vector{T2} = zeros(mag), seed::Vector{T3} = ones(Int, dc),
    ::Type{Val{N}} = Val{7})
    if !(length(mag) == length(dc) == length(seed))
      warn("PRBS(mag, dc, seed, Val{N}): `mag`, `dc` and `seed` must all have the same lengths")
      throw(DomainError())
    elseif N ∉ [7, 9, 11, 15, 20, 23, 31]
      warn("PRBS(mag, dc, seed, Val{N}): `N` must be either of {7, 9, 11, 15, 20, 23, 31}")
      throw(DomainError())
    elseif any(seed .≤ zero(T3))
      warn("PRBS(mag, dc, seed, Val{N}): `seed` must have positive values")
      throw(DomainError())
    end

    S1  = N ≤ 8 ? UInt8 : N ≤ 16 ? UInt16 : UInt32
    T   = promote_type(T1, T2)

    internal = [convert(Tuple{S1,S1,S1,S1}, (_coeffs(Val{N})..., S1(s), S1(0))) for s in seed]

    new{T,length(mag),S1,Val{N}}(convert(Vector{T}, mag), convert(Vector{T}, dc),
      internal)
  end

  function (::Type{PRBS}){T1<:Real,T2<:Real,S1<:Integer,N}(mag::Vector{T1},
    dc::Vector{T2}, internal::Vector{Tuple{S1,S1,S1,S1}}, ::Type{Val{N}})
    if !(length(mag) == length(dc) == length(internal))
      warn("PRBS(mag, dc, internal, Val{N}): `mag`, `dc` and `internal` must all have the same lengths")
      throw(DomainError())
    elseif N ∉ [7, 9, 11, 15, 20, 23, 31]
      warn("PRBS(mag, dc, internal, Val{N}): `N` must be either of {7, 9, 11, 15, 20, 23, 31}")
      throw(DomainError())
    else
      for idx in 1:length(internal)
        a, b, seed, _ = internal[idx]
        if (a, b) ≠ _coeffs(Val{N})
          warn("PRBS(mag, dc, internal, Val{N}): coefficient mismatch in `internal[$(idx)]`")
          throw(DomainError())
        elseif seed == zero(seed)
          warn("PRBS(mag, dc, internal, Val{N}): `seed` is zero in `internal[$(idx)]`")
          throw(DomainError())
        end
      end
    end
    T = promote_type(T1, T2)
    new{T,length(internal),S1,Val{N}}(convert(Vector{T}, mag), convert(Vector{T}, dc), internal)
  end

  (::Type{PRBS}){N}(mag::Real, dc::Real = zero(mag), seed::Integer = 1,
    t::Type{Val{N}} = Val{7}) = PRBS([mag], [dc], [seed], t)
  (::Type{PRBS})() = PRBS(1.)

  function (p::PRBS{T,N,S1}){T,N,S1}()
    temp = zeros(S1, N)
    for idx in 1:length(p.internal)
      a, b, val, _    = p.internal[idx]
      newbit          = ((val >> a) $ (val >> b) & 1)
      val             = ((val << 1) | newbit) & (2 << a - 1)
      p.internal[idx] = (a, b, val, 0)
      temp[idx]       = val
    end
    temp
  end

  function (p::PRBS{T,N,S1,Val{S2}}){T,N,S1,S2}(t::Real, x = nothing)
    temp  = zeros(Bool, N)
    for idx in 1:length(temp)
      a, b, val, k = p.internal[idx]
      if k ≥ S2
        p()
        k   = 0
        val = p.internal[idx][3]
      end
      bit = (val << k) & (1 << (S2-1))
      p.internal[idx] = (a, b, val, k + 1)
      temp[idx]       = bit >> (S2-1)
    end
    p.mag .* temp + p.dc
  end
end

_coeffs(::Type{Val{7}})   = (6  ,  5)
_coeffs(::Type{Val{9}})   = (8  ,  4)
_coeffs(::Type{Val{11}})  = (10 ,  8)
_coeffs(::Type{Val{15}})  = (14 , 13)
_coeffs(::Type{Val{20}})  = (19 ,  2)
_coeffs(::Type{Val{23}})  = (22 , 17)
_coeffs(::Type{Val{31}})  = (30 , 27)

# # Add the below code to tests
# for N in [7, 9, 11, 15, 20, 23, 31]
#   seed    = 2
#   mypr    = Signals.PRBS(1, 0, seed, Val{N})
#   period  = 1
#   while [seed] ≠ mypr()
#     period += 1
#   end
#   println("Period for N = $(N) is $(period)")
# end

# PRBS is discontinuous everywhere. I do not want to pollute discontinuities vector.
# discontinuities{T1<:Real,T2<:Real}(p::PRBS, tspan::Tuple{T1,T2}) =
#   promote_type(T1, T2)[]

convert{T,N,S1,S2}(::Type{AbstractSignal{T,N}}, p::PRBS{T,N,S1,Val{S2}})        = p
convert{T1,T2,N,S1,S2}(::Type{AbstractSignal{T1,N}}, p::PRBS{T2,N,S1,Val{S2}})  =
  PRBS(convert(Vector{T1}, p.mag), convert(Vector{T1}, p.dc), internal, Val{S2})

-{T,N,S1,S2}(p::PRBS{T,N,S1,Val{S2}}) = PRBS(-p.mag, -p.dc, p.internal, Val{S2})

# Relation between `Real`s

+{T1<:Real,T2,N,S1,S2}(x::Union{T1,AbstractVector{T1}}, p::PRBS{T2,N,S1,Val{S2}}) =
  PRBS(p.mag, p.dc + x, p.internal, Val{S2})
+{T<:Real}(p::PRBS, x::Union{T,AbstractVector{T}})  = +(x,  p)
-{T<:Real}(x::Union{T,AbstractVector{T}}, p::PRBS)  = +(x, -p)
-{T<:Real}(p::PRBS, x::Union{T,AbstractVector{T}})  = +(p, -x)

*{T,N,S1,S2}(x::Real, p::PRBS{T,N,S1,Val{S2}}) = PRBS(x*p.mag, x*p.dc, p.internal, Val{S2})
*(p::PRBS, x::Real) = *(x,   p)
/(p::PRBS, x::Real) = *(p, 1/x)
