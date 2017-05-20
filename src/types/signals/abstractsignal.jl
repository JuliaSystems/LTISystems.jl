@compat abstract type AbstractSignal{T<:Real,N} end

discontinuities{T1<:Real,T2<:Real}(u, tspan::Tuple{T1,T2}) = promote_type(T1,T2)[]

_eltype{T,N}(s::AbstractSignal{T,N}) = T
_length{T,N}(s::AbstractSignal{T,N}) = N

+{T1,T2,N}(s1::AbstractSignal{T1,N}, s2::AbstractSignal{T2,N}) = SumOfSignals(s1,  s2)
-{T1,T2,N}(s1::AbstractSignal{T1,N}, s2::AbstractSignal{T2,N}) = SumOfSignals(s1, -s2)
