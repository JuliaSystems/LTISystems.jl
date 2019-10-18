@compat abstract type AbstractSignal{T<:Real,N} end

discontinuities(u, tspan::Tuple{T1,T2}) where {T1<:Real,T2<:Real} = promote_type(T1,T2)[]

_eltype(s::AbstractSignal{T,N}) where {T,N} = T
_length(s::AbstractSignal{T,N}) where {T,N} = N

+(s1::AbstractSignal{T1,N}, s2::AbstractSignal{T2,N}) where {T1,T2,N} = SumOfSignals(s1,  s2)
-(s1::AbstractSignal{T1,N}, s2::AbstractSignal{T2,N}) where {T1,T2,N} = SumOfSignals(s1, -s2)
