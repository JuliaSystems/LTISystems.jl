struct SumOfSignals{T,N} <: AbstractSignal{T,N}
  signals::Vector{AbstractSignal{T,N}}

  function (::Type{SumOfSignals})(signals::Vector{T1}) where {T1<:AbstractSignal}
    if isempty(signals)
      warn("SumOfSignals(signals): `signals` is empty")
      throw(DomainError())
    end
    T = _eltype(signals[1])
    N = _length(signals[1])
    for idx in 2:length(signals)
      if N ≠ _length(signals[idx])
        warn("SumOfSignals(signals): `signals` has signal values with different lengths")
        throw(DomainError())
      end
      T = promote_type(T, _eltype(signals[2]))
    end
    new{T,N}(convert(Vector{AbstractSignal{T,N}}, signals))
  end

  function (::Type{SumOfSignals})(s1::AbstractSignal, s2::AbstractSignal...)
    signals = [s1, s2...]
    SumOfSignals(signals)
  end

  function (sos::SumOfSignals{T,N})(t::Real, x = nothing) where {T,N}
    temp = zeros(T,N)
    for signal in sos.signals
      temp += signal(t, x)
    end
    temp
  end
end

function discontinuities(sos::SumOfSignals{T1},
  tspan::Tuple{T2,T3}) where {T1,T2<:Real,T3<:Real}
  tstops = Set{T1}()
  for signal in sos.signals
    push!(tstops, discontinuities(signal, tspan)...)
  end
  [tstop for tstop in tstops]
end

convert(::Type{AbstractSignal{T,N}}, sos::SumOfSignals{T,N}) where {T,N}       = sos
convert(::Type{AbstractSignal{T1,N}}, sos::SumOfSignals{T2,N}) where {T1,T2,N} =
  SumOfSignals(convert(Vector{AbstractSignal{T1,N}}, sos.signals))

-(sos::SumOfSignals) = SumOfSignals(-signals)

# Relation between `Real`s

function +(x::Union{T,AbstractVector{T}}, sos::SumOfSignals) where {T<:Real}
  temp = [signal + x for signal in sos.signals]
  SumOfSignals(temp)
end
+(sos::SumOfSignals, x::Union{T,AbstractVector{T}}) where {T<:Real} = +(x,  sos)
-(x::Union{T,AbstractVector{T}}, sos::SumOfSignals) where {T<:Real} = +(x, -sos)
-(sos::SumOfSignals, x::Union{T,AbstractVector{T}}) where {T<:Real} = +(sos, -x)

function *(x::Real, sos::SumOfSignals)
  temp = [x*signal for signal in sos.signals]
  SumOfSignals(temp)
end
*(sos::SumOfSignals, x::Real) = *(x,   sos)
/(sos::SumOfSignals, x::Real) = *(sos, 1/x)

# Relationship between `AbstractSignal`s

+(sos::SumOfSignals{T1,N}, s::AbstractSignal{T2,N}) where {T1,T2,N} =
  SumOfSignals(push!(copy(sos.signals), s))
+(s::AbstractSignal{T2,N}, sos::SumOfSignals{T1,N}) where {T1,T2,N} = +(sos,  s)
-(s::AbstractSignal{T2,N}, sos::SumOfSignals{T1,N}) where {T1,T2,N} = +(s, -sos)
-(sos::SumOfSignals{T1,N}, s::AbstractSignal{T2,N}) where {T1,T2,N} = +(sos, -s)

# Relationship among `SumOfSignals`

+(sos1::SumOfSignals{T1,N}, sos2::SumOfSignals{T2,N}) where {T1,T2,N} =
  SumOfSignals(AbstractSignal{promote_type(T1,T2),N}[sos1.signals..., sos2.signals...])
-(sos1::SumOfSignals{T1,N}, sos2::SumOfSignals{T2,N}) where {T1,T2,N} =
  SumOfSignals(AbstractSignal{promote_type(T1,T2),N}[sos1.signals..., -sos2.signals...])
