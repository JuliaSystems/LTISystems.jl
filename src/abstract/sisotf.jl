abstract SisoTf <: SisoSystem

# Printing functions
summary(s::SisoTf) = string("tf(nx=", numstates(s), (isdiscrete(s) ?
  string(",Ts=", samplingtime(s)) : ""), ")")

showcompact(io::IO, s::SisoTf)  = print(io, summary(s))
show(io::IO, s::SisoTf)         = print(io, summary(s))
showall(io::IO, s::SisoTf)      = print(io, summary(s))

function rmgcd{T1<:Number, T2<:Number}(p1::Poly{T1}, p2::Poly{T2})
  R = promote_type(T1,T2)
  gcdp1p2::Poly{R} = gcd(p1,p2)
  gcdp1p2_::Poly{R} = gcdp1p2/gcdp1p2[end]
  p1_::Poly = div(p1,gcdp1p2_)
  p2_::Poly = div(p2,gcdp1p2_)
  return (p1_,p2_,gcdp1p2_)
end

function rmcommon{T1<:AbstractVector, T2<:AbstractVector}(a::T1, b::T2)
  T     = promote_type(eltype(T1),eltype(T2))
  list  = Dict{T,Vector{Int}}()

  l1    = length(a)
  l2    = length(b)
  s1    = l1 > l2 ? b : a
  s2    = l1 > l2 ? a : b

  for elem in s1
    list[elem] = get(list, elem, Int[0, 0]) + [1, 0]
  end

  for elem in s2
    try
      list[elem] += [0, 1]
    end
  end

  temp::Vector{T} = []
  for key in keys(list)
    append!(temp, fill(key, min(list[key]...)))
  end
  ar::T1 = copy(a)
  br::T2 = copy(b)
  for elem in temp
    splice!(ar, findfirst(ar,elem))
    splice!(br, findfirst(br,elem))
  end
  return (ar,br,temp)
end
