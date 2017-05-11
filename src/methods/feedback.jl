# Feedback interconnection
"""
    feedback(s1,s2)

Constructs the (negative) feedback interconnection of systems `s1` and `s2`.

# Examples
```julia
julia> s1 = tf([1],[1,2]);
s2 = tf([3],[1,2]);
s3 = feedback(s1,s2);
showall(s3)
Continuous time rational transfer function model
	y = Gu

  x + 2
------------
x^2 + 4x + 7
```
"""
function feedback{T1<:LtiSystem, T2<:LtiSystem}(s1::T1, s2::T2)
  /(s1, parallel(one(s1), series(s1,s2)))
end

feedback{T1<:LtiSystem, T2<:Real}(s1::T1, n::T2) =
  /(s1, parallel(one(s1), series(s1, n)))
feedback{T1<:LtiSystem, T2<:Real}(n::T2, s1::T1) =
  /(n, parallel(one(s1), series(n, s1)))

feedback{T1<:LtiSystem, T2<:Real}(s1::T1, n::Matrix{T2}) =
  /(s1, parallel(one(s1), series(s1, n)))
feedback{T1<:LtiSystem, T2<:Real}(n::Matrix{T2}, s1::T1) =
  /(n, parallel(one(s1), series(n, s1)))

# # Be smart in specific cases [30x speedup]
# function feedback{T1<:RationalTF{T,S,U,V},T2<:RationalTF}(s1::T1, s2::T2)
#   num = numpoly(s1)*denpoly(s2)
#   den = denpoly(s1)*denpoly(s2) + numpoly(s1)*numpoly(s2)
#   tf(num,den,samplingtime(s1))
# end
#
# function feedback{T1<:RationalTF,T2<:RationalTF}(s1::T1, s2::T2)
#   num = numpoly(s1)*denpoly(s2)
#   den = denpoly(s1)*denpoly(s2) + numpoly(s1)*numpoly(s2)
#   tf(num,den)
# end

# function feedback{T1<:DSisoZpk,T2<:DSisoZpk}(s1::T1, s2::T2)
#   z1,p1,k1 = zpkdata(s1)
#   z2,p2,k2 = zpkdata(s2)
#   z = vcat(z1, p2)
#   P = poly(p1)*poly(p2) + k1*k2*poly(z1)*poly(z2)
#   p = convert(Vector{Complex{Float64}},roots(P))
#   z_,p_,~ = rmcommon(z, p)
#   k = k1/real(P[end])
#   zpk(z_, p_, k, samplingtime(s1))
# end

# function feedback{T1<:DSisoZpk,T2<:DSisoZpk}(s1::T1, s2::T2)
#   s1_ = tf(real(numpoly(s1)), real(denpoly(s1)),samplingtime(s1))
#   s2_ = tf(real(numpoly(s1)), real(denpoly(s1)),samplingtime(s1))
#   feedback(s1_,s2_)
# end
