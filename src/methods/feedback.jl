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

# # TODO: Implement the specific cases for speedup purposes. Take into account the
# #       new type-system
# # Be smart in specific cases [30x speedup]
# function feedback{T1<:TransferFunction{T,S,U,V},T2<:TransferFunction}(s1::T1, s2::T2)
#   n = num(s1)*den(s2)
#   d = den(s1)*den(s2) + num(s1)*num(s2)
#   tf(n,d,samplingtime(s1))
# end
#
# function feedback{T1<:TransferFunction,T2<:TransferFunction}(s1::T1, s2::T2)
#   n = num(s1)*den(s2)
#   d = den(s1)*den(s2) + num(s1)*num(s2)
#   tf(n,d)
# end
