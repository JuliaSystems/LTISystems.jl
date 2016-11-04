"""
`fr, ω = freqresp(sys,ω)`

Evaluate the frequency response of continuous- or discrete-time systems

`G(jω) = C*((jωI -A)^-1)*B + D`
`G(e^jωTs) = C*((e^{jω}I -A)^-1)*B + D`

of system `sys` over the frequency value or vector `ω`.
`freqresp` returns `ω` together with the frequency response `fr` only if `ω` is
a vector. Otherwise, only `fr` is returned.
"""
freqresp{S,M<:Real}(sys::LtiSystem{S,Continuous{true}}, ω::AbstractVector{M}) = (_evalfr(sys,im*ω), ω)
freqresp{S,M<:Real}(sys::LtiSystem{S,Continuous{false}}, ω::AbstractVector{M}) = (_evalfr(sys,exp(ω*im*sys.Ts)), ω)
freqresp{S}(sys::LtiSystem{S,Continuous{true}}, ω::Real) = _evalfr(sys,im*ω)
freqresp{S}(sys::LtiSystem{S,Continuous{false}}, ω::Real) = _evalfr(sys,exp(ω*im*sys.Ts))

"""
Evaluate the frequency response

`H(s) = C*((sI - A)^-1)*B + D`

of system `sys` at the complex value or vector `s`.
"""
evalfr(sys::StateSpace{Siso{true}}, s::Number) = _evalfr(sys,s)[1]
evalfr(sys::StateSpace{Siso{false}}, s::Number) = _evalfr(sys,s)
evalfr{M<:Number}(sys::LtiSystem{Siso{true}}, s::AbstractVector{M}) = _evalfr(sys,s)[1,1,:]
evalfr{M<:Number}(sys::LtiSystem{Siso{false}}, s::AbstractVector{M}) = _evalfr(sys,s)

function _evalfr{M<:Number}(sys::StateSpace, s::AbstractVector{M})
  ny, nu = size(sys)
  nw = length(s)
  sys = preprocess_for_freqresp(sys)
  resp = Array(Complex128, ny, nu, nw)
  for i = 1:nw
    # TODO : This doesn't actually take advantage of Hessenberg structure
    # for statespace version.
    resp[:, :, i] = _evalfr(sys, s[i])
  end
  resp
end

function evalfr{M<:Number}(sys::LtiSystem, s::AbstractVector{M})
  ny, nu = size(sys)
  nw = length(s)
  resp = Array(Complex128, ny, nu, nw)
  for i = 1:nw
    resp[:, :, i] = evalfr(sys, s[i])
  end
  resp
end

# Implements algorithm found in:
# Laub, A.J., "Efficient Multivariable Frequency Response Computations",
# IEEE Transactions on Automatic Control, AC-26 (1981), pp. 407-408.
preprocess_for_freqresp{T}(sys::StateSpace{T,Continuous{true}}) = ss(_preprocess_for_freqresp(sys)...)
preprocess_for_freqresp{T}(sys::StateSpace{T,Continuous{false}}) = ss(_preprocess_for_freqresp(sys)...,sys.Ts)
preprocess_for_freqresp(sys::LtiSystem) = sys

function _preprocess_for_freqresp(sys::StateSpace)
    A, B, C, D = sys.A, sys.B, sys.C, sys.D
    F = hessfact(A)
    H = F[:H]::Matrix{Float64}
    T = full(F[:Q])
    P = C*T
    Q = T\B
    H, Q, P, D
end

function _evalfr(sys::StateSpace, s::Number)
  S = promote_type(typeof(s), Float64)
  R = Diagonal(S[s for i=1:sys.nx]) - sys.A
  sys.D + sys.C*((R\sys.B)::Matrix{S})
end

function evalfr(sys::RationalTF{Siso{true}}, s::Number)
  S = promote_type(typeof(s), Float64)
  den = polyval(sys.den[1], s)
  if den == zero(S)
    convert(S, Inf)
  else
    polyval(sys.num[1], s)/den
  end
end
function evalfr(sys::RationalTF, s::Number)
  S = promote_type(typeof(s), Float64)
  ny, nu = size(sys)
  res = Array(S, ny, nu)
  for j = 1:nu
    for i = 1:ny
      res[i, j] = evalfr(mat[i, j], s)
    end
  end
  return res
end

evalfr(mat::AbstractMatrix, s::Number) = map(sys -> evalfr(sys, s), mat)

"""
`F(s)`
Notation for frequency response evaluation:
F(s) evaluates the LTI system F at the complex value or vector s.
"""
# It is not possible to define this directly for LtiSystem in Julia 0.5, because it is an abstract type.
(sys::RationalTF)(s) = evalfr(sys,s)
(sys::StateSpace)(s) = evalfr(sys,s)
(sys::ZeroPoleGain)(s) = evalfr(sys,s)
(sys::GeneralMimo)(s) = evalfr(sys,s)

# function evalfr(sys::ZeroPoleGain, s::Number)
#     S = promote_type(typeof(s), Float64)
#     ny, nu = size(sys)
#     res = Array(S, ny, nu)
#     for i=1:ny, j=1:nu
#         res[i, j] = evalfr(sys.z[i,j], sys.p[i,j], sys.k[i,j], s)
#     end
#     return res
# end
# function evalfr(z::Vector{Complex128}, p::Vector{Complex128}, k::Float64, s::Number)
#     S = promote_type(typeof(s), Float64)
#     res = k
#     for zi in z
#         res *= s - zi
#     end
#     for pi in p
#         fact = s - pi
#         fact == zero(S) ? (return convert(S, Inf)) : (res /= fact)
#     end
#     return res
# end
