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
evalfr(sys::LtiSystem{Siso{true}}, s::Number) = _evalfr(sys,s)[1]
evalfr(sys::LtiSystem{Siso{false}}, s::Number) = _evalfr(sys,s)
evalfr{M<:Number}(sys::LtiSystem{Siso{true}}, s::AbstractVector{M}) = _evalfr(sys,s)[1,1,:]
evalfr{M<:Number}(sys::LtiSystem{Siso{false}}, s::AbstractVector{M}) = _evalfr(sys,s)
evalfr(mat::AbstractMatrix, s::Number) = map(sys -> evalfr(sys, s), mat)

# Implements algorithm found in:
# Laub, A.J., "Efficient Multivariable Frequency Response Computations",
# IEEE Transactions on Automatic Control, AC-26 (1981), pp. 407-408.
function _evalfr{M<:Number}(sys::StateSpace, s::AbstractVector{M})
  ny, nu = size(sys)
  nw = length(s)

  # Transform system into Hessenberg form
  A, B, C, D = sys.A, sys.B, sys.C, sys.D
  F  = hessfact(A)
  Ah = F[:H]
  T  = full(F[:Q])
  Bh = convert(AbstractMatrix{Complex128},T*B)
  Ch = convert(AbstractMatrix{Complex128},C/T)
  D  = convert(AbstractMatrix{Complex128},D)

  resp = Array(Complex128, ny, nu, nw)
  for i = 1:nw
    R = s[i]*eye(sys.nx) - Ah
    ipiv, info = luhessfact!(R)
    if info > 0
      # s[i] is a pole of the system
      resp[:, :, i] = D + Ch*((s[i]*eye(sys.nx) - Ah)\Bh)
    else
      temp = copy(Bh)
      hesssolve!(R,ipiv,temp)
      resp[:, :, i] = D + Ch*temp
    end
  end
  resp
end

function _evalfr{M<:Number}(sys::LtiSystem, s::AbstractVector{M})
  ny, nu = size(sys)
  nw = length(s)
  resp = Array(Complex128, ny, nu, nw)
  for i = 1:nw
    resp[:, :, i] = evalfr(sys, s[i])
  end
  resp
end

function _evalfr(sys::StateSpace, s::Number)
  S = promote_type(typeof(s), Float64)
  R = Diagonal(S[s for i=1:sys.nx]) - sys.A
  sys.D + sys.C*((R\sys.B)::Matrix{S})
end

function _evalfr(sys::RationalTF{Siso{true}}, s::Number)
  S = promote_type(typeof(s), Float64)
  den = polyval(sys.den[1], s)
  if den == zero(S)
    convert(S, Inf)
  else
    polyval(sys.num[1], s)/den
  end
end
function _evalfr(sys::RationalTF, s::Number)
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

# `F(s)`
# Notation for frequency response evaluation:
# F(s) evaluates the LTI system F at the complex value or vector s.
#
# It is not possible to define this directly for LtiSystem in Julia 0.5, because it is an abstract type.
# WARNING: If one defines a new subtype of LtiSystem in a package requiring ControlCore.jl,
# these lines should be executed again.
for T in subtypes(LtiSystem)
  (sys::T)(s) = evalfr(sys,s)
end

# (sys::RationalTF)(s) = evalfr(sys,s)
# (sys::StateSpace)(s) = evalfr(sys,s)
# (sys::ZeroPoleGain)(s) = evalfr(sys,s)
# (sys::GeneralMimo)(s) = evalfr(sys,s)

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



"""
`ipiv, info = luhessfact!(H)`

To compute an LU factorization of an n x m upper Hessenberg matrix H
using partial pivoting with row interchanges.
Based on the SLICOT routine MB02SZ

Input/Output Parameters

H       (input/output) n x m matrix
        On entry, the n x m upper Hessenberg matrix to be factored.
        On exit, the factors L and U from the factorization H = P*L*U;
        the unit diagonal elements of L are not stored, and L is lower bidiagonal.

ipiv    (output) Integer vector of dimension m
        The pivot indices; for 1 <= i <= m, row i of the matrix
        was interchanged with row ipiv(i).

Error Indicator

info      Integer
= 0:      successful exit;
= i > 0:  U(i,i) is exactly zero. The
          factorization has been completed, but the factor U
          is exactly singular, and division by zero will occur
          if it is used to solve a system of equations.
"""
function luhessfact!{T<:Number}(H::AbstractMatrix{T})
  n, m = size(H)
  ipiv = zeros(Int,m)

  info = 0

  for j = 1:m

    # Find pivot and test for singularity.
    jp = j

    if j < m
      if abs(H[j+1,j]) > abs(H[j,j])
        jp = j + 1
      end
    end
    ipiv[j] = jp

    if H[jp,j] != zero(T)

      # Apply the interchange to columns J:N.
      if jp != j
        for i = j:m
          temp = H[j,i]
          H[j,i] = H[jp,i]
          H[jp,i] = temp
        end
      end

      # Compute element J+1 of J-th column.
      if j < m
        H[j+1, j] /= H[j,j]
      end
    else
      if info == 0
        info = j
      end
    end

    if j < m
      # Update trailing submatrix.
      H[j+1,j+1:m] -= H[j+1,j] * H[j,j+1:m]
    end
  end
  return ipiv, info
end

"""
`hesssolve!(H, ipiv, B)`

To solve the system of linear equations H * X = B
with an upper Hessenberg N-by-N matrix H using the LU
factorization computed by `luhessfact`.
Based on the SLICOT routine MB02RZ.

Input/Output Parameters

H      (input) Matrix LDH x N, containing
       the factors L and U from the factorization H = P*L*U
       as computed by `luhessfact!`.

ipiv   (input) Integer vector of dimension N, containing
       the pivot indices from `luhessfact!`; for 1<=i<=N, row i of the
       matrix was interchanged with row ipiv(i).

B      (input/output) Matrix LDB x NRHS
       On entry, the right hand side matrix B.
       On exit, the solution matrix X.
"""
function hesssolve!{T<:Number}(H::AbstractMatrix{T},
  ipiv::AbstractVector{Int}, B::AbstractMatrix{T})
  LDH, N = size(H)
  LDB, NRHS = size(B)

  # Solve L * X = B, overwriting B with X.
  #
  # L is represented as a product of permutations and unit lower
  # triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
  # where each transformation L(i) is a rank-one modification of
  # the identity matrix.
  for j = 1:N-1
    jp = ipiv[j]

    if jp != j
      for i = 1:NRHS
        temp = B[jp,i]
        B[jp,i] = B[j,i]
        B[j,i] = temp
      end
    end
    B[j+1,:] -= H[j+1,j]*B[j,:]
  end

  # Solve U * X = B, overwriting B with X.
  BLAS.trsm!('L', 'U', 'N', 'N',one(T), H, B)
end
