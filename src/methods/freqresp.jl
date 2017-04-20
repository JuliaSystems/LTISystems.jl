"""
`fr, ω = freqresp(sys,ω)`

Evaluate the frequency response of continuous- or discrete-time systems

`G(jω) = C*((jωI -A)^-1)*B + D`
`G(e^jωTs) = C*((e^{jω}I -A)^-1)*B + D`

of system `sys` over the frequency value or vector `ω`.
`freqresp` returns `ω` together with the frequency response `fr` only if `ω` is
a vector. Otherwise, only `fr` is returned.
"""
freqresp{S}(sys::LtiSystem{Val{S},Val{:cont}}, x::Real)                     =
  sys(1im*x)
freqresp{S,M<:Real}(sys::LtiSystem{Val{S},Val{:cont}}, X::AbstractArray{M}) =
  sys(1im*X)
freqresp{S}(sys::LtiSystem{Val{S},Val{:disc}}, x::Real)                     =
  sys(exp(1im*samplingtime(sys)*x))
freqresp{S,M<:Real}(sys::LtiSystem{Val{S},Val{:disc}}, X::AbstractArray{M}) =
  sys(exp(1im*samplingtime(sys)*X))

# Evaluate system models at given complex numbers
@compat (sys::RationalTF{Val{:siso}})(x::Number)                      =
  (f = sys.mat[1]; convert(Complex128, f(x)))
@compat (sys::RationalTF{Val{:siso}}){M<:Number}(X::AbstractArray{M}) =
  (f = sys.mat[1]; Complex128[f(x) for x in X])
@compat (sys::RationalTF{Val{:mimo}})(x::Number)                      =
  Complex128[f(x) for f in sys.mat]
@compat (sys::RationalTF{Val{:mimo}}){M<:Number}(X::AbstractArray{M}) =
  Complex128[f(x) for f in sys.mat, x in X]

# Implements algorithm found in:
# Laub, A.J., "Efficient Multivariable Frequency Response Computations",
# IEEE Transactions on Automatic Control, AC-26 (1981), pp. 407-408.
function _eval(sys::StateSpace, x::Number)
  R = x*I - sys.A
  convert(AbstractMatrix{Complex128}, sys.D + sys.C*(R\sys.B))
end

function _eval{M<:Number}(sys::StateSpace, X::AbstractArray{M})
  # Transform system into Hessenberg form
  resp = Array(Complex128, size(sys)..., size(X)...)

  sysr        = minreal(sys)
  A, B, C, D  = sysr.A, sysr.B, sysr.C, sysr.D

  if isempty(A)
    for idx in CartesianRange(indices(X))
      resp[:, :, idx] = D + C*((X[idx]*I - A)\B)
    end
    return resp
  end

  F           = hessfact(A)
  Ah          = F[:H]
  T           = full(F[:Q])
  Bh          = T*B
  Ch          = C/T

  for idx in CartesianRange(indices(X))
    R           = X[idx]*I - Ah
    ipiv, info  = luhessfact!(R)
    if info > 0
      # s[i] is a pole of the system
      resp[:, :, idx] = D + Ch*((X[idx]*I - Ah)\Bh)
    else
      temp = convert(AbstractMatrix{eltype(R)}, Bh)
      hesssolve!(R, ipiv, temp)
      resp[:, :, idx] = D + Ch*temp
    end
  end
  resp
end

@compat (sys::StateSpace{Val{:siso}})(x::Number)                      =
  _eval(sys, x)[1]
@compat (sys::StateSpace{Val{:siso}}){M<:Number}(X::AbstractArray{M}) =
  reshape(_eval(sys, X), size(X))
@compat (sys::StateSpace{Val{:mimo}})(x::Number)                      =
  _eval(sys, x)
@compat (sys::StateSpace{Val{:mimo}}){M<:Number}(X::AbstractArray{M}) =
  _eval(sys, X)

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
