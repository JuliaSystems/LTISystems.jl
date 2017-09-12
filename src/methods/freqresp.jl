"""
    freqresp(sys, ω) -> fr

Return the frequency response `fr` of a continuous- or discrete-time system `sys`
over the `Real` frequency value(s) `ω`.

The return type depends on the input pair `sys` and `ω`:

  * When `sys` is a single-input single-output system and `ω` is a scalar, `fr`
    is a scalar,
  * When `sys` is a single-input single-output system and `ω` is an array, `fr`
    is an array with the same dimensions as `ω`,
  * When `sys` is a multiple-input multiple-output system and `ω` is a scalar,
    `fr` is a matrix with the same simensions as `size(sys)`,
  * When `sys` is a multiple-input multiple-output system and `ω` is an array,
    `fr` is an array with the same dimensions as `(size(sys)..., size(ω)...)`.

Shorthand function call notations are defined for frequency-related responses:

  * `sys(x)` evaluates the system `sys` at `Number`s or `Array`s of `Number`s `x`,
    and,
  * `sys(ω = ω₀)` evaluates the system `sys` at `Real` or `Array`s of `Real`
    frequency values `ω₀`, *i.e.*, `sys(ω = ω₀) = freqresp(sys, ω)`.

Please note that `freqresp(sys, ω₀) = sys(1im*ω₀)` for continuous-time systems,
and `freqresp(sys, ω₀) = sys(exp(1im*samplingtime(sys)*ω₀))` for discrete-time
systems.

**See also:** `bode`, `nyquist`, `samplingtime`.
"""
freqresp{S}(sys::LtiSystem{Val{S},Val{:cont}}, x::Real)                     =
  sys(1im*x)
freqresp{S,M<:Real}(sys::LtiSystem{Val{S},Val{:cont}}, X::AbstractArray{M}) =
  sys(1im*X)
freqresp{S}(sys::LtiSystem{Val{S},Val{:disc}}, x::Real)                     =
  sys(exp(1im*samplingtime(sys)*x))
freqresp{S,M<:Real}(sys::LtiSystem{Val{S},Val{:disc}}, X::AbstractArray{M}) =
  sys(exp(1im*samplingtime(sys)*X))

# Implements algorithm found in:
# Laub, A.J., "Efficient Multivariable Frequency Response Computations",
# IEEE Transactions on Automatic Control, AC-26 (1981), pp. 407-408.
function _eval(sys::StateSpace, x::Number)
  R = x*I - sys.A
  convert(AbstractMatrix{Complex128}, sys.D + sys.C*(R\sys.B))
end

function _eval{M<:Number}(sys::StateSpace, X::AbstractArray{M})
  resp = Array{Complex128}(size(sys)..., size(X)...)
  for idx in CartesianRange(indices(X))
    resp[:, :, idx] = sys.D + sys.C*((X[idx]*I - sys.A)\sys.B)
  end
  resp
end
# function _eval{M<:Number}(sys::StateSpace, X::AbstractArray{M})
#   # Transform system into Hessenberg form
#   resp = Array(Complex128, size(sys)..., size(X)...)
#
#   sysr        = minreal(sys)
#   A, B, C, D  = sysr.A, sysr.B, sysr.C, sysr.D
#
#   if isempty(A)
#     for idx in CartesianRange(indices(X))
#       resp[:, :, idx] = D + C*((X[idx]*I - A)\B)
#     end
#     return resp
#   end
#
#   F           = hessfact(A)
#   Ah          = F[:H]
#   T           = full(F[:Q])
#   Bh          = T*B
#   Ch          = C/T
#
#   for idx in CartesianRange(indices(X))
#     R           = X[idx]*I - Ah
#     ipiv, info  = luhessfact!(R)
#     if info > 0
#       # s[i] is a pole of the system
#       resp[:, :, idx] = D + Ch*((X[idx]*I - Ah)\Bh)
#     else
#       temp = convert(AbstractMatrix{eltype(R)}, Bh)
#       hesssolve!(R, ipiv, temp)
#       resp[:, :, idx] = D + Ch*temp
#     end
#   end
#   resp
# end

# """
#     luhessfact!(H) -> ipiv, info
#
# Compute an LU factorization of an `n×m` upper Hessenberg matrix `H` using partial
# pivoting with row interchanges.
#
# Input/Output Parameters
#
#   * `H`: Input/output, `n×m` matrix. At input, the `n×m` upper Hessenberg matrix
#     to be factored. At output, the factors `L` and `U` from the factorization
#     `H = P*L*U`; the unit diagonal elements of `L` are not stored, and `L` is lower
#     bidiagonal.
#   * `ipiv`: Output, Integer vector of length `m`. The pivot indices; for
#     `1 ≤ i ≤ m`, row `i` of the matrix was interchanged with row `ipiv(i)`.
#   * `info`: Output, Integer. `info = 0` implies successful exit. `info = i > 0`
#     implies `U(i,i)` is exactly zero. The factorization has been completed, but
#     the factor `U` is exactly singular, and division by zero will occur if it is
#     used to solve a system of equations.
# """
# function luhessfact!{T<:Number}(H::AbstractMatrix{T})
# end

# """
#     hesssolve!(H, ipiv, B)`
#
# Solve the system of linear equations `H * X = B` with an upper Hessenberg `N×N`
# matrix `H` using the `LU` factorization computed by `luhessfact!`.
#
# Input/Output Parameters
#
#   * `H`: Input, `LDH×N` matrix, containing the factors `L` and `U` from the
#     factorization `H = P*L*U` as computed by `luhessfact!`.
#   * `ipiv`: Input, Integer vector of dimension `N`, containing the pivot indices
#     from `luhessfact!`; for `1 ≤ i ≤ N`, row `i` of the matrix was interchanged
#     with row `ipiv(i)`.
#   * `B`: Input/output, `LDB×NRHS` matrix. At input, the right hand side matrix
#     `B`. At output, the solution matrix `X`.
# """
# Solve L * X = B, overwriting B with X.
#
# L is represented as a product of permutations and unit lower
# triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
# where each transformation L(i) is a rank-one modification of
# the identity matrix.
# function hesssolve!{T<:Number}(H::AbstractMatrix{T},
#   ipiv::AbstractVector{Int}, B::AbstractMatrix{T})
# end
