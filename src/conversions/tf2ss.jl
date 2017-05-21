# TODO use the rfd or lfd conversion depending on dimension of s
ss(s::TransferFunction{Val{:mimo},Val{:cont}}) = minreal(ss(rfd(_tf2rfd(s)...)))
ss(s::TransferFunction{Val{:mimo},Val{:disc}}) = minreal(ss(rfd(_tf2rfd(s)..., s.Ts)))

ss(s::TransferFunction{Val{:siso},Val{:cont}}) = minreal(ss(_tf2ss(s)...))
ss(s::TransferFunction{Val{:siso},Val{:disc}}) = minreal(ss(_tf2ss(s)..., s.Ts))

# TODO num(s) instead of num(s.mat[1])
# controllable canonical form
function _tf2ss(s::TransferFunction{Val{:siso}})
  nump  = RationalFunctions.num(s.mat[1])
  denp  = RationalFunctions.den(s.mat[1])
  nump /= denp[end]
  denp /= denp[end]
  m     = degree(nump)
  n     = degree(denp)

  a     = coeffs(denp)
  b     = zeros(a)

  b[1:m+1] = coeffs(nump) # nump might be of lower order than denp

  A     = diagm(ones(eltype(b), n-1), -1)
  B     = zeros(eltype(b), n, 1)
  C     = zeros(eltype(b), 1, n)

  A[1,:]  = -reverse(a[1:end-1])
  B[1]    = one(eltype(b))
  C[:]    = reverse(b[1:end-1])
  D       = b[end]

  A, B, C, D
end
