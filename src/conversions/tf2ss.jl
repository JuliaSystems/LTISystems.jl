# TODO use the rfd or lfd conversion depending on dimension of s
ss(s::RationalTF{Val{:mimo},Val{:cont}}) = minreal(ss(rfd(_tf2rfd(s)...)))
ss(s::RationalTF{Val{:mimo},Val{:disc}}) = minreal(ss(rfd(_tf2rfd(s)..., s.Ts)))

ss(s::RationalTF{Val{:siso},Val{:cont}}) = minreal(ss(_tf2ss(s)...))
ss(s::RationalTF{Val{:siso},Val{:disc}}) = minreal(ss(_tf2ss(s)..., s.Ts))

# TODO num(s) instead of num(s.mat[1])
# TODO assumes proper tf
function _tf2ss(s::RationalTF{Val{:siso}})
  nump  = RationalFunctions.num(s.mat[1])
  denp  = RationalFunctions.den(s.mat[1])
  nump /= denp[end]
  denp /= denp[end]
  m     = degree(nump)
  n     = degree(denp)
  A     = diagm(ones(eltype(nump), n-1), 1)
  B     = zeros(eltype(nump), n, 1)
  C     = zeros(eltype(nump), 1, n)

  A[end,:] = -coeffs(denp)[1:end-1]
  B[end]   = one(eltype(nump))
  C[1:m+1] = coeffs(nump)
  D        = 0
  A, B, C, D
end
