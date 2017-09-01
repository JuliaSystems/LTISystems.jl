#  Transformation from TransferFunction to state space is based on [1] and the
#  resulting realization is reduced using minreal.
#
#  [1]: R. V. Patel, "Computation of minimal-order state-space realizations
#        and observability indices using orthogonal transforms",
#        International Journal of Automatic Control, vol. 33, no. 2, pp.
#        227-246, Mar. 1981.

# s = R(s)*inv(D(s)) where D(s) is diagonal with the product of all denominators
# of column i in D[i,i]
# R[i,j] is num(s[i,j]) times all denominators in den(s[:,j]) except den(s[i,j])
function _tf2rfd(s::LTISystems.TransferFunction)
  mat = s.mat
  n,m = size(mat)
  dp  = fill(zero(den(mat[1])), m)
  R   = map(num, mat)
  map(num,s.mat)
  for i in 1:m
    dp[i]   = Base.reduce(*, one(den(s.mat[1,i])), map(den,s.mat[:,i]))
    for j in 1:n
      for k in 1:n
        k == j && continue
        R[j,i] *= den(s.mat[k,i])
      end
    end
  end
  R = PolyMatrix(R)
  D = PolyMatrix(diagm(dp))
  R, D
end

# s = inv(D(s))*R(s) where D(s) is diagonal with the product of all denominators
# of column i in D[i,i]
# R[i,j] is num(s[i,j]) times all denominators in den(s[j,:]) except den(s[i,j])
function _tf2lfd(s::LTISystems.TransferFunction)
  mat = s.mat
  n,m = size(mat)
  dp  = fill(zero(den(mat[1])), n)
  R   = map(num, mat)
  map(num,s.mat)
  for i in 1:n
    dp[i] = Base.reduce(*, one(den(s.mat[i,1])), map(den,s.mat[i,:]))
    for j in 1:m
      for k in 1:m
        k == j && continue
        R[i,j] *= den(s.mat[i,k])
      end
    end
  end
  R = PolyMatrix(R)
  D = PolyMatrix(diagm(dp))
  R, D
end

lfd(s::TransferFunction{Val{:mimo},Val{:cont}}) = lfd(rfd(_tf2rfd(s)...))
lfd(s::TransferFunction{Val{:mimo},Val{:disc}}) = lfd(rfd(_tf2rfd(s)...), s.Ts)
lfd(s::TransferFunction{Val{:siso},Val{:cont}}) = ((R,D) = _tf2lfd(s); lfd(R[1], D[1]))
lfd(s::TransferFunction{Val{:siso},Val{:disc}}) = ((R,D) = _tf2lfd(s); lfd(R[1], D[1], s.Ts))

rfd(s::TransferFunction{Val{:mimo},Val{:cont}}) = rfd(lfd(_tf2lfd(s)...))
rfd(s::TransferFunction{Val{:mimo},Val{:disc}}) = rfd(lfd(_tf2lfd(s)...), s.Ts)
rfd(s::TransferFunction{Val{:siso},Val{:cont}}) = ((R,D) = _tf2rfd(s); rfd(R[1], D[1]))
rfd(s::TransferFunction{Val{:siso},Val{:disc}}) = ((R,D) = _tf2rfd(s); rfd(R[1], D[1], s.Ts))
