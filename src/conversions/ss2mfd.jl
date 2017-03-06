# Generates a right MFD whose transfer function coincides with that
# of an observable state-space model.
function _ss2rfd{T,S}(sys::StateSpace{Val{T},Val{S}}, var::Symbol)
  n     = numstates(sys)
  nu    = numinputs(sys)
  ny    = numoutputs(sys)

  # Construct system matrix (Rosenbrock, 1970)
  s     = Poly([0, 1],var)
  c     = vcat(hcat(-sys.A, sys.B), hcat(-eye(A), zeros(B)))
  P     = PolyMatrix(c, (n,n+nu), Val{_pmvaltype(sys)}) # s*eye(n)-sys.A sys.B

  P2, Q = rtriang(P)

  # NOTE: Missing to check if P2[ny,:] is non-zero (this is equivalent to controllability)

  # Write (sI-A)^{-1} B as a left MFD
  Nr = Q[1:n,n+1:end]
  Dr = Q[n+1:end,n+1:end]

  # Include C and D in the final right MFD
  return PolyMatrix(sys.D*Dr - sys.C*Nr), PolyMatrix(Dr)
end
rfd(sys::StateSpace{Val{:siso},Val{:cont}})  =
  rfd(map(x -> x[1], _ss2rfd(sys,:s))...)
  rfd(sys::StateSpace{Val{:siso},Val{:disc}}) =
    rfd(map(x -> x[1], _ss2rfd(sys,:z))...,samplingtime(sys))
rfd(sys::StateSpace{Val{:mimo},Val{:cont}})  =
  rfd(_ss2rfd(sys,:s)...)
rfd(sys::StateSpace{Val{:mimo},Val{:disc}}) =
  rfd(_ss2rfd(sys,:z)...,samplingtime(sys))

_pmvaltype{T,S}(s::LtiSystem{Val{T},Val{S}}) = S
#_pmvaltype{T}(s::LtiSystem{Val{T},Val{:cont}}) = Val{:s}
#_pmvaltype{T}(s::LtiSystem{Val{T},Val{:disc}}) = Val{:z}
