# Generates a right MFD whose transfer function coincides with that
# of an observable state-space model.
function _ss2rfd(sys::StateSpace, var::Symbol)
  n     = numstates(sys)
  nu    = numinputs(sys)
  ny    = numoutputs(sys)

  # Construct system matrix (Rosenbrock, 1970)
  s     = Poly([0, 1],var)
  P     = PolyMatrix([s*eye(n)-sys.A sys.B])

  P2, Q = rtriang(P)

  # NOTE: Missing to check if P2[ny,:] is non-zero (this is equivalent to controllability)

  # Write (sI-A)^{-1} B as a left MFD
  Nr = Q[1:n,n+1:end]
  Dr = Q[n+1:end,n+1:end]

  # Include C and D in the final right MFD
  return PolyMatrix(sys.D*Dr - sys.C*Nr), PolyMatrix(Dr)
end
rfd(sys::StateSpace{Siso{true},Continuous{true}})  =
  rfd(map(x -> x[1], _ss2rfd(sys,:s))...)
  rfd(sys::StateSpace{Siso{true},Continuous{false}}) =
    rfd(map(x -> x[1], _ss2rfd(sys,:z))...,samplingtime(sys))
rfd(sys::StateSpace{Siso{false},Continuous{true}})  =
  rfd(_ss2rfd(sys,:s)...)
rfd(sys::StateSpace{Siso{false},Continuous{false}}) =
  rfd(_ss2rfd(sys,:z)...,samplingtime(sys))
