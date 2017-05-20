# Generates a right MatrixFractionDescription whose transfer function coincides with that
# of an observable state-space model.
function _ss2rfd{T,S}(sys::StateSpace{Val{T},Val{S}}, var::Symbol)
  n     = numstates(sys)
  nu    = numinputs(sys)
  ny    = numoutputs(sys)

  # Construct system matrix (Rosenbrock, 1970)
  Dₗ    = PolyMatrix(vcat(sys.A, -eye(sys.A)), (n,n), Val{_pmvaltype(sys)})
  Nₗ    = PolyMatrix(sys.B, (n,nu), Val{_pmvaltype(sys)})

  # Write (sI-A)^{-1} B as a left MatrixFractionDescription and convert to right MatrixFractionDescription
  L, U = ltriang(hcat(-Dₗ,Nₗ))

  # NOTE: Missing to check if P2[ny,:] is non-zero (this is equivalent to controllability)
  if norm(L[ny,:]) < 1e-12
    warn("rfd: sys might not be controllable")
  end

  Nr = U[1:n,n+1:end]
  Dr = U[n+1:end,n+1:end]

  # Include C and D in the final right MatrixFractionDescription
  return sys.D*Dr - sys.C*Nr, Dr
end

rfd(sys::StateSpace{Val{:siso},Val{:cont}})  =
  rfd(map(x -> x[1], _ss2rfd(sys,:s))...)
  rfd(sys::StateSpace{Val{:siso},Val{:disc}}) =
    rfd(map(x -> x[1], _ss2rfd(sys,:z))...,samplingtime(sys))
rfd(sys::StateSpace{Val{:mimo},Val{:cont}})  =
  rfd(_ss2rfd(sys,:s)...)
rfd(sys::StateSpace{Val{:mimo},Val{:disc}}) =
  rfd(_ss2rfd(sys,:z)...,samplingtime(sys))

function _ss2lfd{T,S}(sys::StateSpace{Val{T},Val{S}}, var::Symbol)
  n     = numstates(sys)
  nu    = numinputs(sys)
  ny    = numoutputs(sys)

  # Construct system matrix (Rosenbrock, 1970)
  Dᵣ    = PolyMatrix(vcat(-sys.A, eye(sys.A)), (n,n), Val{_pmvaltype(sys)})
  Nᵣ    = PolyMatrix(sys.C, (ny,n), Val{_pmvaltype(sys)})
  # Write C(sI-A)^{-1} as a right MatrixFractionDescription and convert to left MatrixFractionDescription
  R,U   = rtriang(vcat(-Dᵣ, Nᵣ))

  # NOTE: Missing to check if R[ny,:] is non-zero (this is equivalent to controllability)
  if norm(R[:,nu]) < 1e-12
    warn("lfd: system might not be observable")
  end
  Nₗ = U[n+1:end,1:n]
  Dₗ = U[n+1:end,n+1:end]

  # Include B and D in the final right MatrixFractionDescription
  return Dₗ*sys.D - Nₗ*sys.B, Dₗ
end

lfd(sys::StateSpace{Val{:siso},Val{:cont}})  =
  lfd(map(x -> x[1], _ss2lfd(sys,:s))...)
lfd(sys::StateSpace{Val{:siso},Val{:disc}}) =
  lfd(map(x -> x[1], _ss2lfd(sys,:z))..., samplingtime(sys))
lfd(sys::StateSpace{Val{:mimo},Val{:cont}})  =
  lfd(_ss2lfd(sys,:s)...)
lfd(sys::StateSpace{Val{:mimo},Val{:disc}}) =
  lfd(_ss2lfd(sys,:z)...,samplingtime(sys))

_pmvaltype{T}(s::LtiSystem{Val{T},Val{:cont}}) = :s
_pmvaltype{T}(s::LtiSystem{Val{T},Val{:disc}}) = :z
