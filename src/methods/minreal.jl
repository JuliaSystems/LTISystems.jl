"""
    minreal(A, B, C, D, tol::Float64 = zero(Float64)) -> Am, Bm, Cm, Dm

Returns the minimal realization `(Am, Bm, Cm, Dm)` of a given state-space
representation `(A, B, C, D)`.

    minreal(sys::LtiSystem, tol::Float64 = zero(Float64))

When applied to LTI systems, it returns their minimal realization.

The method internally uses `rosenbrock` to obtain the minimal realization of the
system in two steps:

1.  obtain the controllable subspace of `(A, B, C, D)`,
2.  obtain the observable subspace of the realization from Step 1.

The method passes the tolerance value `tol` to `rosenbrock`.

See also: `rosenbrock`

# Examples
```julia
julia> A = diagm([1,2,3])
3x3 Array{Int64,2}:
 1  0  0
 0  2  0
 0  0  3

julia> B = reshape([1,0,2],3,1)
3x1 Array{Int64,2}:
 1
 0
 2

julia> C = [1 0 0; 1 1 0]
2x3 Array{Int64,2}:
 1  0  0
 1  1  0

julia> D = zeros(2,1);

julia> Am, Bm, Cm, Dm = minreal(A,B,C,D);

julia> Am
1x1 Array{Float64,2}:
 1.0

julia> Bm
1x1 Array{Float64,2}:
 1.0

julia> Cm
2x1 Array{Float64,2}:
 1.0
 1.0

julia> Dm
2x1 Array{Float64,2}:
 0.0
 0.0
```
"""
function minreal{M1<:AbstractMatrix,M2<:AbstractMatrix,M3<:AbstractMatrix,
  M4<:AbstractMatrix}(A::M1, B::M2, C::M3, D::M4, tol::Float64 = zero(Float64))
  @assert size(A,1) == size(A,2) "minreal: A must be square"
  @assert size(A,1) == size(B,1) "minreal: A and B must have the same row size"
  @assert size(A,1) == size(C,2) "minreal: A and C must have the same column size"
  @assert size(B,2) == size(D,2) "minreal: B and D must have the same column size"
  @assert size(C,1) == size(D,1) "minreal: C and D must have the same row size"

  n     = size(A,1)
  Tc, c = rosenbrock(A, B, tol)
  ATemp = (Tc'*A*Tc)[n-c+1:n,n-c+1:n]
  CTemp = (C*Tc)[:,n-c+1:n]
  To, o = rosenbrock(ATemp.', CTemp.', tol)

  T     = view(Tc, :, n-c+1:n)*To

  Am    = (T'*A*T)[c-o+1:c,c-o+1:c]
  Bm    = (T'*B)[c-o+1:c,:]
  Cm    = (C*T)[:,c-o+1:c]

  return Am, Bm, Cm, D
end

function minreal(s::StateSpace{Val{:siso},Val{:cont}},
  tol::Float64 = zero(Float64))
  Am, Bm, Cm, Dm = minreal(s.A, s.B, s.C, s.D, tol)
  ss(Am, Bm, Cm, Dm[1])
end

function minreal(s::StateSpace{Val{:siso},Val{:disc}},
  tol::Float64 = zero(Float64))
  Am, Bm, Cm, Dm = minreal(s.A, s.B, s.C, s.D, tol)
  ss(Am, Bm, Cm, Dm[1], s.Ts)
end

minreal(s::StateSpace{Val{:mimo},Val{:cont}},
  tol::Float64 = zero(Float64)) = ss(minreal(s.A, s.B, s.C, s.D, tol)...)

minreal(s::StateSpace{Val{:mimo},Val{:disc}},
  tol::Float64 = zero(Float64)) = ss(minreal(s.A, s.B, s.C, s.D, tol)..., s.Ts)
