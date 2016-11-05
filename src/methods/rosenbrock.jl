"""
    rosenbrock(A, B, tol::Float64 = 0.) -> T, c

Stabilized [version][1] of [Rosenbrock][2]'s method to obtain a state
transformation for the state-space model of a system defined by `(A,B,C,D)`
matrices.

If the method is given `A` and `B` matrices, `T` is the unitary state
transformation needed to put the state-space model into its controllable
realization. Then, `c` represents the dimension of the controllable subspace of
the system. If, however, the method is given `A'` and `C'` matrices (i.e., the
dual system), `T` is the similarity transformation needed to obtain the system's
observable realization, with `c` representing the dimension of the observable
subspace.

When a non-positive tolerance value `tol` is given, it will default to
`10*max(size(B,1,2)...)*norm(B,1)*eps(Float64)`.

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

julia> C = zeros(1,3); D = zeros(1,1);

julia> Tc, c = rosenbrock(A, B);

julia> # Controllable realization

julia> Ac = Tc'*A*Tc
3x3 Array{Float64,2}:
 2.0  0.0  0.0
 0.0  1.4  0.8
 0.0  0.8  2.6

julia> Bc = Tc'*B
3x1 Array{Float64,2}:
 0.0
 4.44089e-16
 2.23607

julia> Cc = C*Tc; Dc = D;
```

# References

-  [1]: P.M. Van Dooren, "The generalized eigenstructure problem on linear system
        theory," IEEE Transactions on Automatic Control, vol. 26, no. 1, pp.
        111-129, Feb. 1981.
-  [2]: H.H. Rosenbrock, State-Space and Multivariable Theory. London: Thomas
        Nelson and Sons, 1970.
"""
function rosenbrock{M1<:AbstractMatrix,M2<:AbstractMatrix}(A::M1, B::M2,
  tol::Float64 = zero(Float64))
  n1, n2  = size(A,1,2)
  n3, n4  = size(B,1,2)
  @assert n1 == n2 "rosenbrock: A must be square"
  @assert n1 == n3 "rosenbrock: A and B must have same row sizes"
  c::Int              = 0
  T::Matrix{Float64}  = eye(n1)
  τ₁, τ₂  = (n1,1)
  ρ₁, ρ₂  = (0,1)
  Atemp   = A
  Btemp   = B
  tol     = max(tol, 10*max(n3,n4)*norm(B,1)*eps(Float64) + eps(Float64))
  while ρ₂ != 0 && τ₂ != 0
    n       = size(Btemp,1)
    svdobj  = svdfact(Btemp, thin = false)
    U       = flipdim(svdobj.U, 2)
    ρ₂      = sum(svdobj.S .>= tol)
    τ₂      = n - ρ₂

    if ρ₂ != 0 && τ₂ != 0
      Atemp   = U'*Atemp*U
      Btemp   = view(Atemp, 1:τ₂, τ₂+1:n)
      Atemp   = view(Atemp, 1:τ₂, 1:τ₂)
      T[:]    = T*vcat(hcat(U, zeros(size(U, 1),c)),
                        hcat(zeros(c, size(U, 2)), eye(c)))

      c       = c + ρ₂
      ρ₁, τ₁  = (ρ₂, τ₂)
    end
  end

  c = (ρ₂ == 0) ? size(A,1) - τ₁ : size(A,1)
  return T, c
end
