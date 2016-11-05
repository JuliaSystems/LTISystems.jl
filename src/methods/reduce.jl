function reduce{M1<:AbstractMatrix, M2<:AbstractMatrix, M3<:AbstractMatrix,
  M4<:AbstractMatrix}(A::M1, B::M2, C::M3, D::M4, tol::Float64 = zero(Float64))
  @assert size(A,1) == size(A,2) "reduce: A must be square"
  @assert size(A,1) == size(B,1) "reduce: A and B must have the same row size"
  @assert size(A,1) == size(C,2) "reduce: A and C must have the same column size"
  @assert size(B,2) == size(D,2) "reduce: B and D must have the same column size"
  @assert size(C,1) == size(D,1) "reduce: C and D must have the same row size"

  n1, n2  = size(B,1,2)
  n3, n4  = size(C,1,2)

  ν₁, ν₂  = n1, n1
  σ₁, σ₂  = n3, n3
  δ       = 0
  μ       = n3
  ρ       = n1
  τ       = n3

  Atemp   = convert(Matrix{Float64}, A)
  Btemp   = convert(Matrix{Float64}, B)
  Ctemp   = convert(Matrix{Float64}, C)
  Dtemp   = convert(Matrix{Float64}, D)

  temp1   = [Ctemp Dtemp]
  C̄, D̄    = Ctemp, Dtemp

  tol     = max(tol, 10*max(n3,n2+n4)*norm(temp1,1)*eps(Float64) + eps(Float64))

  while τ ≠ 0 && ρ ≠ 0 && ν₂ ≠ 0
    n3, n4  = size(Ctemp,1,2)
    svdobj  = svdfact(Dtemp, thin = false)
    σ₂      = sum(svdobj.S .>= tol)
    τ       = n3 - σ₂

    U       = svdobj.U
    temp2   = U'*temp1
    C̄       = temp2[1:σ₂, 1:n4]
    D̄       = temp2[1:σ₂, n4+1:end]

    if τ ≠ 0
      svdobj  = svdfact(view(temp2, σ₂+1:n3, 1:n4), thin = false)
      ρ       = sum(svdobj.S .>= tol)
      ν₂      = n4 - ρ

      if ρ ≠ 0 && ν₂ ≠ 0
        μ     = ρ + σ₂
        δ     += ρ
        V     = flipdim(svdobj.Vt', 2)
        temp3 = [V'*Atemp*V ; C̄*V   ]
        temp4 = [V'*Btemp   ; D̄ ]

        Atemp = temp3[1:ν₂, 1:ν₂]
        Ctemp = temp3[ν₂+1:end, 1:ν₂]
        Btemp = temp4[1:ν₂, :]
        Dtemp = temp4[ν₂+1:end,:]

        temp1 = [Ctemp Dtemp]

        ν₁, σ₁ = ν₂, σ₂
      end
    end
  end

  if τ == 0 || ρ == 0 # exit 1
    return (Atemp::Matrix{Float64}, Btemp::Matrix{Float64},
      C̄::Matrix{Float64}, D̄::Matrix{Float64},
      n2::Int, ν₁::Int, σ₂::Int)
  else                # exit 2
    return (Atemp::Matrix{Float64}, Btemp::Matrix{Float64},
      Ctemp::Matrix{Float64}, Dtemp::Matrix{Float64}, n2::Int, 0::Int, σ₂::Int)
  end
end
