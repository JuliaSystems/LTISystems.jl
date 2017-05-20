# Generates a state space system whose transfer matrix coincides
# with the one of a matrix fraction description.
# Based on a column- and row-permuted version of the controller-form
# realization (Kailath, 1980, Section 6.4.1).
# NOTE: Base this function on the SLICOT routine TC04AD.
# NOTE: Is it necessary to make a SISO version of this function?
function _mfd2ss{S,M1,M2}(mfd::MatrixFractionDescription{Val{:mimo},S,Val{:rfd},M1,M2})
  Den, Num  = colred(mfd.D, mfd.N)
  kden, Phc = high_col_deg_matrix(Den)
  knum      = col_degree(Num)

  # Check for properness of the MatrixFractionDescription (Kailath, 1980, Theorem 6.3-12)
  all(knum .≤ kden) || error("mfd2ss: MatrixFractionDescription is not proper")

  # Normalize the coefficient matrices of Den
  Phci = inv(Phc)
  for (k,c) in Den.coeffs
    Den.coeffs[k] = Phci*Den.coeffs[k]
  end
  n  = sum(kden)  # Degree of the state space realization
  nu = mfd.nu
  ny = mfd.ny

  # Initialize state space matrices
  A = zeros(n, n)
  B = zeros(n, nu)
  C = zeros(ny, n)
  D = zeros(ny, nu)

  # Construction of the A and B matrices
  start_i = 0
  for i = 1:nu

    start_j = 0
    for j = 1:nu

      # Add upper diagonal 1's to A
      if i == j
        for ii = 1:kden[i]-1
          A[start_i+ii,start_i+ii+1] = 1
        end
      end

      # Add entries of Den to A
      for l = 0:kden[j]-1
        A[start_i + kden[i],start_j+l+1] = -Den.coeffs[l][i,j]
      end

      # Add entries of 1's to B
      for l = 1:nu
        B[start_i + kden[i],l] = Phci[i,l]
      end

      start_j += kden[j]
    end
    start_i += kden[i]
  end

  # Construction of the C and D matrices
  for i = 1:ny

    start_j = 0
    for j = 1:nu

      # Add entries of Den to A
      for l = 0:min(kden[j]-1,knum[j])
        C[i,start_j+l+1] = Num.coeffs[l][i,j]
      end

      # Modifications in case column j is biproper
      if knum[j] == kden[j]
        for l = 0:kden[j]-1
          C[i,start_j+l+1] -= Num.coeffs[knum[j]][i,j]*Den.coeffs[l][i,j]
        end
        D[i,j] = Num.coeffs[knum[j]][i,j]*Phci[i,j]
      end

      start_j += kden[j]
    end
  end

  return A, B, C, D
end

function _mfd2ss{S,M1,M2}(mfd::MatrixFractionDescription{Val{:mimo},S,Val{:lfd},M1,M2})
  Den, Num  = rowred(mfd.D,mfd.N)
  kden, Phr = high_row_deg_matrix(Den)
  knum      = row_degree(Num)

  # Check for properness of the MatrixFractionDescription (Kailath, 1980, Theorem 6.3-12)
  all(knum .≤ kden) || error("mfd2ss: MatrixFractionDescription is not proper")

  # Normalize the coefficient matrices of Den
  Phri = inv(Phr)
  for (k,c) in Den.coeffs
    Den.coeffs[k] = Den.coeffs[k]*Phri
  end
  n  = sum(kden)  # Degree of the state space realization
  nu = mfd.nu
  ny = mfd.ny

  # Initialize state space matrices
  A = zeros(n,  n )
  B = zeros(n,  nu)
  C = zeros(ny, n )
  D = zeros(ny, nu)

  # Construction of the A and C matrices
  start_i = 0
  for i = 1:nu

    start_j = 0
    for j = 1:nu

      # Add upper diagonal 1's to A
      if i == j
        for ii = 1:kden[i]-1
          A[start_i+ii+1,start_i+ii] = 1
        end
      end

      # Add entries of Den to A
      for l = 0:kden[j]-1
        A[start_j+l+1,start_i + kden[i]] = -Den.coeffs[l][j,i]
      end

      # Add entries of 1's to C
      for l = 1:ny
        C[l, start_i + kden[i]] = Phri[l,i]
      end

      start_j += kden[j]
    end
    start_i += kden[i]
  end

  # Construction of the B and D matrices
  for i = 1:nu

    start_j = 0
    for j = 1:nu

      # Add entries of Den to A
      for l = 0:min(kden[j]-1,knum[j])
        B[start_j+l+1,i] = Num.coeffs[l][j,i]
      end

      # Modifications in case row j is biproper
      if knum[j] == kden[j]
        for l = 0:kden[j]-1
          B[start_j+l+1,i] -= Num.coeffs[knum[j]][j,i]*Den.coeffs[l][j,i]
        end
        D[j,i] = Num.coeffs[knum[j]][j,i]*Phri[j,i]
      end

      start_j += kden[j]
    end
  end

  return A, B, C, D
end

ss(mfd::MatrixFractionDescription{Val{:mimo},Val{:cont}}) = minreal(ss(_mfd2ss(mfd)...))
ss(mfd::MatrixFractionDescription{Val{:mimo},Val{:disc}}) = minreal(ss(_mfd2ss(mfd)..., mfd.Ts))

ss(mfd::MatrixFractionDescription{Val{:siso}})            = ss(tf(mfd))
