immutable NyquistResponse{T,S,U} <: SystemResponse
  freqs::T  # rad/sec
  real::S    # abs, i.e., |G|
  imag::U  # radians

  @compat function (::Type{NyquistResponse}){T<:Real,S<:Real,U<:Real}(
    freqs::AbstractVector{T}, real::AbstractArray{S, 3}, imag::AbstractArray{U, 3})
    @assert size(real) == size(imag) "NyquistResponse: real and imag must have same dimensions"
    @assert length(freqs) == size(real,3) "NyquistResponse: size(real,3) ≠ length(freqs)"

    new{typeof(freqs),typeof(real),typeof(imag)}(freqs, real, imag)
  end
end

# Plot everything
@recipe f(n::NyquistResponse                     ) = (n, 1:size(n.real,1), 1:size(n.real,2))
@recipe f(n::NyquistResponse, ::Colon            ) = (n, 1:size(n.real,1), 1:size(n.real,2))
@recipe f(n::NyquistResponse, ::Colon, ::Colon   ) = (n, 1:size(n.real,1), 1:size(n.real,2))

# Plot a single (sub)system
@recipe f(n::NyquistResponse, row::Int, col::Int ) = (n, row:row, col:col)

# Plot a given output for some outputs
@recipe f(n::NyquistResponse, row::Int, ::Colon              ) = (n, row:row, 1:size(n.real,2))
@recipe f(n::NyquistResponse, row::Int, cols::AbstractVector ) = (n, row:row, cols           )

# Plot some outputs for a given input
@recipe f(n::NyquistResponse, ::Colon, col::Int              ) = (n, 1:size(n.real,1), col:col)
@recipe f(n::NyquistResponse, rows::AbstractVector, col::Int ) = (n, rows,            col:col)

# Plot all outputs for some inputs
@recipe f(n::NyquistResponse, ::Colon, cols::AbstractVector  ) = (n, 1:size(n.real,1), cols   )

# Plot some outputs for all inputs
@recipe f(n::NyquistResponse, rows::AbstractVector, ::Colon  ) = (n, rows, 1:size(n.real,2)   )

# Plot some outputs for some inputs
@recipe function f(n::NyquistResponse, rows::AbstractVector, cols::AbstractVector)
  @assert 1 ≤ minimum(rows) ≤ maximum(rows) ≤ size(n.real,1) "plot(n::NyquistResponse, rows, cols): rows out of bounds"
  @assert 1 ≤ minimum(cols) ≤ maximum(cols) ≤ size(n.real,2) "plot(n::NyquistResponse, rows, cols): cols out of bounds"

  # Maybe we get NTuple and/or some higher dimentional Array for rows and cols
  nrows = length(rows)
  ncols = length(cols)

  # Create labels
  labels      =   ["\$G_{$(row),$(col)}\$" for row in rows, col in cols] |>
                    vec |> x->reshape(x, 1, length(x))

  # Define plotting rules
  layout      :=  (1,1)

  title       --> "Nyquist Plot"
  legend      --> :topright
  grid        --> true

  @series begin
    subplot   :=  1
    label     --> labels

    xlabel    --> "Real Axis (-)"
    ylabel    --> "Imaginary Axis (-)"

    reshape(n.real[rows, cols, :], nrows*ncols, length(n.freqs))',
      reshape(n.imag[rows, cols, :], nrows*ncols, length(n.freqs))'
  end
end

# TODO: Think of how to (elegantly) introduce flexibility in defining some of the
#       desired preferences, here, too.

# TODO: Implement some file I/O functionality such as writedlm, readdlm, etc.
