immutable BodeResponse{T,S,U} <: SystemResponse
  freqs::T  # rad/sec
  mag::S    # abs, i.e., |G|
  phase::U  # radians

  @compat function (::Type{BodeResponse}){T<:Real,S<:Real,U<:Real}(
    freqs::AbstractVector{T}, mag::AbstractArray{S, 3}, phase::AbstractArray{U, 3})
    @assert size(mag) == size(phase) "BodeResponse: mag and phase must have same dimensions"
    @assert length(freqs) == size(mag,3) "BodeResponse: size(mag,3) ≠ length(freqs)"
    # TODO: Do we need an assertion for max ≥ 0 and freqs ≥ 0, too?

    new{typeof(freqs),typeof(mag),typeof(phase)}(freqs, mag, phase)
  end
end

# Plot everything
@recipe f(b::BodeResponse                     ) = (b, 1:size(b.mag,1), 1:size(b.mag,2))
@recipe f(b::BodeResponse, ::Colon            ) = (b, 1:size(b.mag,1), 1:size(b.mag,2))
@recipe f(b::BodeResponse, ::Colon, ::Colon   ) = (b, 1:size(b.mag,1), 1:size(b.mag,2))

# Plot a single (sub)system
@recipe f(b::BodeResponse, row::Int, col::Int ) = (b, row:row, col:col)

# Plot a given output for some outputs
@recipe f(b::BodeResponse, row::Int, ::Colon              ) = (b, row:row, 1:size(b.mag,2))
@recipe f(b::BodeResponse, row::Int, cols::AbstractVector ) = (b, row:row, cols           )

# Plot some outputs for a given input
@recipe f(b::BodeResponse, ::Colon, col::Int              ) = (b, 1:size(b.mag,1), col:col)
@recipe f(b::BodeResponse, rows::AbstractVector, col::Int ) = (b, rows,            col:col)

# Plot all outputs for some inputs
@recipe f(b::BodeResponse, ::Colon, cols::AbstractVector  ) = (b, 1:size(b.mag,1), cols   )

# Plot some outputs for all inputs
@recipe f(b::BodeResponse, rows::AbstractVector, ::Colon  ) = (b, rows, 1:size(b.mag,2)   )

# Plot some outputs for some inputs
@recipe function f(b::BodeResponse, rows::AbstractVector, cols::AbstractVector)
  @assert 1 ≤ minimum(rows) ≤ maximum(rows) ≤ size(b.mag,1) "plot(b::BodeResponse, rows, cols): rows out of bounds"
  @assert 1 ≤ minimum(cols) ≤ maximum(cols) ≤ size(b.mag,2) "plot(b::BodeResponse, rows, cols): cols out of bounds"

  # Maybe we get NTuple and/or some higher dimentional Array for rows and cols
  nrows = length(rows)
  ncols = length(cols)

  # Create labels
  labels      =   ["\$G_{$(row),$(col)}\$" for row in rows, col in cols] |>
                    vec |> x->reshape(x, 1, length(x))

  # Define plotting rules
  layout      :=  (2,1)
  link        :=  :x

  plot_title  --> "Bode Plot"
  legend      --> :topright
  grid        --> true

  @series begin
    subplot   :=  1
    title     --> "Magnitude Response"
    label     --> labels

    xlabel    --> "\$\\omega\$ (rad/sec)"
    xscale    --> :log10
    ylabel    --> "\$\\vert G \\vert\$ (dB)"

    b.freqs, 20*log10(reshape(b.mag[rows, cols, :], nrows*ncols, length(b.freqs))')
  end

  @series begin
    subplot   :=  2
    title     --> "Phase Response"
    label     --> labels

    xlabel    --> "\$\\omega\$ (rad/sec)"
    xscale    --> :log10
    ylabel    --> "\$\\angle G\$ (deg)"

    b.freqs, 180/π*(reshape(b.phase[rows, cols, :], nrows*ncols, length(b.freqs))')
  end
end

# TODO: Think of how to (elegantly) introduce flexibility in defining some of the
#       preferences:
#       -  Magnitude in (-)'s,
#       -  Phase in rad's,
#       -  Frequency in Hz,
#       -  etc.

# TODO: Implement some file I/O functionality such as writedlm, readdlm, etc.
