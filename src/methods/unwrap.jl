"""
  unwrap!(Θ₀[, dim = ndims(Θ)]; tol = π) -> Θ₀
  unwrap(Θ₀[, dim = ndims(Θ)]; tol = π) -> Θ₁

Unwrap the `Array` of phase values `Θ₀` (in radians) along the dimension `dim`
by changing the jump values larger than `tol` to their `2π` complements.

Internally, `unwrap` sets the tolerance to `max(tol, π)` as discountinuities smaller
than `π` should not be unwrapped. `unwrap!` modifies the original input `Θ₀` in
place and returns the modified version, whereas `unwrap` does not modify `Θ₀`.

**See also:* `freqresp`, `bode`, `nyquist`.
"""
function unwrap!{T<:AbstractFloat}(Θ::AbstractArray{T}, dim::Integer = ndims(Θ);
  tol::Real = π)
  size(Θ, dim) < 2 && return Θ
  settol  = max(float(tol), π)
  setrng  = 2settol
  circ    = 2π
  for idx in 2:size(Θ, dim)
    changeddims   = ntuple(n->(n == dim ? (idx:idx) : 1:size(Θ, n)), dim)
    jump          = slicedim(Θ, dim, idx:idx) - slicedim(Θ, dim, idx-1:idx-1)
    map!(j->floor((j+settol)/setrng)*circ, jump)
    Θ[changeddims...] -= jump
  end
  return Θ
end

unwrap{T<:Real}(Θ::AbstractArray{T}, dim::Integer = ndims(Θ);
  tol::Real = π) = unwrap!(map(float,Θ), dim; tol = tol)
