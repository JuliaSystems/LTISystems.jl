module ControlCore

# Import conversion and promotion functions for overloading
import Base: convert, promote_rule

# Import identities for overloading
import Base: one, zero

# Import inv and zeros
import Base: inv, zeros

# Import slicing functions
import Base: ndims, size, getindex

# Import iteration interface functions
import Base: start, next, done, eltype, length, eachindex, endof

# Import printing functions
import Base: showcompact, show, showall, summary

# Import mathematical operations for overloading
import Base: +, .+, -, .-, *, .*, /, ./, ==, !=, isapprox, transpose

# Export only the useful functions
export
  # Interfaces
  numstates,
  numinputs,
  numoutputs,
  poles,
  # zeros = Base.zeros,
  tzeros,
  getmatrix,
  numvec,
  denvec,
  numpoly,
  denpoly,
  isdiscrete,
  samplingtime,
  zpkdata,
  # Constructors
  tf,
  zpk,
  ss,
  mimo,
  # Methods
  series,
  parallel,
  feedback,
  rosenbrock,
  minreal

using Polynomials
using Compat
import Compat.view

include("types/ltisystem.jl")
include("types/rationaltf.jl")
include("types/zeropolegain.jl")
include("types/statespace.jl")
include("types/generalmimo.jl")
include("display.jl")
include("interconnections/series.jl")
include("interconnections/parallel.jl")
include("interconnections/feedback.jl")
include("methods/rosenbrock.jl")
include("methods/minreal.jl")
include("methods/reduce.jl")

end # module
