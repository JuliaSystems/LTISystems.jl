module SystemsBase

# Polynomials-related things
using Polynomials

# RationalFunctions-related things
using RationalFunctions
import RationalFunctions: poles

# PolynomialMatrices-related things
# using PolynomialMatrices

# Plotting recipes
using RecipesBase

# Compatibility among different Julia versions
using Compat
import Compat.view

# Import conversion and promotion functions for overloading
import Base: convert, promote_rule

# Import identities for overloading
import Base: one, zero

# Import inv and zeros
import Base: inv, zeros

# Indexing
import Base: getindex, setindex!, endof

# Import iteration interface functions
import Base: start, next, done, eltype, length, size

# Import printing functions
import Base: summary
@compat import Base.show

# Import mathematical operations for overloading
import Base: +, .+, -, .-, *, .*, /, ./, ==, !=, isapprox, transpose

# Export only the useful functions
export
  # Interfaces
  issiso,
  ismimo,
  iscontinuous,
  isdiscrete,
  samplingtime,
  numstates,
  numinputs,
  numoutputs,
  poles,
  zeros,
  tzeros,
  zpkdata,
  numvec,
  denvec,
  numpoly,
  denpoly,
  # Constructors
  mfd,
  ss,
  tf,
  # Methods
  series,
  parallel,
  feedback,
  rosenbrock,
  minreal,
  freqresp

# System types
include("types/system/ltisystem.jl")
include("types/system/rationaltf.jl")
include("types/system/statespace.jl")
# include("types/system/mfd.jl")

# Conversions
# include("conversions/mfd2ss.jl")
# include("conversions/ss2mfd.jl")

# Response types
# include("types/response/systemresponse.jl")
# include("types/response/boderesponse.jl")
# include("types/response/nyquistresponse.jl")

# Methods
# include("methods/bode.jl")
# include("methods/feedback.jl")
include("methods/freqresp.jl")
include("methods/minreal.jl")
# include("methods/nyquist.jl")
# include("methods/parallel.jl")
# include("methods/reduce.jl")
include("methods/rosenbrock.jl")
# include("methods/series.jl")

end # module
