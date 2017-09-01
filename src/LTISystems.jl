__precompile__(true)

module LTISystems

# Compatibility among different Julia versions
using Compat

# OrdinaryDiffEq
using DiffEqBase, OrdinaryDiffEq
@compat const DEDataVector{T} = isconst(DiffEqBase, :DEDataVector) ?
                                DiffEqBase.DEDataVector{T} :
                                DiffEqBase.DEDataArray{T}

# PolynomialMatrices-related things
using PolynomialMatrices

# Polynomials-related things
using Polynomials

# RationalFunctions-related things
using RationalFunctions
import RationalFunctions: poles

# Plotting recipes
using RecipesBase

# Import conversion and promotion functions for overloading
import Base: convert, promote_rule

# Import identities for overloading
import Base: one, zero

# Import num, den for getting numerator and denominator polynomials
import Base: num, den

# Import inv and zeros
import Base: inv, zeros

# Indexing
import Base: getindex, setindex!, endof

# Import iteration interface functions
import Base: start, next, done, eltype, length, size

# Import printing functions
import Base: summary

# Import step
import Base: step

# Import mathematical operations for overloading
import Base: +, .+, -, .-, *, .*, /, ./, ==, !=, isapprox, transpose

# Export only the useful functions
export
  # Interfaces
  issiso,
  ismimo,
  siso,
  mimo,
  iscontinuous,
  isdiscrete,
  isproper,
  samplingtime,
  numstates,
  numinputs,
  numoutputs,
  poles,
  zeros,
  tzeros,
  zpkdata,
  num,
  den,
  # Constructors
  mfd,
  lfd,
  rfd,
  ss,
  tf,
  #
  islfd,
  isrfd,
  # Methods
  bode,
  feedback,
  freqresp,
  minreal,
  nyquist,
  parallel,
  pzmap,
  # reduce,
  rosenbrock,
  series,
  simulate, step, ramp,
  unwrap, unwrap!,
  Signals

# System types
include("types/system/ltisystem.jl")
include("types/system/transferfunction.jl")
include("types/system/statespace.jl")
include("types/system/mfd.jl")
include("types/system/printing.jl")

# Signal types
module Signals

using Compat

import Base: convert, promote_rule, +, -, *, /

export  discontinuities,
        PRBS,
        Ramp,
        Sinusoid,
        Square,
        Step,
        Triangle

include("types/signals/abstractsignal.jl")
include("types/signals/prbs.jl")
include("types/signals/ramp.jl")
include("types/signals/sinusoid.jl")
include("types/signals/square.jl")
include("types/signals/step.jl")
include("types/signals/sumofsignals.jl")
include("types/signals/triangle.jl")

end
using .Signals

# Conversions
include("conversions/mfd2ss.jl")
include("conversions/mfd2tf.jl")
include("conversions/ss2mfd.jl")
#include("conversions/ss2tf.jl")
include("conversions/tf2mfd.jl")
include("conversions/tf2ss.jl")

# Response types
include("types/response/systemresponse.jl")
include("types/response/boderesponse.jl")
include("types/response/nyquistresponse.jl")
include("types/response/timeresponse.jl")

# Methods
include("methods/feedback.jl")
include("methods/freqresp.jl")
include("methods/minreal.jl")
include("methods/parallel.jl")
include("methods/pzmap.jl")
include("methods/reduce.jl")
include("methods/rosenbrock.jl")
include("methods/series.jl")
include("methods/simulate.jl")
include("methods/unwrap.jl")

end # module
