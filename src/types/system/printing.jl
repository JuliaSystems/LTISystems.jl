_dimensions(s::LtiSystem{Val{:siso}})     = ""
_dimensions(s::LtiSystem{Val{:mimo}})     = "$(numoutputs(s))×$(numinputs(s))"
_time{T}(::LtiSystem{Val{T},Val{:cont}}) = "continuous"
_time{T}(::LtiSystem{Val{T},Val{:disc}}) = "discrete"
_timevar{T}(::LtiSystem{Val{T},Val{:cont}}) = "t"
_timevar{T}(::LtiSystem{Val{T},Val{:disc}}) = "k"
_printsamplingtime{T}(s::LtiSystem{Val{T},Val{:cont}}) = ""
_printsamplingtime{T}(s::LtiSystem{Val{T},Val{:disc}}) = "+$(samplingtime(s))"

import Base: summary

# StateSpace printing
_printxupdate{T}(s::StateSpace{Val{T},Val{:cont}})  = "ẋ(t)"
_printxupdate{T}(s::StateSpace{Val{T},Val{:disc}})  = "x(k+1)"

summary{T}(s::StateSpace{Val{:siso},Val{T}})  = "StateSpace"
summary{T}(s::StateSpace{Val{:mimo},Val{T}})  = "$(_dimensions(s)) StateSpace"

# Compact representations
function _compact{T,S}(stream, ::MIME"text/plain", s::StateSpace{Val{T},Val{S}})
  var = ifelse(S == :cont, "s", "z")
  print(stream, "G($(var))")
end

function _compact{T,S}(stream, ::MIME"text/latex", s::StateSpace{Val{T},Val{S}})
  var = ifelse(S == :cont, "s", "z")
  print(stream, "\$")
  print(stream, "G($(var))")
  print(stream, "\$")
end

# TODO: Think about text/html

# Full representations
function _full{T,S}(stream, m::MIME"text/plain", s::StateSpace{Val{T},Val{S}})
  tvar = _timevar(s)
  println(stream, summary(s))
  println(stream, "$(_printxupdate(s)) = Ax($tvar) + Bu($tvar,x)")
  println(stream, "y($tvar) = Cx($tvar) + Du($tvar,x)")
  println(stream, "with $(numstates(s)) states in $(_time(s)) time.")
end

function _full{T,S}(stream, m::MIME"text/latex", s::StateSpace{Val{T},Val{S}})
  tvar = _timevar(s)
  println(stream, summary(s))
  print(stream, "\\begin{align*}")
  print(stream, "$(_printxupdate(s)) &= Ax($tvar) + Bu($tvar,x)\\\\")
  print(stream, "y($tvar)            &= Cx($tvar) + Du($tvar,x)")
  print(stream, "\\end{align*}")
  println(stream, "with $(numstates(s)) states in $(_time(s)) time.")
end

# TODO: Think about text/html

# `show` function
@compat Base.show(stream::IO, s::StateSpace)                          =
  Base.show(stream, MIME("text/plain"), s)
@compat Base.show(stream::IO, mime::MIME"text/plain", s::StateSpace)  =
  get(stream, :compact, false) ? _compact(stream, mime, s) : _full(stream, mime, s)
@compat Base.show(stream::IO, mime::MIME"text/latex", s::StateSpace)  =
  get(stream, :compact, false) ? _compact(stream, mime, s) : _full(stream, mime, s)

## TransferFunction printing
summary{T}(s::TransferFunction{Val{:siso},Val{T}})  = "TransferFunction"
summary{T}(s::TransferFunction{Val{:mimo},Val{T}})  = "$(_dimensions(s)) TransferFunction"

# Compact representations
function _compact{T,S}(stream, ::MIME"text/plain", s::TransferFunction{Val{T},Val{S}})
  var = ifelse(S == :cont, "s", "z")
  print(stream, "n($(var))/d($(var))")
end

function _compact{T,S}(stream, ::MIME"text/latex", s::TransferFunction{Val{T},Val{S}})
  var = ifelse(S == :cont, "s", "z")
  print(stream, "\$")
  print(stream, "\\tfrac{\\mathrm{n}($(var))}{\\mathrm{d}($(var))}")
  print(stream, "\$")
end

# TODO: Think about text/html

# Full representations
function _full{T,S}(stream, m::MIME"text/plain", s::TransferFunction{Val{T},Val{S}})
  var  = ifelse(S == :cont, "s", "z")
  tvar = _timevar(s)
  println(stream, summary(s))
  println(stream, "y($tvar) = num($(var))/den($(var))u($tvar)")
  println(stream, "in $(_time(s)) time.")
end

function _full{T,S}(stream, m::MIME"text/latex", s::TransferFunction{Val{T},Val{S}})
  var  = ifelse(S == :cont, "s", "z")
  tvar = _timevar(s)
  println(stream, summary(s))
  print(stream, "\$\$")
  print(stream, "y($tvar) = \\tfrac{\\mathrm{num}($(var))}{\\mathrm{den}($(var))} u($tvar)")
  print(stream, "\$\$")
  println(stream, "in $(_time(s)) time.")
end

# TODO: Think about text/html

# `show` function
@compat Base.show(stream::IO, s::TransferFunction)                          =
  Base.show(stream, MIME("text/plain"), s)
@compat Base.show(stream::IO, mime::MIME"text/plain", s::TransferFunction)  =
  get(stream, :compact, false) ? _compact(stream, mime, s) : _full(stream, mime, s)
@compat Base.show(stream::IO, mime::MIME"text/latex", s::TransferFunction)  =
  get(stream, :compact, false) ? _compact(stream, mime, s) : _full(stream, mime, s)

## MatrixFractionDescription printing

_mfdtype{S,T,L}(::MatrixFractionDescription{Val{S},Val{T},Val{L}}) = ifelse(L == :lfd, "Left", "Right")

summary{T,L}(s::MatrixFractionDescription{Val{:siso},Val{T},Val{L}}) =
  "$(_mfdtype(s)) MatrixFractionDescription"
summary{T,L}(s::MatrixFractionDescription{Val{:mimo},Val{T},Val{L}}) =
  "$(_dimensions(s)) $(_mfdtype(s)) MatrixFractionDescription"

# Compact representations
function _compact{T,S}(stream, ::MIME"text/plain", s::MatrixFractionDescription{Val{T},Val{S},Val{:lfd}})
  var = ifelse(S == :cont, "s", "z")
  print(stream, "d($(var))", "\\", "n($(var))")
end

function _compact{T,S}(stream, ::MIME"text/plain", s::MatrixFractionDescription{Val{T},Val{S},Val{:rfd}})
  var = ifelse(S == :cont, "s", "z")
  print(stream, "n($(var))/d($(var))")
end

function _compact{T,S}(stream, ::MIME"text/latex", s::MatrixFractionDescription{Val{T},Val{S},Val{:lfd}})
  var = ifelse(S == :cont, "s", "z")
  print(stream, "\$")
  print(stream, "d($(var))", "\\", "n($(var))")
  print(stream, "\$")
end

function _compact{T,S}(stream, ::MIME"text/latex", s::MatrixFractionDescription{Val{T},Val{S},Val{:rfd}})
  var = ifelse(S == :cont, "s", "z")
  print(stream, "\$")
  print(stream, "n($(var))/d($(var))")
  print(stream, "\$")
end

# TODO: Think about text/html

# Full representations
function _full{T,S}(stream, m::MIME"text/plain", s::MatrixFractionDescription{Val{T},Val{S},Val{:lfd}})
  var  = ifelse(S == :cont, "s", "z")
  tvar = _timevar(s)
  println(stream, summary(s))
  println(stream, "y($tvar) = den($var)", "\\", "num($var) u($tvar)")
  println(stream, "in $(_time(s)) time.")
end

function _full{T,S}(stream, m::MIME"text/plain", s::MatrixFractionDescription{Val{T},Val{S},Val{:rfd}})
  var  = ifelse(S == :cont, "s", "z")
  tvar = _timevar(s)
  println(stream, summary(s))
  println(stream, "y($tvar) = num($var)/den($var) u($tvar)")
  println(stream, "in $(_time(s)) time.")
end

function _full{T,S}(stream, m::MIME"text/latex", s::MatrixFractionDescription{Val{T},Val{S},Val{:lfd}})
  var  = ifelse(S == :cont, "s", "z")
  tvar = _timevar(s)
  println(stream, summary(s))
  print(stream, "\$\$")
  print(stream, "y($tvar) = den^{-1}($var)num($var) u($tvar)")
  print(stream, "\$\$")
  println(stream, "in $(_time(s)) time.")
end

function _full{T,S}(stream, m::MIME"text/latex", s::MatrixFractionDescription{Val{T},Val{S},Val{:rfd}})
  var  = ifelse(S == :cont, "s", "z")
  tvar = _timevar(s)
  println(stream, summary(s))
  print(stream, "\$\$")
  print(stream, "y($tvar) = num($var)den^{-1}($var) u($tvar)")
  print(stream, "\$\$")
  println(stream, "in $(_time(s)) time.")
end

# TODO: Think about text/html

# `show` function
@compat Base.show(stream::IO, s::MatrixFractionDescription)                          =
  Base.show(stream, MIME("text/plain"), s)
@compat Base.show(stream::IO, mime::MIME"text/plain", s::MatrixFractionDescription)  =
  get(stream, :compact, false) ? _compact(stream, mime, s) : _full(stream, mime, s)
@compat Base.show(stream::IO, mime::MIME"text/latex", s::MatrixFractionDescription)  =
  get(stream, :compact, false) ? _compact(stream, mime, s) : _full(stream, mime, s)
