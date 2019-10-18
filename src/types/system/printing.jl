_dimensions(s::LtiSystem{Val{:siso}})              = ""
_dimensions(s::LtiSystem{Val{:mimo}})              = "$(numoutputs(s))×$(numinputs(s))"
_time(::LtiSystem{Val{T},Val{:cont}}) where {T}    = "continuous"
_time(::LtiSystem{Val{T},Val{:disc}}) where {T}    = "discrete"
_timevar(::LtiSystem{Val{T},Val{:cont}}) where {T} = "t"
_timevar(::LtiSystem{Val{T},Val{:disc}}) where {T} = "k"
_printsamplingtime(s::LtiSystem{Val{T},Val{:cont}}) where {T} = ""
_printsamplingtime(s::LtiSystem{Val{T},Val{:disc}}) where {T} = "+$(samplingtime(s))"

import Base: summary

# StateSpace printing
_printxupdate(s::StateSpace{Val{T},Val{:cont}}) where {T} = "ẋ(t)"
_printxupdate(s::StateSpace{Val{T},Val{:disc}}) where {T} = "x(k+1)"

summary(s::StateSpace{Val{:siso},Val{T}}) where {T} = "StateSpace"
summary(s::StateSpace{Val{:mimo},Val{T}}) where {T} = "$(_dimensions(s)) StateSpace"

# Compact representations
function _compact(stream, ::MIME"text/plain", s::StateSpace{Val{T},Val{S}}) where {T,S}
  var = ifelse(S == :cont, "s", "z")
  print(stream, "G($(var))")
end

function _compact(stream, ::MIME"text/latex", s::StateSpace{Val{T},Val{S}}) where {T,S}
  var = ifelse(S == :cont, "s", "z")
  print(stream, "\$")
  print(stream, "G($(var))")
  print(stream, "\$")
end

# TODO: Think about text/html

# Full representations
function _full(stream, m::MIME"text/plain", s::StateSpace{Val{T},Val{S}}) where {T,S}
  tvar = _timevar(s)
  println(stream, summary(s))
  println(stream, "$(_printxupdate(s)) = Ax($tvar) + Bu($tvar,x)")
  println(stream, "y($tvar) = Cx($tvar) + Du($tvar,x)")
  println(stream, "with $(numstates(s)) states in $(_time(s)) time.")
end

function _full(stream, m::MIME"text/latex", s::StateSpace{Val{T},Val{S}}) where {T,S}
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
summary(s::TransferFunction{Val{:siso},Val{T}}) where {T} = "TransferFunction"
summary(s::TransferFunction{Val{:mimo},Val{T}}) where {T} = "$(_dimensions(s)) TransferFunction"

# Compact representations
function _compact(stream, ::MIME"text/plain", s::TransferFunction{Val{T},Val{S}}) where {T,S}
  var = ifelse(S == :cont, "s", "z")
  print(stream, "n($(var))/d($(var))")
end

function _compact(stream, ::MIME"text/latex", s::TransferFunction{Val{T},Val{S}}) where {T,S}
  var = ifelse(S == :cont, "s", "z")
  print(stream, "\$")
  print(stream, "\\tfrac{\\mathrm{n}($(var))}{\\mathrm{d}($(var))}")
  print(stream, "\$")
end

# TODO: Think about text/html

# Full representations
function _full(stream, m::MIME"text/plain", s::TransferFunction{Val{T},Val{S}}) where {T,S}
  var  = ifelse(S == :cont, "s", "z")
  tvar = _timevar(s)
  println(stream, summary(s))
  println(stream, "y($tvar) = num($(var))/den($(var))u($tvar)")
  println(stream, "in $(_time(s)) time.")
end

function _full(stream, m::MIME"text/latex", s::TransferFunction{Val{T},Val{S}}) where {T,S}
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

_mfdtype(::MatrixFractionDescription{Val{S},Val{T},Val{L}}) where {S,T,L} = ifelse(L == :lfd, "Left", "Right")

summary(s::MatrixFractionDescription{Val{:siso},Val{T},Val{L}}) where {T,L} =
  "$(_mfdtype(s)) MatrixFractionDescription"
summary(s::MatrixFractionDescription{Val{:mimo},Val{T},Val{L}}) where {T,L} =
  "$(_dimensions(s)) $(_mfdtype(s)) MatrixFractionDescription"

# Compact representations
function _compact(stream, ::MIME"text/plain", s::MatrixFractionDescription{Val{T},Val{S},Val{:lfd}}) where {T,S}
  var = ifelse(S == :cont, "s", "z")
  print(stream, "d($(var))", "\\", "n($(var))")
end

function _compact(stream, ::MIME"text/plain", s::MatrixFractionDescription{Val{T},Val{S},Val{:rfd}}) where {T,S}
  var = ifelse(S == :cont, "s", "z")
  print(stream, "n($(var))/d($(var))")
end

function _compact(stream, ::MIME"text/latex", s::MatrixFractionDescription{Val{T},Val{S},Val{:lfd}}) where {T,S}
  var = ifelse(S == :cont, "s", "z")
  print(stream, "\$")
  print(stream, "d($(var))", "\\", "n($(var))")
  print(stream, "\$")
end

function _compact(stream, ::MIME"text/latex", s::MatrixFractionDescription{Val{T},Val{S},Val{:rfd}}) where {T,S}
  var = ifelse(S == :cont, "s", "z")
  print(stream, "\$")
  print(stream, "n($(var))/d($(var))")
  print(stream, "\$")
end

# TODO: Think about text/html

# Full representations
function _full(stream, m::MIME"text/plain", s::MatrixFractionDescription{Val{T},Val{S},Val{:lfd}}) where {T,S}
  var  = ifelse(S == :cont, "s", "z")
  tvar = _timevar(s)
  println(stream, summary(s))
  println(stream, "y($tvar) = den($var)", "\\", "num($var) u($tvar)")
  println(stream, "in $(_time(s)) time.")
end

function _full(stream, m::MIME"text/plain", s::MatrixFractionDescription{Val{T},Val{S},Val{:rfd}}) where {T,S}
  var  = ifelse(S == :cont, "s", "z")
  tvar = _timevar(s)
  println(stream, summary(s))
  println(stream, "y($tvar) = num($var)/den($var) u($tvar)")
  println(stream, "in $(_time(s)) time.")
end

function _full(stream, m::MIME"text/latex", s::MatrixFractionDescription{Val{T},Val{S},Val{:lfd}}) where {T,S}
  var  = ifelse(S == :cont, "s", "z")
  tvar = _timevar(s)
  println(stream, summary(s))
  print(stream, "\$\$")
  print(stream, "y($tvar) = den^{-1}($var)num($var) u($tvar)")
  print(stream, "\$\$")
  println(stream, "in $(_time(s)) time.")
end

function _full(stream, m::MIME"text/latex", s::MatrixFractionDescription{Val{T},Val{S},Val{:rfd}}) where {T,S}
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
