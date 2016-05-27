# TODO tidy upp the print functions

function printtransferfunction{T1}(io::IO, s::DSisoRational{T1})
  numstr = print_poly_reverse(s.num)
  denstr = print_poly_reverse(s.den)

  # Figure out the length of the separating line
  len_num = length(numstr)
  len_den = length(denstr)
  dashcount = max(len_num, len_den)

  # Center the numerator or denominator
  if len_num < dashcount
    numstr = "$(repeat(" ", div(dashcount - len_num, 2)))$numstr"
  else
    denstr = "$(repeat(" ", div(dashcount - len_den, 2)))$denstr"
  end
  println(io, numstr)
  println(io, repeat("-", dashcount))
  println(io, denstr)
end

function printtransferfunction(io::IO, s::CSisoRational)
  numstr = print_poly_reverse(s.num)
  denstr = print_poly_reverse(s.den)

  # Figure out the length of the separating line
  len_num = length(numstr)
  len_den = length(denstr)
  dashcount = max(len_num, len_den)

  # Center the numerator or denominator
  if len_num < dashcount
    numstr = "$(repeat(" ", div(dashcount - len_num, 2)))$numstr"
  else
    denstr = "$(repeat(" ", div(dashcount - len_den, 2)))$denstr"
  end
  println(io, numstr)
  println(io, repeat("-", dashcount))
  println(io, denstr)
end

function printtransferfunction(io::IO, s::CSisoZpk)
  numstr = sprint(print_polyroots, s.z, "z")
  denstr = sprint(print_polyroots, s.p, "z")
  gainstr = s.k[1]==1.0 ? "" : "$(round(s.k[1], 6))"

  # Figure out the length of the separating line
  len_num = length(numstr)
  len_den = length(denstr)
  len_gain = length(gainstr)
  dashcount = max(len_num, len_den)

  # Center the numerator or denominator
  if len_num < dashcount
    numstr = "$(repeat(" ", div(dashcount - len_num, 2)))$numstr"
  else
    denstr = "$(repeat(" ", div(dashcount - len_den, 2)))$denstr"
  end
  println(io, repeat(" ", len_gain+1), numstr)
  println(io, gainstr, " ", repeat("-", dashcount))
  println(io, repeat(" ", len_gain+1), denstr)
end

function printtransferfunction(io::IO, s::DSisoZpk)
  numstr = sprint(print_polyroots, s.z, "z")
  denstr = sprint(print_polyroots, s.p, "z")
  gainstr = s.k[1]==1.0 ? "" : "$(round(s.k[1], 6))"

  # Figure out the length of the separating line
  len_num = length(numstr)
  len_den = length(denstr)
  len_gain = length(gainstr)
  dashcount = max(len_num, len_den)

  # Center the numerator or denominator
  if len_num < dashcount
    numstr = "$(repeat(" ", div(dashcount - len_num, 2)))$numstr"
  else
    denstr = "$(repeat(" ", div(dashcount - len_den, 2)))$denstr"
  end
  println(io, repeat(" ", len_gain+1), numstr)
  println(io, gainstr, " ", repeat("-", dashcount))
  println(io, repeat(" ", len_gain+1), denstr)
end

function printterm{T}(p::Poly{T},j,first)
    s = ""
    pj = p[j]
    if pj == zero(T)
        return s
    end
    neg = pj < zero(T)
    if first
      if neg
        s = string("-")    #Prepend - if first and negative
      end
    else
        neg ? s = string(" - ") : s = string(" + ")
    end
    pj = abs(pj)
    if pj != one(T) || j == 0
      s = string(s,pj)
    end
    return string(s,printexponent(p.var,j))
end

function printexponent(var,i)
    if i == 0
      return string("")
    elseif i == 1
        return string(var)
    elseif var == symbol("z^-1")  # ugly trick to handle discrete system printing
        return string("z^",-i)
      else
        return string(var,"^",i)
    end
end

function print_poly_reverse{T}(p::Poly{T})
    first = true
    n = length(p)-1
    s = ""
    for i = n:-1:0
        s = string(s,printterm(p,i,first))
        first = false
    end
    if s == ""
      s = string(zero(T))
    end
    return s
end

function print_polyroots{T1<:Number}(io::IO, z::Vector{T1}, var=:x)
    z = z[imag(z) .>= -abs(z)*sqrt(eps(Float64))]
    n = length(z)
    if n == 0
        print(io, "1.0")
    else
        j = 1
        while length(z) != 0
            zj = z[j]
            tol = abs(zj) * sqrt(eps(Float64))
            if abs(zj) >= 2*eps(Float64)
                sgn = real(zj) < 0 ? (:+) : (:-)
                if imag(zj) >= tol
                    tmp = abs(round(2*real(zj), 6))
                    if tmp == 0
                        str = "($(var)^2 + $(round(abs(zj)^2, 6)))"
                    else
                        str = "($(var)^2 $sgn $tmp$var + $(round(abs(zj)^2, 6)))"
                    end
                else
                    str = "($var $sgn $(round(abs(real(zj)), 6)))"
                end
            else
                str = "$var"
            end
            inds = find(x->abs(x-zj) <= 2*tol, z)
            deleteat!(z, inds)
            exp = length(inds)
            if exp == 1
                print(io, str)
            else
                print(io, str, '^', exp)
            end
        end
    end
end
