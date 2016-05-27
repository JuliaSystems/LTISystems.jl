# zpk
zpk(s::SisoSystem) = isdiscrete(s) ? zpk(zpkdata(s)..., samplingtime(s)) :
                                     zpk(zpkdata(s)...)

function zpk(s::MimoSystem)
  if isdiscrete(s)
    mimo(convert(Array{ControlCore.DSisoZpk}, map(zpk,getmatrix(s))))
  else
    mimo(convert(Array{ControlCore.CSisoZpk}, map(zpk,getmatrix(s))))
  end
end

# tf
tf(s::SisoSystem) = isdiscrete(s) ? tf(numpoly(s), denpoly(s), samplingtime(s)) :
                                    tf(numpoly(s), denpoly(s))

function tf(s::MimoSystem)
  if isdiscrete(s)
    mimo(convert(Array{ControlCore.DSisoZpk}, map(tf,getmatrix(s))))
  else
    mimo(convert(Array{ControlCore.CSisoZpk}, map(tf,getmatrix(s))))
  end
end

# ss

convert(::Type{DSisoZpk}, s::SisoSystem)      = zpk(s)
convert(::Type{DSisoRational}, s::SisoSystem) = tf(s)
