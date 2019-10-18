print("Starting tests for `TransferFunction`...")

# Rational transfer function constructions for SISO systems
rc1   = RationalFunction([1], [2, 1], :s)     # 1/(s+2)
syscs1= tf(rc1)
syscs2= tf([1], [1, 2])
syscs3= tf(1, [1, 2])
syscs4= tf(Poly(1, :s), [1, 2])
syscs5= tf(1, Poly([2, 1], :s))
rd1   = RationalFunction([1], [-0.5, 1], :z)  # 1/(z-0.5)
sysds1= tf(rd1, 0.5)
sysds2= tf([1], [1, -0.5], 0.5)
sysds3= tf(1, [1, -0.5], 0.5)
sysds4= tf(Poly(1, :z), [1, -0.5], 0.5)
sysds5= tf(1, Poly([-0.5, 1], :z), 0.5)
sysds6= tf([0, 1], [1, -0.5], 0.5, :z̄)

@test numerator(syscs1) == numerator(syscs2) == numerator(syscs3) == numerator(syscs4) == numerator(syscs5)
@test denominator(syscs1) == denominator(syscs2) == denominator(syscs3) == denominator(syscs4) == denominator(syscs5)
@test numerator(sysds1) == numerator(sysds2) == numerator(sysds3) == numerator(sysds4) == numerator(sysds5) == numerator(sysds6)
@test denominator(sysds1) == denominator(sysds2) == denominator(sysds3) == denominator(sysds4) == denominator(sysds5) == denominator(sysds6)

# Rational transfer function constructions for MIMO systems
rc2   = RationalFunction([1,1], [3,1], :s)
matc  = diagm([rc1, rc2])
syscm = tf(matc)
rd2   = RationalFunction([-0.1,1], [-0.3,1], :z)
matd  = diagm([rd1, rd2])
sysdm = tf(matd, 0.5)

@test numerator(syscm) == [numerator(rc1) zero(numerator(rc1)); zero(numerator(rc1)) numerator(rc2)]
@test denominator(syscm) == [denominator(rc1) one(numerator(rc1)); one(numerator(rc1)) denominator(rc2)]
@test numerator(sysdm) == [numerator(rd1) zero(numerator(rd1)); zero(numerator(rd1)) numerator(rd2)]
@test denominator(sysdm) == [denominator(rd1) one(numerator(rd1)); one(numerator(rd1)) denominator(rd2)]

@test numerator(tf(diagm([1,2]))) == [Poly(1,:s) Poly(0,:s); Poly(0,:s) Poly(2,:s)]
@test denominator(tf(diagm([1,2]))) == [Poly(1,:s) Poly(1,:s); Poly(1,:s) Poly(1,:s)]

@test numerator(tf(diagm([1,2]), 0.5)) == [Poly(1,:z) Poly(0,:z); Poly(0,:z) Poly(2,:z)]
@test denominator(tf(diagm([1,2]), 0.5)) == [Poly(1,:z) Poly(1,:z); Poly(1,:z) Poly(1,:z)]

# Sampling time
@test samplingtime(syscs1) == samplingtime(syscm) == zero(Float64)
@test samplingtime(sysds1) == samplingtime(sysdm) == 0.5

# Input/Output information
@test numinputs(syscs1) == numinputs(sysds1)  == numoutputs(syscs1) == numoutputs(sysds1) == 1
@test numinputs(syscm)  == numinputs(sysdm)   == 2
@test numoutputs(syscm) == numoutputs(sysdm)  == 2

print("Tests for `TransferFunction` finished.")
