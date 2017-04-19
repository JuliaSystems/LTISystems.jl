info("Starting tests for `TransferFunction`...")

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

@test num(syscs1) == num(syscs2) == num(syscs3) == num(syscs4) == num(syscs5)
@test den(syscs1) == den(syscs2) == den(syscs3) == den(syscs4) == den(syscs5)
@test num(sysds1) == num(sysds2) == num(sysds3) == num(sysds4) == num(sysds5) == num(sysds6)
@test den(sysds1) == den(sysds2) == den(sysds3) == den(sysds4) == den(sysds5) == den(sysds6)

# Rational transfer function constructions for MIMO systems
rc2   = RationalFunction([1,1], [3,1], :s)
matc  = diagm([rc1, rc2])
syscm = tf(matc)
rd2   = RationalFunction([-0.1,1], [-0.3,1], :z)
matd  = diagm([rd1, rd2])
sysdm = tf(matd, 0.5)

@test num(syscm) == [num(rc1) zero(num(rc1)); zero(num(rc1)) num(rc2)]
@test den(syscm) == [den(rc1) one(num(rc1)); one(num(rc1)) den(rc2)]
@test num(sysdm) == [num(rd1) zero(num(rd1)); zero(num(rd1)) num(rd2)]
@test den(sysdm) == [den(rd1) one(num(rd1)); one(num(rd1)) den(rd2)]

@test num(tf(diagm([1,2]))) == [Poly(1,:s) Poly(0,:s); Poly(0,:s) Poly(2,:s)]
@test den(tf(diagm([1,2]))) == [Poly(1,:s) Poly(1,:s); Poly(1,:s) Poly(1,:s)]

@test num(tf(diagm([1,2]), 0.5)) == [Poly(1,:z) Poly(0,:z); Poly(0,:z) Poly(2,:z)]
@test den(tf(diagm([1,2]), 0.5)) == [Poly(1,:z) Poly(1,:z); Poly(1,:z) Poly(1,:z)]

# Sampling time
@test samplingtime(syscs1) == samplingtime(syscm) == zero(Float64)
@test samplingtime(sysds1) == samplingtime(sysdm) == 0.5

# Input/Output information
@test numinputs(syscs1) == numinputs(sysds1)  == numoutputs(syscs1) == numoutputs(sysds1) == 1
@test numinputs(syscm)  == numinputs(sysdm)   == 2
@test numoutputs(syscm) == numoutputs(sysdm)  == 2

info("Tests for `TransferFunction` finished.")
