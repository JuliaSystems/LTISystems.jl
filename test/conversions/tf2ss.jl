# testcases taken from [1].
#
#   [1]: R. V. Patel, "Computation of minimal-order state-space realizations
#        and observability indices using orthogonal transforms",
#        International Journal of Automatic Control, vol. 33, no. 2, pp.
#        227-246, Mar. 1981.

r1 = RationalFunction(2.,[3,1.], :s)
r2 = RationalFunction(1.,[1,1.], :s)
Ptf = tf([r1 r2; r2 r2])

truep = [-3.0, -1.0, -1.0]
truez = [1.0]

Pss = ss(Ptf)
Plfd = lfd(LTISystems._tf2lfd(Ptf)...)
Pss = ss(Plfd)
Prfd = rfd(LTISystems._tf2rfd(Ptf)...)

@test sort(poles(Pss), lt = (l,r)->real(l)<real(r)) ≈ sort(truep)
@test zeros(Pss) ≈ truez

# test 2
charpoly = Poly([720., 3684., 8048, 9833, 7398, 3555, 1092, 207, 22, 1.], :s)

r11 = RationalFunction(3*poly([-3., -5.], :s), poly([-1., -2., -4.], :s))
r12 = RationalFunction(6*poly([-1.], :s), poly([-2., -4.], :s))
r13 = RationalFunction(Poly([7., 2.], :s), poly([-3., -4.], :s))
r14 = RationalFunction(Poly([5., 2.], :s), poly([-2., -3.], :s))
r21 = RationalFunction(Poly([2.], :s), poly([-3., -5.], :s))
r22 = RationalFunction(Poly([1.], :s), poly([-3.], :s))
r23 = RationalFunction(2*poly([-5.], :s), poly([-1., -2., -3.], :s))
r24 = RationalFunction(8*poly([-2.], :s), poly([-1., -3., -5.], :s))
r31 = RationalFunction(2*Poly([18., 7., 1.], :s), poly([-1., -3., -5.], :s))
r32 = RationalFunction(-Poly([0,2.], :s), poly([-1., -3.], :s))
r33 = RationalFunction(Poly([1.], :s), poly([-3.], :s))
r34 = RationalFunction(2*Poly([34., 27., 5.], :s), poly([-1., -3., -5.], :s))

Ptf = tf([r11 r12 r13 r14; r21 r22 r23 r24; r31 r32 r33 r34])

Plfd = lfd(LTISystems._tf2lfd(Ptf)...)
# ss(Plfd) # TODO find bug in _tf2lfd

Prfd = rfd(LTISystems._tf2rfd(Ptf)...)
ss(Prfd)

poly(poles(ss(Ptf)), :s)  ≈ charpoly
#poly(poles(ss(Plfd)), :s) ≈ charpoly
poly(poles(ss(Prfd)), :s) ≈ charpoly
