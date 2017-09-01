type TestSystem{M1,M2}
  v::Vector{LTISystems.LtiSystem}
  sol::M1
  input::M2
  function (::Type{TestSystem})(v,sol,input)
    new{typeof(sol),typeof(input)}(v,sol,input)
  end
end

Base.start(s::TestSystem)       = start(s.v)
Base.next(s::TestSystem, state) = (s.v[state], state+1)
Base.done(s::TestSystem, state) = done(s.v, state)

T = 20.
time = (0, T)
simtol = 1e-4

function testsim(sim, sol, simtol)
  simsum = zero(sim.y[1])
  for (i,t) in enumerate(sim.t)
    simsum += sum(abs2.(sol(t)-sim.y[i,:]))
  end
  return simsum < simtol
end

# first continous system
# s1 = 1/(s+2) step response 1/2 - 1/2*exp(-2t)
c1v = Vector{LTISystems.LtiSystem}(0)
push!(c1v, tf([1],[1, 2]))
push!(c1v, lfd([1], [1, 2]))
push!(c1v, rfd([1], [1, 2]))
push!(c1v, ss(-2ones(1,1),ones(1,1),ones(1,1),0))
c1 = TestSystem(c1v, t->1/2-1/2*exp.(-2t), Signals.Step())

# second continuous system
# step response X(s)
# x₄(s) = 1/3/(s+1)-1/3/(s+4) => x₄(t) = 1/3e(-t)-1/3e(-4t)

# tf
d = Poly([4, 5, 1], :s)
c2tf = tf([RationalFunction([5, 1], d) RationalFunction(1, d);
           RationalFunction(-4,d)      RationalFunction([0, 1],d)])

# lfd/rfd
D = PolyMatrix([d zero(d); zero(d) d])
N = PolyMatrix([Poly([5.,1],:s) Poly(1.,:s); Poly(-4.,:s) Poly([0,1.],:s)])

c2lfd = lfd(N,D)
c2rfd = rfd(N,D)

# ss
A = [ 0.  0.  1.  0.;
      0.  0.  0.  1.;
     -4.  0. -5.  0.;
      0. -4.  0. -5.]
B = [ 0.  0.;
      0.  0.;
      1.  0.;
      0.  1.]
C = [ 5.  1.  1.  0.;
     -4.  0.  0.  1.]
D = zeros(2,2)
c2ss = ss(A,B,C,D)

c2sol = t->[1/4+1/12*exp(-4t)-1/3*exp(-t), 1/3*exp(-1t)-1/3*exp(-4t)]
c2inp = Signals.Step(zeros(2), [0, 1], zeros(2))

c2v = [c2tf, c2lfd, c2rfd, c2ss]
c2 = TestSystem(c2v, c2sol, c2inp)

# third continuous system

# tf
c3tf = tf([RationalFunction([0.,1.], poly([-1.,-1.,-2.,-2.],:s)) RationalFunction([0.,1.], [4.,4.,1.],:s);
           RationalFunction([0.,-1.], [4.,4.,1.],:s)    RationalFunction([0, -1.],  [4.,4.,1.],:s)])

# lfd/rfd
Dr = PolyMatrix([Poly([0.],:s) -Poly([2.,5.,4.,1],:s); Poly([4.,4.,1],:s) Poly([2.,1],:s)])
Nr = PolyMatrix([Poly([0.,1.],:s) Poly([0.,],:s); -Poly([0.,1.],:s) Poly([0.,0.,1.],:s)])

Dl = PolyMatrix([poly([-1.,-1.,-2],:s) Poly([2.,1],:s); Poly([0.],:s) Poly([4.,4.,1],:s)])
Nl = PolyMatrix([Poly([0.0],:s) Poly([0.,0.,1],:s); -Poly([0.,1.],:s) -Poly([0,1.],:s)])

c3rfd = rfd(Nr,Dr)
c3lfd = lfd(Nl,Dl)

# ss
A = [-4. 1 0  0 0;
     -5. 0 1 -1 0;
     -2. 0 0 -2 0;
      0  0 0 -4 1;
      0  0 0 -4 0]
B = [ 0. 1;
      0  0;
      0  0;
     -1 -1;
      0  0]
C = [1. 0 0 0 0;
     0  0 0 1 0]
D = zeros(2,2)
c3ss = ss(A,B,C,D)

c3sol = t->[t*exp(-2t), -t*exp(-2t)]
c3inp = Signals.Step(zeros(2), [0, 1.], zeros(2))

c3v = [c3tf, c3lfd, c3rfd, c3ss]
c3 = TestSystem(c3v, c3sol, c3inp)

# tf
c4tf  = tf([1.,6.,9], [1.,3.,2.])

# lfd/rfd
c4lfd = lfd([1.,6.,9], [1.,3.,2.])
c4rfd = rfd([1.,6.,9], [1.,3.,2.])

# ss
A = [ 0.  1.;
     -2. -3.]

B = [0. 1.].'

C = [7. 3.]
D = 1.
c4ss = ss(A,B,C,D)

c4sol = t->9/2+1/2*exp.(-2t)-4*exp.(-t)
c4inp = Signals.Step()

c4v = [c4tf, c4lfd, c4rfd, c4ss]
c4 = TestSystem(c4v, c4sol, c4inp)

@testset "continuous simulation" begin
  @testset "first order SISO system" begin
    for sys in c1
      sim = simulate(sys, time; input = c1.input, initial = zeros(numstates(sys)))
      @test testsim(sim, c1.sol, simtol)
    end
  end

  @testset "second order MIMO system" begin
    for sys in c2
      sim = simulate(sys, time; input = c2.input, initial = zeros(numstates(sys)))
      @test testsim(sim, c2.sol, simtol)
    end
  end

  @testset "MIMO system with different row and col degrees" begin
    for sys in c3
      sim = simulate(sys, time; input = c3.input, initial = zeros(numstates(sys)))
      @test testsim(sim, c3.sol, simtol)
    end
  end

  @testset "second order SISO system with direct term" begin
    for sys in c4
      sim = simulate(sys, time; input = c4.input, initial = zeros(numstates(sys)))
      @test testsim(sim, c4.sol, simtol)
    end
  end
end

# Discrete simulation
time = (0., T)
# first-order analytic response
fores(t,b,a) = t > 0 ? b*(1-a^t)/(1-a) : 0.

# first discrete system
# s1 = 0.4/(s-0.2) step response b*(1-a^(t))/(1-a) for t>0
b = 0.4
a = 0.2
d1sol = t->fores(t,b,a)
d1inp = Signals.Step([0.], [1.], [0.])

d1v = Vector{LTISystems.LtiSystem}(0)
push!(d1v, tf([b],[1, -a], 1))
push!(d1v, lfd([b],[1, -a], 1))
push!(d1v, rfd([b],[1, -a], 1))
push!(d1v, ss(a*ones(1,1),b*ones(1,1),ones(1,1),0,1))

d1 = TestSystem(d1v, d1sol, d1inp)

# second discrete system
b1  = 1/0.3
a1  = -0.2
a2  = -0.5
b21 = -2/3
b22 = 5/3

d2sol = t->[fores(t,b1,a1)+fores(t,-b1,a2), fores(t,b21,a1)+fores(t,b22,a2)]
d2inp = Signals.Step(zeros(2), [0, 1], zeros(2))

# tf
d = Poly([0.1, 0.7, 1], :z)
d2tf = tf([RationalFunction([5, 1], d) RationalFunction(1, d);
           RationalFunction(-4,d)      RationalFunction([0, 1],d)], 1.0)

# lfd/rfd
D = PolyMatrix([d zero(d); zero(d) d])
N = PolyMatrix([Poly([0.7,1],:z) Poly(1.,:z); Poly(-0.1,:z) Poly([0,1.],:z)])

d2lfd = lfd(N,D,1.0)
d2rfd = rfd(N,D,1.0)

# ss
A = [ 0.  0.   1.   0.;
      0.  0.   0.   1.;
    -0.1  0.  -0.7  0.;
      0. -0.1  0.  -0.7]
B = [ 0.  0.;
      0.  0.;
      1.  0.;
      0.  1.]
C = [ 0.7  1.  1.  0.;
     -0.1  0.  0.  1.]
D = zeros(2,2)
d2ss = ss(A,B,C,D,1.0)

d2v = [d2tf, d2lfd, d2rfd, d2ss]
d2 = TestSystem(d2v, d2sol, d2inp)

@testset "discrete simulation" begin
  @testset "first order SISO system" begin
    for sys in d1
      sim = simulate(sys, time; input = d1.input, initial = zeros(numstates(sys)))
      @test testsim(sim, d1.sol, simtol)
    end
  end

  @testset "second order MIMO system" begin
    for sys in d2
      sim = simulate(sys, time; input = d2.input, initial = zeros(numstates(sys)))
      @test testsim(sim, d2.sol, simtol)
    end
  end
end
