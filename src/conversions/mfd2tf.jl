tf(s::MatrixFractionDescription{Val{:mimo}}) = tf(ss(s))

tf(s::MatrixFractionDescription{Val{:siso},Val{:cont}}) = tf(s.N, s.D)
tf(s::MatrixFractionDescription{Val{:siso},Val{:disc}}) = tf(s.N, s.D, s.Ts)
