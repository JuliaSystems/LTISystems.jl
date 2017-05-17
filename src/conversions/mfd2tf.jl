tf(s::MFD{Val{:mimo}}) = tf(ss(s))

tf(s::MFD{Val{:siso},Val{:cont}}) = tf(s.N[1], s.D[1])
tf(s::MFD{Val{:siso},Val{:disc}}) = tf(s.N[1], s.D[1], s.Ts)
