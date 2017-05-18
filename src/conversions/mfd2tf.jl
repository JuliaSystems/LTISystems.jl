tf(s::MFD{Val{:mimo}}) = tf(ss(s))

tf(s::MFD{Val{:siso},Val{:cont}}) = tf(s.N, s.D)
tf(s::MFD{Val{:siso},Val{:disc}}) = tf(s.N, s.D, s.Ts)
