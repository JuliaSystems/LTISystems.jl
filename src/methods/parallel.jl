# Parallel interconnection
"""
    parallel(s1,s2)

Constructs the parallel interconnection of systems `s1` and `s2`.

# Examples
```julia
julia> s1 = tf([1],[1,2]);
s2 = tf([3],[1,2]);
s3 = parallel(s1,s2);
showall(s3)
Continuous time rational transfer function model
	y = Gu

  4
-----
x + 2
```
"""
parallel{T1<:LtiSystem, T2<:LtiSystem}(s1::T1, s2::T2)     = +(s1, s2)

parallel{T1<:LtiSystem, T2<:Real}(s1::T1, n::T2)           = +(s1, n)
parallel{T1<:LtiSystem, T2<:Real}(n::T2, s1::T1)           = +(n, s1)

parallel{T1<:LtiSystem, M2<:AbstractMatrix}(s1::T1, n::M2) = +(s1, n)
parallel{T1<:LtiSystem, M2<:AbstractMatrix}(n::M2, s1::T1) = +(n, s1)
