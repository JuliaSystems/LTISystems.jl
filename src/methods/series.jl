# Series interconnection
"""
    series(s1,s2)

Constructs the series interconnection of systems `s1` and `s2`.

# Examples
```julia
julia> s1 = tf([1],[1,2]);
s2 = tf([3],[1,2]);
s3 = series(s1,s2);
showall(s3)
Continuous time rational transfer function model
	y = Gu

     3
------------
x^2 + 4x + 4
```
"""
series{T1<:LtiSystem, T2<:LtiSystem}(s1::T1, s2::T2)     = *(s1,s2)

series{T1<:LtiSystem, T2<:Real}(s1::T1, n::T2)           = *(s1, n)
series{T1<:LtiSystem, T2<:Real}(n::T2, s1::T1)           = *(n ,s1)

series{T1<:LtiSystem, M2<:AbstractMatrix}(s1::T1, n::M2) = *(s1, n)
series{T1<:LtiSystem, M2<:AbstractMatrix}(n::M2, s1::T1) = *(n ,s1)
