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
series(s1::T1, s2::T2) where {T1<:LtiSystem, T2<:LtiSystem}     = *(s1,s2)

series(s1::T1, n::T2) where {T1<:LtiSystem, T2<:Real}           = *(s1, n)
series(n::T2, s1::T1) where {T1<:LtiSystem, T2<:Real}           = *(n ,s1)

series(s1::T1, n::M2) where {T1<:LtiSystem, M2<:AbstractMatrix} = *(s1, n)
series(n::M2, s1::T1) where {T1<:LtiSystem, M2<:AbstractMatrix} = *(n ,s1)
