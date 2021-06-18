module Utils

import LightGraphs

"""
    kron_alt(mats...)

Same as `kron`, but, at least for Pauli matrices, it is faster by up to a factor of
two. Also, `kron_alt(one_matrix)` returns `one_matrix` rather than throwing an error.
"""
function kron_alt(mats...)
    if length(mats) == 1
        return only(mats)
    end
    if length(mats) < 6
        return kron(mats...)
    else
        return _kron(mats...)
    end
end

"""
    _kron(mats...)

This is at most 2x faster than `kron`. It is often slower, especially for n < 5 matrices.
`_kron` is called by `kron_alt` when the number of matrices in `mats` is greater than
six.
"""
function _kron(mats...)
    isempty(mats) && error("_kron requires at least one argument")
    n = length(mats)
    if n == 1
        return only(mats)
    elseif n <= 3
        return kron(mats...)
    end
    half_n = n ÷ 2
    krons = Vector{Any}(undef, half_n)
    for i in 1:half_n
        krons[i] = kron(mats[2*i-1], mats[2*i])
    end
    if ! iseven(n)
        krons[half_n] = kron(krons[half_n], last(mats))
    end
    return _kron(krons...)
end

"""
    isapprox_zero(x::Number)

Return `true` if `x` is approximately zero.
This function exists because we define methods for special cases, such as symbolic
libraries.
"""
function isapprox_zero(x::Number)
    return isapprox(x, zero(x), atol=1e-16)
end

## MIME{Symbol("text/input")} is meant to print objects in input
## form, that is Julia code that constructs the object
Base.show(m::MIME{Symbol("text/input")}, item) = show(stdout, m, item)
Base.show(io::IO, ::MIME{Symbol("text/input")}, item) = show(io, item)

"""
    pow_of_minus_one(n::Integer)

Returns `-1` to the power `n`.
"""
pow_of_minus_one(n::Integer) = iseven(n) ? 1 : -1

"""
    property_graph(vector, property_func)

Return a `LightGraphs.SimpleGraph` with `length(vector)` integer-valued vertices, and edges
between exactly those indices of `vector` for which `property_func(vector[i], vector[j])` is
`true`. No self-edges are included.
"""
function property_graph(vector, property_func)
    n_verts = length(vector)
    graph = LightGraphs.SimpleGraph(n_verts)
    @inbounds for i in 1:n_verts, j in i+1:n_verts
        if property_func(vector[i], vector[j])
            LightGraphs.add_edge!(graph, i, j)
        end
    end
    return graph
end

"""
    n_cliffords(n)

Returns the number of elements in the Clifford group on `n` qubits.
"""
n_cliffords(n) = 8 * prod(i -> 2*(4^i - 1) * 4^i, big(1):n)

"""
    triprod(iter)

Iterate over only the upper triangular part of the Cartesian product of `iter` with `iter`,
returning `Tuple`s of elements

# Examples
```jldocstring
julia> collect(triprod('a':'d'))
6-element Vector{Tuple{Char, Char}}:
 ('a', 'b')
 ('a', 'c')
 ('a', 'd')
 ('b', 'c')
 ('b', 'd')
 ('c', 'd')
```
"""
triprod(iter) = ((a, b) for (n, a) in enumerate(iter) for b in Iterators.drop(iter, n))

end # module Utils
