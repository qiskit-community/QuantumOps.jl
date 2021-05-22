"""
    _kron(mats...)

Same as `kron`, but, at least for Pauli matrices, it is faster by up to a factor of
two. Also, `kron(one_matrix)` returns `one_matrix` rather than throwing an error.
"""
function _kron(mats...)
    if length(mats) == 1
        return only(mats)
    end
    if length(mats) < 6
        return kron(mats...)
    else
        return mykron(mats...)
    end
end

"""
    mykron(mats...)

This is at most 2x faster than `kron`. It is often slower, especially for n < 5 matrices.
"""
function mykron(mats...)
    isempty(mats) && error("mykron requires at least one argument")
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
    return mykron(krons...)
end

"""
    isapprox_zero(x::Number)

This function exists because we define methods for special cases, such as symbolic
libraries.
"""
function isapprox_zero(x::Number)
    return isapprox(x, zero(x))
end
