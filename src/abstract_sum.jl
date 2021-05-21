abstract type AbstractSum end

####
#### IO
####

function Base.show(io::IO, asum::AbstractSum)
    for i in eachindex(asum)
        show(io, asum[i])
        if i != lastindex(asum)
            print(io, "\n")
        end
    end
end

## Factor this out
function sort_and_sum_duplicates!(terms, coeffs)
    sort_sums!(terms, coeffs)
    sum_duplicates!(terms, coeffs)
    remove_zeros!(terms, coeffs)
    return nothing
end

function sort_sums!(strings, coeffs)
    p = sortperm(strings; alg=MergeSort)
    permute!(strings, p)
    permute!(coeffs, p)
    return nothing
end

function remove_zeros!(terms, coeffs)
    # ThreadsX is very slow for small arrays. We need to discriminate
    # inds = ThreadsX.findall(isapprox_zero, coeffs)
    # The following is 500ns for two non-zero floats. What is wrong?
    # Appears to be this: iszero.(array) is taking almost all the time.
    # The following is what Base does, but writing it out is faster. A bug.
    inds = findall(isapprox_zero.(coeffs))
    deleteat!(coeffs, inds)
    deleteat!(terms, inds)
    return nothing
end

## Modeled on code in unique!
"""
    sum_duplicates!(op_strings, coeffs)

Find groups of terms whose members differ only in the coefficient.
Replace each group by one term carrying the sum of the coefficients
in that group. This routine assumes that `op_strings` are sorted.
"""
function sum_duplicates!(op_strings, coeffs)
    last_pauli::eltype(op_strings) = first(op_strings)
    coeff = first(coeffs)
    k = 2
    @inbounds for j in 2:length(op_strings)
        if op_strings[j] != last_pauli
            last_pauli = op_strings[k] = op_strings[j]
            coeffs[k] = coeffs[j]
            k += 1
        else
            coeffs[k-1] += coeffs[j]
        end
    end
    resize!(op_strings, k-1)
    resize!(coeffs, k-1)
    return nothing
end
