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
