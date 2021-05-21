abstract type AbstractSum end

function _abstract_sum_inner_constructor_helper!(strings, coeffs, already_sorted=false)
    if length(strings) != length(coeffs)
        throw(DimensionMismatch("bad dims"))
    end
    if ! isempty(strings)
        n = length(first(strings))
        if ! all(x -> length(x) == n, strings)
            throw(DimensionMismatch("Fermi strings are of differing lengths."))
        end
        if ! already_sorted
            sort_and_sum_duplicates!(strings, coeffs)
        end
    end
    return nothing
end

## TODO: find a good way to abstract this
# function _abstract_sum_from_terms(v::AbstractVector, already_sorted=false)
#     strings = [x.ops for x in v]
#     coeffs = [x.coeff for x in v]
#     return FermiSum(strings, coeffs, already_sorted)
# end

function Base.copy(as::AbstractSum)
    (new_strings, new_coeffs) = (copy.(as.strings), copy(as.coeffs))
    already_sorted = true
    return typeof(as)(new_strings, new_coeffs, already_sorted)
end

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

####
#### Canonicalization / sorting
####

function sort_and_sum_duplicates!(asum::AbstractSum)
    sort_and_sum_duplicates!(fsum.strings, fsum.coeffs)
    return fsum
end

sum_duplicates!(asum::AbstractSum) = sum_duplicates!(fsum.strings, fsum.coeffs)

function sort_and_sum_duplicates!(terms, coeffs)
    sort_sums!(terms, coeffs)
    sum_duplicates!(terms, coeffs)
    remove_zeros!(terms, coeffs)
    return nothing
end

## This is expensive. Most time is spent in sortperm.
## There is no ThreadsX.sortperm, only sort.
Base.sort!(asum::AbstractSum) = (sort_sums!(asum.strings, asum.coeffs); asum)

function sort_sums!(strings, coeffs)
    p = sortperm(strings; alg=MergeSort)
    permute!(strings, p)
    permute!(coeffs, p)
    return nothing
end

"""
    remove_zeros!(asum::AbstractSum)
    remove_zeros!(terms, coeffs)

Remove terms from `asum` with coefficient (approximately) equal to zero.
If `terms` and `coeffs` are supplied, then elements are deleted from both `terms`
and `coeffs` at indices corresponding to vanishing elements of `coeff`.
"""
function remove_zeros!(asum::AbstractSum)
    remove_zeros!(asum.strings, asum.coeffs)
    return asum
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

####
#### Container interface
####

## Fails for empty psum
Base.size(asum::AbstractSum) = (length(asum), length(first(asum)))
Base.size(asum::AbstractSum, i::Integer) = size(asum)[i]

# Enables using `findall`, for instance.
# Fallback methods for `values` and `pairs` are OK.
Base.keys(asum::AbstractSum) = eachindex(asum)

# Iterate uses getindex to return `AbstractTerm`s.
function Base.iterate(asum::AbstractSum, state=1)
    state > lastindex(asum) && return nothing
    return (asum[state], state + 1)
end

for func in (:length, :eachindex, :lastindex, :firstindex)
    @eval begin
        function Base.$func(as::AbstractSum, args...)
            return $func(as.coeffs, args...)
        end
    end
end

Base.getindex(asum::AbstractSum, j::Integer) =
    term_type(typeof(asum))(asum.strings[j], asum.coeffs[j])

Base.getindex(asum::AbstractSum, j::Integer, k::Integer) = asum.strings[j][k]
# TODO: Use already_sorted flag ?

Base.getindex(asum::AbstractSum, inds) = typeof(asum)(asum.strings[inds], asum.coeffs[inds])

"""
    reverse(ps::AbstractSum)

Reverse qubit order in `ps` and sort terms.
"""
Base.reverse(as::AbstractSum) = reverse!(copy(as))

"""
    reverse!(ps::AbstractSum)

Reverse qubit order in `ps` in place and sort terms.
"""
function Base.reverse!(as::AbstractSum)
    strings = as.strings
    @inbounds for i in eachindex(strings)
        reverse!(strings[i])
    end
    return Base.sort!(as)
end

"""
    add!(psum::AbstractSum, pt::AbstractTerm...)

Add `AbstractTerm`s to `psum` in place, assuming `psum` is sorted and has no repeated
strings. Either a new term is inserted, or the coefficient is added to an existing
term. After adding the `pt`, `psum` will be left sorted, with no duplicates, and
no zero coefficients. Use `push!` to insert a term at the end of `psum` with no
simplification performed.
"""
add!(asum::AbstractSum, term::AbstractTerm) = add!(asum, op_string(term), term.coeff)

function add!(asum::AbstractSum, op_string, coeff)
    inds = searchsorted(asum.strings, op_string)
    if length(inds) == 0 # op_string not found, add a new term
        insert!(asum, first(inds), (op_string, coeff))
    elseif length(inds) == 1 # one element equal to op_string
        i = first(inds) # get the (single) index
        @inbounds asum.coeffs[i] += coeff # add p to existing term
        @inbounds if isapprox_zero(asum.coeffs[i])
            @inbounds deleteat!(asum, [i])
        end
    else
        throw(ErrorException("Duplicate terms found in operator sum."))
    end
    return asum
end

"""
    insert!(ps::AbstractSum, ind, p::AbstractTerm)

Insert `p` into `ps` without sorting resulting `ps`.
"""
Base.insert!(ps::AbstractSum, ind, p::AbstractTerm) = insert!(ps, ind, (op_string(p), p.coeff))

@inline function Base.insert!(ps::AbstractSum, ind, (paulis, coeff))
    insert!(ps.strings, ind, paulis)
    insert!(ps.coeffs, ind, coeff)
    return ps
end
