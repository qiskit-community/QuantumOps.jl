struct FermiSum{StringT, CoeffT}
    strings::StringT
    coeffs::CoeffT

    function FermiSum(strings, coeffs, already_sorted=false)
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
        return new{typeof(strings), typeof(coeffs)}(strings, coeffs)
    end
end

## Factor this out of here and PauliSum
function FermiSum(v::AbstractVector{<:FermiTerm}, already_sorted=false)
    strings = [x.ops for x in v]
    coeffs = [x.coeff for x in v]
    return FermiSum(strings, coeffs, already_sorted)
end

## TODO: Factor this out
function Base.show(io::IO, fsum::FermiSum)
    for i in eachindex(fsum)
        show(io, fsum[i])
        if i != lastindex(fsum)
            print(io, "\n")
        end
    end
end

## Factor this out
function sort_and_sum_duplicates!(fsum::FermiSum)
    sort_and_sum_duplicates!(fsum.strings, fsum.coeffs)
    return fsum
end

sum_duplicates!(fsum::FermiSum) = sum_duplicates!(fsum.strings, fsum.coeffs)

## Modeled on code in unique!
"""
    sum_duplicates!(op_string, coeffs)

Find groups of terms whose members differ only in the coefficient.
Replace each group by one term carrying the sum of the coefficients
in that group.
"""
function sum_duplicates!(op_string, coeffs)
    last_pauli::eltype(op_string) = first(op_string)
    coeff = first(coeffs)
    k = 2
    @inbounds for j in 2:length(op_string)
        if op_string[j] != last_pauli
            last_pauli = op_string[k] = op_string[j]
            coeffs[k] = coeffs[j]
            k += 1
        else
            coeffs[k-1] += coeffs[j]
        end
    end
    resize!(op_string, k-1)
    resize!(coeffs, k-1)
    return nothing
end

## TODO: Factor this out
Base.sort!(fsum::FermiSum) = (sort_sums!(fsum.strings, fsum.coeffs); fsum)

## Factor all this otu
"""
    remove_zeros!(psum::PauliSum)
    remove_zeros!(terms, coeffs)

Remove terms from `psum` with coefficient (approximately) equal to zero.
If `terms` and `coeffs` are supplied, then elements are deleted from both `terms`
and `coeffs` at indices corresponding to vanishing elements of `coeff`.
"""
function remove_zeros!(psum::FermiSum)
    remove_zeros!(fsum.strings, fsum.coeffs)
    return fsum
end



## Factor this out
for func in (:length, :eachindex, :lastindex, :firstindex)
    @eval begin
        function Base.$func(fs::FermiSum, args...)
            return $func(fs.coeffs, args...)
        end
    end
end

# Factor out
Base.getindex(fsum::FermiSum, j::Integer) = FermiTerm(fsum.strings[j], fsum.coeffs[j])


## Factor out
"""
    isapprox_zero(x::Number)

This function exists because we define
methods for special cases, such as symbolic
libraries.
"""
function isapprox_zero(x::Number)
    return isapprox(x, zero(x))
end
