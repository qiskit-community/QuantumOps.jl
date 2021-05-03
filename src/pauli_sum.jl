export PauliSum, add!

struct PauliSum{StringT, CoeffT}
    strings::StringT
    coeffs::CoeffT

    function PauliSum(strings, coeffs)
        if length(strings) != length(coeffs)
            throw(DimensionMismatch("bad dims"))
        end
        if ! isempty(strings)
            n = length(first(strings))
            if ! all(x -> length(x) == n, strings)
                throw(DimensionMismatch("Pauli strings are of differing lengths."))
            end
            sort_and_sum_duplicates!(strings, coeffs)
        end
        return new{typeof(strings), typeof(coeffs)}(strings, coeffs)
    end
end

PauliSum(strings) = PauliSum(strings, fill(_default_coeff, length(strings)))

function PauliSum(v::AbstractArray{<:PauliTerm})
    strings = getproperty.(v, :paulis)
    coeffs = getproperty.(v, :coeff)
    return PauliSum(strings, coeffs)
end

sort_and_sum_duplicates!(psum::PauliSum) = sort_and_sum_duplicates!(psum.paulis, psum.coeffs)
function sort_and_sum_duplicates!(paulis, coeffs)
    sort_pauli_sum!(paulis, coeffs)
    sum_duplicates!(paulis, coeffs)
    return nothing
end

sum_duplicates!(psum::PauliSum) = sum_duplicates!(psum.paulis, psum.coeffs)
function sum_duplicates!(paulis, coeffs)
    last_pauli = first(paulis)
    coeff = first(coeffs)
    k = 2
    for j in 2:length(paulis)
        if paulis[j] != last_pauli
            last_pauli = paulis[k] = paulis[j]
            coeffs[k] = coeffs[j]
            k += 1
        else
            coeffs[k-1] += coeffs[j]
        end
    end
    resize!(paulis, k-1)
    resize!(coeffs, k-1)
    return nothing
end

sort_pauli_sum!(psum::PauliSum) = sort_pauli_sum!(psum.paulis, psum.coeffs)
function sort_pauli_sum!(strings, coeffs)
    p = sortperm(strings)
    permute!(strings, p)
    permute!(coeffs, p)
    return nothing
end

for func in (:length, :eachindex)
    @eval begin
        function Base.$func(ps::PauliSum, args...)
            return $func(ps.coeffs, args...)
        end
    end
end

Base.getindex(psum::PauliSum, i::Integer) = PauliTerm(psum.strings[i], psum.coeffs[i])
Base.getindex(psum::PauliSum, inds) = PauliSum(psum.strings[inds], psum.coeffs[inds])

# TODO: Should we check that the length of the PauliTerm is correct ?
"""
    push!(psum::PauliSum, ps::PauliTerm...)

Push `ps` to the end of `psum` without regard to order
or possible duplication.

See `sort_and_sum_duplicates!`.
"""
function Base.push!(psum::PauliSum, ps::PauliTerm...)
    for p in ps
        push!(psum.strings, p.paulis)
        push!(psum.coeffs, p.coeff)
    end
    return ps
end

"""
    add!(psum::PauliSum, ps::PauliTerm...)

Add `PauliTerm`s to `psum`, assuming `psum` is sorted and
has no repeated strings. Either a new term is inserted, or
the coefficient is added to an existing term. After adding
the `ps`, `psum` will be left sorted and with no duplicates.
"""
function add!(psum::PauliSum, ps::PauliTerm...)
    for p in ps
        inds = searchsorted(psum.strings, p.paulis)
        if length(inds) == 0
            insert!(psum.strings, first(inds), p.paulis)
            insert!(psum.coeffs, first(inds), p.coeff)
        elseif length(inds) == 1
            psum.coeffs[first(inds)] += p.coeff
        else
            throw(ErrorException("Duplicate terms found in `PauliSum`."))
        end
    end
    return nothing
end

Base.:+(ps::PauliTerm...) = PauliSum(collect(ps))

function Base.show(io::IO, psum::PauliSum)
    for i in eachindex(psum)
        show(io, psum[i])
        print(io, "\n")
    end
end
