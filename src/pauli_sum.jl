export PauliSum

struct PauliSum{StringT, CoeffT}
    strings::StringT
    coeffs::CoeffT

    function PauliSum(strings, coeffs)
        if length(strings) != length(coeffs)
            throw(DimensionMismatch("bad dims"))
        end
        sort_and_sum_duplicates!(strings, coeffs)
        return new{typeof(strings), typeof(coeffs)}(strings, coeffs)
    end
end

sort_and_sum_duplicates!(psum::PauliSum) = sort_and_sum_duplicates!(psum.paulis, psum.coeffs)
function sort_and_sum_duplicates!(paulis, coeffs)
    sort_pauli_sum!(paulis, coeffs)
    sum_duplicates!(paulis, coeffs)
    return (paulis, coeffs)
end

sum_duplicates!(psum::PauliSum) = sum_duplicates!(psum.paulis, psum.coeffs)
function sum_duplicates!(paulis, coeffs)
    a = 2:length(paulis)
    k = 2
    y = first(paulis)
    coeff = first(coeffs)
    for j in a
        if paulis[j] != y
            y = paulis[k] = paulis[j]
            coeffs[k] = coeffs[j]
            k += 1
        else
            coeffs[k-1] += coeffs[j]
        end
    end
    resize!(paulis, k-1)
    resize!(coeffs, k-1)
    return paulis
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

Base.getindex(psum::PauliSum, i::Integer) = PauliString(psum.strings[i], psum.coeffs[i])
Base.getindex(psum::PauliSum, inds) = PauliSum(psum.strings[inds], psum.coeffs[inds])

# TODO: Should we check that the length of the PauliString is correct ?
# Probably not.
function Base.push!(psum::PauliSum, ps::PauliString)
    push!(psum.strings, ps.paulis)
    push!(psum.coeffs, ps.coeff)
    return ps
end

function Base.show(io::IO, psum::PauliSum)
    for i in eachindex(psum)
        show(io, psum[i])
        print(io, "\n")
    end
end
