struct FermiSum{StringT, CoeffT} <: AbstractSum
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

function term_type(::Type{T}) where T <: FermiSum
    return FermiTerm
end

## Factor this out of here and PauliSum
function FermiSum(v::AbstractVector{<:FermiTerm}, already_sorted=false)
    strings = [x.ops for x in v]
    coeffs = [x.coeff for x in v]
    return FermiSum(strings, coeffs, already_sorted)
end

