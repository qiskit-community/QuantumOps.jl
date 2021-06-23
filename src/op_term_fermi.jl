import ..Utils: pow_of_minus_one

####
#### OpTerm{AbstractFermiOp}
####

const FermiTerm = OpTerm{FermiOp}
const FermiSum = OpSum{FermiOp}
const AFermiTerm = OpTerm{<:AbstractFermiOp}
const DenseFermiTerm = DenseOpTerm{<:AbstractFermiOp}
const AFermiSum = OpSum{<:AbstractFermiOp}

macro fermi_str(str)
    return _op_term_macro_helper(FermiTerm, str)
end

####
#### Predicates
####

function LinearAlgebra.ishermitian(ft::AFermiTerm)
    for x in op_string(ft)
        if FermiOps.is_raise_lower(x)
            return false
        end
    end
    if LinearAlgebra.ishermitian(ft.coeff)
        return true
    end
    return false
end

####
#### Algebra / mathematical operations
####

"""
    count_ladder_ops(v::AbstractVector)

Return the number of raising and lowering operators in `v`.
"""
count_ladder_ops(v::AbstractVector) = count(FermiOps.is_raise_lower, v)

function fermi_phase_count(v1::AbstractVector, v2::AbstractVector)
    length(v1) == length(v2) || throw(DimensionMismatch("v1 an v2 differ in length"))
    _count = count_ladder_ops(v1)
    phase_count = 0
    @inbounds for i in eachindex(v1)
        if is_raise_lower(v1[i])
            _count -= 1
        end
        if !iseven(_count) && is_raise_lower(v2[i])
            phase_count += 1
        end
    end
    return phase_count
end

function AbstractOps.compute_phase(_::AbstractOps.DummyPhaseData, s1::AbstractArray{<:FermiOp}, s2::AbstractArray{<:FermiOp})
    phase_count = fermi_phase_count(s1, s2)
    return pow_of_minus_one(phase_count)
end

function Base.:^(ft::AFermiTerm, n::Integer)
    n < 0 && throw(DomainError(n))
    n == 0 && return one(ft)
    n == 1 && return ft
    ## Any +,-,0 sends the entire string to zero
    if any(x -> x === Raise || x === Lower || x === Zero, op_string(ft))
        return zero(ft)
    end
    ## E, N, I, are idempotent
    return strip_typeof(ft)(copy(op_string(ft)), ft.coeff^n)
end

count_bodies(ft::OpTerm) = count_bodies(op_string(ft))

function count_bodies(v::Vector)
    n = sum(FermiOps.count_raise_lower, v)
    return n ÷ 2
end

Base.adjoint(ft::AFermiTerm) = strip_typeof(ft)(adjoint.(op_string(ft)), conj(ft.coeff))

function Base.adjoint(ft::AFermiSum)
    _op_strings = [adjoint.(x) for x in ft.strings]
    coeffs = adjoint.(ft.coeffs)
    return strip_typeof(ft)(_op_strings, coeffs)
end
