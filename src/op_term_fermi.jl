import .Utils: pow_of_minus_one

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

function IsApprox.ishermitian(ft::AFermiTerm, approx_test::AbstractApprox=Equal())
    for x in op_string(ft)
        if FermiOps.is_raise_lower(x)
            return false
        end
    end
    return IsApprox.ishermitian(ft.coeff, approx_test)
end

####
#### Algebra / mathematical operations
####

"""
    count_ladder_ops(v::AbstractVector)

Return the number of explicit raising and lowering operators in `v`.
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

#count_bodies() = count_bodies(op_string(ft))

function count_bodies(ft::OpTerm)
    n = sum(FermiOps.count_raise_lower, op_string(ft))
    return n ÷ 2
end

# This may account correctly for phase
# Looks like each excitation gives a sign flip. That is, each
# -,+ pair gives a sign flip. So odd number of excitations gives sign flip
# No flip: number op, double excitation, coulomb interact.
# a flip: excitation op, number-excitation
function Base.adjoint(ft::AFermiTerm)
    n_excitations = count_ladder_ops(op_string(ft)) ÷ 2
    new_coeff = iseven(n_excitations) ? conj(ft.coeff) : -conj(ft.coeff)
    return strip_typeof(ft)(adjoint.(op_string(ft)), new_coeff)
end

## TODO: improve efficiency
function Base.adjoint(ft::AFermiSum)
    terms = [adjoint(x) for x in ft]
    return strip_typeof(ft)(terms)
end
