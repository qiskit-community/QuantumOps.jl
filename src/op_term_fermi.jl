####
#### OpTerm{AbstractFermiOp}
####

const FermiTerm = OpTerm{FermiOp}
const FermiSum = OpSum{FermiOp}
const AFermiTerm = OpTerm{<:AbstractFermiOp}
const DenseFermiTerm = DenseOpTerm{<:AbstractFermiOp}
#const FermiSums = OpSum{<:AbstractFermiOp}
const AFermiSum = OpSum{<:AbstractFermiOp}

## TODO: We renamed and moved this. Hope nothing untested uses it.
## TODO: This is too specialized to be a constructor method. Better to call it
## one_two_body_fermi_term or s.t.
# function OpTerm{T}(inds::NTuple{N, Int}, coeff, n_modes::Integer) where {T<:AbstractFermiOp, N}
#     (ops, phase) = index_to_ops_phase(inds)
#     (factors, new_coeff) = _fermi_term(ops, phase, coeff, n_modes)
#     return OpTerm{T}(factors, new_coeff)
# end

## TODO: again, this function should be named something more specific.
function OpSum{OpT}(tensor::AbstractArray{<:Number}, tensors::AbstractArray{<:Number}...) where {OpT<:AbstractFermiOp}
    fsum = OpSum(Vector{OpT}[], eltype(tensor)[])
    for tensor in (tensor, tensors...)
        tensor_to_fermi_sum!(fsum, tensor)
    end
    return fsum
end

macro fermi_str(str)
    return _op_term_macro_helper(FermiTerm, str)
end

"""
    FermiSumA(iop::ElectronicStructure.InteractionOperator)

Convert `iop` to a `FermiSum`.
"""
function OpSum{OpT}(iop::ElectronicStructure.InteractionOperator) where {OpT<:AbstractFermiOp}
    fsum = OpSum{OpT}(iop.one_body_tensor, iop.two_body_tensor)
    add!(fsum, iop.nuclear_repulsion * one(fsum[1]))
    return fsum
end

index_to_ops(inds) = index_to_ops(inds...)

function index_to_ops(i, j, k, l)
    ## We filter this elsewhere, but it's dangerous to remove this
    # if i == j || k == l
    #     return ()  # unstable, use a barrier or move this out
    # end
    if i == k || i == l
        op1 = (op=NumberOp, ind=i)
    else
        op1 = (op=Raise, ind=i)
    end
    if j == k || j == l
        op2 = (op=NumberOp, ind=j)
    else
        op2 = (op=Raise, ind=j)
    end
    if k != i && k != j
        op3 = (op=Lower, ind=k)
    else
        op3 = (op=FermiOps.NoOp, ind=k)
    end
    if l != i && l != j
        op4 = (op=Lower, ind=l)
    else
        op4 = (op=FermiOps.NoOp, ind=k)
    end
    return (op1, op2, op3, op4)
end

function index_to_ops(i, j)
    if i == j
        return ((op=NumberOp, ind=i), )
    end
    if i < j
        return ((op=Raise, ind=i), (op=Lower, ind=j))
    else
        return ((op=Lower, ind=i), (op=Raise, ind=j))
    end
end

function count_permutations(inds::NTuple{2})
    if inds[2] < inds[1]
        return 1
    else
        return 0
    end
end

index_phase(inds::NTuple) = pow_of_minus_one(count_permutations(inds))

function count_permutations(inds::NTuple{4})
    count = 0
    @inbounds i = inds[1]
    @inbounds for j in (2,3,4)
        if i > inds[j]
            count += 1
        end
    end
    @inbounds i = inds[2]
    @inbounds for j in (3,4)
        if i > inds[j]
            count += 1
        end
    end
    @inbounds i = inds[3]
    @inbounds if i > inds[4]
        count += 1
    end
    return count
end

index_to_ops_phase(inds) = (ops=index_to_ops(inds), phase=index_phase(inds))

## Embed `ops` in a string of width `n_modes`
function _fermi_term(ops::NTuple, phase::Integer, coeff, n_modes::Integer)
    factors = fill(I, n_modes)
    for (op, ind) in ops
        if op == FermiOps.NoOp
            continue
        end
        factors[ind] = op
    end
    return (factors, phase * coeff)
end

####
#### Algebra / mathematical operations
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

"""
    count_ladder_ops(v::AbstractVector)

Return the number of raising and lowering operators in `v`.
"""
count_ladder_ops(v::AbstractVector) = count(FermiOps.is_raise_lower, v)

function fermi_phase_count(v1::AbstractVector, v2::AbstractVector)
    _count = count_ladder_ops(v1)
    phase_count = 0
    for i in eachindex(v1)
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

Base.adjoint(ft::AFermiTerm) = strip_typeof(ft)(adjoint.(op_string(ft)), conj(ft.coeff))

count_bodies(ft::OpTerm) = count_bodies(op_string(ft))

function count_bodies(v::Vector)
    n = sum(FermiOps.count_raise_lower, v)
    return n ÷ 2
end

####
#### OpSum{<:AbstractFermiOp}
####

_skip_inds(i, j, k, l) = (i == j || k == l) # skip if two raising or two lowering ops on one index
_skip_inds(i, j) = false # always one raising and one lowering, so never skip

function tensor_term_to_fermi_term(::Type{T}, inds::NTuple{N, Int}, coeff, n_modes::Integer) where {T<:AbstractFermiOp, N}
    (ops, phase) = index_to_ops_phase(inds)
    (factors, new_coeff) = _fermi_term(ops, phase, coeff, n_modes)
    return OpTerm{T}(factors, new_coeff)
end

function tensor_to_fermi_sum!(fsum::OpSum{T}, tensor) where T
    n_modes = first(size(tensor))
    for ind in CartesianIndices(tensor)
        if _skip_inds((ind.I)...) || iszero(tensor[ind])
            continue
        end
        add!(fsum, tensor_term_to_fermi_term(T, ind.I, tensor[ind], n_modes))
    end
    return fsum
end

function Base.adjoint(ft::AFermiSum)
    _op_strings = [adjoint.(x) for x in ft.strings]
    coeffs = adjoint.(ft.coeffs)
    return strip_typeof(ft)(_op_strings, coeffs)
end
