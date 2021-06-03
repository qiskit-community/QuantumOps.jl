const PauliSum = OpSum{Pauli}
const PauliTerm = OpTerm{Pauli}
#const APauliTerm = OpTerm{T} where {T <: AbstractPauli}
#const PauliTerm = APauliTerm{Pauli}

## This ruins printing PauliTerm as an alias if `APauliTerm` is in the export list
const APauliTerm = OpTerm{T} where {T <:AbstractPauli}
#const APauliTerm = OpTerm{<:AbstractPauli}

"""
    PauliSumA(::Type{PauliT}, matrix::AbstractMatrix{<:Number}; threads=true)

Construct a Pauli decomposition of `matrix`, that is, a `PauliSumA` representing `matrix`.
If `thread` is `true`, use a multi-threaded algorithm for increased performance.
"""
function OpSum{PauliT}(matrix::AbstractMatrix{<:Number}; threads=true) where PauliT <: AbstractPauli
    if threads
        return pauli_sum_from_matrix_threaded(PauliT, matrix)
    else
        return pauli_sum_from_matrix_one_thread(PauliT, matrix)
    end
end

function pauli_sum_from_matrix_one_thread(::Type{PauliT}, matrix::AbstractMatrix{<:Number}) where PauliT
    nside = LinearAlgebra.checksquare(matrix)
    n_qubits = ILog2.checkispow2(nside)
    denom = 2^n_qubits  # == nside
    s = OpSum{PauliT}()
    for pauli in pauli_basis(PauliT, n_qubits)
        mp = SparseArrays.sparse(pauli)  # Much faster than dense
        coeff = LinearAlgebra.dot(mp, matrix)
        if ! isapprox_zero(coeff)
            push!(s, (op_string(pauli), coeff / denom))  # a bit faster than PauliTermA for small `matrix` (eg 2x2)
        end
    end
    return s
end

function pauli_sum_from_matrix_threaded(::Type{PauliT}, matrix::AbstractMatrix{<:Number}) where PauliT
    nside = LinearAlgebra.checksquare(matrix)
    n_qubits = ILog2.checkispow2(nside)
    denom = 2^n_qubits  # == nside
    ## Create a PauliSumA for each thread, for accumulation.
    sums = [OpSum{PauliT}() for i in 1:Threads.nthreads()]
    Threads.@threads for j in 0:(4^n_qubits - 1)
        pauli = OpTerm{PauliT}(j, n_qubits)
        mp = SparseArrays.sparse(pauli)  # Much faster than dense
        coeff = LinearAlgebra.dot(mp, matrix)
        if ! isapprox_zero(coeff)
            push!(sums[Threads.threadid()], (op_string(pauli), coeff / denom))
        end
    end
    for ind in 2:length(sums)  # Collate results from all threads.
        add!(sums[1], sums[ind])
    end
    return sums[1]
end

function OpTerm{T}(index::Integer, n_paulis::Integer, coeff=_DEFAULT_COEFF) where T <: AbstractPauli
    return OpTerm(pauli_vector(T, index, n_paulis), coeff)
end

####
#### Math
####

"""
    cis(p::PauliTerm)::PauliSum

Compute ``\\exp(i p)``
"""
function Base.cis(pt::OpTerm{<:AbstractPauli})
    return OpSum([op_string(one(pt)), copy(op_string(pt))], [cos(pt.coeff), im * sin(pt.coeff)], true)
end

"""
    exp(p::PauliTerm)::PauliSum

Compute ``\\exp(p)``
"""
function Base.exp(pt::OpTerm{<:AbstractPauli})
    coeff = im * pt.coeff
    return OpSum([op_string(one(pt)), copy(op_string(pt))], [cos(coeff), -im * sin(coeff)], true)
end

"""
    numeric_function(pt::PauliTerm, f)::PauliSum

Compute `f(pt)` by decomposing `f` into odd and even functions.
"""
function numeric_function(pt::OpTerm{<:AbstractPauli}, f)
    c = pt.coeff
    fe = (f(c) + f(-c)) / 2  # Even term
    fo = (f(c) - f(-c)) / 2  # Odd term
    if isapprox_zero(fe)
        strings = [copy(op_string(pt))]
        coeffs = [fo]
    elseif isapprox_zero(fo)
        strings = [op_string(one(pt))]
        coeffs = [fe]
    else
        strings = [op_string(one(pt)), copy(op_string(pt))]
        coeffs = [fe, fo]
    end
    # sorting would take 30x longer
    return OpSum(strings, coeffs; already_sorted=true)
end

# Julia 1.5 does not have cispi
for f in (:cos, :sin, :tan, :sqrt, :sind, :sinpi, :cospi, :sinh, :tanh,
          :acos, :asin, :atan, :sec, :csc, :cot, :log, :log2, :log10,
          :log1p)
    @eval Base.$f(pt::OpTerm{<:AbstractPauli}) = numeric_function(pt, $f)
end

####
#### Conversion
####

function _multiply_coefficient(coeff, matrix)
    if isbool_and_equal(coeff, one(eltype(matrix)))
        return matrix
    end
    if coeff isa eltype(matrix)
        return LinearAlgebra.lmul!(coeff, matrix)
    end
    coeff1 = convert(promote_type(typeof(coeff), eltype(matrix)), coeff)
    return coeff1 .* matrix
end

Base.Matrix(pt::APauliTerm) = Matrix(Float64, pt)

# FIXME: Not type stable.
function Base.Matrix(::Type{Z4Group0}, pt::APauliTerm)
    matrix = _kron((Z4Group0.(m) for m in Matrix.(op_string(pt)))...)
    return _multiply_coefficient(pt.coeff, matrix)
end

function Base.Matrix(::Type{Float64}, pt::APauliTerm)
    matrix = _kron(Matrix.(op_string(pt))...)
    return _multiply_coefficient(pt.coeff, matrix)
end

SparseArrays.sparse(pt::APauliTerm) = SparseArrays.sparse(Float64, pt)

function SparseArrays.sparse(::Type{Float64}, pt::APauliTerm)
    matrix = _kron(SparseArrays.sparse.(op_string(pt))...)
    return _multiply_coefficient(pt.coeff, matrix)
end

# FIXME: broken, should throw inexact error earlier rather than return wrong type
function SparseArrays.sparse(::Type{Z4Group0}, pt::APauliTerm)
    matrix = _kron((Z4Group0.(m) for m in SparseArrays.sparse.(op_string(pt)))...)
    return _multiply_coefficient(pt.coeff, matrix)
end

####
#### Compare / predicates
####

"""
    isunitary(pt::PauliTermA)

Return `true` if `pt` is a unitary operator.
"""
IsApprox.isunitary(pt::APauliTerm) = IsApprox.isunitary(pt.coeff)

"""
    ishermitian(pt::PauliTermA)

Return `true` if `pt` is a Hermitian operator.
"""
LinearAlgebra.ishermitian(pt::APauliTerm) = isreal(pt.coeff)

####
#### Algebra
####

function mul!(target::AbstractArray{<:AbstractPauli}, ps1::APauliTerm, ps2::APauliTerm)
    s_new, phase = multiply_keeping_phase!(target, op_string(ps1), op_string(ps2))
    return strip_typeof(ps1)(s_new, ps1.coeff * ps2.coeff * phase)
end

#const DensePauliTerm = PauliTerm{<:Vector}
const DensePauliTerm = OpTerm{<:AbstractPauli, <:Vector}
function Base.:*(ps1::DensePauliTerm, ps2::DensePauliTerm)
    return mul!(similar(op_string(ps1)), ps1, ps2)
end

Base.inv(p::APauliTerm) = strip_typeof(p)(op_string(p), inv(p.coeff))

function Base.:^(p::APauliTerm, n::Integer)
    new_coeff = p.coeff^n
    if iseven(n)
        return strip_typeof(p)(fill(Pauli(:I), length(p)), new_coeff)
    else
        return strip_typeof(p)(op_string(p), new_coeff)
    end
end

function Base.conj(pt::APauliTerm)
    num_ys = count(p -> op_index(p) == 2, op_string(pt))
    fac = iseven(num_ys) ? 1 : -1
    return strip_typeof(pt)(op_string(pt), conj(pt.coeff) * fac)
end

function Base.transpose(pt::APauliTerm)
    num_ys = count(p -> op_index(p) == 2, op_string(pt))
    fac = iseven(num_ys) ? 1 : -1
    return strip_typeof(pt)(op_string(pt), pt.coeff * fac)
end

function Base.adjoint(pt::APauliTerm)
    return strip_typeof(pt)(op_string(pt), conj(pt.coeff))
end

function LinearAlgebra.eigvals(pt::APauliTerm)
    vals = Vector{promote_type(Float64, typeof(pt.coeff))}(undef, 2 * length(pt))
    pos_eigval = pt.coeff
    neg_eigval = -pos_eigval
    if real(neg_eigval) > real(pos_eigval)
        (pos_eigval, neg_eigval) = (neg_eigval, pos_eigval)
    end
    @inbounds for i in eachindex(pt)
        vals[i] = neg_eigval  # follow the usual order
        vals[i+length(pt)] = pos_eigval
    end
    return vals
end

Base.kron(ps1::APauliTerm, ps2::APauliTerm) = strip_typeof(ps1)(vcat(op_string(ps1), op_string(ps2)), ps1.coeff * ps2.coeff)

Base.kron(paulis::AbstractPauli...) = PauliTerm([paulis...])

## TODO: PauliTerm is hardcoded. Should be OpTerm{T} where T is some AbstractPauli
# TODO: @code_warntype shows red here
function Base.kron(ps::Union{APauliTerm, AbstractPauli}...)
    if ps[1] isa AbstractPauli
        v = typeof(ps[1])[]
    else
        v = eltype(ps[1])[]
    end
    coeffs = []
    for x in ps
        if x isa APauliTerm
            push!(coeffs, x.coeff)
            append!(v, x)
        else
            push!(v, x)
        end
    end
    if isempty(coeffs)
        return PauliTerm(v)
    else
        return PauliTerm(v, reduce(*, coeffs))
    end
end

####
#### Other ?
####

"""
    pauli_basis(::Type{PauliT}, n_qubits; coeff=_DEFAULT_COEFF)

Return a `Generator` over all `PauliTerm`s of `n_qubits`.
"""
function pauli_basis(::Type{PauliT}, n_qubits; coeff=_DEFAULT_COEFF) where PauliT
    return (OpTerm(PauliT, i, n_qubits, coeff) for i in 0:(4^n_qubits - 1))
end
