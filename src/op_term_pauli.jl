import ..AbstractPaulis
import ..AbstractPaulis: AbstractPauli
import .._op_term_macro_helper
import LightGraphs
import IsApprox: commutes

####
#### Constructors
####

const PauliTerm = OpTerm{Pauli}
const DensePauliTerm = DenseOpTerm{<:AbstractPauli}
const PauliSum = OpSum{Pauli}

## This ruins printing PauliTerm as an alias if `APauliTerm` is in the export list
"""
    APauliTerm

A `UnionAll` type representing a term of any type `<:AbstractPauli`.
"""
const APauliTerm = OpTerm{T} where {T <:AbstractPauli}

"""
    APauliSum

A `UnionAll` type representing a sum of terms of factors of any type `<:AbstractPauli`.
"""
const APauliSum = OpSum{T} where {T <:AbstractPauli}

# There is a printing bug in Julia that makes this a bad option
#const PauliTerm = APauliTerm{Pauli}
"""
    @pauli_str(str::String) -> PauliTerm

Construct a `PauliTerm` from `str`.

# Examples
```jldoctest
julia> pauli"XXY"
3-factor PauliTerm{Vector{Pauli}, Complex{Int64}}:
XXY * (1 + 0im)

julia> pauli"3 * XXY"
3-factor PauliTerm{Vector{Pauli}, ComplexF64}:
XXY * (3.0 + 0.0im)

julia> pauli"3.1 + 2.0im * XXY"
3-factor PauliTerm{Vector{Pauli}, ComplexF64}:
XXY * (3.1 + 2.0im)
```
"""
macro pauli_str(str)
    return _op_term_macro_helper(PauliTerm, str)
end

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
    return OpTerm(AbstractPaulis.pauli_vector(T, index, n_paulis), coeff)
end

####
####
####

function Base.zero(pterm::OpTerm{T}) where T <: AbstractPauli
    new_string = similar(op_string(pterm))
    fill!(new_string, one(T))
    return strip_typeof(pterm)(new_string, zero(pterm.coeff))
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
    matrix = kron_alt((Z4Group0.(m) for m in Matrix.(op_string(pt)))...)
    return _multiply_coefficient(pt.coeff, matrix)
end

function Base.Matrix(::Type{Float64}, pt::APauliTerm)
    matrix = kron_alt(Matrix.(op_string(pt))...)
    return _multiply_coefficient(pt.coeff, matrix)
end

SparseArrays.sparse(pt::APauliTerm) = SparseArrays.sparse(Float64, pt)

function SparseArrays.sparse(::Type{Float64}, pt::APauliTerm)
    matrix = kron_alt(SparseArrays.sparse.(op_string(pt))...)
    return _multiply_coefficient(pt.coeff, matrix)
end

# FIXME: broken, should throw inexact error earlier rather than return wrong type
function SparseArrays.sparse(::Type{Z4Group0}, pt::APauliTerm)
    matrix = kron_alt((Z4Group0.(m) for m in SparseArrays.sparse.(op_string(pt)))...)
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

## TODO: anticommutes. for FermiOp as well

commutes(p1::APauliTerm, p2::APauliTerm) = commutes(op_string(p1), op_string(p2))

## Each non-commuting factor introduces a phase factor of -1.
commutes(v1::AbstractArray{<:AbstractPauli}, v2::AbstractArray{<:AbstractPauli}) =
    iseven(count(x -> !commutes(x[1], x[2]), zip(v1, v2)))

####
#### Algebra
####

function Base.conj(pt::APauliTerm)
    return Base.transpose(pt, conj)
end

function Base.transpose(pt::APauliTerm, conj_func=identity)
    num_ys = count(p -> op_index(p) == 2, op_string(pt))
    fac = pow_of_minus_one(num_ys)
    return strip_typeof(pt)(op_string(pt), conj_func(pt.coeff) * fac)
end

function Base.adjoint(pt::APauliTerm)
    return strip_typeof(pt)(op_string(pt), conj(pt.coeff))
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

function LinearAlgebra.eigvals(pt::APauliTerm)
    if all(isone, op_string(pt)) # pt is propto identity
        return fill(float(pt.coeff), 2^length(pt)) # All eigenvalues are one
    end
    vals = Vector{promote_type(Float64, typeof(pt.coeff))}(undef, 2^length(pt))
    pos_eigval = pt.coeff
    neg_eigval = -pos_eigval
    if real(neg_eigval) > real(pos_eigval)
        (pos_eigval, neg_eigval) = (neg_eigval, pos_eigval)
    end
    half_length = 2^(length(pt)-1)
    @inbounds for i in 1:half_length
        vals[i] = neg_eigval  # follow the usual order
        vals[i+half_length] = pos_eigval
    end
    return vals
end

####
#### Other ?
####

"""
    pauli_basis(::Type{PauliT}, n_qubits; coeff=_DEFAULT_COEFF)

Return a `Generator` over all `PauliTerm`s of `n_qubits`.
"""
function pauli_basis(::Type{PauliT}, n_qubits; coeff=_DEFAULT_COEFF) where PauliT
    return (OpTerm{PauliT}(i, n_qubits, coeff) for i in 0:(4^n_qubits - 1))
end

pauli_basis(n_qubits; coeff=_DEFAULT_COEFF) = pauli_basis(PauliDefault, n_qubits; coeff=coeff)

"""
    group_paulis(ps::APauliSum)
    group_paulis(strings::AbstractVector{T}) where {T<:AbstractVector{<:AbstractPauli}}

Return a partition of `ps` or `strings` into groups of mutually commuting strings or terms.
The return type is `Vector{typeof(ps)}` or `Vector{typeof(strings)}`.

# Examples
```julia
julia> ps = rand_op_sum(Pauli, 5, 5; coeff_func=randn)
5x5 PauliSum{Vector{Vector{Pauli}}, Vector{Float64}}:
IXXZZ * 0.7287656446540858
IXZIZ * 0.6497331352636981
YIZZX * -0.7581101586820902
YZIXI * -1.9067637258086527
ZYIII * 1.3523639501623812

julia> group_paulis(ps)
4-element Vector{PauliSum{Vector{Vector{Pauli}}, Vector{Float64}}}:
 2x5 PauliSum{Vector{Vector{Pauli}}, Vector{Float64}}:
IXXZZ * 0.7287656446540858
YZIXI * -1.9067637258086527
 1x5 PauliSum{Vector{Vector{Pauli}}, Vector{Float64}}:
IXZIZ * 0.6497331352636981
 1x5 PauliSum{Vector{Vector{Pauli}}, Vector{Float64}}:
YIZZX * -0.7581101586820902
 1x5 PauliSum{Vector{Vector{Pauli}}, Vector{Float64}}:
ZYIII * 1.3523639501623812

julia> group_paulis(ps.strings)
4-element Vector{Vector{Vector{Pauli}}}:
 [IXXZZ, YZIXI]
 [IXZIZ]
 [ZYIII]
 [YIZZX]
```
"""
group_paulis(ps::APauliSum) = _group_paulis(ps, ps.strings)

group_paulis(strings::AbstractVector{T}) where {T<:AbstractVector{<:AbstractPauli}} =
    _group_paulis(strings, strings)

"""
    _group_paulis(_sum, strings)

Partition indices of `strings` into groups corresponding to mutually commuting terms
in `strings`. Return a partition of `_sum`, which must have the same indices as `strings`
according to partitioned indices. Note that `group_paulis(strings, strings)` simply
partitions `strings`.
"""
function _group_paulis(_sum, strings)
    graph = property_graph(strings, !commutes) # edges for non-commuting pairs
    coloring = LightGraphs.greedy_color(graph)
    groups = [empty(_sum) for i in 1:coloring.num_colors]
    for (i, c) in enumerate(coloring.colors)
        push!(groups[c], _sum[i])
    end
    return groups
end
