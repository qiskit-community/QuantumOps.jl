abstract type AbstractSum{StringT, CoeffT} end

####
#### Constructors
####

## Should we use this ?
#op_strings(as::AbstractSum) = as.strings

function _abstract_sum_inner_constructor_helper!(strings, coeffs; already_sorted=false)
    if length(strings) != length(coeffs)
        throw(DimensionMismatch("bad dims"))
    end
    if ! isempty(strings)
        if ! already_sorted  # Slightly dangerous to only do length check if not sorted
            n = length(first(strings))
            if ! all(x -> length(x) == n, strings)
                throw(DimensionMismatch("Fermi strings are of differing lengths."))
            end
            sort_and_sum_duplicates!(strings, coeffs)
        end
    end
    return nothing
end

function Base.similar(ps::AbstractSum{T, V}, n=0) where {W, C, V <:Vector{C}, T <: Vector{Vector{W}}}
    m = size(ps, 2)
    strings = [Vector{W}(undef, m) for i in 1:n]
    coeffs = Vector{C}(undef, n)
    return typeof(ps)(strings, coeffs; already_sorted=true)
end

## TODO: find a good way to abstract this
# function _abstract_sum_from_terms(v::AbstractVector, already_sorted=false)
#     strings = [x.ops for x in v]
#     coeffs = [x.coeff for x in v]
#     return FermiSum(strings, coeffs, already_sorted)
# end

function Base.copy(as::AbstractSum)
    (new_strings, new_coeffs) = (copy.(as.strings), copy(as.coeffs))
    return typeof(as)(new_strings, new_coeffs; already_sorted=true)
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
    sort_and_sum_duplicates!(asum.strings, asum.coeffs)
    return asum
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
Base.sort!(asum::AbstractSum; alg=MergeSort) = (sort_sums!(asum.strings, asum.coeffs; alg=alg); asum)

function sort_sums!(strings, coeffs; alg=MergeSort)
    p = sortperm(strings; alg=alg)
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
    # ThreadsX is very slow for small arrays
    # The findall(iszero, coeffs) is 500ns for two non-zero floats. What is wrong?
    # Appears to be this: iszero.(array) is taking almost all the time.
    # The following is what Base does, but writing it out is faster. A bug?
    # if length(coeffs) > 10^10
    #     inds = ThreadsX.findall(isapprox_zero.(coeffs))
    # else
    inds = findall(isapprox_zero.(coeffs))
    #    end
    if ! isempty(inds)
        deleteat!(coeffs, inds)
        deleteat!(terms, inds)
    end
    return nothing
end

## Modeled on code in unique! for sorted input
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
function Base.size(asum::AbstractSum)
    n = isempty(asum) ? 0 : length(first(asum))
    (length(asum), n)
end
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

####
#### Compare / predicates
####

# This will fail for an empty `psum`. Use type info instead.
# There no well-defined `one` for `PauliSum`. It depends on the
# width of the string.
function Base.one(asum::AbstractSum)
    t = one(first(asum))
    typeof(asum)([op_string(t)], [t.coeff]; already_sorted=true)
end

function Base.:(==)(asum1::AbstractSum, asum2::AbstractSum)
    if length(asum1) != length(asum2)
        return false
    end
    if length(asum1) > 10^3
        # This is 8x-12x faster for large arrays, 10^4 or 5.
        return ThreadsX.all(i -> asum1[i] == asum2[i], eachindex(asum1))
    end
    return all(i -> asum1[i] == asum2[i], eachindex(asum1))
end

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
    add!(to::AbstractSum, from::AbstractSum)

Adds the terms in `from` to `to` in place. `to` is mutated. `from` is not.
"""
function add!(to::AbstractSum, from::AbstractSum)
    for p in from
        add!(to, p)
    end
    return to
end

####
#### Updating / adding elements
####

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

function Base.deleteat!(ps::AbstractSum, args...)
    deleteat!(ps.coeffs, args...)
    deleteat!(ps.strings, args...)
    return ps
end

"""
    push!(psum::AbstractSum, ps::AbstractTerm...)

Push `ps` to the end of `psum` without regard to order or possible duplication.

See `sort_and_sum_duplicates!`.
"""
function Base.push!(psum::AbstractSum, ps::AbstractTerm...)
    for p in ps
        push!(psum.strings, op_string(p))
        push!(psum.coeffs, p.coeff)
    end
    return psum
end

## TODO: This should probably not be a method of push!
function Base.push!(psum::AbstractSum, (string, coeff))
    push!(psum.strings, string)
    push!(psum.coeffs, coeff)
end

function Base.append!(to::AbstractSum, from::AbstractSum)
    for t in from
        push!(to, t)
    end
    return to
end

####
#### Algebra / mathematical operations
####

Base.:+(terms::T...) where {T <: AbstractTerm} = sum_type(T)([terms...])

function Base.:+(ps0::AbstractSum, pss::AbstractSum...)
    ps_out = copy(ps0)
    for ps in pss
        add!(ps_out, ps)
    end
    return ps_out
end

## TODO: Do something more efficient here.
function Base.:-(pt1::T, pt2::T) where T <: AbstractTerm
    return sum_type(T)([pt1, -one(pt2.coeff) * pt2])
end

function Base.:-(psum::AbstractSum)
    typeof(psum)(psum.strings, -one(eltype(psum.coeffs)) .* psum.coeffs; already_sorted=true)
end

function Base.:*(n::Number, psum::AbstractSum)
    typeof(psum)(psum.strings, n .* psum.coeffs; already_sorted=true)
end

function Base.:/(psum::AbstractSum, n)
    typeof(psum)(psum.strings, psum.coeffs ./ n; already_sorted=true)
end

function mul!(asum_out::T, term::AbstractTerm, asum::T) where T <: AbstractSum
    @inbounds for j in eachindex(asum)
        new_term = term * asum[j]
        asum_out.coeffs[j] = new_term.coeff
#        copy!(asum_out.strings[j], op_string(new_term)) # hmm this is not defined
        asum_out.strings[j] = op_string(new_term)
    end
    return asum_out
end

## TODO: new_coeffs may be a Vector{<:Real}, but because of phase
## we need to set a coefficient to a complex type, which errors out.
## Currently, the user needs to make the coefficients of psum complex.
## We could also widen the type here somehow, say through promotion.
## A similar situation arises in other places in this library.
# function Base.:*(term::AbstractTerm, asum::AbstractSum)
#     new_coeffs = similar(asum.coeffs)
#     new_strings = similar(asum.strings)
#     ## Threading does not help; PauliSum() call is bottleneck
#     ## And it is slower anyway for small sums
#     @inbounds for j in eachindex(asum)
#         new_term = term * asum[j]
#         new_coeffs[j] = new_term.coeff
#         new_strings[j] = op_string(new_term)
#     end
#     return typeof(asum)(new_strings, new_coeffs)
# end

function Base.:*(term::AbstractTerm, asum::AbstractSum)
    asum_out = typeof(asum)(similar(asum.strings), similar(asum.coeffs); already_sorted=true)
    return sort!(mul!(asum_out, term, asum))
end

"""
    *(as1::AbstractSum, as2::AbstractSum)

Multiply `as1` and `as2` returning another `AbstractSum`

# Examples
```jldoctest
julia> a = rand_pauli_sum(5, 3);

julia> b = rand_pauli_sum(5, 3);

julia> a = rand_pauli_sum(5, 3); ma = Matrix(a);

julia> b = rand_pauli_sum(5, 3); mb = Matrix(b);

julia> Matrix(a * b * a) == ma * mb * ma
true
```
"""
function Base.:*(as1::AbstractSum, as2::AbstractSum)
    ## Using  a buffer and mul! ought to be faster, but it is not.
    ## because mul! takes no time here compared to add!
    asum_out = similar(as1)
    asum_cum = typeof(as2)(similar(as2.strings), similar(as2.coeffs); already_sorted=true)
    @inbounds for i in 1:length(as1)
        mul!(asum_cum, as1[i], as2)
#        add!(asum_out, asum_cum)
        append!(asum_out, asum_cum)
    end
#    return asum_out
    return sort_and_sum_duplicates!(asum_out)
end

# function Base.:*(as1::AbstractSum, others::AbstractSum...)
#     asum_out = similar(as1)
#     for as in others
#         @inbounds for i in 1:length(as1)
#             add!(asum_out, as1[i] * as)
#         end
#     end
#     return asum_out
# end

"""
    filter(f, as::AbstractSum)

Return a copy of `as` keeping only terms for which `f` is true.

# Examples
```julia
julia> filter(x -> count_bodies(x) == 1, fermi_sum)
IIIN * -0.47594871522096355
IINI * -0.47594871522096355
INII * -1.2524635735648981
NIII * -1.2524635735648981
```
"""
function Base.filter(f, as::AbstractSum)
    inds = findall(f, as)
    return typeof(as)(as.strings[inds], as.coeffs[inds])
end
