abstract type AbstractTerm{W, CoeffT} end

Base.copy(pt::AbstractTerm) = typeof(pt)(copy(op_string(pt)), copy(pt.coeff))

## type params only to get correct dispatch. There must be a better way
function Base.show(io::IO, term::AbstractTerm) # where { T, V, CoeffT, AbstractTerm{T,V,CoeffT}}
    print(io, op_string(term))
    print(io, " * ")
    if term.coeff isa Real  # could use CoeffT here.
        print(io, term.coeff)
    else
        print(io, "(", term.coeff, ")")
    end
end

function Base.show(io::IO, ps::AbstractTerm{T,Z4Group}) where {T}
    print(io, ps.coeff, " ")
    print(io, op_string(ps))
end

Base.show(m::MIME{Symbol("text/input")}, term::AbstractTerm) = show(stdout, m, term)
function Base.show(io::IO, mime::MIME{Symbol("text/input")}, term::AbstractTerm)
    print(io, strip_typeof(term), "(")
    print(io, "\"", op_string(term), "\"")
    print(io, ", ")
    show(io, mime, term.coeff)
    print(io, ")")
end

Base.one(term::AbstractTerm{W}) where {W} = typeof(term)(fill(one(W), length(term)), one(term.coeff))

Base.:(==)(op1::AbstractTerm, op2::AbstractTerm) = op1.coeff == op2.coeff && op_string(op1) == op_string(op2)

_isless(x, y) = isless(x, y)
_isless(x::Complex, y::Complex) = isless(abs2(x), abs2(y))

function Base.isless(op1::AbstractTerm, op2::AbstractTerm)
    if op_string(op1) == op_string(op2)
        return _isless(op1.coeff, op2.coeff)
    end
    return isless(op_string(op1), op_string(op2))
end

####
#### Container interface
####

# :popat!
for func in (:length, :size, :eltype, :eachindex, :axes, :splice!, :getindex,
             :setindex!, :iterate, :pop!, :popfirst!, :lastindex, :firstindex, :first, :last)
    @eval begin
        Base.$func(ps::AbstractTerm, args...) = $func(op_string(ps), args...)
    end
end

## TODO: Is there a way to combine this with methods for Base functions above ?
for func in (:weight, )
    @eval begin
        $func(ps::AbstractTerm, args...) = $func(op_string(ps), args...)
    end
end

for func in (:push!, :pushfirst!, :insert!)
    @eval begin
        Base.$func(ps::AbstractTerm, args...) = ($func(op_string(ps), args...); ps)
    end
end

Base.iszero(term::AbstractTerm) = iszero(term.coeff) || any(iszero, op_string(term))

"""
    reverse!(pt::AbstractTerm)

Reverse the order of the factors in `pt` in place.
See `reverse`.
"""
Base.reverse!(pt::AbstractTerm) = (reverse!(op_string(pt)); pt)

"""
    reverse(pt::AbstractTerm)

Reverse the order of the factors in `pt`.
See `reverse!`.
"""
Base.reverse(pt::AbstractTerm) = reverse!(copy(pt))

####
#### Algebra
####

Base.:*(z::Number, term::AbstractTerm) = strip_typeof(term)(op_string(term), term.coeff * z)

Base.:*(ps::AbstractTerm, z::Number) = z * ps

Base.:*(z::Number, p::AbstractOp) = term_type(typeof(p))([p], z)
Base.:*(p::AbstractOp, z::Number) = z * p

Base.:-(p::AbstractTerm) = typeof(p)(op_string(p), -p.coeff)

Base.:/(p::AbstractOp, z::Number) = term_type(typeof(p))([p], inv(z))

Base.:/(ps::AbstractTerm, z::Number) = typeof(ps)(op_string(ps), ps.coeff / z)
