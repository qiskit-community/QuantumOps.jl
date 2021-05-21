abstract type AbstractTerm end

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
