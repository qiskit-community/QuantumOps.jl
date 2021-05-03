Base.:*(z::PyCall.PyObject, ps::PauliTerm) = PauliTerm(ps.paulis, ps.coeff * z)
Base.:*(ps::PauliTerm, z::PyCall.PyObject) = *(z, ps)

function Base.zero(x::PyCall.PyObject)
    return sympy.core.numbers.Zero
end

Base.isapprox(x::PyCall.PyObject, y::PyCall.PyObject) = x == y # a hack
