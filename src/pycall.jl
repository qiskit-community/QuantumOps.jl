Base.:*(z::PyCall.PyObject, ps::PauliTerm) = PauliTerm(ps.paulis, ps.coeff * z)
Base.:*(ps::PauliTerm, z::PyCall.PyObject) = *(z, ps)

function isapprox_zero(x::PyCall.PyObject)
    return x == 0
end
