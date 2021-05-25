import QuantumOps.FermiOps

@testset "Jordan-Wigner" begin
    s1 = PauliTerm("ZZZXIIII", complex(-1/2)) +  PauliTerm("ZZZYIIII", im/2)
    @test jordan_wigner(FermiOps.raise_op, 4, 8) == s1

    s2 = PauliTerm("ZZZXIIII", complex(-1/2)) +  PauliTerm("ZZZYIIII", -im/2)
    @test jordan_wigner(FermiOps.lower_op, 4, 8) == s2

    s3 = PauliTerm{Pauli, Vector{Pauli}, ComplexF64}("IIIIIIII", 0.5 + 0.0im) +
        PauliTerm{Pauli, Vector{Pauli}, ComplexF64}("IIIZIIII", -0.5 + 0.0im)
    @test jordan_wigner(FermiOps.number_op, 4, 8) == s3
    @test s3 == s1 * s2

    s4 = PauliTerm{Pauli, Vector{Pauli}, ComplexF64}("IIIIIIII", 0.5 + 0.0im) +
        PauliTerm{Pauli, Vector{Pauli}, ComplexF64}("IIIZIIII", 0.5 + 0.0im)
    @test jordan_wigner(FermiOps.empty_op, 4, 8) == s4
    @test s4 == s2 * s1

    @test jordan_wigner_fermi(FermiOps.raise_op, 4, 8) == FermiTerm("ZZZ+IIII", 1.0 + 0.0im)
    @test jordan_wigner_fermi(FermiOps.lower_op, 4, 8) == FermiTerm("ZZZ-IIII", 1.0 + 0.0im)
    @test jordan_wigner_fermi(FermiOps.number_op, 4, 8) == FermiTerm("IIINIIII", 1.0 + 0.0im)

    # @test jordan_wigner_fermi(FermiOps.raise_op, 4, 8) == FermiTerm{FermiOp, Vector{FermiOp}, ComplexF64}("ZZZ+IIII", 1.0 + 0.0im)
    # @test jordan_wigner_fermi(FermiOps.lower_op, 4, 8) == FermiTerm{FermiOp, Vector{FermiOp}, ComplexF64}("ZZZ-IIII", 1.0 + 0.0im)
    # @test jordan_wigner_fermi(FermiOps.number_op, 4, 8) == FermiTerm{FermiOp, Vector{FermiOp}, ComplexF64}("IIINIIII", 1.0 + 0.0im)

end
