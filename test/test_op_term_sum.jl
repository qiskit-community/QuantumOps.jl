## FIXME: random OpSum{FermiOp} sometimes gets zeros and the size of the sum is wrong
@testset "OpSum" begin
    for _type = (Pauli, PauliI, FermiOp)
        m = rand(_type, 5, 10)
        @test size(OpSum(m)) == (5, 10)
    end
    for _type = (Pauli, PauliI, FermiOp)
        @test length(OpSum{_type}()) == 0
    end
end

@testset "OpTerm" begin
    for term in (PauliTerm("XYZ"), FermiTerm("--++"))
        @test reverse(term) != term
        @test reverse(reverse(term)) == term
        term2 = reverse(term)
        reverse!(term2)
        @test term2 == term
    end
    @test OpTerm(Paulis.X, Paulis.Y, Paulis.Z)  == OpTerm{Pauli}("XYZ")
    @test OpTerm(FermiOps.raise_op, FermiOps.lower_op, FermiOps.number_op) == OpTerm{FermiOp}("+-N")

    pt = PauliTerm("XXIIYY")
    @test weight(pt) == 4
    @test length(pt) == 6
    ft = FermiTerm("++II--")
    @test weight(ft) == 4
    @test length(ft) == 6
end
