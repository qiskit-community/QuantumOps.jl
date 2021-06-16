## FIXME: random OpSum{FermiOp} sometimes gets zeros and the size of the sum is wrong
@testset "OpSum" begin
    for _type = (Pauli, PauliI, FermiOp)
        m = rand(_type, 5, 10)
        @test size(OpSum(m)) == (5, 10)
    end
    for _type = (Pauli, PauliI, FermiOp)
        @test length(OpSum{_type}()) == 0
    end
    for _type in (Pauli, FermiOp)
        _sum = rand_op_sum(_type, 5, 5)
        @test 3 * _sum == _sum * 3
    end

    ps = PauliSum(["XXZ", "YYX"])
    @test length(ps - ps) == 0
    @test iszero(ps - ps)
    @test zero(ps) == (ps - ps)
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
    @test OpTerm(FermiOps.Raise, FermiOps.Lower, FermiOps.NumberOp) == OpTerm{FermiOp}("+-N")

    pt = PauliTerm("XXIIYY")
    @test weight(pt) == 4
    @test length(pt) == 6
    ft = FermiTerm("++II--")
    @test weight(ft) == 4
    @test length(ft) == 6

    @test OpTerm{Pauli}([Paulis.X, Paulis.Y]) == OpTerm([Paulis.X, Paulis.Y])
    (X, Y, Z) = (Paulis.X, Paulis.Y, Paulis.Z)
    @test PauliTerm(@SVector [X, Y, Z, Y ,Z]) ==  PauliTerm([X, Y, Z, Y ,Z])
end
