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
end

@testset "plus, minus, zero, one" begin
    for op in (PauliSum(["XXZ", "YYX"]), FermiSum(["NNE", "++-"]))
        @test length(op - op) == 0
        @test iszero(op - op)
        @test zero(op) == (op - op)
        @test iszero(0 * op)
        @test isone(one(op))
        @test !isone(op)
        @test !iszero(op)
    end
    for term in (PauliTerm("XXX"), FermiTerm("---"))
        @test iszero(0 * term)
        @test iszero(term - term)
        @test iszero(zero(term))
        @test isone(one(term))
        @test !isone(term)
        @test !iszero(term)
        # TODO test that isone is not the expensive fallback method
    end
    ft = FermiTerm("---")
    @test iszero(ft * ft)
    pt = PauliTerm("XYZ")
    @test isone(pt * pt)
end

@testset "indexing" begin
    for _sum = (rand_op_sum(Pauli, 10, 10), rand_op_sum(FermiOp, 10, 10))
        @test _sum[2][3:6] == _sum[2, 3:6]
        @test _sum[2, 3:6] isa Vector
        @test length(_sum[2:5]) == 4
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
