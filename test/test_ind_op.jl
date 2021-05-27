@testset "IndOp" begin
    op = IndOp(Paulis.X, 1)
    @test copy(op) == op

    t1 = PauliTerm("XIIYI", complex(2.0))
    t2 = PauliTerm("IZIYZ", complex(3.0))
    ot1 = OpTerm{IndOp}(t1)
    ot2 = OpTerm{IndOp}(t2)
    @test ot1 * ot2 == OpTerm{IndOp}(t1 * t2)

    ps = PauliSum([PauliTerm("IZYX", 1 + 0im),PauliTerm("XZXX", 1 + 0im),PauliTerm("ZXXX", 1 + 0im),])
    ops = OpSum{IndOp}(ps)
    @test map(weight, ps) == map(weight, ops)
    @test ops isa OpSum{<:IndOp}
    @test typeof(ops) == typeof(OpSum{IndOp}(ops))
end
