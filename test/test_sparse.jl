@testset "sparse terms" begin
    t = PauliTerm("IXIYIZ", 2.0 + im)
    @test t == sparse_op(t)
    @test dense_op(sparse_op(t)) == t
    @test typeof(dense_op(sparse_op(t))) == typeof(t)
    @test typeof(sparse_op(t)) != typeof(t)

    t1 = PauliTerm("IXIYIZ", 2.0 + im)
    t2 = PauliTerm("XIYYZY", 1.0 + 0im)
    f1 = FermiTerm("+I+I-I-", 2.0 + im)
    f2 = FermiTerm("--I-+++", 2.0 + im)
    (st1, st2, sf1, sf2) = sparse_op.((t1, t2, f1, f2))
    @test t1 * t2 == st1 * st2
    ## Make sure that sparse and dense multiplication use different methods
    @test which(*, (typeof(t1), typeof(t2))) !=  which(*, (typeof(st1), typeof(st2)))
    @test f1 * f2 == sf1 * sf2
    @test which(*, (typeof(f1), typeof(f2))) !=  which(*, (typeof(sf1), typeof(sf2)))

    ps = rand_op_sum(Pauli, 100, 100)
    sps = sparse_op(ps)
    @test sort!(sps) == sparse_op(ps)
end
