@testset "Z4Group, Z4Group0 use" begin
    t1 = PauliTerm("XX", Z4Group0(1))
    t2 = PauliTerm("XX", 1)
    (m1, m2) = Matrix.((t1, t2))
    @test m1 == m2
    @test eltype(m1) == eltype(m2) == Float64
end

@testset "Pauli" begin
    for (n, s) in ((0, :I), (1, :X), (2, :Y), (3, :Z))
        @test Pauli(n) == Pauli(s)
        st = string(s)
        @test Pauli(n) == Pauli(st)
        @test Pauli(n) == Pauli(first(st))
    end

    for (a, b, c) in ((:I, :I, :I), (:I, :X, :X), (:I, :Y, :Y), (:I, :Z, :Z),
                      (:X, :I, :X), (:X, :X, :I), (:X, :Y, :Z), (:X, :Z, :Y),
                      (:Y, :I, :Y), (:Y, :X, :Z), (:Y, :Y, :I), (:Y, :Z, :X),
                      (:Z, :I, :Z), (:Z, :X, :Y), (:Z, :Y, :X), (:Z, :Z, :I))
        @test Pauli(a) * Pauli(b) == Pauli(c)
    end

    m = rand(Pauli, (2, 2))
    @test typeof(m) == Matrix{Pauli}
    @test size(m) == (2, 2)
end

@testset "PauliTerm" begin
    (I, X, Y, Z) = Pauli.(0:3)
    @test X ⊗ Y == PauliTerm("XY")
    @test PauliTerm("XY", 3).coeff == 3

    for (a, b, c, phase) in ((:I, :I, :I, 1), (:I, :X, :X, 1), (:I, :Y, :Y, 1), (:I, :Z, :Z, 1),
                             (:X, :I, :X, 1), (:X, :X, :I, 1), (:X, :Y, :Z, im), (:X, :Z, :Y, -im),
                             (:Y, :I, :Y, 1), (:Y, :X, :Z, -im), (:Y, :Y, :I, 1), (:Y, :Z, :X, im),
                             (:Z, :I, :Z, 1), (:Z, :X, :Y, im), (:Z, :Y, :X, -im), (:Z, :Z, :I, 1))
        (as, bs, cs) = String.((a, b, c))
        @test PauliTerm(as) * PauliTerm(bs) == phase * PauliTerm(cs)
    end
    @test PauliTerm("XYZ") * PauliTerm("ZYX") == PauliTerm("YIY")
    @test PauliTerm("XYZX") * PauliTerm("ZYXZ") == -im * PauliTerm("YIYY")
    @test PauliTerm("ZYXZ") * PauliTerm("XYZX") == im * PauliTerm("YIYY")
    t = PauliTerm("ZYXZ", 3)
    @test t^3 == t.coeff^2 * t
    @test t^4 == t.coeff^4 * one(t)
    @test PauliTerm("IXYZ") == PauliTerm([0, 1, 2, 3])
    @test PauliTerm("IXYZ", 2) == PauliTerm([0, 1, 2, 3], 2)
    @test PauliTerm(X, Y, Z) == PauliTerm("XYZ")
    @test PauliTerm(:XYZ) == PauliTerm("XYZ")

    @test sqrt(PauliTerm(X,Y)) == PauliSum(["II", "XY"], [0.5 + 0.5im, 0.5 - 0.5im])

    @test adjoint(PauliTerm("XYZ")) == PauliTerm("XYZ")
    @test conj(PauliTerm("XYZ")) ==  PauliTerm("XYZ", -1)
    @test transpose(PauliTerm("XYZ")) == PauliTerm("XYZ", -1)
    @test transpose(PauliTerm("XYZY")) == PauliTerm("XYZY")
end

@testset "PauliSum" begin
    ps = PauliSum(["YYY", "XXX"], [2, 3])
    @test length(ps) == 2
    @test [t for t in ps] == [PauliTerm("XXX", 3), PauliTerm("YYY", 2)]
    add!(ps, PauliTerm("III", -1))
    @test ps[1] == PauliTerm("III", -1)
    @test length(ps) == 3
    add!(ps, PauliTerm("XXX", -3))
    @test length(ps) == 2

    m = rand(8, 8)
    s = PauliSum(m)
    @test m ≈ Matrix(s)

    term = PauliTerm("XYZ")
    tan_term = tan(term)
    @test first(tan_term)  == tan(1) * term
    @test length(tan_term) == 1
    @test numeric_function(PauliTerm("XX"), x -> x^2 + 2.0 * x^3) == PauliSum([PauliTerm("II"), PauliTerm("XX", 2)])

    a = rand_op_sum(Pauli, 10, 3)
    b = rand_op_sum(Pauli, 10, 3)
    sa = SparseArrays.sparse(a)
    sb = SparseArrays.sparse(b)
    @test SparseArrays.sparse(a * b * a) == sa * sb * sa


    @test all(commutes, group_paulis(rand_op_sum(Pauli, 10, 100)))

end

@testset "Z4Group0" begin
    ms = Matrix.(Pauli.((0, 1, 2, 3)))
    zms = ((Z4Group0.(m) for m in ms)...,)
    @test zms == ms
    @test kron(zms...) == kron(ms...)
end

@testset "kron" begin
    t1 = PauliTerm("XYZ")
    t2 = PauliTerm("ZYX")
    @test Matrix(t1) ⊗ Matrix(t2) ≈ Matrix(t1 ⊗ t2)

    t1 = 3.0 * Paulis.X ⊗ Paulis.Y
    t2 = Paulis.X ⊗ Paulis.Y * 3.0
    @test t1 == t2
    @test_broken typeof(t1.coeff) == typeof(t2.coeff)
end

using LinearAlgebra: eigvals

@testset "eigvals" begin
    @test eigvals(PauliTerm("III")) == ones(2^3)
    for t in (PauliTerm("XYZX"), PauliTerm("XYI"))
        @test eigvals(t) ≈ eigvals(Matrix(t))
    end
end
