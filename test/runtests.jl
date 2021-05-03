using PauliStrings
using Test

@testset "Pauli" begin
    @test Pauli(0) == Pauli(:I)
    @test Pauli(1) == Pauli(:X)
    @test Pauli(2) == Pauli(:Y)
    @test Pauli(3) == Pauli(:Z)
    @test Pauli(0) == Pauli('I')
    @test Pauli(1) == Pauli('X')
    @test Pauli(2) == Pauli('Y')
    @test Pauli(3) == Pauli('Z')
    @test Pauli(1) == Pauli("X")

    for (a, b, c) in ((:I, :I, :I), (:I, :X, :X), (:I, :Y, :Y), (:I, :Z, :Z),
                      (:X, :I, :X), (:X, :X, :I), (:X, :Y, :Z), (:X, :Z, :Y),
                      (:Y, :I, :Y), (:Y, :X, :Z), (:Y, :Y, :I), (:Y, :Z, :X),
                      (:Z, :I, :Z), (:Z, :X, :Y), (:Z, :Y, :X), (:Z, :Z, :I))
        @test Pauli(a) * Pauli(b) == Pauli(c)
    end

    m = rand(Pauli, (2, 2))
    @test typeof(m) == Matrix{Pauli}
    size(m) == (2, 2)
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
end
