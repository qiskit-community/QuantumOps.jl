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
end
