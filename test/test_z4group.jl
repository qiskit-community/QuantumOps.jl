@testset "Z4Group" begin
    for (n, x) in zip(1:4, (im, -1, -im, 1))
        @test z4group(n) == Z4Group(x)
    end
end

@testset "Z4Group0" begin
    for x in (im, -1, -im, 1, 0)
        @test Z4Group0(x) == x
    end
    for (x, y, z) in ((0, 0, 0), (0, 1, 0), (0, -1, 0), (0, im, 0), (0, -im, 0),
                      (1, 0, 0), (1, 1, 1), (1, -1, -1), (1, im, im), (1, -im, -im),
                      (-1, 0, 0), (-1, 1, -1), (-1, -1, 1), (-1, im, -im), (-1, -im, im),
                      (im, 0, 0), (im, 1, im), (im, -1, -im), (im, im, -1), (im, -im, 1),
                      (-im, 0, 0), (-im, 1, -im), (-im, -1, im), (-im, im, 1), (-im, -im, -1))
        @test Z4Group0(x) * Z4Group0(y) == Z4Group0(z)
    end
end
