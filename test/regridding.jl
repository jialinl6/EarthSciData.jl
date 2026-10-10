using EarthSciData
using Test
using Logging

@testset "regridding" begin
    @testset "data2vecormat" begin
        d = zeros(1, 2, 3)
        d2 = EarthSciData.data2vecormat(d, 2, 3)
        @test size(d2) == (6, 1)

        d = zeros(2, 3)
        d2 = EarthSciData.data2vecormat(d, 2, 1)
        @test size(d2) == (6,)
    end

    @testset "planar_regridder" begin
        square(x, y, d) = [(x, y), (x + d, y), (x + d, y + d), (x, y + d), (x, y)]
        dst = [square(x, y, 1.0) for y in 0:2 for x in 0:2]
        src = [square(x, y, 0.7) for y in 0.0:0.7:2.1 for x in 0.0:0.7:2.1]
        # No `crstrait` deprecation warning (logged at Warn under `--depwarn=yes`).
        rg = @test_logs min_level=Logging.Warn EarthSciData.planar_regridder(dst, src)
        ref = with_logger(NullLogger()) do
            EarthSciData.ConservativeRegridding.Regridder(dst, src)
        end
        @test rg.intersections == ref.intersections
        @test rg.dst_areas == ref.dst_areas
        @test rg.src_areas == ref.src_areas
    end
end
