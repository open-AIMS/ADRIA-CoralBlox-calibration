using Test
using Statistics

include("../src/common/constants.jl")
include("../src/common/dist_mean_penalty.jl")

# tabular, corymbose Acropora, corymbose non-Acropora, small massives, large massives
legacy_anchors = [3.751612251, 4.081622683, 4.487465256, 6.165751937, 7.153507902]

@testset "dist_mean_order_penalty" begin
    @test dist_mean_order_penalty(legacy_anchors) == 0.0

    @testset "free order within tiers" begin
        @test dist_mean_order_penalty([3.0, 5.0, 4.0, 7.0, 6.0]) == 0.0
    end

    @testset "violations" begin
        @test dist_mean_order_penalty([5.0, 4.0, 4.5, 6.0, 7.0]) > 0.0   # tabular above corymbose
        @test dist_mean_order_penalty([3.0, 4.0, 6.5, 6.0, 7.0]) > 0.0   # corymbose above small massives
        @test dist_mean_order_penalty([3.0, 4.0, 4.5, 3.5, 7.0]) > 0.0   # massive below corymbose
    end

    @testset "scales with violation size" begin
        small = dist_mean_order_penalty([4.5, 4.0, 4.5, 6.0, 7.0])
        large = dist_mean_order_penalty([6.0, 4.0, 4.5, 6.0, 7.0])
        @test 0.0 < small < large
        @test isfinite(large)
    end

    @testset "margin" begin
        @test dist_mean_order_penalty([4.0, 4.0, 4.5, 6.0, 7.0]) == 0.0
        @test dist_mean_order_penalty([4.0, 4.0, 4.5, 6.0, 7.0]; margin=0.1) > 0.0
    end

    @test_throws ArgumentError dist_mean_order_penalty([1.0, 2.0, 3.0])
end
