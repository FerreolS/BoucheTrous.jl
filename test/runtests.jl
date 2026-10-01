using Test
using Aqua
using CartesianBoxes: CartesianBox
using ExplicitImports
using LinOps: LinOpGrad
using OffsetArrays
using BoucheTrous

@testset "BoucheTrous" begin
    Aqua.test_all(BoucheTrous)
    ExplicitImports.test_explicit_imports(BoucheTrous)
end

@testset "generic indexing" begin
    img = OffsetArray(zeros(Float64, 4, 4, 1), 0:3, 0:3, 1:1)
    badpix = OffsetArray(falses(4, 4, 1), 0:3, 0:3, 1:1)

    @test_throws ArgumentError bouchetrous!(img, badpix)
end

@testset "function construction" begin
    @test BoucheTrous.make_function(Float64, (4, 4), nothing) isa LinOpGrad
    @test typeof(BoucheTrous.make_function(Float64, (4, 4), 0.25)) != LinOpGrad
end

@testset "bad pixel box" begin
    badpix = trues(4, 5, 2)
    badpix[2:3, 3:4, :] .= false

    @test BoucheTrous.get_box(badpix) == CartesianBox(2:3, 3:4)
end

@testset "reconstruction" begin
    img = reshape(Float64[1, 2, 3, 4], 2, 2, 1)
    badpix = falses(2, 2, 1)
    badpix[1, 1, 1] = true

    @test isnothing(bouchetrous!(img, badpix; maxiter = 2))
    @test isfinite(img[1, 1, 1])
end
