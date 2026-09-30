using Test
using Aqua
using ExplicitImports
using BoucheTrous

@testset "BoucheTrous" begin
    Aqua.test_all(BoucheTrous)
    ExplicitImports.test_explicit_imports(BoucheTrous)
end
