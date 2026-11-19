using Test
using CommunityMatrices
using LinearAlgebra: Diagonal, I, eigvals, norm
using Random: Xoshiro

@testset "CommunityMatrices" begin
    include("test_matrices.jl")
    include("test_timeseries.jl")
    include("test_regression.jl")
    include("test_community.jl")
    include("test_docstrings.jl")
end
