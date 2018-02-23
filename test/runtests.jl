using ObjectTracker
using Base.Test

@time @testset "Test blob creation" begin include("test_blobs.jl") end
