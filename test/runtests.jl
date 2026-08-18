using QuantumCircuits
using Test
using LinearAlgebra
using Random

Random.seed!(20260817)

@testset "QuantumCircuits.jl" begin
    include("test_graycode.jl")
    include("test_gates.jl")
    include("test_circuit.jl")
    include("test_decompose.jl")
end
