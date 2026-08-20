using QuantumCircuits
using Test
using LinearAlgebra
using Random

Random.seed!(20260817)

"""Haar random unitary, for checking decompositions."""
const randu = rand_unitary

"""Block-diagonal reference for a uniformly controlled one-qubit gate."""
function ref_multiplexed_gate(Us, controls, target, n)
    N = 1 << n
    k = length(controls)
    U = zeros(ComplexF64, N, N)
    tshift = n - target
    for col in 0:N-1
        j = 0
        for (i, q) in enumerate(controls)
            j |= ((col >> (n - q)) & 1) << (k - i)
        end
        m = Us[j+1]
        b = (col >> tshift) & 1
        for bp in 0:1
            row = (col & ~(1 << tshift)) | (bp << tshift)
            U[row+1, col+1] = m[bp+1, b+1]
        end
    end
    U
end

@testset "QuantumCircuits.jl" begin
    include("test_graycode.jl")
    include("test_gates.jl")
    include("test_circuit.jl")
    include("test_decompose.jl")
    include("test_mathkit.jl")
    include("test_matrixdecomp.jl")
    include("test_phasepoly.jl")
    include("test_applications.jl")
    include("test_shannon.jl")
    include("test_plots.jl")
end
