"""Reference DFT: column x is Σ_y exp(2πi·xy/N)|y⟩/√N."""
ref_dft(n) = (N = 1 << n; [cis(2π * x * y / N) for y in 0:N-1, x in 0:N-1] ./ sqrt(N))

@testset "quantum Fourier transform" begin
    @testset "exact" begin
        for n in 1:5
            F = ref_dft(n)
            c = qft(n)
            @test matrix(c) ≈ F
            @test matrix(qft(n; inverse=true)) ≈ F'
            @test count_gates(c, :H) == n
            @test count_gates(c, :CP) == n * (n - 1) ÷ 2
            @test count_gates(c, :SWAP) == n ÷ 2
            # lowered: 2 CNOTs per controlled phase, 3 per swap, same unitary
            cl = qft(n; lower=true)
            @test matrix(cl) ≈ F
            @test count_cnots(cl) == n * (n - 1) + 3 * (n ÷ 2)
            @test count_gates(cl, :CP) == 0 && count_gates(cl, :SWAP) == 0
            # without swaps the output is bit-reversed
            rev = [parse(Int, reverse(string(y; base=2, pad=n)); base=2) for y in 0:(1 << n)-1]
            @test matrix(qft(n; swaps=false))[rev .+ 1, :] ≈ F
        end
        # on arbitrary wires inside a larger register
        c = Circuit(4); qft!(c, [3, 1, 4])
        @test c.ops[1].qubits == [3]
        @test matrix(c) ≈ embed(ref_dft(3), [3, 1, 4], 4)
        @test_throws ArgumentError qft!(Circuit(2), Int[])
        @test_throws ArgumentError qft!(Circuit(2), [1, 1])
        @test_throws ArgumentError qft(3; cutoff=0)
    end

    @testset "approximate" begin
        n = 6
        F = ref_dft(n)
        errs = Float64[]
        for m in 1:n
            c = qft(n; cutoff=m)
            @test count_gates(c, :CP) == sum(min(n - j, m - 1) for j in 1:n)
            e = opnorm(matrix(c) - F)
            # bounded by the sum of the dropped rotation angles
            dropped = sum(2π / 2^d for j in 1:n for d in 2:n-j+1 if d > m; init=0.0)
            @test e <= dropped + 1e-12
            push!(errs, e)
        end
        @test errs[end] < 1e-12                       # cutoff ≥ n is exact
        @test issorted(errs; rev=true)
        @test matrix(qft(n; cutoff=10)) ≈ F
        # an approximate inverse undoes the matching approximate forward
        @test matrix(qft(n; cutoff=3, inverse=true)) ≈ matrix(qft(n; cutoff=3))'
    end

    @testset "centred DFT as a circuit (Eq. 47)" begin
        for k in 1:5
            c = centered_dft!(Circuit(k), 1:k)
            @test matrix(c) ≈ centered_dft(2^k)        # exact, global phase included
            @test count_gates(c, :P) == 2k
        end
        cl = centered_dft!(Circuit(4), 1:4; lower=true)
        @test matrix(cl) ≈ centered_dft(16)
        # the three-query algorithm's Fourier step, F̂ on every register, as gates
        c = Circuit(4); centered_dft!(c, 1:2); centered_dft!(c, 3:4)
        @test matrix(c) ≈ kron(centered_dft(4), centered_dft(4))
    end
end
