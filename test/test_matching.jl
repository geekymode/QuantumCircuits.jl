@testset "matching decoder" begin
    rng = Random.Xoshiro(41)

    @testset "the correction explains the syndrome" begin
        for d in (3, 5, 7, 9)
            code = rotated_surface_code(d)
            dec = matching_decoder(code)
            n = length(code)
            for _ in 1:100
                E = PauliOp(rand(rng, n) .< 0.15, rand(rng, n) .< 0.15)
                s = syndrome(code, E)
                @test syndrome(code, decode(dec, s)) == s
            end
            @test decode(dec, falses(n - 1)) == PauliOp(falses(n), falses(n))
        end
    end

    @testset "corrects every error of weight ≤ ⌊(d-1)/2⌋" begin
        for d in (3, 5)
            code = rotated_surface_code(d)
            dec = matching_decoder(code)
            n = length(code)
            L = vcat(logical_operators(code)...)
            for w in 1:(d-1)÷2, pos in QuantumCircuits._combinations(n, w),
                kinds in Iterators.product(ntuple(_ -> 1:3, w)...)
                x = falses(n); z = falses(n)
                for (q, k) in zip(pos, kinds); x[q] = k != 3; z[q] = k != 1; end
                E = PauliOp(x, z)
                R = decode(dec, syndrome(code, E)) * E
                @test all(l -> commutes(R, l), L)
            end
        end
    end

    @testset "the matching is minimum weight" begin
        function brute(g, defects)                 # all perfect matchings with boundary
            B = size(g.dist, 1)
            isempty(defects) && return 0
            i, rest = defects[1], defects[2:end]
            best = g.dist[i, B] + brute(g, rest)
            for (k, j) in enumerate(rest)
                best = min(best, g.dist[i, j] + brute(g, deleteat!(copy(rest), k)))
            end
            best
        end
        g = matching_decoder(rotated_surface_code(7)).gx
        for _ in 1:100
            defects = sort(randperm(rng, length(g.checks))[1:rand(rng, 1:7)])
            @test QuantumCircuits._min_matching(g, defects)[2] == brute(g, defects)
        end
        # past maxdefects the greedy fallback still explains the syndrome
        code = rotated_surface_code(7)
        dec = matching_decoder(code; maxdefects=2)
        E = PauliOp(rand(rng, 49) .< 0.3, falses(49))
        s = syndrome(code, E)
        @test syndrome(code, decode(dec, s)) == s
    end

    @testset "other codes" begin
        # repetition code: matching is majority vote
        code = repetition_code(5)
        dec = matching_decoder(code)
        L = vcat(logical_operators(code)...)
        for b in 0:31
            x = BitVector([(b >> (j - 1)) & 1 == 1 for j in 1:5])
            R = decode(dec, syndrome(code, PauliOp(x, falses(5)))) * PauliOp(x, falses(5))
            @test (count(x) <= 2) == all(l -> commutes(R, l), L)
        end
        @test_throws ArgumentError matching_decoder(steane_code())     # qubit in 3 checks
        @test_throws ArgumentError matching_decoder(five_qubit_code()) # not CSS
        # lookup tables go through decode too
        table = lookup_decoder(steane_code())
        @test decode(table, falses(6)) == PauliOp(falses(7), falses(7))
        s1 = syndrome(steane_code(), PauliOp("IIIYIII"))
        @test decode(table, s1) == table[s1]
        # X₁Z₂ leaves different columns in the two halves: no single-qubit error explains it
        @test decode(table, syndrome(steane_code(), PauliOp("XZIIIII"))) === nothing
    end

    @testset "beats the lookup table" begin
        code = rotated_surface_code(5)
        rm = logical_error_rate(code, 0.05; shots=10_000, decoder=matching_decoder(code),
                                rng=Random.Xoshiro(1))[1]
        rl = logical_error_rate(code, 0.05; shots=10_000, decoder=lookup_decoder(code; maxweight=2),
                                rng=Random.Xoshiro(1))[1]
        @test rm < rl / 2                       # measured ~0.02 against ~0.13
        # and larger codes win at p = 3%, where the lookup decoder already fails
        r3 = logical_error_rate(rotated_surface_code(3), 0.03; shots=20_000,
                                decoder=matching_decoder(rotated_surface_code(3)), rng=Random.Xoshiro(2))[1]
        r7 = logical_error_rate(rotated_surface_code(7), 0.03; shots=20_000,
                                decoder=matching_decoder(rotated_surface_code(7)), rng=Random.Xoshiro(2))[1]
        @test r7 < r3
    end
end
