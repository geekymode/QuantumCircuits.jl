@testset "T-count optimisation" begin
    QC = QuantumCircuits
    rng = Random.Xoshiro(31)
    function randct(n, m)
        c = Circuit(n)
        for _ in 1:m
            if n >= 2 && rand(rng) < 0.5
                a, b = randperm(rng, n)[1:2]
                push!(c, rand(rng, [CNOT(), SWAP()]), a, b)
            else
                push!(c, rand(rng, [T(), T(), Tdg(), S(), Sdg(), Z(), RZ(π / 4), PHASE(3π / 4), RZ(-π / 2)]),
                      rand(rng, 1:n))
            end
        end
        c
    end
    parity_gadget!(c, ws, g) = (for q in ws[2:end]; push!(c, CNOT(), q, ws[1]); end;
                                push!(c, g, ws[1]);
                                for q in reverse(ws[2:end]); push!(c, CNOT(), q, ws[1]); end; c)

    @testset "phase polynomial form" begin
        for n in 1:4, _ in 1:10
            c = randct(n, 15n)
            κ, masks, phase = z8_phase_polynomial(c)
            # rebuild the unitary from (κ, masks, phase) by definition
            U = zeros(ComplexF64, 1 << n, 1 << n)
            for x in 0:(1 << n)-1
                f = sum(κ[y+1] * isodd(count_ones(y & x)) for y in 1:(1 << n)-1; init=0)
                out = sum((isodd(count_ones(masks[q] & x)) << (n - q) for q in 1:n); init=0)
                U[out+1, x+1] = cis(phase + π / 4 * f)
            end
            @test U ≈ matrix(c)
        end
        @test t_count(randct(3, 0)) == 0
        c = Circuit(2); push!(c, T(), 1); push!(c, S(), 2); push!(c, RZ(3π / 4), 1); push!(c, PHASE(π / 2), 2)
        @test t_count(c) == 2
        @test_throws ArgumentError z8_phase_polynomial(push!(Circuit(1), H(), 1))
        @test_throws ArgumentError z8_phase_polynomial(push!(Circuit(1), RZ(0.3), 1))
    end

    @testset "every Reed–Muller generator is free" begin
        for n in 4:7, (M, _) in QC._rm_generators(n), x in 0:(1 << n)-1
            @test mod(sum(isodd(count_ones(y & x)) for y in 1:(1 << n)-1 if y & M == M), 8) == 0
        end
        @test length.(QC._rm_generators.(3:7)) == [0, 1, 6, 22, 64]
    end

    @testset "optimisation is exact and never worse" begin
        for n in 2:6, _ in 1:(n <= 5 ? 6 : 3)
            c = randct(n, 25n)
            o = optimize_t_count(c)
            @test matrix(o) ≈ matrix(c)                         # global phase included
            κ, _, _ = z8_phase_polynomial(c)
            @test t_count(o) <= count(isodd, κ) <= t_count(c)
            @test all(op -> op.gate.name in (:CNOT, :SWAP, :T, :Tdg, :S, :Sdg, :Z), o.ops)
        end
        # exact search is optimal over the coset: compare with brute force at n = 5
        gens = QC._rm_generators(5)
        for _ in 1:20
            odd = rand(rng, UInt128) & ((UInt128(1) << 31) - 1)
            chosen = QC._rm_nearest_exact(odd, gens)
            w = odd; for i in chosen; w ⊻= gens[i][2]; end
            brute = minimum(count_ones(odd ⊻ reduce(⊻, (gens[i][2] for i in 1:6 if (s >> (i - 1)) & 1 == 1);
                                                   init=UInt128(0))) for s in 0:63)
            @test count_ones(w) == brute
        end
        # the local search (n = 7's fallback) is never better than exact, and close
        gens6 = QC._rm_generators(6)
        for _ in 1:5
            odd = rand(rng, UInt128) & ((UInt128(1) << 63) - 1)
            we = (w = odd; for i in QC._rm_nearest_exact(odd, gens6); w ⊻= gens6[i][2]; end; count_ones(w))
            ws = (w = odd; for i in QC._rm_nearest_search(odd, gens6); w ⊻= gens6[i][2]; end; count_ones(w))
            @test we <= ws <= we + 6
        end
        c7 = randct(7, 300)
        o7 = optimize_t_count(c7)
        @test matrix(o7) ≈ matrix(c7) && t_count(o7) < t_count(c7)
    end

    @testset "known cases" begin
        # T on all 15 parities of 4 qubits is the identity: 15 T → 0
        c = Circuit(4)
        for y in 1:15
            parity_gadget!(c, [q for q in 1:4 if (y >> (4 - q)) & 1 == 1], T())
        end
        @test matrix(c) ≈ I(16)
        @test t_count(c) == 15 && t_count(optimize_t_count(c)) == 0
        # CCZ needs 7 T, and nothing improves it at n = 3
        ccz = Circuit(3)
        for (ws, g) in (([1], T()), ([2], T()), ([3], T()), ([1, 2], Tdg()), ([1, 3], Tdg()),
                        ([2, 3], Tdg()), ([1, 2, 3], T()))
            parity_gadget!(ccz, ws, g)
        end
        @test matrix(ccz) ≈ Diagonal([1, 1, 1, 1, 1, 1, 1, -1])
        @test t_count(optimize_t_count(ccz)) == 7
        # T then T† merges away; a CNOT ladder alone stays T-free
        c = Circuit(2); push!(c, T(), 1); push!(c, CNOT(), 1, 2); push!(c, CNOT(), 1, 2); push!(c, Tdg(), 1)
        @test t_count(optimize_t_count(c)) == 0 && matrix(optimize_t_count(c)) ≈ I(4)
    end
end
