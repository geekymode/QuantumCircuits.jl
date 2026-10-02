@testset "GF(2) and classical codes" begin
    @testset "GF(2) linear algebra" begin
        for _ in 1:20
            A = rand(Bool, rand(1:6), rand(1:9))
            R, piv = gf2_rref(A)
            N = gf2_nullspace(A)
            @test gf2_rank(A) + size(N, 1) == size(A, 2)           # rank–nullity
            @test all(iszero, (Int.(A) * Int.(N)') .% 2)
            @test gf2_rank(R) == length(piv) == gf2_rank(A)
            for (i, p) in enumerate(piv)                            # pivots are lone 1s
                @test count(R[:, p]) == 1 && R[i, p]
            end
        end
        @test gf2_rank([1 1; 1 1]) == 1 && gf2_rank([1 0; 1 1]) == 2
        @test gf2_rank([2 4; 6 8]) == 0                             # read mod 2
    end

    @testset "Hamming codes" begin
        for r in 2:4
            h = hamming_code(r)
            @test (length(h), dimension(h), minimum_distance(h)) == (2^r - 1, 2^r - 1 - r, 3)
            # the syndrome of a single flip is its position in binary
            for pos in (1, 2^r - 1)
                e = falses(length(h)); e[pos] = true
                @test syndrome(h, e) == BitVector(reverse(digits(pos; base=2, pad=r)))
            end
        end
        h = hamming_code(3)
        for m in Iterators.product(ntuple(_ -> 0:1, 4)...)
            w = encode(h, collect(m))
            @test iscodeword(h, w)
            for pos in 1:7
                y = copy(w); y[pos] = !y[pos]
                @test !iscodeword(h, y) && syndrome_decode(h, y) == w
            end
        end
        @test length(codewords(h)) == 16
        @test_throws ArgumentError encode(h, [1, 0])
    end

    @testset "Reed–Muller codes" begin
        for (r, m) in ((0, 3), (1, 3), (2, 3), (1, 4), (2, 4), (1, 5), (2, 5), (3, 5), (2, 6), (1, 7))
            c = reed_muller(r, m)
            @test length(c) == 2^m
            @test dimension(c) == sum(binomial(m, i) for i in 0:r)
            @test minimum_distance(c) == 2^(m - r)
            if r <= m - 1                                           # RM(r,m)⊥ = RM(m-r-1,m)
                d = reed_muller(m - r - 1, m)
                @test gf2_rank(vcat(dual(c).G, d.G)) == dimension(d) == dimension(dual(c))
            end
        end
        # rows are monomials: encode evaluates the polynomial
        c = reed_muller(1, 3)
        @test Int.(generator_matrix(c)) == [1 1 1 1 1 1 1 1; 0 0 0 0 1 1 1 1; 0 0 1 1 0 0 1 1; 0 1 0 1 0 1 0 1]
        @test encode(c, [1, 1, 0, 0]) == BitVector([x < 4 for x in 0:7])   # 1 + x₁
        # punctured RM(1,3) is the [7,4,3] Hamming code; RM(1,3) is self-dual
        p = puncture(c, [1])
        @test (length(p), dimension(p), minimum_distance(p)) == (7, 4, 3)
        @test gf2_rank(vcat(p.G, hamming_code(3).G)) == 4
        @test gf2_rank(vcat(c.G, dual(c).G)) == 4
        @test_throws ArgumentError reed_muller(4, 3)
    end

    @testset "Reed's local decoding" begin
        rng = Random.Xoshiro(7)
        c = reed_muller(2, 6)
        w = encode(c, rand(rng, Bool, dimension(c)))
        @test all(rm_local_decode(w, 2, 6, x; rng=rng) == w[x+1] for x in 0:63)   # clean: exact
        y = copy(w)
        for i in randperm(rng, 64)[1:2]; y[i] = !y[i]; end
        right = count(rm_local_decode(y, 2, 6, x; trials=31, rng=rng) == w[x+1] for x in 0:63)
        @test right >= 60                                           # majority vote, 2 errors
        @test_throws ArgumentError rm_local_decode(w, 2, 6, 64)
        @test_throws ArgumentError rm_local_decode(w[1:10], 2, 6, 0)
    end
end
