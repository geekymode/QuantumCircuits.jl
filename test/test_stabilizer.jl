@testset "Pauli operators and stabilizer codes" begin
    @testset "Pauli algebra" begin
        X_, Y_, Z_ = PauliOp("X"), PauliOp("Y"), PauliOp("Z")
        @test X_ * Y_ == PauliOp("iZ") && Y_ * X_ == PauliOp("-iZ")
        @test Y_ * Z_ == PauliOp("iX") && Z_ * X_ == PauliOp("iY")
        @test string(PauliOp("-XIZ")) == "-XIZ" && weight(PauliOp("XIZY")) == 3
        @test -PauliOp("XY") == PauliOp("-XY")
        for _ in 1:100
            n = rand(1:3)
            a = PauliOp(rand(Bool, n), rand(Bool, n), rand(0:3))
            b = PauliOp(rand(Bool, n), rand(Bool, n), rand(0:3))
            @test matrix(a * b) ≈ matrix(a) * matrix(b)
            @test commutes(a, b) == (matrix(a) * matrix(b) ≈ matrix(b) * matrix(a))
            ψ = randn(ComplexF64, 2^n)
            @test a * ψ ≈ matrix(a) * ψ
            c = pauli!(Circuit(n), a)
            @test matrix(c) ≈ matrix(a)
        end
        @test_throws ArgumentError PauliOp("XQ")
        @test_throws ArgumentError commutes(PauliOp("X"), PauliOp("XX"))
    end

    @testset "code catalogue" begin
        for (code, n, k, d) in ((repetition_code(3), 3, 1, 1), (five_qubit_code(), 5, 1, 3),
                                (shor_code(), 9, 1, 3), (steane_code(), 7, 1, 3),
                                (quantum_reed_muller(3), 7, 1, 3), (quantum_reed_muller(4), 15, 1, 3),
                                (rotated_surface_code(3), 9, 1, 3), (rotated_surface_code(5), 25, 1, 5))
            @test (length(code), dimension(code)) == (n, k)
            @test code_distance(code) == d
            X̄, Z̄ = logical_operators(code)
            @test !commutes(X̄[1], Z̄[1])
            @test all(commutes(l, s) for l in vcat(X̄, Z̄) for s in stabilizers(code))
        end
        @test !is_css(five_qubit_code()) && is_css(steane_code()) && is_css(rotated_surface_code(3))
        # the Steane code is the m = 3 quantum Reed–Muller code
        st, qr = steane_code(), quantum_reed_muller(3)
        S1 = QuantumCircuits._symp_matrix(stabilizers(st)); S2 = QuantumCircuits._symp_matrix(stabilizers(qr))
        @test gf2_rank(vcat(S1, S2)) == gf2_rank(S1) == gf2_rank(S2)
        # surface-code checks: (d-1)² of weight 4, 2(d-1) of weight 2
        for d in (3, 5)
            ws = weight.(stabilizers(rotated_surface_code(d)))
            @test count(==(4), ws) == (d - 1)^2 && count(==(2), ws) == 2(d - 1)
        end
    end

    @testset "logicals computed for any code" begin
        # strip the given logicals and let Gram–Schmidt find a valid pair
        five = five_qubit_code()
        c = StabilizerCode(stabilizers(five))
        X̄, Z̄ = logical_operators(c)
        @test length(X̄) == 1 && !commutes(X̄[1], Z̄[1])
        # a two-logical-qubit CSS code
        c = css_code([1 1 1 1 0 0; 0 0 1 1 1 1], [1 1 1 1 0 0; 0 0 1 1 1 1])
        @test dimension(c) == 2
        X̄, Z̄ = logical_operators(c)
        @test [commutes(X̄[i], Z̄[j]) for i in 1:2, j in 1:2] == [false true; true false]
    end

    @testset "validation" begin
        @test_throws ArgumentError StabilizerCode(["XX", "ZI"])            # anticommute
        @test_throws ArgumentError StabilizerCode(["ZZ", "ZZ"])            # dependent
        @test_throws ArgumentError StabilizerCode(["iZZ"])                 # not Hermitian
        @test_throws ArgumentError css_code([1 1 0], [0 1 1])              # Hx Hzᵀ ≠ 0
        @test_throws ArgumentError StabilizerCode(["ZZ"]; logical_x=["XI"], logical_z=["ZI"])
        @test_throws ArgumentError rotated_surface_code(4)
    end

    @testset "syndromes and decoding" begin
        for code in (steane_code(), five_qubit_code(), shor_code(), rotated_surface_code(3))
            n = length(code)
            table = lookup_decoder(code)
            for q in 1:n, P in ("X", "Y", "Z")
                E = PauliOp(join(j == q ? P : "I" for j in 1:n))
                s = syndrome(code, E)
                @test any(s)                                   # every single error is seen
                C = table[s]
                @test weight(C) <= 1
                # C·E commutes with everything and is no logical error
                CE = C * E
                @test !any(syndrome(code, CE))
                @test QuantumCircuits._gf2_inrowspace(QuantumCircuits._symp_matrix(stabilizers(code)),
                                                      QuantumCircuits._symp(CE))
            end
        end
    end
end
