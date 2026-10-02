@testset "ancillas" begin
    @testset "bookkeeping" begin
        c = Circuit(3)
        @test isempty(ancillas(c)) && data_qubits(c) == [1, 2, 3]
        @test isometry(c) ≈ matrix(c)                 # no ancillas: nothing changes
        @test logical_matrix(c) ≈ matrix(c)

        c = Circuit(2; ancillas=2)
        @test nqubits(c) == 4
        @test ancillas(c) == [3, 4] && data_qubits(c) == [1, 2]
        @test add_ancillas!(c, 3) == [5, 6, 7]
        @test nqubits(c) == 7 && ancillas(c) == [3, 4, 5, 6, 7]
        @test add_ancillas!(c, 0) == Int[]
        @test_throws ArgumentError Circuit(2; ancillas=-1)
        @test_throws ArgumentError add_ancillas!(c, -1)

        @test sprint(show, Circuit(2; ancillas=1)) == "Circuit(2 qubits + 1 ancilla, 0 gates)"
        @test sprint(show, Circuit(2; ancillas=3)) == "Circuit(2 qubits + 3 ancillas, 0 gates)"
        @test sprint(show, Circuit(2)) == "Circuit(2 qubits, 0 gates)"
        s = sprint(draw, push!(Circuit(1; ancillas=1), CNOT(), 1, 2))
        @test occursin("q1:", s) && occursin("a2:", s)

        a = Circuit(3); b = Circuit(3); b.ancillas = [3]
        append!(a, b)
        @test ancillas(a) == [3]
    end

    @testset "isometry and logical matrix" begin
        # ancilla in the MIDDLE of the register: data wires 1 and 3
        c = Circuit(3); c.ancillas = [2]
        for (g, qs) in ((H(), (1,)), (CNOT(), (1, 3)), (RY(0.4), (3,)), (CZ(), (3, 1)))
            push!(c, g, qs...)
        end
        ref = Circuit(2)
        for (g, qs) in ((H(), (1,)), (CNOT(), (1, 2)), (RY(0.4), (2,)), (CZ(), (2, 1)))
            push!(ref, g, qs...)
        end
        @test logical_matrix(c) ≈ matrix(ref)
        @test implementation_error(c, matrix(ref)) < 1e-12
        @test is_clean(c)
        V = isometry(c)
        @test size(V) == (8, 4) && V' * V ≈ I(4)

        # data-sized input to statevector
        ψ = normalize!(randn(ComplexF64, 4))
        out = statevector(c, ψ)
        @test length(out) == 8 && out[[1, 2, 5, 6]] ≈ matrix(ref) * ψ

        # a dirty ancilla: on data |1⟩ it is rotated away from |0⟩ by θ
        θ = 0.9
        d = Circuit(1; ancillas=1); push!(d, controlled(RY(θ)), 1, 2)
        @test leakage(d) ≈ sin(θ / 2)
        @test !is_clean(d)
        @test logical_matrix(d) ≈ Diagonal([1, cos(θ / 2)])    # a contraction
        @test implementation_error(d, Matrix(I(2))) ≈ 2sin(θ / 4)
        # copying data out and never uncomputing leaks completely on |−⟩
        e = Circuit(1; ancillas=1); push!(e, H(), 1); push!(e, CNOT(), 1, 2)
        @test leakage(e) ≈ 1
        push!(e, CNOT(), 1, 2)                                  # uncompute
        @test is_clean(e) && implementation_error(e, matrix(H())) < 1e-12

        @test_throws ArgumentError implementation_error(d, Matrix(I(4)))
    end

    @testset "AND tree" begin
        for k in 1:7
            c = Circuit(k + 1); and!(c, 1:k, k + 1)
            ref = multicontrolled(matrix(X()), 1:k, k + 1)
            @test implementation_error(c, matrix(ref)) < 1e-10
            @test is_clean(c)
            @test length(ancillas(c)) == max(k - 2, 0)
            if k >= 2
                @test count_gates(c, :CCX) == 2k - 3
                @test depth(c) == 2clog2(k) - 1
            end
        end
        # lowered Toffolis implement the same thing
        c = Circuit(5); and!(c, 1:4, 5; lower=true)
        @test count_gates(c, :CCX) == 0
        @test implementation_error(c, matrix(multicontrolled(matrix(X()), 1:4, 5))) < 1e-10
        # reusing clean ancillas: two ANDs, one ancilla pool
        c = Circuit(6; ancillas=2)
        pool = ancillas(c)
        and!(c, 1:4, 5; ancillas=pool); and!(c, [2, 3, 4, 5], 6; ancillas=pool)
        @test nqubits(c) == 8 && is_clean(c)
        @test_throws ArgumentError and!(Circuit(3), Int[], 3)
        @test_throws ArgumentError and!(Circuit(3), [1, 2], 2)
        @test_throws ArgumentError and!(Circuit(5), 1:4, 5; ancillas=[6])
        @test_throws ArgumentError and!(Circuit(5), 1:3, 5; ancillas=[3])
        c = Circuit(3); try and!(c, [1, 3], 3) catch end
        @test nqubits(c) == 3                          # no stray wires on failure
    end
end
