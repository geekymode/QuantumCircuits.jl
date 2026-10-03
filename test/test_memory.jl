@testset "memory experiment" begin
    QC = QuantumCircuits
    rng = Random.Xoshiro(23)
    code = rotated_surface_code(5)
    zch = [i for (i, s) in enumerate(stabilizers(code)) if any(s.z)]
    g = QC._spacetime_graph(code, 5)

    @testset "detectors" begin
        for _ in 1:5                                      # noiseless: nothing fires
            det, lg = QC._sample_memory(code, NoiseModel(), 5, rng)
            @test !any(det) && !lg
        end
        for q in 1:25                                     # one data X error in round 2
            E = PauliOp(BitVector(j == q for j in 1:25), falses(25))
            det, lg = QC._sample_memory(code, NoiseModel(), 5, rng; inject=(2, E))
            @test det[:, 2] == [stabilizers(code)[i].z[q] for i in zch]
            @test !any(det[:, [1, 3, 4, 5, 6]])
            @test QC._decode_memory(g, det) == lg           # decoded correctly
        end
        for (row, i) in enumerate(zch), τ in 1:4           # one measurement error
            det, lg = QC._sample_memory(code, NoiseModel(), 5, rng; flip=(τ, i))
            @test findall(det) == [CartesianIndex(row, τ), CartesianIndex(row, τ + 1)]
            @test !lg && QC._decode_memory(g, det) == false
        end
        # a Z error leaves the Z-check detectors (and Z̄) alone
        det, lg = QC._sample_memory(code, NoiseModel(), 5, rng; inject=(3, PauliOp("Z" * "I"^24)))
        @test !any(det) && !lg
    end

    @testset "noise models" begin
        @test circuit_noise(0.01) == NoiseModel(p1=0.01, p2=0.01, pmeas=0.01, pdata=0.01)
        @test phenomenological_noise(0.02) == NoiseModel(pmeas=0.02, pdata=0.02)
        @test memory_experiment(code, NoiseModel(); shots=20, rng=rng) == (0.0, 0)
        @test_throws ArgumentError memory_experiment(five_qubit_code(), circuit_noise(0.01))
        @test_throws ArgumentError memory_experiment(steane_code(), circuit_noise(0.01))
        @test_throws ArgumentError memory_experiment(code, circuit_noise(0.01); rounds=0)
    end

    @testset "below threshold, bigger is better" begin
        # phenomenological noise at 1%: measured ~0.008 (d=3) against ~0.003 (d=5)
        r3 = memory_experiment(rotated_surface_code(3), phenomenological_noise(0.01);
                               shots=5000, rng=Random.Xoshiro(1))[1]
        r5 = memory_experiment(code, phenomenological_noise(0.01); shots=5000, rng=Random.Xoshiro(1))[1]
        @test r5 < r3
    end
end
