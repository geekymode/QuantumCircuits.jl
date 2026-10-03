@testset "detector error models" begin
    QC = QuantumCircuits
    code = rotated_surface_code(3)
    rounds = 3
    n, r = 9, 8
    N = n + r
    steps = QC._memory_steps(code, rounds)
    nrec = rounds * r + n

    # Ground truth: run the circuit on the tableau, insert Pauli P after step s,
    # and return which measurement records differ from the fault-free run.
    function records(s, P; rng)
        t = Tableau(N)
        out = falses(nrec)
        s == 0 && P !== nothing && pauli!(t, P)
        for (k, st) in enumerate(steps)
            if st.kind === :gate
                apply!(t, st.gate, st.qubits...)
            elseif st.kind === :meas
                a = st.qubits[1]
                out[st.record] = measure!(t, a; rng=rng)
                out[st.record] && apply!(t, X(), a)
            elseif st.kind === :final
                out[st.record] = measure!(t, st.qubits[1]; rng=rng)
            end
            k == s && P !== nothing && pauli!(t, P)
        end
        out
    end
    # Individual X-check records and final data bits are random; detectors (Z-check
    # records, and Z-check parities of the final data) and Z̄ are not.  Compare those.
    zch = [i for (i, s) in enumerate(stabilizers(code)) if any(s.z)]
    stabs = stabilizers(code)
    zl = logical_operators(code)[2][1].z
    function observables(rec)
        dets = Bool[]
        for i in zch, τ in 1:rounds+1
            push!(dets, τ <= rounds ? rec[(τ-1)*r+i] ⊻ (τ > 1 && rec[(τ-2)*r+i]) :
                        isodd(count(q -> stabs[i].z[q] && rec[rounds*r+q], 1:n)) ⊻ rec[(rounds-1)*r+i])
        end
        push!(dets, isodd(count(q -> zl[q] && rec[rounds*r+q], 1:n)))
        dets
    end
    clean = observables(records(0, nothing; rng=Random.Xoshiro(1)))

    @testset "frame propagation = tableau, every fault location" begin
        rng = Random.Xoshiro(5)
        checked = 0
        for (s, st) in enumerate(steps)
            st.kind in (:noise1, :noise2) || continue
            for q in st.qubits, (xb, zb) in ((true, false), (false, true), (true, true))
                P = PauliOp(BitVector(j == q && xb for j in 1:N), BitVector(j == q && zb for j in 1:N))
                truth = observables(records(s, P; rng=rng)) .⊻ clean
                # the fault-free observables are all zero, so the prediction is
                # the observables of the propagated record flips
                pred = observables(QC._propagate(steps, s, copy(P.x), copy(P.z), nrec))
                @test truth == pred
                checked += 1
            end
        end
        @test checked > 500
    end

    @testset "model structure" begin
        for noise in (phenomenological_noise(0.01), circuit_noise(0.01))
            dem = detector_error_model(code, noise; rounds=rounds)
            @test dem.ndetectors == 4 * (rounds + 1)
            @test all(m -> 1 <= length(m[1]) <= 2, dem.mechanisms)    # graph-like, nothing empty
            @test all(m -> 0 < m[3] < 0.5, dem.mechanisms)
            @test dem.dropped == 0
        end
        # gate faults add diagonal edges the uniform graph lacks
        @test length(detector_error_model(code, circuit_noise(0.01); rounds=rounds).mechanisms) >
              length(detector_error_model(code, phenomenological_noise(0.01); rounds=rounds).mechanisms)
        # no noise, no mechanisms
        @test isempty(detector_error_model(code, NoiseModel(); rounds=rounds).mechanisms)
        @test_throws ArgumentError detector_error_model(five_qubit_code(), circuit_noise(0.01))
    end

    @testset "decoding" begin
        dem = detector_error_model(rotated_surface_code(5), circuit_noise(0.004); rounds=5)
        g = QC._dem_graph(dem)
        # every single fault mechanism is decoded correctly on its own
        m = 12
        for (dets, lg, _) in dem.mechanisms
            det = falses(m, 6)
            for d in dets
                det[(d - 1) % m + 1, (d - 1) ÷ m + 1] = true
            end
            @test QC._decode_dem(g, det, m) == lg
        end
        # and it beats the uniform graph at circuit-level noise (measured
        # ~0.0045 against ~0.022 at d = 5, p = 0.004)
        c5 = rotated_surface_code(5)
        rd = memory_experiment(c5, circuit_noise(0.004); shots=2000, decoder=:dem, rng=Random.Xoshiro(1))[1]
        ru = memory_experiment(c5, circuit_noise(0.004); shots=2000, decoder=:uniform, rng=Random.Xoshiro(1))[1]
        @test rd < ru / 2
        @test_throws ArgumentError memory_experiment(c5, circuit_noise(0.004); decoder=:magic)
    end
end
