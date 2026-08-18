# Plotting lives in a package extension.  These tests run only when a Makie
# backend is installed:  QC_TEST_PLOTS=true julia --project=... test/runtests.jl
# CI exercises the same code path when the documentation builds its figures.
@testset "plotting extension" begin
    if get(ENV, "QC_TEST_PLOTS", "false") != "true"
        @info "skipping plot tests (set QC_TEST_PLOTS=true with CairoMakie installed to run them)"
        @test true
    else
        @eval using CairoMakie
        @test Base.get_extension(QuantumCircuits, :QuantumCircuitsMakieExt) !== nothing

        circuits = [multiplexed_rz([0.1, 0.5, -1.2, 2.0], [1, 2], 3),
                    prepare_state([1, 2, 3im, 4]),
                    diagonal(randn(8)),
                    synthesize_unitary(randu(4)),
                    Circuit(2)]
        for c in circuits, theme in (:light, :dark)
            @test circuitfigure(c; theme=theme) isa Makie.Figure
        end
        U = matrix(multiplexed_1q([randu(2) for _ in 1:4], [1, 2], 3))
        for part in (:abs, :real, :imag, :phase)
            @test matrixfigure(U; part=part) isa Makie.Figure
        end
        @test matrixfigure(circuits[1]) isa Makie.Figure
        @test graycodefigure(3) isa Makie.Figure
        @test costfigure(2:6) isa Makie.Figure
        @test_throws ArgumentError matrixfigure(U; part=:nope)
        @test_throws ArgumentError circuitfigure(circuits[1]; theme=:neon)
    end
end

@testset "plotting stubs without a backend" begin
    # the error must name the fix, not just fail
    if get(ENV, "QC_TEST_PLOTS", "false") != "true"
        for f in (circuitfigure, matrixfigure, graycodefigure, costfigure)
            e = try; f(Circuit(1)); catch err; err; end
            @test e isa ErrorException
            @test occursin("CairoMakie", e.msg)
        end
    else
        @test true
    end
end
