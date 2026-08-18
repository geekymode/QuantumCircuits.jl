# Run me:  julia --project=. examples/demo.jl
using QuantumCircuits
using LinearAlgebra

hdr(s) = println("\n", "="^70, "\n", s, "\n", "="^70)

hdr("1. The Gray code: consecutive entries differ in ONE bit")
println(" i   binary   gray    bit flipped")
for i in 0:7
    flip = i == 0 ? "-" : string(gray_flip_position(i, 3))
    println(lpad(i, 2), "   ", bits(i, 3), "     ", bits(gray(i), 3), "     ", flip)
end
println("\nfull cyclic walk for k=3: ", gray_flip_positions(3))
println("the walk closes (all CNOTs cancel): ",
        foldl(⊻, (1 << p for p in gray_flip_positions(3))) == 0)

hdr("2. A Bell state")
c = Circuit(2)
push!(c, H(), 1)
push!(c, CNOT(), 1, 2)
display(c)
println("\nstatevector: ", round.(statevector(c), digits=4))

hdr("3. Uniformly controlled rotation: 2^k CNOTs, no multi-controlled gates")
α = [0.1, 0.5, -1.2, 2.0]          # rotate target by α[j+1] when controls read j
c = multiplexed_rz(α, [1, 2], 3)
display(c)
println("\nCNOTs: ", count_cnots(c), "  (= 2^2, flat in the number of branches)")
println("branch angles recovered from the circuit's block-diagonal:")
U = matrix(c)
for j in 0:3
    blk = U[2j+1:2j+2, 2j+1:2j+2]
    println("  j=", bits(j, 2), "  α = ", round(2 * angle(blk[2, 2]), digits=4),
            "   (asked for ", α[j+1], ")")
end

hdr("4. Diagonal unitary in 2^n - 2 CNOTs")
φ = [0.0, 0.3, 1.1, -0.7, 2.2, 0.9, -1.5, 0.4]
c = diagonal(φ)
println(c, "   CNOTs = ", count_cnots(c), "  (2^3 - 2 = 6)")
println("exact match incl. global phase: ", matrix(c) ≈ Diagonal(cis.(φ)))

hdr("5. Arbitrary state preparation from |00>")
a = [1, 2, 3im, 4]
c = prepare_state(a)
display(c)
println("\nprepared: ", round.(statevector(c), digits=4))
println("target:   ", round.(normalize(ComplexF64.(a)), digits=4))
println("exact:    ", statevector(c) ≈ normalize(ComplexF64.(a)))

hdr("6. Cost scaling vs. the naive compilation")
println(" n  | CNOTs: state prep (2^(n+1)-4) | naive O(n*2^n)")
for n in 1:8
    circ = prepare_state(normalize(randn(ComplexF64, 1 << n)))
    println(lpad(n, 2), "  | ", lpad(count_cnots(circ), 27), " | ", lpad(n * (1 << n), 14))
end

hdr("7. Matrix decompositions")
U = matrix(prepare_state([1, 2, 3im, 4]))
factors = two_level_decompose(U)
println("two-level (Givens) factors: ", length(factors))
println("  product reproduces U:      ", reduce(*, matrix(t, 4) for t in factors) ≈ U)
println("  synthesised circuit ≈ U:   ", matrix(synthesize_unitary(U)) ≈ U)
α, β, γ, δ = zyz(matrix(H()))
println("ZYZ of H: α=", round(α, digits=3), " β=", round(β, digits=3),
        " γ=", round(γ, digits=3), " δ=", round(δ, digits=3))
V, θ, W = demultiplex(matrix(RY(0.4)), matrix(RZ(1.1)))
mid = Diagonal(vcat(cis.(-θ ./ 2), cis.(θ ./ 2)))
println("demultiplex identity holds: ",
        kron(I(2), V) * mid * kron(I(2), W) ≈ cat(matrix(RY(0.4)), matrix(RZ(1.1)); dims=(1, 2)))

hdr("8. Phase polynomials: Gray ordering vs one gadget per term")
println(" n  | terms | one gadget each | parity network (gray) | optimal 2^n-2")
for n in 2:6
    local pp = phase_polynomial(randn(1 << n))
    local naive = synthesize(pp; order=:gadgets)
    local gray = synthesize(pp; order=:gray)
    println(lpad(n, 2), "  | ", lpad(nterms(pp), 5), " | ", lpad(count_cnots(naive), 15),
            " | ", lpad(count_cnots(gray), 21), " | ", lpad((1 << n) - 2, 13))
end

hdr("9. Hamiltonian simulation from a matrix")
Hm = randn(4, 4); Hm = Hm + Hm'
terms = pauli_decompose(Hm)
println("Pauli terms: ", join([string(s, "→", round(real(c), digits=3)) for (s, c) in terms], "  "))
c = Circuit(2); trotter_step!(c, terms, 0.05)
println("Trotter step: ", length(c), " gates, ", count_cnots(c), " CNOTs")
println("fidelity vs exp(-iH·dt): ", round(gate_fidelity(matrix(c), exp(-im * 0.05 * Hm)), digits=6))

hdr("10. Illustrations")
println("Plotting is a package extension. With a Makie backend installed:")
println("    using CairoMakie")
println("    circuitfigure(prepare_state([1, 2, 3im, 4]))")
println("    graycodefigure(4);  matrixfigure(U; part=:phase);  costfigure(2:9)")
println("    save(\"circuit.pdf\", circuitfigure(c))")
