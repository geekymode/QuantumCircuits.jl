# ---------------------------------------------------------------------------
# Embedding any unitary in a real symmetric involution
#
# Synthesising arbitrary unitaries is no harder than synthesising real,
# symmetric, traceless ones with S² = I (Nehoran & Yuen, Lemma 2.1).  Two
# embeddings, each doubling the dimension and costing one ancilla:
#
#   realify             U  ↦  R_U = [Re U  -Im U; Im U  Re U]     real orthogonal
#   hermitian_dilation  R  ↦  S_R = [0  R; R†  0]                 involution, traceless
#
# R_U acts as U on |−Y⟩ ⊗ |ψ⟩ with |−Y⟩ = (|0⟩ - i|1⟩)/√2, and S_R maps |1⟩ ⊗ |ψ⟩
# to |0⟩ ⊗ R|ψ⟩ — so the second ancilla even resets itself.  A construction that
# only handles the restricted class therefore handles everything, at the cost
# of two wires and three one-qubit gates: see `embedding_circuit`.
# ---------------------------------------------------------------------------

"""
    realify(U) -> Matrix{Float64}

The real representation `[Re U  -Im U; Im U  Re U]` of an `N × N` matrix, of
size `2N`.  It is a homomorphism — `realify(U * V) == realify(U) * realify(V)`
— and maps unitaries to real orthogonal matrices.

On a qubit placed as the most significant bit, it acts as `U` on
`|−Y⟩ ⊗ |ψ⟩`, where `|−Y⟩ = (|0⟩ - i|1⟩)/√2 = RX(π/2)|0⟩`.
"""
realify(U::AbstractMatrix) = [real(U) -imag(U); imag(U) real(U)]

"""
    hermitian_dilation(R) -> Matrix

`[0 R; R† 0]`.  For unitary `R` this is a Hermitian, traceless involution
(`S² = I`); for real `R` it is real symmetric.

On a qubit placed as the most significant bit, it maps `|1⟩ ⊗ |ψ⟩` to
`|0⟩ ⊗ R|ψ⟩`: prepare the ancilla in `|1⟩` and it comes back in `|0⟩`.
"""
function hermitian_dilation(R::AbstractMatrix)
    size(R, 1) == size(R, 2) || throw(ArgumentError("expected a square matrix, got $(size(R))"))
    Z = zero(R)
    [Z R; R' Z]
end

"""
    symmetric_embedding(U) -> Matrix{Float64}

`hermitian_dilation(realify(U))`: a real, symmetric, traceless involution of
four times the dimension, from which `U` is recovered with two ancillas — see
[`embedding_circuit`](@ref).  This is the reduction of Nehoran and Yuen's
Lemma 2.1: synthesising this class is as hard as synthesising every unitary.
"""
symmetric_embedding(U::AbstractMatrix) = hermitian_dilation(realify(U))

"""
    embedding_circuit(U) -> Circuit
    embedding_circuit(sc) -> Circuit

A circuit implementing `U` on `n` data wires through its
[`symmetric_embedding`](@ref) `S`, with two ancillas returned clean:

1. `X` on the first ancilla (`|1⟩`, selecting the off-diagonal block of `S`),
   and `RX(π/2)` on the second (`|−Y⟩`, selecting `U` inside `realify(U)`);
2. `S` on (ancilla 1, ancilla 2, data), ancilla 1 being most significant;
3. `RX(-π/2)` on the second ancilla.  The first has already reset itself.

Given a matrix, `S` is one dense gate.  Given a circuit `sc` on `n + 2` wires
that implements `S` (wires 1 and 2 standing for the two ancillas), `sc` is
spliced in instead — so anything able to synthesise real symmetric involutions
synthesises `U`.  Check with [`implementation_error`](@ref).
"""
function embedding_circuit(U::AbstractMatrix)
    N = size(U, 1)
    size(U) == (N, N) && ispow2(N) ||
        throw(ArgumentError("expected a 2ⁿ × 2ⁿ matrix, got $(size(U))"))
    n = trailing_zeros(N)
    sc = Circuit(n + 2)
    push!(sc, Gate(:Sym, symmetric_embedding(U)), 1:n+2...)
    embedding_circuit(sc)
end

function embedding_circuit(sc::Circuit)
    n = nqubits(sc) - 2
    n >= 0 || throw(ArgumentError("an embedding circuit needs at least two wires"))
    a1, a2 = n + 1, n + 2
    wire(q) = q == 1 ? a1 : q == 2 ? a2 : q - 2
    c = Circuit(n; ancillas=2)
    push!(c, X(), a1)
    push!(c, RX(π / 2), a2)
    for op in sc.ops
        push!(c, op.gate, map(wire, op.qubits)...)
    end
    c.global_phase += sc.global_phase
    push!(c, RX(-π / 2), a2)
    c
end
