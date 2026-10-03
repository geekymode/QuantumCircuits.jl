# ---------------------------------------------------------------------------
# Circuits for stabilizer codes
#
# Encoding a CSS code is a CNOT network — the same kind of object as the Gray
# encoder.  For checks Hx in reduced row echelon form, the encoded basis state
#
#     |x̄⟩ = Σ_{s ∈ rowspace(Hx)} |x·Lx + s⟩
#
# is built by copying each logical input along its X̄ row, then putting each
# pivot qubit of Hx in |+⟩ and copying it along its row.
#
# Syndrome extraction measures each stabilizer onto its own ancilla: Hadamard,
# the stabilizer controlled on the ancilla, Hadamard.  Circuits here stay
# unitary — the ancillas end in the syndrome basis state rather than being
# measured — so everything can be checked exactly against the statevector.
# ---------------------------------------------------------------------------

# CNOTs (and swaps) on `wires` implementing |x⟩ ↦ |x M⟩ for invertible M over
# GF(2), x a row vector: reduce M to I by column operations M·C₁⋯C_t = I, so
# x M = x C_t ⋯ C₁, and each C is its own inverse.
function _linear_map!(c::Circuit, M::AbstractMatrix, wires::Vector{Int})
    A = _gf2(M)
    k = size(A, 1)
    ops = Tuple{Symbol,Int,Int}[]
    for col in 1:k
        if !A[col, col]
            p = findfirst(j -> A[col, j], col+1:k)
            p === nothing && throw(ArgumentError("map is singular over GF(2)"))
            p += col
            A[:, [col, p]] = A[:, [p, col]]
            push!(ops, (:swap, col, p))
        end
        for j in 1:k
            if j != col && A[col, j]
                A[:, j] .⊻= A[:, col]                 # col_j += col_col
                push!(ops, (:cnot, col, j))           # x_j ⊕= x_col
            end
        end
    end
    for (kind, a, b) in reverse(ops)
        kind === :swap ? push!(c, SWAP(), wires[a], wires[b]) : push!(c, CNOT(), wires[a], wires[b])
    end
    c
end

"""
    encoding_circuit(code) -> (circuit, inputs)

An encoder for a CSS code with `X`-type logical `X̄` (as [`css_code`](@ref)
and the catalogue produce), made of Hadamards and CNOTs only.  Place logical
qubit `j` on wire `inputs[j]` and every other wire in `|0⟩`; the circuit maps
`|x⟩` to the encoded `|x̄⟩` of the code's own logical operators
([`logical_operators`](@ref)) — `Z̄ⱼ` reads back `xⱼ`, `X̄ⱼ` flips it.

Non-CSS codes (such as [`five_qubit_code`](@ref)) are not supported.
"""
function encoding_circuit(code::StabilizerCode)
    is_css(code) || throw(ArgumentError("encoding_circuit supports CSS codes only"))
    all(l -> !any(l.z), code.logical_x) ||
        throw(ArgumentError("encoding_circuit needs X-type logical X̄ operators, as css_code computes"))
    n, k = length(code), dimension(code)
    Xs = [s.x for s in code.stabilizers if any(s.x)]
    R, P = isempty(Xs) ? (falses(0, n), Int[]) : gf2_rref(reduce(vcat, permutedims.(Xs)))
    R = R[1:length(P), :]

    Lx = BitMatrix(reduce(vcat, permutedims(l.x) for l in code.logical_x; init=falses(0, n)))
    for j in 1:k, (i, p) in enumerate(P)                # clear the pivot columns
        Lx[j, p] && (Lx[j, :] .⊻= R[i, :])
    end
    # rref of Lx, tracking the row transformation A (rref = A·Lx)
    aug, Q = gf2_rref(hcat(Lx, BitMatrix(LinearAlgebra.I(k))))
    length(Q) == k && all(q -> q <= n, Q) || throw(ArgumentError("logical operators are not independent of the stabilizers"))
    Lr, A = aug[:, 1:n], aug[:, n+1:end]

    c = Circuit(n)
    # inputs x in the code's basis become x A⁻¹ in the reduced basis
    k > 1 && _linear_map!(c, _gf2_inv(A), Q)
    for j in 1:k, col in 1:n
        col != Q[j] && Lr[j, col] && push!(c, CNOT(), Q[j], col)
    end
    for (i, p) in enumerate(P)
        push!(c, H(), p)
        for col in 1:n
            col != p && R[i, col] && push!(c, CNOT(), p, col)
        end
    end
    c, Q
end

"""
    syndrome_circuit(code) -> Circuit

Syndrome extraction for any stabilizer code, on `n` data wires and one
ancilla per stabilizer (wires `n+1 …`, marked as ancillas).  For each
stabilizer: `H` on its ancilla, the stabilizer's Paulis controlled on the
ancilla, `H` again.  On an encoded state hit by a Pauli error `E`, the
ancillas finish in the basis state [`syndrome`](@ref)`(code, E)` and the data
in `E|ψ̄⟩`: nothing is learned about `|ψ̄⟩` itself.
"""
function syndrome_circuit(code::StabilizerCode)
    n, r = length(code), length(code.stabilizers)
    c = Circuit(n; ancillas=r)
    gates = Dict((true, false) => controlled(X()), (true, true) => controlled(Y()),
                 (false, true) => controlled(Z()))
    for (i, s) in enumerate(code.stabilizers)
        a = n + i
        push!(c, H(), a)
        for j in 1:n
            (s.x[j] || s.z[j]) && push!(c, gates[(s.x[j], s.z[j])], a, j)
        end
        s.phase == 2 && push!(c, Z(), a)            # -P: flip the outcome
        push!(c, H(), a)
    end
    c
end
