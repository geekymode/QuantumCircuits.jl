# ---------------------------------------------------------------------------
# Two-qubit unitaries: the KAK (Cartan) decomposition
#
# Every two-qubit unitary factors as
#
#     U = e^{iφ} (A₁ ⊗ A₂) · exp(i(a XX + b YY + c ZZ)) · (B₁ ⊗ B₂)
#
# — local gates around a *canonical gate* fixed by three angles.  The local
# parts are free; the entangling content is (a, b, c), and the canonical gate
# needs at most three CNOTs.  So does every two-qubit unitary, which is optimal.
#
# How: in the magic basis, SU(2) ⊗ SU(2) becomes the real group SO(4), and XX,
# YY, ZZ become diagonal.  For U′ = M† U M, the matrix U′ᵀU′ is symmetric and
# unitary; a real orthogonal P diagonalises it, and U′ = K · D^{1/2} · Pᵀ with
# K real orthogonal.  K and Pᵀ map back to local gates; D^{1/2} is the
# canonical gate.
# ---------------------------------------------------------------------------

# The magic basis: columns are Bell states with phases chosen so that local
# gates become real orthogonal matrices.
const _MAGIC = ComplexF64[1 0 0 im; 0 im 1 0; 0 im -1 0; 1 0 0 -im] ./ sqrt(2)

# Signs of XX, YY, ZZ on the magic basis vectors: in that basis
# exp(i(aXX + bYY + cZZ)) = diag(exp(i·_MAGIC_SIGNS·(a, b, c))).
const _MAGIC_SIGNS = [1 -1 1; 1 1 -1; -1 -1 -1; -1 1 1]

"""
    canonical_gate(a, b, c) -> Matrix{ComplexF64}

`exp(i(a XX + b YY + c ZZ))`: the entangling core of every two-qubit unitary.
The three terms commute, so the angles are independent; each is defined
modulo `π/2` up to local gates.  `canonical_gate(π/4, 0, 0)` is locally
equivalent to a CNOT.
"""
function canonical_gate(a::Real, b::Real, c::Real)
    d = cis.(_MAGIC_SIGNS * [a, b, c])
    _MAGIC * LinearAlgebra.Diagonal(d) * _MAGIC'
end

"""
    KAK

The factors of a two-qubit unitary,
`U = e^{i·phase} (A1 ⊗ A2) · canonical_gate(a, b, c) · (B1 ⊗ B2)`,
with `A1, A2, B1, B2 ∈ SU(2)` up to sign and each angle in `(-π/4, π/4]`.
`Matrix(k)` multiplies it back out.  Produced by [`kak`](@ref).
"""
struct KAK
    A1::Matrix{ComplexF64}
    A2::Matrix{ComplexF64}
    B1::Matrix{ComplexF64}
    B2::Matrix{ComplexF64}
    a::Float64
    b::Float64
    c::Float64
    phase::Float64
end

Base.Matrix(k::KAK) = cis(k.phase) .* (kron(k.A1, k.A2) * canonical_gate(k.a, k.b, k.c) *
                                       kron(k.B1, k.B2))

# Split a 4×4 Kronecker product A = A₁ ⊗ A₂ (unitary up to scale): take the
# largest 2×2 block, which is A₁[p,q]·A₂, normalise it to SU(2), and read A₁
# off by projecting every block onto it.
function _kron_factor(A::AbstractMatrix)
    blocks = [A[2i-1:2i, 2j-1:2j] for i in 1:2, j in 1:2]
    _, k = findmax(LinearAlgebra.norm, blocks)
    A2 = blocks[k] ./ sqrt(LinearAlgebra.det(blocks[k]))
    A1 = [LinearAlgebra.tr(A2' * blocks[i, j]) / 2 for i in 1:2, j in 1:2]
    A1, A2
end

"""
    kak(U; atol=1e-9) -> KAK

The KAK (Cartan) decomposition of a two-qubit unitary:

    U = e^{iφ} (A₁ ⊗ A₂) · exp(i(a XX + b YY + c ZZ)) · (B₁ ⊗ B₂)

with each angle reduced to `(-π/4, π/4]` — a whole `π/2` turn of `XX` is
`i·X⊗X`, a local gate, so it moves into `B₁ ⊗ B₂`.  Exact, global phase
included.  [`two_qubit!`](@ref) turns it into at most three CNOTs.
"""
function kak(U::AbstractMatrix; atol::Real=1e-9)
    size(U) == (4, 4) || throw(ArgumentError("kak expects a 4×4 matrix, got $(size(U))"))
    Uc = ComplexF64.(U)
    is_unitary(Uc; atol=1e-8) || throw(ArgumentError("kak expects a unitary matrix"))

    φ = angle(LinearAlgebra.det(Uc)) / 4
    Um = _MAGIC' * (Uc .* cis(-φ)) * _MAGIC            # SU(4), magic basis
    M2 = transpose(Um) * Um                            # symmetric unitary

    # Re and Im of M2 commute (M2 is unitary and symmetric), so one real
    # orthogonal P diagonalises both; a generic real combination finds it.
    # Within a degenerate eigenspace both are scalar, so any basis works.
    P = LinearAlgebra.eigen(LinearAlgebra.Symmetric(real(M2) .+ 0.4142135623730951 .* imag(M2))).vectors
    LinearAlgebra.det(P) < 0 && (P[:, 1] .*= -1)       # P ∈ SO(4)
    d = sqrt.(LinearAlgebra.diag(transpose(P) * M2 * P))
    K = Um * P * LinearAlgebra.Diagonal(1 ./ d)        # real orthogonal
    if real(LinearAlgebra.det(K)) < 0                  # move to SO(4)
        d[1] = -d[1]
        K = Um * P * LinearAlgebra.Diagonal(1 ./ d)
    end

    # angle(d) = _MAGIC_SIGNS·(a, b, c) + g
    sol = [_MAGIC_SIGNS ones(4)] \ angle.(d)
    abc, g = sol[1:3], sol[4]
    A1, A2 = _kron_factor(_MAGIC * real(K) * _MAGIC')
    B1, B2 = _kron_factor(_MAGIC * transpose(P) * _MAGIC')
    phase = φ + g

    # reduce each angle to (-π/4, π/4]: exp(i(θ + kπ/2)P⊗P) = exp(iθ P⊗P)·(i P⊗P)^k
    paulis = (ComplexF64[0 1; 1 0], ComplexF64[0 -im; im 0], ComplexF64[1 0; 0 -1])
    for (j, σ) in enumerate(paulis)
        k = -floor(Int, (π / 4 - abc[j]) / (π / 2))    # abc[j] - kπ/2 ∈ (-π/4, π/4]
        abc[j] -= k * π / 2
        if isodd(k)
            B1 = σ * B1
            B2 = σ * B2
        end
        phase += k * π / 2                             # i^k
    end
    clean(x) = abs(x) < atol ? 0.0 : x
    KAK(A1, A2, B1, B2, clean(abc[1]), clean(abc[2]), clean(abc[3]), phase)
end

"""
    two_qubit!(c, U, qubits; atol=1e-9) -> c

Append an arbitrary two-qubit unitary on `qubits = [q1, q2]` (`q1` most
significant), through its [`kak`](@ref) decomposition, using as few CNOTs as
the canonical angles allow:

| canonical angles `(a, b, c)` | CNOTs |
|:--|:--|
| all zero — a product of one-qubit gates | 0 |
| two zero, the third `±π/4` — locally a CNOT | 1 |
| one zero | 2 |
| otherwise | 3 |

Each count is optimal for its class; three is what a generic two-qubit
unitary needs.  Exact, global phase included.
"""
function two_qubit!(c::Circuit, U::AbstractMatrix, qubits::AbstractVector{<:Integer}; atol::Real=1e-9)
    qs = collect(Int, qubits)
    length(qs) == 2 && qs[1] != qs[2] || throw(ArgumentError("two_qubit! needs two distinct qubits, got $qs"))
    q1, q2 = qs
    k = kak(U; atol=atol)
    c.global_phase += k.phase
    a, b, cc = k.a, k.b, k.c
    zero_ = (abs(a) < atol, abs(b) < atol, abs(cc) < atol)
    I2 = Matrix{ComplexF64}(LinearAlgebra.I, 2, 2)

    if all(zero_)                                      # local
        decompose_1q!(c, k.A1 * k.B1, q1)
        decompose_1q!(c, k.A2 * k.B2, q2)
    elseif count(zero_) == 2 && isapprox(abs(a + b + cc), π / 4; atol=atol)
        # CNOT class: Can = L·exp(isπ/4·ZZ)·L†, and exp(iπ/4·ZZ) is
        # e^{-iπ/4}·CZ·(RZ(-π/2) ⊗ RZ(-π/2)), with CZ = (I⊗H)·CNOT·(I⊗H)
        θ = a + b + cc
        s = sign(θ)
        L = !zero_[3] ? I2 : !zero_[1] ? matrix(H()) : matrix(RX(π / 2))
        Zfix = s > 0 ? I2 : ComplexF64[1 0; 0 -1]       # exp(-iπ/4 ZZ) = exp(iπ/4 ZZ)·(-i Z⊗Z)
        Hm, Rm = matrix(H()), matrix(RZ(-π / 2))
        decompose_1q!(c, Rm * Zfix * L' * k.B1, q1)
        decompose_1q!(c, Hm * Rm * Zfix * L' * k.B2, q2)
        push!(c, CNOT(), q1, q2)
        decompose_1q!(c, k.A1 * L, q1)
        decompose_1q!(c, k.A2 * L * Hm, q2)
        c.global_phase += -π / 4 + (s > 0 ? 0.0 : -π / 2)
    elseif any(zero_)                                  # Can = L·Can(x, 0, z)·L†
        if zero_[2]
            x, z, L = a, cc, I2
        elseif zero_[1]                                # S⊗S swaps the XX, YY angles
            x, z, L = b, cc, matrix(S())
        else                                           # RX(π/2)⊗2 swaps YY, ZZ
            x, z, L = a, b, matrix(RX(π / 2))
        end
        decompose_1q!(c, L' * k.B1, q1)
        decompose_1q!(c, L' * k.B2, q2)
        push!(c, CNOT(), q1, q2)                       # CNOT·exp(i(x X₁ + z Z₂))·CNOT
        push!(c, RX(-2x), q1)
        push!(c, RZ(-2z), q2)
        push!(c, CNOT(), q1, q2)
        decompose_1q!(c, k.A1 * L, q1)
        decompose_1q!(c, k.A2 * L, q2)
    else                                               # Vatan–Williams, 3 CNOTs
        decompose_1q!(c, k.B1, q1)
        decompose_1q!(c, matrix(RZ(π / 2)) * k.B2, q2)
        push!(c, CNOT(), q2, q1)
        push!(c, RZ(-2cc - π / 2), q1)
        push!(c, RY(-2a - π / 2), q2)
        push!(c, CNOT(), q1, q2)
        push!(c, RY(2b + π / 2), q2)
        push!(c, CNOT(), q2, q1)
        decompose_1q!(c, k.A1 * matrix(RZ(-π / 2)), q1)
        decompose_1q!(c, k.A2, q2)
        c.global_phase -= π / 4                        # the circuit is e^{iπ/4}·Can
    end
    c
end

"""
    two_qubit(U) -> Circuit

Standalone circuit for [`two_qubit!`](@ref).
"""
two_qubit(U::AbstractMatrix; kwargs...) = two_qubit!(Circuit(2), U, [1, 2]; kwargs...)
