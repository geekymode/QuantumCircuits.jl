# ---------------------------------------------------------------------------
# Stabilizer (Clifford tableau) simulation
#
# A stabilizer state on n qubits is fixed by n commuting Paulis, so it fits in
# O(n²) bits instead of 2ⁿ amplitudes, and Clifford gates update those bits in
# O(n) time (Gottesman–Knill).  This is the Aaronson–Gottesman tableau: rows
# 1…n are destabilizers, n+1…2n stabilizers, 2n+1 scratch, each a packed Pauli
# with a sign bit.  Destabilizers make measurement O(n²) instead of O(n³).
#
# Circuits stay unitary everywhere else in the package; measurement lives
# here.  `apply!` runs any Clifford `Circuit` on a tableau.
# ---------------------------------------------------------------------------

"""
    Tableau(n)

A stabilizer-state simulator on `n` qubits, starting in `|0…0⟩`, after
Aaronson and Gottesman.  Memory and gate cost grow as `n²` and `n`, so
hundreds of qubits are cheap — but only Clifford operations are allowed.

Operations: [`apply!`](@ref) (gates or a whole `Circuit`), [`pauli!`](@ref),
[`measure!`](@ref), [`expectation`](@ref), [`stabilizers`](@ref).  Global
phase is not tracked.
"""
mutable struct Tableau
    n::Int
    W::Int                       # 64-bit words per row
    x::Matrix{UInt64}            # W × (2n+1): column = row of the tableau
    z::Matrix{UInt64}
    r::BitVector                 # sign bit per row (-1 when set)
end

function Tableau(n::Integer)
    n >= 1 || throw(ArgumentError("need at least one qubit"))
    W = cld(n, 64)
    t = Tableau(n, W, zeros(UInt64, W, 2n + 1), zeros(UInt64, W, 2n + 1), falses(2n + 1))
    for q in 1:n
        _setbit!(t.x, q, q, true)                 # destabilizer q = X_q
        _setbit!(t.z, q, n + q, true)             # stabilizer q = Z_q
    end
    t
end

nqubits(t::Tableau) = t.n

@inline _word(q) = (q - 1) >> 6 + 1
@inline _mask(q) = UInt64(1) << ((q - 1) & 63)
@inline _getbit(M, q, row) = (M[_word(q), row] & _mask(q)) != 0
@inline function _setbit!(M, q, row, v::Bool)
    w, m = _word(q), _mask(q)
    M[w, row] = v ? (M[w, row] | m) : (M[w, row] & ~m)
end

function _check_qubit(t::Tableau, q::Integer)
    1 <= q <= t.n || throw(ArgumentError("qubit $q out of range 1:$(t.n)"))
end

# --- the three generators of the Clifford group -----------------------------

function _h!(t::Tableau, a::Int)
    w, m = _word(a), _mask(a)
    @inbounds for i in 1:2t.n
        xa = t.x[w, i] & m; za = t.z[w, i] & m
        (xa != 0 && za != 0) && (t.r[i] = !t.r[i])
        t.x[w, i] = (t.x[w, i] & ~m) | za
        t.z[w, i] = (t.z[w, i] & ~m) | xa
    end
end

function _s!(t::Tableau, a::Int)
    w, m = _word(a), _mask(a)
    @inbounds for i in 1:2t.n
        xa = t.x[w, i] & m; za = t.z[w, i] & m
        (xa != 0 && za != 0) && (t.r[i] = !t.r[i])
        t.z[w, i] ⊻= xa
    end
end

function _cnot!(t::Tableau, a::Int, b::Int)
    wa, ma, wb, mb = _word(a), _mask(a), _word(b), _mask(b)
    @inbounds for i in 1:2t.n
        xa = (t.x[wa, i] & ma) != 0; za = (t.z[wa, i] & ma) != 0
        xb = (t.x[wb, i] & mb) != 0; zb = (t.z[wb, i] & mb) != 0
        (xa && zb && (xb == za)) && (t.r[i] = !t.r[i])
        xa && (t.x[wb, i] ⊻= mb)
        zb && (t.z[wa, i] ⊻= ma)
    end
end

# Pauli on one qubit: flips the sign of every row it anticommutes with.
function _pauli1!(t::Tableau, a::Int, px::Bool, pz::Bool)
    w, m = _word(a), _mask(a)
    @inbounds for i in 1:2t.n
        xa = (t.x[w, i] & m) != 0; za = (t.z[w, i] & m) != 0
        ((px && za) ⊻ (pz && xa)) && (t.r[i] = !t.r[i])
    end
end

# Exponent of i, mod 4, in (row i)·(row h) — +i for XY, YZ, ZX; -i reversed.
function _rowphase(t::Tableau, i::Int, h::Int)
    s = 0
    @inbounds for w in 1:t.W
        x1, z1, x2, z2 = t.x[w, i], t.z[w, i], t.x[w, h], t.z[w, h]
        plus = (x1 & ~z1 & x2 & z2) | (x1 & z1 & ~x2 & z2) | (~x1 & z1 & x2 & ~z2)
        minus = (x1 & z1 & x2 & ~z2) | (~x1 & z1 & x2 & z2) | (x1 & ~z1 & ~x2 & z2)
        s += count_ones(plus) - count_ones(minus)
    end
    s
end

# Row h ← row i · row h.
function _rowsum!(t::Tableau, h::Int, i::Int)
    ph = 2 * t.r[h] + 2 * t.r[i] + _rowphase(t, i, h)
    t.r[h] = mod(ph, 4) == 2
    @inbounds for w in 1:t.W
        t.x[w, h] ⊻= t.x[w, i]
        t.z[w, h] ⊻= t.z[w, i]
    end
end

# --- gates and circuits -------------------------------------------------------

# A rotation angle that is a whole number of quarter turns, as 0…3; or nothing.
function _quarter_turns(θ::Real)
    k = θ / (π / 2)
    isapprox(k, round(k); atol=1e-9) ? mod(round(Int, k), 4) : nothing
end

"""
    apply!(t::Tableau, g::Gate, qubits...) -> t
    apply!(t::Tableau, c::Circuit) -> t

Apply a Clifford gate, or every instruction of a circuit, to a tableau.
Accepted: `Id X Y Z H S Sdg CNOT CZ SWAP`, controlled `X Y Z`, and `RX RY RZ
PHASE` at multiples of `π/2`.  Anything else (a `T`, a generic rotation) is
not Clifford and throws.  Global phase is ignored.
"""
function apply!(t::Tableau, g::Gate, qubits::Integer...)
    qs = collect(Int, qubits)
    foreach(q -> _check_qubit(t, q), qs)
    nm = g.name
    if nm === :I
    elseif nm === :X; _pauli1!(t, qs[1], true, false)
    elseif nm === :Y; _pauli1!(t, qs[1], true, true)
    elseif nm === :Z; _pauli1!(t, qs[1], false, true)
    elseif nm === :H; _h!(t, qs[1])
    elseif nm === :S; _s!(t, qs[1])
    elseif nm === :Sdg; _s!(t, qs[1]); _s!(t, qs[1]); _s!(t, qs[1])
    elseif nm === :CNOT || nm === :CX; _cnot!(t, qs[1], qs[2])
    elseif nm === :CZ; _h!(t, qs[2]); _cnot!(t, qs[1], qs[2]); _h!(t, qs[2])
    elseif nm === :CY; _s!(t, qs[2]); _s!(t, qs[2]); _s!(t, qs[2])        # S†
                       _cnot!(t, qs[1], qs[2]); _s!(t, qs[2])               # S·CX·S† = CY
    elseif nm === :SWAP; _cnot!(t, qs[1], qs[2]); _cnot!(t, qs[2], qs[1]); _cnot!(t, qs[1], qs[2])
    elseif nm in (:RZ, :P, :RX, :RY)
        k = _quarter_turns(g.params[1])
        k === nothing && throw(ArgumentError("$(label(g)) is not a Clifford gate"))
        q = qs[1]
        nm === :RX && _h!(t, q)
        nm === :RY && (_s!(t, q); _s!(t, q); _s!(t, q); _h!(t, q))     # S† then H: Y → Z
        for _ in 1:k; _s!(t, q); end
        nm === :RX && _h!(t, q)
        nm === :RY && (_h!(t, q); _s!(t, q))
    else
        throw(ArgumentError("$(label(g)) is not a supported Clifford gate"))
    end
    t
end

function apply!(t::Tableau, c::Circuit)
    c.nqubits == t.n || throw(ArgumentError("circuit has $(c.nqubits) qubits, tableau $(t.n)"))
    for op in c.ops
        apply!(t, op.gate, op.qubits...)
    end
    t
end

"""
    pauli!(t::Tableau, p::PauliOp, qubits=1:length(p)) -> t

Apply a Pauli (an error, or a correction) to a tableau.
"""
function pauli!(t::Tableau, p::PauliOp, qubits::AbstractVector{<:Integer}=1:length(p))
    length(qubits) == length(p) || throw(ArgumentError("need $(length(p)) qubits"))
    for (j, q) in enumerate(qubits)
        (p.x[j] || p.z[j]) && _pauli1!(t, Int(q), p.x[j], p.z[j])
    end
    t
end

# --- Pauli rows ↔ PauliOp -----------------------------------------------------

function _row_pauli(t::Tableau, row::Int)
    PauliOp([_getbit(t.x, q, row) for q in 1:t.n], [_getbit(t.z, q, row) for q in 1:t.n],
            t.r[row] ? 2 : 0)
end

function _load_row!(t::Tableau, row::Int, p::PauliOp)
    t.x[:, row] .= 0; t.z[:, row] .= 0
    for q in 1:t.n
        p.x[q] && _setbit!(t.x, q, row, true)
        p.z[q] && _setbit!(t.z, q, row, true)
    end
    t.r[row] = false
end

# Whether the Pauli in `row` anticommutes with p (packed into a scratch row).
function _anticommutes_row(t::Tableau, row::Int, prow::Int)
    s = 0
    @inbounds for w in 1:t.W
        s += count_ones(t.x[w, row] & t.z[w, prow]) + count_ones(t.z[w, row] & t.x[w, prow])
    end
    isodd(s)
end

"""
    stabilizers(t::Tableau) -> Vector{PauliOp}

The current stabilizer generators, with signs: the state is their joint
`+1` eigenstate.
"""
stabilizers(t::Tableau) = [_row_pauli(t, t.n + i) for i in 1:t.n]

# --- measurement --------------------------------------------------------------

function _check_pauli(t::Tableau, P::PauliOp)
    length(P) == t.n || throw(ArgumentError("Pauli acts on $(length(P)) qubits, tableau has $(t.n)"))
    iseven(P.phase) || throw(ArgumentError("can only measure Hermitian Paulis (phase ±1)"))
end

"""
    measure!(t, q; rng) -> Bool
    measure!(t, P::PauliOp; rng) -> Bool

Measure qubit `q` in the computational basis, or any Hermitian Pauli `P`, and
collapse the state.  Returns `true` for outcome `1` (eigenvalue `-1`).  The
outcome is random exactly when the state is not an eigenstate — when `P`
anticommutes with some stabilizer — and then each value has probability ½.
"""
function measure!(t::Tableau, P::PauliOp; rng=Random.default_rng())
    _check_pauli(t, P)
    n, scr = t.n, 2t.n + 1
    _load_row!(t, scr, P)
    p = findfirst(i -> _anticommutes_row(t, n + i, scr), 1:n)
    if p !== nothing                                       # random outcome
        p += n
        for i in 1:2n
            i != p && _anticommutes_row(t, i, scr) && _rowsum!(t, i, p)
        end
        t.x[:, p-n] .= t.x[:, p]; t.z[:, p-n] .= t.z[:, p]; t.r[p-n] = t.r[p]
        t.x[:, p] .= t.x[:, scr]; t.z[:, p] .= t.z[:, scr]
        outcome = rand(rng, Bool)
        t.r[p] = outcome ⊻ (P.phase == 2)                  # the row is (-1)^outcome · P
        return outcome
    end
    # deterministic: P = ± product of the stabilizers whose destabilizer anticommutes
    t.x[:, scr] .= 0; t.z[:, scr] .= 0; t.r[scr] = false
    for i in 1:n
        _anticommutes_with(t, i, P) && _rowsum!(t, scr, n + i)
    end
    t.r[scr] ⊻ (P.phase == 2)
end

measure!(t::Tableau, q::Integer; rng=Random.default_rng()) =
    (_check_qubit(t, q); measure!(t, PauliOp(falses(t.n), BitVector(j == q for j in 1:t.n)); rng=rng))

function _anticommutes_with(t::Tableau, row::Int, P::PauliOp)
    s = 0
    for q in 1:t.n
        s += (_getbit(t.x, q, row) & P.z[q]) + (_getbit(t.z, q, row) & P.x[q])
    end
    isodd(s)
end

"""
    expectation(t, P::PauliOp) -> Int

`⟨P⟩` for a Hermitian Pauli on a stabilizer state: `+1` or `-1` when the
state is an eigenstate of `P`, otherwise `0`.  Does not disturb the state.
"""
function expectation(t::Tableau, P::PauliOp)
    _check_pauli(t, P)
    n = t.n
    any(i -> _anticommutes_with(t, n + i, P), 1:n) && return 0
    scr = 2n + 1
    t.x[:, scr] .= 0; t.z[:, scr] .= 0; t.r[scr] = false
    for i in 1:n
        _anticommutes_with(t, i, P) && _rowsum!(t, scr, n + i)
    end
    (t.r[scr] ⊻ (P.phase == 2)) ? -1 : 1
end

"""
    statevector(t::Tableau) -> Vector{ComplexF64}

The state as a dense vector, up to global phase — for checking against the
statevector simulator on small `n`.
"""
function statevector(t::Tableau)
    t.n <= 16 || throw(ArgumentError("statevector(::Tableau) is for n ≤ 16"))
    S = stabilizers(t)
    # project a basis state with nonzero overlap onto the joint +1 eigenspace
    for b in 0:(1 << t.n)-1
        ψ = zeros(ComplexF64, 1 << t.n); ψ[b+1] = 1
        for s in S
            ψ = (ψ .+ s * ψ) ./ 2
        end
        nrm = LinearAlgebra.norm(ψ)
        nrm > 1e-6 && return ψ ./ nrm
    end
    error("empty stabilizer space — inconsistent tableau")
end
