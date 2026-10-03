# ---------------------------------------------------------------------------
# Circuits and a dense statevector simulator
#
# Convention: qubit 1 is the MOST significant bit of a basis index.  Basis
# state |q1 q2 … qn⟩ has index  Σ_q  bit_q * 2^(n-q).  This is the "reading
# order" convention: the top wire in a drawing is the leftmost bit in a ket.
# ---------------------------------------------------------------------------

"""
    Instruction(gate, qubits)

A gate together with the (1-based) wires it acts on.  `qubits[1]` maps to the
most significant index bit of `gate.mat`.
"""
struct Instruction
    gate::Gate
    qubits::Vector{Int}
end

"""
    Circuit(n; ancillas=0)

An `n`-qubit circuit: an ordered list of [`Instruction`](@ref)s plus an
accumulated `global_phase` (in radians).  Global phase is unobservable on its
own but is tracked so that [`matrix`](@ref) and [`statevector`](@ref) are exact
rather than correct-up-to-phase, which makes tests much sharper.

`ancillas` adds that many wires after the `n` data wires, promised to start in
`|0⟩`.  More can be added later with [`add_ancillas!`](@ref).  Ancillas are
ordinary wires to every gate; they only change what [`isometry`](@ref),
[`logical_matrix`](@ref) and [`leakage`](@ref) report.
"""
mutable struct Circuit
    nqubits::Int
    ops::Vector{Instruction}
    global_phase::Float64
    ancillas::Vector{Int}
end
Circuit(n::Integer, ops::Vector{Instruction}, global_phase::Real) =
    Circuit(Int(n), ops, Float64(global_phase), Int[])
function Circuit(n::Integer; ancillas::Integer=0)
    ancillas >= 0 || throw(ArgumentError("ancilla count must be non-negative"))
    Circuit(Int(n + ancillas), Instruction[], 0.0, collect(n+1:n+ancillas))
end

nqubits(c::Circuit) = c.nqubits
Base.length(c::Circuit) = length(c.ops)

"""
    push!(c, gate, qubits...)

Append `gate` acting on `qubits`.
"""
function Base.push!(c::Circuit, g::Gate, qubits::Integer...)
    qs = collect(Int, qubits)
    length(qs) == nqubits(g) || throw(ArgumentError("$(label(g)) acts on $(nqubits(g)) qubit(s), got $(length(qs))"))
    allunique(qs) || throw(ArgumentError("repeated qubit in $qs"))
    all(1 .<= qs .<= c.nqubits) || throw(ArgumentError("qubit index out of range 1:$(c.nqubits) in $qs"))
    push!(c.ops, Instruction(g, qs))
    c
end

"""
    append!(c, other)

Append every instruction of `other` (and its global phase) onto `c`.
"""
function Base.append!(c::Circuit, other::Circuit)
    c.nqubits == other.nqubits || throw(ArgumentError("qubit count mismatch"))
    append!(c.ops, other.ops)
    c.global_phase += other.global_phase
    union!(c.ancillas, other.ancillas)
    sort!(c.ancillas)
    c
end

# --- ancillas --------------------------------------------------------------

"""
    add_ancillas!(c, k=1) -> Vector{Int}

Append `k` fresh wires to `c`, marked as ancillas that start in `|0⟩`, and
return their indices.  Existing wire numbers never change, so this is safe in
the middle of building a circuit.
"""
function add_ancillas!(c::Circuit, k::Integer=1)
    k >= 0 || throw(ArgumentError("ancilla count must be non-negative"))
    new = collect(c.nqubits+1:c.nqubits+k)
    c.nqubits += k
    append!(c.ancillas, new)
    new
end

"""
    ancillas(c) -> Vector{Int}

The wires of `c` marked as ancillas, in increasing order.
"""
ancillas(c::Circuit) = copy(c.ancillas)

"""
    data_qubits(c) -> Vector{Int}

The wires of `c` that are not ancillas, in increasing order.  The first is the
most significant bit of the index of [`logical_matrix`](@ref) and
[`isometry`](@ref).
"""
data_qubits(c::Circuit) = setdiff(1:c.nqubits, c.ancillas)

# Full basis index of data value x (MSB = first data wire) with ancillas |0⟩.
function _data_indices(c::Circuit)
    dq = data_qubits(c)
    d, n = length(dq), c.nqubits
    [sum((((x >> (d - i)) & 1) << (n - q) for (i, q) in enumerate(dq)); init=0)
     for x in 0:(1 << d)-1]
end

"""
    isometry(c) -> Matrix{ComplexF64}

`C (I ⊗ |0⟩)`: the circuit's action on data inputs with every ancilla in
`|0⟩`, a `2ⁿ × 2ᵈ` isometry for `n` wires and `d` data wires.  Column `x+1` is
the output for data input `|x⟩`; rows use the full wire ordering.

Without ancillas this is [`matrix`](@ref).
"""
function isometry(c::Circuit)
    idx = _data_indices(c)
    N = 1 << c.nqubits
    V = Matrix{ComplexF64}(undef, N, length(idx))
    ψ = Vector{ComplexF64}(undef, N)
    for (j, i) in enumerate(idx)
        fill!(ψ, 0)
        ψ[i+1] = 1
        for op in c.ops
            apply!(ψ, op.gate.mat, op.qubits, c.nqubits)
        end
        @views V[:, j] .= ψ
    end
    c.global_phase == 0 ? V : V .* cis(c.global_phase)
end

"""
    logical_matrix(c) -> Matrix{ComplexF64}

`(I ⊗ ⟨0|) C (I ⊗ |0⟩)`: the `2ᵈ × 2ᵈ` block of the circuit seen from the data
wires, with ancillas prepared in `|0⟩` and post-selected on `|0⟩`.  When the
ancillas come back clean ([`is_clean`](@ref)) this is the unitary the circuit
implements; otherwise it is a contraction, and the missing weight is the
[`leakage`](@ref).
"""
logical_matrix(c::Circuit) = isometry(c)[_data_indices(c) .+ 1, :]

"""
    leakage(c) -> Float64

How far the ancillas can end away from `|0⟩`, over all data inputs: the
operator norm of the part of [`isometry`](@ref) outside the ancilla-`|0⟩`
subspace.  Zero exactly when every ancilla is returned clean.
"""
function leakage(c::Circuit)
    V = isometry(c)
    V[_data_indices(c) .+ 1, :] .= 0
    opnorm(V)
end

"""
    is_clean(c; atol=1e-9) -> Bool

Whether every ancilla of `c` returns to `|0⟩` for every data input, i.e.
`leakage(c) ≤ atol`.
"""
is_clean(c::Circuit; atol::Real=1e-9) = leakage(c) <= atol

"""
    implementation_error(c, U) -> Float64

`‖C (I ⊗ |0⟩) - U ⊗ |0⟩‖` in operator norm: how far the circuit, run with
clean ancillas, is from applying `U` to the data wires *and* returning the
ancillas clean.  This is the error measure of Nehoran and Yuen's Theorem 5.1.
Exact in global phase, like everything else here.
"""
function implementation_error(c::Circuit, U::AbstractMatrix)
    idx = _data_indices(c)
    size(U) == (length(idx), length(idx)) ||
        throw(ArgumentError("U is $(size(U)), the circuit has $(length(idx)) data states"))
    V = isometry(c)
    V[idx .+ 1, :] .-= U
    opnorm(V)
end

"""
    count_gates(c, name) -> Int

Number of instructions whose gate has the given `name`, e.g.
`count_gates(c, :CNOT)`.  Handy for checking synthesis cost.
"""
count_gates(c::Circuit, name::Symbol) = count(op -> op.gate.name === name, c.ops)

"""
    count_cnots(c) -> Int

Number of CNOTs in the circuit — the usual proxy for cost on hardware, and the
figure of merit every decomposition in this package is trying to minimise.
"""
count_cnots(c::Circuit) = count_gates(c, :CNOT)

"""
    depth(c) -> Int
    depth(c, name) -> Int

Number of time steps when every gate runs as early as its wires allow (ASAP
scheduling), i.e. the longest chain of gates sharing wires.

With a gate `name`, only those gates cost a step and the rest are free:
`depth(c, :CNOT)` is the CNOT depth.  A [`FANOUT`](@ref) or [`PARITY`](@ref)
gate counts as one step however many wires it touches.
"""
function depth(c::Circuit, name::Union{Symbol,Nothing}=nothing)
    level = zeros(Int, c.nqubits)
    for op in c.ops
        l = maximum(level[q] for q in op.qubits) + (name === nothing || op.gate.name === name)
        level[op.qubits] .= l
    end
    maximum(level; init=0)
end

"""
    layers(c) -> Vector{Vector{Instruction}}

The ASAP schedule behind [`depth`](@ref): layer `t` holds the instructions
that run at step `t`.  Instructions within a layer act on disjoint wires, so
they commute and can run in parallel.
"""
function layers(c::Circuit)
    level = zeros(Int, c.nqubits)
    out = Vector{Instruction}[]
    for op in c.ops
        l = maximum(level[q] for q in op.qubits) + 1
        level[op.qubits] .= l
        l > length(out) && push!(out, Instruction[])
        push!(out[l], op)
    end
    out
end

# --- simulation ------------------------------------------------------------

# Index of the basis state obtained from `base` by setting the gate-local
# pattern `a` (MSB = qubits[1]) into the positions given by `shifts`.
@inline function _scatter_index(base::Int, a::Int, shifts::Vector{Int}, k::Int)
    idx = base
    @inbounds for i in 1:k
        if (a >> (k - i)) & 1 == 1
            idx |= (1 << shifts[i])
        end
    end
    idx
end

"""
    apply!(ψ, U, qubits, n) -> ψ

Apply the `2^k × 2^k` matrix `U` to wires `qubits` of the `n`-qubit
statevector `ψ`, in place.
"""
function apply!(ψ::AbstractVector{ComplexF64}, U::AbstractMatrix{ComplexF64}, qubits::Vector{Int}, n::Int)
    k = length(qubits)
    shifts = [n - q for q in qubits]
    mask = 0
    for s in shifts
        mask |= (1 << s)
    end
    buf = Vector{ComplexF64}(undef, 1 << k)
    out = Vector{ComplexF64}(undef, 1 << k)
    idxs = Vector{Int}(undef, 1 << k)
    for base in 0:(1 << n)-1
        (base & mask) == 0 || continue
        @inbounds for a in 0:(1 << k)-1
            idxs[a+1] = _scatter_index(base, a, shifts, k)
            buf[a+1] = ψ[idxs[a+1]+1]
        end
        LinearAlgebra.mul!(out, U, buf)
        @inbounds for a in 0:(1 << k)-1
            ψ[idxs[a+1]+1] = out[a+1]
        end
    end
    ψ
end

"""
    statevector(c, ψ0 = |0…0⟩) -> Vector{ComplexF64}

Run the circuit on `ψ0` and return the resulting statevector, including the
tracked global phase.

If the circuit has ancillas, `ψ0` may instead be a state of the data wires
alone; the ancillas are then prepared in `|0⟩`.
"""
function statevector(c::Circuit, ψ0::AbstractVector=zero_state(c.nqubits))
    ψ = ComplexF64.(collect(ψ0))
    if !isempty(c.ancillas) && length(ψ) == 1 << (c.nqubits - length(c.ancillas))
        full = zeros(ComplexF64, 1 << c.nqubits)
        full[_data_indices(c) .+ 1] .= ψ
        ψ = full
    end
    length(ψ) == 1 << c.nqubits || throw(ArgumentError("state has wrong length"))
    for op in c.ops
        apply!(ψ, op.gate.mat, op.qubits, c.nqubits)
    end
    c.global_phase == 0 ? ψ : ψ .* cis(c.global_phase)
end

"""
    zero_state(n) -> Vector{ComplexF64}

The `|0…0⟩` statevector on `n` qubits.
"""
function zero_state(n::Integer)
    ψ = zeros(ComplexF64, 1 << n)
    ψ[1] = 1
    ψ
end

"""
    matrix(c) -> Matrix{ComplexF64}

The full `2^n × 2^n` unitary of the circuit (column `j` is the circuit applied
to basis state `j-1`).  Exponential in `n` — for inspection and testing.
"""
function matrix(c::Circuit)
    N = 1 << c.nqubits
    U = Matrix{ComplexF64}(undef, N, N)
    ψ = Vector{ComplexF64}(undef, N)
    for j in 1:N
        fill!(ψ, 0)
        ψ[j] = 1
        for op in c.ops
            apply!(ψ, op.gate.mat, op.qubits, c.nqubits)
        end
        @views U[:, j] .= ψ
    end
    c.global_phase == 0 ? U : U .* cis(c.global_phase)
end

# --- drawing ---------------------------------------------------------------

function _center(s::AbstractString, w::Int, fill::Char='─')
    pad = w - length(s)
    left = pad ÷ 2
    string(repeat(fill, left), s, repeat(fill, pad - left))
end

"""
    draw([io], c)

Print an ASCII diagram of the circuit, one column per instruction.  Ancilla
wires are labelled `a` instead of `q`, keeping their wire number.
"""
function draw(io::IO, c::Circuit)
    n = c.nqubits
    cols = [String[] for _ in 1:n]
    for op in c.ops
        labels = fill("", n)
        q = op.qubits
        g = op.gate
        if g.name === :CNOT
            labels[q[1]] = "●"; labels[q[2]] = "⊕"
        elseif g.name === :CZ
            labels[q[1]] = "●"; labels[q[2]] = "●"
        elseif g.name === :CCX
            labels[q[1]] = "●"; labels[q[2]] = "●"; labels[q[3]] = "⊕"
        elseif g.name === :SWAP
            labels[q[1]] = "×"; labels[q[2]] = "×"
        elseif g.name === :FANOUT
            labels[q[1]] = "●"
            for qq in q[2:end]; labels[qq] = "⊕"; end
        elseif g.name === :PARITY
            for qq in q[1:end-1]; labels[qq] = "●"; end
            labels[q[end]] = "⊕"
        elseif length(q) == 1
            labels[q[1]] = label(g)
        else
            for (i, qq) in enumerate(q)
                labels[qq] = string(label(g), "[", i, "]")
            end
        end
        w = maximum(length, labels)
        lo, hi = extrema(q)
        for r in 1:n
            if !isempty(labels[r])
                push!(cols[r], _center(labels[r], w))
            elseif lo < r < hi
                push!(cols[r], _center("│", w))
            else
                push!(cols[r], repeat("─", w))
            end
        end
    end
    pad = length(string(n))
    for r in 1:n
        print(io, r in c.ancillas ? "a" : "q", lpad(r, pad), ": ─")
        print(io, join(cols[r], "─"))
        println(io, "─")
    end
    if c.global_phase != 0
        println(io, "global phase: ", @sprintf("%.4f", c.global_phase), " rad")
    end
    nothing
end
draw(c::Circuit) = draw(stdout, c)

function Base.show(io::IO, ::MIME"text/plain", c::Circuit)
    println(io, "Circuit(", _width_string(c), ", ", length(c.ops), " gates, ",
            count_cnots(c), " CNOTs)")
    length(c.ops) <= 64 && draw(io, c)
end
Base.show(io::IO, c::Circuit) = print(io, "Circuit(", _width_string(c), ", ", length(c.ops), " gates)")

function _width_string(c::Circuit)
    a = length(c.ancillas)
    a == 0 ? string(c.nqubits, " qubits") :
             string(c.nqubits - a, " qubits + ", a, a == 1 ? " ancilla" : " ancillas")
end
