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
    Circuit(n)

An `n`-qubit circuit: an ordered list of [`Instruction`](@ref)s plus an
accumulated `global_phase` (in radians).  Global phase is unobservable on its
own but is tracked so that [`matrix`](@ref) and [`statevector`](@ref) are exact
rather than correct-up-to-phase, which makes tests much sharper.
"""
mutable struct Circuit
    nqubits::Int
    ops::Vector{Instruction}
    global_phase::Float64
end
Circuit(n::Integer) = Circuit(Int(n), Instruction[], 0.0)

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
    c
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
"""
function statevector(c::Circuit, ψ0::AbstractVector=zero_state(c.nqubits))
    ψ = ComplexF64.(collect(ψ0))
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

Print an ASCII diagram of the circuit, one column per instruction.
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
        elseif g.name === :SWAP
            labels[q[1]] = "×"; labels[q[2]] = "×"
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
        print(io, "q", lpad(r, pad), ": ─")
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
    println(io, "Circuit(", c.nqubits, " qubits, ", length(c.ops), " gates, ",
            count_cnots(c), " CNOTs)")
    length(c.ops) <= 64 && draw(io, c)
end
Base.show(io::IO, c::Circuit) = print(io, "Circuit(", c.nqubits, " qubits, ", length(c.ops), " gates)")
