# ---------------------------------------------------------------------------
# Gates
#
# A gate is a name, its parameters, and its dense matrix.  Multi-qubit matrices
# use the big-endian convention: for a gate acting on `qubits = [q1, q2, …]`,
# `q1` is the most significant bit of the matrix index.
# ---------------------------------------------------------------------------

"""
    Gate(name, params, mat)

A named unitary.  `mat` is `2^k × 2^k` for a `k`-qubit gate.
"""
struct Gate
    name::Symbol
    params::Vector{Float64}
    mat::Matrix{ComplexF64}
end

Gate(name::Symbol, mat::AbstractMatrix) = Gate(name, Float64[], ComplexF64.(mat))

nqubits(g::Gate) = trailing_zeros(size(g.mat, 1))

"""
    matrix(g::Gate)
    matrix(c::Circuit)

The dense unitary of a gate, or of a whole circuit.
"""
matrix(g::Gate) = g.mat

const _SQ2 = 1 / sqrt(2)

"""
    nqubits(g::Gate)
    nqubits(c::Circuit)

Number of qubits a gate acts on, or the width of a circuit's register.
"""
function nqubits end

"""
    Id()

The one-qubit identity, `[1 0; 0 1]`.
"""
Id()   = Gate(:I, ComplexF64[1 0; 0 1])

"""
    X()

Pauli-X (bit flip), `[0 1; 1 0]`.
"""
X()    = Gate(:X, ComplexF64[0 1; 1 0])

"""
    Y()

Pauli-Y, `[0 -im; im 0]`.
"""
Y()    = Gate(:Y, ComplexF64[0 -im; im 0])

"""
    Z()

Pauli-Z (phase flip), `[1 0; 0 -1]`.
"""
Z()    = Gate(:Z, ComplexF64[1 0; 0 -1])

"""
    H()

Hadamard, `[1 1; 1 -1]/√2`.
"""
H()    = Gate(:H, ComplexF64[_SQ2 _SQ2; _SQ2 -_SQ2])

"""
    S()

Phase gate `diag(1, im)`; `S == PHASE(π/2)`.
"""
S()    = Gate(:S, ComplexF64[1 0; 0 im])

"""
    Sdg()

Adjoint of [`S`](@ref), `diag(1, -im)`.
"""
Sdg()  = Gate(:Sdg, ComplexF64[1 0; 0 -im])

"""
    T()

`π/8` gate `diag(1, exp(iπ/4))`; `T^2 == S`.
"""
T()    = Gate(:T, ComplexF64[1 0; 0 cis(pi / 4)])

"""
    Tdg()

Adjoint of [`T`](@ref), `diag(1, exp(-iπ/4))`.
"""
Tdg()  = Gate(:Tdg, ComplexF64[1 0; 0 cis(-pi / 4)])

"""
    RX(θ)

Rotation `exp(-i θ X / 2)` about the x-axis.  See also [`RY`](@ref), [`RZ`](@ref).
"""
RX(θ::Real) = Gate(:RX, [float(θ)], ComplexF64[cos(θ/2) -im*sin(θ/2); -im*sin(θ/2) cos(θ/2)])

"""
    RY(θ)

Rotation `exp(-i θ Y / 2)` about the y-axis.

`X RY(θ) X == RY(-θ)` — the sign flip that the Gray-code multiplexor
([`multiplexed_rotation!`](@ref)) exploits.
"""
RY(θ::Real) = Gate(:RY, [float(θ)], ComplexF64[cos(θ/2) -sin(θ/2); sin(θ/2) cos(θ/2)])

"""
    RZ(θ)

Rotation `exp(-i θ Z / 2) = diag(exp(-iθ/2), exp(iθ/2))` about the z-axis.

`X RZ(θ) X == RZ(-θ)` — the sign flip that the Gray-code multiplexor
([`multiplexed_rotation!`](@ref)) exploits.
"""
RZ(θ::Real) = Gate(:RZ, [float(θ)], ComplexF64[cis(-θ/2) 0; 0 cis(θ/2)])

"""
    PHASE(λ)

`diag(1, exp(im*λ))`.
"""
PHASE(λ::Real) = Gate(:P, [float(λ)], ComplexF64[1 0; 0 cis(λ)])

"""
    CNOT()

Controlled-NOT.  Used as `push!(c, CNOT(), control, target)`; the first wire is
the control (the most significant bit of the matrix index).
"""
CNOT() = Gate(:CNOT, ComplexF64[1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0])

"""
    CZ()

Controlled-Z, `diag(1, 1, 1, -1)`.  Symmetric in its two wires.
"""
CZ()   = Gate(:CZ,   ComplexF64[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1])

"""
    SWAP()

Exchange the states of two qubits.
"""
SWAP() = Gate(:SWAP, ComplexF64[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1])

# FANOUT and PARITY are stored densely like every other gate, so their matrix
# grows as 4^(r+1): r = 10 is already 64 MB.  Past that, build the gate from
# CNOTs (`fanout!` / `parity!` with `style=:ladder` or `:tree`) instead.
const _MAX_DENSE_FANOUT = 10

function _check_fanout_size(r::Integer, what::String)
    r >= 1 || throw(ArgumentError("$what needs at least one wire besides the control/target"))
    r <= _MAX_DENSE_FANOUT || throw(ArgumentError(
        "a $what gate on $(r + 1) wires would need a dense $(1 << (r + 1))² matrix; " *
        "the limit is $(_MAX_DENSE_FANOUT + 1) wires — use style=:tree to build it from CNOTs"))
    nothing
end

"""
    FANOUT(r)

Quantum fan-out on `r + 1` wires: `|b, t₁ … t_r⟩ ↦ |b, t₁⊕b … t_r⊕b⟩`.  The
first wire is the control.  Equal to `r` CNOTs sharing a control, but kept as
one gate so that [`depth`](@ref) can count it as a single step — the
unbounded fan-out model of constant-depth circuits.  See [`fanout!`](@ref).

Stored as a dense matrix like every gate, so `r` is capped at 10 (an 11-wire
gate, 64 MB).
"""
function FANOUT(r::Integer)
    _check_fanout_size(r, "fan-out")
    D = 1 << (r + 1)
    flip = (1 << r) - 1
    m = zeros(ComplexF64, D, D)
    for i in 0:D-1
        j = i >> r == 1 ? i ⊻ flip : i
        m[j+1, i+1] = 1
    end
    Gate(:FANOUT, m)
end

"""
    PARITY(r)

Parity on `r + 1` wires: `|x₁ … x_r, t⟩ ↦ |x₁ … x_r, t ⊕ x₁ ⊕ … ⊕ x_r⟩`.  The
last wire is the target.  It is [`FANOUT`](@ref) conjugated by Hadamards on
every wire, with the control becoming the target.  See [`parity!`](@ref).
Capped at `r = 10`, like [`FANOUT`](@ref).
"""
function PARITY(r::Integer)
    _check_fanout_size(r, "parity")
    D = 1 << (r + 1)
    m = zeros(ComplexF64, D, D)
    for i in 0:D-1
        j = i ⊻ (count_ones(i >> 1) & 1)
        m[j+1, i+1] = 1
    end
    Gate(:PARITY, m)
end

"""
    controlled(g, n=1) -> Gate

Add `n` control qubits to `g`.  In the resulting matrix the controls are the
most significant bits, so `controlled(X())` is `CNOT()`.
"""
function controlled(g::Gate, n::Integer=1)
    d = size(g.mat, 1)
    D = d * (1 << n)
    m = Matrix{ComplexF64}(LinearAlgebra.I, D, D)
    m[D-d+1:D, D-d+1:D] = g.mat
    Gate(Symbol("C"^n, g.name), g.params, m)
end

"""
    label(g) -> String

Short display label, e.g. `"RZ(1.57)"`.
"""
function label(g::Gate)
    isempty(g.params) && return string(g.name)
    string(g.name, "(", join((@sprintf("%.2f", p) for p in g.params), ","), ")")
end

Base.show(io::IO, g::Gate) = print(io, label(g))
