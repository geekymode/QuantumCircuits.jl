# ---------------------------------------------------------------------------
# Phase polynomials, phase gadgets, Pauli exponentials
#
# A diagonal unitary is the exponential of a real function on bitstrings, and
# any such function has a unique expansion over parity functions:
#
#     φ(x) = Σ_S  a_S · (-1)^(S · x)        S ⊆ {1 … n}
#
# — the *phase polynomial*.  Each term is one `exp(i a_S Z_S)`, and each of
# those is a *phase gadget*: a CNOT ladder that computes the parity of `S` onto
# one wire, an `RZ`, and the ladder undone.  Because a CNOT ladder is the
# expensive part, ordering the terms so that consecutive supports differ in one
# bit — Gray order — lets neighbouring ladders cancel.
# ---------------------------------------------------------------------------

"""
    phase_gadget!(c, θ, qubits) -> c

Append `exp(-im*θ/2 * Z⊗Z⊗…)` over `qubits`: a CNOT ladder collecting the
parity onto the last wire, `RZ(θ)` there, and the ladder reversed.

Uses `2(k-1)` CNOTs for `k` qubits — and none at all for `k == 1`, where it is
just a rotation.
"""
function phase_gadget!(c::Circuit, θ::Real, qubits::AbstractVector{<:Integer})
    qs = collect(Int, qubits)
    isempty(qs) && throw(ArgumentError("a phase gadget needs at least one qubit"))
    allunique(qs) || throw(ArgumentError("repeated qubit in $qs"))
    for i in 1:length(qs)-1
        push!(c, CNOT(), qs[i], qs[i+1])
    end
    push!(c, RZ(θ), qs[end])
    for i in length(qs)-1:-1:1
        push!(c, CNOT(), qs[i], qs[i+1])
    end
    c
end

"""
    phase_gadget(θ, qubits; n) -> Circuit

Standalone circuit for [`phase_gadget!`](@ref).
"""
phase_gadget(θ::Real, qubits::AbstractVector{<:Integer}; n::Integer=maximum(qubits)) =
    phase_gadget!(Circuit(n), θ, qubits)

# basis changes mapping each Pauli onto Z: U P U† = Z, applied as U … U†
_basis_in(::Val{'X'}) = (H(),)
_basis_in(::Val{'Y'}) = (Sdg(), H())          # H·S† conjugates Y to Z
_basis_out(::Val{'X'}) = (H(),)
_basis_out(::Val{'Y'}) = (H(), S())

"""
    pauli_rotation!(c, θ, s; qubits=1:length(s)) -> c

Append `exp(-im*θ/2 * P)` for the Pauli string `s` (e.g. `"XIZY"`) acting on
`qubits`.

Identity factors are skipped; `X` and `Y` factors are conjugated onto the `Z`
axis, and what remains is one [`phase_gadget!`](@ref). This is the elementary
term of a Trotter step — see [`trotter_step!`](@ref).
"""
function pauli_rotation!(c::Circuit, θ::Real, s::AbstractString;
                         qubits::AbstractVector{<:Integer}=1:length(s))
    qs = collect(Int, qubits)
    length(qs) == length(s) || throw(ArgumentError("Pauli string \"$s\" needs $(length(s)) qubits, got $(length(qs))"))
    support = Int[]
    pre = Tuple{Gate,Int}[]
    post = Tuple{Gate,Int}[]
    for (ch, q) in zip(s, qs)
        ch == 'I' && continue
        ch in ('X', 'Y', 'Z') || throw(ArgumentError("bad Pauli character '$ch'"))
        push!(support, q)
        if ch != 'Z'
            for g in _basis_in(Val(ch)); push!(pre, (g, q)); end
            for g in _basis_out(Val(ch)); push!(post, (g, q)); end
        end
    end
    if isempty(support)                      # the identity string is a global phase
        c.global_phase -= θ / 2
        return c
    end
    for (g, q) in pre;  push!(c, g, q); end
    phase_gadget!(c, θ, support)
    for (g, q) in post; push!(c, g, q); end
    c
end

"""
    trotter_step!(c, terms, dt; qubits=nothing, parity_network=true) -> c

One first-order Trotter step `∏_P exp(-im*dt*c_P*P)` for a Hamiltonian given as
Pauli terms — exactly the format [`pauli_decompose`](@ref) returns.

`parity_network` chooses between two orderings, and the choice is not just
cosmetic:

* `false` — emit every term in the order given, one gadget or rotation each.
  The circuit is exactly `∏` over the input list, in that order.
* `true` (default) — pool the diagonal (`I`/`Z`-only) terms, which all commute
  with one another, and synthesise them together through the Gray-ordered
  parity network before the rest. Pooling commuting terms is a *different*
  ordering from the input list — a standard and generally better one — so the
  circuit is a different (equally valid) first-order step, not a cheaper
  rendering of the same product.

How much the parity network saves depends on how dense the diagonal part is.
A Hamiltonian whose `Z` terms are a few low-weight couplings (transverse-field
Ising, say) sees no change — a two-body gadget is already two CNOTs. A dense
diagonal — a QAOA cost operator, a Hamiltonian with high-weight `Z` terms —
drops from `Σ 2(|S|-1)` to about `2ⁿ`. Measure with [`count_cnots`](@ref).

Terms carrying `X` or `Y` need their own basis change and are emitted
individually either way.

```julia
H = randn(4, 4); H = H + H'                 # a two-qubit Hamiltonian
terms = pauli_decompose(H)
c = Circuit(2); trotter_step!(c, terms, 0.05)
```
"""
function trotter_step!(c::Circuit, terms, dt::Real; qubits=nothing,
                       parity_network::Bool=true)
    isempty(terms) && return c
    len = length(first(first(terms)))
    qs = qubits === nothing ? collect(1:len) : collect(Int, qubits)
    length(qs) == len || throw(ArgumentError("Pauli strings of length $len need $len qubits"))
    _isdiag(s) = all(ch -> ch == 'I' || ch == 'Z', s) && any(ch -> ch == 'Z', s)

    function checked(s, coeff)
        isapprox(imag(coeff), 0; atol=1e-9) ||
            throw(ArgumentError("term $s has a complex coefficient; the Hamiltonian is not Hermitian"))
        real(coeff)
    end

    if !parity_network                        # strictly the given order
        for (s, coeff) in terms
            abs(coeff) < 1e-15 && continue
            pauli_rotation!(c, 2 * checked(s, coeff) * dt, s; qubits=qs)
        end
        return c
    end

    diag = Pair{Int,Float64}[]                # commuting, so poolable
    rest = Any[]
    for (s, coeff) in terms
        abs(coeff) < 1e-15 && continue
        r = checked(s, coeff)
        # exp(-i·dt·c·Z_S) is exp(i·a_S·Z_S) with a_S = -dt·c
        _isdiag(s) ? push!(diag, support_mask(s) => -dt * r) : push!(rest, s => coeff)
    end
    if length(diag) == 1
        S, a = diag[1]
        phase_gadget!(c, -2a, [qs[len-b] for b in 0:len-1 if (S >> b) & 1 == 1])
    elseif !isempty(diag)
        _parity_network!(c, reverse(qs), diag) do circ, a, anchor
            push!(circ, RZ(-2a), anchor)
        end
    end
    for (s, coeff) in rest
        pauli_rotation!(c, 2 * checked(s, coeff) * dt, s; qubits=qs)
    end
    c
end

# --- phase polynomials -----------------------------------------------------

"""
    PhasePolynomial(n, coeffs)

The expansion `φ(x) = Σ_S coeffs[S+1] · (-1)^(S·x)` of a real function on
`n`-bit strings over the parity basis. `coeffs[1]` (the `S = 0` term) is a
constant, i.e. a global phase.
"""
struct PhasePolynomial
    n::Int
    coeffs::Vector{Float64}
end

"""
    phase_polynomial(φ) -> PhasePolynomial

Expand a vector of `2ⁿ` phases over the parity basis, via [`fwht!`](@ref) in
`O(n 2ⁿ)`. Inverse of [`phases`](@ref).
"""
function phase_polynomial(φ::AbstractVector{<:Real})
    N = length(φ)
    ispow2(N) || throw(ArgumentError("need a power-of-two number of phases, got $N"))
    PhasePolynomial(trailing_zeros(N), fwht(φ) ./ N)
end

"""
    phases(pp) -> Vector{Float64}

Evaluate a phase polynomial at every bitstring, recovering the phase vector.
"""
phases(pp::PhasePolynomial) = fwht(pp.coeffs)

"""
    support(pp; atol=1e-12) -> Vector{Int}

The non-constant terms with a non-negligible coefficient, as parity masks.
Sparse support is what makes a diagonal unitary cheap: `nterms` gadgets rather
than `2ⁿ`.
"""
support(pp::PhasePolynomial; atol::Real=1e-12) =
    [S for S in 1:(1 << pp.n)-1 if abs(pp.coeffs[S+1]) > atol]

"""
    nterms(pp; atol=1e-12) -> Int

Number of non-constant terms — the number of phase gadgets a naive synthesis
would emit.
"""
nterms(pp::PhasePolynomial; atol::Real=1e-12) = length(support(pp; atol=atol))

Base.show(io::IO, pp::PhasePolynomial) =
    print(io, "PhasePolynomial(", pp.n, " qubits, ", nterms(pp), " terms)")

"""
    synthesize(pp; order=:gray, atol=1e-12) -> Circuit

Realise `diag(exp(im * phases(pp)))` from its phase polynomial.

`order` picks the construction:

* `:gadgets` — one independent [`phase_gadget!`](@ref) per term, `2(|S|-1)`
  CNOTs each. The obvious baseline.
* `:gray` — a *parity network*. Pick an anchor wire, and keep the running
  parity of the current term on it; moving to the next term then costs one CNOT
  per bit of difference, so visiting the terms in Gray order costs about one
  CNOT per term instead of a whole ladder. Terms missing the anchor are handled
  by recursing on what is left, greedily choosing the wire that covers the most
  remaining terms.

Both are followed by [`cancel_adjacent_cnots!`](@ref). For a *dense* polynomial
[`diagonal`](@ref) is better still — its recursive multiplexors hit `2ⁿ - 2`.
"""
function synthesize(pp::PhasePolynomial; order::Symbol=:gray, atol::Real=1e-12)
    order in (:gray, :gadgets, :natural) ||
        throw(ArgumentError("order must be :gray or :gadgets, got $order"))
    n = pp.n
    c = Circuit(n)
    c.global_phase += pp.coeffs[1]
    terms = [S => pp.coeffs[S+1] for S in support(pp; atol=atol)]
    if order === :gray
        _parity_network!(c, _maskwires(n), terms) do circ, a, anchor
            push!(circ, RZ(-2a), anchor)
        end
    else
        for (S, a) in terms
            phase_gadget!(c, -2a, _wires(S, n))
        end
    end
    cancel_adjacent_cnots!(c)
end

# wires (ascending) carrying the bits of a parity mask
_wires(S::Integer, n::Integer) = sort!([n - b for b in 0:n-1 if (S >> b) & 1 == 1])

# wire carrying bit b of a mask, under the package's big-endian convention
_maskwires(n::Integer) = [n - b for b in 0:n-1]

"""
    _parity_network!(emit, c, wires, terms)

Walk a set of parity masks in Gray order, carrying the running parity on an
anchor wire, and call `emit(c, value, anchor)` once per term with the anchor
holding the parity of that term's mask.

Invariant inside a chain: `anchor` holds the parity of `acc`, every other wire
still holds its input value, and the anchor's own bit is in `acc` — so moving
`acc → S` costs one CNOT per bit of `acc ⊻ S`, and Gray order keeps that at one
bit for as long as possible.  Masks missing the anchor are handled by recursing
on what is left, greedily choosing the wire covering the most remaining terms.
Bit `b` of a mask lives on `wires[b+1]`.

This is the engine behind both [`synthesize`](@ref) (emit an `RZ`) and
[`multicontrolled!`](@ref) (emit a controlled-`V^{±1}`).
"""
function _parity_network!(emit, c::Circuit, wires::AbstractVector{<:Integer},
                          terms::AbstractVector{<:Pair{Int,<:Any}})
    ws = collect(Int, wires)
    nb = length(ws)
    remaining = collect(terms)
    while !isempty(remaining)
        counts = zeros(Int, nb)
        for (S, _) in remaining, b in 0:nb-1
            ((S >> b) & 1) == 1 && (counts[b+1] += 1)
        end
        b0 = argmax(counts) - 1              # anchor bit: covers the most terms
        t = ws[b0+1]
        chain = [q for q in remaining if ((q.first >> b0) & 1) == 1]
        sort!(chain; by = q -> ungray(q.first & ~(1 << b0)))
        acc = 1 << b0                        # the anchor starts holding its own value
        for (S, val) in chain
            for b in 0:nb-1
                (b == b0 || ((acc ⊻ S) >> b) & 1 == 0) && continue
                push!(c, CNOT(), ws[b+1], t)
            end
            acc = S
            emit(c, val, t)
        end
        for b in 0:nb-1                      # unwind the anchor back to its input
            (b == b0 || (acc >> b) & 1 == 0) && continue
            push!(c, CNOT(), ws[b+1], t)
        end
        remaining = [q for q in remaining if ((q.first >> b0) & 1) == 0]
    end
    c
end

# --- peephole optimisation -------------------------------------------------

"""
    cancel_adjacent_cnots!(c) -> c

Remove pairs of identical CNOTs that meet with nothing between them touching
either wire — `CNOT² = I`, so the pair is a no-op. Repeats until no pair is
left. Conservative: it never moves a gate past another that shares a wire, so
the circuit's unitary is unchanged.
"""
function cancel_adjacent_cnots!(c::Circuit)
    ops = c.ops
    changed = true
    while changed
        changed = false
        i = 1
        while i <= length(ops)
            op = ops[i]
            if op.gate.name === :CNOT
                j = i + 1
                while j <= length(ops)
                    o2 = ops[j]
                    if o2.gate.name === :CNOT && o2.qubits == op.qubits
                        deleteat!(ops, j)
                        deleteat!(ops, i)
                        changed = true
                        break
                    elseif any(q -> q in op.qubits, o2.qubits)
                        break                # blocked: something acts on our wires
                    end
                    j += 1
                end
                changed && break
            end
            i += 1
        end
    end
    c
end
