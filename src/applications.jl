# ---------------------------------------------------------------------------
# Applications of Gray coding
#
# Five more places the one-bit-at-a-time property pays, beyond the multiplexors
# and parity networks:
#
#   multicontrolled!   C^n(U) with no ancillas, from inclusion–exclusion over
#                      parities walked in Gray order
#   gray_encoder       |x⟩ ↦ |gray(x)⟩ in n-1 CNOTs — the code itself is a
#                      linear reversible circuit
#   increment!         a counter, built from multi-controlled X
#   select!            LCU/QROM address decoding: sweeping the address in Gray
#                      order costs one X per step instead of a Hamming distance
#   gray_order         Trotter term scheduling
# ---------------------------------------------------------------------------

"""
    matrix_root(U, m) -> Matrix{ComplexF64}

Principal `2^m`-th root of a unitary: the `V` with `V^(2^m) = U` whose
eigenvalue phases are each divided by `2^m`.  Computed from a Schur
factorisation, which stays unitary even when the spectrum is degenerate.
"""
function matrix_root(U::AbstractMatrix, m::Integer)
    m >= 0 || throw(ArgumentError("m must be non-negative"))
    Uc = ComplexF64.(U)
    m == 0 && return Uc
    F = LinearAlgebra.schur(Uc)
    Z = Matrix{ComplexF64}(F.Z)
    λ = LinearAlgebra.diag(F.T)
    Z * LinearAlgebra.Diagonal(cis.(angle.(λ) ./ (1 << m))) * Z'
end

"""
    multicontrolled!(c, U, controls, target) -> c

Append `C^n(U)`: apply the one-qubit unitary `U` to `target` only when every
control is `|1⟩`.  No ancillas.

The construction is inclusion–exclusion over parities. Writing the AND of the
control bits in the parity basis,

    x₁ x₂ ⋯ xₙ  =  2^{-(n-1)} · Σ_{S ≠ ∅} (-1)^{|S|+1} · (⊕_{i ∈ S} xᵢ)

so `U^{AND(x)}` factors into one controlled-`V^{±1}` per non-empty subset `S`,
with `V = U^{1/2^{n-1}}` and the sign set by the parity of `|S|`, each
controlled on the *parity* of `S` rather than on individual wires. Those `2ⁿ-1`
parities are then walked in Gray order by the same parity network that drives
[`synthesize`](@ref) — one CNOT per step.

Costs `2ⁿ-1` controlled-`V` gates and about `2ⁿ` CNOTs: exponential in the
number of controls, which is the price of using no ancillas. Fine for the
handful of controls that show up inside [`two_level!`](@ref).
"""
function multicontrolled!(c::Circuit, U::AbstractMatrix, controls::AbstractVector, target::Integer)
    ctrl = collect(Int, controls)
    n = length(ctrl)
    size(U) == (2, 2) || throw(ArgumentError("multicontrolled! expects a one-qubit unitary"))
    target in ctrl && throw(ArgumentError("target $target is also a control"))
    Uc = ComplexF64.(U)
    n == 0 && (push!(c, Gate(:U, Uc), target); return c)
    n == 1 && (push!(c, controlled(Gate(:U, Uc)), ctrl[1], target); return c)

    V = matrix_root(Uc, n - 1)
    gV, gVd = Gate(:V, V), Gate(:Vdg, Matrix(V'))
    terms = [S => (isodd(count_ones(S)) ? 1 : -1) for S in 1:(1 << n)-1]
    _parity_network!(c, ctrl, terms) do circ, sgn, anchor
        push!(circ, controlled(sgn == 1 ? gV : gVd), anchor, target)
    end
    c
end

"""
    multicontrolled(U, controls, target; n) -> Circuit

Standalone circuit for [`multicontrolled!`](@ref).
"""
multicontrolled(U::AbstractMatrix, controls::AbstractVector, target::Integer;
                n::Integer=max(Int(target), maximum(Int, controls; init=0))) =
    multicontrolled!(Circuit(n), U, controls, target)

# --- the code as a circuit -------------------------------------------------

"""
    gray_encoder(n) -> Circuit

The map `|x⟩ ↦ |gray(x)⟩`, in `n-1` CNOTs.

`gray(x) = x ⊻ (x >> 1)` is linear over GF(2), so the Gray code *is* a CNOT
circuit: bit `b` of the output is `x_b ⊻ x_{b+1}`, one CNOT each. Sweeping from
the least significant bit upwards means every source bit is still untouched
when it is read.

See [`gray_decoder`](@ref) for the inverse and [`gray_increment`](@ref) for what
they are good for.
"""
function gray_encoder(n::Integer)
    c = Circuit(n)
    for b in 0:n-2
        push!(c, CNOT(), n - b - 1, n - b)     # bit b+1 controls bit b
    end
    c
end

"""
    gray_decoder(n) -> Circuit

Inverse of [`gray_encoder`](@ref): `|gray(x)⟩ ↦ |x⟩`, also `n-1` CNOTs.
Recovers `x_b = ⊕_{j ≥ b} g_j` by running the same ladder in reverse.
"""
function gray_decoder(n::Integer)
    c = Circuit(n)
    for b in n-2:-1:0
        push!(c, CNOT(), n - b - 1, n - b)
    end
    c
end

"""
    increment!(c, wires) -> c

Append `|x⟩ ↦ |x+1 mod 2ᵏ⟩` on `wires` (`wires[1]` most significant).

A carry chain: flip bit `b` when every lower bit is `1`. Sweeping from the top
bit down means each multi-controlled `X` reads bits that have not been touched
yet. Uses [`multicontrolled!`](@ref) for the controls.
"""
function increment!(c::Circuit, wires::AbstractVector{<:Integer})
    ws = collect(Int, wires)
    k = length(ws)
    for b in k-1:-1:0
        tgt = ws[k-b]                          # bit b
        lower = [ws[k-j] for j in 0:b-1]       # bits 0 … b-1
        multicontrolled!(c, matrix(X()), lower, tgt)
    end
    c
end

"""
    increment(k) -> Circuit

Standalone `k`-qubit incrementer, `|x⟩ ↦ |x+1 mod 2ᵏ⟩`.
"""
increment(k::Integer) = increment!(Circuit(k), 1:k)

"""
    gray_increment(n) -> Circuit

A Gray-code counter: `|gray(i)⟩ ↦ |gray(i+1)⟩`.

Decode to binary, add one, re-encode — `n-1` CNOTs each side around an
ordinary [`increment`](@ref). Successive states of the register differ in
exactly one bit, which is why Gray counters are used where a partly-updated
readout must still be a valid code word.
"""
function gray_increment(n::Integer)
    c = gray_decoder(n)
    increment!(c, 1:n)
    append!(c, gray_encoder(n))
    c
end

# --- address decoding: LCU / QROM SELECT -----------------------------------

"""
    select!(c, Us, controls, targets; order=:gray) -> c

Append a SELECT operation: apply `Us[j+1]` to `targets` when the address
register `controls` reads `|j⟩`. The building block of linear-combination-of-
unitaries and of QROM table lookup.

Each branch is a multi-controlled gate, and controls fire on `|1⟩`, so an
address pattern `j` has to be reached by `X`-ing every control where `j` has a
zero. Sweeping the addresses in **Gray order** means consecutive patterns differ
in one bit, so each step costs a single `X` instead of the Hamming distance
between neighbouring addresses:

* `order = :gray` — `2ᵏ + O(1)` X gates;
* `order = :natural` — `Σᵢ hamming(i, i+1)` of them, about `2ᵏ⁺¹`.

The multi-controlled gates themselves are emitted as single instructions; run
them through [`multicontrolled!`](@ref) to lower one-qubit targets to CNOTs.
"""
function select!(c::Circuit, Us::AbstractVector{<:AbstractMatrix},
                 controls::AbstractVector, targets::AbstractVector;
                 order::Symbol=:gray)
    order in (:gray, :natural) || throw(ArgumentError("order must be :gray or :natural, got $order"))
    ctrl = collect(Int, controls)
    tgts = collect(Int, targets)
    k = length(ctrl)
    length(Us) == 1 << k || throw(ArgumentError("expected $(1 << k) unitaries for $k address qubits, got $(length(Us))"))
    isempty(intersect(ctrl, tgts)) || throw(ArgumentError("address and target wires overlap"))
    mask = (1 << k) - 1
    addresses = order === :gray ? [gray(i) for i in 0:mask] : collect(0:mask)
    flipped = 0                                  # which controls currently carry an X
    for j in addresses
        want = ~j & mask                         # a zero in the address needs an X
        for b in 0:k-1
            ((want ⊻ flipped) >> b) & 1 == 1 && push!(c, X(), ctrl[k-b])
        end
        flipped = want
        push!(c, controlled(Gate(:U, ComplexF64.(Us[j+1])), k), ctrl..., tgts...)
    end
    for b in 0:k-1                               # restore the address register
        (flipped >> b) & 1 == 1 && push!(c, X(), ctrl[k-b])
    end
    c
end

"""
    select(Us, controls, targets; n, order=:gray) -> Circuit

Standalone circuit for [`select!`](@ref).
"""
function select(Us::AbstractVector{<:AbstractMatrix}, controls::AbstractVector,
                targets::AbstractVector;
                n::Integer=max(maximum(Int, controls; init=0), maximum(Int, targets; init=0)),
                order::Symbol=:gray)
    select!(Circuit(n), Us, controls, targets; order=order)
end

# --- Trotter term scheduling -----------------------------------------------

"""
    support_mask(s) -> Int

Bit mask of the non-identity positions of a Pauli string, with the leftmost
character as the most significant bit.
"""
support_mask(s::AbstractString) =
    sum(Int[1 << (length(s) - i) for i in 1:length(s) if s[i] != 'I']; init=0)

"""
    gray_order(terms) -> Vector

Reorder Pauli terms (as returned by [`pauli_decompose`](@ref)) by the Gray
position of their support, so that consecutive terms overlap as much as
possible.

**This does not by itself reduce the CNOT count**, and it is worth being clear
about why, because the intuition that it should is a natural one. Each term's
gadget is bracketed by its own basis-change gates, and those sit on the same
wires as the ladder ends, so [`cancel_adjacent_cnots!`](@ref) can never see two
CNOTs meet. Measured on dense random Hamiltonians the reordering is neutral to
slightly worse (3 qubits: 162 CNOTs as given, 158 reordered; the difference is
noise, not a trend).

What *does* pay is pooling terms that commute and synthesising them together —
`trotter_step!(c, terms, dt; parity_network = true)`, which carries one running
parity instead of rebuilding a ladder per term. This function is kept as an
ordering utility, and for inspecting how a Hamiltonian's supports are
distributed.
"""
gray_order(terms) = sort(collect(terms); by = q -> ungray(support_mask(first(q))))

# --- Walsh series ----------------------------------------------------------

"""
    truncate_terms(pp, k) -> PhasePolynomial

Keep only the `k` largest non-constant coefficients of a phase polynomial.

Truncating the Walsh series of a smooth phase function is the standard way to
trade accuracy for gates: a diagonal `exp(i f(x))` that would cost `2ⁿ-2` CNOTs
exactly is often good to several digits with a few dozen terms.
"""
function truncate_terms(pp::PhasePolynomial, k::Integer)
    idx = sortperm(abs.(pp.coeffs[2:end]); rev=true)
    keep = Set(idx[1:min(k, length(idx))])
    co = zeros(Float64, length(pp.coeffs))
    co[1] = pp.coeffs[1]
    for i in keep
        co[i+1] = pp.coeffs[i+1]
    end
    PhasePolynomial(pp.n, co)
end
