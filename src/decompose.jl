# ---------------------------------------------------------------------------
# Gray-code based decompositions
#
# The workhorse is the *uniformly controlled rotation* (a "multiplexed" or
# "multiplexor" gate):
#
#     UCR(α) = Σ_j |j⟩⟨j| ⊗ R(α_j)          j = 0 … 2^k-1 over k control qubits
#
# Naively you would implement each of the 2^k branches as a k-controlled
# rotation — O(k) CNOTs each, O(k 2^k) total.  The Gray-code construction does
# the whole multiplexor in exactly 2^k CNOTs and 2^k single-qubit rotations,
# with no ancillas and no multi-controlled gates at all.
#
# Why it works (see docs/graycode.md for the full derivation):
#
#   Alternate R(θ_i) on the target with CNOT(c_i → target).  With the controls
#   in basis state |j⟩ the CNOTs are just X^{bit} on the target, and since
#   X R(θ) X = R(-θ), each rotation picks up a sign equal to the parity of the
#   control bits toggled so far.  If the CNOT controls follow the Gray-code
#   flip sequence, the accumulated toggle mask after i-1 steps is exactly
#   gray(i-1), so
#
#       α_j = Σ_i (-1)^{gray(i-1) · j} θ_i
#
#   which is an orthogonal (Walsh–Hadamard) transform: θ = Mᵀα / 2^k.
#   Walking the Gray code cyclically returns the mask to 0, so all the CNOTs
#   cancel on the control register.
# ---------------------------------------------------------------------------

"""
    multiplex_angles(α) -> Vector{Float64}

Solve `α = M θ` where `M[j+1, i] = (-1)^(gray(i-1) · j)`, i.e. return the
rotation angles `θ` that a Gray-code multiplexor must apply in order to realise
the branch angles `α`.

`M` is a column-permuted Walsh–Hadamard matrix, so `M Mᵀ = 2^k I` and the
inverse is just the transpose over `2^k`:

    θ_i = 2^{-k} Σ_j (-1)^{gray(i-1) · j} α_j
"""
function multiplex_angles(α::AbstractVector{<:Real})
    m = length(α)
    ispow2(m) || throw(ArgumentError("need a power-of-two number of angles, got $m"))
    # The sum over j is a Walsh-Hadamard transform, so take it in O(m log m)
    # rather than O(m^2) and read the result off in Gray order.
    h = fwht(α)
    [h[gray(i - 1) + 1] / m for i in 1:m]
end

"""
    multiplex_matrix(k) -> Matrix{Float64}

The `2^k × 2^k` transform `M[j+1, i] = (-1)^(gray(i-1) · j)` used by
[`multiplex_angles`](@ref).  Exposed for inspection/tests; the solver never
builds it.
"""
multiplex_matrix(k::Integer) = [(-1.0)^parity(gray(i - 1), j) for j in 0:(1 << k)-1, i in 1:(1 << k)]

"""
    multiplexed_rotation!(c, kind, α, controls, target) -> c

Append a uniformly controlled `:RY` or `:RZ` rotation to `c`: when the control
register is in state `|j⟩`, `R(α[j+1])` is applied to `target`.

`controls[1]` is the **most significant** bit of `j`, matching the ket reading
order used everywhere else in this package.

Costs exactly `2^k` rotations and `2^k` CNOTs for `k = length(controls)`
(`k == 0` degenerates to a single rotation, no CNOTs).
"""
function multiplexed_rotation!(c::Circuit, kind::Symbol, α::AbstractVector{<:Real},
                               controls::AbstractVector, target::Integer)
    kind in (:RY, :RZ) || throw(ArgumentError("kind must be :RY or :RZ, got $kind"))
    ctrl = collect(Int, controls)
    k = length(ctrl)
    length(α) == 1 << k || throw(ArgumentError("expected $(1 << k) angles for $k controls, got $(length(α))"))
    target in ctrl && throw(ArgumentError("target $target is also a control"))
    allunique(ctrl) || throw(ArgumentError("repeated control qubit in $ctrl"))
    R = kind === :RY ? RY : RZ
    if k == 0
        push!(c, R(α[1]), target)
        return c
    end
    θ = multiplex_angles(α)
    for i in 1:(1 << k)
        push!(c, R(θ[i]), target)
        b = gray_flip_position(i, k)          # bit of j that flips at this step (LSB = 0)
        push!(c, CNOT(), ctrl[k-b], target)   # bit b ↔ ctrl[k-b] (controls are MSB-first)
    end
    c
end

"""
    multiplexed_ry(α, controls, target; n)

Standalone circuit for a uniformly controlled `RY` on `n` qubits (default: the
smallest register that fits the given wires).  See
[`multiplexed_rotation!`](@ref).
"""
function multiplexed_ry(α::AbstractVector{<:Real}, controls::AbstractVector, target::Integer;
                        n::Integer=max(Int(target), maximum(Int, controls; init=0)))
    multiplexed_rotation!(Circuit(n), :RY, α, controls, target)
end
"""
    multiplexed_rz(α, controls, target; n)

Standalone circuit for a uniformly controlled `RZ` on `n` qubits (default: the
smallest register that fits the given wires).  See
[`multiplexed_rotation!`](@ref).
"""
function multiplexed_rz(α::AbstractVector{<:Real}, controls::AbstractVector, target::Integer;
                        n::Integer=max(Int(target), maximum(Int, controls; init=0)))
    multiplexed_rotation!(Circuit(n), :RZ, α, controls, target)
end

"""
    diagonal!(c, φ, qubits) -> global_phase

Append a circuit implementing `diag(exp(im .* φ))` on `qubits`, returning the
leftover global phase (the caller decides what to do with it; the
[`diagonal`](@ref) wrapper folds it into the circuit).

Recursive peeling: split the index `j = (j', b)` on the last qubit.  The
*difference* of the two phases in each pair is a uniformly controlled `RZ`
(one Gray-code multiplexor), the *mean* is a diagonal on one fewer qubit.

    CNOTs(n) = CNOTs(n-1) + 2^(n-1)   ⟹   CNOTs(n) = 2^n - 2
"""
function diagonal!(c::Circuit, φ::AbstractVector{<:Real}, qubits::AbstractVector)
    n = length(qubits)
    n == 0 && return float(φ[1])
    length(φ) == 1 << n || throw(ArgumentError("expected $(1 << n) phases, got $(length(φ))"))
    half = length(φ) ÷ 2
    δ = Vector{Float64}(undef, half)   # φ(j',1) - φ(j',0)  → relative phase = RZ
    μ = Vector{Float64}(undef, half)   # mean               → diagonal on the rest
    for jp in 0:half-1
        lo, hi = φ[2jp+1], φ[2jp+2]
        δ[jp+1] = hi - lo
        μ[jp+1] = (hi + lo) / 2
    end
    gp = diagonal!(c, μ, @view qubits[1:n-1])
    multiplexed_rotation!(c, :RZ, δ, @view(qubits[1:n-1]), qubits[n])
    gp
end

"""
    diagonal(φ) -> Circuit

Circuit for the diagonal unitary `diag(exp(im .* φ))` on `log2(length(φ))`
qubits, using `2^n - 2` CNOTs.
"""
function diagonal(φ::AbstractVector{<:Real})
    N = length(φ)
    ispow2(N) || throw(ArgumentError("need a power-of-two number of phases, got $N"))
    n = trailing_zeros(N)
    c = Circuit(n)
    c.global_phase += diagonal!(c, φ, 1:n)
    c
end

"""
    prepare_state(a) -> Circuit

Möttönen-style state preparation: a circuit mapping `|0…0⟩` to the (normalised)
amplitude vector `a`.

Two Gray-code passes:

1. a cascade of uniformly controlled `RY` rotations, one level per qubit,
   distributing the *magnitudes* down the binary tree of amplitudes;
2. a [`diagonal`](@ref) unitary stamping on the *phases*.

Total CNOT count: `(2^n - 2)` from the RY cascade plus `(2^n - 2)` from the
diagonal, i.e. `2^{n+1} - 4` — versus the O(n 2^n) you would pay with
multi-controlled rotations.
"""
function prepare_state(a::AbstractVector{<:Number})
    N = length(a)
    ispow2(N) || throw(ArgumentError("need a power-of-two number of amplitudes, got $N"))
    n = trailing_zeros(N)
    nrm = LinearAlgebra.norm(a)
    nrm > 0 || throw(ArgumentError("cannot prepare the zero vector"))
    v = collect(ComplexF64, a) ./ nrm
    r = abs.(v)

    c = Circuit(n)
    for level in 1:n
        blk = 1 << (n - level + 1)          # amplitudes sharing a prefix at this level
        m = 1 << (level - 1)                # number of prefixes = branches of the multiplexor
        α = Vector{Float64}(undef, m)
        for p in 0:m-1
            off = p * blk
            n0 = LinearAlgebra.norm(@view r[off+1:off+blk÷2])
            n1 = LinearAlgebra.norm(@view r[off+blk÷2+1:off+blk])
            α[p+1] = 2 * atan(n1, n0)       # atan(0, 0) == 0: empty subtree ⟹ no rotation
        end
        multiplexed_rotation!(c, :RY, α, 1:level-1, level)
    end

    ω = angle.(v)
    c.global_phase += diagonal!(c, ω, 1:n)
    c
end
