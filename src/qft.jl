# ---------------------------------------------------------------------------
# Quantum Fourier transform
#
#     QFT |x⟩ = 2^(-n/2) Σ_y exp(2πi·xy/2ⁿ) |y⟩
#
# Textbook circuit: on each wire a Hadamard, then controlled phases 2π/2^k
# from every later wire, then swaps to undo the bit reversal — n(n-1)/2
# controlled phases.  Most of them are tiny: dropping those below 2π/2^m
# (Coppersmith's approximate QFT) keeps O(n·m) gates at an error under
# n·2π/2^m, the same accuracy-for-gates trade as `truncate_terms`.
#
# The centred DFT of the three-query algorithm is this QFT, inverted, between
# two diagonal phases that are linear in the index (Nehoran–Yuen Eq. 47), so it
# is one phase gate per wire on each side: `centered_dft!`.
# ---------------------------------------------------------------------------

# Controlled phase diag(1, 1, 1, e^{iλ}); lowered, two CNOTs and three phases.
function _cphase!(c::Circuit, λ::Real, a::Int, b::Int, lower::Bool)
    if lower
        push!(c, PHASE(λ / 2), a)
        push!(c, PHASE(λ / 2), b)
        push!(c, CNOT(), a, b)
        push!(c, PHASE(-λ / 2), b)
        push!(c, CNOT(), a, b)
    else
        push!(c, controlled(PHASE(λ)), a, b)
    end
    c
end

function _swap!(c::Circuit, a::Int, b::Int, lower::Bool)
    if lower
        push!(c, CNOT(), a, b); push!(c, CNOT(), b, a); push!(c, CNOT(), a, b)
    else
        push!(c, SWAP(), a, b)
    end
    c
end

"""
    qft!(c, qubits; inverse=false, swaps=true, cutoff=nothing, lower=false) -> c

Append the quantum Fourier transform on `qubits` (`qubits[1]` most
significant): `|x⟩ ↦ 2^(-n/2) Σ_y exp(2πi·xy/2ⁿ) |y⟩`, or its inverse.

* `swaps=false` leaves the output bit-reversed and saves `⌊n/2⌋` swaps.
* `cutoff=m` is the approximate QFT: controlled phases smaller than `2π/2^m`
  are dropped, keeping `O(n·m)` of the `n(n-1)/2`.  The operator-norm error
  is at most the sum of the dropped angles, under `n·2π/2^m`, and in practice
  well below it (`0.098` against `1.18` at `n = 6, m = 5`).  `cutoff ≥ n` is
  exact.
* `lower=true` writes each controlled phase as 2 CNOTs and 3 phase gates, and
  each swap as 3 CNOTs, so that [`count_cnots`](@ref) is meaningful.
"""
function qft!(c::Circuit, qubits::AbstractVector{<:Integer}; inverse::Bool=false,
              swaps::Bool=true, cutoff::Union{Integer,Nothing}=nothing, lower::Bool=false)
    qs = collect(Int, qubits)
    isempty(qs) && throw(ArgumentError("a QFT needs at least one qubit"))
    allunique(qs) || throw(ArgumentError("repeated qubit in $qs"))
    cutoff === nothing || cutoff >= 1 || throw(ArgumentError("cutoff must be at least 1"))
    n = length(qs)
    m = cutoff === nothing ? n : cutoff

    # forward circuit as a list of steps, so the inverse is its reverse
    steps = Tuple{Symbol,Float64,Int,Int}[]
    for j in 1:n
        push!(steps, (:H, 0.0, qs[j], 0))
        for k in j+1:n
            d = k - j + 1                          # rotation 2π/2^d
            d <= m && push!(steps, (:CP, 2π / 2^d, qs[k], qs[j]))
        end
    end
    if swaps
        for j in 1:n÷2
            push!(steps, (:SWAP, 0.0, qs[j], qs[n+1-j]))
        end
    end

    for (kind, λ, a, b) in (inverse ? reverse(steps) : steps)
        if kind === :H
            push!(c, H(), a)
        elseif kind === :CP
            _cphase!(c, inverse ? -λ : λ, a, b, lower)
        else
            _swap!(c, a, b, lower)
        end
    end
    c
end

"""
    qft(n; kwargs...) -> Circuit

Standalone `n`-qubit circuit for [`qft!`](@ref), same keywords.
"""
qft(n::Integer; kwargs...) = qft!(Circuit(n), 1:n; kwargs...)

"""
    centered_dft!(c, qubits; kwargs...) -> c

Append the centred discrete Fourier transform [`centered_dft`](@ref) of the
three-query algorithm on `k = length(qubits)` wires (`K = 2^k` grid points):

    F̂ = D₁ · QFT† · D₁,    D₁ = diag(exp(2πi/K · (K-1)/2 · (u - (K-1)/4)))

(Nehoran–Yuen Eq. 47).  `D₁` is linear in `u`, so it is one phase gate per
wire plus a global phase, and the whole thing is an inverse [`qft!`](@ref)
with `2k` extra single-qubit gates.  Keywords are passed to `qft!`.
"""
function centered_dft!(c::Circuit, qubits::AbstractVector{<:Integer}; kwargs...)
    qs = collect(Int, qubits)
    k = length(qs)
    K = 2.0^k
    h = (K - 1) / 2
    d1!() = for (j, q) in enumerate(qs)
        push!(c, PHASE(2π * h / 2^j), q)          # bit j has weight 2^(k-j)
    end
    d1!()
    qft!(c, qs; inverse=true, kwargs...)
    d1!()
    c.global_phase += 2 * (-2π * h * (K - 1) / (4K))
    c
end
