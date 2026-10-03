# ---------------------------------------------------------------------------
# T-count optimisation as Reed–Muller decoding (after Amy and Mosca)
#
# A circuit of CNOTs and phase gates T = diag(1, ω) (ω = e^{iπ/4}), S = T²,
# Z = T⁴ acts as
#
#     |x⟩ ↦ ω^{f(x)} |A x⟩,     f(x) = Σ_{y ≠ 0} κ_y · (y·x)  mod 8,
#
# where (y·x) is the parity of the input bits selected by y and A is the
# linear map of the CNOTs.  Each *odd* coefficient κ_y needs one T gate; even
# ones are S and Z, which are Clifford.  So the T-count of the circuit is at
# least the number of odd κ_y, and merging the gates by parity reaches it.
#
# The Reed–Muller step: for a monomial M of degree k ≤ n - 4, adding 1 to κ_y
# for every y ≠ 0 containing M changes f by
#
#     Σ_{y ⊇ M} (y·x) = 2^{n-k-1}·[x outside M ≠ 0] + 2^{n-k}·(M·x)·[x outside M = 0]
#
# which is 0 mod 8.  The unitary is unchanged — global phase included — while
# the odd pattern moves by the evaluation table of M, a generator of the
# punctured Reed–Muller code RM(n-4, n)*.  The minimum T-count is the distance
# from the odd pattern to that code: decoding.
# ---------------------------------------------------------------------------

# Phase gates as powers of T (mod 8), with the global phase they carry.
function _t_power(g::Gate)
    nm = g.name
    nm === :T && return 1, 0.0
    nm === :Tdg && return 7, 0.0
    nm === :S && return 2, 0.0
    nm === :Sdg && return 6, 0.0
    nm === :Z && return 4, 0.0
    nm === :I && return 0, 0.0
    if nm === :P || nm === :RZ
        θ = g.params[1]
        k = θ / (π / 4)
        isapprox(k, round(k); atol=1e-9) ||
            throw(ArgumentError("$(label(g)) is not a multiple of π/4"))
        return mod(round(Int, k), 8), nm === :RZ ? -θ / 2 : 0.0   # RZ(θ) = e^{-iθ/2} P(θ)
    end
    nothing
end

"""
    t_count(c) -> Int

Number of `T`-type gates in a circuit: `T`, `T†`, and `PHASE` or `RZ` by an
odd multiple of `π/4`.
"""
function t_count(c::Circuit)
    count(c.ops) do op
        tp = _t_power(op.gate)
        tp !== nothing && isodd(tp[1])
    end
end

"""
    z8_phase_polynomial(c) -> (κ, masks, phase)

The sum-over-paths form of a circuit of `CNOT`, `SWAP` and phase gates
(`T`, `T†`, `S`, `S†`, `Z`, and `PHASE`/`RZ` by multiples of `π/4`):

    |x⟩ ↦ e^{i·phase} ω^{f(x)} |x′⟩,    f(x) = Σ_y κ[y+1]·(y·x) mod 8,

with `ω = e^{iπ/4}` and output wire `q` holding the parity `masks[q]` of the
input.  Parity masks use the package's convention: bit `n - q` is wire `q`.
"""
function z8_phase_polynomial(c::Circuit)
    n = c.nqubits
    n <= 20 || throw(ArgumentError("phase polynomials are tabulated for n ≤ 20"))
    masks = [1 << (n - q) for q in 1:n]
    κ = zeros(Int, 1 << n)
    phase = c.global_phase
    for op in c.ops
        g, qs = op.gate, op.qubits
        if g.name === :CNOT || g.name === :CX
            masks[qs[2]] ⊻= masks[qs[1]]
        elseif g.name === :SWAP
            masks[qs[1]], masks[qs[2]] = masks[qs[2]], masks[qs[1]]
        else
            tp = _t_power(g)
            tp === nothing && throw(ArgumentError(
                "$(label(g)) is not a CNOT, SWAP or π/4 phase gate; T-count optimisation " *
                "works on {CNOT, T} circuits"))
            k, φ = tp
            y = masks[qs[1]]
            κ[y+1] = mod(κ[y+1] + k, 8)
            phase += φ
        end
    end
    κ, masks, phase
end

# Phase gate T^k on one wire, as at most one T-type gate and one Clifford.
function _emit_t_power!(c::Circuit, k::Int, q::Int)
    k = mod(k, 8)
    k == 7 && return push!(c, Tdg(), q)
    isodd(k) && push!(c, T(), q)
    rest = isodd(k) ? k - 1 : k
    rest == 2 && push!(c, S(), q)
    rest == 4 && push!(c, Z(), q)
    rest == 6 && push!(c, Sdg(), q)
    c
end

# Rebuild a circuit from its phase polynomial: a parity network for the
# diagonal part, then a CNOT network for the linear part.
function _synthesize_z8(n::Int, κ::Vector{Int}, masks::Vector{Int}, phase::Float64)
    c = Circuit(n)
    c.global_phase = phase
    terms = [y => κ[y+1] for y in 1:(1 << n)-1 if mod(κ[y+1], 8) != 0]
    _parity_network!(c, _maskwires(n), terms) do circ, k, anchor
        _emit_t_power!(circ, k, anchor)
    end
    # wire q ends with parity masks[q] of the inputs: |x⟩ ↦ |x M⟩ with
    # M[j, q] = bit of masks[q] belonging to wire j
    M = BitMatrix([(masks[q] >> (n - j)) & 1 == 1 for j in 1:n, q in 1:n])
    M == BitMatrix(LinearAlgebra.I(n)) || _linear_map!(c, M, collect(1:n))
    cancel_adjacent_cnots!(c)
end

# Generators of RM(n-4, n)* on the nonzero points: one per monomial of degree
# ≤ n - 4, as (mask, support word) with bit y-1 set when y ⊇ mask.
function _rm_generators(n::Int)
    n >= 4 || return Tuple{Int,UInt128}[]
    [(M, sum((UInt128(1) << (y - 1) for y in 1:(1 << n)-1 if y & M == M); init=UInt128(0)))
     for M in 0:(1 << n)-1 if count_ones(M) <= n - 4]
end

"""
    optimize_t_count(c; exact=nothing) -> Circuit

An equivalent circuit (global phase included) with fewer `T` gates, for
circuits of `CNOT`, `SWAP` and `π/4` phase gates (see
[`z8_phase_polynomial`](@ref)).

1. Merge every phase gate into its parity's coefficient: the T-count drops to
   the number of odd coefficients.
2. Move the odd pattern by the closest codeword of the punctured Reed–Muller
   code `RM(n-4, n)*` — each generator changes nothing about the unitary —
   to reach the fewest odd coefficients.  For `n ≤ 6` the search is exact, a
   Gray-code walk over all `2ᵏ` codewords (`k = 22` at `n = 6`).  At `n = 7`
   (`k = 64`) it is a local search over one- and two-generator moves with
   restarts — on `n = 6` patterns that lands within 4 of the optimum, mean
   gap about 2 — so it is good but not guaranteed minimal.
3. Resynthesise with a parity network and a CNOT network.

For `n ≤ 3` the code is empty and step 2 does nothing.
"""
function optimize_t_count(c::Circuit; exact::Union{Bool,Nothing}=nothing)
    n = c.nqubits
    κ, masks, phase = z8_phase_polynomial(c)
    gens = _rm_generators(n)
    if !isempty(gens)
        n <= 7 || throw(ArgumentError("optimize_t_count supports n ≤ 7"))
        odd = sum((UInt128(1) << (y - 1) for y in 1:(1 << n)-1 if isodd(κ[y+1])); init=UInt128(0))
        chosen = (exact === nothing ? n <= 6 : exact) ? _rm_nearest_exact(odd, gens) :
                                                        _rm_nearest_search(odd, gens)
        for i in chosen                               # subtract each generator: f unchanged
            M = gens[i][1]
            for y in 1:(1 << n)-1
                y & M == M && (κ[y+1] = mod(κ[y+1] - 1, 8))
            end
        end
    end
    _synthesize_z8(n, κ, masks, phase)
end

# Exact: walk all combinations of generators in Gray order (one XOR per step).
function _rm_nearest_exact(odd::UInt128, gens)
    k = length(gens)
    k <= 26 || throw(ArgumentError("exact search over 2^$k codewords is too large"))
    w, best, bestt = odd, count_ones(odd), 0
    for t in 1:(1 << k)-1
        w ⊻= gens[gray_flip_position(t)+1][2]
        cw = count_ones(w)
        cw < best && ((best, bestt) = (cw, t))
    end
    g = gray(bestt)
    [i for i in 1:k if (g >> (i - 1)) & 1 == 1]
end

# Heuristic for when 2ᵏ is too many: local search that toggles one or two
# generators at a time, from the empty set and from random starts, keeping
# the best.  On random n = 6 patterns it lands within 4 of the exact optimum
# (mean gap ≈ 2); plain single-toggle descent was ≈ 12 off.
function _rm_nearest_search(odd::UInt128, gens; restarts::Int=30, rng=Random.Xoshiro(0))
    k = length(gens)
    words = [g[2] for g in gens]
    bestw, bestset = count_ones(odd), falses(k)
    for r in 1:restarts
        set = r == 1 ? falses(k) : BitVector(rand(rng, Bool, k))
        w = odd
        for i in findall(set)
            w ⊻= words[i]
        end
        while true
            cur, bi, bj = count_ones(w), 0, 0
            for i in 1:k
                wi = w ⊻ words[i]
                ci = count_ones(wi)
                ci < cur && ((cur, bi, bj) = (ci, i, 0))
                for j in i+1:k
                    cj = count_ones(wi ⊻ words[j])
                    cj < cur && ((cur, bi, bj) = (cj, i, j))
                end
            end
            bi == 0 && break
            w ⊻= words[bi]; set[bi] = !set[bi]
            bj != 0 && (w ⊻= words[bj]; set[bj] = !set[bj])
        end
        count_ones(w) < bestw && ((bestw, bestset) = (count_ones(w), copy(set)))
    end
    findall(bestset)
end
