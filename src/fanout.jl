# ---------------------------------------------------------------------------
# Fan-out and parity: the same CNOTs, arranged for depth
#
# Fan-out copies one bit onto r wires; parity collects r bits onto one.  Each
# is r CNOTs, and the obvious ladder runs them one after another in depth r.
# A doubling tree runs them in depth ~2 log₂ r.  If an r-qubit fan-out is
# allowed as a single gate, depth stops depending on r at all — the model in
# which Nehoran and Yuen show every unitary has a constant-depth circuit.
#
# Every routine takes `style`:
#   :ladder  r CNOTs in sequence             CNOTs r       depth r
#   :tree    doubling tree, exact on all     CNOTs 2r-1    depth 2⌈log₂ r⌉+1
#            inputs (not just clean targets)
#   :gate    one FANOUT / PARITY gate        CNOTs 0       depth 1
# ---------------------------------------------------------------------------

const _STYLES = (:ladder, :tree, :gate)

_check_style(style::Symbol) = style in _STYLES ||
    throw(ArgumentError("style must be one of $(_STYLES), got :$style"))

# Rounds of a doubling tree over positions 1:r rooted at 1: round h pairs
# (i, i+h) for i ≤ h.  Pairs within a round are disjoint, so each round is one
# layer, and after ⌈log₂ r⌉ rounds every position has been reached from 1.
function _tree_rounds(r::Int)
    rounds = Vector{Tuple{Int,Int}}[]
    h = 1
    while h < r
        push!(rounds, [(i, i + h) for i in 1:h if i + h <= r])
        h <<= 1
    end
    rounds
end

# Copy wire w[1] onto every wire of w (XOR) by the doubling tree, or undo it.
function _copy_tree!(c::Circuit, w::Vector{Int}; inverse::Bool=false)
    rounds = _tree_rounds(length(w))
    for round in (inverse ? reverse(rounds) : rounds), (i, j) in round
        push!(c, CNOT(), w[i], w[j])
    end
    c
end

# Collect the parity of every wire of w onto w[1], in place, or undo it.  The
# transpose of `_copy_tree!`: same pairs, CNOTs pointing the other way.
function _collect_tree!(c::Circuit, w::Vector{Int}; inverse::Bool=false)
    rounds = _tree_rounds(length(w))
    for round in (inverse ? rounds : reverse(rounds)), (i, j) in round
        push!(c, CNOT(), w[j], w[i])
    end
    c
end

function _check_wires(control::Int, ts::Vector{Int}, what::String)
    isempty(ts) && throw(ArgumentError("$what needs at least one wire"))
    allunique(ts) || throw(ArgumentError("repeated qubit in $ts"))
    control in ts && throw(ArgumentError("qubit $control appears on both sides"))
    nothing
end

"""
    fanout!(c, control, targets; style=:tree) -> c

Append a fan-out: XOR `control` into every wire of `targets`,
`|b, t₁ … t_r⟩ ↦ |b, t₁⊕b … t_r⊕b⟩`.

| `style`   | CNOTs    | depth             |
|:----------|:---------|:------------------|
| `:ladder` | `r`      | `r`               |
| `:tree`   | `2r - 1` | `2⌈log₂ r⌉ + 1`   |
| `:gate`   | `0`      | `1` (one [`FANOUT`](@ref)) |

The tree is `T · CNOT(control, t₁) · T⁻¹`, where `T` is a doubling CNOT tree
among the targets copying `t₁` to all of them.  It is exact for every input,
not only for targets prepared in `|0⟩`, at the price of `r - 1` extra CNOTs.
"""
function fanout!(c::Circuit, control::Integer, targets::AbstractVector{<:Integer};
                 style::Symbol=:tree)
    _check_style(style)
    ts = collect(Int, targets)
    _check_wires(Int(control), ts, "fan-out")
    if style === :gate
        push!(c, FANOUT(length(ts)), control, ts...)
    elseif style === :ladder
        for t in ts; push!(c, CNOT(), control, t); end
    else
        _copy_tree!(c, ts; inverse=true)
        push!(c, CNOT(), control, ts[1])
        _copy_tree!(c, ts)
    end
    c
end

"""
    parity!(c, sources, target; style=:tree) -> c

Append a parity gate: XOR the parity of `sources` into `target`,
`|x₁ … x_r, t⟩ ↦ |x₁ … x_r, t ⊕ x₁ ⊕ … ⊕ x_r⟩`.  The mirror image of
[`fanout!`](@ref), with the same costs per `style`.

(Not to be confused with [`parity`](@ref), the parity of an integer's bits.)
"""
function parity!(c::Circuit, sources::AbstractVector{<:Integer}, target::Integer;
                 style::Symbol=:tree)
    _check_style(style)
    ss = collect(Int, sources)
    _check_wires(Int(target), ss, "parity")
    if style === :gate
        push!(c, PARITY(length(ss)), ss..., target)
    elseif style === :ladder
        for s in ss; push!(c, CNOT(), s, target); end
    else
        _collect_tree!(c, ss)
        push!(c, CNOT(), ss[1], target)
        _collect_tree!(c, ss; inverse=true)
    end
    c
end
