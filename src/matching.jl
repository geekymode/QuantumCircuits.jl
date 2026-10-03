# ---------------------------------------------------------------------------
# Minimum-weight matching decoding
#
# In a surface code every qubit lies in at most two checks of each type, so an
# error is a set of edges in a graph whose nodes are checks, and its syndrome
# is the set of nodes with odd degree — the endpoints of error chains.  (A
# qubit in only one check gets an edge to a virtual *boundary* node.)
# Decoding is then: pair up the flipped checks, or send them to the boundary,
# so that the total chain length is smallest — the most likely error under
# independent noise — and flip the qubits along those chains.
#
# X and Z errors decode separately: Z checks find X errors, X checks find Z
# errors.  A Y error is one of each.
#
# The matching is exact, by dynamic programming over subsets of defects:
# dp[S] = min( dp[S - i] + d(i, boundary),  min_j dp[S - i - j] + d(i, j) )
# for i the lowest defect in S.  That is O(2ᵐ m) for m defects, so beyond
# `maxdefects` the decoder falls back to greedy pairing and says so.
# ---------------------------------------------------------------------------

# One matching graph: checks 1…m, boundary node m+1, and per node a BFS tree
# for shortest paths (distance and the qubit/edge leading back to the root).
struct _MatchGraph
    checks::Vector{Int}                  # stabilizer indices, in node order
    dist::Matrix{Int}                    # (m+1) × (m+1) shortest-path lengths
    parent::Array{Tuple{Int,Int},2}      # parent[root, v] = (previous node, qubit)
end

function _match_graph(code::StabilizerCode, idx::Vector{Int}, onqubits::Function)
    n = length(code)
    m = length(idx)
    B = m + 1
    # which checks contain each qubit
    owners = [Int[] for _ in 1:n]
    for (node, i) in enumerate(idx), q in 1:n
        onqubits(code.stabilizers[i], q) && push!(owners[q], node)
    end
    adj = [Tuple{Int,Int}[] for _ in 1:B]           # (neighbour, qubit)
    for q in 1:n
        o = owners[q]
        length(o) > 2 && throw(ArgumentError(
            "qubit $q lies in $(length(o)) checks of one type; matching needs at most 2 " *
            "(surface-code-like codes)"))
        length(o) == 0 && continue
        u, v = length(o) == 2 ? (o[1], o[2]) : (o[1], B)
        push!(adj[u], (v, q)); push!(adj[v], (u, q))
    end
    dist = fill(typemax(Int) ÷ 4, B, B)
    parent = fill((0, 0), B, B)
    for root in 1:B                                    # BFS from every node
        dist[root, root] = 0
        queue = [root]
        while !isempty(queue)
            u = popfirst!(queue)
            for (v, q) in adj[u]
                if dist[root, v] > dist[root, u] + 1
                    dist[root, v] = dist[root, u] + 1
                    parent[root, v] = (u, q)
                    push!(queue, v)
                end
            end
        end
    end
    _MatchGraph(idx, dist, parent)
end

# Qubits on a shortest path between nodes a and b.
function _path_qubits(g::_MatchGraph, a::Int, b::Int)
    qs = Int[]
    v = b
    while v != a
        u, q = g.parent[a, v]
        u == 0 && throw(ArgumentError("checks $a and $b are not connected"))
        push!(qs, q)
        v = u
    end
    qs
end

"""
    MatchingDecoder

A minimum-weight matching decoder for a CSS code in which every qubit lies in
at most two checks of each type — surface codes, repetition codes.  Build it
with [`matching_decoder`](@ref); use it with [`decode`](@ref) or pass it to
[`logical_error_rate`](@ref).
"""
struct MatchingDecoder
    n::Int
    gx::_MatchGraph          # Z checks → X errors
    gz::_MatchGraph          # X checks → Z errors
    maxdefects::Int
end

"""
    matching_decoder(code; maxdefects=20) -> MatchingDecoder

Build the two matching graphs of a CSS code — `Z` checks locate `X` errors,
`X` checks locate `Z` errors — with all shortest paths precomputed.
Throws for codes where a qubit lies in three or more checks of one type
(such as the Steane code), which are not graph-like.

Matching is exact for up to `maxdefects` flipped checks of a type; beyond
that it falls back to greedy pairing (see [`decode`](@ref)).
"""
function matching_decoder(code::StabilizerCode; maxdefects::Integer=20)
    is_css(code) || throw(ArgumentError("matching decoding needs a CSS code"))
    zidx = [i for (i, s) in enumerate(code.stabilizers) if any(s.z)]
    xidx = [i for (i, s) in enumerate(code.stabilizers) if any(s.x)]
    gx = _match_graph(code, zidx, (s, q) -> s.z[q])
    gz = _match_graph(code, xidx, (s, q) -> s.x[q])
    MatchingDecoder(length(code), gx, gz, maxdefects)
end

# Exact minimum-weight matching of `defects` (node numbers), each either paired
# with another or sent to the boundary.  Returns pairs (a, b), b may be boundary.
function _min_matching(g::_MatchGraph, defects::Vector{Int})
    m = length(defects)
    B = size(g.dist, 1)
    full = (1 << m) - 1
    dp = fill(typemax(Int) ÷ 4, full + 1)
    choice = zeros(Int, full + 1)                     # 0 = boundary, else partner
    dp[1] = 0
    for S in 1:full
        i = trailing_zeros(S) + 1
        rest = S & ~(1 << (i - 1))
        best = dp[rest+1] + g.dist[defects[i], B]
        ch = 0
        r = rest
        while r != 0
            j = trailing_zeros(r) + 1
            r &= r - 1
            c = dp[(rest & ~(1 << (j - 1)))+1] + g.dist[defects[i], defects[j]]
            if c < best
                best, ch = c, j
            end
        end
        dp[S+1], choice[S+1] = best, ch
    end
    pairs = Tuple{Int,Int}[]
    S = full
    while S != 0
        i = trailing_zeros(S) + 1
        j = choice[S+1]
        if j == 0
            push!(pairs, (defects[i], B))
            S &= ~(1 << (i - 1))
        else
            push!(pairs, (defects[i], defects[j]))
            S &= ~(1 << (i - 1)) & ~(1 << (j - 1))
        end
    end
    pairs, dp[full+1]
end

# Greedy fallback: repeatedly take the cheapest remaining pair or boundary hop.
function _greedy_matching(g::_MatchGraph, defects::Vector{Int})
    B = size(g.dist, 1)
    left = copy(defects)
    pairs = Tuple{Int,Int}[]
    while !isempty(left)
        best, bi, bj = typemax(Int), 0, 0
        for (a, u) in enumerate(left)
            g.dist[u, B] < best && ((best, bi, bj) = (g.dist[u, B], a, 0))
            for b in a+1:length(left)
                g.dist[u, left[b]] < best && ((best, bi, bj) = (g.dist[u, left[b]], a, b))
            end
        end
        if bj == 0
            push!(pairs, (left[bi], B)); deleteat!(left, bi)
        else
            push!(pairs, (left[bi], left[bj])); deleteat!(left, sort([bi, bj]))
        end
    end
    pairs
end

function _decode_graph(g::_MatchGraph, s::AbstractVector, maxdefects::Int)
    defects = [node for (node, i) in enumerate(g.checks) if s[i]]
    flip = Set{Int}()
    isempty(defects) && return flip, true
    exact = length(defects) <= maxdefects
    pairs = exact ? _min_matching(g, defects)[1] : _greedy_matching(g, defects)
    for (a, b) in pairs, q in _path_qubits(g, a, b)
        q in flip ? delete!(flip, q) : push!(flip, q)
    end
    flip, exact
end

"""
    decode(decoder, s) -> Union{PauliOp, Nothing}

A correction for syndrome `s`.  For a [`MatchingDecoder`](@ref): the
`X` and `Z` parts from minimum-weight matchings of the flipped `Z` and `X`
checks; its syndrome always equals `s`.  For a lookup table
([`lookup_decoder`](@ref)): the table entry, or `nothing` if the syndrome is
not in it.
"""
function decode(d::MatchingDecoder, s::AbstractVector)
    fx, _ = _decode_graph(d.gx, s, d.maxdefects)
    fz, _ = _decode_graph(d.gz, s, d.maxdefects)
    PauliOp(BitVector(q in fx for q in 1:d.n), BitVector(q in fz for q in 1:d.n))
end

decode(table::AbstractDict, s::AbstractVector) = get(table, s, nothing)
