# ---------------------------------------------------------------------------
# Circuit-level noise and the memory experiment
#
# Real syndrome extraction is itself noisy: gates fail, measurements lie.  So
# syndromes are measured over many rounds and the decoder looks at
# *detectors* — the change in a check's value from one round to the next.  A
# data error flips one check's history from some round on (one detector per
# check it touches, in that round); a measurement error flips a single
# outcome (two detectors, the same check in consecutive rounds).  Decoding is
# matching on a space-time graph: space edges for data errors, time edges for
# measurement errors.
#
# The memory experiment: prepare |0̄⟩, run `rounds` noisy rounds, measure every
# data qubit, and ask whether the decoded logical value is still 0.
# ---------------------------------------------------------------------------

"""
    NoiseModel(; p1=0, p2=0, pmeas=0, pdata=0)

Pauli noise for [`memory_experiment`](@ref):

* `p1` — after each one-qubit gate, a random `X`, `Y` or `Z` with probability `p1`;
* `p2` — after each two-qubit gate, one of the 15 non-identity two-qubit
  Paulis with probability `p2`;
* `pmeas` — each measurement result is flipped with probability `pmeas`;
* `pdata` — at the start of every round, each data qubit is depolarized with
  probability `pdata`.

[`circuit_noise`](@ref)`(p)` and [`phenomenological_noise`](@ref)`(p)` are
the two standard settings.
"""
Base.@kwdef struct NoiseModel
    p1::Float64 = 0.0
    p2::Float64 = 0.0
    pmeas::Float64 = 0.0
    pdata::Float64 = 0.0
end

"""
    circuit_noise(p) -> NoiseModel

Every gate, measurement and idle round fails with probability `p` — the usual
"circuit-level" benchmark.
"""
circuit_noise(p::Real) = NoiseModel(p1=p, p2=p, pmeas=p, pdata=p)

"""
    phenomenological_noise(p) -> NoiseModel

Data qubits depolarize and measurements flip, each with probability `p`, but
the gates themselves are perfect — the textbook model in which matching on the
space-time graph is exactly the right decoder.
"""
phenomenological_noise(p::Real) = NoiseModel(pmeas=p, pdata=p)

function _depolarize!(t::Tableau, q::Int, p::Float64, rng)
    p > 0 && rand(rng) < p || return
    k = rand(rng, 1:3)
    _pauli1!(t, q, k != 3, k != 1)
end

function _depolarize2!(t::Tableau, a::Int, b::Int, p::Float64, rng)
    p > 0 && rand(rng) < p || return
    k = rand(rng, 1:15)                       # 1…15: pairs (ka, kb) ≠ (0, 0)
    ka, kb = divrem(k, 4)
    ka != 0 && _pauli1!(t, a, ka != 3, ka != 1)
    kb != 0 && _pauli1!(t, b, kb != 3, kb != 1)
end

# One shot.  Returns the Z-check detector matrix (checks × rounds+1) and the
# measured logical Z̄ (true = flipped).  For tests: `inject` = (round, PauliOp
# on data) adds an error at the start of that round, and `flip` = (round, i)
# flips the reported outcome of stabilizer i in that round.
function _sample_memory(code::StabilizerCode, noise::NoiseModel, rounds::Int, rng;
                        inject=nothing, flip=nothing)
    n = length(code)
    stabs = code.stabilizers
    r = length(stabs)
    zchecks = [i for (i, s) in enumerate(stabs) if any(s.z)]
    t = Tableau(n + r)                         # data |0…0⟩: every Z check is +1
    sc = syndrome_circuit(code)
    prev = falses(r)
    det = falses(length(zchecks), rounds + 1)
    for round in 1:rounds
        for q in 1:n
            _depolarize!(t, q, noise.pdata, rng)
        end
        if inject !== nothing && inject[1] == round
            pauli!(t, inject[2], 1:n)
        end
        for op in sc.ops
            apply!(t, op.gate, op.qubits...)
            if length(op.qubits) == 1
                _depolarize!(t, op.qubits[1], noise.p1, rng)
            else
                _depolarize2!(t, op.qubits[1], op.qubits[2], noise.p2, rng)
            end
        end
        outcome = falses(r)
        for i in 1:r
            a = n + i
            actual = measure!(t, a; rng=rng)
            actual && _pauli1!(t, a, true, false)          # reset to |0⟩
            outcome[i] = actual ⊻ (noise.pmeas > 0 && rand(rng) < noise.pmeas) ⊻ (flip == (round, i))
        end
        for (row, i) in enumerate(zchecks)
            det[row, round] = outcome[i] ⊻ prev[i]
        end
        prev = outcome
    end
    # final transversal Z measurement of the data
    bits = BitVector([measure!(t, q; rng=rng) ⊻ (noise.pmeas > 0 && rand(rng) < noise.pmeas) for q in 1:n])
    for (row, i) in enumerate(zchecks)
        det[row, rounds+1] = isodd(count(stabs[i].z .& bits)) ⊻ prev[i]
    end
    logical = isodd(count(code.logical_z[1].z .& bits))
    det, logical
end

# Space-time matching graph for the Z checks: node (c, τ) for check c and
# layer τ = 1…L, numbered (τ-1)·m + c, plus one boundary node.  Space edges
# carry the qubit (to track the logical), time edges carry 0.
struct _SpaceTimeGraph
    m::Int
    L::Int
    adj::Vector{Vector{Tuple{Int,Int}}}          # (neighbour, qubit or 0)
    onlogical::BitVector                          # qubit lies on the support of Z̄
end

function _spacetime_graph(code::StabilizerCode, rounds::Int)
    n = length(code)
    zchecks = [i for (i, s) in enumerate(code.stabilizers) if any(s.z)]
    m, L = length(zchecks), rounds + 1
    owners = [Int[] for _ in 1:n]
    for (c, i) in enumerate(zchecks), q in 1:n
        code.stabilizers[i].z[q] && push!(owners[q], c)
    end
    any(o -> length(o) > 2, owners) &&
        throw(ArgumentError("memory decoding needs every qubit in at most two Z checks"))
    B = m * L + 1
    adj = [Tuple{Int,Int}[] for _ in 1:B]
    node(c, τ) = (τ - 1) * m + c
    for τ in 1:L
        for q in 1:n
            o = owners[q]
            isempty(o) && continue
            u = node(o[1], τ)
            v = length(o) == 2 ? node(o[2], τ) : B
            push!(adj[u], (v, q)); push!(adj[v], (u, q))
        end
        if τ < L
            for c in 1:m
                push!(adj[node(c, τ)], (node(c, τ + 1), 0))
                push!(adj[node(c, τ + 1)], (node(c, τ), 0))
            end
        end
    end
    _SpaceTimeGraph(m, L, adj, BitVector(code.logical_z[1].z))
end

# BFS from `src`: distance to every node, and whether the shortest path to it
# crosses the logical support an odd number of times.
function _bfs_parity(g::_SpaceTimeGraph, src::Int)
    N = length(g.adj)
    dist = fill(-1, N)
    par = falses(N)
    dist[src] = 0
    queue = [src]
    head = 1
    while head <= length(queue)
        u = queue[head]; head += 1
        for (v, q) in g.adj[u]
            if dist[v] < 0
                dist[v] = dist[u] + 1
                par[v] = par[u] ⊻ (q > 0 && g.onlogical[q])
                push!(queue, v)
            end
        end
    end
    dist, par
end

# Decode one shot's detectors: does the minimum-weight explanation flip Z̄?
function _decode_memory(g::_SpaceTimeGraph, det::BitMatrix)
    defects = [(τ - 1) * g.m + c for τ in 1:g.L for c in 1:g.m if det[c, τ]]
    k = length(defects)
    k == 0 && return false
    B = length(g.adj)
    bfs = [_bfs_parity(g, d) for d in defects]
    edges = Tuple{Int,Int,Int}[]
    for a in 1:k
        push!(edges, (a, k + a, bfs[a][1][B]))
        for b in a+1:k
            push!(edges, (a, b, bfs[a][1][defects[b]]))
            push!(edges, (k + a, k + b, 0))
        end
    end
    mate = min_weight_perfect_matching(2k, edges)
    flip = false
    for a in 1:k
        b = mate[a]
        if b > k
            flip ⊻= bfs[a][2][B]
        elseif a < b
            flip ⊻= bfs[a][2][defects[b]]
        end
    end
    flip
end

"""
    memory_experiment(code, noise; rounds=code distance, shots=1000, rng) -> (rate, failures)

Estimate the logical error rate of storing `|0̄⟩` under circuit-level `noise`
([`NoiseModel`](@ref)), on the tableau simulator: start in `|0…0⟩`, run
`rounds` rounds of the [`syndrome_circuit`](@ref) — every gate followed by its
noise, every ancilla measured (with flips) and reset — then measure every data
qubit.  Detectors are the round-to-round changes of the `Z` checks, decoded by
minimum-weight matching on the space-time graph; a shot fails when the decoded
`Z̄` disagrees with `0`.

For CSS codes whose qubits lie in at most two `Z` checks (surface codes).

Under [`phenomenological_noise`](@ref) this decoder is the standard one and
the surface code's crossover sits near the known threshold.  Under
[`circuit_noise`](@ref) it is deliberately simple, and the crossover (about
0.3% here) sits below the 0.5–1% that tuned decoders reach, for three reasons:
the graph lacks the diagonal edges that a gate fault mid-round creates, all
edges weigh the same although data and measurement faults have different
probabilities, and the checks are measured one after another rather than in
parallel.  A detector error model — every single fault traced to its detector
signature and probability — would remove the first two.
"""
function memory_experiment(code::StabilizerCode, noise::NoiseModel;
                           rounds::Integer=something(code.distance, 3),
                           shots::Integer=1000, rng=Random.default_rng())
    is_css(code) || throw(ArgumentError("memory_experiment needs a CSS code"))
    dimension(code) == 1 || throw(ArgumentError("memory_experiment needs one logical qubit"))
    rounds >= 1 || throw(ArgumentError("need at least one round"))
    g = _spacetime_graph(code, rounds)
    failures = 0
    for _ in 1:shots
        det, logical = _sample_memory(code, noise, rounds, rng)
        (_decode_memory(g, det) != logical) && (failures += 1)
    end
    failures / shots, failures
end
