# ---------------------------------------------------------------------------
# Detector error models
#
# Under Pauli noise in a Clifford circuit, each single fault flips a fixed set
# of detectors and possibly the logical observable — deterministically, and
# linearly over GF(2).  Listing every fault with its probability and effect
# gives the *detector error model*; merging faults with identical effect and
# weighting each by log((1-p)/p) gives the decoding graph the noise actually
# induces, diagonal edges and all.
#
# Effects are found by Pauli-frame propagation: push the fault's Pauli forward
# through the rest of the circuit (conjugating by each gate), and record which
# measurements it flips.  By linearity only X and Z at each location need
# propagating; Y and every two-qubit Pauli are XORs of those.
# ---------------------------------------------------------------------------

# The memory-experiment circuit as a flat list of steps, mirroring
# `_sample_memory` exactly.  Each step is (kind, a, b, gate):
#   :noise1  (qubit, prob-field)        :noise2 (qubit, qubit)
#   :gate    (qubits…) with gate        :meas   (ancilla, record)
#   :mflip   (record)                   :final  (data qubit, record)
struct _Step
    kind::Symbol
    qubits::Vector{Int}
    gate::Union{Gate,Nothing}
    record::Int
    prob::Symbol            # which NoiseModel field governs a noise step
end

function _memory_steps(code::StabilizerCode, rounds::Int)
    n, r = length(code), length(code.stabilizers)
    sc = syndrome_circuit(code)
    steps = _Step[]
    rec(τ, i) = (τ - 1) * r + i                  # ancilla records first
    for τ in 1:rounds
        for q in 1:n
            push!(steps, _Step(:noise1, [q], nothing, 0, :pdata))
        end
        for op in sc.ops
            push!(steps, _Step(:gate, op.qubits, op.gate, 0, :none))
            if length(op.qubits) == 1
                push!(steps, _Step(:noise1, op.qubits, nothing, 0, :p1))
            else
                push!(steps, _Step(:noise2, op.qubits, nothing, 0, :p2))
            end
        end
        for i in 1:r
            push!(steps, _Step(:meas, [n + i], nothing, rec(τ, i), :none))
            push!(steps, _Step(:mflip, Int[], nothing, rec(τ, i), :pmeas))
        end
    end
    for q in 1:n
        push!(steps, _Step(:final, [q], nothing, rounds * r + q, :none))
        push!(steps, _Step(:mflip, Int[], nothing, rounds * r + q, :pmeas))
    end
    steps
end

# Conjugate a Pauli frame (x, z bits) by a Clifford gate.
function _frame_gate!(x::BitVector, z::BitVector, g::Gate, qs::Vector{Int})
    nm = g.name
    if nm === :H
        a = qs[1]; x[a], z[a] = z[a], x[a]
    elseif nm === :S || nm === :Sdg
        a = qs[1]; z[a] ⊻= x[a]
    elseif nm === :CNOT || nm === :CX
        a, b = qs; x[b] ⊻= x[a]; z[a] ⊻= z[b]
    elseif nm === :CZ
        a, b = qs; z[a] ⊻= x[b]; z[b] ⊻= x[a]
    elseif nm === :CY                                  # CY = S_b · CX · S_b†
        a, b = qs
        z[b] ⊻= x[b]; x[b] ⊻= x[a]; z[a] ⊻= z[b]; z[b] ⊻= x[b]
    elseif nm === :SWAP
        a, b = qs; x[a], x[b] = x[b], x[a]; z[a], z[b] = z[b], z[a]
    elseif nm in (:I, :X, :Y, :Z)
    else
        throw(ArgumentError("$(label(g)) is not supported in a detector error model"))
    end
end

# Records flipped by a Pauli (x, z) inserted just after step `s`.
function _propagate(steps::Vector{_Step}, s::Int, x::BitVector, z::BitVector, nrec::Int)
    flipped = falses(nrec)
    for k in s+1:length(steps)
        st = steps[k]
        if st.kind === :gate
            _frame_gate!(x, z, st.gate, st.qubits)
        elseif st.kind === :meas
            a = st.qubits[1]
            flipped[st.record] = x[a]
            x[a] = false; z[a] = false                 # reset clears the frame
        elseif st.kind === :final
            flipped[st.record] = x[st.qubits[1]]
        end
    end
    flipped
end

"""
    DetectorErrorModel

The independent fault mechanisms of a noisy memory experiment: for each, the
detectors it flips, whether it flips the logical `Z̄`, and its probability.
Built by [`detector_error_model`](@ref); decoded by
[`memory_experiment`](@ref) with `decoder = :dem`.
"""
struct DetectorErrorModel
    ndetectors::Int
    mechanisms::Vector{Tuple{Vector{Int},Bool,Float64}}
    decomposed::Int          # mechanisms with > 2 detectors split into edges
    dropped::Int             # … that could not be split
end

function Base.show(io::IO, dem::DetectorErrorModel)
    print(io, "DetectorErrorModel(", dem.ndetectors, " detectors, ",
          length(dem.mechanisms), " mechanisms)")
end

"""
    detector_error_model(code, noise; rounds) -> DetectorErrorModel

Every single fault of the memory experiment under `noise`, traced by
Pauli-frame propagation to the detectors it flips and its effect on `Z̄`;
faults with the same effect are merged (`p = p₁ + p₂ - 2p₁p₂`).  Mechanisms
flipping more than two detectors — a gate fault spreading to several data
qubits — are split into existing edges where the split is consistent.
"""
function detector_error_model(code::StabilizerCode, noise::NoiseModel;
                              rounds::Integer=something(code.distance, 3))
    is_css(code) || throw(ArgumentError("detector_error_model needs a CSS code"))
    n, r = length(code), length(code.stabilizers)
    stabs = code.stabilizers
    zch = [i for (i, s) in enumerate(stabs) if any(s.z)]
    m, L = length(zch), rounds + 1
    steps = _memory_steps(code, rounds)
    nrec = rounds * r + n
    N = n + r
    rec(τ, i) = (τ - 1) * r + i
    zl = code.logical_z[1].z

    # records → (detectors, logical)
    function effect(flipped::BitVector)
        dets = Int[]
        for (c, i) in enumerate(zch), τ in 1:L
            v = if τ <= rounds
                flipped[rec(τ, i)] ⊻ (τ > 1 && flipped[rec(τ - 1, i)])
            else
                isodd(count(q -> stabs[i].z[q] && flipped[rounds*r+q], 1:n)) ⊻ flipped[rec(rounds, i)]
            end
            v && push!(dets, (τ - 1) * m + c)
        end
        dets, isodd(count(q -> zl[q] && flipped[rounds*r+q], 1:n))
    end

    acc = Dict{Tuple{Vector{Int},Bool},Float64}()
    add!(key, p) = p > 0 && (acc[key] = (q = get(acc, key, 0.0); q + p - 2q * p))
    unit(q, xb, zb) = (x = falses(N); z = falses(N); x[q] = xb; z[q] = zb; (x, z))
    for (s, st) in enumerate(steps)
        if st.kind === :noise1 || st.kind === :noise2
            p = getfield(noise, st.prob)
            p > 0 || continue
            qs = st.qubits
            # X and Z basis effects per qubit, combined by XOR
            base = Dict{Tuple{Int,Int},BitVector}()
            for (j, q) in enumerate(qs), kd in 1:2
                x, z = unit(q, kd == 1, kd == 2)
                base[(j, kd)] = _propagate(steps, s, x, z, nrec)
            end
            paulis = st.kind === :noise1 ? [(k,) for k in 1:3] :
                     [(ka, kb) for ka in 0:3, kb in 0:3 if (ka, kb) != (0, 0)]
            for P in paulis
                f = falses(nrec)
                for (j, k) in enumerate(P)
                    k == 0 && continue
                    k != 3 && (f .⊻= base[(j, 1)])       # X part (X or Y)
                    k != 1 && (f .⊻= base[(j, 2)])       # Z part (Z or Y)
                end
                key = effect(f)
                isempty(key[1]) && !key[2] && continue
                add!(key, p / length(paulis))
            end
        elseif st.kind === :mflip
            noise.pmeas > 0 || continue
            f = falses(nrec); f[st.record] = true
            key = effect(f)
            isempty(key[1]) && !key[2] && continue        # e.g. an X-check record
            add!(key, noise.pmeas)
        end
    end

    # split hyperedges into existing ≤2-detector mechanisms, logical-consistently
    simple = Dict(k => v for (k, v) in acc if length(k[1]) <= 2)
    edgekeys = Dict{Vector{Int},Vector{Bool}}()
    for (k, _) in simple
        push!(get!(edgekeys, k[1], Bool[]), k[2])
    end
    mech = [(k[1], k[2], v) for (k, v) in simple]
    decomposed = dropped = 0
    for (k, v) in acc
        length(k[1]) <= 2 && continue
        split = _split_hyperedge(k[1], k[2], edgekeys)
        if split === nothing
            dropped += 1
        else
            decomposed += 1
            for (part, lg) in split
                push!(mech, (part, lg, v))
            end
        end
    end
    DetectorErrorModel(m * L, mech, decomposed, dropped)
end

# Write a detector set as a disjoint union of existing edges (pairs or single
# boundary hops) whose logical flips XOR to `logical`; nothing if impossible.
function _split_hyperedge(dets::Vector{Int}, logical::Bool, edges::Dict{Vector{Int},Vector{Bool}})
    isempty(dets) && return logical ? nothing : Tuple{Vector{Int},Bool}[]
    d = dets[1]
    rest = dets[2:end]
    for k in 0:length(rest)
        part = k == 0 ? [d] : sort([d, rest[k]])
        haskey(edges, part) || continue
        remaining = k == 0 ? rest : deleteat!(copy(rest), k)
        for lg in unique(edges[part])
            sub = _split_hyperedge(remaining, logical ⊻ lg, edges)
            sub !== nothing && return vcat([(part, lg)], sub)
        end
    end
    nothing
end

# The weighted decoding graph: detectors 1…D and boundary D+1; parallel edges
# keep the most likely mechanism.  Weights log((1-p)/p), scaled to integers.
struct _DEMGraph
    D::Int
    adj::Vector{Vector{Tuple{Int,Int,Bool}}}     # (neighbour, weight, logical)
end

function _dem_graph(dem::DetectorErrorModel)
    D = dem.ndetectors
    B = D + 1
    best = Dict{Tuple{Int,Int},Tuple{Float64,Bool}}()
    for (dets, lg, p) in dem.mechanisms
        isempty(dets) && continue                  # undetectable: nothing to match
        u, v = length(dets) == 1 ? (dets[1], B) : (dets[1], dets[2])
        key = (min(u, v), max(u, v))
        (!haskey(best, key) || best[key][1] < p) && (best[key] = (p, lg))
    end
    adj = [Tuple{Int,Int,Bool}[] for _ in 1:B]
    for ((u, v), (p, lg)) in best
        w = max(1, round(Int, 1000 * log((1 - p) / p)))
        push!(adj[u], (v, w, lg)); push!(adj[v], (u, w, lg))
    end
    _DEMGraph(D, adj)
end

# Dijkstra from `src`: distances and logical parity of each shortest path.
function _dijkstra_parity(g::_DEMGraph, src::Int)
    N = length(g.adj)
    dist = fill(typemax(Int) ÷ 4, N)
    par = falses(N)
    done = falses(N)
    dist[src] = 0
    pq = [(0, src)]
    while !isempty(pq)
        i = argmin(first.(pq))
        d, u = pq[i]
        deleteat!(pq, i)
        done[u] && continue
        done[u] = true
        for (v, w, lg) in g.adj[u]
            if d + w < dist[v]
                dist[v] = d + w
                par[v] = par[u] ⊻ lg
                push!(pq, (dist[v], v))
            end
        end
    end
    dist, par
end

function _decode_dem(g::_DEMGraph, det::BitMatrix, m::Int)
    defects = [(τ - 1) * m + c for τ in 1:size(det, 2) for c in 1:m if det[c, τ]]
    k = length(defects)
    k == 0 && return false
    B = g.D + 1
    sp = [_dijkstra_parity(g, d) for d in defects]
    edges = Tuple{Int,Int,Int}[]
    for a in 1:k
        push!(edges, (a, k + a, sp[a][1][B]))
        for b in a+1:k
            push!(edges, (a, b, sp[a][1][defects[b]]))
            push!(edges, (k + a, k + b, 0))
        end
    end
    mate = min_weight_perfect_matching(2k, edges)
    flip = false
    for a in 1:k
        b = mate[a]
        if b > k
            flip ⊻= sp[a][2][B]
        elseif a < b
            flip ⊻= sp[a][2][defects[b]]
        end
    end
    flip
end
