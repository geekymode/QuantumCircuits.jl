# ---------------------------------------------------------------------------
# Minimum-weight perfect matching — Edmonds' blossom algorithm
#
# A port of Joris van Rantwijk's maximum-weight matching (the O(n³) primal–dual
# blossom algorithm after Galil, "Efficient algorithms for finding maximum
# matching in graphs", 1986).  Internally everything is 0-based, as in the
# original; the public wrapper is 1-based.  Weights are integers and doubled,
# so the dual variables stay integral.
#
# Minimum-weight perfect matching is maximum-weight maximum-cardinality
# matching on weights W - w: every perfect matching has the same size, so
# maximising Σ(W - w) minimises Σ w.
# ---------------------------------------------------------------------------

# Maximum-weight matching; returns mate (0-based, -1 = unmatched).
function _max_weight_matching(edges::Vector{NTuple{3,Int}}, maxcardinality::Bool)
    isempty(edges) && return Int[]
    nedge = length(edges)
    nvertex = 0
    for (i, j, _) in edges
        nvertex = max(nvertex, i + 1, j + 1)
    end
    maxweight = max(0, maximum(e[3] for e in edges))
    endpoint = [p % 2 == 0 ? edges[p÷2+1][1] : edges[p÷2+1][2] for p in 0:2nedge-1]
    neighbend = [Int[] for _ in 1:nvertex]
    for k in 0:nedge-1
        i, j, _ = edges[k+1]
        push!(neighbend[i+1], 2k + 1)
        push!(neighbend[j+1], 2k)
    end
    mate = fill(-1, nvertex)
    label = zeros(Int, 2nvertex)
    labelend = fill(-1, 2nvertex)
    inblossom = collect(0:nvertex-1)
    blossomparent = fill(-1, 2nvertex)
    blossomchilds = Vector{Union{Nothing,Vector{Int}}}(nothing, 2nvertex)
    blossombase = vcat(collect(0:nvertex-1), fill(-1, nvertex))
    blossomendps = Vector{Union{Nothing,Vector{Int}}}(nothing, 2nvertex)
    bestedge = fill(-1, 2nvertex)
    blossombestedges = Vector{Union{Nothing,Vector{Int}}}(nothing, 2nvertex)
    unusedblossoms = collect(nvertex:2nvertex-1)
    dualvar = vcat(fill(maxweight, nvertex), zeros(Int, nvertex))
    allowedge = falses(nedge)
    queue = Int[]

    slack(k) = (e = edges[k+1]; dualvar[e[1]+1] + dualvar[e[2]+1] - 2e[3])

    function leaves(b)
        out = Int[]
        stack = [b]
        while !isempty(stack)
            t = pop!(stack)
            if t < nvertex
                push!(out, t)
            else
                append!(stack, blossomchilds[t+1])
            end
        end
        out
    end

    function assignlabel(w, t, p)
        b = inblossom[w+1]
        label[w+1] = label[b+1] = t
        labelend[w+1] = labelend[b+1] = p
        bestedge[w+1] = bestedge[b+1] = -1
        if t == 1
            append!(queue, leaves(b))
        elseif t == 2
            base = blossombase[b+1]
            assignlabel(endpoint[mate[base+1]+1], 1, mate[base+1] ⊻ 1)
        end
    end

    function scanblossom(v, w)
        path = Int[]
        base = -1
        while v != -1 || w != -1
            b = inblossom[v+1]
            if label[b+1] & 4 != 0
                base = blossombase[b+1]
                break
            end
            push!(path, b)
            label[b+1] = 5
            if labelend[b+1] == -1
                v = -1
            else
                v = endpoint[labelend[b+1]+1]
                b = inblossom[v+1]
                v = endpoint[labelend[b+1]+1]
            end
            if w != -1
                v, w = w, v
            end
        end
        for b in path
            label[b+1] = 1
        end
        base
    end

    function addblossom(base, k)
        v, w, _ = edges[k+1]
        bb = inblossom[base+1]
        bv = inblossom[v+1]
        bw = inblossom[w+1]
        b = pop!(unusedblossoms)
        blossombase[b+1] = base
        blossomparent[b+1] = -1
        blossomparent[bb+1] = b
        path = Int[]; endps = Int[]
        while bv != bb
            blossomparent[bv+1] = b
            push!(path, bv)
            push!(endps, labelend[bv+1])
            v = endpoint[labelend[bv+1]+1]
            bv = inblossom[v+1]
        end
        push!(path, bb)
        reverse!(path); reverse!(endps)
        push!(endps, 2k)
        while bw != bb
            blossomparent[bw+1] = b
            push!(path, bw)
            push!(endps, labelend[bw+1] ⊻ 1)
            w = endpoint[labelend[bw+1]+1]
            bw = inblossom[w+1]
        end
        blossomchilds[b+1] = path
        blossomendps[b+1] = endps
        label[b+1] = 1
        labelend[b+1] = labelend[bb+1]
        dualvar[b+1] = 0
        for v in leaves(b)
            label[inblossom[v+1]+1] == 2 && push!(queue, v)
            inblossom[v+1] = b
        end
        bestedgeto = fill(-1, 2nvertex)
        for bv in path
            nblists = blossombestedges[bv+1] === nothing ?
                      [[p ÷ 2 for p in neighbend[v+1]] for v in leaves(bv)] :
                      [blossombestedges[bv+1]]
            for nblist in nblists, k in nblist
                i, j, _ = edges[k+1]
                inblossom[j+1] == b && ((i, j) = (j, i))
                bj = inblossom[j+1]
                if bj != b && label[bj+1] == 1 &&
                   (bestedgeto[bj+1] == -1 || slack(k) < slack(bestedgeto[bj+1]))
                    bestedgeto[bj+1] = k
                end
            end
            blossombestedges[bv+1] = nothing
            bestedge[bv+1] = -1
        end
        blossombestedges[b+1] = [k for k in bestedgeto if k != -1]
        bestedge[b+1] = -1
        for k in blossombestedges[b+1]
            (bestedge[b+1] == -1 || slack(k) < slack(bestedge[b+1])) && (bestedge[b+1] = k)
        end
    end

    function expandblossom(b, endstage)
        for s in blossomchilds[b+1]
            blossomparent[s+1] = -1
            if s < nvertex
                inblossom[s+1] = s
            elseif endstage && dualvar[s+1] == 0
                expandblossom(s, endstage)
            else
                for v in leaves(s)
                    inblossom[v+1] = s
                end
            end
        end
        if !endstage && label[b+1] == 2
            childs = blossomchilds[b+1]; endps = blossomendps[b+1]
            L = length(childs)
            entrychild = inblossom[endpoint[(labelend[b+1] ⊻ 1)+1]+1]
            j = findfirst(==(entrychild), childs) - 1
            if isodd(j)
                j -= L; jstep = 1; endptrick = 0
            else
                jstep = -1; endptrick = 1
            end
            at(arr, i) = arr[mod(i, L)+1]                    # Python negative indexing
            p = labelend[b+1]
            while j != 0
                label[endpoint[(p ⊻ 1)+1]+1] = 0
                label[endpoint[(at(endps, j - endptrick) ⊻ endptrick ⊻ 1)+1]+1] = 0
                assignlabel(endpoint[(p ⊻ 1)+1], 2, p)
                allowedge[at(endps, j - endptrick)÷2+1] = true
                j += jstep
                p = at(endps, j - endptrick) ⊻ endptrick
                allowedge[p÷2+1] = true
                j += jstep
            end
            bv = at(childs, j)
            label[endpoint[(p ⊻ 1)+1]+1] = label[bv+1] = 2
            labelend[endpoint[(p ⊻ 1)+1]+1] = labelend[bv+1] = p
            bestedge[bv+1] = -1
            j += jstep
            while at(childs, j) != entrychild
                bv = at(childs, j)
                if label[bv+1] == 1
                    j += jstep
                    continue
                end
                found = -1
                for v in leaves(bv)
                    if label[v+1] != 0
                        found = v
                        break
                    end
                end
                if found != -1
                    v = found
                    label[v+1] = 0
                    label[endpoint[mate[blossombase[bv+1]+1]+1]+1] = 0
                    assignlabel(v, 2, labelend[v+1])
                end
                j += jstep
            end
        end
        label[b+1] = labelend[b+1] = -1
        blossomchilds[b+1] = blossomendps[b+1] = nothing
        blossombase[b+1] = -1
        blossombestedges[b+1] = nothing
        bestedge[b+1] = -1
        push!(unusedblossoms, b)
    end

    function augmentblossom(b, v)
        t = v
        while blossomparent[t+1] != b
            t = blossomparent[t+1]
        end
        t >= nvertex && augmentblossom(t, v)
        childs = blossomchilds[b+1]; endps = blossomendps[b+1]
        L = length(childs)
        at(arr, i) = arr[mod(i, L)+1]
        i = j = findfirst(==(t), childs) - 1
        if isodd(i)
            j -= L; jstep = 1; endptrick = 0
        else
            jstep = -1; endptrick = 1
        end
        while j != 0
            j += jstep
            t = at(childs, j)
            p = at(endps, j - endptrick) ⊻ endptrick
            t >= nvertex && augmentblossom(t, endpoint[p+1])
            j += jstep
            t = at(childs, j)
            t >= nvertex && augmentblossom(t, endpoint[(p ⊻ 1)+1])
            mate[endpoint[p+1]+1] = p ⊻ 1
            mate[endpoint[(p ⊻ 1)+1]+1] = p
        end
        blossomchilds[b+1] = vcat(childs[i+1:end], childs[1:i])
        blossomendps[b+1] = vcat(endps[i+1:end], endps[1:i])
        blossombase[b+1] = blossombase[blossomchilds[b+1][1]+1]
    end

    function augmentmatching(k)
        v, w, _ = edges[k+1]
        for (s, p) in ((v, 2k + 1), (w, 2k))
            while true
                bs = inblossom[s+1]
                bs >= nvertex && augmentblossom(bs, s)
                mate[s+1] = p
                labelend[bs+1] == -1 && break
                t = endpoint[labelend[bs+1]+1]
                bt = inblossom[t+1]
                s = endpoint[labelend[bt+1]+1]
                j = endpoint[(labelend[bt+1] ⊻ 1)+1]
                bt >= nvertex && augmentblossom(bt, j)
                mate[j+1] = labelend[bt+1]
                p = labelend[bt+1] ⊻ 1
            end
        end
    end

    for _ in 1:nvertex
        fill!(label, 0)
        fill!(bestedge, -1)
        for b in nvertex+1:2nvertex
            blossombestedges[b] = nothing
        end
        fill!(allowedge, false)
        empty!(queue)
        for v in 0:nvertex-1
            (mate[v+1] == -1 && label[inblossom[v+1]+1] == 0) && assignlabel(v, 1, -1)
        end
        augmented = false
        while true
            while !isempty(queue) && !augmented
                v = pop!(queue)
                for p in neighbend[v+1]
                    k = p ÷ 2
                    w = endpoint[p+1]
                    inblossom[v+1] == inblossom[w+1] && continue
                    kslack = 0
                    if !allowedge[k+1]
                        kslack = slack(k)
                        kslack <= 0 && (allowedge[k+1] = true)
                    end
                    if allowedge[k+1]
                        if label[inblossom[w+1]+1] == 0
                            assignlabel(w, 2, p ⊻ 1)
                        elseif label[inblossom[w+1]+1] == 1
                            base = scanblossom(v, w)
                            if base >= 0
                                addblossom(base, k)
                            else
                                augmentmatching(k)
                                augmented = true
                                break
                            end
                        elseif label[w+1] == 0
                            label[w+1] = 2
                            labelend[w+1] = p ⊻ 1
                        end
                    elseif label[inblossom[w+1]+1] == 1
                        b = inblossom[v+1]
                        (bestedge[b+1] == -1 || kslack < slack(bestedge[b+1])) && (bestedge[b+1] = k)
                    elseif label[w+1] == 0
                        (bestedge[w+1] == -1 || kslack < slack(bestedge[w+1])) && (bestedge[w+1] = k)
                    end
                end
            end
            augmented && break

            deltatype = -1
            delta = 0; deltaedge = -1; deltablossom = -1
            if !maxcardinality
                deltatype = 1
                delta = minimum(dualvar[1:nvertex])
            end
            for v in 0:nvertex-1
                if label[inblossom[v+1]+1] == 0 && bestedge[v+1] != -1
                    d = slack(bestedge[v+1])
                    if deltatype == -1 || d < delta
                        delta, deltatype, deltaedge = d, 2, bestedge[v+1]
                    end
                end
            end
            for b in 0:2nvertex-1
                if blossomparent[b+1] == -1 && label[b+1] == 1 && bestedge[b+1] != -1
                    d = slack(bestedge[b+1]) ÷ 2
                    if deltatype == -1 || d < delta
                        delta, deltatype, deltaedge = d, 3, bestedge[b+1]
                    end
                end
            end
            for b in nvertex:2nvertex-1
                if blossombase[b+1] >= 0 && blossomparent[b+1] == -1 && label[b+1] == 2 &&
                   (deltatype == -1 || dualvar[b+1] < delta)
                    delta, deltatype, deltablossom = dualvar[b+1], 4, b
                end
            end
            if deltatype == -1
                deltatype = 1
                delta = max(0, minimum(dualvar[1:nvertex]))
            end
            for v in 0:nvertex-1
                l = label[inblossom[v+1]+1]
                l == 1 ? (dualvar[v+1] -= delta) : l == 2 ? (dualvar[v+1] += delta) : nothing
            end
            for b in nvertex:2nvertex-1
                if blossombase[b+1] >= 0 && blossomparent[b+1] == -1
                    label[b+1] == 1 ? (dualvar[b+1] += delta) :
                    label[b+1] == 2 ? (dualvar[b+1] -= delta) : nothing
                end
            end
            if deltatype == 1
                break
            elseif deltatype == 2
                allowedge[deltaedge+1] = true
                i, j, _ = edges[deltaedge+1]
                label[inblossom[i+1]+1] == 0 && ((i, j) = (j, i))
                push!(queue, i)
            elseif deltatype == 3
                allowedge[deltaedge+1] = true
                i, _, _ = edges[deltaedge+1]
                push!(queue, i)
            else
                expandblossom(deltablossom, false)
            end
        end
        augmented || break
        for b in nvertex:2nvertex-1
            if blossomparent[b+1] == -1 && blossombase[b+1] >= 0 && label[b+1] == 1 && dualvar[b+1] == 0
                expandblossom(b, true)
            end
        end
    end
    [m >= 0 ? endpoint[m+1] : -1 for m in mate]
end

"""
    min_weight_perfect_matching(n, edges) -> Vector{Int}

A minimum-weight perfect matching of the graph on vertices `1:n` with
`edges = [(u, v, w), …]` (integer weights), by Edmonds' blossom algorithm in
`O(n³)`.  Returns `mate` with `mate[u] == v` for each matched pair.  Throws if
no perfect matching exists.

This is the engine behind [`matching_decoder`](@ref) once a syndrome has too
many defects for exhaustive search.
"""
function min_weight_perfect_matching(n::Integer, edges)
    iseven(n) || throw(ArgumentError("a perfect matching needs an even number of vertices"))
    n == 0 && return Int[]
    W = maximum(e[3] for e in edges; init=0) + 1
    E = NTuple{3,Int}[(Int(u) - 1, Int(v) - 1, 2 * (W - Int(w))) for (u, v, w) in edges]
    # make sure every vertex appears, even if isolated (it then has no match)
    mate = _max_weight_matching(E, true)
    length(mate) < n && append!(mate, fill(-1, n - length(mate)))
    any(<(0), mate) && throw(ArgumentError("the graph has no perfect matching"))
    mate .+ 1
end
