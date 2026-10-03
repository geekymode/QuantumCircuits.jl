@testset "minimum-weight perfect matching" begin
    rng = Random.Xoshiro(17)
    function brute(n, W)
        best = typemax(Int)
        function rec(left, acc)
            isempty(left) && (best = min(best, acc); return)
            for k in 2:length(left)
                W[left[1], left[k]] == typemax(Int) && continue
                rec(deleteat!(copy(left), [1, k]), acc + W[left[1], left[k]])
            end
        end
        rec(collect(1:n), 0)
        best
    end
    tested = 0
    for _ in 1:800
        n = 2rand(rng, 1:5)
        dense = rand(rng, Bool)
        W = fill(typemax(Int), n, n)
        edges = Tuple{Int,Int,Int}[]
        for i in 1:n, j in i+1:n
            (dense || rand(rng) < 0.5) || continue
            w = rand(rng, 0:(rand(rng, Bool) ? 3 : 20))          # ties are where bugs hide
            W[i, j] = W[j, i] = w
            push!(edges, (i, j, w))
        end
        opt = brute(n, W)
        if opt == typemax(Int)
            @test_throws ArgumentError min_weight_perfect_matching(n, edges)
            continue
        end
        mate = min_weight_perfect_matching(n, edges)
        tested += 1
        @test all(mate[mate[v]] == v != mate[v] for v in 1:n)
        @test sum(W[v, mate[v]] for v in 1:n if v < mate[v]) == opt
    end
    @test tested > 500
    @test min_weight_perfect_matching(0, Tuple{Int,Int,Int}[]) == Int[]
    @test_throws ArgumentError min_weight_perfect_matching(3, [(1, 2, 1)])
    # a 200-vertex complete graph is quick
    pts = rand(rng, 1:100, 200, 2)
    edges = [(i, j, sum(abs.(pts[i, :] .- pts[j, :]))) for i in 1:200 for j in i+1:200]
    mate = min_weight_perfect_matching(200, edges)
    @test all(mate[mate[v]] == v for v in 1:200)
end
