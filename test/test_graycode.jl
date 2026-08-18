@testset "gray code" begin
    @test gray.(0:7) == [0, 1, 3, 2, 6, 7, 5, 4]

    for n in 1:10
        seq = graycode(n)
        @test length(seq) == 1 << n
        @test sort(seq) == collect(0:(1 << n)-1)              # bijection
        @test all(gray_adjacent(seq[i], seq[i+1]) for i in 1:length(seq)-1)
        @test gray_adjacent(seq[end], seq[1])                 # cyclic
        @test all(ungray(gray(i)) == i for i in 0:(1 << n)-1) # inverse
    end

    @testset "flip positions" begin
        for k in 1:8
            fl = gray_flip_positions(k)
            @test length(fl) == 1 << k
            # cumulative XOR of flip masks reproduces the Gray sequence
            acc = 0
            for i in 1:(1 << k)-1
                acc ⊻= 1 << fl[i]
                @test acc == gray(i)
            end
            acc ⊻= 1 << fl[end]
            @test acc == 0                       # the cycle closes: CNOTs cancel
            @test fl[end] == k - 1               # last flip is the top control
            @test all(0 .<= fl .<= k - 1)
        end
    end

    @testset "walks and helpers" begin
        @test hamming(0b1011, 0b1110) == 2
        @test sort(gray_walk(0b1011, 0b1110)) == [0, 2]
        @test isempty(gray_walk(5, 5))
        @test bits(5, 4) == "0101"
        @test parity(0b101, 0b111) == 0
        @test parity(0b100, 0b111) == 1
    end
end
