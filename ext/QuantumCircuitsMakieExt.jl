module QuantumCircuitsMakieExt

using Makie
using QuantumCircuits
using QuantumCircuits: Circuit, Gate, Instruction, label, nqubits, count_cnots,
                       gray, gray_flip_position, gray_flip_positions

# ---------------------------------------------------------------------------
# Palette
#
# Two categorical hues (indigo / amber), each stepped separately for the light
# and dark surface rather than flipped, and validated for colour-vision
# deficiency separation on both.  Magnitude gets a single-hue sequential ramp;
# signed quantities get a diverging ramp with a neutral grey midpoint; phase
# gets a cyclic ramp, because phase wraps and so must its colour.
# ---------------------------------------------------------------------------

struct Palette
    surface::String
    ink::String
    muted::String
    line::String
    accent::String       # categorical 1
    accent2::String      # categorical 2
    fill::String         # gate box fill
    sequential::Vector{String}
    diverging::Vector{String}
    cyclic::Vector{String}
end

const LIGHT = Palette(
    "#FCFCFB", "#16161F", "#5C5C6E", "#3A3A48", "#4B41B8", "#C87A0E", "#EFEEFA",
    ["#F7F6FD", "#D3CFF3", "#9A93E6", "#4B41B8", "#1F1A57"],
    ["#8A5407", "#DDA24A", "#E9E8E5", "#7C74D8", "#2E2680"],
    ["#241E63", "#8A82E0", "#F4F3FC", "#DDA24A", "#8A5407", "#241E63"],
)

const DARK = Palette(
    "#1A1A19", "#E8E8F0", "#A0A0B4", "#C6C6D4", "#8A82E0", "#B8842F", "#26263A",
    ["#14132B", "#2C2670", "#4B41B8", "#8A82E0", "#DAD7F8"],
    ["#B8842F", "#7A5A22", "#3A3A44", "#5D55C0", "#A9A2F0"],
    ["#DAD7F8", "#8A82E0", "#2C2670", "#7A5A22", "#B8842F", "#DAD7F8"],
)

_palette(theme::Symbol) = theme === :dark ? DARK :
    theme === :light ? LIGHT : throw(ArgumentError("theme must be :light or :dark, got $theme"))

# ---------------------------------------------------------------------------
# Circuit diagrams
# ---------------------------------------------------------------------------

# Greedy column packing: an instruction takes the first column in which every
# wire it spans (including the ones its vertical connector crosses) is free.
function _columns(c::Circuit)
    nfree = ones(Int, c.nqubits)
    cols = Int[]
    for op in c.ops
        lo, hi = extrema(op.qubits)
        col = maximum(view(nfree, lo:hi))
        push!(cols, col)
        for q in lo:hi
            nfree[q] = col + 1
        end
    end
    cols
end

# is this a controlled gate emitted by `controlled(...)`?  name is C…C<base>
function _control_split(g::Gate, nq::Int)
    s = String(g.name)
    k = something(findfirst(ch -> ch != 'C', s), length(s) + 1) - 1
    (k >= 1 && k == nq - 1) ? (k, s[k+1:end]) : (0, "")
end

_boxwidth(txt::AbstractString) = 0.42 + 0.155 * length(txt)

function _opwidth(op::Instruction)
    g = op.gate
    g.name in (:CNOT, :CZ, :SWAP) && return 0.5
    nq = length(op.qubits)
    nq == 1 && return _boxwidth(label(g))
    k, base = _control_split(g, nq)
    k > 0 ? _boxwidth(base) : _boxwidth(label(g))
end

"""
    circuitplot!(ax, c; kwargs...)

Draw circuit `c` into `ax`.  See `QuantumCircuits.circuitplot!`.
"""
function QuantumCircuits.circuitplot!(ax, c::Circuit;
                                      theme::Symbol=:light,
                                      wirelabels::Bool=true,
                                      fontsize::Real=13,
                                      p::Palette=_palette(theme))
    n = c.nqubits
    cols = _columns(c)
    ncol = isempty(cols) ? 0 : maximum(cols)

    # column geometry: each column is as wide as its widest instruction
    widths = fill(0.5, max(ncol, 1))
    for (op, col) in zip(c.ops, cols)
        widths[col] = max(widths[col], _opwidth(op))
    end
    pad = 0.30
    edges = Vector{Float64}(undef, max(ncol, 1) + 1)
    edges[1] = 0.0
    for i in 1:max(ncol, 1)
        edges[i+1] = edges[i] + widths[i] + pad
    end
    centre(col) = (edges[col] + edges[col+1] - pad) / 2
    x0, x1 = -0.35, edges[end] - pad + 0.35

    ink, muted, accent = parse.(Makie.Colorant, (p.ink, p.muted, p.accent))
    linec, fillc = parse.(Makie.Colorant, (p.line, p.fill))
    surf = parse(Makie.Colorant, p.surface)

    for q in 1:n                                            # the wires
        lines!(ax, [x0, x1], [q, q]; color=linec, linewidth=1.1)
    end
    if wirelabels
        for q in 1:n
            text!(ax, Point2f(x0 - 0.12, q); text="q$q", align=(:right, :center),
                  fontsize=fontsize, color=muted)
        end
    end

    dot!(x, q) = scatter!(ax, [x], [q]; color=ink, markersize=11, marker=:circle)
    function target!(x, q)                                  # the ⊕ of a CNOT
        scatter!(ax, [x], [q]; color=surf, strokecolor=ink, strokewidth=1.6,
                 markersize=17, marker=:circle)
        r = 0.085
        lines!(ax, [x - r, x + r], [q, q]; color=ink, linewidth=1.6)
        lines!(ax, [x, x], [q - r * 2.2, q + r * 2.2]; color=ink, linewidth=1.6)
    end
    function box!(x, q, txt; w=_boxwidth(txt))
        poly!(ax, Rect2f(x - w/2, q - 0.30, w, 0.60);
              color=fillc, strokecolor=accent, strokewidth=1.4)
        text!(ax, Point2f(x, q); text=txt, align=(:center, :center),
              fontsize=fontsize, color=ink)
    end

    for (op, col) in zip(c.ops, cols)
        x = centre(col)
        g, qs = op.gate, op.qubits
        lo, hi = extrema(qs)
        length(qs) > 1 && lines!(ax, [x, x], [lo, hi]; color=ink, linewidth=1.3)
        if g.name === :CNOT
            dot!(x, qs[1]); target!(x, qs[2])
        elseif g.name === :CZ
            dot!(x, qs[1]); dot!(x, qs[2])
        elseif g.name === :SWAP
            for q in qs
                scatter!(ax, [x], [q]; color=ink, markersize=13, marker=:xcross)
            end
        elseif length(qs) == 1
            box!(x, qs[1], label(g))
        else
            k, base = _control_split(g, length(qs))
            if k > 0
                for q in qs[1:k]; dot!(x, q); end
                box!(x, qs[end], base)
            else
                for (i, q) in enumerate(qs)
                    box!(x, q, string(label(g), "[", i, "]"))
                end
            end
        end
    end

    limits!(ax, x0 - (wirelabels ? 0.75 : 0.15), x1 + 0.15, n + 0.75, 0.25)
    ax
end

"""
    circuitfigure(c; kwargs...)

Publication-quality circuit diagram.  See `QuantumCircuits.circuitfigure`.
"""
function QuantumCircuits.circuitfigure(c::Circuit;
                                       theme::Symbol=:light,
                                       title=nothing,
                                       scale::Real=1,
                                       fontsize::Real=13,
                                       wirelabels::Bool=true)
    p = _palette(theme)
    cols = _columns(c)
    ncol = isempty(cols) ? 1 : maximum(cols)
    w = scale * (140 + 62 * ncol)
    h = scale * (70 + 46 * c.nqubits)
    ttl = if title !== nothing
        title
    else
        base = "$(c.nqubits) qubits · $(length(c.ops)) gates · $(count_cnots(c)) CNOTs"
        c.global_phase == 0 ? base :
            base * " · global phase " * string(round(c.global_phase; digits=4)) * " rad"
    end
    fig = Figure(size=(w, h), backgroundcolor=parse(Makie.Colorant, p.surface))
    ax = Axis(fig[1, 1]; title=String(ttl), titlesize=fontsize + 2,
              titlecolor=parse(Makie.Colorant, p.muted), titlealign=:left,
              backgroundcolor=parse(Makie.Colorant, p.surface))
    hidedecorations!(ax); hidespines!(ax)
    QuantumCircuits.circuitplot!(ax, c; theme=theme, fontsize=fontsize,
                                 wirelabels=wirelabels, p=p)
    fig
end

# ---------------------------------------------------------------------------
# Matrix heatmaps
# ---------------------------------------------------------------------------

function QuantumCircuits.matrixfigure(U::AbstractMatrix;
                                      part::Symbol=:abs,
                                      theme::Symbol=:light,
                                      title=nothing,
                                      atol::Real=1e-9,
                                      scale::Real=1,
                                      labels::Bool=true)
    p = _palette(theme)
    N = size(U, 1)
    A = ComplexF64.(U)
    if part === :abs
        vals = abs.(A); cmap = Makie.cgrad(p.sequential); rng = (0.0, maximum(vals))
        cblabel = "|Uᵢⱼ|"
    elseif part === :real || part === :imag
        vals = part === :real ? real.(A) : imag.(A)
        m = max(maximum(abs, vals), eps())
        cmap = Makie.cgrad(p.diverging); rng = (-m, m)
        cblabel = part === :real ? "Re Uᵢⱼ" : "Im Uᵢⱼ"
    elseif part === :phase
        vals = [abs(z) > atol ? angle(z) : NaN for z in A]
        cmap = Makie.cgrad(p.cyclic); rng = (-Float64(π), Float64(π))
        cblabel = "arg Uᵢⱼ"
    else
        throw(ArgumentError("part must be :abs, :real, :imag or :phase, got $part"))
    end

    surf = parse(Makie.Colorant, p.surface)
    fig = Figure(size=(scale * (300 + 15N), scale * (250 + 15N)), backgroundcolor=surf)
    ax = Axis(fig[1, 1]; title=title === nothing ? String(part) : String(title),
              titlealign=:left, titlecolor=parse(Makie.Colorant, p.muted),
              aspect=DataAspect(), yreversed=true, backgroundcolor=surf,
              xlabel="column (input basis state)", ylabel="row (output basis state)",
              xlabelcolor=parse(Makie.Colorant, p.muted),
              ylabelcolor=parse(Makie.Colorant, p.muted))
    hm = heatmap!(ax, 1:N, 1:N, permutedims(vals); colormap=cmap, colorrange=rng,
                  nan_color=parse(Makie.Colorant, p.fill))
    if labels && N <= 16
        n = trailing_zeros(N)
        ticks = (1:N, [string(j; base=2, pad=n) for j in 0:N-1])
        ax.xticks = ticks; ax.yticks = ticks
        ax.xticklabelrotation = π/2
        ax.xticklabelsize = 9; ax.yticklabelsize = 9
        ax.xticklabelcolor = parse(Makie.Colorant, p.muted)
        ax.yticklabelcolor = parse(Makie.Colorant, p.muted)
    else
        hidedecorations!(ax; label=false)
    end
    Colorbar(fig[1, 2], hm; label=cblabel, labelcolor=parse(Makie.Colorant, p.muted),
             ticklabelcolor=parse(Makie.Colorant, p.muted))
    fig
end

QuantumCircuits.matrixfigure(c::Circuit; kwargs...) =
    QuantumCircuits.matrixfigure(QuantumCircuits.matrix(c); kwargs...)

# ---------------------------------------------------------------------------
# Gray code walk
# ---------------------------------------------------------------------------

function QuantumCircuits.graycodefigure(n::Integer=4; theme::Symbol=:light, scale::Real=1)
    p = _palette(theme)
    N = 1 << n
    seq = [gray(i) for i in 0:N-1]
    flips = gray_flip_positions(n)          # bit flipped entering step i (cyclic)

    surf = parse(Makie.Colorant, p.surface)
    ink = parse(Makie.Colorant, p.ink)
    muted = parse(Makie.Colorant, p.muted)
    accent = parse(Makie.Colorant, p.accent)
    accent2 = parse(Makie.Colorant, p.accent2)
    fillc = parse(Makie.Colorant, p.fill)

    fig = Figure(size=(scale * (160 + 34N), scale * (130 + 34n)), backgroundcolor=surf)
    ax = Axis(fig[1, 1]; title="$n-bit Gray code — one bit changes per step",
              titlealign=:left, titlecolor=muted, backgroundcolor=surf, yreversed=true)
    hidedecorations!(ax); hidespines!(ax)

    for i in 1:N, b in 0:n-1
        set = (seq[i] >> b) & 1 == 1
        y = n - b                            # bit n-1 on top, bit 0 at the bottom
        poly!(ax, Rect2f(i - 0.42, y - 0.42, 0.84, 0.84);
              color=set ? accent : fillc, strokecolor=muted, strokewidth=0.4)
        text!(ax, Point2f(i, y); text=set ? "1" : "0", align=(:center, :center),
              fontsize=11, color=set ? surf : muted)
    end
    # ring the bit that changed on entry to each step
    for i in 2:N
        y = n - flips[i-1]
        poly!(ax, Rect2f(i - 0.46, y - 0.46, 0.92, 0.92);
              color=(:white, 0.0), strokecolor=accent2, strokewidth=2.4)
    end
    for b in 0:n-1
        text!(ax, Point2f(0.25, n - b); text="bit $b", align=(:right, :center),
              fontsize=11, color=muted)
    end
    for i in 1:N
        text!(ax, Point2f(i, n + 0.85); text=string(i - 1), align=(:center, :center),
              fontsize=10, color=muted)
    end
    text!(ax, Point2f(0.25, n + 0.85); text="i", align=(:right, :center), fontsize=11, color=muted)
    limits!(ax, -1.6, N + 0.8, n + 1.4, 0.2)

    ax2 = Axis(fig[2, 1]; backgroundcolor=surf, xlabel="step i   (i = 2ⁿ closes the cycle)",
               ylabel="flipped bit", xlabelcolor=muted, ylabelcolor=muted,
               xticklabelcolor=muted, yticklabelcolor=muted,
               yticks=(0:n-1, ["$b" for b in 0:n-1]))
    hidespines!(ax2, :t, :r)
    xs = 1:N
    for (x, f) in zip(xs, flips)
        lines!(ax2, [x, x], [-0.35, f]; color=ink, linewidth=1.0)
    end
    scatter!(ax2, collect(xs), flips; color=accent2, markersize=9)
    limits!(ax2, 0.2, N + 0.8, -0.5, n - 0.5)
    rowsize!(fig.layout, 2, Relative(0.32))
    fig
end

# ---------------------------------------------------------------------------
# Cost comparison
# ---------------------------------------------------------------------------

function QuantumCircuits.costfigure(ns=2:9; theme::Symbol=:light, scale::Real=1)
    p = _palette(theme)
    xs = collect(ns)
    naive = [n * (1 << n) for n in xs]                       # a multi-controlled
    prep = [2 * (1 << n) - 4 for n in xs]                    # gate per branch
    diag = [(1 << n) - 2 for n in xs]
    mux = [1 << (n - 1) for n in xs]

    surf = parse(Makie.Colorant, p.surface)
    muted = parse(Makie.Colorant, p.muted)
    accent = parse(Makie.Colorant, p.accent)
    accent2 = parse(Makie.Colorant, p.accent2)

    fig = Figure(size=(scale * 620, scale * 400), backgroundcolor=surf)
    ax = Axis(fig[1, 1]; title="CNOT count vs register size", titlealign=:left,
              titlecolor=muted, backgroundcolor=surf, yscale=log10,
              xlabel="qubits n", ylabel="CNOTs (log scale)",
              xlabelcolor=muted, ylabelcolor=muted,
              xticklabelcolor=muted, yticklabelcolor=muted,
              xticks=xs, ygridcolor=(muted, 0.18), xgridvisible=false)
    hidespines!(ax, :t, :r)

    series = [("naive, one multi-controlled gate per branch", naive, accent2, :dash),
              ("state preparation, 2ⁿ⁺¹-4", prep, accent, :solid),
              ("diagonal unitary, 2ⁿ-2", diag, accent, :dot),
              ("multiplexed rotation, 2ⁿ⁻¹", mux, accent, :dashdot)]
    for (lbl, ys, col, ls) in series
        lines!(ax, xs, ys; color=col, linewidth=2, linestyle=ls, label=lbl)
        scatter!(ax, xs, ys; color=col, markersize=8)
    end
    axislegend(ax; position=:lt, framevisible=false, labelcolor=muted, labelsize=11)
    fig
end

end # module
