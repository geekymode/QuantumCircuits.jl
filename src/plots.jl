# ---------------------------------------------------------------------------
# Plotting API
#
# The implementations live in `ext/QuantumCircuitsMakieExt.jl`, a package
# extension loaded automatically as soon as a Makie backend is available.  The
# stubs here exist so the names are documented and so the error you get without
# a backend tells you what to do.
# ---------------------------------------------------------------------------

const _MAKIE_HINT = """
requires a Makie backend. Run

    using Pkg; Pkg.add("CairoMakie")   # once
    using CairoMakie                   # loads the QuantumCircuits plotting extension

CairoMakie renders publication-quality vector output (`save("fig.pdf", fig)`);
GLMakie or WGLMakie work too if you want an interactive window."""

"""
    circuitfigure(c; kwargs...) -> Figure

Draw a circuit as a publication-quality diagram.

Instructions are packed into columns greedily, so gates on disjoint wires share
a column the way a hand-drawn diagram would. Keyword arguments:

* `title` — figure title (default: a gate/CNOT count summary)
* `theme` — `:light` (default) or `:dark`
* `wirelabels` — `true` to label the wires `q1 … qn`
* `gatecolor`, `wirecolor`, `accent` — override individual colours
* `scale` — overall size multiplier

$_MAKIE_HINT
"""
circuitfigure(args...; kwargs...) = error("`circuitfigure` ", _MAKIE_HINT)

"""
    circuitplot!(ax, c; kwargs...) -> ax

Draw a circuit into an existing `Makie.Axis`, for composing several circuits
into one figure. Same keywords as [`circuitfigure`](@ref).

$_MAKIE_HINT
"""
circuitplot!(args...; kwargs...) = error("`circuitplot!` ", _MAKIE_HINT)

"""
    matrixfigure(U; part=:abs, kwargs...) -> Figure

Heatmap of a matrix (or of `matrix(c)` for a circuit), with a colourbar.

`part` selects what is shown and, with it, the kind of colour scale:

* `:abs` — magnitude, on a single-hue sequential ramp (light → dark);
* `:real`, `:imag` — signed, on a diverging ramp with a neutral grey midpoint;
* `:phase` — the argument, on a cyclic ramp (phase wraps, so the colour must
  wrap with it), masked to entries whose magnitude clears `atol`.

Block structure is the point: a multiplexor is block-diagonal, a diagonal
unitary is a single diagonal line, and a permutation is a scatter of ones.

$_MAKIE_HINT
"""
matrixfigure(args...; kwargs...) = error("`matrixfigure` ", _MAKIE_HINT)

"""
    graycodefigure(n; kwargs...) -> Figure

Illustrate the `n`-bit Gray code: the bit grid over the whole walk, with the
bit that changes at each step highlighted, and the flip-position sequence
beneath it.

$_MAKIE_HINT
"""
graycodefigure(args...; kwargs...) = error("`graycodefigure` ", _MAKIE_HINT)

"""
    costfigure(ns=1:8; kwargs...) -> Figure

CNOT count versus register size for the package's decompositions against the
naive `O(n 2ⁿ)` compilation, on a log scale.

$_MAKIE_HINT
"""
costfigure(args...; kwargs...) = error("`costfigure` ", _MAKIE_HINT)
