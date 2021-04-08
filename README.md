# Koebe Realizations

Computes a Koebe realization of a 3-dimensional polytope following [Variational principles for circle patterns and Koebe's theorem](https://arxiv.org/abs/math/0203250). The input is the vertex-facet description of a polytope.

## Installation

```julia
using Pkg
pkg"add https://github.com/saschatimme/KoebeRealizations.jl.git"
```


## Example 1: The Cube
```julia
using KoebeRealizations

cube = Bool[
    1 1 1 1 0 0 0 0
    1 1 0 0 1 1 0 0
    0 1 1 0 0 1 1 0
    0 0 1 1 0 0 1 1
    1 0 0 1 1 0 0 1
    0 0 0 0 1 1 1 1
]
koebe_realization(cube)
```

```
8×3 Matrix{Float64}:
  0.316228  -0.707107   0.948683
 -0.948684  -0.707106   0.316228
 -0.316228  -0.707107  -0.948683
  0.948683  -0.707107  -0.316228
  0.316228   0.707107   0.948683
 -0.948683   0.707107   0.316228
 -0.316228   0.707107  -0.948683
  0.948683   0.707107  -0.316228
```

## Example 2: The Icosahedron

The input can also be a polytope from `Polymake.jl`. In this case the function also returns a polymake polytope.

```julia
using KoebeRealizations, Polymake

ico = polytope.icosahedron()
Q = koebe_realization(ico)
Polymake.visual(Q)
```

## Example 3: The Circle Packing

Underneath, the algorithm to compute a Koebe realization computes a particular circle packing. This can be visualized as follows:

```julia
using KoebeRealizations
cube = Bool[
    1 1 1 1 0 0 0 0
    1 1 0 0 1 1 0 0
    0 1 1 0 0 1 1 0
    0 0 1 1 0 0 1 1
    1 0 0 1 1 0 0 1
    0 0 0 0 1 1 1 1
]
plot_circle_packing(cube)
```
