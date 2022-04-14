# Reordering

Reorder the elemental loop or nodes.
Most of these assume a uniform grid of squares or cubes such as
those produced by the internal mesh creation utilities.

```@index
Pages = ["reorder.md"]
```

```@docs
mortonElements(griddim)
hilbertElements(griddim)
tiledElements(griddim, tiledim)
randomElements(seed = 17)
mortonNodes(griddim)
hilbertNodes(griddim)
tiledNodes(griddim, tiledim)
elementFirstNodes()
randomNodes(seed = 17)
```
