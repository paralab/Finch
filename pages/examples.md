---
title: Examples
---

## Examples

These examples start simple and demonstrate how to use the various aspects of finch.

<style>
img {float: right;}
</style>
<h3>Poisson</h3> 
<p> <img src="images/poisson1d.png" alt="poisson1d" width="200">
The simplest example possible. 1D Poisson with homogeneous Dirichlet boundary. It demonstrates the basics of setting up and solving a problem directly in finch.
<br>Page: <a href="https://paralab.github.io/finch/pages/poisson1d.html">poisson 1D page</a>
<br>Script: <a href="https://github.com/paralab/finch/blob/master/src/examples/example-poisson1d.jl">example-poisson1d.jl</a>
</p>

<h3>Linear Elasticity</h3>
<p> <img src="images/elasticity.png" alt="elasticity" width="200">
The linear elasticity example demonstrates vector entities and mixed boundary conditions.
<br>Page: <a href="https://paralab.github.io/finch/pages/elasticity.html">elasticity page</a>
<br><a href="https://github.com/paralab/finch/blob/master/src/examples/example-elasticity.jl">example-elasticity.jl</a>
</p>

<h3>Heat</h3>
<p> <img src="images/heat.png" alt="heat" width="200">
The heat equation demonstrates support for time dependent problems.
<br>Page: <a href="https://paralab.github.io/finch/pages/heat.html">heat page</a>
<br><a href="https://github.com/paralab/finch/blob/master/src/examples/example-heat2d.jl">example-heat2d.jl</a>
</p>

<h3>Navier Stokes</h3>
<p> <img src="images/NS.png" alt="NS" width="150">
This 2D lid-driven NS problem demonstrates some of the nonlinear capabilities and the use of parameter entities.
<br>Page: <a href="https://paralab.github.io/finch/pages/NS.html">NS page</a>
<br><a href="https://github.com/paralab/finch/blob/master/src/examples/example-NS.jl">example-NS.jl</a>
</p>

<h3>Unstructured Meshes</h3>
<p><img src="images/umeshes.png" alt="unstructured" width="200">
Unstructured meshes in 2D and 3D made of triangles, quadrilaterals, hexahedra and tetrahedra are imported from .MSH files and demonstrated on some simple problems.
<br>Page: <a href="https://paralab.github.io/finch/pages/unstructured.html">Unstructured Meshes page</a>
<br><a href="https://github.com/paralab/finch/blob/master/src/examples/example-unstructured2d.jl">example-unstructured2d.jl</a>
<br><a href="https://github.com/paralab/finch/blob/master/src/examples/example-unstructured3d.jl">example-unstructured3d.jl</a>
</p>

<h3>Finite Volume: Advection</h3>
<p><img src="images/fvadvec2d.png" alt="fvadvec2d" width="200">
Finite volume method for a basic advection problem. This example demonstrates the finite volume capability as well as exporting/importing generated code.
<br>Page: <a href="https://paralab.github.io/finch/pages/FVadvection.html">FV Advection page</a>
<br><a href="https://github.com/paralab/finch/blob/master/src/examples/example-advection1d-fv.jl">example-advection1d-fv.jl</a>
<br><a href="https://github.com/paralab/finch/blob/master/src/examples/example-advection2d-fv.jl">example-advection2d-fv.jl</a>
</p>

<h3>Indexed Variables</h3>
<p><img src="images/addiff1dindexed.png" alt="addiff1dindexed" width="200">
Finite volume method for an advection-diffusion problem using indexed variables to compute the same equation over many values of advection speed and diffusion rate. This example demonstrates indexed variables and customizable assembly code generation.
<br>Page: <a href="https://paralab.github.io/finch/pages/indexed.html">Indexed Variables page</a>
<br><a href="https://github.com/paralab/finch/blob/master/src/examples/example-addiff1d-indexed.jl">example-addiff1d-indexed.jl</a>
</p>