# LRSLib

| **Build Status** | **References to cite** |
|:----------------:|:----------------------:|
| [![Build Status][build-img]][build-url] | [![DOI][zenodo-img]][zenodo-url] |
| [![Codecov branch][codecov-img]][codecov-url] | |

LRSLib.jl is a wrapper for [lrs](http://cgm.cs.mcgill.ca/~avis/C/lrs.html). This module can either be used in a "lower level" using the [API of lrslib](http://cgm.cs.mcgill.ca/~avis/C/lrslib/lrslib.html) or using the higher level interface of [Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl).

As written in the [user guide of lrs](http://cgm.cs.mcgill.ca/~avis/C/lrslib/USERGUIDE.html#Introduction):
> A polyhedron can be described by a list of inequalities (H-representation) or as by a list of its vertices and extreme rays (V-representation). lrs is a C program that converts a H-representation of a polyhedron to its V-representation, and vice versa.  These problems are known respectively at the vertex enumeration and convex hull problems.

I have [forked lrs](https://github.com/blegat/lrslib) to add a few functions to help doing the wrapper.
These changes are not upstream yet so this version is used instead of the upstream version.

[build-img]: https://github.com/JuliaPolyhedra/LRSLib.jl/workflows/CI/badge.svg?branch=master
[build-url]: https://github.com/JuliaPolyhedra/LRSLib.jl/actions?query=workflow%3ACI
[codecov-img]: http://codecov.io/github/JuliaPolyhedra/LRSLib.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/JuliaPolyhedra/LRSLib.jl?branch=master

[zenodo-url]: https://doi.org/10.5281/zenodo.1214579
[zenodo-img]: https://zenodo.org/badge/DOI/10.5281/zenodo.1214579.svg
