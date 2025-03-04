# 4k_SlidedPolygon4COCH
Efficient Real-Time and Parallel Algorithm for Connected Orthogonal Convex Hulls on Large Point Sets Using $(4k+4)$-Sided Orthogonal Polygons


1. 2ksidedPolygon4OCH.jl: It is a sequential algorithm that uses $(4k+4)$-Sided Orthogonal Polygons to find the Connected Orthogonal Convex Hulls.
2. 2ksidedPolygon4OCH_parallel_ver1_1.jl: It is the Parallel Algorithm for Connected Orthogonal Convex Hulls Using $(4k+4)$-Sided Orthogonal Polygons. This version is used as the result of the paper "Efficient Real-Time and Parallel Algorithm for Connected Orthogonal Convex Hulls on Large Point Sets Using $(4k+4)$-Sided Orthogonal Polygons".
3. O_Quickhull.jl, O_Graham.jl: The algorithms are O-Quickhull (introduced in 2022 by Linh et al. in Applied Mathematics and Computation, 429) and O-Graham (introduced in 2021 by An et al. in Applied Mathematics and Computation, 397).
2ksidedPolygon4OCH_parallel_ver1_2.jl, 2ksidedPolygon4OCH_parallel_ver2_1.jl, 2ksidedPolygon4OCH_parallel_ver2_2.jl: These are several other parallel versions of 2ksidedPolygon4OCH.jl.
utils.jl: This file generates several types of test data and includes some other preparation functions.
main.jl: These files are included in main.jl to run the algorithms.

### Install libraries
- Open a terminal and go to the directory containing the code
>        julia install.jl

**Run the programs**
- Open a terminal and go to the directory containing the code.
> julia main.jl

- Run Algorithms in parallel mode.
> julia -t numberOfThreads main.jl

**Note**
Creat a file "result" in the directory containing the codes.

**Setting**
Benchmarking mode
Set benchmarking = true in the main functions.

Export the convex hull to file
Set benchmarking = false and exportResult = true in the main functions.
