rectangular_conformal_map: Conformally map a simply-connected open triangle mesh to a rectangle

This code computes the rectangular conformal parameterizations (i.e. angle-preserving mappings onto a rectangle) of triangle meshes with disk topology using the fast method in [1], which has been applied for texture mapping, surface registration, shape analysis and so on.
Any comments and suggestions are welcome. 

If you use this code in your own work, please cite the following paper:
[1] T. W. Meng, G. P.-T. Choi and L. M. Lui, 
    "TEMPO: Feature-Endowed Teichmüller Extremal Mappings of Point Clouds."
    SIAM Journal on Imaging Sciences, 9(4), pp. 1922-1962, 2016.

Copyright (c) 2015-2018, Gary Pui-Tung Choi
https://math.mit.edu/~ptchoi

===============================================================



Usage:
map = rectangular_conformal_map(v,f,corner)


Input:
v: nv x 3 vertex coordinates of a simply-connected open triangle mesh
f: nf x 3 triangulations of a simply-connected open triangle mesh
(optional) corner: 4 x 1 vertex indices for the four corners of the rectangle, with anti-clockwise orientation
 4 - 3
 |   |
 1 - 2

Output:
map: nv x 2 vertex coordinates of the rectangular conformal map

Remarks:
1. Please make sure that the input mesh does not contain any unreferenced vertices/non-manifold vertices/non-manifold edges.
2. Please remove all valence 1 boundary vertices (i.e. vertices with only 1 face attached to them) before running the program.
3. Please make sure that the input triangulations f are with anti-clockwise orientation.
4. The output rectangular domain will always have width = 1, while the height depends on the choice of the corners and may not be 1.
   (The Riemann mapping theorem guarantees that there exists a conformal map from any simple-connected open surface to the unit square, but if four vertices on the surface boundary are specified to be the four corners of the planar domain, the theorem is no longer applicable.)
