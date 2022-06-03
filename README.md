# SpectralClustering
Symmetric matrices' eigenvalues, eigenvectors.


I made this spectral clustering algorithm for a school project to draw a  graph with 50 nodes and 61 edges with the least intersections but it calculates too slow for that.
It has a jacobi algorithm for symmetric matrices (that gives back a diagonal matrix which has all it's eigenvalues in the diagonal), and a Gauss elimination to get the eigenvectors. 

If you want to use it for clustering you can get the second smallest eigenvalue and with that you can get the Fiedler vector. In the Fiedler vector the positive components make a cluster and the negative components make a cluster. (The indexes in the vector are the same as the indexes in the matrix therefor they index the nodes)

It could use some more math but I currently have no time for it.
