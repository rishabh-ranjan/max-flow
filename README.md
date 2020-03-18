# Fast approximate algorithms for the Maximum Flow problem

Algorithms are described in https://people.csail.mit.edu/madry/docs/maxflow.pdf  
Laplacian solver is here: https://github.com/danspielman/Laplacians.jl/tree/master/src  
SparseMatrixCSC documentation: https://docs.julialang.org/en/v1/stdlib/SparseArrays/  

Points to keep in mind:  
* julia indices are 1-indexed
* Accessing the sparse matrix by rows is considerably slower.
Operations such as insertion of previously unstored entries one at a time in the CSC structure tend to be slow.

Done by Rishabh Ranjan and Kalash Gupta as a semester project under Prof. Amit Kumar.
