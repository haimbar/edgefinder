# edgefinder
The edgefinder package is used to find edges in gene networks using co-expression
data. The input to the program is a normalized expression matrix, with genes (nodes)
in the rows, and samples in the columns.
The program calculates the pair-wise correlations, performs Fisher's Z 
transformation, and fits the L2N model to the transformed data. L2N is a mixture
model with three components: the uncorrelated pairs belong to the null component
which is assumed to be normally distributed, and the correlated pairs belong to one
of the two non-null components which are assumed to follow lognormal distributions.

Typical datasets consist of hundreds, or thousands of genes, and hence a very
large number of pairs. Therefore, edgefinder randomly selects a subset of the pairs (the
default number of pairs is 20,000), fits the L2N model to the subset, and calculates
the component probabilities for *all* possible pairs.
Using the posterior probabilities, edgefinder determines which pairs are
highly correlated while controlling the false discovery rate.
Note that edgefinder makes no assumptions about the structure of the network.

The edgefinder package depends on the 'Matrix' package, to allow for efficient
storage and computation of large co-occurrence matrices. For simulating datasets
we used the 'huge' and 'MASS' packages, but they are not required when
using edgefinder.

See vignettes/edgefinder.md for details.
