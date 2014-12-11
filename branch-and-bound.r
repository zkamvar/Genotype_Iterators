# Way to manipulate the output of the C permutation to be able to use it for
# branch and bound search

library(poppr)

# Interesting thing: if the permutation matrix is ordered, then you would expect
# the slowest changing columns to change every (n*n!)/(n^2) rows.

sapply(1:10, function(x) length(.Call("permuto", x))/(x^2))

# This lends itself to branch and bound search.
# Unfortunately, the output of .Call("permuto", n) is in a reversed order.

.Call("permuto", 4)

# So, we write a little function to return an ordered matrix

permuto <- function(x) matrix(rev(.Call("permuto", x)), nrow = x)
permuto(4)

# That's much better, but now the indices are still reversed. We can fix this by
# using the apply function:

t(apply(permuto(4), 1, rev))

# This, will take more time, however. A lot more time.