
library("rhdf5")
f <- "example.hdf5"

trees <- h5read(f, "trees")
L <- order(trees$left, trees$time)
R <- order(trees$right, -trees$time)
p <- rep(0, max(trees$node))
N <- length(trees$left)

l <- 0
k <- 1
j <- 1
while (k <= N) {
    while (j <= N && trees$right[R[j]] == l) {
        u <- trees$node[R[j]]
        p[trees$children[1, R[j]]] <- 0
        p[trees$children[2, R[j]]] <- 0
        j <- j + 1
    }
    while (k <= N && trees$left[L[k]] == l) {
        u <- trees$node[L[k]]
        p[trees$children[1, L[k]]] <- u
        p[trees$children[2, L[k]]] <- u
        k <- k + 1
    }
    r <- trees$right[R[j]]
    # The new tree is ready -- this is where we'd do something with
    # it.
    print(paste("New tree: left = ", l, "right = ", r, "p = (next line)"))
    print(p)
    l <- r 
}


