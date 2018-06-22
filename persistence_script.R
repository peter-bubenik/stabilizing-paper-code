# Example from Section 1.2.1 of 
# Stabilizing the unstable output of persistent homology computations
# by Paul Bendich, Peter Bubenik, and Alexander Wagner

# Copyright 2018, Peter Bubenik, All rights reserved.

FIND <- function(i){
  if (P[i] != i){
    P[i] <<- FIND(P[i]) # superassignment to change global variable
  }
  return(P[i])
}

UNION <- function(i,j,birth_vertex){
  k <- FIND(i)
  l <- FIND(j)
  if (k != l){
    if (S[k] < S[l]){
      m <- k
      k <- l
      l <- m
    }
    P[l] <<- k
    B[k] <<- birth_vertex
    S[k] <<- S[k] + S[l]
  }
}

# calculate the degree 0 persistence diagram for one of the brain artery tree
# input t is the subject number from 1 to 97
# output is a matrix, each row is birth value, death value, birth vertex number, and 
# death vertex number
# followed by persistence (death-birth), and the x,y,z coordinates of the birth vertex
# sorted by persistence in decreasing order

pd_list <- vector("list",num_repeats)
bar_of_interest_list <- vector("list",num_repeats)

for (c in 1:num_repeats){

V <- dim(vertices_list[[t]])[1]
E <- dim(edges_list[[t]])[1]
vertex_values <- vertices_list[[t]][,1:3] + cbind(rnorm(V,sd=sd),rnorm(V,sd=sd),rnorm(V,sd=sd))

edge_values <- matrix(0,E,3)
for (i in 1:E){ 
  edge_values[i,1] <- i
  edge_values[i,2] <- max(vertex_values[edges_list[[t]][i,1],3], vertex_values[edges_list[[t]][i,2],3])
  edge_values[i,3] <- min(vertex_values[edges_list[[t]][i,1],3], vertex_values[edges_list[[t]][i,2],3])
}
edge_values <- edge_values[order(edge_values[,2],edge_values[,3]),] # sort edges

P <- seq(1,V) # parent of vertex (vertex number of root of connected component)
B <- seq(1,V) # birth vertex (vertex number of oldest vertex in connected component - only correct for root vertex)
S <- rep(1,V) # size (number of vertices in connected component - only correct for root vertex)

pd <- vector("list", length = 8*V)
bd_pairs <- 0

for (i in 1:E){
  e <- edge_values[i,1]
  i <- edges_list[[t]][e,1]
  j <- edges_list[[t]][e,2]
  k <- FIND(i)
  l <- FIND(j)
  height_i <- vertices_list[[t]][i,3]
  height_j <- vertices_list[[t]][j,3]
  height_k <- vertices_list[[t]][k,3]
  height_l <- vertices_list[[t]][l,3]
  death_val <- max(height_i,height_j)
  birth_val <- max(height_k,height_l)
  if (height_k < height_l){
    birth_vertex <- l
  } else{
    birth_vertex <- k
  }
  if (birth_val != death_val){
    if (height_i < height_j){
      death_vertex <- j
    } else{
      death_vertex <- i
    }
    x <- vertices_list[[t]][birth_vertex,1]
    y <- vertices_list[[t]][birth_vertex,2]
    z <- vertices_list[[t]][birth_vertex,3]
    pd[(bd_pairs*8+1):(bd_pairs*8+8)] <- c(birth_val,death_val,birth_vertex,death_vertex,death_val-birth_val,x,y,z)
    bd_pairs <- bd_pairs+1
  }
  UNION(i,j,birth_vertex)
}

pd <- unlist(pd)
pd <- matrix(pd, ncol = 8, byrow = TRUE)
pd <- pd[order(pd[,5], decreasing=TRUE),]
bar_of_interest_list[[c]] <- pd[28,]
pd_list[[c]] <- pd
}

bl <- unlist(bar_of_interest_list)
bars_of_interest <- matrix(bl, ncol = 8, byrow = TRUE)