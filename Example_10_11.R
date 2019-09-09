# Aux functions
sample.annulus <- function(n, R_in = 20, R_out = 30){
  ptcloud <- matrix(0, ncol = 2, nrow = n)
  for(t in 1:n){
    while(T){
      temp <- runif(2, -R_out, R_out)
      d1 <- sqrt(temp[1]^2 + temp[2]^2)
      if(d1 <= R_out & d1 >= R_in){
        ptcloud[t,] <- temp
        break
      }
    }
  }
  return(ptcloud)
}
pers_hom_filtered_graph <- function(edges,vertex_values){
  # Peter Bubenik, June 11, 2019
  
  # calculate the degree 0 persistence diagram for a filtered graph
  # using the union-find algorithm
  # input: 
  #   edges: an E by 2 matrix, each row is an edge whose entries are the vertex numbers of its boundary
  #   vertex_values: a vector of length V, entries are the scalar values (birth times) of each vertex
  # output: an n by 2 matrix, each row is the birth value and death value or a persistence pair
  
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
      if (vertex_values[k] > vertex_values[l]){
        m <- k
        k <- l
        l <- m
      }
      P[l] <<- k
      B[k] <<- birth_vertex
      S[k] <<- S[k] + S[l]
    }
  }
  
  V <- length(vertex_values)
  E <- dim(edges)[1]
  
  edge_values <- matrix(0,E,3)
  for (i in 1:E){ 
    edge_values[i,1] <- i
    edge_values[i,2] <- max(vertex_values[edges[i,1]], vertex_values[edges[i,2]])
    edge_values[i,3] <- min(vertex_values[edges[i,1]], vertex_values[edges[i,2]])
  }
  edge_values <- edge_values[order(edge_values[,2],-edge_values[,3]),] # sort edges
  
  P <- seq(1,V) # parent of vertex (vertex number of root of connected component)
  B <- seq(1,V) # birth vertex (vertex number of oldest vertex in connected component - only correct for root vertex)
  S <- rep(1,V) # size (number of vertices in connected component - only correct for root vertex)
  
  pd <- vector("list", length = 2*V)
  bd_pairs <- 0
  
  c <- 1
  for (c in 1:E){
    e <- edge_values[c,1]
    i <- edges[e,1]
    j <- edges[e,2]
    k <- FIND(i)
    l <- FIND(j)
    death_val <- max(vertex_values[i],vertex_values[j])
    birth_val <- max(vertex_values[k],vertex_values[l])
    if (vertex_values[k] < vertex_values[l]){
      birth_vertex <- k
    } else{
      birth_vertex <- l
    }
    if (birth_val != death_val){
      if (vertex_values[i] < vertex_values[j]){
        death_vertex <- j
      } else{
        death_vertex <- i
      }
      pd[(bd_pairs*2+1):(bd_pairs*2+2)] <- c(birth_val,death_val) #,birth_vertex,death_vertex,death_val-birth_val,x,y,z)
      bd_pairs <- bd_pairs+1
    }
    UNION(i,j,birth_vertex)
  }
  
  pd <- unlist(pd)
  if(is.null(pd)){return()}
  pd <- matrix(pd, ncol = 2, byrow = TRUE)
}

# Examples
outlier.problem <- function(n = 150){
  library('tdatools')
  X <- sample.annulus(n)
  X_outlier <- rbind(X, matrix(runif(6, -10, 10), nrow = 3))
  layout(matrix(1:4,ncol=2))
  plot(X, xlim = c(-30, 30), ylim = c(-30, 30), xlab = "", ylab = "", pch = 19, cex.axis = 2)
  plot(X_outlier, xlim = c(-30, 30), ylim = c(-30, 30), xlab = "", ylab = "", pch = 19, cex.axis = 2)
  pd <- diagram(X, "point-cloud", 1)
  pd_outlier <- diagram(X_outlier, "point-cloud", 1)
  plot(pd$pairs[[2]],  xlim = c(0, 40), ylim = c(0, 40), xlab = "", ylab = "", pch = 19, cex.axis = 2, cex = 1.5)
  lines(c(0, 40), c(0, 40))
  plot(pd_outlier$pairs[[2]],  xlim = c(0, 40), ylim = c(0, 40), xlab = "", ylab = "", pch = 19, cex.axis = 2, cex = 1.5)
  lines(c(0, 40), c(0, 40))
}

first.line.graph <- function(num_repeats = 10000){
  layout(matrix(1:3, ncol = 3), width = c(1, 1, 1))
  
  # Plot data
  edges <- matrix(c(1:6, 2:7), ncol = 2)
  vertices <- c(10, 11, 12.5, 13, 9.9, 20, 1)
  plot(NULL, xlim = c(1, 7), ylim = c(0, 20), xlab = "", ylab = "", cex.axis = 2)
  lines(1:7, vertices, lwd = 2)
  
  # Plot pd
  pd0 <- pers_hom_filtered_graph(edges, vertices)
  # Append global min/max pair
  pd0 <- rbind(pd0, range(vertices))
  plot(pd0,  xlim = c(0, 20), ylim = c(0, 20), xlab = "", ylab = "", pch = 19, cex.axis = 2, cex = 1.5)
  lines(c(0, 20), c(0, 20))

  sd <- seq(0.01, 0.1, length.out = 100)
  g5 <- numeric(100)
  g1 <- numeric(100)

  for(t in 1:100){
    h5 <- 0
    h1 <- 0
    for(i in 1:num_repeats){
      vertices_perturb <- vertices + rnorm(7, sd = sd[t])
      pd0_perturb <- pers_hom_filtered_graph(edges, vertices_perturb)
      pd0_perturb <- rbind(pd0_perturb, range(vertices_perturb))
      if(any(pd0_perturb[, 1] == vertices_perturb[1])){
        birth_at_1 <- which(pd0_perturb[, 1] == vertices_perturb[1])
        h1 <- h1 + pd0_perturb[birth_at_1, 2] - pd0_perturb[birth_at_1, 1]
      }
      if(any(pd0_perturb[, 1] == vertices_perturb[5])){
        birth_at_5 <- which(pd0_perturb[, 1] == vertices_perturb[5])
        h5 <- h5 + pd0_perturb[birth_at_5, 2] - pd0_perturb[birth_at_5, 1]
      }
    }
    g5[t] <- h5/num_repeats
    g1[t] <- h1/num_repeats
  }
  plot(sd, g5, type = 'l', col = 'red', xlab = "Alpha", ylab = "Approximation of g", xlim = c(0.01, 0.1), ylim = c(0, 14), cex.axis = 2, cex.lab = 1.6, lwd = 2)
  lines(sd, g1, col = 'blue', lwd = 2)
  lines(sd, g1 + g5, col = 'orange', lwd = 2)
  return(list(g5 = g5, g1 = g1))
}

second.line.graph <- function(num_repeats = 10000){
  layout(matrix(1:3, ncol = 3))
  
  # Plot data
  edges <- matrix(c(1:4, 2:5), ncol = 2)
  vertices <- c(5, 1.1, 1, 1.05, 15)
  plot(NULL, xlim = c(1, 5), ylim = c(0, 15), xlab = "", ylab = "", cex.axis = 2)
  lines(1:5, vertices, lwd = 2)
  
  # Plot pd
  pd0 <- pers_hom_filtered_graph(edges, vertices)
  # Append global min/max pair
  pd0 <- rbind(pd0, range(vertices))
  plot(pd0,  xlim = c(0, 15), ylim = c(0, 15), xlab = "", ylab = "", pch = 19, cex.axis = 2, cex = 1.5)
  lines(c(0, 20), c(0, 15))
  
  sd <- seq(0.01, 0.5, length.out = 100)
  g <- matrix(0, nrow = 5, ncol = 100)
  
  for(t in 1:100){
    h <- numeric(5)
    for(i in 1:num_repeats){
      vertices_perturb <- vertices + rnorm(5, sd = sd[t])
      pd0_perturb <- pers_hom_filtered_graph(edges, vertices_perturb)
      pd0_perturb <- rbind(pd0_perturb, range(vertices_perturb))
      for(b in 1:dim(pd0_perturb)[1]){
        birth_vertex <- which(vertices_perturb == pd0_perturb[b, 1])
        h[birth_vertex] <- h[birth_vertex] + pd0_perturb[b, 2] - pd0_perturb[b, 1]
      }
    }
    g[, t] <- h/num_repeats
  }
  
  plot(sd, g[2,], type = 'l', col = 'red', xlab = "Alpha", ylab = "Approximation of g", xlim = range(sd), ylim = c(0, 15), cex.axis = 2, cex.lab = 1.6, lwd = 2)
  lines(sd, g[4,], col = 'orange', lwd = 2)
  lines(sd, g[3,], col = 'blue', lwd = 2)
  lines(sd, g[2,] + g[4,] + g[3,], col = 'purple', lwd = 2)
  return(g)
}

distance.to.curve <- function(num_repeats = 10000){
  library('tdatools')
  layout(matrix(1:3, ncol = 3))
  # Original data
  ptcloud <- matrix(c(0, 1, 2, 7, 12, 7, 2, 1, 0, 0.1, 1, 0.12, 5, 0, -5, -0.12, -1, -0.1),
                    ncol = 2)
  plot(NULL, xlim = c(0, 12), ylim = c(-5, 5), xlab = "", ylab = "", cex.axis = 2)
  lines(ptcloud, lwd = 2)
  
  # Original pd
  D <- as.matrix(dist(ptcloud))
  for(i in 1:8){
    D[i+1, i] <- 0
  }
  pd <- diagram(D, "distance-matrix", 1)
  pd1 <- pd$pairs[[2]]
  plot(pd1,  xlim = c(0, 12), ylim = c(0, 12), xlab = "", ylab = "", pch = 19, cex.axis = 2, cex = 1.5)
  lines(c(0, 12), c(0, 12))
  
  sd <- seq(0.01, 0.1, length.out = 100)
  g <- array(0, c(9, 9, 100))
  pb <- txtProgressBar(0, 100, style = 3)
  for(t in 1:100){
    h <- matrix(0, nrow = 9, ncol = 9)
    for(i in 1:num_repeats){
      ptcloud_perturb <- ptcloud + matrix(rnorm(18, sd = sd[t]), ncol = 2)
      D <- as.matrix(dist(ptcloud_perturb))
      for(i in 1:8){
        D[i+1, i] <- 0
      }
      D <- round(D, 5)
      pd <- diagram(D, "distance-matrix", 1)
      pd1 <- pd$pairs[[2]]
      pd1 <- round(pd1, 5)
      for(b in dim(pd1)[1]){
        if(!any(D == pd1[b, 1])){
          print(D)
          print(pd1[b, 1])
          stop("Didn't find a generating edge")
          }
        h[D == pd1[b, 1]] <- h[D == pd1[b, 1]] + pd1[b, 2] - pd1[b, 1]
      }
    }
    g[,, t] <- h/num_repeats
    setTxtProgressBar(pb, t)
  }
  plot(sd, g[1,9,], type = 'l', col = 'red', xlab = "Alpha", ylab = "Approximation of g", xlim = range(sd), ylim = c(0, 10), cex.axis = 2, cex.lab = 1.6, lwd = 2)
  lines(sd, g[3,7,], col = 'orange', lwd = 2)
  return(g)
}
