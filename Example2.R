# Function to computed 0-pd using union find
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
  pd <- matrix(pd, ncol = 2, byrow = TRUE)
}
get_birth_pt <- function(bar_of_interest, pd, vertices){
  pd <- pd[order(pd[, 2] - pd[, 1], decreasing = T),]
  birth_vertex <- which(vertices[, 3] == pd[bar_of_interest, 1])
  return(vertices[birth_vertex, 1:3])
}

main <- function(bar_of_interest = 28, 
                 dist_to_birth_point = 5, 
                 num_repeats = 1000,
                 sd = 0.1){
  # Load dataset
  # load("~/R/BullittTreesSubject1.Rdata")
  
  # Get persistence diagram
  pd <- pers_hom_filtered_graph(edges, vertices[, 3])
  birth_pt <- get_birth_pt(bar_of_interest, pd, vertices[, 1:3])
  
  D <- numeric(num_repeats)
  pb <- txtProgressBar(0, num_repeats, style = 3)
  for(i in 1:num_repeats){
    vertices_perturb <- vertices[, 1:3] + 
      matrix(rnorm(dim(vertices)[1]*3, sd = sd), ncol = 3)
    pd_perturb <- pers_hom_filtered_graph(edges, vertices_perturb[, 3])
    birth_pt_perturb <- get_birth_pt(bar_of_interest,
                                     pd_perturb,
                                     vertices_perturb)
    D[i] <- dist(rbind(birth_pt, birth_pt_perturb))
    setTxtProgressBar(pb, i)
  }
  print(sum(D <= dist_to_birth_point) / num_repeats)
  
  library(ggplot2)
  df <- data.frame(dist = D)
  p <- ggplot(df, aes(dist)) + 
    stat_ecdf(geom = "step") + 
    labs(title="Empirical cumulatative distribution function", x="Distance to birth point", y="Percentage of points") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  # Plot dataset with birth pt and neighborhood
  x <- vertices[,1]
  y <- vertices[,2]
  z <- vertices[,3]
  library(plotrix)
  E <- dim(edges)[1]
  edge_segments <- matrix(0,E,8) # (x1,y1,z1,r1,x2,y2,z2,r2)
  theta <- pi/24 # amount to rotate brain image around the z axis
  for (i in 1:E)
    edge_segments[i,] <- c(t(vertices[edges[i,],]))
  plot(c(min(x*cos(theta)+y*sin(theta)),max(x*cos(theta)+y*sin(theta))), c(min(z),max(z)), type="n", xlab="", ylab="", axes=FALSE, asp=1)
  segments(edge_segments[,1]*cos(theta)+edge_segments[,2]*sin(theta),edge_segments[,3],edge_segments[,5]*cos(theta)+edge_segments[,6]*sin(theta),edge_segments[,7],col="blue",lwd=(edge_segments[,4]+edge_segments[,8]))
  points(birth_pt[1]*cos(theta)+birth_pt[2]*sin(theta), birth_pt[3], type="p", pch=20, col="black", lwd=4)
  draw.circle(birth_pt[1]*cos(theta)+birth_pt[2]*sin(theta), birth_pt[3],radius=dist_to_birth_point,border="black",lwd=4)
  
  return(p)
}