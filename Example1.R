# Run the following to see a demo
demo <- function(){
  dependencies()
  
  grid <- choose.grid(-75, 75, -75, 200, 1)
  X <- sample.double.annulus(250)
  double.annulus.output <- double.annulus.experiment(X, grid$xs, grid$ys)
  cols <- colorRampPalette(c("white", "red"))(100)
  maxcol <- max(max(double.annulus.output$grid1),
                max(double.annulus.output$grid2))
  cols_range <- seq(0, maxcol, length.out = 100)
  dev.new()
  double.annulus.plot(X,
                      grid$xs,
                      grid$ys,
                      double.annulus.output$grid1,
                      double.annulus.output$grid2,
                      cols,
                      cols_range,
                      maxcol)
  
  X <- sample.annulus(250)
  grid <- choose.grid(-200, 200, -200, 200, 1)
  grid_bw0.2 <- bandwidth.experiment(X, grid$xs, grid$ys, bw = 0.2)
  grid_bw1 <- bandwidth.experiment(X, grid$xs, grid$ys, bw = 1)
  grid_bw3 <- bandwidth.experiment(X, grid$xs, grid$ys, bw = 3)
  grid_bw10 <- bandwidth.experiment(X, grid$xs, grid$ys, bw = 10)
  dev.new()
  bandwidth.plot(X,
                 grid$xs,
                 grid$ys,
                 grid_bw0.2,
                 grid_bw1,
                 grid_bw3,
                 grid_bw10,
                 cols,
                 cols_range)
}

# Run this to load required packages
dependencies <- function(){
  library("TDA")
  library("Matrix")
}

# Sampling functions
sample.annulus <- function(n, R1 = 20){
  ptcloud <- matrix(0, ncol = 2, nrow = n)
  for(t in 1:n){
    while(T){
      temp <- cbind(runif(1, -50, 50), runif(1, -50, 50))
      d1 <- sqrt(temp[1]^2 + temp[2]^2)
      if(d1 <= 50 & d1 >= R1){
        ptcloud[t,] <- temp
        break
      }
    }
  }
  return(ptcloud)
}
sample.double.annulus <- function(n, R1 = 20, R2 = 40){
  ptcloud <- matrix(0, ncol = 2, nrow = n)
  for(t in 1:n){
    while(T){
      temp <- cbind(runif(1, -50, 50), runif(1, -50, 140))
      d1 <- sqrt(temp[1]^2 + temp[2]^2)
      if(d1 <= 50 & d1 >= R1){
        ptcloud[t,] <- temp
        break
      }
      d2 <- sqrt(temp[1]^2 + (temp[2]- 90)^2)
      if(d2 <= 50 & d2 >= R2){
        ptcloud[t,] <- temp
        break
      }
    }
  }
  return(ptcloud)
}

# Auxiliary functions
do.segments.intersect <- function(a1, a2, b1, b2){
  # a1, a2, b1, b2 are endpoints of two segments
  x00 <- a1[1]
  y00 <- a1[2]
  x10 <- b1[1]
  y10 <- b1[2]
  x01 <- a2[1]-a1[1]
  y01 <- a2[2]-a1[2]
  x11 <- b2[1]-b1[1]
  y11 <- b2[2]-b1[2]
  
  d <- x11*y01 - x01*y11
  if(d == 0){return(F)}
  
  s = (1/d)*((x00 - x10)*y01 - (y00 - y10)* x01)
  t = (1/d)*-(-(x00 - x10)*y11 + (y00 - y10)*x11)
  return(0 <= s & s <= 1 & 0 <= t & t <= 1)
}
does.segment.intersect.square <- function(x, y, v1, v2) {
  if (x[1] <= v1[1] &
      v1[1] <= x[2] &
      y[1] <= v1[2] &
      v1[2] <= y[2]) {
    return(T)
  }
  if (x[1] <= v2[1] &
      v2[1] <= x[2] &
      y[1] <= v2[2] &
      v2[2] <= y[2]) {
    return(T)
  }
  return(
    do.segments.intersect(c(x[1], y[1]), c(x[2], y[1]), v1, v2) |
      do.segments.intersect(c(x[1], y[1]), c(x[1], y[2]), v1, v2) |
      do.segments.intersect(c(x[2], y[1]), c(x[2], y[2]), v1, v2) |
      do.segments.intersect(c(x[1], y[2]), c(x[2], y[2]), v1, v2)
  )
}
find.segment.grid.intersection <- function(xs, ys, grid, v1, v2){
  # grid is a sparse matrix with nrow = length(xs) - 1, ncol = length(ys) - 1
  
  # Find bounding box
  for(i in 1:length(xs)){
    xstart <- i
    if(xs[i + 1] > min(v1[1], v2[1])){
      break
    }
  }
  for(i in length(xs):1){
    xend <- i
    if(xs[i-1] < max(v1[1], v2[1])){
      break
    }
  }
  for(i in 1:length(ys)){
    ystart <- i
    if(ys[i + 1] > min(v1[2], v2[2])){
      break
    }
  }
  for(i in length(ys):1){
    yend <- i
    if(ys[i-1] < max(v1[2], v2[2])){
      break
    }
  }
  
  # Check each square in bounding box
  for(s in xstart:(xend - 1)){
    for(t in ystart:(yend - 1)){
      if(grid[s,t] | does.segment.intersect.square(c(xs[s], xs[s+1]),
                                                   c(ys[t], ys[t+1]),
                                                   v1, v2)){
        grid[s,t] <- 1
      }
    }
  }
  return(grid)
}
choose.grid <- function(x_low, x_high, y_low, y_high, grid_width){
  xs <- seq(x_low, x_high, by = grid_width)
  ys <- seq(y_low, y_high, by = grid_width)
  return(list(xs = xs, ys = ys))
}

# Main functions
double.annulus.experiment <- function(X, xs, ys, num_repeats = 50, bw = 3){
  grid1 <- Matrix(0, nrow = length(xs) - 1, ncol = length(ys) - 1, sparse = T)
  grid2 <- Matrix(0, nrow = length(xs) - 1, ncol = length(ys) - 1, sparse = T)
  n <- dim(X)[1]
  
  pb <- txtProgressBar(0, num_repeats, style = 3)
  for(t in 1:num_repeats){
    # Perturb data
    temp_X <- X + matrix(rnorm(2*n, mean = 0, sd = bw), ncol = 2)
    
    # persistence diagram of alpha complex
    DiagAlphaCmplx <- alphaComplexDiag(
      X = temp_X, library = c("GUDHI", "Dionysus"), location = TRUE)
    one <- which(DiagAlphaCmplx[["diagram"]][, 1] == 1)
    lifetimes <- order(DiagAlphaCmplx[["diagram"]][one, 3] - DiagAlphaCmplx[["diagram"]][one, 2], decreasing = T)
    first <- one[lifetimes[1]]
    second <- one[lifetimes[2]]
    
    temp_grid1 <- Matrix(0, nrow = length(xs) - 1, ncol = length(ys) - 1, sparse = T)
    for (j in seq_len(dim(DiagAlphaCmplx[["cycleLocation"]][[first]])[1])) {
      v1 <- DiagAlphaCmplx[["cycleLocation"]][[first]][j, 1, ]
      v2 <- DiagAlphaCmplx[["cycleLocation"]][[first]][j, 2, ]
      temp_grid1 <- find.segment.grid.intersection(xs, ys, temp_grid1, v1, v2)
    }
    grid1 <- grid1 + temp_grid1
    
    temp_grid2 <- Matrix(0, nrow = length(xs) - 1, ncol = length(ys) - 1, sparse = T)
    for (j in seq_len(dim(DiagAlphaCmplx[["cycleLocation"]][[second]])[1])) {
      v1 <- DiagAlphaCmplx[["cycleLocation"]][[second]][j, 1, ]
      v2 <- DiagAlphaCmplx[["cycleLocation"]][[second]][j, 2, ]
      temp_grid2 <- find.segment.grid.intersection(xs, ys, temp_grid2, v1, v2)
    }
    grid2 <- grid2 + temp_grid2
    setTxtProgressBar(pb, t)
  }
  return(list(grid1 = grid1/num_repeats, grid2 = grid2/num_repeats))
}
bandwidth.experiment <- function(X, xs, ys, num_repeats = 50, bw){
  grid1 <- Matrix(0, nrow = length(xs) - 1, ncol = length(ys) - 1, sparse = T)
  n <- dim(X)[1]
  pb <- txtProgressBar(0, num_repeats, style = 3)
  for(t in 1:num_repeats){
    # Perturb data
    temp_X <- X + matrix(rnorm(2*n, mean = 0, sd = bw), ncol = 2)
    
    # persistence diagram of alpha complex
    DiagAlphaCmplx <- alphaComplexDiag(
      X = temp_X, library = c("GUDHI", "Dionysus"), location = TRUE)
    one <- which(DiagAlphaCmplx[["diagram"]][, 1] == 1)
    lifetimes <- order(DiagAlphaCmplx[["diagram"]][one, 3] - DiagAlphaCmplx[["diagram"]][one, 2], decreasing = T)
    first <- one[lifetimes[1]]
    
    temp_grid1 <- Matrix(0, nrow = length(xs) - 1, ncol = length(ys) - 1, sparse = T)
    for (j in seq_len(dim(DiagAlphaCmplx[["cycleLocation"]][[first]])[1])) {
      v1 <- DiagAlphaCmplx[["cycleLocation"]][[first]][j, 1, ]
      v2 <- DiagAlphaCmplx[["cycleLocation"]][[first]][j, 2, ]
      temp_grid1 <- find.segment.grid.intersection(xs, ys, temp_grid1, v1, v2)
    }
    grid1 <- grid1 + temp_grid1
    
    setTxtProgressBar(pb, t)
  }
  return(grid1/num_repeats)
}

# Plotting functions
double.annulus.plot <- function(X, xs, ys, grid1, grid2, cols, cols_range, maxcol){
  
  layout(matrix(1:4,ncol=4), width = c(2, 2, 2, 1),height = c(1,1,1,1))

  DiagAlphaCmplx <- alphaComplexDiag(
    X = X, library = c("GUDHI", "Dionysus"), location = TRUE,
    printProgress = TRUE)
  one <- which(DiagAlphaCmplx[["diagram"]][, 1] == 1)
  lifetimes <- order(DiagAlphaCmplx[["diagram"]][one, 3] - DiagAlphaCmplx[["diagram"]][one, 2], decreasing = T)
  first <- one[lifetimes[1]]
  second <- one[lifetimes[2]]
  
  # Plot dataset
  plot(X, xlim=c(-50, 50), ylim=c(-50, 150), ylab="", xlab="", xaxt = 'n', yaxt = 'n', pch = 19, cex = 0.35)
  
  # Draw first cycle
  for (j in seq_len(dim(DiagAlphaCmplx[["cycleLocation"]][[first]])[1])) {
    lines(
      DiagAlphaCmplx[["cycleLocation"]][[first]][j, , ], pch = 19, lwd = 3, col = "#228b22")
  }
  # Draw second cycle
  for (j in seq_len(dim(DiagAlphaCmplx[["cycleLocation"]][[second]])[1])) {
    lines(
      DiagAlphaCmplx[["cycleLocation"]][[second]][j, , ], pch = 19, lwd = 3, col = "blue")
  }
  
  # Plot first grid
  plot(NULL, xlim=c(-50, 50), ylim=c(-50, 150), ylab="", xlab="", xaxt = 'n', yaxt = 'n')
  for(i in 1:dim(grid1)[1]){
    for(j in 1:dim(grid1)[2]){
      col = cols[findInterval(grid1[i,j], cols_range, rightmost.closed = T)]
      rect(xs[i], ys[j], xs[i+1], ys[j+1], col = col, border = NA)
    }
  }
  points(X, pch = 19, cex = 0.35)
  
  # Plot second grid
  plot(NULL, xlim=c(-50, 50), ylim=c(-50, 150), ylab="", xlab="", xaxt = 'n', yaxt = 'n')
  #cols <- heat.colors(100, alpha = 0.5)
  #cols <- rgb(1, 0, 0, seq(0, 1, length.out = 100))
  for(i in 1:dim(grid2)[1]){
    for(j in 1:dim(grid2)[2]){
      col = cols[findInterval(grid2[i,j], cols_range, rightmost.closed = T)]
      rect(xs[i], ys[j], xs[i+1], ys[j+1], col = col, border = NA)
    }
  }
  points(X, pch = 19, cex = 0.35)
  
  # Plot legend
  legend_image <- as.raster(matrix(rev(cols), ncol=1))
  plot(c(0,2),c(0,2),type = 'n', axes = F,xlab = '', ylab = '', main = '')
  text(x=1.4, y = seq(0,1,l=5), labels = round(seq(0, maxcol, l=5), 2),  cex = 2.5)
  rasterImage(legend_image, -0.3, 0, 0.5,1)
}
bandwidth.plot <- function(X, xs, ys, grid1, grid2, grid3, grid4, cols, cols_range){
  plot.single.grid <- function(X, xs, ys, grid, cols, cols_range){
    # Plot second grid
    plot(NULL, xlim=c(-50, 50), ylim=c(-50, 50), ylab="", xlab="", xaxt = 'n', yaxt = 'n', asp = 1)
    for(i in 1:dim(grid)[1]){
      for(j in 1:dim(grid)[2]){
        if(xs[i] >= -50 & xs[i+1] <= 50 & ys[j] >= -50 & ys[j+1] <= 50){
          col = cols[findInterval(grid[i,j], cols_range, rightmost.closed = T)]
          rect(xs[i], ys[j], xs[i+1], ys[j+1], col = col, border = NA)
        }
      }
    }
    points(X, pch = 19, cex = 0.7)
  }
  
  layout(matrix(1:4,ncol=4), width = c(2,2, 2, 2),height = c(1,1,1,1))
  for(grid in c(grid1, grid2, grid3, grid4)){
    plot.single.grid(X, xs, ys, grid, cols, cols_range)
  }
}