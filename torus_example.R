# Example from Section 1.2.2 of 
# Stabilizing the unstable output of persistent homology computations
# by Paul Bendich, Peter Bubenik, and Alexander Wagner

# Copyright 2018, Alexander Wagner, All rights reserved.

dependencies <- function(){
  library("rgl")
  library("deldir")
  library("scatterplot3d")
  # library("perseus4r")
  library("viridisLite")
  library("ggplot2")
  library("plotly")
}

run.me <- function(N=100, M=10, a=0.2){
  dependencies()
  
  # Take sample
  pts <- mfld.sample(N)
  
  # Plot underlying function on square
  print(demo.function())
  
  # Show flat plots
  mesh <- make.mesh(pts, torus = F)
  
  plot.on.cube(N, mesh, plot.type = "full")
  par3d("windowRect"=c(0,0,800,800))
  
  open3d()
  plot.on.cube(N, mesh, plot.type = "sample")
  par3d("windowRect"=c(0,0,800,800))
  
  print(plot.on.square(N, mesh))
  
  # Show torus plots
  mesh <- make.mesh(pts)
  
  open3d()
  plot.on.torus(mesh, plot.type = "full", R=2)
  par3d("windowRect"=c(0,0,800,800))
  
  open3d()
  plot.on.torus(mesh, plot.type = "flat", R=2)
  par3d("windowRect"=c(0,0,800,800))
  
  open3d()
  plot.on.torus(mesh, plot.type = "sample", R=2)
  par3d("windowRect"=c(0,0,800,800))
  
  # # Compute 0-dim barcode of induced simplicial complex
  # ph0 <- compute.mesh.pd(mesh$p,mesh$e)
  # 
  # # Plot 0-dim persistence diagram
  # bd <- data.frame(b = mesh$p$z[ph0[,1]], d = mesh$p$z[ph0[,2]])
  # threshold <- mean(order(bd$d-bd$b, decreasing=TRUE)[1:2])
  # print(plot.persistence.diagram(bd, t = threshold))
  # 
  # # Approximate convolution
  # big.bars <- compute.conv(pts = pts, M = M, alpha = a, t = threshold)
  # plot.conv(big.bars = big.bars)
  # inquad2 <- get.quad(big.bars$bx, big.bars$by)==2
  # g_approx <- sum(big.bars$size[inquad2])/M
  # print(paste("Approximate value of convolution:", g_approx))
  # return()
}

demo.function <- function(){
  x <- seq(from=-pi,to=pi,length.out = 500)
  y <- seq(from = -pi,to = pi,length.out = 500)
  z <- matrix(0,nrow = 500,ncol = 500)
  for(i in 1:500){
    for(j in 1:500){
      if(x[i]<0 & y[j]<0){
        z[i,j]=.1*sin(x[i])*sin(y[j])
        next
      }
      z[i,j]=sin(x[i])*sin(y[j])
    }
  }
  scene = list(aspectratio = list(x = 1, y = 1, z = 1.25))
  p <- plot_ly(x=x, y=y, z=z, showscale = F, type = "surface", reversescale = T) %>%
    layout(scene = scene)
  return(p)
}

mfld.sample <- function(N){
  # Sample x,y-points
  x <- runif(N, min=-pi, max=pi)
  y <- runif(N, min=-pi, max=pi)
  # Apply heights
  z <- sin(x)*sin(y)
  # Shrink 3rd quadrant
  z <- z*(1-.9*(x<0 & y<0))
  pts <- data.frame(x,y,z)
  return(pts)
}

make.mesh <- function(pts, alpha=0, torus = T){
  # Add alpha-bandwidth noise
  N <- dim(pts)[1]
  pts <- pts + alpha*matrix(rnorm(3*N),N,3)
  # Put points back on torus
  pts[,c("x","y")] <- (pts[,c("x","y")]+pi)%%(2*pi)-pi
  # Sort by height
  pts <- pts[order(pts$z),]
  # Form buffer for Delaunay on torus
  buffer <- vector("list", 9)
  buffer[[1]] <- pts
  k <- 2
  for(i in -1:1){
    for(j in -1:1){
      if(i==0 & j==0){next}
      temp.pts <- pts
      temp.pts$x <- temp.pts$x + i*2*pi
      temp.pts$y <- temp.pts$y + j*2*pi
      buffer[[k]] <- temp.pts
      k <- k+1
    }
  }
  pts <- do.call(rbind, buffer)
  # Compute Delaunay triangulation
  deltri <- deldir(pts$x, pts$y, z=pts$z)
  # Store edges
  edges <- cbind(deltri$delsgs$ind1, deltri$delsgs$ind2)
  edges <- edges[(edges[,1] <= N) | (edges[,2] <= N),]
  if(torus){
    edges[edges>N] <- edges[edges>N]%%N
    pts <- pts[1:N,]}
  # Sort edges by appearance
  edges <- t(apply(edges,1,sort))
  edges <- edges[order(edges[,2]),]
  return(list("p"=pts,"e"=edges))
}

plot.on.square <- function(N, mesh){
  pts <- mesh$p[1:N,]
  min.pt <- pts[which.min(pts$z),]
  pts <- pts[-which.min(pts$z),]
  nbcol <- 200
  color <- rev(viridis(nbcol))
  color <- color[cut(pts$z, nbcol)]
  pts$label <- 1:length(pts$x)
  
  p <- ggplot(pts, aes(x, y)) + labs(x = "X", y = "Y")
  edges <- data.frame(x1 = mesh$p$x[mesh$e[,1]], x2 = mesh$p$x[mesh$e[,2]],
                      y1 = mesh$p$y[mesh$e[,1]], y2 = mesh$p$y[mesh$e[,2]])
  p <- p + geom_segment(aes(x=x1, y=y1, xend = x2, yend = y2), alpha = 0.2,  data = edges)
  p <- p + geom_point(aes(color = factor(label)), show.legend = F) + scale_colour_manual(values = color)
  p <- p + geom_point(data = min.pt, col = "Red", size = 2.5) 
  t <- .1
  p <- p + coord_cartesian(xlim = c(-pi+t,pi-t), ylim = c(-pi+t,pi-t))
  p <- p + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.title.y=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())
  return(p)
}

plot.on.cube <- function(N, mesh, plot.type, quad.col = F){
  if(!plot.type %in% c("sample", "full")){
    stop("plot.type not sample or full")
  }
  pts.in <- mesh$p[1:N,]
  gmin <- which.min(pts.in[,3])
  nbcol <- 200
  color <- rev(viridis(nbcol))
  zcol <- cut(pts.in[-gmin,3], nbcol)
  color <- color[zcol]
  # Plotting
  t <- 0
  plot3d(pts.in[-gmin,], size = 6, col = color, alpha = 0.8,
         axes = F, xlab = "", ylab = "", zlab = "")
  points3d(x=pts.in[gmin,1], y=pts.in[gmin,2], z=pts.in[gmin,3],
           size = 8, col = "Red", add = T)
  if(plot.type == "full"){
    inside.idx <- (mesh$e[,1] <= N & mesh$e[,2] <= N)
    edges.in <- mesh$p[as.vector(t(mesh$e[inside.idx,])),]
    segments3d(edges.in, alpha = 0.3, add = T)
    edges.out <- clip.edges(edges = mesh$e[!inside.idx,], pts = mesh$p)
    segments3d(edges.out, alpha = 0.3, add = T)
    border = matrix(c(pi,-pi,-pi,-pi,-pi,pi,pi,pi,
                      pi,pi,pi,-pi,-pi,-pi,-pi,pi,
                      0,0,0,0,0,0,0,0), ncol = 3)
    segments3d(border, add = T)
  }
}

plot.on.torus <- function(mesh, R = 2, plot.type){
  if(!plot.type %in% c("sample", "flat", "full")){
    stop("plot.type not sample, flat, or full")
  }
  r <- 1.1
  theta <- mesh$p$x
  phi <- mesh$p$y
  # Form flat embedding
  pts <- cbind((R+r*cos(theta))*cos(phi),
               (R+r*cos(theta))*sin(phi),
               r*sin(theta))
  if(plot.type == "full"){
    # Compute tangent vectors at each sample pt
    bx <- -sin(phi)
    by <- cos(phi)
    sx <- -cos(phi)*sin(theta)
    sy <- -sin(phi)*sin(theta)
    sz <- cos(theta)
    # Compute normal vectors and normalize
    n <- cbind(by*sz, -bx*sz, bx*sy - by*sx)
    n <- n/sqrt(n[,1]*n[,1] + n[,2]*n[,2] + n[,3]*n[,3])
    # Scale by sample heights
    n <- n*mesh$p$z
    pts <- pts + n}
  
  gmin <- which.min(mesh$p$z)
  nbcol <- 200
  color <- rev(viridis(nbcol))
  zcol <- cut(mesh$p$z[-gmin], nbcol)
  color <- color[zcol]
  
  # Plotting
  plot3d(pts[-gmin,], size = 6, col = color, alpha = 0.8,
         axes = F, xlab = "", ylab = "", zlab = "")
  aspect3d(x=(R+r), y=(R+r), z=1)
  points3d(x=pts[gmin,1], y=pts[gmin,2], z=pts[gmin,3],
           size = 8, col = "Red", add = T)
  if(plot.type != "sample"){
    edges <- pts[as.vector(t(mesh$e)),]
    segments3d(edges, alpha = 0.2, add = T)}
}

compute.mesh.pd <- function(pts,edges){
  #Create input for Perseus
  N <- dim(pts)[1]
  vtces <- cbind(rep(0,N),1:N,1:N)
  edges <- cbind(rep(1,dim(edges)[1]),edges,edges[,2])
  cplx <- list(vtces,edges)
  # Apply perseus
  barcode <- perseus.nmfsimtop(n=1,cplx)
  ph0 <- barcode[barcode[,1]==0,2:3]
  # Replace infinite bars with global max a la extended persistence
  ph0[ph0 == -1] <- N
  return(ph0)
}

plot.persistence.diagram <- function(bd,t){
  # bd = Nx2 matrix of birth/death times
  # t = threshold
  
  # Determine big bars
  barsize <- bd[,2]-bd[,1]
  barsize[barsize<=t] <- 1
  barsize[barsize>t] <- 2
  bd$barsize <- barsize
  # Plot pd
  p <- ggplot(bd,aes(b,d)) + geom_point(aes(color=factor(barsize)),size=3) + scale_color_manual(values = c("Black", "Red"))
  p <- p + geom_abline(intercept=0,slope=1) + geom_abline(intercept = 1.5,slope = 1,linetype = "dotted",color = "Black")
  p <- p + labs(title = "Persistence Diagram")
  p <- p + xlim(-1.5,.5) + ylim(-.5,1.5) + xlab("Birth") + ylab("Death") + theme(legend.position="none")
  return(p)
}

compute.conv <- function(pts, M=1000, alpha=0.2, t){
  bar.list <- vector("list", M)
  pb <- txtProgressBar(min = 0, max = M, style = 3)
  for(trial in 1:M){
    mesh <- make.mesh(pts = pts, alpha = alpha)
    ph0 <- compute.mesh.pd(mesh$p,mesh$e)
    ph0 <- data.frame(bx = mesh$p$x[ph0[,1]],
                      by = mesh$p$y[ph0[,1]],
                      size = mesh$p$z[ph0[,2]] - mesh$p$z[ph0[,1]])
    ph0 <- ph0[ph0$size>t,]
    bar.list[[trial]] <- ph0
    setTxtProgressBar(pb,trial)
  }
  big.bars <- do.call(rbind,bar.list)
  return(big.bars)
}

plot.conv <- function(big.bars){
  quad <- get.quad(big.bars$bx,big.bars$by)
  quad[quad!=2] <- "Blue"
  quad[quad==2] <- "Red"
  quad <- adjustcolor(quad,alpha.f=0.7)
  scatterplot3d(big.bars, color = quad, pch = 19, type = "h",
                xlab = "X,Y-plane", ylab = "", zlab = "Size",
                xlim = c(-3,3), ylim = c(-3,3), box = F)
}

clip.edges <- function(edges, pts){
  num.edges <- dim(edges)[1]
  edge.list <- vector("list", num.edges)
  for(i in 1:num.edges){
    new.pair <- pts[edges[i,],]
    x1 <- new.pair[1,1]
    x2 <- new.pair[2,1]
    y1 <- new.pair[1,2]
    y2 <- new.pair[2,2]
    z1 <- new.pair[1,3]
    z2 <- new.pair[2,3]
    #Bottom
    t <- (y1+pi)/(y1-y2)
    xnew <- (1-t)*x1+t*x2
    if(xnew >= -pi & xnew <= pi & t>=0 & t<=1){
      new.pair[2,1] <- xnew
      new.pair[2,2] <- (1-t)*y1+t*y2
      new.pair[2,3] <- 0 #(1-t)*z1+t*z2
      edge.list[[i]] <- new.pair
      next}
    #Left
    t <- (x1+pi)/(x1-x2)
    ynew <- (1-t)*y1+t*y2
    if(ynew >= -pi & ynew <= pi & t>=0 & t<=1){
      new.pair[2,1] <- (1-t)*x1+t*x2
      new.pair[2,2] <- ynew
      new.pair[2,3] <- 0 #(1-t)*z1+t*z2
      edge.list[[i]] <- new.pair
      next}
    #Top
    t <- (y1-pi)/(y1-y2)
    xnew <- (1-t)*x1+t*x2
    if(xnew >= -pi & xnew <= pi & t>=0 & t<=1){
      new.pair[2,1] <- xnew
      new.pair[2,2] <- (1-t)*y1+t*y2
      new.pair[2,3] <- 0 #(1-t)*z1+t*z2
      edge.list[[i]] <- new.pair
      next}
    #Right
    t <- (x1-pi)/(x1-x2)
    ynew <- (1-t)*y1+t*y2
    if(ynew >= -pi & ynew <= pi & t>=0 & t<=1){
      new.pair[2,1] <- (1-t)*x1+t*x2
      new.pair[2,2] <- ynew
      new.pair[2,3] <- 0 #(1-t)*z1+t*z2
      edge.list[[i]] <- new.pair
      next}
  }
  edge.list <- do.call(rbind,edge.list)
  return(edge.list)
}

get.quad <- function(x,y){
  x <- (x>0)
  y <- (y>0)
  return (-2*x*y+x-y+3)
}