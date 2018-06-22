# Example from Section 1.2.1 of 
# Stabilizing the unstable output of persistent homology computations
# by Paul Bendich, Peter Bubenik, and Alexander Wagner

# Copyright 2018, Peter Bubenik, All rights reserved.

main_directory <- "~/R/brain_arteries" # sample directory for Mac/linux users
setwd(main_directory) # change to the main directory

load("BullittTrees.RData")

bar_of_interest <- 28
t <- 1  # subject number from 1 to 97
dist_to_birth_point <- 5 # h equals the proportion of birth points within this distance of the observed birth point
theta <- pi/24 # amount to rotate brain image around the z axis

num_repeats <- 1
sd <- 0 # standard deviation for Gaussian noise
source("persistence_script.R")
pd_list[[1]][26:30,]
diffs <- abs(pd_list[[1]][26:30,5] - pd_list[[1]][28,5])
rank_diffs <- sort(diffs)[c(2,3)]
c(pd_list[[1]][28,5],rank_diffs,rank_diffs[2]/rank_diffs[1])
birth_point <- pd_list[[1]][28,6:8]
closest_point <- pd_list[[1]][25+which(rank(diffs)==2),6:8]
second_closest_point <- pd_list[[1]][25+which(rank(diffs)==3),6:8]
dist_closest <- dist(rbind(birth_point,closest_point))
dist_second_closest <- dist(rbind(birth_point,second_closest_point))
c(dist_closest,dist_second_closest)

num_repeats <- 100
sd <- 0.1 # standard deviation for Gaussian noise
source("persistence_script.R")
D <- numeric(num_repeats)
for (i in 1:num_repeats)
  D[i] <- dist(rbind(bars_of_interest[i,6:8],birth_point))
hist(D,breaks=200,main="Histogram",xlab="Distance to birth point")
P <- ecdf(D)
plot(P, xlab="Distance to birth point", ylab="Percentage of points",main="Empirical cumulatative distribution function")
h <- sum(D <= dist_to_birth_point) / num_repeats
h

library(ggplot2)
df <- data.frame(dist = D)
ggplot(df, aes(dist)) + stat_ecdf(geom = "step") + labs(title="Empirical cumulatative distribution function",x="Distance to birth point", y="Percentage of points") + theme(plot.title = element_text(hjust = 0.5))

vertices <- vertices_list[[t]]
edges <- edges_list[[t]]

x <- vertices[,1]
y <- vertices[,2]
z <- vertices[,3]

library(plotrix)
edge_segments <- matrix(0,E,8) # (x1,y1,z1,r1,x2,y2,z2,r2)
for (i in 1:E)
  edge_segments[i,] <- c(t(vertices[edges[i,],]))
plot(c(min(x*cos(theta)+y*sin(theta)),max(x*cos(theta)+y*sin(theta))), c(min(z),max(z)), type="n", xlab="", ylab="", axes=FALSE, asp=1)
segments(edge_segments[,1]*cos(theta)+edge_segments[,2]*sin(theta),edge_segments[,3],edge_segments[,5]*cos(theta)+edge_segments[,6]*sin(theta),edge_segments[,7],col="blue",lwd=(edge_segments[,4]+edge_segments[,8]))
#points(bars_of_interest[,6]*cos(theta)+bars_of_interest[,7]*sin(theta), bars_of_interest[,8], type="p", pch=1, col="red", lwd=4)
points(birth_point[1]*cos(theta)+birth_point[2]*sin(theta), birth_point[3], type="p", pch=20, col="black", lwd=4)
draw.circle(birth_point[1]*cos(theta)+birth_point[2]*sin(theta), birth_point[3],radius=5,border="black",lwd=4)

