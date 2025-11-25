# Diffusion in Networks with R

# 0. Set up working directory

setwd("c:/Users/lengyel.balazs/Desktop/LengyelB_2022June/BCE/SDS/Applied Network Science")

# 0. Install and call in package
#install.packages("igraph") # In case you have not installed the package yet. 
#It is worth to repeat installation on a -say- yearly basis. Igraph R syntax is changing frequently.

#install.packages("sand")

library(igraph) # Call this package in case you work with network objects.
library(sand) # Sand is important for teaching purposes only: it contains built-in data.


set.seed(42) # set the seed for random number generation - this makes 
#data(ppi.CC) # call in data from the SAND package

# 1. The SIR model on toy networks

# Create toy networks
  #The parameters have been chosen so as to guarantee graphs of roughly the same
  #average degree (i.e., around 10), and a number of vertices of similar order of magnitude,
  #since both are basic characteristics expected to fundamentally affect the
  #progression of an epidemic in ways that do not reflect interesting differences in topology.
set.seed(42)
gl <- list() # we are going to store the networks in a list
gl$ba <- barabasi.game(250, m=5, directed=FALSE) # a scale-free network 
gl$er <- erdos.renyi.game(250, 1250, type=c("gnm")) # random network alternative: "gnp" p= probability
gl$ws <- watts.strogatz.game(1, 250, 12, 0.01) # small-world network


V(gl$ws)

png(filename="ba_net.png", width=1800, height=600, units="px")
par(mfrow=c(1,3))
plot(gl$ba, vertex.size=2)
plot(gl$er, vertex.size=2)
plot(gl$ws, vertex.size=2)
dev.off()
# Set the model parameters of the SIR model
b <- 0.5 # The rate of infection of an individual that is susceptible and has a single infected neighbor. 
            # The infection rate of a susceptible individual with n infected neighbors is n times beta.
ga <- 1  # The rate of recovery of an infected individual.

ntrials <- 100 # Number of simulations.

# Run the simulation
sim <- lapply(gl, sir, beta=b, gamma=ga,
   no.sim=ntrials) # We run the diffusion simulations 

sim_ba<-sir(gl$ba, beta=b, gamma=ga,
           no.sim=ntrials)
sim_er<-sir(gl$er, beta=b, gamma=ga,
             no.sim=ntrials)
sim_ws<-sir(gl$ws, beta=b, gamma=ga,
            no.sim=ntrials)

  # The output from each simulation is an sir object, containing information about
  # the times at which changes of states occurred and the values of the processes
  # NS(t),NI (t), and NR(t) at those times. 

# Present simulation results
png(filename="simulations.png", width=900, height=300, units="px")
par(mfrow=c(1,3))
plot(sim$er)
plot(sim$ba, color="palegoldenrod",
   median_color="blue", quantile_color="red")
plot(sim$ws, color="pink", median_color="red",
   quantile_color="red")
dev.off()
  # The results of plotting the total number of infectives NI(t) for each of these three networks is shown.

# Plot the median of the curves NI(t) together to see the differences

png(filename="sim_compare.png", width=600, height=600, units="px")
par(mfrow=c(1,1))
x.max <- max(sapply(sapply(sim, time_bins), max)) # set the maximum value on the x axis
y.max <- 1.05 * max(sapply(sapply(sim, function(x)
   median(x)[["NI"]]), max, na.rm=TRUE)) # set the maximum value on the y axis
plot(time_bins(sim$er), median(sim$er)[["NI"]],
   type="l", lwd=2, col="blue", xlim=c(0, x.max),
   ylim=c(0, y.max), xlab="Time",
   ylab=expression(N[I](t)))
lines(time_bins(sim$ba), median(sim$ba)[["NI"]],
   lwd=2, col="gold")
lines(time_bins(sim$ws), median(sim$ws)[["NI"]],
   lwd=2, col="red")
legend("topright", c("ER", "BA", "WS"),
   col=c("blue", "gold", "red"), lty=1)
dev.off()


# 2. Simulate complex diffusion manually - a Toy iWiW network of 1.000 nodes

# 2.1. Read Data
ver <- read.csv("vertices_sample.csv")
edg <- read.csv("edges_sample.csv")

# 2.2. Create the graph
network_fs <- graph_from_data_frame(edg, directed = FALSE, vertices = ver)

# 2.3. Initialize nodes and parameters
nodes_fs <- V(network_fs)
diffusers_fs <- which(ver$month == 3 )
susceptibles_fs <- setdiff(nodes_fs, diffusers_fs)

time_fs <- numeric(length(nodes_fs))
time_fs[diffusers_fs] <- 1

p_fs <- 0.000104
q_fs <- 0.12

# 2.4. Precompute neighborhood information

library(tictoc)

tic("Precompute neighborhoods")
neigh <- ego(network_fs, 1, nodes = nodes_fs)
toc()

# 2.5. Simulate diffusion
set.seed(1042)
time_adoption <- 2

tic("Diffusion process")
while (time_adoption < 129) {
  tic(paste("Time step", time_adoption))
  
  # Compute the fraction of neighbors that are diffusers
  a <- vapply(neigh, function(x) sum(x %in% diffusers_fs) / (length(x) - 1), numeric(1))
  
  # Determine new adopters
  prob <- p_fs + a * q_fs
  adopters <- susceptibles_fs[runif(length(susceptibles_fs)) < prob[susceptibles_fs]]
  
  # Update states
  time_fs[adopters] <- time_adoption
  diffusers_fs <- c(diffusers_fs, adopters)
  susceptibles_fs <- setdiff(nodes_fs, diffusers_fs)
  
  # Update neighborhoods
  neigh <- neigh[setdiff(nodes_fs, diffusers_fs)]
  
  cat(sum(time_fs != 0), "-->", time_adoption, "\n")
  time_adoption <- time_adoption + 1
  toc()
}
toc()

# 2.6. Plot Results
results <- data.frame(
  Time = 1:max(time_fs),
  CDF = cumsum(tabulate(time_fs)) / vcount(network_fs)
)

plot(results$Time, results$CDF, type = "b",
     ylab = "CDF", xlab = "Time")



# Homework: 
# H1. Simulate the virus diffusion on the flight network with beta=0.5 and gamma=1. When do we expect the peak of infections?
# H2. Exclude top diffuser airports from the network. How does this modify the diffusion curve?
# H3. What is the best network centrality metric that should be applied to close down airports in order to slow down diffusion?

# Visualization of the network can be found at the bottom of the code


# The global airline network

# Read, clean the data and plot the network
air_nodes=read.table("airlines_nodes.csv", header=T, sep=";")
air_nodes$lng=as.numeric(air_nodes$lng)
air_nodes=air_nodes[!is.na(air_nodes$lng),]
air_edges=read.table("airlines_edges.csv", header=T, sep=";")

air_g=graph_from_data_frame(air_edges, vertices=air_nodes, directed=F) # network object from edgelist and nodelist
air_g=simplify(air_g) # remove self-loops and multiple edges

V(air_g)$lng # vertices of g
E(air_g) # edges of g

coordinates<-matrix(c(V(air_g)$lng, V(air_g)$lat),nrow=length(V(air_g)$lat),ncol=2) # coordinate matrix to locate nodes in geographical networks

# plot the network into a PNG file
png(filename="air_net.png", width=6000, height=2000, units="px")
plot.igraph(air_g,
            layout=coordinates, 
            # vertex.color=V(air_g)$col, 
            # vertex.frame.color=V(air_g)$fc, 
            vertex.size=0.5,
            edge.color="gray50", edge.curved=F,
            edge.width=1,
            axes=F, frame=F, rescale=T,
            vertex.label=NA
)
dev.off()

# Color network after community finding algorithm 
#install.packages("viridis")
library(viridis)

set.seed(1024)
air_c=cluster_louvain(air_g) # the Louvain algorithm is a very popular method that is fast on large graphs
length(air_c)
length(V(air_g))
length(air_c[sizes(air_c)>10])

c=as.data.frame(sizes(air_c)) # create community table
names(c)=c("comm", "size")

v=data.frame(as.integer(V(air_g)$name), as.vector(membership(air_c))) # membership table
names(v)=c("Id", "comm")

col=c[c$size>10,] # colors into community and membership tables
#col$col=magma(length(col$comm))
col$col=c(magma(length(air_c[sizes(air_c)>10])-11),
          "olivedrab1", "navajowhite3", "khaki2", "darkblue",
          "orange","cyan", "deeppink", "red",
          "green4", "brown", "dodgerblue")
col=col[,c("comm", "col")]
c=merge(c,col, by="comm", all.x=T)
v=merge(v,c, by="comm", all.x=T)

air_nodes_new=merge(air_nodes,v, by="Id", all.x=T) # colors into nodes table

# regenerate network and coordinates
air_g=graph_from_data_frame(air_edges, vertices=air_nodes_new, directed=F) # network object from edgelist and nodelist
coordinates<-matrix(c(V(air_g)$lng, V(air_g)$lat),nrow=length(V(air_g)$lat),ncol=2) # coordinate matrix to locate nodes in geographical networks


# plot the network with the chosen colors of communities
png(filename="air_comm.png", width=6000, height=2000, units="px")
plot(air_g,
     layout=coordinates, 
     vertex.color=V(air_g)$col,
     vertex.frame.color=V(air_g)$col, 
     vertex.size=1,
     edge.color="grey",
     edge.width=0.1,edge.curved=F,
     axes=F, frame=F, rescale=T,
     vertex.label=NA
)
dev.off()

# Evaluate the community finding algorithm

set.seed(24)
c1=cluster_louvain(air_g, resolution = 2)
c2=cluster_louvain(air_g, resolution = 2)

compare(membership(c1), membership(c2), method="nmi")