library(readr)
library(igraph)
library(ggplot2)
library(poweRlaw)

# source('/afs/inf.ed.ac.uk/user/s18/s1882216/MSCPROJECT/MScProject/PPI network/setUp.R');
# all_proteins <- read_csv("C:/Users/homoe/OneDrive - University of Edinburgh/PROJECT/MScProject/Networks/Combined.csv")
# PPI_List_Unique <- read_csv("C:/Users/homoe/OneDrive - University of Edinburgh/PROJECT/MScProject/Networks/PPI_List_Unique.csv")
# PPI_List <- read_csv("C:/Users/homoe/OneDrive - University of Edinburgh/PROJECT/MScProject/Networks/PPI_List.csv")

all_proteins <- read_csv("/afs/inf.ed.ac.uk/user/s18/s1882216/MSCPROJECT/MScProject/Networks/Combined.csv")
PPI_List_Unique <- read_csv("/afs/inf.ed.ac.uk/user/s18/s1882216/MSCPROJECT/MScProject/Networks/PPI_List_Unique.csv")
PPI_List <- read_csv("/afs/inf.ed.ac.uk/user/s18/s1882216/MSCPROJECT/MScProject/Networks/PPI_List.csv")


graph <- graph_from_data_frame(PPI_List,directed=FALSE)

plot(graph,vertex.size=5, vertex.label=NA, edge.arrow.mode=0,layout=layout_nicely)

#simplifying the graph ( removeing duplicates etc)/
agg <- function(.x)toString(unique(.x))
sgraph <- simplify(graph,remove.multiple=TRUE,remove.loops=TRUE,edge.attr.comb=agg)

#checking the components (subgraphs)
comps <- components(sgraph)

#taking only biggest subgraph
i <- which.max(comps$csize)
vg <- groups(comps)
csgraph <- induced_subgraph(sgraph,vg[[i]])
write_graph(csgraph, "PPI_Network.csv","ncol")
c(vcount(csgraph), max(comps$csize)) #making sure we got the biggest one by looking at vertice count


#graph metrics

#Degree
degrees <- igraph::degree(csgraph)
max_degree <- which.max(degrees)
write.table(degrees,'C:/Users/homoe/OneDrive - University of Edinburgh/PROJECT/MScProject/Supplimentary/degrees.csv',sep = ',')
#V(csgraph)[max_degree] incident_edges(csgraph,max_degree)

#Scale free network
deg_dist <- degree_distribution(csgraph)
xx <- 1:length(deg_dist) -1
ind <- which(deg_dist > 0)
deg_dist <- deg_dist[ind]
xx <- xx[ind]
qplot(xx, deg_dist, log='xy', xlab='degree', ylab='distribution/fraction of nodes')
ggsave('C:/Users/homoe/OneDrive - University of Edinburgh/PROJECT/MScProject/supplimentary/scale_free_graph.pdf')

fit1 <- fit_power_law(degrees)

#checking the power law 
degree_displ = displ$new(degrees)
est = estimate_xmin(degree_displ)
checking = bootstrap_p(degree_displ)
#checking$ p #gives the p value that degrees is a power law distrubution
degree_displ$setXmin(est)
bs = bootstrap(degree_displ, no_of_sims = 5000, threads = 2) # try more threads for faster processing
alpha = mean(bs$bootstraps[,3])
alpha_sd = sd(bs$bootstraps[,3])


#path length
distances(csgraph, v=max_degree) -> path_dists
qplot(path_dists[1,], binwidth=0.25,xlab = 'path distances',ylab='number of proteins')
ggsave('C:/Users/homoe/OneDrive - University of Edinburgh/PROJECT/MScProject/supplimentary/path_distances.pdf')
#table( path_dists )  
summary(path_dists[1])

#centrality / closeness
closns <- closeness(csgraph)
max_closns <- which.max(closns)
qplot(closns,xlab='closeness') 
summary(closns)
ggsave('C:/Users/homoe/OneDrive - University of Edinburgh/PROJECT/MScProject/supplimentary/closeness.pdf')

#betweenness
betwness <- betweenness(csgraph, directed = FALSE)
qplot(betwness, log = 'x',xlab='betweenness')
ggsave('C:/Users/homoe/OneDrive - University of Edinburgh/PROJECT/MScProject/supplimentary/betweenness.pdf')

#cluster coefficient / transivity
trans <- transitivity(csgraph, type='local',vids=V(csgraph))
qplot(trans,xlab='transitivity')

#cluster coefficient vs other stuff

#assortativity degree
assortativity_degree(csgraph,directed = FALSE)

#page rank
pagerank <- page_rank(csgraph,directed = FALSE)$vector
qplot(pagerank,log='x')


#clustering
