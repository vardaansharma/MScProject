library(igraph)
library(CINNA )
PPI_Network <- V(gg)$name
degrees <- igraph::degree(gg)
betweenness <- betweenness(gg, directed = FALSE)
CC <- transitivity(gg, type='local', vids = V(gg))
cl <- calculate_centralities(gg, include = 'Semi Local Centrality')
topD <- which.max(degrees)
dists <- distances(gg, v=topD)
Pr <- page_rank(gg, directed = FALSE)$vector
PPI_Network_rds <- data.frame(PPI_Network, degrees, betweenness, CC, cl, as.vector(dists), Pr)
saveRDS(PPI_Network_rds, file = sprintf("%s/originalNtwrkMeasures.rds",OUT[2]))