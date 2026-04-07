df <- read_csv("/Users/clairedonnat/Documents/CCAR3_code/experiments/data/fMRI-data/activation_groups.xlsx")
df <- read_xlsx("/Users/clairedonnat/Documents/CCAR3_code/experiments/data/fMRI-data/data/")
positions = read_csv("/Users/clairedonnat/Documents/CCAR3_code/experiments/data/fMRI-data/parcellation_coordinates.csv")
real_labels <- readxl::read_xlsx("/Users/clairedonnat/Documents/CCAR3_code/experiments/data/fMRI-data/BNA_subregions.xlsx") %>% 
  filter(`Label ID.L`<211)
real_labels["name"] = sapply(real_labels$...6, function(x){str_split(x, ",")[[1]][1]})
real_labels["full_label"] = sapply(real_labels$...6, function(x){str_split(x, ",")[[1]][2]})

positions["name"] = sapply(positions$Region, function(x){return(str_split(x, "_")[[1]][2])})
positions  = positions %>% left_join(real_labels %>% dplyr::select(`...6`,
                                                                   name,
                                                                   `Gyrus`),
                                     by = c("name"))


positions = positions %>%
  dplyr::select(-c("...6"))
positions$Gyrus[1:19] = positions$name[1:19] 
positions$Gyrus[15] = "Brain Stem" 
positions$name[15] = "Brain Stem" 
positions$Gyrus[which(is.na(positions$Gyrus))] = "STG, Superior Temporal Gyrus"

write_csv(positions, "~/Documents/CCAR3_code/experiments/data/fMRI-data/groups.csv")
unique(positions$Gyrus)
distanceMatrix <- as.matrix(dist(positions[, c("x", "y", "z")]))

# Find the 6 nearest neighbors for each point
numNeighbors <- 4
edges <- c()
for (i in 1:nrow(distanceMatrix)){
  x = distanceMatrix[i,] 
  edges = rbind(edges,
                t(rbind(rep(i, numNeighbors), order(x)[1:numNeighbors+1])))
}
edges = data.frame(edges)
colnames(edges) <- c("Source", "Target")
# Create a graph from these edges
g <- graph_from_edgelist(as.matrix(edges), directed = FALSE)

# Optionally, plot the graph
#plot(g)
# Find the connected components
comp <- components(g)
# Get the number of connected components
num_components <- comp$no
print(num_components)
Gamma <- get_edge_incidence(g, weight = 1)
Gamma_dagger = pinv(Gamma)


edges2 <- c()
for (i in 1:nrow(positions)){
  neighbours= which(positions$Gyrus == positions$Gyrus[i])
  edges2 = rbind(edges2,
                 t(rbind(rep(i, length(neighbours)), neighbours)))
}
edges2 = data.frame(edges2)
colnames(edges2) <- c("Source", "Target")
edges2 = edges2 %>%
  filter(Source < Target)
# Create a graph from these edges
g2 <- graph_from_edgelist(as.matrix(edges2), directed = FALSE)

# Optionally, plot the graph
#plot(g2)
# Find the connected components
comp <- components(g2)
# Get the number of connected components
num_components <- comp$no
print(num_components)
Gamma2 <- get_edge_incidence(g2, weight = 1)
Gamma_dagger2 = pinv(Gamma2)

unique_regions = unique(positions$Gyrus)
groups <- sapply(1:length(unique_regions),
                 function(i){which(positions$Gyrus == unique_regions[i])}
)
###