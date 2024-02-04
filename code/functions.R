
  #Libraries
library(igraph)
library(Matrix) #Regularized Spectral Clustering
library(randnet) #Regularized Spectral Clustering
library(clustAnalytics) #library to compute weighted_clustering_coefficient
#library(fpc) # Evaluation and comparison of different clusterings
#library(dplyr)
#require(xtable)
library(ggraph) #nice plots for modularity optimization


# ###################### SOCIAL NETWORK CONSTRUCTION ######################




################################### PREPROCESSING ################################### 


#1. Giant Component
#The "giant component" of a network refers to the largest connected component within that network. In the context of social network analysis, this usually means the largest subgraph where most nodes are interconnected

preprocess_graph<- function(G){
    #get the connected components
    largest_component_subgraph <- largest_component(G)
    
    # cat("Number of nodes in G:", vcount(G), "\n")
    # cat("Number of edges in G:", ecount(G), "\n")
    # 
    # cat("Number of nodes in the largest connected component:", vcount(largest_component_subgraph), "\n")
    # cat("Number of edges in the largest connected component:", ecount(largest_component_subgraph), "\n")
    # 
    #2. Simplify graphs
    simple_G <- simplify(largest_component_subgraph, remove.multiple = TRUE, remove.loops = TRUE)
    
    
    
    ggraph(simple_G, layout = "fr") +  geom_edge_link() +
      theme_void() +
      labs(title = "Simplified Graph",
           color = "Clusters") +
      theme(legend.position = "topright")

    return(simple_G)

}





# ###################### GRAPH PARTITIONING ######################


################################### CONSENSUS CLUSTERING ################################### 

# Define a function to compute the consensus matrix
consensus_matrix <- function(clusterings) {
  #
  # Get the number of nodes and clusterings
  n <- length(clusterings[[1]])
  k <- length(clusterings)
  # Initialize an empty matrix
  mat <- matrix(0, nrow = n, ncol = n)
  # Loop over all pairs of nodes
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Count how many times they are in the same cluster
      count <- 0
      for (l in 1:k) {
        if (clusterings[[l]][i] == clusterings[[l]][j]) {
          count <- count + 1
        }
      }
      # Compute the proportion and store it in the matrix
      mat[i, j] <- count / k
      mat[j, i] <- mat[i, j]
    }
  }
  # Return the matrix
  return(mat)
}



redirect_subgraph <- function(original_graph, sub_graph) {
  # Get the edge list of the subgraph
  sub_g_edges <- as_edgelist(sub_graph)
  original_g_edges <- as_edgelist(original_graph)
  
  # Find indices of edges in the subgraph that exist in the original graph
  edge_indices <- unique(match(as.character(sub_g_edges), as.character(original_g_edges)))
  # I apply unique since the result is a duplicated set of edges, since each undirected edge is counted twice
  
  # Filter out NA indices (edges not found in the original graph)
  valid_indices <- !is.na(edge_indices)
  
  # Extract the valid edges from the subgraph
  valid_edges <- sub_g_edges[valid_indices, ]
  
  # Initialize an empty graph to store the directed subgraph
  directed_sub_g <- make_empty_graph(n = vcount(sub_graph), directed = TRUE)
  
  # Set the vertex names to be the same as the subgraph
  V(directed_sub_g)$name <- V(sub_graph)$name
  V(directed_sub_g)$cluster <- V(sub_graph)$cluster
  
  # Add the valid edges to the directed subgraph
  directed_sub_g <- add_edges(directed_sub_g, valid_edges)
  
  return(directed_sub_g) 
}



###################################  MODULARITY OPTIMIZATION CLUSTERING ###################################
# Function to compute and plot clusterings given the input  data

#Tested also Label Propagation, Walktrap, Edge Betweeness, Spinglass, Infomap but performed very badly
plot_clusterings <- function(graph, k_values) {
  set.seed(123) 
  
  ######LOUVAIN########
  # Louvain algorithm
  louvain_clust <- multilevel.community(graph, resolution = 0.1 ) #the lower the resolution the bigger the clusters
  memb_louvain <- membership(louvain_clust)
  
  #plot
  g<-ggraph(graph, layout = "fr") +
    geom_edge_link() +
    geom_node_point(aes(color = factor(memb_louvain))) +
    theme_void() +
    labs(title = "Louvain Spectral Clustrering",
         color = "Clusters") +
    theme(legend.position = "topright")
  print(g)
  
  
  #######LEIDEN########
  #Leiden algorithm (Modularity Optimization)
  #Algorithm to maximize modularity is Louvain. 
  #A more recently developed community detection algorithm which improves on the Louvain algorithm is the Leiden algorithm.
  leiden_clust <- cluster_leiden( graph,
                                  objective_function = "modularity",
                                  weights = E(graph)$weight,
                                  resolution_parameter = 0.2,
                                  n_iterations = 100
  )
  memb_leiden <- membership(leiden_clust)
  
  g<-ggraph(graph, layout = "fr") +
    geom_edge_link() +
    geom_node_point(aes(color = factor(memb_leiden))) +
    theme_void() +
    labs(title = paste0("Leiden Clustering"), color = "Clusters") +
    theme(legend.position = "topright")
  
  print(g)
  
  
  
  #######RSC#######
  adj_matrix <- as_adjacency_matrix(graph, sparse = FALSE)
  for (k in k_values) {
    rsc <- reg.SP(adj_matrix, K = k, tau = 10, lap = FALSE, nstart = 30, iter.max = 100)
    
    g<-ggraph(graph, layout = "fr") +
      geom_edge_link() +
      geom_node_point(aes(color = factor(rsc$cluster))) +
      theme_void() +
      labs(title = paste("Regularized Spectral Clustering (k =", k, ")"),
           color = "Clusters") +
      theme(legend.position = "topright")
    
    print(g)
    
  }

  
  
}




# CHOOSE BEST CLUSTERING AMONG DIFFERENT ONES
# function to determine the best-ranked clustering among a set of clusterings for a dataset with no ground truth
choose_best_algo <- function(G, k_values) {
  best_scores <- list()
  best_clusterings <- list()
  adj_matrix <- as_adjacency_matrix(G, sparse = FALSE)
  
  for (k in k_values) {
    # RSC
    rsc <- reg.SP(adj_matrix, K = k, tau = 10, lap = FALSE, nstart = 30, iter.max = 100)
    # Calculate metrics for RSC
    rsc_metrics <- c(modularity = modularity(G, rsc$cluster))
    
    best_scores[[paste0("rsc_", k)]] <- rsc_metrics["modularity"]
    
    best_clusterings[[paste0("rsc_", k)]]<-rsc$cluster
  }
  
  #Louvain
  clust_louvain <- membership(multilevel.community(G, resolution = 0.1 )) #the lower the resolution the bigger the clusters
  
  louvain_metrics <- c(modularity = modularity(G, clust_louvain))
  
  best_scores[["Louvain"]] <- louvain_metrics["modularity"]
  
  best_clusterings[["Louvain"]] <- clust_louvain
  
  
  #Leiden
  clust_leiden <- membership(cluster_leiden(G, objective_function = "modularity", weights = E(G)$weight,  
                                            resolution_parameter = 0.2,  n_iterations = 100))
  
  leiden_metrics <- c(modularity = modularity(G, clust_leiden))
  
  best_scores[["Leiden"]] <- leiden_metrics["modularity"] 
  
  
  best_clusterings[["Leiden"]]<-clust_leiden
  
  # Find the algorithm with the maximum score
  best_algorithm <- names(which.max(unlist(best_scores)))
  best_clustering<- best_clusterings$best_algorithm
  
  # Print the scores in a table
  best_scores_df <- as.data.frame(do.call(rbind, best_scores))
  colnames(best_scores_df) <- "Maximized Total Score"
  print(best_scores_df)
  
  return(list(name = best_algorithm, clusters = best_clusterings))
}



###################### CREATE RANDOMIZED MODELS #########################
#d=0 fix the number of links
generate_ER <- function(graph){
  set.seed(123)
  er_graph <- sample_gnm(n = vcount(graph), m = ecount(graph), directed = TRUE)
  V(er_graph)$name <- V(graph)$name
  return(er_graph)
}


#d=1 fix the degree sequence
generate_config_model <- function(graph){
  #The “simple” method connects the out-stubs of the edges (undirected graphs) or the out-stubs and in-stubs (directed graphs) together. 
  #This way loop edges and also multiple edges may be generated. This method is not adequate if one needs to generate simple graphs with a
  #given degree sequence. The multiple and loop edges can be deleted, but then the degree sequence is distorted and there is nothing
  #to ensure that the graphs are sampled uniformly.
  
  
  #set.seed(123) 
  configuration_model <- sample_degseq(degree(graph, mode = "out"), degree(graph, mode = "in"), method = "simple")
  
  conf_graph <- graph_from_edgelist(as_edgelist(configuration_model), directed = TRUE)
  V(conf_graph)$name <- V(graph)$name
         
  return(conf_graph)
  
}





###################### COMPUTE POLARIZATION SCORES #########################

#Logistic regression

compute_accuracy <- function(labels, scores){

  #set.seed(123) #for reproducibility
  
  scores$label <- labels

  # Drop the 'net_name' from the dataframe
  scores$net_name <- NULL
  #scores$smi <- NULL
  
  #replace NA with column average
  scores[sapply(scores, is.numeric)] <- lapply(scores[sapply(scores, is.numeric)], function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
  
  #compute_roc(labels, scores)
    
    
  train_indices <- sample(1:length(labels), 0.8 * length(labels))  # 80% training data
  # Split the data into training and testing sets
  train_data <- scores[train_indices, ]
  test_data <- scores[-train_indices, ]

  print(train_data)
  
  # Train a logistic regression model
  model <- glm(label ~ ., data = train_data, family = binomial)


  print("--------Model summary--------")
  print(summary(model))
  print("----------------")
  
  threshold <- 0.5
  
  print("Training accuracy: ")
  predictions <- predict(model, newdata = train_data, type = "response")
  predicted_labels <- ifelse(predictions > threshold, 1, 0)
  accuracy <- sum(predicted_labels == train_data$label) / length(train_data$label)
  print(accuracy)
  
   
  print("Test accuracy:")
  # Make predictions on the test set
  predictions <- predict(model, newdata = test_data, type = "response")
  
  # Assign labels based on  threshold
  predicted_labels <- ifelse(predictions > threshold, 1, 0)
  # Evaluate model accuracy
  accuracy <- sum(predicted_labels == test_data$label) / length(test_data$label)
  print(accuracy)
  
  
  return(list(model = model, accuracy = accuracy))
}




compute_ensemble_accuracy<-function(labels, scores){
  
  scores$ei <- -scores$ei
  scores$net_name<-NULL
  scores$label <- labels
  
  
  #replace NA with column average
  scores[sapply(scores, is.numeric)] <- lapply(scores[sapply(scores, is.numeric)], function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
  
  
  #compute_roc(labels,scores)
  
  
  
  train_indices <- sample(1:length(labels), 0.8 * length(labels))  # 80% training data
  # Split the data into training and testing sets
  train_data <- scores[train_indices, ]
  test_data <- scores[-train_indices, ]
  
  
  
  # Fit a decision tree model
  tree_model <- rpart(label ~ ., data = train_data, method = "class")
  

  print("--------DEcision tree training--------")
  

  print("--------Model summary--------")
  print(summary(tree_model))
  print("----------------")
  rpart.plot(tree_model)
  

  print("Training accuracy: ")
  predictions <- predict(tree_model, newdata = train_data, type = "class")
  accuracy <- sum(predictions == train_data$label) / length(train_data$label)
  print(accuracy)
  
  
  print("Test accuracy:")
  # Make predictions on the test set
  predictions <- predict(tree_model, test_data, type = "class")
  accuracy <- sum(predictions == test_data$label) / length(test_data$label)
  print(accuracy)
  
  
  return(list(model = tree_model, accuracy = accuracy))
}


compute_roc<-function(labels,scores){
  fit1<-glm(label~ei, data = scores, family="binomial")
  fit2<-glm(label~modularity, data = scores, family="binomial")
  fit3<-glm(label~smi, data = scores, family="binomial")
  fit4<-glm(label~rwc, data = scores, family="binomial")
  fit5<-glm(label~bp, data = scores, family="binomial")
  fit6<-glm(label~aei, data = scores, family="binomial")
  
  preds=predict(fit1)
  roc1=roc(labels~preds)
  preds=predict(fit2)
  roc2=roc(labels~preds)
  preds=predict(fit3)
  roc3=roc(labels~preds)
  preds=predict(fit4)
  roc4=roc(labels~preds)
  preds=predict(fit5)
  roc5=roc(labels~preds)
  preds=predict(fit6)
  roc6=roc(labels~preds)
  # Create a list to store ROC curve data
  plot(roc1, col = "red", main = "ROC Curves for Scores")
  plot(roc2, add = TRUE, col = "blue")
  plot(roc3, add = TRUE, col = "green")
  plot(roc4, add = TRUE, col = "orange")
  plot(roc5, add = TRUE, col = "purple")
  plot(roc6, add = TRUE, col = "brown")
  legend("bottomright", 
        legend = c(paste0("ei (AUC = ",round(auc(roc1), 3),")"), paste0("modularity (AUC = ",round(auc(roc2), 3),")" ), paste0("smi (AUC = ",round(auc(roc3), 3),")"), paste0("rwc (AUC = ",round(auc(roc4), 3),")" ), paste0("bp (AUC = ",round(auc(roc5), 3),")" ), paste0("aei (AUC = ",round(auc(roc6), 3),")")), 
        col = c("red", "blue", "green", "orange", "purple", "brown"), lty = 1,
        text.col = "black",
        bty = "n") # Setting bty = "n" removes the legend box
  
}