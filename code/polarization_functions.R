library(igraph)
library(Matrix) #Regularized Spectral Clustering
library(dplyr)
require(clustAnalytics)
require(xtable)
library(greed)
#library(NetIndices)
library(randnet)
library(modMax)
library(ggraph)
library(netseg)
library(ks)
library(philentropy)

## Random Walk Controvery (RWC)
## Adaptive Random Walk Controversy (ARWC)

get_influencer_nodes <- function(G, left_nodes, right_nodes, n_influencers) {
  # Validate inputs
  if(length(left_nodes) == 0 | length(right_nodes) == 0 | n_influencers < 0) {
    stop("Invalid input: nodes lists must be non-empty and n_influencers must be non-negative.")
  }
  
  left_degrees <- degree(G, v=left_nodes)
  right_degrees <- degree(G, v=right_nodes)
  
  # Adjusting n_influencers if it's a proportion
  k_left <- k_right <- n_influencers
  if(n_influencers < 1) {
    k_left <- max(1, as.integer(n_influencers * length(left_nodes)))
    k_right <- max(1, as.integer(n_influencers * length(right_nodes)))
  }
  
  # Sort and select top influencers, keeping their degrees
  left_influencers <- sort(left_degrees, decreasing = TRUE)[1:k_left]
  right_influencers <- sort(right_degrees, decreasing = TRUE)[1:k_right]

  return(list(left_influencers, right_influencers))
}

rwc_perform_random_walk <- function(G, left_influencers, right_influencers, starting_node) {
  # Determine the side of the starting node
  starting_side <- if(starting_node %in% left_influencers) "left" else "right"
  
  current_node <- starting_node
  
  while(TRUE) {
    neighbors <- neighbors(G, current_node)
    
    if(length(neighbors) == 0) {
      stop("The current node has no neighbors to walk to.")
    }
    
    next_node <- sample(V(G)$name[neighbors], 1)
    
    if(next_node %in% left_influencers) {
      return(list(starting_side, "left"))
    } else if(next_node %in% right_influencers) {
      return(list(starting_side, "right"))
    } else {
      current_node <- next_node
    }
  }
}

RWC_polarization <- function(G, ms, n_influencers, n_sim, n_walks) {
  left_nodes <- names(ms[ms == 1])
  right_nodes <- names(ms[ms == 2])
  
  influencers <- get_influencer_nodes(G, left_nodes, right_nodes, n_influencers)
  left_influencers <- names(influencers[[1]])
  right_influencers <- names(influencers[[2]])
  
  rwc_dist <- numeric(n_sim)
  
  
  for (i in 1:n_sim) {
    counts <- integer(4)  # To store counts for left_left, right_left, left_right, right_right
    names(counts) <- c("left_left", "right_left", "left_right", "right_right")
    
    # Inside the RWC_polarization function
    for (j in 1:n_walks) {
      starting_side <- sample(c("left", "right"), 1)
      which_random_starting_node <- sample(if(starting_side == "left") left_nodes else right_nodes, 1)
      result <- rwc_perform_random_walk(G, left_influencers, right_influencers, which_random_starting_node)
      ending_side <- result[[2]]
      
      counts[paste(starting_side, ending_side, sep = "_")] <- counts[paste(starting_side, ending_side, sep = "_")] + 1
    }
    
    
    # Calculating RWC for each simulation
    e1 <- counts["left_left"] / max(1, counts["left_left"] + counts["right_left"])
    e2 <- counts["right_left"] / max(1, counts["left_left"] + counts["right_left"])
    e3 <- counts["left_right"] / max(1, counts["right_right"] + counts["left_right"])
    e4 <- counts["right_right"] / max(1, counts["right_right"] + counts["left_right"])
    
    rwc <- e1 * e4 - e2 * e3
    rwc_dist[i] <- rwc
  }
  
  rwc_ave <- mean(rwc_dist)
  return(rwc_ave)
}


## Boundary Polarization
# Function to calculate boundary polarization
calculate_boundary_polarization <- function(graph, community_membership) {
  boundary_nodes <- which(sapply(1:vcount(graph), function(node) {
    # Check if the node has neighbors in different communities
    any(sapply(neighbors(graph, node), function(neighbor) {
      community_membership[node] != community_membership[neighbor]
    }))
  }))
  
  polarization_scores <- sapply(boundary_nodes, function(node) {
    intra_community_edges <- length(which(sapply(neighbors(graph, node), function(neighbor) {
      community_membership[node] == community_membership[neighbor]
    })))
    inter_community_edges <- length(which(sapply(neighbors(graph, node), function(neighbor) {
      community_membership[node] != community_membership[neighbor]
    })))
    
    # Calculate the ratio of inter-community to total edges
    if (intra_community_edges + inter_community_edges > 0) {
      inter_community_edges / (intra_community_edges + inter_community_edges)
    } else {
      0
    }
  })
  
  # Calculate the mean polarization across all boundary nodes
  mean_polarization <- mean(polarization_scores)
  return(mean_polarization)
}


## Adaptive EI Index (AEI)
AEI_polarization <- function(G, ms) {
  # Create two blocks of nodes based on the community assignments
  block_a <- names(ms[ms == 1])
  block_b <- names(ms[ms == 2])
  
  n_a <- length(block_a)
  n_b <- length(block_b)
  
  # Count internal connections within each block
  c_a <- sum(sapply(E(G)[from(block_a) & to(block_a)], function(e) is.element(ends(G, e)[2], block_a)))
  c_b <- sum(sapply(E(G)[from(block_b) & to(block_b)], function(e) is.element(ends(G, e)[2], block_b)))
  
  # Count external connections between the blocks
  c_ab <- sum(sapply(E(G), function(e) {
    ends_e <- ends(G, e)
    (is.element(ends_e[1], block_a) && is.element(ends_e[2], block_b)) ||
      (is.element(ends_e[1], block_b) && is.element(ends_e[2], block_a))
  }))
  
  # Density calculations with checks to prevent division by zero
  B_aa <- if (n_a > 1) c_a / (n_a * (n_a - 1) / 2) else 0
  B_bb <- if (n_b > 1) c_b / (n_b * (n_b - 1) / 2) else 0
  B_ab <- if (n_a > 0 && n_b > 0) c_ab / (n_a * n_b) else 0
  
  # B_ba is equivalent to B_ab, so it's not calculated separately
  B_ba <- B_ab
  
  # Calculate the AEI score
  aei_score <- (B_aa + B_bb - B_ab - B_ba) / (B_aa + B_bb + B_ab + B_ba)
  
  return(aei_score)
}


