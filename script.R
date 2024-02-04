
  #Libraries
library(rpart)#decision tree
library(rpart.plot)

library(pROC) # to compute ROC
library(reshape2)
library(beanplot) #to create beanplot
library(igraph)
library(Matrix) #Regularized Spectral Clustering
library(randnet) #Regularized Spectral Clustering
library(clustAnalytics) #library to compute weighted_clustering_coefficient
#library(fpc) # Evaluation and comparison of different clusterings
#library(dplyr)
#require(xtable)
library(ggraph) #nice plots for modularity optimization

library(greed)
library(randnet)
library(modMax)
library(netseg)
library(ks)
library(philentropy)

library(openxlsx) #to export excel files




source("functions.R")
source("polarization_functions.R")



set.seed(123)


###################### LOAD FILES ######################

#load file with labels of all networks
network_info <- read.csv("network_data/network_info.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

#print the structure of the dataframe
str(network_info)

#df for labels (polarized or not)
network_labels <- data.frame(net_name=character(), label = numeric(), stringsAsFactors = FALSE)

#df for polarization scores
polarization_scores <- data.frame(net_name = character(), ei = numeric(), modularity = numeric(), smi = numeric(), rwc = numeric(),  bp = numeric(), aei = numeric(), stringsAsFactors = FALSE)


#net_name = "translaki_p1.edgelist" #network_info[1,1], [1,2..]
#edgelist_file = paste0("network_data/training_set/", net_name) 


#set variables for random walk -> n_influencers must be integer >= 2
n_influencers <-10
adaptive_n_influencers <- 0.01
n_sim <- 25 #50
n_walks <- 200#300

threshold <- 0.5


#loop over the rows of network_info
for (i in 1:nrow(network_info)) {
  #get the file name from the first column
  net_name <- network_info[i, 1]
  
  
  print(paste("Processing network",i,":", net_name))
  edgelist_file <- paste0("network_data/training_data/", net_name) 
  
  
  ####VERSION 2 OF DATA ###########
  # Read the data from the .edgelist file
  data <- read.table(edgelist_file, header = TRUE, sep = ",", col.names = c("retweeter", "retweeted", "weight"))
  
  # Convert the dataframe  to  numeric
  data <- as.data.frame(lapply(data, as.numeric))
  
  # Create a weighted directed graph from the data (# undirected graph to be used for algorithms that do not support directed)
  G <- graph_from_data_frame(data, directed = TRUE)#FALSE)
  
  #summary of the graph
  summary(G)
  
  if(vcount(G)>5000){
    print("Network size >5000. Discarded.")
  }else{
  
    
    
    #get network label
    label <- as.numeric(network_info[i, 3])
  
    # Append the network name and label to the network_labels dataframe
    #original
    network_labels <- rbind(network_labels, data.frame(net_name = net_name, label = label, stringsAsFactors = FALSE))
    
    #randomized ER model (label=0)
    network_labels <- rbind(network_labels, data.frame(net_name = paste0(net_name,"_ER"), label = 0, stringsAsFactors = FALSE))
    
    #randomized configuration model (label=0)
    network_labels <- rbind(network_labels, data.frame(net_name = paste0(net_name,"_conf1"), label = 0, stringsAsFactors = FALSE))
    network_labels <- rbind(network_labels, data.frame(net_name = paste0(net_name,"_conf2"), label = 0, stringsAsFactors = FALSE))
    #network_labels <- rbind(network_labels, data.frame(net_name = paste0(net_name,"_conf3"), label = 0, stringsAsFactors = FALSE))
    
    
    
  
    
    #plot graph
    # plot(
    #   G,
    #   layout = layout_with_fr(G), 
    #   vertex.size = 10,
    #   edge.arrow.size = 0.5,
    #   edge.width = E(G)$weight,
    #   main = "Weighted Undirected Graph",
    #   vertex.label = NA  
    # )
  
  
    ################################ PREPROCESSING ################################### 
  
    preprocessed_g <- preprocess_graph(G)
    #CREATE RANDOMIZED MODELS
    #d=0
    erdos_renyi_model <- generate_ER(preprocessed_g)
    preprocessed_er_g <-preprocess_graph(erdos_renyi_model)
    #d=1
    configuration_model1 <- generate_config_model(preprocessed_g)
    preprocessed_conf_g1 <-preprocess_graph(configuration_model1)
    configuration_model2 <- generate_config_model(preprocessed_g)
    preprocessed_conf_g2 <-preprocess_graph(configuration_model2)
    #configuration_model3 <- generate_config_model(preprocessed_g)
    #preprocessed_conf_g3 <-preprocess_graph(configuration_model3)
  
  
    ################################### GRAPH PARTITIONING################################### 
    #CITE: https://arxiv.org/abs/1703.09307
    
    clust_fluid <- cluster_fluid_communities(preprocessed_g, no.of.communities=2)
    clust_fluid_er <- cluster_fluid_communities(preprocessed_er_g, no.of.communities=2)
    clust_fluid_conf1 <- cluster_fluid_communities(preprocessed_conf_g1, no.of.communities=2)
    clust_fluid_conf2 <- cluster_fluid_communities(preprocessed_conf_g2, no.of.communities=2)
    #clust_fluid_conf3 <- cluster_fluid_communities(preprocessed_conf_g3, no.of.communities=2)
    
    ################## POLARIZATION SCORES #################
    
    ms <- membership(clust_fluid)
    ms_er <- membership(clust_fluid_er)
    ms_conf1 <- membership(clust_fluid_conf1)
    ms_conf2 <- membership(clust_fluid_conf2)
    #ms_conf3 <- membership(clust_fluid_conf3)
    
    
  
    
    print("----------------------")
    print("Start computing scores....")
    
     ### Adaptive EI Index (AEI)
     aei_score <- AEI_polarization(preprocessed_g, ms)
     aei_score_er <- AEI_polarization(preprocessed_er_g, ms_er)
     aei_score_conf1 <- AEI_polarization(preprocessed_conf_g1, ms_conf1)
     aei_score_conf2 <- AEI_polarization(preprocessed_conf_g2, ms_conf2)
     #aei_score_conf3 <- AEI_polarization(preprocessed_conf_g3, ms_conf3)
     
     
     ### Segregation Matrix Index (SMI) -> needs directed graph
     # Assuming your communities are labeled with numbers and you choose communities 1 and 2
    
     V(preprocessed_g)$community <- ms
     V(preprocessed_er_g)$community <- ms_er
     V(preprocessed_conf_g1)$community <- ms_conf1
     V(preprocessed_conf_g2)$community <- ms_conf2
     #V(preprocessed_conf_g3)$community <- ms_conf3
     
     V(preprocessed_g)$selected_community <- ifelse(ms == 1, 1, ifelse(ms == 2, 2, NA))  # Apply the smi function to the directed graph using the new attribute
     smi_values <- smi(preprocessed_g,"community", normalize = TRUE)
     cluster_sizes<-table(clust_fluid$membership)
     smi_values <- sum(smi_values * cluster_sizes) / sum(cluster_sizes)
     
     
     V(preprocessed_er_g)$selected_community <- ifelse(ms_er == 1, 1, ifelse(ms_er == 2, 2, NA))  # Apply the smi function to the directed graph using the new attribute
     smi_values_er<- smi(preprocessed_er_g,"community", normalize = TRUE)
     cluster_sizes<-table(clust_fluid_er$membership)
     smi_values_er <- sum(smi_values_er * cluster_sizes) / sum(cluster_sizes)
     # 
     # 
     V(preprocessed_conf_g1)$selected_community <- ifelse(ms_conf1 == 1, 1, ifelse(ms_conf1 == 2, 2, NA))  # Apply the smi function to the directed graph using the new attribute
     smi_values_conf1 <- smi(preprocessed_conf_g1,"community", normalize = TRUE)
     cluster_sizes<-table(clust_fluid_conf1$membership)
     smi_values_conf1 <- sum(smi_values_conf1 * cluster_sizes) / sum(cluster_sizes)
     
     
     V(preprocessed_conf_g2)$selected_community <- ifelse(ms_conf2 == 1, 1, ifelse(ms_conf2 == 2, 2, NA))  # Apply the smi function to the directed graph using the new attribute
     smi_values_conf2 <- smi(preprocessed_conf_g2,"community", normalize = TRUE)
     cluster_sizes<-table(clust_fluid_conf2$membership)
     smi_values_conf2 <- sum(smi_values_conf2 * cluster_sizes) / sum(cluster_sizes)
     
     # V(preprocessed_conf_g3)$selected_community <- ifelse(ms_conf3 == 1, 1, ifelse(ms_conf2 == 2, 2, NA))  # Apply the smi function to the directed graph using the new attribute
     # smi_values_conf3 <- smi(preprocessed_conf_g3,"community", normalize = TRUE)
     # cluster_sizes<-table(clust_fluid_conf3$membership)
     # smi_values_conf3 <- sum(smi_values_conf3 * cluster_sizes) / sum(cluster_sizes)
     # 

     
    
    print("Compute scores for original graph")
     
    # Append the polarization scores to the dataframe
    polarization_scores <- rbind(polarization_scores, data.frame(
      net_name = net_name,
      ei = ei(preprocessed_g, ms, directed = TRUE, loops = F),
      modularity = modularity(preprocessed_g, ms, directed = TRUE),
      smi = smi_values,
      rwc =   RWC_polarization(as.undirected(preprocessed_g), ms, n_influencers, n_sim, n_walks),
      #arwc =   RWC_polarization(as.undirected(preprocessed_g), ms, adaptive_n_influencers, n_sim, n_walks),
      bp= calculate_boundary_polarization(preprocessed_g, ms),
      aei = aei_score,
      stringsAsFactors = FALSE
    ))
    
    # print("Compute scores for ER randomized graph")

    polarization_scores <- rbind(polarization_scores, data.frame(
      net_name = paste0(net_name,"_ER"),
      ei = ei(preprocessed_er_g, ms_er, directed = TRUE, loops = F),
      smi = smi_values_er,
      rwc = RWC_polarization(as.undirected(preprocessed_er_g), ms_er, n_influencers, n_sim, n_walks),
      arwc =   RWC_polarization(as.undirected(preprocessed_er_g), ms_er, adaptive_n_influencers, n_sim, n_walks),
      modularity = modularity(preprocessed_er_g, ms_er, directed = TRUE),
      bp= calculate_boundary_polarization(preprocessed_er_g, ms_er),
      aei = aei_score_er,
      stringsAsFactors = FALSE
    ))
    
    
    print("Compute scores for configuration model graphs")
    
    polarization_scores <- rbind(polarization_scores, data.frame(
      net_name = paste0(net_name,"_conf_1"),
      ei = ei(preprocessed_conf_g1, ms_conf1, directed = TRUE, loops = F),
      modularity =  modularity(preprocessed_conf_g1, ms_conf1, directed = TRUE),
      smi = smi_values_conf1,
      rwc =   RWC_polarization(as.undirected(preprocessed_conf_g1), ms_conf1, n_influencers, n_sim, n_walks),
      #arwc =   RWC_polarization(as.undirected(preprocessed_conf_g1), ms_conf1, adaptive_n_influencers, n_sim, n_walks),
      bp= calculate_boundary_polarization(preprocessed_conf_g1, ms_conf1),
      aei = aei_score_conf1,
      stringsAsFactors = FALSE
    ))
    
    polarization_scores <- rbind(polarization_scores, data.frame(
      net_name = paste0(net_name,"_conf_2"),
      ei = ei(preprocessed_conf_g2, ms_conf2, directed = TRUE, loops = F),
      modularity =  modularity(preprocessed_conf_g2, ms_conf2, directed = TRUE),
      smi = smi_values_conf2,
      rwc =   RWC_polarization(as.undirected(preprocessed_conf_g2), ms_conf2, n_influencers, n_sim, n_walks),
      #arwc =   RWC_polarization(as.undirected(preprocessed_conf_g2), ms_conf2, adaptive_n_influencers, n_sim, n_walks),
      bp= calculate_boundary_polarization(preprocessed_conf_g2, ms_conf2),
      aei = aei_score_conf2,
      stringsAsFactors = FALSE
    ))
    
    # polarization_scores <- rbind(polarization_scores, data.frame(
    #   net_name = paste0(net_name,"_conf_3"),
    #   ei = ei(preprocessed_conf_g3, ms_conf3, directed = TRUE, loops = F),
    #   modularity =  modularity(preprocessed_conf_g3, ms_conf3, directed = TRUE),
    #   smi = smi_values_conf3,
    #   rwc =   RWC_polarization(as.undirected(preprocessed_conf_g3), ms_conf3, n_influencers, n_sim, n_walks),
    #   arwc =   RWC_polarization(as.undirected(preprocessed_conf_g3), ms_conf3, adaptive_n_influencers, n_sim, n_walks),
    #   bp= calculate_boundary_polarization(preprocessed_conf_g3, ms_conf3),
    #   aei = aei_score_conf3,
    #   stringsAsFactors = FALSE
    # ))
    # 
    
    print("End of scores computation")
    print("----------------------")
    
    
  
  }
  
}

# Save the DataFrame as an Excel file
output_file <- 'total_scores.csv'
write.csv(polarization_scores, file = output_file, row.names = FALSE)


output_file <- 'total_labels.csv'
write.csv(network_labels, file = output_file, row.names = FALSE)


# Print a message indicating successful save
cat("DataFrame saved as CSV file:", output_file, "\n")

print(polarization_scores)


############ TRAINING MODELS ############ 
print("----------Training before normalization ----------------")


#LOGISTIC REGRESSION
  result_ml_model <- compute_accuracy(network_labels$label,polarization_scores)
lr_model <- result_ml_model$model
print("Model: Logistic Regression")
print(paste("Accuracy:",result_ml_model$accuracy))

#DECISION TREE
result_dt_model <- compute_ensemble_accuracy(network_labels$label,polarization_scores)
dt_model <- result_dt_model$model
print("Model: Decision Tree")
print(paste("Accuracy:",result_dt_model$accuracy))





################ NORMALIZE SCORES############

#1. Get scores of configuration models
# Filter rows referred to configuration models
conf_scores <- polarization_scores %>%
  filter(grepl("_conf_", net_name))

print(conf_scores)

#2. Compute mean and stdev for each score of each configuration model
#Create an empty data frame to store the results
conf_stats_df <- data.frame(
  net_name = character(),
  mean_ei = numeric(), sd_ei = numeric(),
  mean_modularity = numeric(), sd_modularity = numeric(),
  mean_smi = numeric(), sd_smi = numeric(),
  mean_rwc = numeric(), sd_rwc = numeric(),
  mean_bp = numeric(), sd_bp = numeric(),
  mean_aei = numeric(), sd_aei = numeric(),
  stringsAsFactors = FALSE
)

# Iterate over each pair of rows
for(i in seq(1, nrow(conf_scores), by = 2))  {
  # Extract the two rows
    row1 <- conf_scores[i, ]
    #assign to each row the original network name
    name <- paste0(strsplit(row1$net_name, split = ".")[1], "edgelist")
    # Combine results into a single row
    result_row <- data.frame(
      net_name = name,
      mean_ei = mean(conf_scores[c(i, i+1), "ei"]),
      sd_ei = sd(conf_scores[c(i, i+1), "ei"]),
      mean_modularity = mean(conf_scores[c(i, i+1), "modularity"]),
      sd_modularity = sd(conf_scores[c(i, i+1), "modularity"]),
      mean_smi = mean(conf_scores[c(i, i+1), "smi"]),
      sd_smi = sd(conf_scores[c(i, i+1), "smi"]),
      mean_rwc = mean(conf_scores[c(i, i+1), "rwc"]),
      sd_rwc = sd(conf_scores[c(i, i+1), "rwc"]),
      mean_bp = mean(conf_scores[c(i, i+1), "bp"]),
      sd_bp = sd(conf_scores[c(i, i+1), "bp"]),
      mean_aei = mean(conf_scores[c(i, i+1), "aei"]),
      sd_aei = sd(conf_scores[c(i, i+1), "aei"]),
      stringsAsFactors = FALSE
      
    )
    # Append the result row to the result data frame
    conf_stats_df <- rbind(conf_stats_df, result_row)
  }



#3. Replace Na with mean value of the column.
conf_stats_df[sapply(conf_stats_df, is.numeric)] <- lapply(conf_stats_df[sapply(conf_stats_df, is.numeric)], function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))

#4. subtract the std and mean from original graphs to normalize them.

#4.a. Select original network scores
original_scores_df <- polarization_scores %>%
  filter(!grepl("_conf", net_name) & !grepl("_ER", net_name))


original_labels_df <- network_labels %>%
  filter(!grepl("_conf", net_name) & !grepl("_ER", net_name))

#4.b.
mean_df <- conf_stats_df[, seq(2, ncol(conf_stats_df), by = 2)]
sd_df <- conf_stats_df[, seq(3, ncol(conf_stats_df), by = 2)]

#4.c normalize
normalize_original_df <- (original_scores_df[,2:ncol(original_scores_df)]-mean_df)/sd_df
normalize_original_df$net_name<-original_scores_df[,1]



remaining_scores <- polarization_scores %>%
  filter(grepl("_conf", net_name) | grepl("_ER", net_name))

train_df <- rbind(normalize_original_df, remaining_scores)
random_models_labels_df <- network_labels %>%
  filter(grepl("_conf", net_name) | grepl("_ER", net_name))
label_train_df <- rbind(original_labels_df,random_models_labels_df)



############ TRAINING MODELS on NORMALIZED DATA ############ 
print("----------Training after normalization ----------------")
#LOGISTIC REGRESSION
result_ml_model_n <- compute_accuracy(label_train_df$label,train_df)
lr_model_n<- result_ml_model_n$model
print("Model: Logistic Regression")
print(paste("Accuracy:",result_ml_model$accuracy))

#DECISION TREE
result_en_model <- compute_ensemble_accuracy(label_train_df$label,train_df)
dt_model_n<- result_en_model$model
print("Model: Decision Tree")
print(paste("Accuracy:",result_en_model$accuracy))









################## TEST on UNSEEN NETWORKS#################

# Set path to the folder containing your files
test_folder_path <- "network_data/test_data"

# Get list of all file names in the folder
file_list <- list.files(test_folder_path, full.names = TRUE)







print("------- Evaluate on Test Set -------")
data<-data.frame()
# Loop through the files and read them into R
for (net_file in file_list) {
  test_data <- data.frame()
  print("----------------------------")
  print(paste("Processing", net_file))
  data <- read.table(net_file, header = FALSE, sep = ",", col.names = c("retweeter", "retweeted", "Weight")) 
  
  # Convert the dataframe  to  numeric
  data <- as.data.frame(lapply(data, as.numeric))
  
  # Create a directed graph from the data
  G <- graph_from_data_frame(data, directed = TRUE)
  
  #summary of the graph
  summary(G)
  
  preprocessed_g <- preprocess_graph(G)
  
  #plot_clusterings(as.undirected(preprocessed_g), 2)
  #CREATE RANDOMIZED MODELS
  #d=0
  erdos_renyi_model <- generate_ER(preprocessed_g)
  preprocessed_er_g <-preprocess_graph(erdos_renyi_model)

  #d=1
  configuration_model <- generate_config_model(preprocessed_g)
  preprocessed_conf_g <-preprocess_graph(configuration_model)
  
  
  clust_fluid <- cluster_fluid_communities(preprocessed_g, no.of.communities=2) 
  clust_fluid_er <- cluster_fluid_communities(preprocessed_er_g, no.of.communities=2)
  clust_fluid_conf <- cluster_fluid_communities(preprocessed_conf_g, no.of.communities=2)
  

  
  
  # g<-ggraph(preprocessed_er_g, layout = "fr") +
  #   geom_edge_link() +
  #   geom_node_point(aes(color = factor(clust_fluid_er$membership))) +
  #   theme_void() +
  #   labs(title = "ER Graph",
  #        color = "Clusters") +
  #   theme(legend.position = "topright")
  # ggsave(paste0("img/",net_file,"_er.png"), plot = g, width = 6, height = 4, units = "in")
  # 
  # 
  # g<-ggraph(preprocessed_g, layout = "fr") +
  #   geom_edge_link() +
  #   geom_node_point(aes(color = factor(clust_fluid$membership))) +
  #   theme_void() +
  #   labs(title = "Initial Graph",
  #        color = "Clusters") +
  #   theme(legend.position = "topright")
  # ggsave(paste0("img/",net_file,".png"), plot = g, width = 6, height = 4, units = "in")
  # 
  # g<-ggraph(preprocessed_conf_g, layout = "fr") +
  #   geom_edge_link() +
  #   geom_node_point(aes(color = factor(clust_fluid_conf$membership))) +
  #   theme_void() +
  #   labs(title = "Configuration Model Graph",
  #        color = "Clusters") +
  #   theme(legend.position = "topright")
  # ggsave(paste0("img/",net_file,"_conf.png"), plot = g, width = 6, height = 4, units = "in")
  # 

  
  ################## POLARIZATION SCORES #################
  
  ms <- membership(clust_fluid)
  ms_er <- membership(clust_fluid_er)
  ms_conf <- membership(clust_fluid_conf)
  
  
  V(preprocessed_g)$community <- ms
  V(preprocessed_er_g)$community <- ms_er
  V(preprocessed_conf_g)$community <- ms_conf
  
  V(preprocessed_g)$selected_community <- ifelse(ms == 1, 1, ifelse(ms == 2, 2, NA))  # Apply the smi function to the directed graph using the new attribute
  smi_values <- smi(preprocessed_g,"community", normalize = TRUE)
  cluster_sizes<-table(clust_fluid$membership)
  smi_values <- sum(smi_values * cluster_sizes) / sum(cluster_sizes)
  
  
  V(preprocessed_er_g)$selected_community <- ifelse(ms_er == 1, 1, ifelse(ms_er == 2, 2, NA))  # Apply the smi function to the directed graph using the new attribute
  smi_values_er<- smi(preprocessed_er_g,"community", normalize = TRUE)
  cluster_sizes<-table(clust_fluid_er$membership)
  smi_values_er <- sum(smi_values_er * cluster_sizes) / sum(cluster_sizes)
  # 
  # 
  V(preprocessed_conf_g)$selected_community <- ifelse(ms_conf == 1, 1, ifelse(ms_conf == 2, 2, NA))  # Apply the smi function to the directed graph using the new attribute
  smi_values_conf <- smi(preprocessed_conf_g,"community", normalize = TRUE)
  cluster_sizes<-table(clust_fluid_conf$membership)
  smi_values_conf <- sum(smi_values_conf * cluster_sizes) / sum(cluster_sizes)
  
  
  
  aei_score <- AEI_polarization(preprocessed_g, ms)
  aei_score_er <- AEI_polarization(preprocessed_er_g, ms_er)
  aei_score_conf <- AEI_polarization(preprocessed_conf_g, ms_conf)
  
  
  #compute scores
  test_data <- rbind(test_data, data.frame(
    ei = ei(preprocessed_g, ms, directed = TRUE, loops = F),
    modularity =  modularity(preprocessed_g,ms, directed = TRUE),
    smi = smi_values_conf,
    rwc =   RWC_polarization(as.undirected(preprocessed_g), ms, n_influencers, n_sim, n_walks),
    bp= calculate_boundary_polarization(preprocessed_g, ms),
    aei = aei_score_conf,
    stringsAsFactors = FALSE
  ))
  #compute scores
  test_data <- rbind(test_data, data.frame(
    ei = ei(preprocessed_er_g, ms_er, directed = TRUE, loops = F),
    modularity =  modularity(preprocessed_er_g, ms_er, directed = TRUE),
    smi = smi_values_conf,
    rwc =   RWC_polarization(as.undirected(preprocessed_er_g), ms_er, n_influencers, n_sim, n_walks),
    bp= calculate_boundary_polarization(preprocessed_er_g, ms_er),
    aei = aei_score_conf,
    stringsAsFactors = FALSE
  ))
  #compute scores
  test_data <- rbind(test_data, data.frame(
    ei = ei(preprocessed_conf_g, ms_conf, directed = TRUE, loops = F),
    modularity =  modularity(preprocessed_conf_g, ms_conf, directed = TRUE),
    smi = smi_values_conf,
    rwc =   RWC_polarization(as.undirected(preprocessed_conf_g), ms_conf, n_influencers, n_sim, n_walks),
    bp= calculate_boundary_polarization(preprocessed_conf_g, ms_conf),
    aei = aei_score_conf,
    stringsAsFactors = FALSE
  ))
  # 
  data<-rbind(data,test_data)
  print("-----------Logistic regression -----------")
  y_pred <- predict(lr_model, newdata = test_data, type = "response")
  predicted_labels <- ifelse(y_pred > threshold, 1, 0)
  
  print(paste("Prediction: ", predicted_labels[1]))
  print(paste("Final prediction: ", y_pred[1]))
  print(paste("Prediction ER: ", predicted_labels[2]))
  print(paste("Final prediction ER: ", y_pred[2]))
  print(paste("Prediction configuration model: ", predicted_labels[3]))
  print(paste("Final prediction configuration model: ", y_pred[3]))
  

  
  print("-----------decision tree -----------")
  y_pred <- predict(dt_model, newdata = test_data, type = "class")
  print(paste("Prediction: ", y_pred[1]))
  print(paste("Prediction ER network: ", y_pred[2]))
  print(paste("Prediction configuration network: ", y_pred[3]))

  
  print("-----------Logistic regression after normalization -----------")
  y_pred <- predict(lr_model_n, newdata = test_data, type = "response")
  predicted_labels<- ifelse(y_pred > threshold, 1, 0)
  
  print(paste("Prediction: ", predicted_labels[1]))
  print(paste("Final prediction: ", y_pred[1]))
  print(paste("Prediction ER: ", predicted_labels[2]))
  print(paste("Final prediction ER: ", y_pred[2]))
  print(paste("Prediction configuration model: ", predicted_labels[3]))
  print(paste("Final prediction configuration model: ", y_pred[3]))
  
  
  
  
  
  print("-----------decision tree after normalization -----------")
  y_pred <- predict(dt_model_n, newdata = test_data, type = "class")

  print(paste("Prediction: ", y_pred[1]))
  print(paste("Prediction ER network: ", y_pred[2]))
  print(paste("Prediction configuration network: ", y_pred[3]))
  # 
  # 
  print("-------------------")
}


print(data)
xtable(data, type = "latex")








################## CONSENSUS CLUSTERING#################
# 
# #Louvain
# clust_louvain <- membership(multilevel.community(preprocessed_g, resolution = 0.1 )) #the lower the resolution the bigger the clustersclust_louvain<-cut_at(cluster_edge_betweenness(simple_G), 2)
# #Leiden
# clust_leiden <- membership(cluster_leiden(preprocessed_g, objective_function = "modularity", weights = E(preprocessed_g)$weight,  
#                                           resolution_parameter = 0.2,  n_iterations = 100))
# #clust_leiden <-  cut_at(cluster_walktrap(simple_G), 2)
# # Apply the function to the list of clusterings
# best_clusterings <- list(Louvain = clust_louvain, Leiden = clust_leiden)#, rsc = rsc$cluster)
# cons_mat <- consensus_matrix(best_clusterings)
# 
# cluster_assignments <- apply(cons_mat, 1, which.max)
# 
# V(preprocessed_g)$cluster <- cluster_assignments
# 
# # Calculate the total number of nodes in the graph
# #total_nodes <- vcount(simple_G)
# 
# #n. of nodes in each cluster
# cluster_counts <- table(V(preprocessed_g)$cluster)
# 
# # Plot only the nodes in selected clusters
# ggraph(preprocessed_g, layout = "fr") +
#   geom_edge_link() +
#   geom_node_point(aes(color = factor(V(preprocessed_g)$cluster))) +
#   theme_void() +
#   labs(title = paste0("Consensus Clustering"), color = "Clusters") +
#   theme(legend.position = "topright")
# 
# # # Visualize the consensus matrix using heatmap
# # heatmap(cons_mat, symm = TRUE, col = heat.colors(10))
# # 
# # # Visualize the consensus matrix using ggplot2
# # library(ggplot2)
# # library(reshape2)
# # cons_mat_melt <- melt(cons_mat)
# # ggplot(cons_mat_melt, aes(x = Var1, y = Var2, fill = value)) +
# #   geom_tile() +
# #   scale_fill_gradient(low = "white", high = "red") +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# 
# 
# 
# #Now, consider only the subgraph induce by the 2 biggest communities and plot it.
# 
# largest_clusters <- names(sort(cluster_counts, decreasing = TRUE)[1:2])
# 
# # Filter nodes belonging to the two largest clusters
# selected_nodes <- V(preprocessed_g)$name[V(preprocessed_g)$cluster %in% largest_clusters]
# 
# # Create a subgraph with only the selected nodes
# final_graph <- induced_subgraph(preprocessed_g, selected_nodes)
# 
# # Plot the final graph
# ggraph(final_graph, layout = "fr") +
#   geom_edge_link() +
#   geom_node_point(aes(color = factor(V(final_graph)$cluster))) +
#   theme_void() +
#   labs(title = paste0("Consensus Clustering"), color = "Clusters") +
#   theme(legend.position = "topright")






