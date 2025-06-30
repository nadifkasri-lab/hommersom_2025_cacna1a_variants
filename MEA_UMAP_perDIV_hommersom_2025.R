########## Code for making UMAPs using preprocessed data ##########
#please read through the comments to change/check the necessary things before starting the code

#Clean Global Environment
rm(list = ls())

#install pacakges and load them
install.packages("umap")
install.packages("plotly")
install.packages("readxl")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("writexl")
install.packages("Rtsne")
install.packages("ggplot2")
install.packages("janitor")
install.packages("factoextra")
install.packages("tidyr")
install.packages("RColorBrewer")
install.packages("dendextend")

library(tidyr)
library(plotly) 
library(umap)
library(readxl)
library(tidyverse)
library(dplyr)
library(writexl)
library(Rtsne)
library(ggplot2)
library(janitor)
library(writexl)
library(factoextra)
library(RColorBrewer)
library(dendextend)

#set working directory
#correct the "/" to the right direction or go to "Files > More > Set as working directory" in the right bottom panel
setwd("insert path")

################################# DATA PREP #################################
#Load the data and filter out patient cell lines
data <- read_excel("insert/folder/name.xlsx")
data <- data %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))
data <- data %>%
  filter(Treatment == "Basal") %>%
  filter(Cell_line %in% c("PAT4-Rescue", "PAT2-Rescue",
                          "PAT1", "PAT2", "PAT3", "PAT4",
                          "VUS1", "VUS2", "VUS3"
                          ) | 
        Cell_line == "WTC" & str_detect(PT_all_path, "HDR"))
print(unique(data$Phenotype))
#clean up the data 
#remove wells with no spikes, remove empty columns, and remove double/irrelevant parameters
data <- data[data$`Number_of_spikes`!= 0, , drop = FALSE]
data <- data %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))
#check if there are any NA values left and remove them if that is the case
any(is.na(data))


################################# DATA BLINDING #################################

#remove all defining information, e.g., cell line, div, batch
#note: make sure to check your data to change the column from which you start data blinding
#in this case, the the blinded data will include everything from column 7 on
datablinded <- data[ , c(11:ncol(data))]

#make data numeric
datablinded <- apply(datablinded, 2, as.numeric)
datablinded <- as.data.frame(datablinded)

#scale data to make sure all parameters have the same influence on the final UMAP
datablinded <- scale(datablinded)
datablinded <- as.data.frame(datablinded)


#export scaled data set 
#correct the "/" to the right direction or go to "Files > More > Set as working directory" in the right bottom panel
#write_xlsx(datablinded,"insert folder/name.xlsx") 

################################# RUNNING THE PCA #################################

#run PCA on data first
pca <- prcomp(datablinded) 
summary(pca)

#make a scree plot to see which number of PCs define your data variance and save it
pdf("insert name.pdf")
plot(1:length(pca$sdev), pca$sdev^2, type = "b", 
     main = "Scree Plot", xlab = "Principal Component", ylab = "Eigenvalue",
     col = "purple")
abline(h = 1, col = "red", lty = 2)
dev.off()

#make a graph to see which parameters are contributing to top 2 PC and save it
pdf("insert name.pdf")
fviz_pca_var(pca, axes = c(1, 2), col.var = "contrib", repel = TRUE, labelsize = 2)
dev.off()

# Specify the PCs of interest
selected_pcs <- c(1, 2, 38, 39)

# Initialize an empty list to store results
top_contributions_list <- list()

# Loop through the selected PCs
for (pc_index in selected_pcs) {
  # Extract loadings for the current PC
  loadings_pca <- pca$rotation[, pc_index]
  
  # Create a data frame to display variable contributions for the current PC
  contributions_df <- data.frame(
    Variable = colnames(datablinded),
    Contribution = loadings_pca
  )
  
  # Sort by absolute value of contributions and take the top 10
  top_contributions <- head(contributions_df[order(abs(contributions_df$Contribution), decreasing = TRUE), ], 10)
  
  # Add a column indicating the current PC
  top_contributions$PC <- paste0("PCA", pc_index)
  
  # Save the results for the current PC into the list
  top_contributions_list[[pc_index]] <- top_contributions
}

# Combine the results for all selected PCs into a single data frame
PCA_parameters <- do.call(rbind, top_contributions_list)
write_xlsx(PCA_parameters,"insert path/name.xlsx")

# Extract top 10 parameters for PCs 1 and 2
top_params_pc1 <- head(contributions_df[order(abs(pca$rotation[, 1]), decreasing = TRUE), "Variable"], 10)
top_params_pc2 <- head(contributions_df[order(abs(pca$rotation[, 2]), decreasing = TRUE), "Variable"], 10)
top_params_12 <- unique(c(top_params_pc1, top_params_pc2))

# Extract top 10 parameters for PCs 38 and 39
top_params_pc38 <- head(contributions_df[order(abs(pca$rotation[, 38]), decreasing = TRUE), "Variable"], 10)
top_params_pc39 <- head(contributions_df[order(abs(pca$rotation[, 39]), decreasing = TRUE), "Variable"], 10)
top_params_3839 <- unique(c(top_params_pc38, top_params_pc39))

# Find unique parameters in PCs 1 and 2 that are not in PCs 38 and 39
unique_in_12 <- setdiff(top_params_12, top_params_3839)

# Find unique parameters in PCs 38 and 39 that are not in PCs 1 and 2
unique_in_3839 <- setdiff(top_params_3839, top_params_12)

# Print results
cat("Unique parameters in PCs 1 and 2, but not in PCs 38 and 39:\n")
print(unique_in_12)

cat("\nUnique parameters in PCs 38 and 39, but not in PCs 1 and 2:\n")
print(unique_in_3839)


################################# UMAP TIME #################################

#create a directory to store the UMAPs if it doesn't exist
umap_results_path <- "insert foldername/UMAPs/"
if (!file.exists(umap_results_path)) {
  dir.create(umap_results_path, recursive = TRUE)
}

#define the desired orders and colours
#your desired order cell line should be based on Treatment/ID
desired_order_cell_line <- c("WTC",  
                             "PAT2-Rescue", 
                             "PAT4-Rescue", "PAT2", 
                             "PAT3", "PAT1",
                             #"PAT4-Mimic", 
                             "PAT4", "VUS1", "VUS2", "VUS3")
desired_order_DIV <- c("DIV35", "DIV42", "DIV49")
color_order <- c("#606060", 
                 "#a0a0a4", 
                 "#fbfbfb","#8cc63f", 
                 "#2c6a31",  "#56B4e9",
                 #"#f66d9b", 
                 "#9b59b6", "#e69f00","#d55e00", "#b20e00")

#define the number of PCs you want to use in the UMAP analysis
#this number is based on the PCA that was ran before, check the results to determine the pc_num
pc_nums <- c(7)

for (pc_num in pc_nums){

#create a dataframe where the PCs can be found
scores_pcs <- predict(pca, newdata = datablinded)[, 1:pc_num]

#if you want to find the PC score again later, you can combine it with your data
data <- cbind(data, scores_pcs)

#define the columns where the PCs can be found
pcadata <- scores_pcs
pcadata <- na.omit(pcadata)
pcadata <- apply(pcadata, 2, as.numeric)
pcadata <- as.data.frame(pcadata)

#calculate umap and merge coordinates with the data set
data.umap = umap(pcadata, n_components = 2, mind_dist = 0.005, random_state = 2)
layout <- data.umap[["layout"]] 
layout <- data.frame(layout) 

#combine the UMAP1 and UMAP2 into the data so that we can correlate the UMAP to other info for each wells
final_pca <- cbind(layout, data) 

#export the data with PCA and UMAP, 
#correct the "/" to the right direction or go to "Files > More > Set as working directory" in the right bottom panel
#write_xlsx(final_pca,"~/OneDrive - Radboudumc/Test data set MEA workshop/Exported data/20240909_Test data MEA workshop_data-pca-umap.xlsx")

#find unique DIVs in your data frame
unique_divs <- unique(final_pca$DIV_range)

#loop through each level of DIVs and create a plot
for (div_level in unique_divs) {
  #filter data for the current div level
  filtered_data <- subset(final_pca, DIV_range == div_level)
  #arrange data to desired_order_cell_line
  filtered_data$Cell_line <- factor(filtered_data$Cell_line, levels = desired_order_cell_line)
  filtered_data <- filtered_data[order(filtered_data$Cell_line), ]
  
  #create the UMAP plot for the filtered data
  umap_plot <- ggplot(filtered_data, aes(x = X1, y = X2, color = factor(`Cell_line`, levels = desired_order_cell_line),
                                         fill = factor(Cell_line, levels = desired_order_cell_line))) + 
    geom_point(alpha = 0.9, show.legend = TRUE, shape = 21, size = 3, stroke = 0.2, color = 'black') +
    #geom_label(aes(label = `Number_of_spikes`, fill = factor(`Cell_line`, levels = desired_order_cell_line)), 
               #color = "black", size = 2, show.legend = FALSE) +
    scale_color_manual(values = color_order) +
    scale_fill_manual(values = color_order) + 
    labs(title = paste('Dimension reduction of MEA recordings for', div_level), x = 'UMAP1', y = 'UMAP2') +
    theme_bw() +
    #expand_limits(x = c(min(filtered_data$X1) - 0.5, max(filtered_data$X1) + 0.5), 
    #              y = c(min(filtered_data$X2) - 0.5, max(filtered_data$X2) + 0.5)) +
    scale_x_continuous(
      breaks = seq(-6, 6, by = 2), # Set tick positions
      limits = c(-6.2, 6.2),          # Set axis limits
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq(-10, 6, by = 2), # Set tick positions
      limits = c(-10.2, 7.2),            # Set axis limits
      expand = c(0, 0)
    )
  
  #display the plot
  print(umap_plot)
  
  #save the plot to a file
  umap_plot_file <- paste("UMAP_", div_level, "_PC", pc_num, ".pdf", sep = "")
  umap_plot_file_path <- file.path(umap_results_path, umap_plot_file)
  ggsave(umap_plot_file_path, umap_plot, width = 8, height = 5)
  
  cat("UMAP analysis completed for", div_level, "with PC =", pc_num, "\n")
}
}
