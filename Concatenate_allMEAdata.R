########## Code for making UMAPs using preprocessed data ##########
#please read through the comments to change/check the necessary things before starting the code

#Clean Global Environment
rm(list = ls())

#install pacakges and load them
install.packages("readxl")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("writexl")
install.packages("janitor")
install.packages("tidyr")

library(tidyr)
library(readxl)
library(tidyverse)
library(dplyr)
library(writexl)
library(janitor)
library(writexl)


#set working directory
#correct the "/" to the right direction or go to "Files > More > Set as working directory" in the right bottom panel
setwd("insert path")

################################# DATA PREP #################################

#get list of all files in the directory
file_list <- list.files(path = getwd(), pattern = "Well_Summary_Parameters_hfNB_s70\\.csv$", recursive = TRUE, full.names = TRUE)
#include certain files
pattern_to_include <- c('adjusted')
file_list <- file_list[grep(paste(pattern_to_include, collapse = '|'), file_list)]
patterns_to_include <- c("DIV35", "DIV42", "DIV49")
file_list <- file_list[grep(paste(patterns_to_include, collapse = '|'), file_list)]
patterns_to_include <- c("101-4964", "101-4908", "101-4947", "101-4960", "100-2110", "118-2257", "118-2253", "118-2251")
file_list <- file_list[grep(paste(patterns_to_include, collapse = '|'), file_list)]

#create an empty list to store the data frames
data_list <- list()

#loop through each file and read in the data
for (i in seq_along(file_list)) {
  data_list[[i]] <- read.csv(file_list[i], header = TRUE)

}

#combine all data frames into one big data frame, remove any rows in which the Treatment/ID is NA
#note: if this is the case, go back to your raw recording and re-record with the plate map, this should resolve the issue
data <- bind_rows(data_list)
#make parameter names the column header and not the first row
#colnames(data) <- data[1,]
#data <- data[-1, ]
data <- data[!is.na(data$`Phenotype`), , drop = FALSE]

#remove the exclusion wells and rename Phenotypes, check with print if everything is correct
data <- data %>%
  filter(Phenotype != "Exclude") %>%
  mutate(Phenotype = if_else(Phenotype == 'HET1', 'CACNA1A+/- (1)', Phenotype)) %>%
  mutate(Phenotype = if_else(Phenotype == 'HET2', 'CACNA1A+/- (2)', Phenotype)) %>%
  mutate(Phenotype = if_else(Phenotype == 'PAT2-Rescue', 'PAT2 Rescue', Phenotype)) %>%
  mutate(Phenotype = if_else(Phenotype == 'PAT2 ', 'PAT2', Phenotype)) #%>%
  #filter(!grepl('-[A-Za-z]', Phenotype))

print(unique(data$Phenotype))

#clean up the data 
#remove wells with no spikes and/or network bursts, remove empty columns, and remove double/irrelevant parameters
data <- data[data$`Number_of_spikes`!= 0, , drop = FALSE]
data <- data %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))
#data <- data[!is.na(data$`Network Burst Frequency`), , drop = FALSE ]
#data <- remove_empty(data, which = c("cols"), cutoff = 1, quiet = FALSE)
#col_remove <- c("Network Bursts Ignored Flag", "Start Electrode", "Percent Bursts with Start Electrode", 
                # "Network Normalized Duration IQR", "Number of Spikes", "Number of Active Electrodes",
                # "Number of Bursting Electrodes", "Number of Bursts", "Number of Network Bursts",
                # "Number of Elecs Participating in Burst - Avg", "Number of Elecs Participating in Burst - Std",
                # "Number of Elecs Participating in Burst - Median", "Number of Elecs Participating in Burst - MAD",
                # "Number of Spikes per Network Burst per Channel - Avg", "Number of Spikes per Network Burst per Channel - Std",
                # "Number of Spikes per Network Burst per Channel - Median", "Number of Spikes per Network Burst per Channel - MAD",
                # "Median/Mean ISI within Burst - Avg", "Median/Mean ISI within Burst - Std", 
                # "Median/Mean ISI within Network Burst - Avg", "Median/Mean ISI within Network Burst - Std")
#data <- data[, !names(data) %in% col_remove]

#check if there are any NA values left and remove them if that is the case
any(is.na(data))
#data <- na.omit(data)

#create new columns that contains the cell lines
#you can change and add to this to match your data set
data <- data %>%
  mutate(
    Cell_line = case_when(
      str_detect(`Phenotype`, "WTC") ~ "WTC", #finds "Cell line 1" in Treatment/ID and adds CTRL to the new Cell_line column
      str_detect(`Phenotype`, "VUS1") ~ "VUS1",
      str_detect(`Phenotype`, "VUS2") ~ "VUS2",
      str_detect(`Phenotype`, "CACNA1A\\+/- \\(1\\)") ~ "HET1",
      str_detect(`Phenotype`, "VUS3") ~ "VUS3",
      str_detect(`Phenotype`, "CACNA1A\\+/- \\(2\\)") ~ "HET2",
      str_detect(`Phenotype`, "PAT2 Rescue") ~ "PAT2-Rescue",      
      str_detect(`Phenotype`, "PAT2") ~ "PAT2",
      str_detect(`Phenotype`, "PAT1") ~ "PAT1",
      str_detect(`Phenotype`, "PAT3") ~ "PAT3",
      str_detect(`Phenotype`, "V1393M Mimic") ~ "PAT4-Mimic",
      str_detect(`Phenotype`, "V1393M PAT") ~ "PAT4",
      str_detect(`Phenotype`, "V1393M Rescue") ~ "PAT4-Rescue",
      TRUE ~ NA_character_),
    
    #uncomment and change the following lines if you want to add a Treatment column
    Treatment = "Basal" #Set default to "Basal" for the entire column 
  ) %>%
  mutate(Treatment = if_else(str_detect(`Phenotype`, "DMSO"), "DMSO", Treatment)) %>%
  mutate(Treatment = if_else(str_detect(`Phenotype`, "PD212"), "PD212", Treatment)) %>%
  mutate(Treatment = if_else(str_detect(`Phenotype`, "NS309"), "NS309", Treatment)) %>%
  mutate(Treatment = if_else(str_detect(`Phenotype`, "CNV944"), "CNV944", Treatment))
    
any(is.na(data))

#re-ordering the column so that all defining information is at the front
data <- data[, c(1:3, 43:49, 4:42)]

#export the final data set
#correct the "/" to the right direction or go to "Files > More > Set as working directory" in the right bottom panel
write_xlsx(data,"insert path/name.xlsx")

#Calculate Avg number of fragments per NB
data$Avg_Fragments_per_NB <- ((data$Number_of_Network_Bursts + data$Number_of_hf_Network_Bursts)/data$Number_of_Network_Bursts)
write_xlsx(data,"insert path/name_avg fragments per NB.xlsx")

