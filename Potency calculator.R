#Potency Calculator

library(reshape)
library(ggplot2)
library(manipulate)
library(dplyr)
library(pheatmap)
library(tidyr)


setwd("C:/Users/Patrick Niekamp/Desktop/Paper/Additional")
#read in file with average values for IDO1, CSF1, CD63 and CCL2 for each donor
data <- read.csv("Potency_Markers.csv", na.strings = c("", "NA"))

# Pre-processing
data <- data[, colSums(is.na(data)) < nrow(data)]
data <- data[rowSums(is.na(data)) < ncol(data), ]
data <- data.frame(data)
data$Donor_ID <- factor(data$Donor_ID, levels = c(257, 320, 43, 127, 17, 74, 75, 63))

# Determine quartiles
score_quartile <- function(column) {
  q <- quantile(column, probs = c(0.25, 0.5, 0.75))  # Calculate quartiles
  score <- rep(0, length(column))  # Initialize scores with 1 (1st quartile)
  score[column > q[1]] <- 1  # 2nd quartile
  score[column > q[2]] <- 2  # 3rd quartile
  score[column > q[3]] <- 3  # 4th quartile
  return(score)
}


# Unweighted matrix
# Add scores
data$IDO1_score <- score_quartile(data$IDO1) 
data$CSF1_score <- score_quartile(data$CSF1) 
data$CD63_score <- score_quartile(data$CD63) 
data$CCL2_score <- score_quartile(data$CCL2)

# Calculate the total score by summing individual scores
data$total_score <- data$IDO1_score + data$CSF1_score + data$CD63_score + data$CCL2_score

# View the scored data
print(data[, c("Donor_ID", "IDO1_score", "CSF1_score", "CD63_score", "CCL2_score", "total_score")])
write.csv(data, file = "unweighted_matrix.csv")

data_long <- data %>%
  pivot_longer(cols = c(IDO1_score, CSF1_score, CD63_score, CCL2_score),
               names_to = "Score_Type",
               values_to = "Score_Value")

# Fig 7E / Unweighted Matrix
ggplot(data_long, aes(x = Donor_ID, y = Score_Value, fill = Score_Type)) +
  geom_bar(stat = "identity") +
  labs(title = "Scores by Donor ID",
       x = "Donor ID",
       y = "Score Value") +
  theme_minimal() +
  scale_x_discrete(drop = FALSE)


# Weighted matrix 

data <- read.csv("Potency_Markers.csv", na.strings = c("", "NA"))

# Pre-processing
data <- data[, colSums(is.na(data)) < nrow(data)]
data <- data[rowSums(is.na(data)) < ncol(data), ]
data <- data.frame(data)
data$Donor_ID <- factor(data$Donor_ID, levels = c(257, 320, 43, 127, 17, 74, 75, 63))

# Determine quartiles
score_quartile <- function(column) {
  q <- quantile(column, probs = c(0.25, 0.5, 0.75))  # Calculate quartiles
  score <- rep(0, length(column))  # Initialize scores with 1 (1st quartile)
  score[column > q[1]] <- 1  # 2nd quartile
  score[column > q[2]] <- 2  # 3rd quartile
  score[column > q[3]] <- 3  # 4th quartile
  return(score)
}

# Apply the scoring function with additional multipliers 
data$IDO1_score <- score_quartile(data$IDO1) *2
data$CSF1_score <- score_quartile(data$CSF1) 
data$CD63_score <- score_quartile(data$CD63) *2
data$CCL2_score <- score_quartile(data$CCL2) *1/2

# Calculate the total score by summing individual scores
data$total_score <- data$IDO1_score + data$CSF1_score + data$CD63_score + data$CCL2_score

# View the scored data
print(data[, c("Donor_ID", "IDO1_score", "CSF1_score", "CD63_score", "CCL2_score", "total_score")])
write.csv(data, file = "weighted_matrix.csv")

data_long <- data %>%
  pivot_longer(cols = c(IDO1_score, CSF1_score, CD63_score, CCL2_score),
               names_to = "Score_Type",
               values_to = "Score_Value")

# Fig 7F / Weighted Matrix
ggplot(data_long, aes(x = Donor_ID, y = Score_Value, fill = Score_Type)) +
  geom_bar(stat = "identity") +
  labs(title = "Scores by Donor ID",
       x = "Donor ID",
       y = "Score Value") +
  theme_minimal() +
  scale_x_discrete(drop = FALSE)


