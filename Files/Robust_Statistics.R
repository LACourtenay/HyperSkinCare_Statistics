
# Code for Hyperspectral Band Analysis Statistics - HYPERSKINCARE Project
#
# Code by Lloyd A. Courtenay - ladc1995@gmail.com
# TIDOP Research Group (http://tidop.usal.es/) - University of Salamanca
#
# Last Update: 09/04/2021
#

# Libraries ---------------------------

# check to see whether the required packages are installed on your system. If they are not, they will be installed

for (i in c("ggplot2", "tibble", "dplyr", "tidyverse", "e1071",
            "philentropy", "gridExtra", "RVAideMemoire", "asbio",
            "car", "caret", "pROC")) {
  if(i %in% rownames(installed.packages()) == TRUE) {
    print(paste("The package named", i, "has already been installed", sep = " "))
  } else {install.packages(i)}
}; rm(i)

# load packages

library(ggplot2)
library(tibble)
library(dplyr)
library(tidyverse)
library(e1071)
library(philentropy)
library(gridExtra)
library(RVAideMemoire)
library(asbio)
library(car)
library(caret)
library(pROC)

#

# functions ----------------------------

# statistics function

# P-Value calibration functions

# Bayes factor bounds

BFB<-function(p){
  return((1/(-exp(1)*p*log(p))))
}

# posterior odds calculation

post_odds<-function(p, priors = 0.5){
  return(priors * BFB(p))
}

# false positive risk

FPR<-function(p, priors = 0.5){
  pH = priors/(1-priors)
  return(1/(1+(pH*BFB(p))))
}

# pU(Ha|p)

p_BFB<-function(p){
  return((1/(-exp(1)*p*log(p)))/(1+(1/(-exp(1)*p*log(p)))))
}

# p(H0)

p_H0<-function(p, priors = 0.5){
  if(p <= 0.3681) {
    value = FPR(p, priors = priors)
  } else {
    value = 1-FPR(p, priors = 1-priors)
  }
  return(value)
}; p_H0<-Vectorize(p_H0)

# calculate quantile

quantile_function<-function(a, p = 0.5) {
  answer = sort(a)[ceiling(p * length(a))]
  return(answer)
}

# calculate normality data

normality_tests <- function(data) {
  normality_results <- tibble(
    Band_name = as.character(),
    Shapiro_w = as.numeric(),
    Shapiro_p = as.numeric(),
    Shapiro_Lower_pH0 = as.numeric(),
    Shapiro_pH0 = as.numeric(),
    Shapiro_Upper_pH0 = as.numeric(),
    Skew = as.numeric(),
    Kurtosis = as.numeric()
  )
  BAND_data <- data[,3:ncol(data)]
  for (i in 1:ncol(BAND_data)){
    target <- BAND_data[[paste("BAND_", i, sep = "")]]
    normality_results <- normality_results %>%
      add_row(
        Band_name = paste("BAND_", i, sep = ""),
        Shapiro_w = shapiro.test(target)$statistic[[1]],
        Shapiro_p = shapiro.test(target)$p.value[[1]],
        Shapiro_Lower_pH0 = p_H0(shapiro.test(target)$p.value[[1]], priors = 0.8),
        Shapiro_pH0 = p_H0(shapiro.test(target)$p.value[[1]], priors = 0.5),
        Shapiro_Upper_pH0 = p_H0(shapiro.test(target)$p.value[[1]], priors = 0.2),
        Skew = skewness(target),
        Kurtosis = kurtosis(target)
      )
  }
  normality_results <- normality_results %>%
    mutate(Band_name = factor(1:n()))
  return(normality_results)
}

# calculation for residuals

residual_calculation <- function(data){
  residual_results <- tibble(
    Band_name = as.character(),
    Res_count = as.numeric()
  )
  BAND_data <- data[,3:ncol(data)]
  for (i in 1:ncol(BAND_data)) {
    target <- BAND_data[[paste("BAND_", i, sep = "")]]
    linear <- lm(target ~ data$Sample)
    res <- 1 - sum(linear$residuals ^ 2) / sum((target - mean(target)) ^ 2)
    residual_results <- residual_results %>%
      add_row(
        Band_name = paste("BAND_", i, sep = ""),
        Res_count = res
      )
  }
  residual_results <- residual_results %>%
    mutate(Band_name = factor(1:n()))
  return(residual_results)
}

# calculate signature

signature_data <- function(data, type) {
  BAND_data <- data[,3:ncol(data)]
  if (type == "gaussian") {
    signature <- tibble(
      Band_name = as.character(),
      Lower_CI = as.numeric(),
      Central = as.numeric(),
      Upper_CI = as.numeric(),
      Deviation = as.numeric()
    )
    for (i in 1:ncol(BAND_data)){
      target <- BAND_data[[paste("BAND_", i, sep = "")]]
      signature <- signature %>%
        add_row(
          Band_name = paste("BAND_", i, sep = ""),
          Lower_CI = quantile_function(target, p = 0.05),
          Central = mean(target),
          Upper_CI = quantile_function(target, p = 0.95),
          Deviation = sd(target)
        )
    }
    signature <- signature %>%
      mutate(Band_name = factor(1:n()))
    return(signature)
  }
  else if (type == "robust") {
    signature <- tibble(
      Band_name = as.character(),
      Lower_CI = as.numeric(),
      Central = as.numeric(),
      Upper_CI = as.numeric(),
      Deviation = as.numeric()
    )
    for(i in 1:ncol(BAND_data)){
      target <- BAND_data[[paste("BAND_", i, sep = "")]]
      signature <- signature %>%
        add_row(
          Band_name = paste("BAND_", i, sep = ""),
          Lower_CI = quantile_function(target, p = 0.05),
          Central = median(target),
          Upper_CI = quantile_function(target, p = 0.95),
          Deviation = sqrt(r.bw(target)$`S.xx`[[1]])
        )
    }
    signature <- signature %>%
      mutate(Band_name = factor(1:n()))
    return(signature)
  }
  else {"Choose type as either 'robust' or 'gaussian"}
}

# homoscedasticity tests

homoscedasticity_test <- function(data, method) {
  BAND_data <- data[, 3:ncol(data)]
  if (method == "levene") {
    homosc_results <- tibble(
      Band_name = as.character(),
      Test_Statistic = as.numeric(),
      p_Value = as.numeric()
    )
    for (i in 1:ncol(BAND_data)) {
      target <- BAND_data[[paste("BAND_", i, sep = "")]]
      homosc_results <- homosc_results %>%
        add_row(
          Band_name = paste("BAND_", i, sep = ""),
          Test_Statistic = leveneTest(target, data$Sample, location = "median")$`F value`[[1]],
          p_Value = leveneTest(target, data$Sample, location = "median")$`Pr(>F)`[[1]]
        )
    }
    homosc_results <- homosc_results %>%
      mutate(Band_name = factor(1:n()),
             pBFB = if_else(
               p_Value <= 0.3681,
               p_BFB(p_Value),
               NA_real_
             ),
             lower_post_odds = if_else(
               p_Value <= 0.3681,
               post_odds(p_Value, priors = 0.8),
               NA_real_
             ),
             post_odds = if_else(
               p_Value <= 0.3681,
               post_odds(p_Value, priors = 0.5),
               NA_real_
             ),
             upper_post_odds = if_else(
               p_Value <= 0.3681,
               post_odds(p_Value, priors = 0.2),
               NA_real_
             ),
             lower_FPR = if_else(
               p_Value <= 0.3681,
               FPR(p_Value, priors = 0.8),
               NA_real_
             ),
             FPR = if_else(
               p_Value <= 0.3681,
               FPR(p_Value, priors = 0.5),
               NA_real_
             ),
             upper_FPR = if_else(
               p_Value <= 0.3681,
               FPR(p_Value, priors = 0.2),
               NA_real_
             ),
             lower_pH0 = p_H0(p_Value, priors = 0.8),
             pH0 = p_H0(p_Value, priors = 0.5),
             upper_pH0 = p_H0(p_Value, priors = 0.2)
      )
    return(homosc_results)
  }
  else if (method == "bartlett") {
    homosc_results <- tibble(
      Band_name = as.character(),
      Test_Statistic = as.numeric(),
      p_Value = as.numeric()
    )
    for (i in 1:ncol(BAND_data)) {
      target <- BAND_data[[paste("BAND_", i, sep = "")]]
      homosc_results <- homosc_results %>%
        add_row(
          Band_name = paste("BAND_", i, sep = ""),
          Test_Statistic = bartlett.test(target, data$Sample)$statistic[[1]],
          p_Value = bartlett.test(target, data$Sample)$p.value[[1]]
        )
    }
    homosc_results <- homosc_results %>%
      mutate(Band_name = factor(1:n()),
             pBFB = if_else(
               p_Value <= 0.3681,
               p_BFB(p_Value),
               NA_real_
             ),
             lower_post_odds = if_else(
               p_Value <= 0.3681,
               post_odds(p_Value, priors = 0.8),
               NA_real_
             ),
             post_odds = if_else(
               p_Value <= 0.3681,
               post_odds(p_Value, priors = 0.5),
               NA_real_
             ),
             upper_post_odds = if_else(
               p_Value <= 0.3681,
               post_odds(p_Value, priors = 0.2),
               NA_real_
             ),
             lower_FPR = if_else(
               p_Value <= 0.3681,
               FPR(p_Value, priors = 0.8),
               NA_real_
             ),
             FPR = if_else(
               p_Value <= 0.3681,
               FPR(p_Value, priors = 0.5),
               NA_real_
             ),
             upper_FPR = if_else(
               p_Value <= 0.3681,
               FPR(p_Value, priors = 0.2),
               NA_real_
             ),
             lower_pH0 = p_H0(p_Value, priors = 0.8),
             pH0 = p_H0(p_Value, priors = 0.5),
             upper_pH0 = p_H0(p_Value, priors = 0.2)
      )
    return(homosc_results)
  }
  else {print("Choose a method between 'levene' and 'bartlett'")}
}

# jensen-shannon distance calculations

JSD_Calculations <- function(data) {
  split_data <- split(data, data$Sample)
  data_1 <- split_data[[1]]
  data_1 <- data_1[,3:ncol(data_1)]
  data_2 <- split_data[[2]]
  data_2 <- data_2[,3:ncol(data_2)]
  results <- c()
  for (i in 1:ncol(data_1)) {
    dist_1 <- data_1[[paste("BAND_", i, sep = "")]]
    dist_2 <- data_2[[paste("BAND_", i, sep = "")]]
    prob_1 <- hist(dist_1)$counts / nrow(data_1)
    prob_2 <- hist(dist_2)$counts / nrow(data_2)
    dist <- distance(rbind(prob_1, prob_2), method = "jensen-shannon")[[1]]
    results <- c(results, dist)
  }
  distance_values <- tibble(
    Band_name = paste("BAND_", i, sep = ""),
    Similarity = results
  ) %>% 
    mutate(Band_name = factor(1:n()))
  return(distance_values)
}

# train-test machine learning split function

split.data = function(data, p = 0.7, s = 666) {
  set.seed(2)
  index = sample (1:dim(data)[1])
  train = data [index[1: floor(dim(data)[1]*p)],]
  test = data [index[((ceiling(dim(data)[1]*p)) + 1) :dim(data)[1]],]
  return(list(train = train, test = test))
}

# machine learning univariate filter

ml_feature_selection <- function(data, min_band, max_band, method){
  ctrl <- trainControl(method = "repeatedcv",
                       number = 10, repeats = 10,
                       classProbs = TRUE, summaryFunction = twoClassSummary)
  ml_results <- c()
  raw_data <- data[(min_band+2):(max_band+2)]
  samples <- data$Sample
  for (i in min_band:max_band){
    target_raw_data <- raw_data[[paste("BAND_", i, sep = "")]]
    target_data <- data.frame(Data = target_raw_data, Sample = samples)
    allset <- split.data(target_data, p = 0.7)
    trainset <- allset$train
    testset <- allset$test
    glm.model <- train(Sample ~ Data, data = target_data, method = method,
                       metric = "ROC", trControl = ctrl)
    glm.probs <- predict(glm.model, testset, type = "prob")
    glm.ROC <- roc(response = testset$Sample, predictor = glm.probs[,2], levels = levels(testset$Sample))
    ml_results <- c(ml_results, glm.ROC$auc[[1]])
  }
  ml_results <- tibble(AUC = ml_results) %>%
    mutate(Band_name = factor(min_band:max_band))
  return(ml_results)
}

# multivariate tests

multivariate_tests <- function(data, min_band, max_band, method) {
  if (method == "Wilcox") {
    pairwise.wilcox.test(as.matrix(data[,(min_band+2):(max_band+2)]), data$Sample, p.adjust.method = "BH")
  } else if (method == "MANOVA") {
    pairwise.perm.manova(as.matrix(data[,(min_band+2):(max_band+2)]), data$Sample, test = c("Hotelling-Lawley"), nperm = 999, 
                         progress = FALSE, p.method = "none")
  } else {print("Select either 'Wilcox' or 'MANOVA'")}
}

#

# load data ---------------------------

a <- read.csv(url("https://raw.githubusercontent.com/LACourtenay/HyperSkinCare_Statistics/main/Files/Hyperspectral_Data.csv"))
a$Sample <- factor(a$Sample); a$Patient_ID <- factor(a$Patient_ID)

# seperate data according to sample

split <- split(a, a$Sample)

SCC<-split$SCC; SCC<-droplevels(SCC, exclude = c("H","BCC"))
BCC<-split$BCC; BCC<-droplevels(BCC, exclude = c("H","SCC"))
H<-split$Healthy; H<-droplevels(H, exclude = c("SCC","BCC"))

H_SCC<-rbind(split$Healthy, split$SCC); H_SCC<-droplevels(H_SCC,exclude = "BCC")
H_BCC<-rbind(split$Healthy, split$BCC); H_BCC<-droplevels(H_BCC,exclude = "SCC")
SCC_BCC<-rbind(split$SCC, split$BCC); SCC_BCC <- droplevels(SCC_BCC, exclude = "H")

rm(split)

#

# Check data normality ----------------

# Calculate Shapiro results, skewness and kurtosis

H_normality_results <- normality_tests(H)
BCC_normality_results <- normality_tests(BCC)
SCC_normality_results <- normality_tests(SCC)

# View numeric results

View(H_normality_results)
View(BCC_normality_results)
View(SCC_normality_results)

#

# Prepare and view normality plots ----------------

# Shapiro W and P value Plots

shap_pvalues_SCC<-ggplot(data = SCC_normality_results,
                    aes(x = Band_name,
                        y = log(Shapiro_p))) +
  geom_bar(stat = "identity", fill = "#999999") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        plot.title = element_text(face = "bold")) +
  xlab("Band") +
  labs(title = "SCC - Shapiro Test (p Values)") +
  geom_hline(yintercept = log(0.003), colour = "black", size = 1)
shap_pvalues_BCC<-ggplot(data = BCC_normality_results,
                    aes(x = Band_name,
                        y = log(Shapiro_p))) +
  geom_bar(stat = "identity", fill = "#999999") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        plot.title = element_text(face = "bold")) +
  xlab("Band") +
  labs(title = "BCC - Shapiro Test (p Values)") +
  geom_hline(yintercept = log(0.003), colour = "black", size = 1)
shap_pvalues_H<-ggplot(data = H_normality_results,
                  aes(x = Band_name,
                      y = log(Shapiro_p))) +
  geom_bar(stat = "identity", fill = "#999999") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        plot.title = element_text(face = "bold")) +
  xlab("Band") +
  labs(title = "Healthy - Shapiro Test (p Values)") +
  geom_hline(yintercept = log(0.003), colour = "black", size = 1)
shap_wvalues_H<-ggplot(data = H_normality_results,
                  aes(x = Band_name,
                      y = Shapiro_w)) +
  geom_bar(stat = "identity", fill = "#999999") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        plot.title = element_text(face = "bold")) +
  xlab("Band") +
  labs(title = "Healthy - Shapiro Test (w Values)") +
  coord_cartesian(ylim = c(0.91,1))
shap_wvalues_SCC<-ggplot(data = SCC_normality_results,
                    aes(x = Band_name,
                        y = Shapiro_w)) +
  geom_bar(stat = "identity", fill = "#999999") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        plot.title = element_text(face = "bold")) +
  xlab("Band") +
  labs(title = "SCC - Shapiro Test (w Values)") +
  coord_cartesian(ylim = c(0.91,1))
shap_wvalues_BCC<-ggplot(data = BCC_normality_results,
                    aes(x = Band_name,
                        y = Shapiro_w)) +
  geom_bar(stat = "identity", fill = "#999999") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        plot.title = element_text(face = "bold")) +
  xlab("Band") +
  labs(title = "BCC - Shapiro Test (w Values)") +
  coord_cartesian(ylim = c(0.91,1))
skew_values_BCC<-ggplot(data = BCC_normality_results,
                        aes(x = Band_name,
                            y = Skew)) +
  geom_bar(stat = "identity", fill = "#999999") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        plot.title = element_text(face = "bold")) +
  xlab("Band") +
  labs(title = "BCC - Skewness")
skew_values_SCC<-ggplot(data = SCC_normality_results,
                        aes(x = Band_name,
                            y = Skew)) +
  geom_bar(stat = "identity", fill = "#999999") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        plot.title = element_text(face = "bold")) +
  xlab("Band") +
  labs(title = "SCC - Skewness")
skew_values_H<-ggplot(data = H_normality_results,
                      aes(x = Band_name,
                          y = Skew)) +
  geom_bar(stat = "identity", fill = "#999999") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        plot.title = element_text(face = "bold")) +
  xlab("Band") +
  labs(title = "Healthy - Skewness")
kurt_values_H<-ggplot(data = H_normality_results,
                      aes(x = Band_name,
                          y = Kurtosis)) +
  geom_bar(stat = "identity", fill = "#999999") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        plot.title = element_text(face = "bold")) +
  xlab("Band") +
  labs(title = "Healthy - Kurtosis") +
  geom_hline(yintercept = 0, colour = "black", size = 1)
kurt_values_SCC<-ggplot(data = SCC_normality_results,
                        aes(x = Band_name,
                            y = Kurtosis)) +
  geom_bar(stat = "identity", fill = "#999999") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        plot.title = element_text(face = "bold")) +
  xlab("Band") +
  labs(title = "SCC - Kurtosis") +
  geom_hline(yintercept = 0, colour = "black", size = 1)
kurt_values_BCC<-ggplot(data = BCC_normality_results,
                        aes(x = Band_name,
                            y = Kurtosis)) +
  geom_bar(stat = "identity", fill = "#999999") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        plot.title = element_text(face = "bold")) +
  xlab("Band") +
  labs(title = "BCC - Kurtosis") +
  geom_hline(yintercept = 0, colour = "black", size = 1)

#

H_shap_ph0_plot<-ggplot(data = H_normality_results,
                        aes(x = Band_name, y = Shapiro_pH0, group = 1)) +
  geom_ribbon(aes(ymin = Shapiro_Lower_pH0, ymax = Shapiro_Upper_pH0),
              alpha = 0.2) +
  geom_line(size = 0.75) +
  theme_bw() +
  xlab("Band") +
  ylab("Probability of Null Hypothesis") +
  ggtitle("Shapiro Test Results - Healthy") +
  theme(plot.margin = unit(c(1,1,0.5,1), "cm"),
        axis.title = element_text(face = "bold", size = 15),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.text.x = element_text(margin = margin(t = 5))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(breaks = seq(0, 260, 50))
BCC_shap_ph0_plot<-ggplot(data = BCC_normality_results,
                          aes(x = Band_name, y = Shapiro_pH0, group = 1)) +
  geom_ribbon(aes(ymin = Shapiro_Lower_pH0, ymax = Shapiro_Upper_pH0),
              alpha = 0.2) +
  geom_line(size = 0.75) +
  theme_bw() +
  xlab("Band") +
  ylab("Probability of Null Hypothesis") +
  ggtitle("Shapiro Test Results - BCC") +
  theme(plot.margin = unit(c(0.5,1,1,1), "cm"),
        axis.title = element_text(face = "bold", size = 15),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.text.x = element_text(margin = margin(t = 5))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(breaks = seq(0, 260, 50))
SCC_shap_ph0_plot<-ggplot(data = SCC_normality_results,
                          aes(x = Band_name, y = Shapiro_pH0, group = 1)) +
  geom_ribbon(aes(ymin = Shapiro_Lower_pH0, ymax = Shapiro_Upper_pH0),
              alpha = 0.2) +
  geom_line(size = 0.75) +
  theme_bw() +
  xlab("Band") +
  ylab("Probability of Null Hypothesis") +
  ggtitle("Shapiro Test Results - SCC") +
  theme(plot.margin = unit(c(0.5,1,1,1), "cm"),
        axis.title = element_text(face = "bold", size = 15),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.text.x = element_text(margin = margin(t = 5))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(breaks = seq(0, 260, 50))

#

X11(); grid.arrange(shap_pvalues_H, shap_wvalues_H,
                    shap_pvalues_SCC, shap_wvalues_SCC,
                    shap_pvalues_BCC, shap_wvalues_BCC,
                    nrow = 3,
                    ncol = 2)

rm(shap_pvalues_H, shap_wvalues_H, shap_pvalues_SCC, shap_wvalues_SCC, shap_pvalues_BCC, shap_wvalues_BCC) # clean memory

X11(); grid.arrange(skew_values_H, kurt_values_H,
                    skew_values_SCC, kurt_values_SCC,
                    skew_values_BCC, kurt_values_BCC,
                    nrow = 3,
                    ncol = 2)

rm(skew_values_H, kurt_values_H, skew_values_SCC, kurt_values_SCC, skew_values_BCC, kurt_values_BCC) # clean memory

X11(); grid.arrange(H_shap_ph0_plot, BCC_shap_ph0_plot, SCC_shap_ph0_plot,
                    ncol = 2, nrow = 2)

rm(H_shap_ph0_plot, BCC_shap_ph0_plot, SCC_shap_ph0_plot)

#

# Calculate and plot residuals -----------------------

C<-rbind(SCC, BCC); C$Sample <- fct_collapse(C$Sample, C = c("SCC", "BCC"))

ggplot(data = residual_calculation(rbind(H, C)), aes(x = Band_name, y = Res_count, group = 1)) +
  geom_line(stat = "identity", colour = "black", size = 1) +
  theme_bw() +
  theme(plot.margin = unit(c(1,0.5,1,1), "cm"),
        axis.title = element_text(face = "bold", size = 15),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", size = 20),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.text.x = element_text(margin = margin(t = 5))
  ) +
  xlab("Band") +
  scale_x_discrete(breaks = seq(0, 260, 50)) +
  labs(title = "Residuals from Linear Gaussian Models")

#

# Calculate signatures ---------------

H_signature <- signature_data(H, type = "robust") # select either robust signature or gaussian signature
SCC_signature <- signature_data(SCC, type = "robust") # select either robust signature or gaussian signature
BCC_signature <- signature_data(BCC, type = "robust") # select either robust signature or gaussian signature

# view numeric results

View(H_signature)
View(SCC_signature)
View(BCC_signature)

# Visualise plot

grid.arrange(ggplot() +
               geom_line(data = SCC_signature, aes(x = Band_name, y = Central, group = 1),
                         size = 0.75, colour = "#FF0000") +
               geom_line(data = SCC_signature, aes(x = Band_name, y = Central + Upper_CI, group = 1),
                         size = 0.5, colour = "#FF0000", linetype = "solid") +
               geom_line(data = SCC_signature, aes(x = Band_name, y = Central - Lower_CI, group = 1),
                         size = 0.5, colour = "#FF0000", linetype = "solid") +
               geom_line(data = BCC_signature, aes(x = Band_name, y = Central, group = 1),
                         size = 0.75, colour = "#000000") +
               geom_line(data = BCC_signature, aes(x = Band_name, y = Central + Upper_CI, group = 1),
                         size = 0.5, colour = "#000000", linetype = "solid") +
               geom_line(data = BCC_signature, aes(x = Band_name, y = Central - Lower_CI, group = 1),
                         size = 0.5, colour = "#000000", linetype = "solid") +
               geom_line(data = H_signature, aes(x = Band_name, y = Central, group = 1),
                         size = 0.75, colour = "#3399FF") +
               geom_line(data = H_signature, aes(x = Band_name, y = Central + Upper_CI, group = 1),
                         size = 0.5, colour = "#3399FF", linetype = "solid") +
               geom_line(data = H_signature, aes(x = Band_name, y = Central - Lower_CI, group = 1),
                         size = 0.5, colour = "#3399FF", linetype = "solid") +
               theme_bw() +
               theme(plot.margin = unit(c(1,0.5,1,1), "cm"),
                     axis.title = element_text(face = "bold", size = 15),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.y = element_text(margin = margin(r = 10)),
                     axis.text = element_text(size = 12, face = "bold"),
                     plot.title = element_text(face = "bold", size = 20),
                     axis.text.y = element_text(margin = margin(r = 5)),
                     axis.text.x = element_text(margin = margin(t = 5))
               ) +
               ggtitle("Robust Signature") +
               xlab("Band_name") +
               ylab("Intensity") +
               coord_cartesian(ylim = c(0, 80)) +
               scale_x_discrete(breaks = seq(0, 260, 50)), ggplot() +
               geom_line(data = SCC_signature, aes(x = Band_name, y = Deviation, group = 1),
                         size = 0.75, colour = "#FF0000") +
               geom_line(data = BCC_signature, aes(x = Band_name, y = Deviation, group = 1),
                         size = 0.75, colour = "#000000") +
               geom_line(data = H_signature, aes(x = Band_name, y = Deviation, group = 1),
                         size = 0.75, colour = "#3399FF") +
               theme_bw() +
               theme(plot.margin = unit(c(1,1,1,0.5), "cm"),
                     axis.title = element_text(face = "bold", size = 15),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.y = element_text(margin = margin(r = 10)),
                     axis.text = element_text(size = 12, face = "bold"),
                     plot.title = element_text(face = "bold", size = 20),
                     axis.text.y = element_text(margin = margin(r = 5)),
                     axis.text.x = element_text(margin = margin(t = 5))
               ) +
               ggtitle("Robust Variance") +
               xlab("Band_name") +
               ylab("Intensity") +
               coord_cartesian(ylim = c(4, 12.1)) +
               scale_x_discrete(breaks = seq(0, 260, 50)),
             ncol = 2, nrow = 1)

# clear memory

rm(H_signature, SCC_signature, BCC_signature)


#

# Hypothesis Testing ---------------

# homoscedasticity

H_BCC_homosc <- homoscedasticity_test(H_BCC, method = "levene") # select either levene or bartlett test
H_SCC_homosc <- homoscedasticity_test(H_SCC, method = "levene") # select either levene or bartlett test
SCC_BCC_homosc <- homoscedasticity_test(SCC_BCC, method = "levene") # select either levene or bartlett test

# view numeric homoscedasticity results

View(H_BCC_homosc)
View(H_SCC_homosc)
View(SCC_BCC_homosc)

# prepare homoscedasticity p-value plots

homosc_H_BCC_plot<-ggplot(data = H_BCC_homosc,
                          aes(x = Band_name, y = pH0, group = 1)) +
  geom_ribbon(aes(ymin = lower_pH0, ymax = upper_pH0),
              alpha = 0.2) +
  geom_line(size = 0.75) +
  theme_bw() +
  xlab("Band") +
  ylab("P(H0)") +
  ggtitle("H vs BCC") +
  theme(plot.margin = unit(c(1,1,0.5,1), "cm"),
        axis.title = element_text(face = "bold", size = 15),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.text.x = element_text(margin = margin(t = 5))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(breaks = seq(0, 260, 50))
homosc_H_SCC_plot<-ggplot(data = H_SCC_homosc,
                          aes(x = Band_name, y = pH0, group = 1)) +
  geom_ribbon(aes(ymin = lower_pH0, ymax = upper_pH0),
              alpha = 0.2) +
  geom_line(size = 0.75) +
  theme_bw() +
  xlab("Band") +
  ylab("P(H0)") +
  ggtitle("H vs SCC") +
  theme(plot.margin = unit(c(1,1,0.5,1), "cm"),
        axis.title = element_text(face = "bold", size = 15),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.text.x = element_text(margin = margin(t = 5))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(breaks = seq(0, 260, 50))
homosc_SCC_BCC_plot<-ggplot(data = SCC_BCC_homosc,
                            aes(x = Band_name, y = pH0, group = 1)) +
  geom_ribbon(aes(ymin = lower_pH0, ymax = upper_pH0),
              alpha = 0.2) +
  geom_line(size = 0.75) +
  theme_bw() +
  xlab("Band") +
  ylab("P(H0)") +
  ggtitle("SCC vs BCC") +
  theme(plot.margin = unit(c(1,1,0.5,1), "cm"),
        axis.title = element_text(face = "bold", size = 15),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.text.x = element_text(margin = margin(t = 5))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(breaks = seq(0, 260, 50))

# prepare homoscedasticity test statistic plots

homosc_H_SCC_test_stat_plot<-ggplot(data = H_SCC_homosc, 
                                    aes(x = Band_name,
                                        group = 1,
                                        y = Test_Statistic)) +
  geom_line(stat = "identity", size = 0.75) +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,0.5,1), "cm"),
        axis.title = element_text(face = "bold", size = 15),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.text.x = element_text(margin = margin(t = 5))) +
  labs(title = "H vs SCC") +
  theme(plot.title = element_text(face = "bold")) +
  scale_x_discrete(breaks = seq(0, 260, 50))
homosc_H_BCC_test_stat_plot<-ggplot(data = H_BCC_homosc, 
                                    aes(x = Band_name,
                                        group = 1,
                                        y = Test_Statistic)) +
  geom_line(stat = "identity", size = 0.75) +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,0.5,1), "cm"),
        axis.title = element_text(face = "bold", size = 15),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.text.x = element_text(margin = margin(t = 5))) +
  labs(title = "H vs BCC") +
  theme(plot.title = element_text(face = "bold")) +
  scale_x_discrete(breaks = seq(0, 260, 50))
homosc_SCC_BCC_test_stat_plot<-ggplot(data = SCC_BCC_homosc, 
                                      aes(x = Band_name,
                                          group = 1,
                                          y = Test_Statistic)) +
  geom_line(stat = "identity", size = 0.75) +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,0.5,1), "cm"),
        axis.title = element_text(face = "bold", size = 15),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.text.x = element_text(margin = margin(t = 5))) +
  labs(title = "SCC vs BCC") +
  theme(plot.title = element_text(face = "bold")) +
  scale_x_discrete(breaks = seq(0, 260, 50))

# View plots

grid.arrange(homosc_H_SCC_plot, homosc_H_SCC_test_stat_plot,
             homosc_H_BCC_plot, homosc_H_BCC_test_stat_plot,
             homosc_SCC_BCC_plot, homosc_SCC_BCC_test_stat_plot,
             nrow = 3, ncol = 2)

# clear memory

rm(homosc_H_SCC_plot, homosc_H_SCC_test_stat_plot,
   homosc_H_BCC_plot, homosc_H_BCC_test_stat_plot,
   homosc_SCC_BCC_plot, homosc_SCC_BCC_test_stat_plot,
   H_BCC_homosc, H_SCC_homosc, SCC_BCC_homosc)

#

# Jensen-Shannon Divergence ------------------------------

H_BCC_JSD <- suppressMessages(suppressWarnings(JSD_Calculations(H_BCC)))
H_SCC_JSD <- suppressMessages(suppressWarnings(JSD_Calculations(H_SCC)))
SCC_BCC_JSD <- suppressMessages(suppressWarnings(JSD_Calculations(SCC_BCC)))

ggplot() +
  geom_line(data = H_BCC_JSD,
            aes(x = Band_name, y = Similarity, group = 1), size = 0.75,
            color = "black") +
  geom_line(data = H_SCC_JSD,
            aes(x = Band_name, y = Similarity, group = 1), size = 0.75,
            color = "red") +
  geom_line(data = SCC_BCC_JSD,
            aes(x = Band_name, y = Similarity, group = 1), size = 0.75,
            color = "blue") +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(face = "bold", size = 15),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", size = 20),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.text.x = element_text(margin = margin(t = 5))
  ) +
  ggtitle("Jenson-Shannon Divergence") +
  xlab("Band") +
  ylab("Similarity") +
  scale_x_discrete(breaks = seq(0, 260, 50))

# clear memory

rm(H_BCC_JSD, H_SCC_JSD, SCC_BCC_JSD)

#

# Univariate Machine Learning Feature Selection ---------------------------------------

H_SCC_ml_results <- suppressMessages(ml_feature_selection(H_SCC, min_band = 53, max_band = 205,
                                                          method = "glm") # choose either method = "svmRadial" or "glm"
)
H_BCC_ml_results <- suppressMessages(ml_feature_selection(H_BCC, min_band = 53, max_band = 205,
                                                          method = "glm") # choose either method = "svmRadial" or "glm"
)
SCC_BCC_ml_results <- suppressMessages(ml_feature_selection(SCC_BCC, min_band = 53, max_band = 205,
                                                          method = "glm") # choose either method = "svmRadial" or "glm"
)

ggplot() +
  geom_line(data = H_BCC_ml_results,
            aes(x = Band_name, y = AUC, group = 1), size = 0.75,
            color = "black") +
  geom_line(data = H_SCC_ml_results,
            aes(x = Band_name, y = AUC, group = 1), size = 0.75,
            color = "red") +
  geom_line(data = SCC_BCC_ml_results,
            aes(x = Band_name, y = AUC, group = 1), size = 0.75,
            color = "blue") +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(face = "bold", size = 15),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", size = 20),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.text.x = element_text(margin = margin(t = 5))
  ) +
  ggtitle("Receiver Operating Characteristics") +
  xlab("Band") +
  ylab("Area Under Curve") +
  coord_cartesian(ylim = c(0.45, 0.65)) +
  scale_x_discrete(breaks = seq(50, 210, 20))

# clear memory

rm(H_SCC_ml_results, H_BCC_ml_results, SCC_BCC_ml_results)

#

# Multivariate Hypothesis entire spectrum ---------------------------------------

mult_diff <- multivariate_tests(a, min_band = 1, max_band = 260,
                   method = "Wilcox") # select either Wilcox or MANOVA

# p-value table

mult_diff

# H vs BCC p-value calibrations

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("posterior odds with priors of ", i, " calculated at: ",
              post_odds(mult_diff$p.value[1,1], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("False Positive Risk with priors of ", i, " calculated at: ",
              FPR(mult_diff$p.value[1,1], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("Probability of Null Hypothesis with priors of ", i, " calculated at: ",
              p_H0(mult_diff$p.value[1,1], priors = i),
              sep = ""))
}

# H vs SCC p-value calibrations

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("posterior odds with priors of ", i, " calculated at: ",
              post_odds(mult_diff$p.value[2,2], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("False Positive Risk with priors of ", i, " calculated at: ",
              FPR(mult_diff$p.value[2,2], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("Probability of Null Hypothesis with priors of ", i, " calculated at: ",
              p_H0(mult_diff$p.value[2,2], priors = i),
              sep = ""))
}

# BCC vs SCC p-value calibrations

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("posterior odds with priors of ", i, " calculated at: ",
              post_odds(mult_diff$p.value[2,1], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("False Positive Risk with priors of ", i, " calculated at: ",
              FPR(mult_diff$p.value[2,1], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("Probability of Null Hypothesis with priors of ", i, " calculated at: ",
              p_H0(mult_diff$p.value[2,1], priors = i),
              sep = ""))
}

# clear memory

rm(mult_diff)

#

# Multivariate Hypothesis between 573 and 779nm ---------------------------------------

mult_diff <- multivariate_tests(a, min_band = 75, max_band = 168,
                   method = "Wilcox") # select either Wilcox or MANOVA

# p-value table

mult_diff

# H vs BCC p-value calibrations

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("posterior odds with priors of ", i, " calculated at: ",
              post_odds(mult_diff$p.value[1,1], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("False Positive Risk with priors of ", i, " calculated at: ",
              FPR(mult_diff$p.value[1,1], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("Probability of Null Hypothesis with priors of ", i, " calculated at: ",
              p_H0(mult_diff$p.value[1,1], priors = i),
              sep = ""))
}

# H vs SCC p-value calibrations

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("posterior odds with priors of ", i, " calculated at: ",
              post_odds(mult_diff$p.value[2,2], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("False Positive Risk with priors of ", i, " calculated at: ",
              FPR(mult_diff$p.value[2,2], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("Probability of Null Hypothesis with priors of ", i, " calculated at: ",
              p_H0(mult_diff$p.value[2,2], priors = i),
              sep = ""))
}

# BCC vs SCC p-value calibrations

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("posterior odds with priors of ", i, " calculated at: ",
              post_odds(mult_diff$p.value[2,1], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("False Positive Risk with priors of ", i, " calculated at: ",
              FPR(mult_diff$p.value[2,1], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("Probability of Null Hypothesis with priors of ", i, " calculated at: ",
              p_H0(mult_diff$p.value[2,1], priors = i),
              sep = ""))
}

# clear memory

rm(mult_diff)

#

# Multivariate Hypothesis between 429 and 520 nm ---------------------------------------

mult_diff <- multivariate_tests(a, min_band = 10, max_band = 51,
                   method = "Wilcox") # select either Wilcox or MANOVA

# p-value table

mult_diff

# H vs BCC p-value calibrations

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("posterior odds with priors of ", i, " calculated at: ",
              post_odds(mult_diff$p.value[1,1], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("False Positive Risk with priors of ", i, " calculated at: ",
              FPR(mult_diff$p.value[1,1], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("Probability of Null Hypothesis with priors of ", i, " calculated at: ",
              p_H0(mult_diff$p.value[1,1], priors = i),
              sep = ""))
}

# H vs SCC p-value calibrations

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("posterior odds with priors of ", i, " calculated at: ",
              post_odds(mult_diff$p.value[2,2], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("False Positive Risk with priors of ", i, " calculated at: ",
              FPR(mult_diff$p.value[2,2], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("Probability of Null Hypothesis with priors of ", i, " calculated at: ",
              p_H0(mult_diff$p.value[2,2], priors = i),
              sep = ""))
}

# BCC vs SCC p-value calibrations

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("posterior odds with priors of ", i, " calculated at: ",
              post_odds(mult_diff$p.value[2,1], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("False Positive Risk with priors of ", i, " calculated at: ",
              FPR(mult_diff$p.value[2,1], priors = i),
              sep = ""))
}

for (i in c(0.2, 0.5, 0.8)) {
  print(paste("Probability of Null Hypothesis with priors of ", i, " calculated at: ",
              p_H0(mult_diff$p.value[2,1], priors = i),
              sep = ""))
}

# clear memory

rm(mult_diff)

#


