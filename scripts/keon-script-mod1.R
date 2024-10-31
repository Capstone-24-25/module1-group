# Choose a larger number (more than ten) of top predictive proteins using each selection method
# We will try choosing 25 proteins for each method!

setwd("C:/Users/Keon School/OneDrive/Documents/GitHub/module1-group")

library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')

## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

#ttests_out

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 25) %>%
  pull(protein)

proteins_s1

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422) # Set the same seed to reproduce the same analysis
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 25) %>%
  pull(protein)

proteins_s2


corr_data <- biomarker_clean %>% 
  filter(!is.na(ados))

corr_data

proteins_s3 <- corr_data %>%
  pivot_longer(cols = -c(ados, group),
               names_to = 'protein',
               values_to = 'level') %>%
  group_by(protein) %>%
  summarize(correlation = cor(ados, level)) %>%
  arrange(desc(abs(correlation))) %>%
  slice_head(n = 25)  # Select top 25 proteins

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

proteins_sstar

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422) # Set the same seed to reproduce the same analysis
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')



# Original intersected proteins: "DERM"  "RELT"  "IgD"   "FSTL1"

# Model metrics on original analysis (on testing set): 
# Sensitivity = 0.875
# Specificity = 0.8
# Accuracy = 0.839
# Roc_auc = 0.908


# New intersected proteins: "DERM", "RELT", "Calcineurin", "IgD", "PTN", 
# "FSTL1", "MAPK2", "TGF-b R III", "MMP-2", "gp130, soluble", "Notch 1", "ALCAM", "MATN2"

# Model metrics on larger model: 
# Sensitivity = 0.812
# Specificity = 0.867
# Accuracy = 0.839
# Roc_auc = 0.946

# Accuracy is the same, and ROC_AUC is higher, which suggest that the large model performed better. 
# Sensitivity decreased while specificity increased, suggesting more false positives and fewer false negatives. 





# Fuzzy intersection instead of hard intersection

#fuzzy_proteins_s1
#fuzzy_proteins_s2


rf_importance <- importance(rf_out)
rf_importance <- as.data.frame(rf_importance)


# Normalize random forest importance scores
rf_importance$norm_importance <- (rf_importance$MeanDecreaseGini - min(rf_importance$MeanDecreaseGini)) / 
  (max(rf_importance$MeanDecreaseGini) - min(rf_importance$MeanDecreaseGini))

# Convert p-values to a score (1 - p_value), we then compare with normalized rf importance scores.
ttests_out$p_score <- 1 - ttests_out$p_value

# Merge the two data frames
merged_results <- merge(rf_importance, ttests_out, by.x = "row.names", by.y = "protein", all = TRUE)
colnames(merged_results)[1] <- "protein"
#merged_results

# Calculate fuzzy intersection score
merged_results$fuzzy_intersection <- pmin(merged_results$norm_importance, merged_results$p_score)

# Sort by fuzzy intersection score
merged_results <- merged_results[order(-merged_results$fuzzy_intersection), ]
#merged_results

# Select variables with fuzzy intersection score above 0.5
fuzzy_proteins <- merged_results[merged_results$fuzzy_intersection > 0.5, "protein"]

fuzzy_proteins # This fuzzy intersection gives us a panel of5 proteins: "DERM", "IgD", "TGF-b R III", "MAPK14", "FSTL1"


biomarker_sstar_2 <- biomarker_clean %>%
  select(group, any_of(fuzzy_proteins)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422) # Set the same seed to reproduce the same analysis
biomarker_split2 <- biomarker_sstar_2 %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit2 <- glm(class ~ ., 
           data = training(biomarker_split2), 
           family = 'binomial')

# evaluate errors on test set
class_metrics2 <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

testing(biomarker_split2) %>%
  add_predictions(fit2, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')


# With this method of fuzzy intersection we found a panel of five proteins. We got metrics of: 
# sensitivity: 0.75 
# specificity: 0.8  
# accuracy: 0.774
# roc_auc: 0.888