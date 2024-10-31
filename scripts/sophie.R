# problem 3
# install.packages("randomForest")

library(tidyverse)
library(tidymodels)
library(modelr)
library(rsample)
library(yardstick)
library(purrr)
library(dplyr)
library(randomForest)

load("./data/biomarker-clean.RData")

# partition
set.seed(3435)

partitions <- biomarker_clean %>%
  initial_split(prop = 0.8)

train_bio <- training(partitions)
test_bio <- testing(partitions)

asd_clean <- train_bio %>% 
  select(-ados)

# t tests
asd_nested <- asd_clean %>%
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  nest(data = c(level, group))

asd_nested %>% head(4)

# compute for groups
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

tt_out <- asd_nested %>%
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value)

# multiple testing corrections
m <- nrow(tt_out)
hm <- log(m) + 1/(2*m) - digamma(1)

tt_corrected <- tt_out %>%
  select(data, protein, p_value) %>%
  mutate(rank = row_number()) %>%
  mutate(p_bh = p_value*m/rank,
         p_by = p_value*m*hm/rank,
         p_bonf = p_value*m)

# top 10 for multiple testing
s1 <- tt_corrected %>%
  select(protein, p_by) %>%
  slice_min(order_by = p_by, n = 10) %>%
  pull(protein)

print("Top 10 proteins for multiple testing:")
s1

# random forest
asd_response <- asd_clean$group
asd_preds <- asd_clean[, -1] 
rf_out <- randomForest(x = asd_preds, # predictors
                       y = as.factor(asd_response), # response
                       ntree = 100, # number of trees
                       importance = T) # compute importance

# select most important predictors
s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

print("Top 10 proteins for random forest:")
s2

## correlation method
corr_data <- train_bio %>% 
  filter(!is.na(ados))

s3 <- corr_data %>%
  pivot_longer(cols = -c(ados, group),
               names_to = 'protein',
               values_to = 'level') %>%
  group_by(protein) %>%
  summarize(correlation = cor(as.numeric(ados), level)) %>%
  arrange(desc(correlation)) %>%
  slice_head(n = 10)  # Select top 10 proteins

print("Top 10 proteins for correlation method:")
s3$protein


proteins_sstar <- intersect(s1, s2)

proteins_sstar

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

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