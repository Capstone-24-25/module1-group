library(tidyverse)


# get names
var_names <- read_csv('data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

# function for trimming outliers (good idea??)
trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}

# read in data
biomarker_clean <- read_csv('data/biomarker-raw.csv', 
                            skip = 2,
                            col_select = -2L,
                            col_names = c('group', 
                                          'empty',
                                          pull(var_names, abbreviation),
                                          'ados'),
                            na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # log transform, center and scale, and trim
  mutate(across(.cols = -c(group, ados), 
                ~ trim(scale(log10(.x))[, 1], .at = 3))) %>%
  # reorder columns
  select(group, ados, everything())

# THIS IS THE ACTUAL CODE BELOW, THE ABOVE IS JUST SETUP

set.seed(333) # Set seed for reproducibility
sample_proteins <- sample(names(biomarker_clean)[-c(1, ncol(biomarker_clean))], 5)  # Exclude 'group' and 'ados'

# Pivot data to a long format for easier plotting
biomarker_long <- biomarker_clean %>%
  select(all_of(sample_proteins)) %>%
  pivot_longer(cols = everything(), names_to = "protein", values_to = "value")

# Plot histograms of raw values for the sampled proteins
ggplot(biomarker_long, aes(x = value)) +
  geom_histogram(bins = 30, color = "black", fill = "skyblue", alpha = 0.7) +
  facet_wrap(~ protein, scales = "free") +
  labs(title = "Distribution of Raw Protein Levels",
       x = "Protein Level",
       y = "Frequency") +
  theme_minimal()

# The data seems to be skewed to the right, which would be helped by a log transform. 

library(glmnet)
set.seed(234)


biomarker_LASSO = biomarker_clean %>% 
  select(-ados) %>% 
  mutate(class = as.numeric(group == 'ASD'))
partitions <- biomarker_LASSO %>%
  initial_split(prop = 0.8)

x_train <- training(partitions) %>%
  select(-group, -class) %>%
  as.matrix()
y_train <- training(partitions) %>%
  pull(class)

plot(cv.glmnet(x_train, y_train, family = 'binomial', nfolds = 5, type.measure = 'deviance'))

lambda = exp(-1.8)

lambda_test <- cv.glmnet(x_train, y_train, family = 'binomial', nfolds=5)
lambda_min <- lambda_test$lambda.1se
lambda_min = exp(-3.5)
final_fit <- glmnet(x_train, y_train, family = 'binomial', lambda = lambda_min)
final_fit_df = tidy(final_fit)


proteins_panel = final_fit_df %>%
  filter(term != "(Intercept)") %>% 
  pull(term)

biomarker_panel <- biomarker_clean %>%
  select(group, any_of(proteins_panel)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_panel

biomarker_split <- biomarker_panel %>%
  initial_split(prop = 0.8)

fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

class_metrics <- metric_set(accuracy)

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')

