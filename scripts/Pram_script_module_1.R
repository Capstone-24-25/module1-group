#Pram's script

getwd()
setwd("C:/Users/pramu/OneDrive/Documents/GitHub/module1-group")
getwd()


library(tidyverse)
library(ggplot2)

# Read variable names as before
var_names <- read_csv('data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

# Read in data without trimming
biomarker_clean_no_trim <- read_csv('data/biomarker-raw.csv', 
                                    skip = 2,
                                    col_select = -2L,
                                    col_names = c('group', 
                                                  'empty',
                                                  pull(var_names, abbreviation),
                                                  'ados'),
                                    na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # log transform, center and scale without trimming
  mutate(across(.cols = -c(group, ados), 
                ~ scale(log10(.x))[, 1])) %>%
  # reorder columns
  select(group, ados, everything())







# Define threshold for identifying outliers
outlier_threshold <- 3

# Count outliers per subject
outlier_counts <- biomarker_clean_no_trim %>%
  rowwise() %>%
  mutate(outlier_count = sum(abs(c_across(-c(group, ados))) > outlier_threshold, na.rm = TRUE)) %>%
  ungroup()

############ EDA

# Display the summary of outliers per group
outlier_summary <- outlier_counts %>%
  group_by(group) %>%
  summarize(mean_outliers = mean(outlier_count),
            max_outliers = max(outlier_count),
            subjects_with_outliers = sum(outlier_count > 0))

# See that the average amount of outliers for ASD is 13.2 and 17.6 for TD







biomarker_data_outliers <- read_csv('data/biomarker-raw.csv', 
                           skip = 2,
                           col_select = -2L,
                           col_names = c('group', 'empty', pull(var_names, abbreviation), 'ados'),
                           na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  mutate(across(-c(group, ados), ~ abs(scale(log10(.x))[, 1]) > 3, .names = "outlier_{.col}"))









# Count outliers for each subject
outlier_summary_2 <- biomarker_data_outliers %>%
  mutate(total_outliers = rowSums(select(., starts_with("outlier_")))) %>%
  select(group, total_outliers)

# Display the table
outlier_summary_2
View(outlier_summary_2)


#Every subject had at least 1 outlier except for 2
#5 subjects had a dispropportionate amount of outliers (>100)
## 2 of these subjects were ASD and 4 were TD
#Outliers did not appear more frequent in one group in regards to another

# Create the scatter plot
ggplot(outlier_summary_2, aes(x = group, y = total_outliers, color = group)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
  labs(
    title = "Distribution of Outliers by Group",
    x = "Group",
    y = "Total Outliers"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("ASD" = "skyblue", "TD" = "salmon"))






#4

library(tidyverse)
library(tidymodels)
library(glmnet)
library(broom)

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

################# Code starts here
set.seed(42239)
biomarker_LASSO <- biomarker_clean %>% 
  select(-ados) %>% 
  mutate(class = as.numeric(group == 'ASD'))

partitions <- biomarker_LASSO %>%
  initial_split(prop = 0.8)

x_train <- training(partitions) %>%
  select(-group, -class) %>%
  as.matrix()
y_train <- training(partitions) %>%
  pull(class)

#lambda_test <- cv.glmnet(x_train, y_train, family = 'binomial', nfolds = 10) # Increased nfolds for more robust CV
#lambda_min <- lambda_test$lambda.min  # Use lambda.min to ensure lower bias if accuracy improves


#### For less proteins
lambda_test <- cv.glmnet(x_train, y_train, family = 'binomial', nfolds = 5)
lambda_min <- lambda_test$lambda.1se
# **Adjust lambda to control number of proteins in the panel (Used 15)
lambda_adjusted <- exp(-1.95) # 


final_fit <- glmnet(x_train, y_train, family = 'binomial', lambda = lambda_min)
final_fit_df <- tidy(final_fit)

proteins_panel <- final_fit_df %>%
  filter(term != "(Intercept)") %>% 
  pull(term)

biomarker_panel <- biomarker_clean %>%
  select(group, any_of(proteins_panel)) %>%
  mutate(class = as.factor(group == 'ASD')) %>%
  select(-group)

biomarker_split <- biomarker_panel %>%
  initial_split(prop = 0.8)

fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

class_metrics <- metric_set(accuracy, roc_auc)

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')


#Accuracy is 0.871

