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

# Count outliers for each subject
outlier_summary_2 <- biomarker_data %>%
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

library(glmnet)
install.packages("glmnet")
library(tidymodels)

# read in data
biomarker_clean <- read_csv('data/biomarker-raw.csv')
# partition
set.seed(101622)
partitions <- biomarker_clean %>%
  initial_split(prop = 0.8)

x_train <- training(partitions) %>%
  as.matrix()
y_train <- training(partitions)



# reproducibility
set.seed(102022)

# multiple partitioning for lambda selection
cv_out <- cv.glmnet(x_train, 
                    y_train, 
                    family = 'binomial', 
                    nfolds = 5, 
                    type.measure = 'deviance')

cvout_df <- tidy(cv_out) 















