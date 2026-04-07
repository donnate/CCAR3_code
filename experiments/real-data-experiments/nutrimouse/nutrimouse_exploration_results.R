library(tidyverse)

df <- read_csv("/Users/clairedonnat/Documents/CCAR3_code/experiments/real-data-experiments/nutrimouse/nutrimouse_benchmark_all_methods_cv/benchmark_runs.csv")


df_summ <- df %>%
  group_by(repeat_id,method ) %>%
  summarise_if(is.numeric, mean, na.rm=TRUE) %>%
  group_by(method ) %>%
  summarise(mean_test_mse = mean(test_mse),
            sd_test_mse = sd(test_mse),
            mean_test_cor = mean(test_cor),
            sd_test_cor = sd(test_cor),
            mean_time = mean(runtime_sec),
            sd_time = sd(runtime_sec))


df2 <- read_csv("/Users/clairedonnat/Documents/CCAR3_code/experiments/real-data-experiments/nutrimouse/nutrimouse_benchmark_all_methods_scaled_with_fantope/benchmark_runs.csv")


df_summ2 <- df2 %>%
  group_by(repeat_id,method ) %>%
  summarise_if(is.numeric, median, na.rm=TRUE) %>%
  group_by(method ) %>%
  summarise(mean_test_mse = mean(test_mse),
            sd_test_mse = sd(test_mse),
            mean_test_cor = mean(test_cor),
            sd_test_cor = sd(test_cor),
            mean_time = median(runtime_sec),
            sd_time = sd(runtime_sec))
