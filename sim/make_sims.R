library(tidyverse)
library(mrgsolve)
library(mapbayr) 
stopifnot(packageVersion("mapbayr") >= "0.10.0")

source("sim/funs.R")

my_patients <- tibble(
  # fake dataset of n = 5 patients, 
  # should be a table of covariate representative of the population
  BW = c(32.6, 57.2, 60.1, 13.1, 16.5), 
  BSA = c(1.06, 1.6, 1.44, 0.475, 0.675), 
  AGE = c(15.3, 19.8, 7.69, 0.762, 5.79),
  DIAG = c(0, 1, 0, 1, 1), 
  GSTA1 = c(0, 0, 1, 1, 0), 
  SEX = c(1, 0, 1, 0, 1)
)

# We do it once. Only dosing according to the genetics.
treatment(
  target_AUC = 16000,
  method1 = "genet_equation",
  method2 = "genet_equation",
  method3 = "genet_equation", 
  method4 = "genet_equation", 
  patient = slice_sample(my_patients, n = 1000, replace = TRUE) # sampling with replacement
)$cumAUC %>% hist

# Table of the simulated strategies
to_test <- tibble(
  Day = 1:4,
  NoGenet_eq = c("nogenet_equation", "nogenet_equation","nogenet_equation", "nogenet_equation"),
  Genet_eq   = c("genet_equation",   "genet_equation",  "genet_equation",   "genet_equation"),
  TDM_day1   = c("nogenet_equation", "TDM_bayes_1",     "TDM_bayes_1",      "TDM_bayes_1"), 
  TDM_day13  = c("nogenet_equation", "TDM_bayes_1",     "TDM_bayes_1",      "TDM_bayes_3"), 
  TDM_day13_7samp  = c("nogenet_equation", "TDM_bayes_1",     "TDM_bayes_1",      "TDM_bayes_3"), 
  TDM_day123 = c("nogenet_equation", "TDM_bayes_1",     "TDM_bayes_2",      "TDM_bayes_3")
) %>% 
  pivot_longer(-Day) %>% 
  pivot_wider(names_from = "Day", names_prefix = "method") %>% 
  mutate(sampletimes = if_else(str_detect(name, "7samp"), list(protocol_times), list(c(3.25, 5, 8))))

to_test

# Computing (can be quite long)
ans <- to_test[,-1] %>% 
    pmap(.f = treatment,
         target = 16000,
         patient = slice_sample(my_patients, n = 100, replace = TRUE)
         # only n = 100 to spare time but should be more 
    )

# Matching scenarios and results
ans_tab <- ans %>% 
  set_names(to_test$name) %>% 
  bind_rows(.id = "scenar") %>% 
  mutate(scenar = factor(scenar, levels = unique(scenar)))

# Summary
ans_tab %>% 
  group_by(scenar) %>% 
  summarise(
    Q05 = quantile(cumAUC/16000, 0.05) %>% scales::percent(), 
    Q95 = quantile(cumAUC/16000, 0.95) %>% scales::percent(), 
  ) %>% 
  mutate(CI90 = paste0(Q05, " - ", Q95), .keep = "unused")

# Plot
ans_tab %>% 
  ggplot(aes(cumAUC/16000)) + 
  geom_histogram(aes(fill = scenar),col = "black") +
  facet_grid(~scenar) +
  coord_flip() + 
  scale_x_log10() +
  theme(legend.position = "none")
