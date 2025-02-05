library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)

df_myco_t48 <- read_csv("data/chemical_growth_kinetic/myco_T48.csv", col_types = cols(.default = "c"))
df_myco_t48$run <- "T48"

df_myco_t49 <- read_csv("data/chemical_growth_kinetic/myco_T49.csv", col_types = cols(.default = "c"))
df_myco_t49$run <- "T49"

df_myco_exp1 <- bind_rows(df_myco_t48, df_myco_t49)

df_meta_exp1 <- read_csv("data/meta/meta_combined_exp1.csv", col_types = cols(.default = "c"))

df_time <- read_csv("data/meta/hours_to_day.csv", col_types = cols(.default = "c"))

df_myco_comb <- left_join(df_myco_exp1, df_meta_exp1, by = c("run", "reactor")) %>%
    inner_join(df_time, by = "time")

df_myco_summary <- df_myco_comb %>%
    
    gather(compound, concentration, -time, -run, -reactor, -treatment, -day) %>%
    mutate(concentration = as.numeric(concentration)) %>%

    group_by(time, treatment, compound, run) %>%
    mutate(median_concentration = median(concentration)) %>%
    ungroup() %>%
    
    select(day, treatment, compound, run, median_concentration) %>%
    distinct()

df_select_plot <- df_myco_summary %>%
    filter(compound != "DON") %>%
    filter(compound != "ZEN")

df_select_plot$compound <- factor(df_select_plot$compound, levels = c("DOM", "aZEL", "bZEL"))
df_select_plot$day <- as.integer(df_select_plot$day)

ggplot(df_select_plot, aes(x = day, y = median_concentration, colour = treatment)) +
    theme_bw() +
    geom_point(aes(colour = treatment)) +
    facet_wrap(run ~ compound, scales = "free_y") +
    geom_smooth(aes(colour = treatment), se = FALSE) +
    ylab("Median Concentration (uM)") +
    xlab("Time Point (day)")
