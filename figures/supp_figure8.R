library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)

df_vfa_t48 <- read_csv("data/chemical_growth_kinetic/VFA_T48.csv", col_types = cols(.default = "c"))
df_vfa_t48$run <- "T48"

df_vfa_t49 <- read_csv("data/chemical_growth_kinetic/VFA_T49.csv", col_types = cols(.default = "c"))
df_vfa_t49$run <- "T49"

df_vfa_exp1 <- bind_rows(df_vfa_t48, df_vfa_t49) %>%
  
  gather(compound, concentration, -reactor, -run, -time) %>%
  mutate(concentration = as.numeric(concentration)) %>%
  drop_na()

df_meta_exp1 <- read_csv("data/meta/meta_combined_exp1.csv", col_types = cols(.default = "c"))
df_time <- read_csv("data/meta/hours_to_day.csv", col_types = cols(.default = "c"))

df_vfa_summary <- df_vfa_exp1 %>%
  
  left_join(df_meta_exp1, by = c("run", "reactor")) %>%
  left_join(df_time, by = c("time")) %>%
  
  mutate(concentration = concentration / 20) %>%
  
  group_by(run, day, treatment, compound) %>%
  mutate(median_concentration = median(concentration)) %>%
  ungroup() %>%
  select(day, treatment, compound, run, median_concentration) %>%
  distinct() %>%
  drop_na() %>%
  
  filter(compound != "Isobutyric_acid") 


df_vfa_summary$day <- as.integer(df_vfa_summary$day)

df_vfa_summary$compound <- factor(df_vfa_summary$compound, levels = c("Acetic_acid", "Propionic_acid", "Butyric_acid", "Isovaleric_acid", "Valeric_acid"))

ggplot(df_vfa_summary, aes(x = day, y = median_concentration, colour = treatment)) +
  theme_bw() +
  geom_point(aes(colour = treatment)) +
  facet_grid(compound ~ run, scales = "free_y") +
  geom_smooth(aes(colour = treatment), se = FALSE) +
  ylab("Median Concentration (mM)") +
  xlab("Time Point (day)")

####

df_mcfa_t49 <- read_csv("data/chemical_growth_kinetic/MCFA_T49.csv", col_types = cols(.default = "c"))

df_meta <- read_csv("data/meta/meta_combined_exp1.csv", col_types = cols(.default = "c")) %>%
  filter(run == "T49")

df_mcfa_conv <- read_csv("data/chemical_growth_kinetic/MCFA_T49_conversion.csv", col_types = cols(.default = "c"))

df_mcfa_summary <- df_mcfa_t49 %>%
  left_join(df_meta, by = c("reactor")) %>%
  
  gather(compound, concentration, -reactor, -treatment, -run) %>%
  
  inner_join(df_mcfa_conv) %>%
  
  mutate(concentration = as.numeric(concentration)) %>%
  mutate(divide_by = as.numeric(divide_by)) %>%
  mutate(concentration = concentration / divide_by) %>%
  
  group_by(treatment, compound) %>%
  mutate(median_concentration = median(concentration)) %>%
  mutate(sd_concentration = sd(concentration)) %>%
  ungroup() %>%
  select(treatment, compound, median_concentration, sd_concentration, concentration) %>%
  distinct() %>%
  
  mutate(upper = median_concentration + sd_concentration) %>%
  mutate(lower = median_concentration - sd_concentration)

df_mcfa_summary <- df_mcfa_summary %>%
  filter(compound != "2-Methylbutyric") %>%
  filter(compound != "Isobutyric")

df_mcfa_summary$compound <- factor(df_mcfa_summary$compound, levels = c("Acetic", "Propionic", "Butyric", "Isovaleric", "Valeric", "Hexanoic", "Heptanoic"))

ggplot(df_mcfa_summary, aes(x = treatment, y = median_concentration, fill = treatment)) +
  
  theme_bw() +
  geom_bar(stat = "identity", position = "dodge") +
  geom_jitter(aes(y = concentration), width = 0.1, height = 0, alpha = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  
  facet_grid(compound~., scales = "free", ) +
  ylab("Median Concentration (mM)") +
  xlab("Treatment") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))