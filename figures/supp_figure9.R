library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

mycotoxin_inhibition_df <- read_csv("data/chemical_growth_kinetic/mycotoxin_inhibition.csv") %>%
  
  mutate(time_hrs = round(time / 60 / 60)) %>%
  mutate(od_blank = od - blank) %>%
  
  group_by(treatment, time_hrs) %>%
  mutate(med_od_blank = median(od_blank)) %>%
  mutate(sd_od_blank = sd(od_blank)) %>%
  ungroup() %>%
  
  select(treatment, time_hrs, med_od_blank, sd_od_blank) %>%
  distinct() %>%
  
  filter(time_hrs <= 30) %>%
  
  mutate(group = if_else(treatment == "DMSO" | treatment == "AUR", "dmso", "ethanol"))

mycotoxin_inhibition_df$group <- factor(mycotoxin_inhibition_df$group, levels = c("ethanol", "dmso"))

ggplot(mycotoxin_inhibition_df, aes(x=time_hrs, y=med_od_blank, color = treatment)) +
  geom_errorbar(aes(ymin=med_od_blank-sd_od_blank, ymax=med_od_blank+sd_od_blank), width=0.2,  alpha=0.5, size=0.5) +  
  geom_point(stat="identity", fill="skyblue", size = 3) + 
  theme_minimal() + 
  facet_wrap(. ~ group, scales = "free_y") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) 