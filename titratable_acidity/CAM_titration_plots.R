###CAM plots####################################################################
## Anri Chomentowska: last updated March 2025

rm(list = ls()) # clear environment
setwd("~/Dropbox/Yale/Research/CAM") #change your environment
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)


############Load and wrangle data
cam <- read.csv("CAMtritrationsAllExp.csv")
cam2 <- cam %>%
  dplyr::select(PlantID,
                Day,
                Time,
                Treatment,
                State,
                Daylength,
                H,
                Experiment) %>%
  dplyr::filter(!is.na(H))%>%
  pivot_wider(id_cols = c(PlantID, Day, Treatment, State, Daylength, Experiment), names_from = Time , values_from = H)%>%
  mutate(deltaH = AM - PM)

View(cam)
View(cam2)

cam2 %>%
  group_by(State) %>% 
  summarise(number = length(PlantID))

#re-order how I want it
cam2$State <- factor(cam2$State, levels = c("Vegetative","Reproductive","Flowering","Fruiting"))
cam2$Treatment <- factor(cam2$Treatment, levels = c("Watered","Drought"))

#filter out the NAs and extremely negative values
cam3 <- cam2 %>%
  dplyr::filter(!is.na(deltaH))%>%
  dplyr::filter(!is.na(State))

#optionally
write.csv(cam3, "CAMtitrationsALL_wrangled.csv", row.names=F)

#Assign the categories
cam_new <- cam3 %>%
  mutate(New_State = case_when(
    State %in% c("Vegetative", "Reproductive") ~ "Pre_Flowering",
    State %in% c("Flowering", "Fruiting") ~ "Post_Flowering",
    TRUE ~ as.character(State)
  ))

#re-order again
cam_new$New_State <- factor(cam_new$New_State, levels = c("Pre_Flowering","Post_Flowering"))
cam_new$Treatment <- factor(cam_new$Treatment, levels = c("Watered","Drought"))

#split data into experiments, and filter out huge titration outliers (values of less than -10)
#if a value like that exists in your data, probably human error
cam_exp1 <- cam_new %>%
  dplyr::filter(New_State!="Post_Flowering") %>%
  dplyr::filter(Daylength=="Short") %>%
  dplyr::filter(deltaH > -10)

cam_exp2 <- cam_new %>%
  dplyr::filter(Daylength=="Long") %>%
  dplyr::filter(deltaH > -10)

#wrangle data so that it can be plotted on one plot
cam_exp1_new <- cam_exp1 %>%
  mutate(New_State = Treatment)

cam_all <- cam_exp1_new %>%
  bind_rows(cam_exp2)

View(cam_all)


#############Plot box plot
p <- ggplot(cam_all, aes(x = New_State, y = deltaH, fill = New_State))
p + geom_boxplot(alpha=0.5, outlier.shape = 1, linewidth = 1, aes(colour = New_State)) +
  geom_jitter(alpha=0.7, width = 0.1, size = 3, aes(colour = New_State)) + 
  scale_fill_manual(values=c("#619CFF", "#F8766D", "seagreen3", "magenta2")) +
  scale_color_manual(values=c("#619CFF", "#F8766D", "seagreen3", "magenta2")) +
  scale_x_discrete(labels=c("Watered" = "Well-Watered", "Drought" = "Drought", "Pre_Flowering" = "Pre-Flowering", "Post_Flowering" = "Post-Flowering")) +
  scale_y_continuous(limits = c(-10, 43)) +
  ylab("Titratable acidity (∆H µmol/g FW)" ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.position = "none") +
  geom_vline(xintercept = 2.5, linetype = "dotted", color = "black", linewidth = 0.5)


#############Stats: 
####One-sided t-test to compare between treatments/states
res <- wilcox.test(deltaH ~ Treatment, data = cam_exp1,
                   alternative = "less")
res #shows watered samples have significantly smaller deltaH value

res <- wilcox.test(deltaH ~ New_State, data = cam_exp2,
                   alternative = "less")
res #shows pre-flowering samples have slightly (and somewhat significant) smaller deltaH value


####One-sample t-test if means are significantly > than zero for each bar (==CAM)
cam_water <- cam_exp1 %>%
  dplyr::filter(Treatment!="Drought")
res <- wilcox.test(cam_water$deltaH, mu = 0)
res 

cam_drought <- cam_exp1 %>%
  dplyr::filter(Treatment!="Watered")
res <- wilcox.test(cam_drought$deltaH, mu = 0)
res

cam_pre <- cam_exp2 %>%
  dplyr::filter(New_State!="Post_Flowering")
res <- wilcox.test(cam_pre$deltaH, mu = 0)
res

cam_post <- cam_exp2 %>%
  dplyr::filter(New_State!="Pre_Flowering")
res <- wilcox.test(cam_post$deltaH, mu = 0)
res
