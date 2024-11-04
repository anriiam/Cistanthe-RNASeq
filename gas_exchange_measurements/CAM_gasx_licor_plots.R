###LICOR plots##################################################################
rm(list = ls()) # clear environment
setwd("~/Dropbox/Yale/Research/CAM/Licor") # set your 
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(ggtext)
library(RColorBrewer)
set.seed(11)

### reads data
licor <- read.csv("Spotcheck_cropped_adjusted.csv")

### filter data
licor2 <- licor %>%
  dplyr::select(obs,
                samples,
                elapsed,
                A,
                gsw,
                A.SD) %>%
  dplyr::filter(A > -10) %>%
  dplyr::filter(samples!="veg_old") %>%
  dplyr::filter(samples!="fruit_old")

head(licor2)

### plot co2 assimilation rate

p <- ggplot(licor2, aes(x = elapsed, y = A, colour = samples))
p + geom_point(alpha=0.9, stat="identity") + 
  geom_errorbar(aes(ymin=A-A.SD, ymax=A+A.SD), width=.2, size = 2, position=position_dodge(0.05)) + 
  geom_line(size = 2) +
  ylab("Assimilation (µmol m-2 s-1)" ) +
  scale_color_manual(labels = c("Drought", "Well-Watered"), values=c("#F8766D", "#619CFF"), name = "Treatment") +
  xlab("Time elapsed (s)" ) +
  geom_hline(yintercept=0, linetype="dashed") +
  annotate("rect",
           alpha = 0.3,
           fill = "lightgrey",
           xmin = 0, xmax =36000,
           ymin = -Inf, ymax = Inf) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "right")
  
### plot stomatal conductance
p1 <- ggplot(licor2, aes(x = obs, y = gsw))
p1 + geom_bar(alpha=0.9, stat="identity", color = 'steelblue4')

ggplot(licor2, aes(x = obs)) +
  geom_bar(aes(y = gsw),
           fill = 'steelblue4',
           stat="identity",
           alpha=0.8)+
  geom_bar(aes(y = A / 100),
           fill = 'magenta4',
           stat="identity",
           alpha=0.4)+
  scale_y_continuous(name = "Stomatal conductance (gsw)",
                     sec.axis = sec_axis( trans=~.*100, name="CO2 assimilation (A)")
  ) 


