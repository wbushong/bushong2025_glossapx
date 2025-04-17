############################
## Analyze norming study for Experiments 3-4
############################

############################
## Preliminaries
############################

library(tidyverse)
library(lme4)

theme_clean <- theme_bw() + 
  theme(strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        text = element_text(size = 15))

############################
## Load data
############################
load("data/norming_data_preprocessed.RData")

d$respond_t <- ifelse(d$respWord == "tent", 1, 0)
d$VOT <- as.numeric(as.character(d$VOT))

############################
## Plot
############################
d.by.subj <- d %>%
	group_by(subject, VOT) %>%
	summarise(respond_t = mean(respond_t))

p.avg <- ggplot(d.by.subj, aes(x = VOT, y = respond_t)) +
	stat_summary(fun.data = mean_cl_boot, geom = "pointrange") +
  theme_clean +
  ylab("Proportion /t/ responses") +
  xlab("VOT (ms)")
ggsave("figures/FigureS6.pdf",
       p.avg,
       height = 4,
       width = 4)

############################
## Statistics
############################
# Gelman-scale VOT
gs <- function(x) (x - mean(x)) / (sd(x)*2)

d$vot.scaled <- gs(d$VOT)

# fit GLMM
m <- glmer(respond_t ~ vot.scaled + (1 + vot.scaled | subject), 
	d, family = "binomial")

# find category boundary
cat_boundary <- (- summary(m)$coef[1,1]) / summary(m)$coef[2,1]
cat_boundary_VOTspace <- (cat_boundary * 2*(sd(d$VOT))) + mean(d$VOT)

# empirical ceiling & floor of responses
propresps <- d.by.subj %>%
  group_by(VOT) %>%
  summarise(mean_resp = mean(respond_t))


