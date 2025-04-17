########################################
## Demonstrate the qualitative model predictions
## Corresponds to Figures 4-5 in main text & Figure S1 in SI
########################################
## ---- Load libraries & make plot theme ----
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)

theme_clean <- theme_bw() + 
  theme(strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        text = element_text(size = 15))

## ---- Specify models ----
# ideal integration
compute_logodds_ideal <- function(p_t_VOT, p_d_VOT, p_t_context, p_d_context) {
  logodds_t <- log(p_t_VOT * p_t_context) - log(p_d_VOT * p_d_context)
  return(logodds_t)
}

# ambiguity
compute_logodds_ambiguity <- function(p_t_VOT, p_d_VOT, p_t_context, p_d_context, p_context_VOT) {
  lambda <- abs(p_t_VOT - 0.5) * 2
  logodds_t <- log(lambda * p_t_VOT + (1-lambda) * (p_t_VOT * p_t_context)/p_context_VOT) - 
    log(lambda * p_d_VOT + (1-lambda) * (p_d_VOT * p_d_context)/p_context_VOT)
  return(logodds_t)
}

# categorize-&-discard
compute_logodds_catdiscard <- function(p_t_VOT, p_d_VOT, p_t_context, p_d_context) {
  logodds_t <- log(p_t_VOT / p_d_VOT)
  return(logodds_t)
}

# categorize-discard-&-switch
compute_logodds_catswitch <- function(p_t_VOT, p_d_VOT, p_t_context, p_d_context) {
  logodds_t <- log(p_t_VOT + (1 - p_t_VOT) * p_t_context) - log(p_d_VOT + (1 - p_d_VOT) * p_d_context)
  return(logodds_t)
}

# context-only
compute_logodds_co <- function(p_t_context, p_d_context) {
  logodds_t <- log(p_t_context) - log(p_d_context)
  return(logodds_t)
}

## ---- Calculate model predictions ----
# set up template of VOT & context conditions
VOT <- seq(10, 85, 0.01) # same VOT bounds as in the experiments
context <- c("tent", "dent") # two subsequent context conditions

# grid of VOT-context pairings
# convert VOT to scaled polynomial
m.template <- expand.grid(VOT, context)
names(m.template) <- c("VOT", "context")
m.template$vot.1 <- poly(m.template$VOT, degree = 2)[, 1]
m.template$vot.2 <- poly(m.template$VOT, degree = 2)[, 2]
m.template$p_context_VOT <- 0.5

## give some cue weights for purposes of demonstration
# set p(t|VOT) & p(t|context)
# show context effects of 1 and 3 in log-odds space (divided by 2 b/c of how fns work)
# i.e., log(p(t|context)) - log(p(d|context)) = 1 (or 3) 
# in our experiments we typically find context effects between 1-3 so these are plausible points to visualize
ce_lo <- c(.5, 1.5) 
context_effects <- exp(ce_lo) / (1 + exp(ce_lo)) # probability space

# look at 2 vot effect sizes (this often changes qual. results in the models)
# first -- difference of 4 in log-odds between min & max VOT (small VOT effect)
# second -- difference of 8 in log-odds between min & max VOT (large VOT effect)
vot_lo_1 <- c(2 / max(m.template$vot.1), 4 / max(m.template$vot.1))

# don't visualize quadratic vot effects
vot_lo_2 <- 0

# expand out conditions
conditions <- expand.grid(vot_lo_1, vot_lo_2, context_effects)
names(conditions) <- c("vot_lo_1", "vot_lo_2", "c_p")

assign_probabilities <- function(template, vot_lo_1, vot_lo_2, c_p) {
  curr_df <- template
  lo_t_VOT <- curr_df$vot.1 * vot_lo_1 + curr_df$vot.2 * vot_lo_2
  curr_df$p_t_VOT <- exp(lo_t_VOT) / (1 + exp(lo_t_VOT))
  curr_df$p_d_VOT <- rev(curr_df$p_t_VOT)
  curr_df$p_t_context <- c(rep(c_p, length(unique(curr_df$VOT))), rep(1-c_p, length(unique(curr_df$VOT))))
  curr_df$p_d_context <- c(rep(1-c_p, length(unique(curr_df$VOT))), rep(c_p, length(unique(curr_df$VOT))))
  return(curr_df)
}

# compute model predictions
all_predictions <- c()
for (i in 1:nrow(conditions)) {
  curr_df <- assign_probabilities(m.template, conditions$vot_lo_1[i], conditions$vot_lo_2[i], conditions$c_p[i])
  curr_df$ideal_observer <- with(curr_df, compute_logodds_ideal(p_t_VOT, p_d_VOT, p_t_context, p_d_context))
  curr_df$ambiguity <- with(curr_df, compute_logodds_ambiguity(p_t_VOT, p_d_VOT, p_t_context, p_d_context, p_context_VOT))
  curr_df$cat_discard <- with(curr_df, compute_logodds_catdiscard(p_t_VOT, p_d_VOT, p_t_context, p_d_context))
  curr_df$cat_switch <- with(curr_df, compute_logodds_catswitch(p_t_VOT, p_d_VOT, p_t_context, p_d_context))
  curr_df$context_only <- with(curr_df, compute_logodds_co(p_t_context, p_d_context))
  models <- c("ideal_observer", "cat_discard", "cat_switch", "ambiguity")
  curr_df$vot_lo_1 <- conditions$vot_lo_1[i]
  curr_df$vot_lo_2 <- conditions$vot_lo_2[i]
  curr_df$c_lo <- log(conditions$c_p[i] / (1 - conditions$c_p[i]))
  all_predictions <- rbind(all_predictions, curr_df)
}

# label VOT effects as "small" or "large"
all_predictions$vot_lo_1.fac <- factor(all_predictions$vot_lo_1, levels = unique(all_predictions$vot_lo_1),
                                     labels = c("Small VOT effect", 
                                                "Large VOT effect"))
all_predictions$vot_lo_2.fac <- factor(all_predictions$vot_lo_2)
 
# show context effects (subtract dent-biasing from tent-biasing condition) 
all_predictions_agg <- all_predictions %>%
  group_by(VOT, vot_lo_1.fac, c_lo, vot_lo_2.fac) %>%
  summarise_at(vars(ideal_observer:context_only), function(x) x[context=="tent"] - x[context=="dent"])


## ---- VOT & context effects (Figure S1) ----
# make predictions for each individual plot
p.ambiguity.unagg <- ggplot(all_predictions, aes(x = VOT, y = ambiguity, color = context,
                                                 group = interaction(context,c_lo),
                                                 linetype = factor(c_lo))) + 
  geom_line(size = 1) +
  scale_color_manual(values = c("blue", "red"), name = "Context bias",
                     labels = c('t-biasing (e.g., "?ent...forest")', 'd-biasing (e.g., "?ent...fender")')) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Context effect",
                        labels = c("p(t|t-biasing context) = .75",
                                   "p(t|t-biasing context) = .95")) +
  ylab("Predicted log-odds of\n/t/ response") +
  xlab("") +
  ylim(c(-5.5, 5.5)) +
  facet_wrap(~ vot_lo_1.fac) +
  theme_clean

p.cat_switch.unagg <- ggplot(all_predictions, aes(x = VOT, y = cat_switch, color = context,
                                                 group = interaction(context,c_lo),
                                                 linetype = factor(c_lo))) + 
  geom_line(size = 1) +
  scale_color_manual(values = c("blue", "red"), name = "Context bias",
                     labels = c('t-biasing (e.g., "?ent...forest")', 'd-biasing (e.g., "?ent...fender")')) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Context effect",
                        labels = c("p(t|t-biasing context) = .75",
                                   "p(t|t-biasing context) = .95")) +
  ylab("Predicted log-odds of\n/t/ response") +
  xlab("") +
  ylim(c(-5.5, 5.5)) +
  facet_wrap(~ vot_lo_1.fac) +
  theme_clean

p.idealobserver.unagg <- ggplot(all_predictions, aes(x = VOT, y = ideal_observer, color = context,
                                                     group = interaction(context,c_lo),
                                                     linetype = factor(c_lo))) + 
  geom_line(size = 1) +
  scale_color_manual(values = c("blue", "red"), name = "Context bias",
                     labels = c('t-biasing (e.g., "?ent...forest")', 'd-biasing (e.g., "?ent...fender")')) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Context effect",
                        labels = c("p(t|t-biasing context) = .75",
                                   "p(t|t-biasing context) = .95")) +
  ylab("Predicted log-odds of\n/t/ response") +
  xlab("") +
  ylim(c(-5.5, 5.5)) +
  facet_wrap(~ vot_lo_1.fac) +
  theme_clean

p.discard.unagg <- ggplot(all_predictions, aes(x = VOT, y = cat_discard, color = context,
                                               group = interaction(context,c_lo),
                                               linetype = factor(c_lo))) + 
  geom_line(size = 1) +
  scale_color_manual(values = c("blue", "red"), name = "Context bias",
                     labels = c('t-biasing (e.g., "?ent...forest")', 'd-biasing (e.g., "?ent...fender")')) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Context effect",
                        labels = c("p(t|t-biasing context) = .75",
                                   "p(t|t-biasing context) = .95")) +
  ylab("Predicted log-odds of\n/t/ response") +
  xlab("") +
  ylim(c(-5.5, 5.5)) +
  facet_wrap(~ vot_lo_1.fac) +
  theme_clean

p.co.unagg <- ggplot(all_predictions, aes(x = VOT, y = context_only, color = context,
                                               group = interaction(context,c_lo),
                                               linetype = factor(c_lo))) + 
  geom_line(size = 1) +
  scale_color_manual(values = c("blue", "red"), name = "Context bias",
                     labels = c('t-biasing (e.g., "?ent...forest")', 'd-biasing (e.g., "?ent...fender")')) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Context effect",
                        labels = c("p(t|t-biasing context) = .75",
                                   "p(t|t-biasing context) = .95")) +
  ylab("Predicted log-odds of\n/t/ response") +
  xlab("") +
  ylim(c(-5.5, 5.5)) +
  facet_wrap(~ vot_lo_1.fac) +
  theme_clean

## now, context effect only
p.ambiguity <- ggplot(all_predictions_agg, aes(x = VOT, y = ambiguity, 
        linetype = factor(c_lo))) + 
  geom_line(size = 1.3) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Context effect",
                        labels = c("p(t|t-biasing context) = .75",
                                   "p(t|t-biasing context) = .95")) +
  ylab("Predicted context\neffect (log-odds)") +
  ylim(c(0, 3.5)) +
  xlab("VOT (ms)") +
  facet_wrap(~ vot_lo_1.fac) +
  theme_clean +
  theme(strip.text = element_blank())

p.cat_switch <- ggplot(all_predictions_agg, aes(x = VOT, y = cat_switch, 
                                                linetype = factor(c_lo))) + 
  geom_line(size = 1.3) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Context effect",
                        labels = c("p(t|t-biasing context) = .75",
                                   "p(t|t-biasing context) = .95")) +
  ylab("Predicted context effect\n(log-odds)") +
  ylim(c(0, 3.5)) +
  xlab("VOT (ms)") +
  facet_wrap(~ vot_lo_1.fac) +
  theme_clean +
  theme(strip.text = element_blank())

p.idealobserver <- ggplot(all_predictions_agg, aes(x = VOT, y = ideal_observer, 
                                                   linetype = factor(c_lo))) + 
  geom_line(size = 1.3) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Context effect",
                        labels = c("p(t|t-biasing context) = .75",
                                   "p(t|t-biasing context) = .95")) +
  ylab("Predicted context effect\n(log-odds)") +
  ylim(c(0, 3.5)) +
  xlab("VOT (ms)") +
  facet_wrap(~ vot_lo_1.fac) +
  theme_clean +
  theme(strip.text = element_blank())

p.discard <- ggplot(all_predictions_agg, aes(x = VOT, y = cat_discard, 
                                             linetype = factor(c_lo))) + 
  geom_line(size = 1.3) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Context effect",
                        labels = c("p(t|t-biasing context) = .75",
                                   "p(t|t-biasing context) = .95")) +
  ylab("Predicted context effect\n(log-odds)") +
  ylim(c(0, 3.5)) +
  xlab("VOT (ms)") +
  facet_wrap(~ vot_lo_1.fac) +
  theme_clean +
  theme(strip.text = element_blank())

p.co <- ggplot(all_predictions_agg, aes(x = VOT, y = context_only, 
                                             linetype = factor(c_lo))) + 
  geom_line(size = 1.3) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Context effect",
                        labels = c("p(t|t-biasing context) = .75",
                                   "p(t|t-biasing context) = .95")) +
  ylab("Predicted context effect\n(log-odds)") +
  ylim(c(0, 3.5)) +
  xlab("VOT (ms)") +
  facet_wrap(~ vot_lo_1.fac) +
  theme_clean +
  theme(strip.text = element_blank())

## plot grid
shared_legend <- get_legend(p.ambiguity.unagg)

pa1 <- p.ambiguity.unagg + theme(legend.position="none")
pa2 <- p.ambiguity + theme(legend.position="none")
pi1 <- p.idealobserver.unagg + theme(legend.position="none") + ylab("")
pi2 <- p.idealobserver + theme(legend.position="none") + ylab("")
pcd1 <- p.discard.unagg + theme(legend.position="none") + ylab("")
pcd2 <- p.discard + theme(legend.position="none") + ylab("")
pcs1 <- p.cat_switch.unagg + theme(legend.position="none") + ylab("")
pcs2 <- p.cat_switch + theme(legend.position="none") + ylab("")
pco1 <- p.co.unagg + theme(legend.position="none") + ylab("")
pco2 <- p.co + theme(legend.position="none") + ylab("")

p.models.combined <- cowplot::plot_grid(
  cowplot::plot_grid(
    pa1, pi1, pcd1, pcs1, pco1,
    pa2, pi2, pcd2, pcs2, pco2,
    nrow = 2, ncol = 5, labels = c(letters[1:10])
  ),
  cowplot::plot_grid(NULL, shared_legend, NULL, ncol = 1, rel_widths = c(0.5, 1, 0.5)),
  ncol = 2,
  rel_widths = c(1, .2)
)
ggsave("figures/FigureS1.pdf",
       p.models.combined,
       width = 25,
       height = 6,
       units = "in")

## ---- Simple individual model plots for main text (Figure 4 in main text) ----
p.ambiguity.maintext <- ggplot(subset(all_predictions, vot_lo_1.fac == "Large VOT effect" & c_lo == 1.5), 
                               aes(x = VOT, y = ambiguity, color = context, group = context)) + 
  geom_line(size = 1) +
  scale_color_manual(values = c("blue", "red"), name = "Context bias",
                     labels = c('t-biasing (e.g., "?ent...forest")', 'd-biasing (e.g., "?ent...fender")')) +
  ylab("") +
  xlab("") +
  ylim(c(-5.5, 5.5)) +
  theme_clean
p.ideal.maintext <- ggplot(subset(all_predictions, vot_lo_1.fac == "Large VOT effect" & c_lo == 1.5), 
                               aes(x = VOT, y = ideal_observer, color = context, group = context)) + 
  geom_line(size = 1) +
  scale_color_manual(values = c("blue", "red"), name = "Context bias",
                     labels = c('t-biasing (e.g., "?ent...forest")', 'd-biasing (e.g., "?ent...fender")')) +
  ylab("") +
  xlab("") +
  ylim(c(-5.5, 5.5)) +
  theme_clean
p.catdiscard.maintext <- ggplot(subset(all_predictions, vot_lo_1.fac == "Large VOT effect" & c_lo == 1.5), 
                               aes(x = VOT, y = cat_discard, color = context, group = context)) + 
  geom_line(size = 1) +
  scale_color_manual(values = c("blue", "red"), name = "Context bias",
                     labels = c('t-biasing (e.g., "?ent...forest")', 'd-biasing (e.g., "?ent...fender")')) +
  ylab("") +
  xlab("") +
  ylim(c(-5.5, 5.5)) +
  theme_clean
p.catswitch.maintext <- ggplot(subset(all_predictions, vot_lo_1.fac == "Large VOT effect" & c_lo == 1.5), 
                               aes(x = VOT, y = cat_switch, color = context, group = context)) + 
  geom_line(size = 1) +
  scale_color_manual(values = c("blue", "red"), name = "Context bias",
                     labels = c('t-biasing (e.g., "?ent...forest")', 'd-biasing (e.g., "?ent...fender")')) +
  ylab("") +
  xlab("") +
  ylim(c(-5.5, 5.5)) +
  theme_clean
p.co.maintext <- ggplot(subset(all_predictions, vot_lo_1.fac == "Large VOT effect" & c_lo == 1.5), 
                               aes(x = VOT, y = context_only, color = context, group = context)) + 
  geom_line(size = 1) +
  scale_color_manual(values = c("blue", "red"), name = "Context bias",
                     labels = c('t-biasing (e.g., "?ent...forest")', 'd-biasing (e.g., "?ent...fender")')) +
  ylab("") +
  xlab("") +
  ylim(c(-5.5, 5.5)) +
  theme_clean

shared_legend_2 <- get_legend(p.ambiguity.maintext)
pa <- p.ambiguity.maintext + theme(legend.position="none", text = element_text(size=10)) + ggtitle("ambiguity-dependent") + ylab("")
pi <- p.ideal.maintext + theme(legend.position="none", text = element_text(size=10)) + ggtitle("ideal integration")
pcd <- p.catdiscard.maintext + theme(legend.position="none", text = element_text(size=10)) + ggtitle("categorize-&-discard")
pcs <- p.catswitch.maintext + theme(legend.position="none", text = element_text(size=10)) + ylab("") + ggtitle("categorize-discard-&-switch")
pco <- p.co.maintext + theme(legend.position="none", text = element_text(size=10)) + ylab("") + ggtitle("context-only")

p.mods <- cowplot::plot_grid(
    pi, pa, pcd, pcs, pco, shared_legend_2,
    nrow = 2, ncol = 3, labels = c(letters[1:5])
  )
y.grob2 <- textGrob("Predicted /t/-responses (log-odds)", 
                   gp=gpar(fontsize=12), rot=90)
x.grob2 <- textGrob("VOT (ms)", 
                   gp=gpar(fontsize=12))
p.mods2 <- grid.arrange(arrangeGrob(p.mods, left = y.grob2, bottom = x.grob2))
ggsave("figures/Figure4.pdf",
       p.mods2,
       width = 9,
       height = 5.5,
       units = "in")


## ---- Plot all models on same grid (Figure 5 in main text) ----
vot1s <- unique(all_predictions_agg$vot_lo_1.fac)
vot2s <- unique(all_predictions_agg$vot_lo_2.fac)
all_predictions_agg2 <- all_predictions_agg %>%
  pivot_longer(cols = ideal_observer:context_only,
               names_to = "model",
               values_to = "log_odds")
p.allpreds.context <- ggplot(subset(all_predictions_agg2, 
                            vot_lo_1.fac == vot1s[2] & 
                            vot_lo_2.fac == vot2s[1] &
                            c_lo == 1.5),
                     aes(x = VOT, y = log_odds, color = model)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = c("#d95f02", "#e7298a", "#7570b3", 
                                "#1b9e77", "#66A61E"),
                     labels = c("ambiguity-dependent", "categorize-&-discard",
                                "categorize-discard-&-switch", "context-only", "ideal integration"),
                     name = "Model") +
  ylab("Predicted context\neffect (log-odds)") +
  xlab("") +
  theme_clean 
  # theme(legend.position = "left")

all_predictions_agg3 <- all_predictions %>%
  group_by(VOT, vot_lo_1.fac, c_lo, vot_lo_2.fac) %>%
  summarise_at(vars(ideal_observer:context_only), mean)
all_predictions_agg3 <- all_predictions_agg3 %>%
  pivot_longer(cols = ideal_observer:context_only,
               names_to = "model",
               values_to = "log_odds")
p.allpreds.vot <- ggplot(subset(all_predictions_agg3, 
                            vot_lo_1.fac == vot1s[2] & 
                              vot_lo_2.fac == vot2s[1] &
                              c_lo == 1.5),
                     aes(x = VOT, y = log_odds, color = model)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = c("#d95f02", "#e7298a", "#7570b3", 
                                "#1b9e77", "#66A61E"),
                     labels = c("ambiguity-dependent", "categorize-&-discard",
                                "categorize-discard-&-switch", "context-only", "ideal integration")) +
  ylab("Predicted /t/-responses\n(log-odds)") +
  xlab("") +
  theme_clean +
  theme(legend.position = "none")

shared_legend_individ <- get_legend(p.allpreds.context)
p.allpreds.context <- p.allpreds.context + theme(legend.position = "none")

x.grob3 <- textGrob("VOT (ms)", 
                   gp=gpar(fontsize=15))
p.noleg <- cowplot::plot_grid(p.allpreds.vot, p.allpreds.context, labels = letters[1:2])
p.noleg2 <- grid.arrange(arrangeGrob(p.noleg, bottom = x.grob3))

p.votctxt.individ.preds <- cowplot::plot_grid(
  p.noleg2,
  shared_legend_individ,
  ncol = 2,
  rel_widths = c(1, .5)
)

ggsave("figures/Figure5.pdf",
       p.votctxt.individ.preds,
       width = 8,
       height = 4,
       units = "in")
