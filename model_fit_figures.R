######################################
## Comparisons of models + plots of model fits
######################################
## ---- Libraries + other preliminaries ----
library(tidyverse)
library(brms)
library(cowplot)
library(tidybayes)
theme_set(theme_cowplot() + 
            theme(panel.grid.major = element_line(color = "grey", size = 0.2),
                  panel.grid.minor = element_line(color = "grey90", size = 0.2)))

source("preprocess_functions.R")

inv_logit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}

## load datasets
exp.files <- paste0("data/exp", 1:4, ".RDS")
d.combined <- preprocess_and_combine(exp.files, "acousticPoints.csv")

## gather all the model files together
dir.prefix <- c("models/")
model.prefixes <- c("ambiguity_model_all_ranefs", "catdiscard_model_all_ranefs",
                    "catswitch_model_all_ranefs", "idealobserver_model_all_ranefs",
                    "contextonly_model_all_ranefs")
exps <- paste("Experiment", c(1:4))
model.files <- as.vector(sapply(model.prefixes, function(.) paste0(dir.prefix, ., "_", exps, ".RDS")))

## ---- Compute pairwise elpd-waic comparisons ----
## this is the basis of Table 2 in the main text
for (i in 1:4) {
  expmodels <- model.files[grep(exps[i], model.files)]
  models <- list()
  for (j in 1:length(expmodels)) {
    models[[j]] <- readRDS(expmodels[j])
    models[[j]] <- add_criterion(models[[j]], "waic", file = gsub(".RDS", "", expmodels[j]), force_save = TRUE)
  }
  names(models) <- c("ambiguity", "catdiscard", "catswitch", "ideal", "contextonly")
  model.comps <- c()
  # pairwise model comparisons
  for (j in 1:(length(models)-1)) {
    # for each model go 1 further
    for (k in (j+1):length(models)) {
      comparison <- as.data.frame(loo_compare(models[[j]], models[[k]], criterion = "waic")) %>%
        rownames_to_column("model") %>%
        mutate(comp = paste0(names(models)[j], "_vs_", names(models)[k]))
      comparison$model <- ifelse(comparison$model == "models[[j]]",
                                 names(models)[[j]], names(models)[[k]])
      model.comps <- rbind(model.comps, comparison)
    }
  }
  saveRDS(model.comps, file = paste0("model_comparisons/WAICS_and_comparisons_Experiment", i, ".RDS"))
}

## ---- Visualize predictions for each model against empirical data  ----
## in order of best -> worst fit top -> bottom
## Corresponds to Figure 6 in the text + Figures S3-6 in SI
for (i in 1:4) {
  plots = list()
  expmodels <- model.files[grep(exps[i], model.files)]
  models <- list()
  for (j in 1:length(expmodels)) {
    models[[j]] <- readRDS(expmodels[j])
  }
  names(models) <- c("ambiguity", "catdiscard", "catswitch", "ideal", "contextonly")
  fullnames <- c("ambiguity-dependent", "categorize-&-discard", "categorize-discard-&-switch",
                 "ideal integration", "context-only")
  model.order <- c(4, 1, 2:3, 5)
  
  for (m in 1:length(models)) {
    n.samples <- 1000 # increase for accuracy
    model <- models[[m]]
    
    # Get polys based on experiment, as well as their standard deviations 
    # (used below to scale the data for the visualization of fits to the original predictors used during model fitting)
    p = with(subset(d.combined, explabel == exps[i]), poly(VOT, 2))
    sd.1 = with(subset(d.combined, explabel == exps[i]), sd(poly(VOT, 2)[, 1]))
    sd.2 = with(subset(d.combined, explabel == exps[i]), sd(poly(VOT, 2)[, 2]))
    
    # Add predictions (slow, if uncertainty about random effects is considered)
    # Choose bias.numeric / bias.fac, depending on model
    d <- 
      crossing(
        bias.numeric = c(-.5, .5),
        VOT = 10:85) %>%
      mutate(
        context = ifelse(bias.numeric %in% c(-.5, "dent"), "dent", "tent"), # written to be general, so that you switch out bias.fac/numeric
        vot.1.unstandardized = predict(p, VOT)[, 1],
        vot.2.unstandardized = predict(p, VOT)[, 2],
        vot.1 = vot.1.unstandardized / (2 * sd.1),
        vot.2 = vot.2.unstandardized / (2 * sd.2)
      ) %>%
      # Extract predictions in logits for increased accuracy
      # scale = response or linear (in this case, they are identical)
      # re_formula = NULL to add uncertainty about random effects (necessary also for by-subject plots)
      add_fitted_draws(model, n = n.samples, scale = "linear", re_formula = NA) %>% 
      distinct()
    
    # Static plots
    p.vot.ctxt <- 
      d %>%
      ggplot(aes(x = VOT, y = plogis(.value))) +
      # first call creates ribbon
      stat_lineribbon(data = ~ filter(.x, context == "dent"), geom = "ribbon", alpha = .25, fill = "red",
                      point_interval = mean_hdci, .width = c(.95)) +
      stat_lineribbon(data = ~ filter(.x, context == "tent"), geom = "ribbon", alpha = .25, fill = "blue",
                      point_interval = mean_hdci, .width = c(.95)) +
      # second call creates lines
      stat_lineribbon(data = ~ filter(.x, context == "dent"), geom = "line", color = "red",
                      point_interval = mean_hdci, linetype = "dashed", size = 1) +
      stat_lineribbon(data = ~ filter(.x, context == "tent"), geom = "line", color = "blue",
                      point_interval = mean_hdci, linetype = "dashed", size = 1) +
      stat_summary(
        data = subset(d.combined,explabel==exps[i]) %>% 
          group_by(subject, VOT, context) %>%
          summarise(respond_t = mean(respond_t)),
        aes(y = respond_t, color = context), 
        fun.data = mean_cl_boot, 
        geom = "pointrange", size = .5) +
      scale_x_continuous("VOT (ms)") +
      scale_y_continuous("Proportion /t/\nresponses", limits = c(0, 1)) +
      scale_color_manual("Context",
                         breaks = c("tent", "dent"),
                         values = c("blue", "red"),
                         labels = c("tent-biasing", "dent-biasing")) +
      scale_linetype_discrete("Context", guide = FALSE) +
    theme_bw() + theme(legend.position = "bottom")
    ggsave(plot = p.vot.ctxt,
           file = paste0("model_fit_figures/fitted_noRandomEffects_prop_preds_", names(models)[m], "_Experiment", i, ".pdf"),
           height = 5,
           width = 5,
           unit = "in",
           device = "pdf")
    
    p.vot.ctxt.logit <- 
      d %>%
      ggplot(aes(x = VOT, y = .value)) +
      # first call creates ribbon
      stat_lineribbon(data = ~ filter(.x, context == "dent"), geom = "ribbon", alpha = .25, fill = "red",
                      point_interval = mean_hdci, .width = c(.95)) +
      stat_lineribbon(data = ~ filter(.x, context == "tent"), geom = "ribbon", alpha = .25, fill = "blue",
                      point_interval = mean_hdci, .width = c(.95)) +
      # second call creates lines
      stat_lineribbon(data = ~ filter(.x, context == "dent"), geom = "line", color = "red",
                      point_interval = mean_hdci, linetype = "dashed", size = 1) +
      stat_lineribbon(data = ~ filter(.x, context == "tent"), geom = "line", color = "blue",
                      point_interval = mean_hdci, linetype = "dashed", size = 1) +
      scale_x_continuous("VOT (ms)") +
      scale_y_continuous("Log-odds of /t/ responses") +
      scale_color_manual("Context",
                         breaks = c("tent", "dent"),
                         values = c("blue", "red")) +
      scale_linetype_discrete("Context", guide = FALSE) + 
      theme_bw() + theme(legend.position = "bottom")
    ggsave(plot = p.vot.ctxt.logit,
           file = paste0("model_fit_figures/fitted_noRandomEffects_logits_preds_", names(models)[m], "_Experiment", i, ".pdf"),
           height = 5,
           width = 5,
           unit = "in",
           device = "pdf")
    # if categorize-&-discard, manually draw context effect line at 0 b/c I'm too lazy to make it work with the estimates with axis limits etc
    if (m == 2) { 
      p.ctxt <- 
        d %>%
        arrange(.draw, VOT, context) %>%
        group_by(.draw, VOT) %>%
        summarise(.value = last(.value) - first(.value)) %>%
        ggplot(aes(x = VOT, y = .value)) +
        # just make a line at 0
        geom_segment(aes(x = 10, xend = 85, y = 0, yend = 0), color = "black", size = 1) +
      scale_x_continuous("VOT (ms)") +
        scale_y_continuous("Context effect (log-odds)") +
        theme_bw() + theme(legend.position = "bottom")
      
      p.ctxt = p.ctxt + 
        scale_y_continuous("Context effect (log-odds)", breaks = 0:4) +
        coord_cartesian(ylim = c(-.1, 3))  
    } else {    
      p.ctxt <- 
        d %>%
        arrange(.draw, VOT, context) %>%
        group_by(.draw, VOT) %>%
        summarise(.value = last(.value) - first(.value)) %>%
        ggplot(aes(x = VOT, y = .value)) +
        # first call creates ribbon
        stat_lineribbon(geom = "ribbon", alpha = .25, fill = "black",
                        point_interval = mean_hdci, .width = c(.95)) +
        # second call creates lines
        stat_lineribbon(geom = "line", color = "black",
                        point_interval = mean_hdci, size = 1) +
      scale_x_continuous("VOT (ms)") +
        scale_y_continuous("Context effect (log-odds)") +
        theme_bw() + theme(legend.position = "bottom")
      
      p.ctxt = p.ctxt + 
        scale_y_continuous("Context effect (log-odds)", breaks = 0:4) +
        coord_cartesian(ylim = c(-.1, 3))
    }
    ggsave(plot = p.ctxt,
           file = paste0("model_fit_figures/fitted_noRandomEffects_ctxt_logits_preds_", names(models)[m], "_Experiment", i, "_ConstantLimits.pdf"),
           height = 5,
           width = 5,
           unit = "in",
           device = "pdf")
    
    rm(d)
    plots[[m]] = cowplot::plot_grid(
      plotlist = 
        list(NULL,
             cowplot::plot_grid(
               plotlist = list(p.vot.ctxt + 
                                 { if (m != length(names(models))) theme(axis.title.x = element_blank(), axis.text.x = element_blank()) } +
                                 theme(legend.position = "none"), 
                               p.vot.ctxt.logit + 
                                 { if (m != length(names(models))) theme(axis.title.x = element_blank(), axis.text.x = element_blank()) } +
                                 theme(legend.position = "none"), 
                               p.ctxt + 
                                 { if (m != length(names(models))) theme(axis.title.x = element_blank(), axis.text.x = element_blank()) } +
                                 theme(legend.position = "none")),
               align = "hv", axis = "btrl", nrow = 1)),
      rel_heights = c(.05, .95), ncol = 1)
  }
  plots[[m+1]] = get_legend(p.vot.ctxt)
  p.combined = cowplot::plot_grid(
    plotlist = plots[c(model.order, length(model.order) + 1)],
    ncol = 1, hjust = 0, vjust = 1.3, label_size = 10,
    labels = c(paste0(" ", letters[1:m], "    ", fullnames[model.order]), ""),
    rel_heights = c(rep(1 / length(names(models)), length(names(models))), 1/30) +
      c(rep(0, length(names(models)) - 1), 1/35, 0))
  
  # combined plot for Experiment 2 is fig 6 in paper
  # combined plots for all experiments are figs S3-6 in SI
  # name accordingly
  ggsave(plot = p.combined,
         file = paste0("figures/FigureS", i+1, ".pdf"),
         height = 2.5 * length(names(models)) + .75,
         width = 8,
         unit = "in",
         device = "pdf")
  if (i == 2) system("cp figures/FigureS3.pdf figures/Figure6.pdf")
  
  rm(p.combined, p.vot.ctxt, p.vot.ctxt.logit, p.ctxt)
}

