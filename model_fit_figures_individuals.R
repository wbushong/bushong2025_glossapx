######################################
## Comparisons of models + plots of model fits
## within each individual
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
d.combined$unique.subjectID <- with(d.combined, paste(subject, explabel, sep = "_"))
subjs <- unique(d.combined$unique.subjectID)

## gather all the model files together
dir.prefix <- c("models_individuals/")
model.prefixes <- c("ambiguity/ambiguity_model_all_ranefs_subject", "catdiscard/catdiscard_model_all_ranefs_subject",
                    "catswitch/catswitch_new_model_all_ranefs_subject", "ideal/idealobserver_model_all_ranefs_subject",
                    "contextonly/contextonly_model_all_ranefs_subject")
model.files <- as.vector(sapply(model.prefixes, function(.) paste0(dir.prefix, ., "_", subjs, ".RDS")))

## ---- Compute pairwise elpd-waic comparisons ----
## this is the basis of Table 3 in the main text
best.models <- c()
all.model.comps <- c()
ind.removes <- c()
for (i in 1:length(subjs)) {
  # load the models
  expmodels <- model.files[grep(paste0("_", subjs[i]), model.files)]
  models <- list()
  conv <- c()
  for (j in 1:length(expmodels)) {
    models[[j]] <- readRDS(expmodels[j])
    conv.stats <- c(summary(models[[j]])$random$sFrame$Rhat, summary(models[[j]])$fixed$Rhat)
    conv <- c(conv, conv.stats)
    models[[j]] <- add_criterion(models[[j]], "waic", file = gsub(".RDS", "", expmodels[j]), force_save = TRUE)
  }
  # if models didn't converge, move to next subj & add to list of nonconvergence
  conv <- all(conv < 1.01) & all(conv > .99)
  if (conv == FALSE) {
    ind.removes <- c(ind.removes, subjs[i])
    next
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
      comparison$dif_se_scaled <- with(comparison, abs(elpd_diff / se_diff))
      comparison$subject <- subjs[i]
      model.comps <- rbind(model.comps, comparison)
    }
  }
  all.model.comps <- rbind(all.model.comps, model.comps)
  # find best-fitting model & compare to second-best-fitting model
  # each model comparison consists of 2 rows: 
  # better-fitting model's row w/ an elpd_diff of 0 & worse-fitting model's row has the actual elpd difference
  # so the best-fitting overall model will have all elpd_diff's of 0
  best.fit.model <- model.comps %>% group_by(model) %>% mutate(best_model = all(elpd_diff == 0))
  best.fit.model.info <- data.frame(subject = subjs[i],
                                    best.model = unique(best.fit.model$model[best.fit.model$best_model==TRUE]),
                                    lowest.elpd.se.diff = min(model.comps$dif_se_scaled,na.rm=TRUE))
  best.models <- rbind(best.models, best.fit.model.info)
}
saveRDS(best.models, file = "model_comparisons/individual_model_comparisons_summarized.RDS")
saveRDS(all.model.comps, file = "model_comparisons/individual_model_comparisons_all.RDS")

# descriptive statistics
best.models$weak.evidence <- best.models$lowest.elpd.se.diff > 2.5
best.models$strong.evidence <- best.models$lowest.elpd.se.diff > 5

## ---- Visualize model comparisons across individuals  ----
## Corresponds to Figure 7 in the text
fullnames <- c("categorize-&-discard", "ambiguity-dependent",
               "categorize-discard-&-switch", "ideal integration", "context-only")

extract_comps <- function(df) {
  comparison <- strsplit(df$comp[1], split="_")[[1]][c(1,3)]
  df$comparison <- ifelse(df$model==comparison[1], comparison[2], comparison[1])
  df$elpd_diff <- ifelse(df$elpd_diff==0, abs(min(df$elpd_diff)), df$elpd_diff)
  df$se_diff <- ifelse(df$se_diff==0, max(df$se_diff), df$se_diff)
  df$dif_se_scaled <- with(df, elpd_diff / se_diff)
  return(df[,c("model", "comparison", "dif_se_scaled")])
}

pairwise.model.comps <- c()
subjs2 <- unique(all.model.comps$subject)
for (i in 1:length(subjs2)) {
  curr.comp <- subset(all.model.comps, subject == subjs2[i])
  new.comps <- c()
  for (j in unique(curr.comp$comp)) {
    new.comps <- rbind(new.comps, extract_comps(subset(curr.comp,comp==j)))
  }
  new.comps$subject <- subjs[i]
  pairwise.model.comps <- rbind(pairwise.model.comps, new.comps)
}

pairwise.model.comps$model <- factor(pairwise.model.comps$model,
                                     levels = c("ideal", "ambiguity",
                                                "catdiscard", "catswitch",
                                                "contextonly"),
                                     labels = fullnames[c(4, 2, 1, 3, 5)])
pairwise.model.comps$comparison <- factor(pairwise.model.comps$comparison,
                                     levels = rev(c("ideal", "ambiguity",
                                                "catdiscard", "catswitch",
                                                "contextonly")),
                                     labels = rev(fullnames[c(4, 2, 1, 3, 5)]))

p.all_comparisons <- ggplot(pairwise.model.comps, aes(x = comparison,
                                                      y = dif_se_scaled)) +
  geom_point(alpha = .01) +
  ylab("Difference in model fit\n(elpd difference scaled by standard error)") +
  xlab("") +
  coord_cartesian(ylim = c(0, 8), xlim = c(0.5, 5)) +
  coord_flip() +
  theme(legend.position = "none") +
  annotate("rect", xmin = 0, xmax = 6, ymin = -2.5, ymax = 2.5,
           alpha = .1,fill = "black") +
  annotate("rect", xmin = 0, xmax = 6, ymin = 2.5, ymax = 5,
           alpha = .2,fill = "green3") +
  annotate("rect", xmin = 0, xmax = 6, ymin = 5, ymax = 20,
           alpha = .35,fill = "green3") +
  annotate("rect", xmin = 0, xmax = 6, ymin = -5, ymax = -2.5,
           alpha = .2,fill = "red3") +
  annotate("rect", xmin = 0, xmax = 6, ymin = -20, ymax = -5,
           alpha = .35,fill = "red3") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(alpha = .6) +
  # stat_summary(geom = "pointrange", fun.data = mean_cl_boot, color = "red", size = .5) +
  facet_wrap(~ model) +
  theme(strip.background = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(color = "black"))
ggsave("figures/Figure7.pdf",
       p.all_comparisons,
       height = 5,
       width = 10)
