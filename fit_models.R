##################################
## Fit models to datasets
##################################
## --------------------- Libraries and constants ----
library(tidyverse)
library(brms)
library(bayesplot)
library(tidybayes)
library(rstan)
library(lme4)

# Stan settings
max.cores <- min(parallel::detectCores(), 4)
options(mc.cores = max.cores)

## --------------------- Functions ----
source("preprocess_functions.R")

inv_logit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}

## --------------------- Load data ----
exp.files <- paste0("data/exp", 1:4, ".RDS")
d.combined <- preprocess_and_combine(exp.files, "acousticPoints.csv")

## --------------------- Misc. statistics ----
## pull out the empirically observed min/max proportion of /t/-responses --
## how close to ceiling/floor are they?
d.by.subj <- d.combined %>%
  group_by(explabel, subject, VOT) %>%
  summarise(respond_t = mean(respond_t)) %>%
  group_by(explabel, VOT) %>%
  summarise(mean_t = mean(respond_t))

## --------------------- Set priors ----
# priors for ideal integration & ambiguity models
priors <- c(prior(student_t(3, 0, 2.5), class = "b", nlpar = "bvot"),
            prior(student_t(3, 0, 2.5), class = "b", nlpar = "bcontext"),
            prior(cauchy(0,5), class = "sd", nlpar = "bvot"),
            prior(cauchy(0,5), class = "sd", nlpar = "bcontext"))
# priors for categorize-&-discard models
priors.io_discard <- c(prior(student_t(3, 0, 2.5), class = "b"),
            prior(cauchy(0,5), class = "sd"))

## --------------------- Fit models ----
## First, fit to experiment 1
d.combined$trial.log0 <- log(d.combined$Trial+1)
d.exp1 <- subset(d.combined, explabel == "Experiment 1")

# categorize-discard-&-switch model
catswitch.formula <- bf(respond_t ~ log(inv_logit(bvot) + ((1-inv_logit(bvot)) * inv_logit(bcontext)))
                        - log((1-inv_logit(bvot)) + (inv_logit(bvot)) * (1-inv_logit(bcontext))), 
                        bvot ~ 1 + vot.1 + vot.2 + (vot.1 + vot.2 | g | subject) + (vot.1 + vot.2 | g2 | sFrame),
                        bcontext ~ 0 + bias.numeric + (0 + bias.numeric | g | subject) + (0 + bias.numeric | g2 | sFrame),
                        nl = TRUE)
m.catswitch <- brm(catswitch.formula,
                   data = d.exp1,
                   family = bernoulli(link="logit"),
                   prior = priors,
                   chains = 4,
                   thin = 4,
                   warmup = 2000,
                   iter = 6000, 
                   control = list(adapt_delta = .9999, max_treedepth = 15),
                   file = "./models/catswitch_new_model_all_ranefs_Experiment 1.RDS") 

# ambiguity model
ambiguity.formula <- bf(respond_t ~ 
                              log((     fabs(inv_logit(bvot) - 0.5) * 2) *       inv_logit(bvot)  + ((1 - (fabs(inv_logit(bvot) - 0.5) * 2)) *     inv_logit(bvot + bcontext))) -
                              log((     fabs(inv_logit(bvot) - 0.5) * 2) * (1 - inv_logit(bvot)) +       ((1 - (fabs(inv_logit(bvot) - 0.5) * 2)) * (1 - inv_logit(bvot + bcontext)))),
                            bvot ~ 1 + vot.1 + vot.2 + (vot.1 + vot.2 | g | subject) + (vot.1 + vot.2 | g2 | sFrame), # this is essentially p(t|VOT) in the ambiguity model formula in paper
                            bcontext ~ 0 + bias.numeric + (0 + bias.numeric | g | subject) + (0 + bias.numeric | g2 | sFrame), # this is p(t|context) in the ambiguity model formula in paper
                            nl = TRUE)
m.ambiguity <- brm(ambiguity.formula,
                   data = d.exp1,
                   family = bernoulli(link="logit"),
                   prior = priors,
                   chains = 4,
                   thin = 4,
                   warmup = 2000,
                   iter = 6000, 
                   control = list(adapt_delta = .9999, max_treedepth = 15),
                   file = "./models/ambiguity_model_all_ranefs_Experiment 1.RDS") 
                   

# ideal integration
io.formula <- bf(respond_t ~ vot.1 + vot.2 + bias.numeric +
                   (vot.1 + vot.2 + bias.numeric | subject) +
                   (vot.1 + vot.2 + bias.numeric | sFrame))
m.io <- brm(io.formula,
                   data = d.exp1,
                   family = bernoulli(link="logit"),
                   prior = priors.io_discard,
                   chains = 4,
                   thin = 1,
                   warmup = 2000,
                   iter = 6000, 
                   control = list(adapt_delta = .995, max_treedepth = 10),
                   file = "./models/idealobserver_model_all_ranefs_Experiment 1.RDS") 

# vanilla categorize-&-discard model
discard.formula <- bf(respond_t ~ vot.1 + vot.2 +
                        (vot.1 + vot.2 | subject) +
                        (vot.1 + vot.2 | sFrame))
m.discard <- brm(discard.formula,
            data = d.exp1,
            family = bernoulli(link="logit"),
            prior = priors.io_discard,
            chains = 4,
            thin = 1,
            warmup = 2000,
            iter = 6000, 
            control = list(adapt_delta = .9, max_treedepth = 10),
            refresh = 100,
            file = "./models/catdiscard_model_all_ranefs_Experiment 1.RDS") 

# update 2/5/25: new context-only model
contextonly.formula <- bf(respond_t ~ bias.numeric +
                        (bias.numeric | subject) +
                        (bias.numeric | sFrame))
m.contextonly <- brm(contextonly.formula,
                 data = d.exp1,
                 family = bernoulli(link="logit"),
                 prior = priors.io_discard,
                 chains = 4,
                 thin = 1,
                 warmup = 2000,
                 iter = 6000, 
                 control = list(adapt_delta = .9, max_treedepth = 10),
                 refresh = 100,
                 file = "./models/contextonly_model_all_ranefs_Experiment 1.RDS") 

## fit brms w/ trial 
# for simplicity, exclude VOT^2 (convergence issues)
trial.formula <- bf(respond_t ~ vot.1 * trial.log0 + 
                      bias.numeric * trial.log0 +
                   (vot.1 + bias.numeric + trial.log0 | subject) +
                   (vot.1 + bias.numeric + trial.log0 | sFrame))
m.trial <- brm(trial.formula,
            data = d.exp1,
            family = bernoulli(link="logit"),
            prior = priors.io_discard,
            chains = 4,
            thin = 1,
            warmup = 2000,
            iter = 6000, 
            control = list(adapt_delta = .999, max_treedepth = 10),
            file = "./models/trial_model_all_ranefs_Experiment 1.RDS") 

## now, loop through other datasets using update() to avoid recompiling
for (i in paste("Experiment", 2:4)) {
  d.current <- subset(d.combined, explabel == i)
  update(m.discard,
         newdata = d.current,
         recompile = FALSE,
         file = paste0("./models/catdiscard_model_all_ranefs_", i, ".RDS"))
  update(m.catswitch,
         newdata = d.current,
         recompile = FALSE,
         file = paste0("./models/catswitch_new_model_all_ranefs_", i, ".RDS"))
  update(m.io,
         newdata = d.current,
         recompile = FALSE, control = list(adapt_delta = .9999),
         file = paste0("./models/idealobserver_model_all_ranefs_", i, ".RDS"))
  update(m.ambiguity,
         newdata = d.current,
         recompile = FALSE, control = list(adapt_delta = .9999),
         file = paste0("./models/ambiguity_model_all_ranefs_", i, ".RDS"))
  update(m.contextonly,
         newdata = d.current,
         recompile = FALSE, control = list(adapt_delta = .9999),
         file = paste0("./models/contextonly_model_all_ranefs_", i, ".RDS"))
  update(m.trial,
                newdata = d.current,
                recompile = FALSE, control = list(adapt_delta = .9999),
         file = paste0("./models/trial_model_all_ranefs_", i, ".RDS"))
}
