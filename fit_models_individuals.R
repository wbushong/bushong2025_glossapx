##################################
## Fit models to individuals
##################################
## --------------------- Libraries and constants ----
library(tidyverse)
library(brms)
library(bayesplot)
library(tidybayes)
library(rstan)

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
## Fit to first subject
d.combined$unique.subjectID <- with(d.combined, paste(subject, explabel, sep = "_"))
subjs <- unique(d.combined$unique.subjectID)

# categorize-discard-&-switch model
catswitch.formula <- bf(respond_t ~ log(inv_logit(bvot) + ((1-inv_logit(bvot)) * inv_logit(bcontext)))
                        - log((1-inv_logit(bvot)) + (inv_logit(bvot)) * (1-inv_logit(bcontext))), 
                        bvot ~ 1 + vot.1 + vot.2 + (vot.1 + vot.2 | g2 | sFrame),
                        bcontext ~ 0 + bias.numeric + (0 + bias.numeric | g2 | sFrame),
                        nl = TRUE)
m.catswitch <- brm(catswitch.formula,
                   data = subset(d.combined, unique.subjectID == subjs[1]),
                   family = bernoulli(link="logit"),
                   prior = priors,
                   chains = 4,
                   thin = 4,
                   warmup = 2000,
                   iter = 6000, 
                   control = list(adapt_delta = .9999, max_treedepth = 15),
                   file = paste0("./models_individuals/catswitch/catswitch_new_model_all_ranefs_subject_", subjs[1], ".RDS")) 

# ambiguity model
ambiguity.formula <- bf(respond_t ~ 
                              log((     fabs(inv_logit(bvot) - 0.5) * 2) *       inv_logit(bvot)  + ((1 - (fabs(inv_logit(bvot) - 0.5) * 2)) *     inv_logit(bvot + bcontext))) -
                              log((     fabs(inv_logit(bvot) - 0.5) * 2) * (1 - inv_logit(bvot)) +       ((1 - (fabs(inv_logit(bvot) - 0.5) * 2)) * (1 - inv_logit(bvot + bcontext)))),
                            bvot ~ 1 + vot.1 + vot.2 + (vot.1 + vot.2 | g2 | sFrame), # this is essentially p(t|VOT) in the ambiguity model formula in paper
                            bcontext ~ 0 + bias.numeric + (0 + bias.numeric | g2 | sFrame), # this is p(t|context) in the ambiguity model formula in paper
                            nl = TRUE)
m.ambiguity <- brm(ambiguity.formula,
                   data = subset(d.combined, unique.subjectID == subjs[1]),
                   family = bernoulli(link="logit"),
                   prior = priors,
                   chains = 4,
                   thin = 4,
                   warmup = 2000,
                   iter = 6000, 
                   control = list(adapt_delta = .9999, max_treedepth = 15),
                   file = paste0("./models_individuals/ambiguity/ambiguity_model_all_ranefs_subject_", subjs[1], ".RDS")) 
                   

# ideal integration
io.formula <- bf(respond_t ~ vot.1 + vot.2 + bias.numeric +
                   (vot.1 + vot.2 + bias.numeric | sFrame))
m.io <- brm(io.formula,
                   data = subset(d.combined, unique.subjectID == subjs[1]),
                   family = bernoulli(link="logit"),
                   prior = priors.io_discard,
                   chains = 4,
                   thin = 1,
                   warmup = 2000,
                   iter = 6000, 
                   control = list(adapt_delta = .995, max_treedepth = 10),
                   file = paste0("./models_individuals/ideal/idealobserver_model_all_ranefs_subject_", subjs[1], ".RDS")) 

# vanilla categorize-&-discard model
discard.formula <- bf(respond_t ~ vot.1 + vot.2 + 
							(vot.1 + vot.2 | sFrame))
m.discard <- brm(discard.formula,
            data = subset(d.combined, unique.subjectID == subjs[1]),
            family = bernoulli(link="logit"),
            prior = priors.io_discard,
            chains = 4,
            thin = 1,
            warmup = 2000,
            iter = 6000, 
            control = list(adapt_delta = .9, max_treedepth = 10),
            refresh = 600,
            file = paste0("./models_individuals/catdiscard/catdiscard_model_all_ranefs_subject_", subjs[1], ".RDS")) 

contextonly.formula <- bf(respond_t ~ bias.numeric +
                            (bias.numeric | sFrame))
m.contextonly <- brm(contextonly.formula,
                     data = subset(d.combined, unique.subjectID == subjs[1]),
                     family = bernoulli(link="logit"),
                     prior = priors.io_discard,
                     chains = 4,
                     thin = 1,
                     warmup = 2000,
                     iter = 6000, 
                     control = list(adapt_delta = .9, max_treedepth = 10),
                     refresh = 100,
                     file = paste0("./models_individuals/contextonly/contextonly_model_all_ranefs_subject_", subjs[1], ".RDS"))
     
## now, loop through other subjects using update() to avoid recompiling
for (i in 2:length(subjs)) {
  update(m.discard,
  newdata = subset(d.combined, unique.subjectID == subjs[i]),
  recompile = FALSE,
  file = paste0("./models_individuals/catdiscard/catdiscard_model_all_ranefs_subject_", subjs[i], ".RDS"))
  
  update(m.io,
  newdata = subset(d.combined, unique.subjectID == subjs[i]),
  recompile = FALSE, control = list(adapt_delta = .9999),
  file = paste0("./models_individuals/ideal/idealobserver_model_all_ranefs_subject_", subjs[i], ".RDS"))

  update(m.ambiguity,
  newdata = subset(d.combined, unique.subjectID == subjs[i]),
  recompile = FALSE, control = list(adapt_delta = .9999),
  file = paste0("./models_individuals/ambiguity/ambiguity_model_all_ranefs_subject_", subjs[i], ".RDS"))

  update(m.contextonly,
  newdata = subset(d.combined, unique.subjectID == subjs[i]),
  recompile = FALSE, control = list(adapt_delta = .9999),
  file = paste0("./models_individuals/contextonly/contextonly_model_all_ranefs_subject_", subjs[i], ".RDS"))

  update(m.catswitch,
         newdata = subset(d.combined, unique.subjectID == subjs[i]),
         recompile = FALSE,
         file = paste0("./models_individuals/catswitch/catswitch_new_model_all_ranefs_subject_", subjs[i], ".RDS"))
}
