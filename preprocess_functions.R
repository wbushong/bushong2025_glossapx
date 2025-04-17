##################################
## Functions for loading & preprocessing data
##################################
library(tidyverse)

## ---- Exclude people with no VOT effect ----
initial_clean <- function(exp_file, expno) {
  d <- readRDS(exp_file)
  d$respond_t <- ifelse(d$respWord == "tent", 1, 0)
  d$VOT <- as.numeric(as.character(d$VOT))
  # update 2/6/25: for reasons unknown (I suspect an R update), different subjects
  # from Experiments 3 & 4 are removed during the below step than in initial analyses 
  # (from 2024). torn on what to do about this in the long-term, but for now am removing
  # subjects manually based on who was excluded in the initial 2024 analysis
  # and commenting out the old system
  # ultimately, in future work it may be more useful to not exclude any subjects
  # exclude.by.vot <- split(d, f = d$subject) %>%
  #   purrr::map(., ~ glm(respond_t ~ VOT, data = ., family = "binomial")) %>%
  #   purrr::map(summary) %>%
  #   purrr::map("coefficients") %>%
  #   purrr::map(., ~ as.data.frame(.)) %>%
  #   purrr::map(., function(.) .["VOT", "Pr(>|z|)"] < 0.05)
  # excludes <- which(exclude.by.vot == FALSE)
  # nexcl.1 <- length(excludes)
  m <- readRDS(paste0("models/ambiguity_model_all_ranefs_Experiment ", expno, ".RDS"))
  includes <- unique(m$data$subject)
  d <- filter(d, subject %in% as.character(includes))
  return(d)
}

## ---- Convert relevant variables to numeric & label experiments ----
assign_vars <- function(exp_file, audiosync_file) {
  expno <- unlist(strsplit(exp_file,"*"))
  expno <- expno[grep("[0-9]",expno)]
  
  d <- initial_clean(exp_file, expno)
  ## get numeric version of important vars
  d$vot.1.unstandardized <- poly(d$VOT, 2)[, 1]
  d$vot.2.unstandardized <- poly(d$VOT, 2)[, 2]
  # scale VOT by 2 SDs
  d$vot.1 <- with(d, vot.1.unstandardized / (2 * sd(vot.1.unstandardized)))
  d$vot.2 <- with(d, vot.2.unstandardized / (2 * sd(vot.2.unstandardized)))
  d$bias.numeric <- ifelse(d$context == "tent", 0.5, -0.5)
  d$distance.numeric <- ifelse(d$distance == "far", 0.5, -0.5)
  d$bias.fac <- factor(d$context, levels = c("tent", "dent"))
  contrasts(d$bias.fac) <- c(0.5, -0.5)
  d$distance.fac <- factor(d$distance, levels = c("short", "long"))
  contrasts(d$distance.fac) <- c(-0.5, 0.5)
  
  ## assign response type (hard coded -- exps 1/3 are forced response)
  d$response.type <- ifelse(length(grep("2|4", exp_file))!=0,
                            "free response",
                            "forced response")
  d$explabel <- paste("Experiment", expno)
  
  ## exp 4 didn't have annotated RTs yet
  if (expno == "4") {
    annotatedAudio <- read.csv("acousticPoints.csv")
    annotatedAudio$filename <- as.character(annotatedAudio$filename)
    d$filename <- with(d, paste0(paste(sFrame, distance, context, sep = "_"), ".wav")) # not true filename, an edited one
    d <- left_join(d, annotatedAudio, by = "filename")
    d$RT.wordOffset <- with(d, RT - wordOffset)
    d$RT.disambigOffset <- with(d, RT - semanticOffset)
  }

  ## for free-response exps, remove responses before 200ms after biasing context
  if ("RT.disambigOffset" %in% colnames(d)) {
    d <- subset(d, RT.disambigOffset > 200)
  }
  
  return(d)
}

## ---- Preprocess experiment files and combine into one file ----
preprocess_and_combine <- function(exp_files, audiosync_file) {
  d_combined <- c()
  for (i in exp_files) {
    d <- assign_vars(i, audiosync_file)
    rel_columns <- c("VOT", "vot.1", "vot.2", "context", "bias.fac", "bias.numeric", 
                     "distance", "distance.fac", "Trial", 
                     "explabel", "response.type", 
                     "subject", "sFrame", "respond_t")
    d_combined <- rbind(d_combined, d[, rel_columns])
  }
  contrasts(d_combined$bias.fac) <- c(0.5, -0.5)
  return(d_combined)
}
