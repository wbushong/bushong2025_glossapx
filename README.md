# Data, Analyses, & Supplementary Information for Bushong (2025) *Glossa Psycholinguistics*

This repository contains all data, analyses, and figures for Bushong (2025) Glossa Psycholinguistics. 

This repository includes the following directories and files:

## Supplementary Information

Written supplementary information to the paper. This file is entitled "SupplementaryInformation.pdf" in the main directory.

## Stimuli 
- acousticPoints.csv contains temporal information for each stimulus recording (temporal location of target word and subsequent context), which is used for excluding observations where participants responded before hearing biasing subsequent context
- stimuli.csv contains transcriptions of all stimulus sentences and indicates which items were used in which experiments

## Data 
- data/ contains the empirical data for Experiments 1-4 plus the norming study for the stimuli used in Experiments 3-4. Some preprocessing has already been done on these data (anonymizing participants, etc.)
- models/ contains the fitted models for each model & experiment, generated from fit_models.R
- models_individuals/ contains the fitted models for each individual participant, generated from fit_models_individuals.R
- model_comparisons/ contains the pairwise comparisons between each model for each experiment and each individual participant

## R Scripts 
- preprocess_functions.R contains various functions for preprocessing the data (removing participants w/ null VOT slopes, creating relevant standardized variables, etc.)
- norming_analysis.R contains the analysis for the norming study for the stimuli used in  Experiments 3-4
- fit_models.R defines and fits each of the five computational models to Experiments 1-4
- fit_models_individuals.R defines and fits each of the five models to each individual subject
- visualize_model_predictions.R creates the qualitative model predictions in the main text
- model_fit_figures.R generates the quantitative pairwise comparisons for each model, and makes the fit figures (Figure 6 in main text, Figures S3-6 in SI)
- model_fit_figures_individuals.R generates the quantitative pairwise comparisons for each individual subject, and makes corresponding fit figures (Figure 7 in main text)
- If you want to reproduce the analyses of this project, run fit_models.R followed by model_fit_figures.R

## Figures 
- figures/ contains all figures in the main text & SI
- model_fit_figures/ contains more detailed model fit figures for each model & experiment