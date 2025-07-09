library(readr)
library(tidyverse)
library(tmle)
data <- read_csv("C:/Users/tanzh/OneDrive/Desktop/Causal Inference/final project/student_dataset_canue_Version2_49900par_358var.csv")

# outcome: type 2 diabetes
data$type2_diab <- ifelse(data$DIS_DIAB_TYPE==2,1,0)

# exposure: smk_cig_ever (whether the participant has ever smoked at least 100 cigarettes during his lifetime)
data$SMK_CIG_EVER


data2 <- select(data, type2_diab, SMK_CIG_EVER, SDC_SEX, SDC_AGE_CALC, SDC_EDU_LEVEL, SDC_INCOME, ALC_EVER, DIS_DIAB_FAM_EVER)


# complete case for missing data
tmle_data <- data2 %>% na.omit()


# confounders: sex,age,edu level,income, ever consumed alcohol, family history of diabetes 
covariates <- tmle_data %>% select(SDC_SEX,SDC_AGE_CALC,SDC_EDU_LEVEL,SDC_INCOME,ALC_EVER, DIS_DIAB_FAM_EVER)

# TMLE
SL.library = c("SL.glm", 
               "SL.glmnet", 
               "SL.xgboost")

tmle.fit <- tmle(Y = tmle_data$type2_diab, 
                 A = tmle_data$SMK_CIG_EVER, 
                 W = covariates, 
                 family = "binomial",
                 V.Q = 3, #outcome model;
                 V.g = 3, #treatment model;
                 Q.SL.library = SL.library, 
                 g.SL.library = SL.library)


summary(tmle.fit)


# Propensity score weighting using generalized boosted modeling
library(WeightIt)

baselines <- colnames(covariates)
ps.formula <- as.formula(paste("SMK_CIG_EVER~", 
                               paste(baselines, collapse = "+")))

IPTW_gbm <- weightit(ps.formula,
                     data = tmle_data,
                     method = "gbm",
                     stabilize = TRUE)
# saving the model output as a R object to avoid rerunning the same model;
saveRDS(IPTW_gbm, file = "IPTW_gbm")

# reading saved model output;
require(sjPlot)
IPTW_gbm <- readRDS(file = "IPTW_gbm")
summary(IPTW_gbm)

fit2_gbm <- glm(type2_diab ~ SMK_CIG_EVER, 
                family = "binomial",
                weights = IPTW_gbm$weights,
                data = tmle_data)
tab_model(fit2_gbm)



# Propensity score weighting using SuperLearner

IPTW_SL <- weightit(ps.formula,
                    data = tmle_data,
                    method = "super",
                    SL.library=c("SL.randomForest", "SL.glmnet", "SL.nnet"), 
                    stabilize = TRUE)
# saving the model output as a R object to avoid rerunning the same model;
saveRDS(IPTW_SL, file = "IPTW_SL")

# reading saved model output;
IPTW_SL <- readRDS(file = "IPTW_SL")
summary(IPTW_SL)

fit2_SL <- glm(type2_diab ~ SMK_CIG_EVER, 
               family = "binomial",
               weights = IPTW_SL$weights,
               data = tmle_data)
tab_model(fit2_SL)

# PS model with Bayesian additive regression trees

IPTW_bart <- weightit(ps.formula,
                      data = tmle_data,
                      method = "bart",
                      stabilize = TRUE)
# saving the model output as a R object to avoid rerunning the same model;
saveRDS(IPTW_bart, file = "IPTW_bart")

# reading saved model output;
IPTW_bart <- readRDS(file = "IPTW_bart")
summary(IPTW_bart)

fit2_bart <- glm(type2_diab ~ SMK_CIG_EVER, 
                 family = "binomial",
                 weights = IPTW_bart$weights,
                 data = tmle_data)
tab_model(fit2_bart)

# weight diagnostic
library(cobalt)

bal.tab(IPTW_gbm) # For GBM
bal.tab(IPTW_SL) # For SuperLearner
bal.tab(IPTW_bart) # For BART


# descriptive table
library(tableone)

# Define variables to summarize
baseline_vars <- c("SDC_SEX", "SDC_AGE_CALC", "SDC_EDU_LEVEL", 
                   "SDC_INCOME", "ALC_EVER", "DIS_DIAB_FAM_EVER")

# Create factor versions if needed
tmle_data$SMK_CIG_EVER <- factor(tmle_data$SMK_CIG_EVER, labels = c("Never", "Ever"))
tmle_data$SDC_SEX <- factor(tmle_data$SDC_SEX, labels = c("Female", "Male"))
tmle_data$ALC_EVER <- factor(tmle_data$ALC_EVER, labels = c("Never", "Ever"))
tmle_data$DIS_DIAB_FAM_EVER <- factor(tmle_data$DIS_DIAB_FAM_EVER, labels = c("No", "Yes", "Presume No"))

# Create Table 1 stratified by smoking
table1 <- CreateTableOne(vars = baseline_vars, strata = "SMK_CIG_EVER", data = tmle_data)
print(table1, showAllLevels = TRUE)
