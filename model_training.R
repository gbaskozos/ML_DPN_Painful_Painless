library(caret)
library(foreach)
library(doParallel)
library(parallel)
library(mice)
library(plyr)
library(dplyr)

set.seed(5623)

####################################################################
## Train ML algorithms for binary classification Painful Painless ##
####################################################################
PATH_ML <- "/data/ndcn-diabetes_pn_ml/ndcn0569/diabetes_ML/Current_run/best_models"
#PATH_ML <- "/home/george/Desktop/arc_ndcn_data_ML/diabetes_ML/Current_run/best_models"

n <- 16

print(paste0("Number of cores: ", n))

cl <- makePSOCKcluster(n)

registerDoParallel(cl)

load(paste0(PATH_ML,"/data_Painful_Painless.RData"))

load(paste0(PATH_ML, "/MI_dataset.RData"))

#### Model training, no predictor profiling
#caretFuncs$summary <- prSummary

predictors <- predictors[!predictors %in% c("MNSI_score", "DN4_score", "PCS_score")]
predictors_dv <- predictors_dv[!predictors_dv %in% c("MNSI_score", "DN4_score", "PCS_score")]

classWeights <- ifelse(data_Painful_Painless[index,]$Outcome == "Painful_neuropathy",
                        (1/table(data_Painful_Painless[index,]$Outcome)[1]) * 0.5,
                        (1/table(data_Painful_Painless[index,]$Outcome)[2]) * 0.5)

prSummary_mod <- function (data, lev = NULL, model = NULL) {
    require(mltools)
    require(caret)
    caret:::requireNamespaceQuietStop("MLmetrics")
    caret:::requireNamespaceQuietStop("pROC")
    if (length(levels(data$obs)) > 2) 
        stop(paste("Your outcome has", length(levels(data$obs)), 
            "levels. `prSummary`` function isn't appropriate.", 
            call. = FALSE))
    if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
        stop("Levels of observed and predicted data do not match.", 
            call. = FALSE)
    pr_auc <- try(MLmetrics::PRAUC(y_pred = data[, lev[1]], y_true = ifelse(data$obs == 
        lev[1], 1, 0)), silent = TRUE)
    if (inherits(pr_auc, "try-error")) 
        pr_auc <- NA
        sens <- sensitivity(data[, "pred"], data[, "obs"], lev[1])
        spec <- specificity(data[, "pred"], data[, "obs"], lev[2])
    c(AUC = pr_auc, Precision = caret:::precision.default(data = data$pred, 
        reference = data$obs, relevant = lev[1]), Recall = caret:::recall.default(data = data$pred, 
        reference = data$obs, relevant = lev[1]), F = caret:::F_meas.default(data = data$pred, 
        reference = data$obs, relevant = lev[1]), MCC = mcc(preds=data[, "pred"], actuals=data[, "obs"]), Sensitivity = sens, Specificity = spec, Balanced_Accuracy=((sens+spec)/2))
}



ctrl_train <- trainControl(method = "repeatedcv", classProbs = TRUE, number=10, repeats=5, allowParallel = TRUE, verboseIter=TRUE, savePredictions=TRUE, summaryFunction = prSummary_mod)

##############
# Accepts case weights #

print("Train RF model case weights")
model_rf_weighted <- foreach(i=1:imps) %do% {
require(caret)
train(as.formula(paste(outcomeName, "~", paste(predictors, collapse="+"))), data=train_list[[i]], method='ranger', trControl= ctrl_train, metric = "MCC", tuneLength = 60, weights = classWeights, importance = "impurity", preProc = c("center", "scale"))
}
save(file = "/data/ndcn-diabetes_pn_ml/ndcn0569/diabetes_ML/Current_run/best_models/RF_weighted.RData", model_rf_weighted)

print("Train Multivariate Adaptive Regression Splines model with case weights")
model_gcv_weighted <- foreach(i=1:imps) %do% {
require(caret)
train(as.formula(paste(outcomeName, "~", paste(predictors_dv, collapse="+"))), data=train_list_dv[[i]], method='gcvEarth', metric = "MCC", trControl= ctrl_train, tuneLength = 60, weights = classWeights, preProc = c("center", "scale"))
}
save(file = "/data/ndcn-diabetes_pn_ml/ndcn0569/diabetes_ML/Current_run/best_models/GCV_weighted.RData", model_gcv_weighted)

# Unweighted version

print("Train RF model ranger")
model_rf_ranger <- foreach(i=1:imps) %do% {
require(caret)
train(as.formula(paste(outcomeName, "~", paste(predictors, collapse="+"))), data=train_list[[i]], method='ranger', metric = "MCC", trControl= ctrl_train, tuneLength = 60, importance = "impurity", preProc = c("center", "scale"))
}
save(file = "/data/ndcn-diabetes_pn_ml/ndcn0569/diabetes_ML/Current_run/best_models/RF_ranger.RData", model_rf_ranger)


print("Train Multivariate Adaptive Regression Splines model")
model_gcv <- foreach(i=1:imps) %do% {
require(caret)
train(as.formula(paste(outcomeName, "~", paste(predictors_dv, collapse="+"))), data=train_list_dv[[i]], method='gcvEarth', metric = "MCC", trControl= ctrl_train, tuneLength = 60, preProc = c("center", "scale"))
}
save(file = "/data/ndcn-diabetes_pn_ml/ndcn0569/diabetes_ML/Current_run/best_models/GCV.RData", model_gcv)

##############################################################################################################################
#No case weights#

predictors_nb <- names(data_Painful_Painless)[names(data_Painful_Painless) %in% c("EQ5D_Index", "Depression_tscore", "Anxiety_tscore", "Alcohol_consumption_likert", "TIPIAgreeableness",     "TIPIEmotionalStability", "Trauma", "Age", "HBA1C", "BMI", "Gender")]

print("Train naive_bayes model")
model_naive_bayes <- foreach(i=1:imps) %do% {
require(caret)
train(as.formula(paste(outcomeName, "~", paste(predictors_nb, collapse="+"))), data=train_list[[i]], method='naive_bayes', metric = "MCC", trControl= ctrl_train, tuneLength = 60, preProc = c("center", "scale"))
}
save(file = "/data/ndcn-diabetes_pn_ml/ndcn0569/diabetes_ML/Current_run/best_models/naive_bayes.RData", model_naive_bayes)


stopCluster(cl)


