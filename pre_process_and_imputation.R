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

outcomeName <- "Outcome"
predictors <- names(data_Painful_Painless)[!(names(data_Painful_Painless) %in% c("Outcome", "Neuropathy", "MNSI_status", "MNSI_score", "DN4_status", "DN4_score", "Depression_metric", "Anxiety_metric", "Center", "Alcohol_consumption", "Alcohol_status", "Set_index", "PCS_score"))]

index <- which(data_Painful_Painless$Center!="Dundee")

missing_values <- data.frame(Variable = names(data_Painful_Painless), Missing_values = (sapply(data_Painful_Painless[index,], function(x) sum(is.na(x))))/nrow(data_Painful_Painless[index,]) )

missing_values_test <- data.frame(Variable = names(data_Painful_Painless), Missing_values = (sapply(data_Painful_Painless[-index,], function(x) sum(is.na(x))))/nrow(data_Painful_Painless[-index,]) )

intersect(missing_values[missing_values$Missing_values < 0.5,]$Variable, missing_values_test[missing_values_test$Missing_values < 0.5,]$Variable)

predictors <- predictors[predictors %in% intersect(missing_values[missing_values$Missing_values < 0.5,]$Variable, missing_values_test[missing_values_test$Missing_values < 0.5,]$Variable)]

trainSet <- data_Painful_Painless[index, c(predictors, outcomeName)]
testSet <- data_Painful_Painless[-index, c(predictors, outcomeName)]

names(trainSet)

var_check <- nearZeroVar(trainSet[,predictors], saveMetrics = TRUE)

numeric <- unlist(lapply(trainSet[,predictors], is.numeric))  

dataNumeric <- trainSet[,predictors[numeric]]
correlationMatrix <- cor(dataNumeric, use= "pairwise.complete.obs")
write.csv(file = paste0(PATH_ML, "/correlationMatrix.csv"), correlationMatrix)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.8, exact=TRUE, names=TRUE)
highlyCorrelated

predictors <- predictors[!predictors %in% highlyCorrelated & !predictors %in% rownames(var_check[var_check[,"zeroVar"] > 0, ])]
trainSet <- trainSet[, c(predictors, outcomeName)]
write.csv(file = paste0(PATH_ML, "/predictors_final_uncorellated.csv"), predictors)


library(Hmisc)

numeric <- unlist(lapply(trainSet[,predictors], is.numeric))  

pdf(paste0(PATH_ML, "/Predictors_histogram.pdf"))
hist.data.frame(trainSet[, predictors[numeric]])
dev.off()

data_dummies <- data_Painful_Painless[, c(outcomeName, predictors)]

dummies <- dummyVars(Outcome ~ ., data =  data_dummies, fullRank = TRUE)

data_Painful_Painless_dv <- predict(dummies, newdata = data_Painful_Painless)

data_Painful_Painless_dv <- as.data.frame(data_Painful_Painless_dv)

data_Painful_Painless_dv$Gender.Male <- as.factor(data_Painful_Painless_dv$Gender.Male)
data_Painful_Painless_dv$Trauma.TRUE <- as.factor(data_Painful_Painless_dv$Trauma.TRUE)
data_Painful_Painless_dv$Hospital_stay.TRUE <- as.factor(data_Painful_Painless_dv$Hospital_stay.TRUE)
data_Painful_Painless_dv$Ever_smoked_status.TRUE <- as.factor(data_Painful_Painless_dv$Ever_smoked_status.TRUE)

data_Painful_Painless_dv$Outcome <- data_Painful_Painless$Outcome

predictors_dv <- names(data_Painful_Painless_dv)[!names(data_Painful_Painless_dv) %in% "Outcome"]

trainSet_dv <- data_Painful_Painless_dv[index, c(predictors_dv, outcomeName)]
testSet_dv <- data_Painful_Painless_dv[-index, c(predictors_dv, outcomeName)]


trainSet <- data_Painful_Painless[index, c(predictors, outcomeName)]
testSet <- data_Painful_Painless[-index, c(predictors, outcomeName)]


trainSet_dv$Outcome <- factor(trainSet_dv$Outcome, levels = c("Painful_neuropathy", "Painless_neuropathy"))
testSet_dv$Outcome <- factor(testSet_dv$Outcome, levels = c("Painful_neuropathy", "Painless_neuropathy"))

trainSet$Outcome <- factor(trainSet$Outcome, levels = c("Painful_neuropathy", "Painless_neuropathy"))
testSet$Outcome <- factor(testSet$Outcome, levels = c("Painful_neuropathy", "Painless_neuropathy"))

print("Prepare imputed datasets")


imps <- 45

imp <- parlmice(trainSet,maxit=100,meth='pmm', cluster.seed=5623, n.core = 15, n.imp.core=3)
summary(imp)

pdf(paste0(PATH_ML, "/imputed_density_plots_painful_painless.pdf"), width=12, height=12)
densityplot(imp)
dev.off()

train_list <- foreach(i=1:imps) %dopar% {
require(mice)
complete(imp, i)
}

imp_dv <- parlmice(trainSet_dv,maxit=100,meth='pmm', cluster.seed=5623, n.core = 15, n.imp.core=3)

train_list_dv <- foreach(i=1:imps) %dopar% {
require(mice)
complete(imp_dv, i)
}

save(file = "/data/ndcn-diabetes_pn_ml/ndcn0569/diabetes_ML/Current_run/best_models/MI_dataset.RData", index, imp, imp_dv, train_list, train_list_dv, predictors, predictors_dv, outcomeName, trainSet, testSet, trainSet_dv, testSet_dv, imps)

###############################################
# Outcome agnostic imputation of test dataset #
###############################################
dummies <- dummyVars(Outcome ~ ., data =  data_Painful_Painless[, c(outcomeName, predictors)], fullRank = TRUE)

data_Painful_Painless_dv <- predict(dummies, newdata = data_Painful_Painless)

data_Painful_Painless_dv <- as.data.frame(data_Painful_Painless_dv)

data_Painful_Painless_dv$Gender.Male <- as.factor(data_Painful_Painless_dv$Gender.Male)
data_Painful_Painless_dv$Trauma.TRUE <- as.factor(data_Painful_Painless_dv$Trauma.TRUE)
data_Painful_Painless_dv$Hospital_stay.TRUE <- as.factor(data_Painful_Painless_dv$Hospital_stay.TRUE)
data_Painful_Painless_dv$Ever_smoked_status.TRUE <- as.factor(data_Painful_Painless_dv$Ever_smoked_status.TRUE)

test_dataset_dv <- data_Painful_Painless_dv[,predictors_dv]

test_dataset <- data_Painful_Painless[,predictors]
real_outcome <-factor(data_Painful_Painless[-index,outcomeName], levels = c("Painful_neuropathy", "Painless_neuropathy")) 

#Perform a separate multiple imputation on test dataset that is outcome agnostic. Notice that we also use the training dataset to help the multiple imputation but we have removed the outcome!
imp_test <- parlmice(test_dataset,maxit=100,meth='pmm',cluster.seed=5623, n.core = 15, n.imp.core=3)
summary(imp_test)

pdf(paste0(PATH_ML, "/TEST_TIME_imputed_density_plots_painful_painless.pdf"), width=12, height=12)
densityplot(imp_test)
dev.off()

#Create a list holding the completed, imputed datasets. Remove instance indicated by index, so we have a list of test instances only that have been imputed in an outcome agnostic way.
test_list <- foreach(i=1:imps) %dopar% {
require(mice)
data.frame(complete(imp_test, i)[-index,], Outcome=real_outcome)
}

imp_test_dv <- parlmice(test_dataset_dv,maxit=100,meth='pmm',cluster.seed=5623, n.core = 15, n.imp.core=3)

test_list_dv <- foreach(i=1:imps) %dopar% {
require(mice)
data.frame(complete(imp_test_dv, i)[-index,], Outcome=real_outcome)
}

save(file = "/data/ndcn-diabetes_pn_ml/ndcn0569/diabetes_ML/Current_run/best_models/MI_dataset_test_time.RData", imp_test, imp_test_dv, index, test_list, test_list_dv)

stopCluster(cl)
