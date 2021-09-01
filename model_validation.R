library(caret)
library(foreach)
library(doParallel)
library(parallel)
library(mice)
#Use this customised function
source("classMI.R")

set.seed(5623)

n <- 1

print(paste0("Number of cores: ", n))

cl <- makePSOCKcluster(n)

registerDoParallel(cl)

################

# Test time imputation
load("../best_models/data_Painful_Painless.RData")
load("../best_models/MI_dataset.RData")
load("../best_models/MI_dataset_test_time.RData")


##########################################
#Load the saved models and use the function classMI, save results in res_. Open the file "classMI.R" to see how it works.
load("../best_models/RF_weighted.RData")
res_rf_weighted <- classMI(model_rf_weighted, train_list, test_list, outcomeName, predictors, rfe=FALSE)

load("../best_models/GCV_weighted.RData")
res_gcv_weighted <- classMI(model_gcv_weighted, train_list_dv, test_list_dv, outcomeName, predictors_dv, rfe=FALSE)

######

load("../best_models/RF_ranger.RData")
res_rf_ranger <- classMI(model_rf_ranger, train_list, test_list, outcomeName, predictors, rfe=FALSE)

load("../best_models/GCV.RData")
res_gcv <- classMI(model_gcv, train_list_dv, test_list_dv, outcomeName, predictors_dv, rfe=FALSE)

#####
predictors_nb <- names(data_Painful_Painless)[names(data_Painful_Painless) %in% c("EQ5D_Index", "Depression_tscore", "Anxiety_tscore", "Alcohol_consumption_likert", "TIPIAgreeableness",     "TIPIEmotionalStability", "Trauma", "Age", "HBA1C", "BMI", "Gender")]

load("../best_models/naive_bayes.RData")
res_naive_bayes <- classMI(model_naive_bayes, train_list, test_list, outcomeName, predictors_nb, rfe=FALSE)


save(file = "../best_models/model_performance_metrics.RData", res_rf_weighted, res_gcv_weighted, res_rf_ranger, res_gcv, res_naive_bayes, predictors_dv, predictors, predictors_nb, train_list, test_list, train_list_dv, test_list_dv)

sink("../best_models/confusion_matrices.txt")
print("RF weighted\n")
res_rf_weighted$cMatrix
print("GCV weighted\n")
res_gcv_weighted$cMatrix
print("RF ranger\n")
res_rf_ranger$cMatrix
print("GCV\n")
res_gcv$cMatrix
print("NaiveBayes\n")
res_naive_bayes$cMatrix
sink()


write.csv(res_rf_weighted$pooledResampleProfile, file = "../best_models/res_rf_weighted_resamples.csv")
write.csv(res_gcv_weighted$pooledResampleProfile, file = "../best_models/res_gcv_weighted_resamples.csv")
write.csv(res_rf_ranger$pooledResampleProfile, file = "../best_models/res_rf_ranger_resamples.csv")
write.csv(res_gcv$pooledResampleProfile, file = "../best_models/res_gcv_resamples.csv")
write.csv(res_naive_bayes$pooledResampleProfile, file = "../best_models/res_naive_bayes_resamples.csv")


training_estimates <- rbind(res_rf_weighted$pooledResampleProfile[c("AUC", "Balanced_Accuracy", "MCC"),], res_rf_ranger$pooledResampleProfile[c("AUC", "Balanced_Accuracy", "MCC"),], res_gcv_weighted$pooledResampleProfile[c("AUC", "Balanced_Accuracy", "MCC"),], res_gcv$pooledResampleProfile[c("AUC", "Balanced_Accuracy", "MCC"),], res_naive_bayes$pooledResampleProfile[c("AUC", "Balanced_Accuracy", "MCC"),])
training_estimates$Model <- c(rep("Random Forest (weighted)",3), rep("Random Forest",3), rep("Adaptive Regression Splines (weighted)",3), rep("Adaptive Regression Splines",3), rep("NaiveBayes",3))

training_estimates <- droplevels(training_estimates)

levels(training_estimates$variable) <- c("AUPRC", "Balanced_Accuracy", "MCC")

ggplot(training_estimates, aes(x=variable, y=value, fill=Model)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value - sd, ymax=value + sd), width=.2,
                 position=position_dodge(.9), data = training_estimates, color='black', linetype = "dashed") + ggtitle("Performance estimates accross resamplings and imputations") + theme(axis.text.y = element_text(size= 10, face="bold"), axis.title.y = element_text(size=10, face="bold"), axis.title.x = element_blank(), axis.text.x = element_text(size= 12, face="bold"), legend.text=element_text(size=10, face="bold"), plot.title=element_text(size=14, face="bold", hjust = 0.5)) + scale_y_continuous(breaks=seq(0,1,0.05))

pdf(paste0(PATH_ML,"/model_training_estimates.pdf"), width = 12, height = 12)
print(ggplot(training_estimates, aes(x=variable, y=value, fill=Model)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value - sd, ymax=value + sd), width=.2,
                 position=position_dodge(.9), data = training_estimates, color='black', linetype = "solid") + ggtitle("Performance estimates accross resamplings and imputations") + theme(axis.text.y = element_text(size= 10, face="bold", angle = 45, vjust = 0.5, hjust=1), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(size= 10, face="bold", angle = 45, vjust = 0.8, hjust=1), legend.text=element_text(size=10, face="bold"), legend.position="top", plot.title=element_text(size=14, face="bold", hjust = 0.5)) + scale_y_continuous(breaks=seq(0,1,0.05)) + coord_flip()
)
dev.off()

#

#Create a dataframe with the probabilities of the out of fold predictions
classifier_probabilities <- data.frame(RF_weighted = res_rf_weighted$OOF_pred, RF_ranger = res_rf_ranger$OOF_pred, GCV_weighted = res_gcv_weighted$OOF_pred, GCV = res_gcv$OOF_pred, NaiveBayes = res_naive_bayes$OOF_pred)

write.csv(file = "../best_models/OOF_class_probs.csv", classifier_probabilities)


library("precrec")
scores_AUC <- join_scores(res_rf_weighted$testSet_pred$pred_probs, res_rf_ranger$testSet_pred$pred_probs, res_gcv_weighted$testSet_pred$pred_probs, res_gcv$testSet_pred$pred_probs, res_naive_bayes$testSet_pred$pred_probs)

msmdat <- mmdata(scores_AUC, test_list[[1]][,outcomeName], modnames = c("RF weighted", "RF (ranger)", "GCV weighted", "GCV", "NaiveBayes"), posclass="Painful_neuropathy")

mscurves <- evalmod(msmdat, posclass="Painful_neuropathy")

aucs <- auc(mscurves)
aucs_prc <- subset(aucs, curvetypes == "PRC")
aucs_prc$modnames <- factor(aucs_prc$modnames, levels = c("RF weighted", "RF (ranger)", "GCV weighted", "GCV", "NaiveBayes"))
names(aucs_prc)[[4]] <- "AUC"

#Compare performace estimated in training set and evaluated in test set
resample_test_performance <- data.frame(Model = rep(c("Random Forest (weighted)", "Random Forest", "Adaptive Regression Splines (weighted)", "Adaptive Regression Splines", "NaiveBayes"),2), Time = c(rep("Training",5), rep("Validation",5)), AUC=c(res_rf_weighted$pooledResampleProfile["AUC",]$value, res_rf_ranger$pooledResampleProfile["AUC",]$value, res_gcv_weighted$pooledResampleProfile["AUC",]$value, res_gcv$pooledResampleProfile["AUC",]$value, res_naive_bayes$pooledResampleProfile["AUC",]$value, aucs_prc$AUC))

write.csv(file = "../best_models/resample_test_performance.csv", resample_test_performance)

#Create a dataframe holding the achieved performance accross imputed test datasets
performance_table <- data.frame(Model = c("Random Forest (weighted)", "Random Forest", "Adaptive Regression Splines (weighted)", "Adaptive Regression Splines", "NaiveBayes"), Accuracy = c(res_rf_weighted$cMatrix$overall[[1]], res_rf_ranger$cMatrix$overall[[1]], res_gcv_weighted$cMatrix$overall[[1]], res_gcv$cMatrix$overall[[1]], res_naive_bayes$cMatrix$overall[[1]]), P_Value = c(res_rf_weighted$cMatrix$overall[[6]], res_rf_ranger$cMatrix$overall[[6]], res_gcv_weighted$cMatrix$overall[[6]], res_gcv$cMatrix$overall[[6]], res_naive_bayes$cMatrix$overall[[6]]), Balanced_Accuracy = c(res_rf_weighted$cMatrix$byClass[[11]], res_rf_ranger$cMatrix$byClass[[11]], res_gcv_weighted$cMatrix$byClass[[11]], res_gcv$cMatrix$byClass[[11]], res_naive_bayes$cMatrix$byClass[[11]]), MCC = c(res_rf_weighted$MCC, res_rf_ranger$MCC, res_gcv_weighted$MCC, res_gcv$MCC, res_naive_bayes$MCC), AUPRC = aucs_prc$AUC)

#Order the performance table
performance_table <- performance_table[order(-performance_table$MCC),]

write.csv(file = "../best_models/performance_table.csv", performance_table)

library(mltools)

resample_test_performance_MCC <- data.frame(Model = rep(c("Random Forest (weighted)", "Random Forest", "Adaptive Regression Splines (weighted)", "Adaptive Regression Splines", "NaiveBayes"),2), Time = c(rep("Training",5), rep("Validation",5)), MCC=c(res_rf_weighted$pooledResampleProfile["MCC",]$value, res_rf_ranger$pooledResampleProfile["MCC",]$value, res_gcv_weighted$pooledResampleProfile["MCC",]$value, res_gcv$pooledResampleProfile["MCC",]$value, res_naive_bayes$pooledResampleProfile["MCC",]$value, res_rf_weighted$MCC, res_rf_ranger$MCC, res_gcv_weighted$MCC, res_gcv$MCC, res_naive_bayes$MCC))

write.csv(file = "../best_models/resample_test_performance_MCC.csv", resample_test_performance_MCC)

resample_test_performance_Balanced_Accuracy <- data.frame(Model = rep(c("Random Forest (weighted)", "Random Forest", "Adaptive Regression Splines (weighted)", "Adaptive Regression Splines", "NaiveBayes"),2), Time = c(rep("Training",5), rep("Validation",5)), Balanced_Accuracy=c(res_rf_weighted$pooledResampleProfile["Balanced_Accuracy",]$value, res_rf_ranger$pooledResampleProfile["Balanced_Accuracy",]$value, res_gcv_weighted$pooledResampleProfile["Balanced_Accuracy",]$value, res_gcv$pooledResampleProfile["Balanced_Accuracy",]$value, res_naive_bayes$pooledResampleProfile["Balanced_Accuracy",]$value, res_rf_weighted$cMatrix$byClass[11], res_rf_ranger$cMatrix$byClass[11], res_gcv_weighted$cMatrix$byClass[11], res_gcv$cMatrix$byClass[11], res_naive_bayes$cMatrix$byClass[11]))

write.csv(file = "../best_models/resample_test_performance_Balanced_Accuracy.csv", resample_test_performance_Balanced_Accuracy)

MCC_plot <- ggplot(resample_test_performance_MCC) + geom_point(aes(x = MCC, y = Model, color = Time), size=2) + geom_vline(xintercept=0.25, color='black', linetype = "dashed") + ggtitle("Mathew's Correlation Coefficient") + theme(axis.text.y = element_text(size= 10, face="bold"), axis.title.y = element_text(size=10, face="bold"), axis.title.x = element_blank(), axis.text.x = element_text(size= 10, face="bold"), legend.text=element_text(size=10, face="bold"), plot.title=element_text(size=14, face="bold", hjust = 0.5))

b_accuracy_plot <- ggplot(resample_test_performance_Balanced_Accuracy) + geom_point(aes(x = Balanced_Accuracy, y = Model, color = Time), size=2) + geom_vline(xintercept=table(test_list[[1]]$Outcome)[[1]]/(table(test_list[[1]]$Outcome)[[1]]+ table(test_list[[1]]$Outcome)[[2]]), color='black', linetype = "dashed") + ggtitle("Balanced Accuracy") + theme(axis.text.y = element_text(size= 10, face="bold"), axis.title.y = element_text(size=10, face="bold"), axis.title.x = element_blank(), axis.text.x = element_text(size= 10, face="bold"), legend.text=element_text(size=10, face="bold"), plot.title=element_text(size=14, face="bold", hjust = 0.5))

AUPRC_plot <- ggplot(resample_test_performance) + geom_point(aes(x = AUC, y = Model, color = Time), size=2) + geom_vline(xintercept=0.7, color='black', linetype = "dashed") + ggtitle("Area Under the P/R curve") + theme(axis.text.y = element_text(size= 10, face="bold"), axis.title.y = element_text(size=10, face="bold"), axis.title.x = element_blank(), axis.text.x = element_text(size= 10, face="bold"), legend.text=element_text(size=10, face="bold"), plot.title=element_text(size=14, face="bold", hjust = 0.5))

perf_table_plot <- tableGrob(performance_table, rows=NULL, theme = ttheme_default(base_size = 9, base_colour = "black", base_family = "",
  parse = FALSE, padding = unit(c(4, 2), "mm")))

title <- ggdraw() + draw_label("Model validation",fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))

perf_table_grid <- plot_grid(perf_table_plot, labels=LETTERS[1])
model_validation_plot <- plot_grid(MCC_plot, b_accuracy_plot, AUPRC_plot, labels=LETTERS[2:4], ncol=1)

pdf(paste0(PATH_ML,"/model_validation.pdf"), width = 10, height = 10)
print(plot_grid(title,perf_table_grid,model_validation_plot, ncol = 1, rel_heights = c(0.1, 1,3)))
dev.off()


pdf("../best_models/train_test_MCC.pdf")
print(MCC_plot )
dev.off()

pdf("../best_models/train_test_balanced_accuracy.pdf")
print(b_accuracy_plot)
dev.off()

pdf("../best_models/train_test_AUC.pdf")
print(AUPRC_plot)
dev.off()


# Cluster clasifiers according to the probabilities they have assigned on out of fold samples
#Calculate 1- |cor|
dissimilarity <- 1 - abs(cor(as.matrix(classifier_probabilities), use="pairwise.complete.obs"))
#transform it to an euclidean distance matrix
distance <- as.dist(dissimilarity)
#perform hierarchical clustering and plot
pdf("../best_models/Correlation_probs.pdf")
plot(hclust(distance),  main="OOF similarity", xlab=NULL)
dev.off()




