library(caret)
library(foreach)
library(doParallel)
library(parallel)
library(mice)
library(iml)
library(forcats)
#Use this customised function
source("classMI.R")
library(gridExtra)



PATH_ML <- "../best_models"


n <- 3

print(paste0("Number of cores: ", n))

cl <- makePSOCKcluster(n)

registerDoParallel(cl)

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


set.seed(5623)

#load(paste0(PATH_ML,"/data_Painful_Painless.RData"))

load(paste0(PATH_ML,"/MI_dataset.RData"))

load(paste0(PATH_ML,"/MI_dataset_test_time.RData"))

load(paste0(PATH_ML,"/model_performance_metrics.RData"))

# keep GCV, RF, stepped_LDA, naive+_bayes 

p1 <- res_gcv$varImp_plot + ggtitle("Adaptive Regression Splines")
p2 <- res_rf_ranger$varImp_plot + ggtitle("Random Forest")
p3 <- res_naive_bayes$varImp_plot + ggtitle("Naive Bayes")

p4 <- res_gcv_weighted$varImp_plot + ggtitle("Adaptive Regression Splines (weighted)")
p5 <- res_rf_weighted$varImp_plot + ggtitle("Random Forest (weighted)")

title <- ggdraw() + draw_label("Pooled Variable Importance",fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))

p_importance <- plot_grid(p1, p2, p3, p4, p5, labels=LETTERS[1:5])

pdf(paste0(PATH_ML,"/variable_importance.pdf"), width = 13, height = 13)
print(plot_grid(title, p_importance, ncol = 1, rel_heights = c(0.1, 1)))
dev.off()

best_predictors <- unique(c(rownames(res_gcv$varImpROC[c(1:5),]), rownames(res_rf_ranger$varImpROC[c(1:5),]), rownames(res_naive_bayes$varImpROC[c(1:5),])))


#GCV
best_gcv_train = Predictor$new(res_gcv$bestModel, data = train_list_dv[[res_gcv$bestModel_MI_index]], y = outcomeName, type = "prob")

effs_gcv_train <- lapply(gsub("[^a-zA-Z.5_]", "", rownames(varImp(res_gcv$bestModel)$importance)[varImp(res_gcv$bestModel)$importance > 10]), FeatureEffect$new, predictor=best_gcv_train, method = "pdp+ice", grid.size = 100)

effs_gcv_plot <- lapply(effs_gcv_train, plot)

effs_gcv_plot[[length(effs_gcv_plot)+1]] <- tableGrob(data.frame(Variable = rownames(varImp(res_gcv$bestModel)$importance)[varImp(res_gcv$bestModel)$importance > 10], Importance = varImp(res_gcv$bestModel)$importance[rownames(varImp(res_gcv$bestModel)$importance)[varImp(res_gcv$bestModel)$importance > 10],]))

gcv_plots <- plot_grid(plotlist=effs_gcv_plot, labels=LETTERS[1:length(effs_gcv_plot)])

title <- ggdraw() + draw_label("Partial dependence plots, Adaptive Regression Splines (best model)",fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))

pdf(paste0(PATH_ML,"/predictor_effects_gcv.pdf"), width = 12, height = 8)
print(plot_grid(title, gcv_plots, ncol = 1, rel_heights = c(0.1, 1)))
dev.off()

#GCV_weighted
best_gcv_weighted_train = Predictor$new(res_gcv_weighted$bestModel, data = train_list_dv[[res_gcv_weighted$bestModel_MI_index]], y = outcomeName, type = "prob")
effs_gcv_weighted_train <- FeatureEffects$new(best_gcv_weighted_train, feature = gsub("[^a-zA-Z.5_]", "", rownames(varImp(res_gcv_weighted$bestModel)$importance)[varImp(res_gcv_weighted$bestModel)$importance > 5]), method = "pdp+ice", grid.size = 100)

pdf(paste0(PATH_ML,"/predictor_effects_gcv_weighted.pdf"), width = 12, height = 10)
plot(effs_gcv_weighted_train) + ggtitle("Partial dependence plots, Adaptive regression splines (weighted)")
dev.off()

#RF
best_rf_train = Predictor$new(res_rf_ranger$bestModel, data = train_list[[res_rf_ranger$bestModel_MI_index]], y = outcomeName, type = "prob")
effs_rf_train <- lapply(rownames(varImp(res_rf_ranger$bestModel)$importance)[varImp(res_rf_ranger$bestModel)$importance > 10], FeatureEffect$new, predictor=best_rf_train, method = "pdp+ice", grid.size = 100)

gender_BMI <- FeatureEffect$new(best_rf_train, feature = c("BMI", "Gender"), method = "pdp", grid.size = 100) 
gender_depression <- FeatureEffect$new(best_rf_train, feature = c("Depression_tscore", "Gender"), method = "pdp", grid.size = 100) 
gender_anxiety <- FeatureEffect$new(best_rf_train, feature = c("Anxiety_tscore", "Gender"), method = "pdp", grid.size = 100) 
gender_HBA1C <- FeatureEffect$new(best_rf_train, feature = c("HBA1C", "Gender"), method = "pdp", grid.size = 100) 
gender_alcohol <- FeatureEffect$new(best_rf_train, feature = c("Alcohol_consumption_likert", "Gender"), method = "pdp", grid.size = 100) 
gender_smoke <- FeatureEffect$new(best_rf_train, feature = c("Ever_smoked_status", "Gender"), method = "pdp", grid.size = 100)
gender_age <- FeatureEffect$new(best_rf_train, feature = c("Age", "Gender"), method = "pdp", grid.size = 100)

rf_interactions <- plot_grid(plot(gender_BMI), plot(gender_depression), plot(gender_anxiety), plot(gender_HBA1C), plot(gender_alcohol),  labels=LETTERS[1:5])

title <- ggdraw() + draw_label("Partial dependence plots, Interactions with Gender",fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))

pdf(paste0(PATH_ML,"/predictor_interactions_rf.pdf"), width = 14, height = 8)
print(plot_grid(title, rf_interactions, ncol = 1, rel_heights = c(0.1, 2.5)))
dev.off()


effs_rf_plot <- lapply(effs_rf_train, plot)

effs_rf_plot[[length(effs_rf_plot)+1]] <- tableGrob(data.frame(Variable = rownames(varImp(res_rf_ranger$bestModel)$importance)[varImp(res_rf_ranger$bestModel)$importance > 10], Importance = varImp(res_rf_ranger$bestModel)$importance[rownames(varImp(res_rf_ranger$bestModel)$importance)[varImp(res_rf_ranger$bestModel)$importance > 10],]), theme = ttheme_default(base_size = 9, base_colour = "black", base_family = "",
  parse = FALSE, padding = unit(c(4, 2), "mm")))

rf_plots <- plot_grid(plotlist=effs_rf_plot, labels=LETTERS[1:length(effs_rf_plot)])

title <- ggdraw() + draw_label("Partial dependence plots, Radom Forest (best model)",fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))

pdf(paste0(PATH_ML,"/predictor_effects_rf.pdf"), width = 12, height = 12)
print(plot_grid(title, rf_plots, ncol = 1, rel_heights = c(0.1, 2.5)))
dev.off()


#RF_weighted
best_rf_weighted_train = Predictor$new(res_rf_weighted$bestModel, data = train_list[[res_rf_weighted$bestModel_MI_index]], y = outcomeName, type = "prob")
effs_rf_weighted_train <- FeatureEffects$new(best_rf_weighted_train, feature = gsub("TRUE", "", rownames(varImp(res_rf_weighted$bestModel)$importance)[varImp(res_rf_weighted$bestModel)$importance > 5]), method = "pdp+ice", grid.size = 100)

pdf(paste0(PATH_ML,"/predictor_effects_rf_weighted.pdf"), width = 12, height = 10)
plot(effs_rf_weighted_train) + ggtitle("Partial dependence plots, Random Forest (weighted)")
dev.off()

#Naive Bayes
best_naive_bayes_train = Predictor$new(res_naive_bayes$bestModel, data = train_list[[res_naive_bayes$bestModel_MI_index]], y = outcomeName, type = "prob")
effs_naive_bayes_train <- lapply(rownames(varImp(res_naive_bayes$bestModel)$importance)[varImp(res_naive_bayes$bestModel)$importance$Painful_neuropathy > 10], FeatureEffect$new, predictor=best_naive_bayes_train, method = "pdp+ice", grid.size = 100)

effs_naive_bayes_plot <- lapply(effs_naive_bayes_train, plot)

effs_naive_bayes_plot[[length(effs_naive_bayes_plot)]] <- tableGrob(data.frame(Variable = rownames(varImp(res_naive_bayes$bestModel)$importance)[varImp(res_naive_bayes$bestModel)$importance$Painful_neuropathy > 10], Importance = varImp(res_naive_bayes$bestModel)$importance[rownames(varImp(res_naive_bayes$bestModel)$importance)[varImp(res_naive_bayes$bestModel)$importance$Painful_neuropathy > 10],]$Painful_neuropathy), theme = ttheme_default(base_size = 10, base_colour = "black", base_family = "", parse = FALSE, padding = unit(c(4, 2), "mm")))

naive_bayes_plots <- plot_grid(plotlist=effs_naive_bayes_plot, labels=LETTERS[1:length(effs_naive_bayes_plot)])

title <- ggdraw() + draw_label("Partial dependence plots, Naive Bayes (best model)",fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))

pdf(paste0(PATH_ML,"/predictor_effects_naive_bayes.pdf"), width = 13, height = 9)
print(plot_grid(title, naive_bayes_plots, ncol = 1, rel_heights = c(0.1, 1)))
dev.off()

#Lift curves and calibration plots
lift_results <- data.frame(Class = testSet$Outcome)

lift_results$RF <- res_rf_ranger$testSet_pred$pred_probs
lift_results$RF_weighted <- res_rf_weighted$testSet_pred$pred_probs
lift_results$GCV <- res_gcv$testSet_pred$pred_probs
lift_results$GCV_weighted <- res_gcv_weighted$testSet_pred$pred_probs
lift_results$Naive_Bayes <- res_naive_bayes$testSet_pred$pred_probs

trellis.par.set(caretTheme())
lift_obj <- caret::lift(Class ~ RF + RF_weighted + GCV + GCV_weighted + Naive_Bayes, data = lift_results)

pdf(paste0(PATH_ML,"/lift_curves.pdf"))
plot(lift_obj, values = 60, auto.key = list(columns = 5, lines = TRUE, points = FALSE))
dev.off()

cal_obj <- caret::calibration(Class ~ RF + RF_weighted + GCV + GCV_weighted + Naive_Bayes,
                       data = lift_results,
                       cuts = 5)

calib_data <- cal_obj$data

calibration <- data.frame(sapply(by(calib_data, calib_data$calibModelVar, function(z) coef(lm(Percent ~ midpoint, data = z))), cbind), row.names=c("Intercept", "Slope"))

names(calibration) <- c("RF", "RF_weighted", "Adaptive_Regression_Splines", "Adaptive_Regression_Splines_weighted", "Naive_Bayes")

write.csv(file = paste0(PATH_ML, "/calibration_slope_intercept.csv"), calibration)

tableGrob(calibration, theme = ttheme_default(base_size = 10, base_colour = "black", base_family = "", parse = FALSE, padding = unit(c(4, 4), "mm")))

cal_plots <- plot_grid(plot(cal_obj), tableGrob(t(calibration), theme = ttheme_default(base_size = 10, base_colour = "black", base_family = "", parse = FALSE, padding = unit(c(4, 4), "mm"))), labels=LETTERS[1:2])

title <- ggdraw() + draw_label("Model calibration on validation dataset",fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))

pdf(paste0(PATH_ML,"/calibration.pdf"), width = 12, height = 6)
print(plot_grid(title, cal_plots, ncol = 1, rel_heights = c(0.1, 1)))
dev.off()
