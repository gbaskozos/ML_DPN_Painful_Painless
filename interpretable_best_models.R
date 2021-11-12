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
library(cowplot)

aggregate_pdp <- function(model, feature, train_data, y="Outcome", class, method="pdp+ice") {
model_list <- mapply(Predictor$new, model, train_data, y=y, type="prob", class=class)

model_res_1 <- lapply(model_list, FeatureEffect$new, method = method, feature = feature, grid.size=100)
model_res_2 <- lapply(model_res_1, getElement, name="results")
template <- model_res_1[[1]]
rm(model_res_1)

composite_df <- do.call(rbind, model_res_2)
agg_pdp <- aggregate(as.formula(paste(".value ~",eval(feature), "+ .type")), composite_df[composite_df$".type" == "pdp",], mean)
if (method=="pdp+ice") {
agg_ice <- aggregate(as.formula(paste(".value ~",eval(feature), "+ .type + .id")), composite_df[composite_df$".type" == "ice",], mean)
rm(composite_df)

agg_pdp$".id"=NA

template$results <- rbind(agg_pdp, agg_ice)
}
else if (method=="pdp") {
template$results <- agg_pdp
}

p <- plot(template)
return(p)

}




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

top_predictors_gcv <- res_gcv$varImpROC[res_gcv$varImpROC$Importance > 10,]
top_predictors_rf <- res_rf_ranger$varImpROC[res_rf_ranger$varImpROC$Importance > 10,]
top_predictors_naive_bayes <- res_naive_bayes$varImpROC[res_naive_bayes$varImpROC$Importance > 10,]

save(file = paste0(PATH_ML, "/top_predictors.RData"), top_predictors_gcv, top_predictors_rf, top_predictors_naive_bayes)

load(paste0(PATH_ML, "/top_predictors.RData"))

#GCV

load(paste0(PATH_ML, "/GCV.RData"))

effs_gcv_plot <- vector(mode = "list", length = length(top_predictors_gcv$Variable))

for(i in 1:length(top_predictors_gcv$Variable)) {
effs_gcv_plot[[i]] <- aggregate_pdp(model= model_gcv, feature = top_predictors_gcv$Variable[i], train_data = train_list_dv, y="Outcome", class="Painful_neuropathy", method="pdp") + geom_line(size=2, colour = 'orange', alpha=0.5) + geom_smooth(span=1)
}

#effs_gcv_plot[[length(effs_gcv_plot)+1]] <- tableGrob(top_predictors_gcv[,c(1,2)])

gcv_plots <- plot_grid(plotlist=effs_gcv_plot, labels=LETTERS[1:length(effs_gcv_plot)])

title <- ggdraw() + draw_label("Partial dependence plots Painful Neuropathy, Adaptive Regression Splines",fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))

pdf(paste0(PATH_ML,"/predictor_effects_gcv.pdf"), width = 11, height = 8)
print(plot_grid(title, gcv_plots, ncol = 1, rel_heights = c(0.1, 1)))
dev.off()

load(paste0(PATH_ML, "/RF_ranger.RData"))

effs_rf_plot <- vector(mode = "list", length = length(top_predictors_rf$Variable))

for(i in 1:length(top_predictors_rf$Variable)) {
effs_rf_plot[[i]] <- aggregate_pdp(model= model_rf_ranger, feature = top_predictors_rf$Variable[i], train_data = train_list, y="Outcome", class="Painful_neuropathy", method="pdp") + geom_line(size=2, colour = 'orange', alpha=0.5) + geom_smooth(span=1)
}

#effs_rf_plot[[length(effs_rf_plot)+1]] <- tableGrob(top_predictors_rf[,c(1,2)])

rf_plots <- plot_grid(plotlist=effs_rf_plot, labels=LETTERS[1:length(effs_rf_plot)], ncol=3)

title <- ggdraw() + draw_label("Partial dependence plots Painful Neuropathy, Random Forest",fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))

pdf(paste0(PATH_ML,"/predictor_effects_rf.pdf"), width = 11, height = 13)
print(plot_grid(title, rf_plots, ncol = 1, rel_heights = c(0.1, 1)))
dev.off()


load(paste0(PATH_ML, "/naive_bayes.RData"))

effs_naive_bayes_plot <- vector(mode = "list", length = length(top_predictors_naive_bayes$Variable))

for(i in 1:length(top_predictors_naive_bayes$Variable)) {
effs_naive_bayes_plot[[i]] <- aggregate_pdp(model= model_naive_bayes, feature = top_predictors_naive_bayes$Variable[i], train_data = train_list, y="Outcome", class="Painful_neuropathy", method="pdp") + geom_line(size=2, colour = 'orange', alpha=0.5) + geom_smooth(span=1)
}

naive_bayes_plots <- plot_grid(plotlist=effs_naive_bayes_plot, labels=LETTERS[1:length(effs_naive_bayes_plot)], ncol=3)

title <- ggdraw() + draw_label("Partial dependence plots Painful Neuropathy, Naive Bayes",fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))

pdf(paste0(PATH_ML,"/predictor_effects_naive_bayes.pdf"), width = 11, height = 11)
print(plot_grid(title, naive_bayes_plots, ncol = 1, rel_heights = c(0.1, 1)))
dev.off()

save(file = paste0(PATH_ML, "/plots.RData"), naive_bayes_plots, rf_plots, gcv_plots)


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


