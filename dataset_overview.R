library(dplyr)
library(finalfit)
library(kableExtra)
library(rtf)

PATH_res <- "../best_models/"

load("../data_painful_painless.RData")

data_Painful_Painless$Outcome <- factor(data_Painful_Painless$Outcome, levels = c("Painful_neuropathy", "Painless_neuropathy"))

index <- which(data_Painful_Painless$Center!="Dundee")


explanatory = names(data_Painful_Painless)[!names(data_Painful_Painless) %in% c("Neuropathy", "MNSI_status", "DN4_status", "Depression_metric", "Anxiety_metric", "Set_index")]
dependent = "Set_index"
data_Painful_Painless %>%
summary_factorlist(dependent, explanatory, p = TRUE, na_include=TRUE, add_dependent_label=TRUE, column=TRUE, total_col=TRUE) -> t1

rtffile <- RTF(paste0(PATH_res,"train_test.table.doc"))
addTable(rtffile, t1)
done(rtffile)

knitr::kable(t1, row.names=FALSE, align=c("l", "l", "r", "r", "r")) %>% kable_styling("striped", "scale_down", full_width = F) %>% row_spec(which(t1$p == "<0.001"), bold = T, color = "white", background = "red") %>% save_kable(paste0(PATH_res,"data_overview_SetIndex_T1.html"))

knitr::kable(t1, row.names=FALSE, align=c("l", "l", "r", "r", "r")) %>% kable_styling("striped", "scale_down", full_width = F) %>% row_spec(which(t1$p == "<0.001"), bold = T, color = "white", background = "red") %>% save_kable(paste0(PATH_res,"data_overview_SetIndex_T1.pdf"))

#Remove CKD very low complexity factor too few  positives
#Remove outliers with too low HBA1c
data_Painful_Painless <- droplevels(data_Painful_Painless[!rownames(data_Painful_Painless) %in% rownames(data_Painful_Painless[data_Painful_Painless$HBA1C < 5 & !is.na(data_Painful_Painless$HBA1C),]), !names(data_Painful_Painless) %in% "CKD"])


explanatory = names(data_Painful_Painless[index,])[!names(data_Painful_Painless) %in% c("Outcome", "Neuropathy", "MNSI_status", "DN4_status", "Depression_metric", "Anxiety_metric", "Set_index")]
dependent = "Outcome"
data_Painful_Painless[index,] %>%
summary_factorlist(dependent, explanatory, p = TRUE, na_include=TRUE, add_dependent_label=TRUE, column=TRUE, total_col=TRUE) -> t2

rtffile <- RTF(paste0(PATH_res,"train.table.doc"))
addTable(rtffile, t2)
done(rtffile)

knitr::kable(t2, row.names=FALSE, align=c("l", "l", "r", "r", "r")) %>% kable_styling("striped", "scale_down", full_width = F) %>% row_spec(which(t2$p == "<0.001"), bold = T, color = "white", background = "red") %>% save_kable(paste0(PATH_res,"data_overview_Outcome_T2.html"))

knitr::kable(t2, row.names=FALSE, align=c("l", "l", "r", "r", "r")) %>% kable_styling("striped", "scale_down", full_width = F) %>% row_spec(which(t2$p == "<0.001"), bold = T, color = "white", background = "red") %>% save_kable(paste0(PATH_res,"data_overview_Outcome_T2.pdf"))

explanatory = names(data_Painful_Painless[index,])[!names(data_Painful_Painless) %in% c("Outcome", "Neuropathy", "MNSI_status", "DN4_status", "Depression_metric", "Anxiety_metric", "Set_index", "Center")]
dependent = "Outcome"
data_Painful_Painless[-index,] %>%
summary_factorlist(dependent, explanatory, p = TRUE, na_include=TRUE, add_dependent_label=TRUE, column=TRUE, total_col=TRUE) -> t3

rtffile <- RTF(paste0(PATH_res,"test.table.doc"))
addTable(rtffile, t3)
done(rtffile)

knitr::kable(t3, row.names=FALSE, align=c("l", "l", "r", "r", "r")) %>% kable_styling("striped", "scale_down", full_width = F) %>% row_spec(which(t3$p == "<0.001"), bold = T, color = "white", background = "red") %>% save_kable(paste0(PATH_res,"data_overview_Outcome_test_T3.html"))

knitr::kable(t3, row.names=FALSE, align=c("l", "l", "r", "r", "r")) %>% kable_styling("striped", "scale_down", full_width = F) %>% row_spec(which(t3$p == "<0.001"), bold = T, color = "white", background = "red") %>% save_kable(paste0(PATH_res,"data_overview_Outcome_test_T3.pdf"))



data_Painful_Painless[index,] %>%
  group_by(Gender, Age, HBA1C, Outcome) %>%
  summarise(n = n())

save(file = paste0(PATH_res, "data_Painful_Painless.RData"), data_Painful_Painless)  
  
  
