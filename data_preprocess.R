library(FactoMineR)
library(missMDA)
library("factoextra")
library(VIM)
library(plotrix)
library(finalfit)
library(MissMech)
library(naniar)
library(eq5d)
library(ggplot2)
library(finalfit)
library(plyr)
library(dplyr)
#library(MissMech)
library(naniar)
library(kableExtra)


PATH="/home/george/Desktop/arc_ndcn_data_ML/diabetes_ML/Current_run/"

options(stringsAsFactors=FALSE)

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

Q_format <- function(x) {

x[is.na(x)] <- 0

x[x=="N/A"] <- 0

x[x==""] <- 0

x[x!=0] <- 1

return(as.numeric(x))

}

pMiss <- function(x){sum(is.na(x))/length(x)*100}

#dundee <- read.table(paste0(PATH,"raw_data/DOLORisk_GoDARTS.csv"), sep = "|", na.strings=c("","na", "NA"), stringsAsFactors=FALSE, header=TRUE, strip.white=TRUE, quote = "", fill = TRUE)

DOLORISK <- read.csv(paste0(PATH,"raw_data/DataExport19_02_2020.csv"), na.strings=c("","na", "NA", "Not calculated – Missing data", "ND", " Not calculated – Missing data", "Not calculated – Not enough data"), stringsAsFactors=FALSE) 
oxford_imperial <- DOLORISK[(DOLORISK$Centre == 1 | DOLORISK$Centre == 2 | DOLORISK$Centre == 8) & DOLORISK$Q75a == "Diabetic polyneuropathy" & !is.na(DOLORISK$Q75a), ]
oxford_imperial <- oxford_imperial[-which(duplicated(oxford_imperial$DOLORiskNum)), ]
oxford_imperial_data <- data.frame(matrix(ncol=0,nrow=nrow(oxford_imperial)))
oxford_imperial_data$DOLORiskNum <- oxford_imperial$DOLORiskNum

rownames(oxford_imperial_data) <- oxford_imperial_data$DOLORiskNum

oxford_imperial_data$Center <- oxford_imperial$Centre

#oxford_imperial feature extraction
#EQ5D
eq5d.scores.df <- data.frame(MO=oxford_imperial$Q49a, SC=oxford_imperial$Q49b, UA=oxford_imperial$Q49c, PD=oxford_imperial$Q49d, AD=oxford_imperial$Q49e)

eq5d.scores.df$MO[grep("I have no", eq5d.scores.df$MO)] <- 1 
eq5d.scores.df$MO[grep("slight", eq5d.scores.df$MO)] <- 2 
eq5d.scores.df$MO[grep("moderate", eq5d.scores.df$MO)] <- 3 
eq5d.scores.df$MO[grep("severe", eq5d.scores.df$MO)] <- 4 
eq5d.scores.df$MO[grep("unable", eq5d.scores.df$MO)] <- 5 
eq5d.scores.df$MO <- as.numeric(eq5d.scores.df$MO)

eq5d.scores.df$SC[grep("I have no", eq5d.scores.df$SC)] <- 1 
eq5d.scores.df$SC[grep("slight", eq5d.scores.df$SC)] <- 2 
eq5d.scores.df$SC[grep("moderate", eq5d.scores.df$SC)] <- 3 
eq5d.scores.df$SC[grep("severe", eq5d.scores.df$SC)] <- 4 
eq5d.scores.df$SC[grep("unable", eq5d.scores.df$SC)] <- 5 
eq5d.scores.df$SC <- as.numeric(eq5d.scores.df$SC)

eq5d.scores.df$UA[grep("I have no", eq5d.scores.df$UA)] <- 1 
eq5d.scores.df$UA[grep("slight", eq5d.scores.df$UA)] <- 2 
eq5d.scores.df$UA[grep("moderate", eq5d.scores.df$UA)] <- 3 
eq5d.scores.df$UA[grep("severe", eq5d.scores.df$UA)] <- 4 
eq5d.scores.df$UA[grep("unable", eq5d.scores.df$UA)] <- 5 
eq5d.scores.df$UA <- as.numeric(eq5d.scores.df$UA)

eq5d.scores.df$PD[grep("I have no", eq5d.scores.df$PD)] <- 1 
eq5d.scores.df$PD[grep("slight", eq5d.scores.df$PD)] <- 2 
eq5d.scores.df$PD[grep("moderate", eq5d.scores.df$PD)] <- 3 
eq5d.scores.df$PD[grep("severe", eq5d.scores.df$PD)] <- 4 
eq5d.scores.df$PD[grep("extreme", eq5d.scores.df$PD)] <- 5 
eq5d.scores.df$PD <- as.numeric(eq5d.scores.df$PD)

eq5d.scores.df$AD[grep("I am not", eq5d.scores.df$AD)] <- 1 
eq5d.scores.df$AD[grep("slightly", eq5d.scores.df$AD)] <- 2 
eq5d.scores.df$AD[grep("moderately", eq5d.scores.df$AD)] <- 3 
eq5d.scores.df$AD[grep("severely", eq5d.scores.df$AD)] <- 4 
eq5d.scores.df$AD[grep("extremely", eq5d.scores.df$AD)] <- 5 
eq5d.scores.df$AD <- as.numeric(eq5d.scores.df$AD)


oxford_imperial$EQ5D_Index <- NA 
oxford_imperial$EQ5D_Index[complete.cases(eq5d.scores.df)] <- eq5d(na.omit(eq5d.scores.df), country = "UK", version = "5L", type = "CW")

oxford_imperial_data$EQ5D_Index <- oxford_imperial$EQ5D_Index

#PROMIS
oxford_imperial_data$Depression_tscore <- as.numeric(oxford_imperial$PROMISDepressionTScore)
oxford_imperial_data$Anxiety_tscore <- as.numeric(oxford_imperial$PROMISAnxietyTScore)
oxford_imperial_data$Sleep_Disturbance_tscore <- as.numeric(oxford_imperial$PROMISSleepTScore)

#Trauma
oxford_imperial_data$Trauma <- NA
oxford_imperial_data$Trauma[grep("Yes:", oxford_imperial$Q54)] <- 1 
oxford_imperial_data$Trauma[grep("None", oxford_imperial$Q54)] <- 0 
oxford_imperial_data$Trauma <- as.factor(oxford_imperial_data$Trauma)
levels(oxford_imperial_data$Trauma)[levels(oxford_imperial_data$Trauma)=="1"] <- TRUE 
levels(oxford_imperial_data$Trauma)[levels(oxford_imperial_data$Trauma)=="0"] <- FALSE 

#Hospital stay
oxford_imperial_data$Hospital_stay <- NA
oxford_imperial_data$Hospital_stay[grep("Yes",oxford_imperial$Q55)] <- 1 
oxford_imperial_data$Hospital_stay[grep("No",oxford_imperial$Q55)] <- 0 
oxford_imperial_data$Hospital_stay <- as.factor(oxford_imperial_data$Hospital_stay)
levels(oxford_imperial_data$Hospital_stay)[levels(oxford_imperial_data$Hospital_stay)=="1"] <- TRUE 
levels(oxford_imperial_data$Hospital_stay)[levels(oxford_imperial_data$Hospital_stay)=="0"] <- FALSE 

#TIPI
oxford_imperial_data$TIPIExtraversion <- oxford_imperial$TIPIExtraversion
oxford_imperial_data$TIPIAgreeableness <- oxford_imperial$TIPIAgreeableness
oxford_imperial_data$TIPIConscientiousness <- oxford_imperial$TIPIConscientiousness
oxford_imperial_data$TIPIEmotionalStability <- oxford_imperial$TIPIEmotionalStability
oxford_imperial_data$TIPIOpenness <- oxford_imperial$TIPIOpenness

#Smoking
oxford_imperial_data$Current_smoker <- oxford_imperial$Q60Smoker
oxford_imperial_data$Ever_smoked_status <- NA
oxford_imperial_data$Ever_smoked_status[oxford_imperial_data$Current_smoker == 1 | grepl("Yes", oxford_imperial$Q58)] <- 1
oxford_imperial_data$Ever_smoked_status[oxford_imperial_data$Current_smoker == 0 & grepl("No", oxford_imperial$Q58)] <- 0 
oxford_imperial_data$Ever_smoked_status <- as.factor(oxford_imperial_data$Ever_smoked_status)
oxford_imperial_data$Current_smoker <- as.factor(oxford_imperial_data$Current_smoker)
levels(oxford_imperial_data$Ever_smoked_status)[levels(oxford_imperial_data$Ever_smoked_status)=="1"] <- TRUE 
levels(oxford_imperial_data$Ever_smoked_status)[levels(oxford_imperial_data$Ever_smoked_status)=="0"] <- FALSE 

#Alcohol
#Keep full alcohol consumption information
oxford_imperial_data$Alcohol_consumption <- ordered(oxford_imperial$Q62)
oxford_imperial_data$Alcohol_consumption <- ordered(oxford_imperial_data$Alcohol_consumption, levels = c("Never", "Less than 1 day per month", "1 to 3 days per month", "1 or 2 days per week", "3 or 4 days per week", "Daily or almost daily"))
#Alcohol_status[Alcohol_status=="Never" & !is.na(Alcohol_status)] <- 0
#Alcohol_status[Alcohol_status!=0 & !is.na(Alcohol_status)] <- 1
oxford_imperial_data$Alcohol_consumption_likert <- oxford_imperial_data$Alcohol_consumption
levels(oxford_imperial_data$Alcohol_consumption_likert)[levels(oxford_imperial_data$Alcohol_consumption_likert)=="Never"] <- 0
levels(oxford_imperial_data$Alcohol_consumption_likert)[levels(oxford_imperial_data$Alcohol_consumption_likert)=="Less than 1 day per month"] <- 1
levels(oxford_imperial_data$Alcohol_consumption_likert)[levels(oxford_imperial_data$Alcohol_consumption_likert)=="1 to 3 days per month"] <- 2
levels(oxford_imperial_data$Alcohol_consumption_likert)[levels(oxford_imperial_data$Alcohol_consumption_likert)=="1 or 2 days per week"] <- 3
levels(oxford_imperial_data$Alcohol_consumption_likert)[levels(oxford_imperial_data$Alcohol_consumption_likert)=="3 or 4 days per week"] <- 4
levels(oxford_imperial_data$Alcohol_consumption_likert)[levels(oxford_imperial_data$Alcohol_consumption_likert)=="Daily or almost daily"] <- 5
oxford_imperial_data$Alcohol_consumption_likert <- as.numeric(as.character(oxford_imperial_data$Alcohol_consumption_likert))

oxford_imperial_data$Alcohol_status <- NA
oxford_imperial_data$Alcohol_status[oxford_imperial_data$Alcohol_consumption_likert == 0] <- FALSE
oxford_imperial_data$Alcohol_status[oxford_imperial_data$Alcohol_consumption_likert > 0] <- TRUE
oxford_imperial_data$Alcohol_status <- as.factor(oxford_imperial_data$Alcohol_status)

#PCS score
oxford_imperial_data$PCS_score <- oxford_imperial$PCSTotal

#MNSI_score
oxford_imperial_data$MNSI_score <- oxford_imperial$MNSITotalScore
oxford_imperial_data$MNSI_status <- NA
oxford_imperial_data$MNSI_status[oxford_imperial_data$MNSI_score > 3 & !is.na(oxford_imperial_data$MNSI_score)] <- 1
oxford_imperial_data$MNSI_status[oxford_imperial_data$MNSI_score < 4 & !is.na(oxford_imperial_data$MNSI_score)] <- 0

#DN4
oxford_imperial_data$DN4_score <- as.numeric(substr(oxford_imperial$DN4Score, 2, 2))
oxford_imperial_data$DN4_status <- NA
oxford_imperial_data$DN4_status[oxford_imperial_data$DN4_score > 3 & !is.na(oxford_imperial_data$DN4_score)] <- 1
oxford_imperial_data$DN4_status[oxford_imperial_data$DN4_score < 4 & !is.na(oxford_imperial_data$DN4_score)] <- 0

#Age
oxford_imperial_data$Age <- oxford_imperial$Q1a

#Gender
oxford_imperial_data$Gender <- as.factor(oxford_imperial$Q1b)

#BMI
oxford_imperial_data$Weight <- oxford_imperial$Q1e
oxford_imperial_data$Height <- oxford_imperial$Q1f
oxford_imperial_data$Height[oxford_imperial_data$Height < 0] <- NA
oxford_imperial_data$Weight[oxford_imperial_data$Weight < 0] <- NA

oxford_imperial_data$BMI <- oxford_imperial_data$Weight / (oxford_imperial_data$Height/100)^2

#Diagnosis
oxford_imperial_data$Diagnosis <- as.factor(oxford_imperial$Q75a)
oxford_imperial_data$Diagnosis_Other <- as.factor(oxford_imperial$Q75aOther)
oxford_imperial_data$Diabetes_type <- as.factor(oxford_imperial$Q75bDiabetes)

#HBA1C
HBA1c_oxford_imperial <- data.frame(HBA1C = as.numeric(oxford_imperial$Q75bPerc), HBA1C_mmol_mol = oxford_imperial$Q75bmmol)
HBA1C_model <- lm(HBA1C ~ HBA1C_mmol_mol, data = HBA1c_oxford_imperial)
HBA1c_oxford_imperial$predicted <- round(predict(HBA1C_model, HBA1c_oxford_imperial), digits = 1)
HBA1c_oxford_imperial[is.na(HBA1c_oxford_imperial$HBA1C),]$HBA1C <- HBA1c_oxford_imperial[is.na(HBA1c_oxford_imperial$HBA1C),]$predicted
oxford_imperial_data$HBA1C <- HBA1c_oxford_imperial$HBA1C
sink(paste0(PATH,"HBA1C_model.txt"))
summary(HBA1C_model)
sink()

#Outcome
oxford_imperial_data$Neuropathy <- oxford_imperial$NeuropathyGrading
oxford_imperial_data$Pain <- oxford_imperial$GradingOutcome
oxford_imperial_data$Pain[is.na(oxford_imperial_data$Pain)] <- "No pain"
oxford_imperial_data$Pain <- as.factor(oxford_imperial_data$Pain)
#####################
# Data only in PINS #
#####################
PINS_DOLORISK <- read.csv(paste0(PATH,"raw_data/DOLORISK_PINS_ID.csv"), stringsAsFactors=FALSE, header=FALSE)
load(paste0(PATH, "PINS_DOLORISK_Oxford.RData"))
PINS_DOLORISK <- PINS_DOLORISK[PINS_DOLORISK[,1] %in% oxford_imperial_data$DOLORiskNum,]
PINS_DOLORISK[,2] <- paste0("X", PINS_DOLORISK[,2])
PINS_DOLORISK <- PINS_DOLORISK[PINS_DOLORISK[,2] %in% rownames(diabetes),]

diabetes <- diabetes[PINS_DOLORISK[,2], ]

mapping <- cbind(rownames(oxford_imperial_data[PINS_DOLORISK[,1],]), rownames(diabetes))

oxford_imperial_data$Diabetes_Duration <- NA
oxford_imperial_data$Waist <- NA
oxford_imperial_data$Cholesterol <- NA
oxford_imperial_data$LDL <- NA
oxford_imperial_data$HDL <- NA
oxford_imperial_data$Creatinine <- NA
oxford_imperial_data$TRIGLYCERIDES <- NA          
oxford_imperial_data$DAPOS_ANXIETY <- NA
oxford_imperial_data$DAPOS_DEPRESSION <- NA
#oxford_imperial_data$Depression_quantile <- NA
oxford_imperial_data$Depression_metric <- NA
#oxford_imperial_data$Anxiety_quantile <- NA
oxford_imperial_data$Anxiety_metric <- NA

#Duration of diabetes
date_Interview <- strptime(diabetes$DATEOFINTERVIEW_1, format = "%d/%m/%Y")
date_Diagnosis <- strptime(diabetes$DIABDIAGDATE, format = "%d/%m/%Y")

diff_in_days = difftime(date_Interview, date_Diagnosis, units = "days")
Diabetes_Duration <- round(as.double(diff_in_days)/365)

oxford_imperial_data[PINS_DOLORISK[,1],]$Diabetes_Duration <- Diabetes_Duration

#Waist
oxford_imperial_data[PINS_DOLORISK[,1],]$Waist <- as.numeric(diabetes$WAISTCIRCUMFERENCE) 

#LDL
oxford_imperial_data[PINS_DOLORISK[,1],]$LDL <- as.numeric(diabetes$LPD_LDL_CONC)

#HDL
oxford_imperial_data[PINS_DOLORISK[,1],]$HDL <- as.numeric(diabetes$LPD_HDL_CONC)

#Creatinine
oxford_imperial_data[PINS_DOLORISK[,1],]$Creatinine <- as.numeric(diabetes$BLOOD_CR_VAL)

#TRIGLYCERIDES
oxford_imperial_data[PINS_DOLORISK[,1],]$TRIGLYCERIDES <- as.numeric(diabetes$LPD_TRIG_CONC)

#Cholesterol
oxford_imperial_data[PINS_DOLORISK[,1],]$Cholesterol <- as.numeric(diabetes$LPD_CHOL_CONC)

#DAPOS
oxford_imperial_data[PINS_DOLORISK[,1],]$DAPOS_ANXIETY <- as.numeric(diabetes$"DAPOS ANXIETY")
oxford_imperial_data[PINS_DOLORISK[,1],]$DAPOS_DEPRESSION <- as.numeric(diabetes$"DAPOS DEPRESSION")
oxford_imperial_data$DAPOS_ANXIETY <- as.numeric(oxford_imperial_data$DAPOS_ANXIETY)
oxford_imperial_data$DAPOS_DEPRESSION <- as.numeric(oxford_imperial_data$DAPOS_DEPRESSION)


#CKD
oxford_imperial_data$CKD <- 0
oxford_imperial_data$CKD[grep("kidney", oxford_imperial_data$Diagnosis_Other)] <- 1
oxford_imperial_data$CKD <- as.factor(oxford_imperial_data$CKD)

levels(oxford_imperial_data$CKD)[levels(oxford_imperial_data$CKD)=="1"] <- TRUE 
levels(oxford_imperial_data$CKD)[levels(oxford_imperial_data$CKD)=="0"] <- FALSE 


#Create outcomes
Neuropathy <- oxford_imperial_data$Neuropathy == "Confirmed peripheral neuropathy" | oxford_imperial_data$Neuropathy == "Probable peripheral neuropathy" | oxford_imperial_data$Neuropathy == "Possible peripheral neuropathy" | !is.na(oxford_imperial_data$Neuropathy)
#Neuropathy[oxford_imperial_data$Neuropathy == "Possible peripheral neuropathy"]

Painful <- oxford_imperial_data$Pain == "Definite neuropathic pain" | oxford_imperial_data$Pain == "Probable neuropathic pain" | !is.na(oxford_imperial_data$Pain)
Painless <- oxford_imperial_data$Pain == "No pain" | oxford_imperial_data$Pain == "Unlikely to be neuropathic pain"

oxford_imperial_data$Outcome <- NA
oxford_imperial_data$Outcome[Neuropathy & Painful] <- "Painful_neuropathy"
oxford_imperial_data$Outcome[Neuropathy & Painless] <- "Painless_neuropathy"
oxford_imperial_data$Outcome[!Neuropathy | is.na(Neuropathy)] <- "No_neuropathy"

oxford_imperial_data$Outcome <- ordered(oxford_imperial_data$Outcome, levels= c("No_neuropathy", "Painless_neuropathy", "Painful_neuropathy"))
oxford_imperial_data$Neuropathy <- ordered(oxford_imperial_data$Neuropathy, levels= c("No_neuropathy", "Possible peripheral neuropathy", "Probable peripheral neuropathy", "Confirmed peripheral neuropathy"))
#######################################################
oxford_imperial_data$Center[oxford_imperial_data$Center == 1] <- "Oxford"
oxford_imperial_data$Center[oxford_imperial_data$Center == 2] <- "Imperial"
oxford_imperial_data$Center[oxford_imperial_data$Center == 8] <- "Technion"

###############################################
###### Harmonise anxiety and depression #######
###############################################

explanatory <- c("Outcome",                                            
"Age",                                         
"Gender",                                                                               
"BMI",                                 
"Diabetes_type",                               
"HBA1C",
"PCS_score")

data <- oxford_imperial_data[,c("Anxiety_tscore", "DAPOS_ANXIETY", "Depression_tscore", "DAPOS_DEPRESSION", explanatory)]
data$Outcome <- droplevels(data$Outcome)
data$DAPOS_ANXIETY <- scale(data$DAPOS_ANXIETY)
data$Anxiety_tscore <- scale(data$Anxiety_tscore)
data$DAPOS_DEPRESSION <- scale(data$DAPOS_DEPRESSION)
data$Depression_tscore <- scale(data$Depression_tscore)

fitlm_dapos_anxiety <- lm(DAPOS_ANXIETY ~ 0 + Outcome + Gender + PCS_score, data=data)
summary(fitlm_dapos_anxiety)

fitlm_promis_anxiety <- lm(Anxiety_tscore ~ 0 + Outcome + Gender + PCS_score, data=data)
summary(fitlm_promis_anxiety)

fit_dapos_promis_coef_anxiety <- fix.coef(fitlm_dapos_anxiety, beta = coef(fitlm_promis_anxiety))

summary(fit_dapos_promis_coef_anxiety)

data$DAPOS_ANXIETY_reg <- data$DAPOS_ANXIETY - coef(fit_dapos_promis_coef_anxiety)

fitlm_dapos_anxiety_reg <- lm(DAPOS_ANXIETY_reg ~ 0 + Outcome + Gender + PCS_score, data=data)

summary(fitlm_dapos_anxiety_reg)

# compare proportion explained variance
cor(predict(fitlm_dapos_anxiety), predict(fitlm_dapos_anxiety) + residuals(fitlm_dapos_anxiety))^2
cor(predict(fit_dapos_promis_coef_anxiety), predict(fit_dapos_promis_coef_anxiety) + residuals(fit_dapos_promis_coef_anxiety))^2

fitlm_dapos_depression <- lm(DAPOS_DEPRESSION ~ 0 + Outcome + Gender + PCS_score, data=data)
summary(fitlm_dapos_depression)

fitlm_promis_depression <- lm(Depression_tscore ~ 0 + Outcome + Gender + PCS_score, data=data)
summary(fitlm_promis_depression)

fit_dapos_promis_coef_depression <- fix.coef(fitlm_dapos_depression, beta = coef(fitlm_promis_depression))

summary(fit_dapos_promis_coef_depression)

data$DAPOS_DEPRESSION_reg <- data$DAPOS_DEPRESSION - coef(fit_dapos_promis_coef_depression)

fitlm_dapos_depression_reg <- lm(DAPOS_DEPRESSION_reg ~ 0 + Outcome + Gender + PCS_score, data=data)

summary(fitlm_dapos_depression_reg)

# compare proportion explained variance
cor(predict(fitlm_dapos_depression), predict(fitlm_dapos_depression) + residuals(fitlm_dapos_depression))^2
cor(predict(fit_dapos_promis_coef_depression), predict(fit_dapos_promis_coef_depression) + residuals(fit_dapos_promis_coef_depression))^2

save(file = paste0(PATH, "depression_anxiety_models.RData"), data, fitlm_dapos_anxiety, fitlm_promis_anxiety, fit_dapos_promis_coef_anxiety, fitlm_dapos_depression, fitlm_promis_depression, fit_dapos_promis_coef_depression)

#Anxiety
Anxiety_scores <- cbind(as.numeric(oxford_imperial_data$Anxiety_tscore), as.numeric(oxford_imperial_data$DAPOS_ANXIETY))
oxford_imperial_data$Anxiety_metric[!is.na(Anxiety_scores[,1])] <- "PROMIS_tscore"
oxford_imperial_data$Anxiety_metric[is.na(Anxiety_scores[,1]) & !is.na(Anxiety_scores[,2])] <- "DAPOS"

Anxiety_scores <- scale(Anxiety_scores)
Anxiety_scores[,2] <- Anxiety_scores[,2] - coef(fit_dapos_promis_coef_anxiety)

Anxiety_compare <- data.frame(Score=factor(rep(c("PROMIS_tscore", "DAPOS"), each=nrow(oxford_imperial_data))), Anxiety=as.numeric(c(Anxiety_scores[,1], Anxiety_scores[,2])))
library(plyr)
mu <- ddply(Anxiety_compare, "Score", summarise, grp.mean=mean(Anxiety, na.rm=TRUE))

pdf(paste0(PATH, "Anxiety_scores_density_plot.pdf"))
print(ggplot(Anxiety_compare, aes(x=Anxiety, fill=Score)) + geom_density(alpha=0.4) + geom_vline(data=mu, aes(xintercept=grp.mean, color=Score), linetype="dashed"))
#ggplot(Anxiety_compare, aes(x=Score, y=Anxiety)) + geom_violin(trim=TRUE, scale = "width", aes(fill = Score)) + stat_summary(fun.y=median, geom="point", size=2, color="black") 
dev.off()

Anxiety_scores[is.na(Anxiety_scores[,1]),1] <- Anxiety_scores[is.na(Anxiety_scores[,1]),2]

oxford_imperial_data$Anxiety_tscore <- Anxiety_scores[,1]*attr(Anxiety_scores, 'scaled:scale')[[1]] + attr(Anxiety_scores, 'scaled:center')[[1]]
#Depression
Depression_scores <- cbind(as.numeric(oxford_imperial_data$Depression_tscore), as.numeric(oxford_imperial_data$DAPOS_DEPRESSION))
oxford_imperial_data$Depression_metric[!is.na(Depression_scores[,1])] <- "PROMIS_tscore"
oxford_imperial_data$Depression_metric[is.na(Depression_scores[,1]) & !is.na(Depression_scores[,2])] <- "DAPOS"
Depression_scores <- scale(Depression_scores)
Depression_scores[,2] <- Depression_scores[,2] - coef(fit_dapos_promis_coef_Depression)

Depression_compare <- data.frame(Score=factor(rep(c("PROMIS_tscore", "DAPOS"), each=nrow(oxford_imperial_data))), Depression=as.numeric(c(Depression_scores[,1], Depression_scores[,2])))
library(plyr)
mu <- ddply(Depression_compare, "Score", summarise, grp.mean=mean(Depression, na.rm=TRUE))

pdf(paste0(PATH,"Depression_scores_density_plot.pdf"))
print(ggplot(Depression_compare, aes(x=Depression, fill=Score)) + geom_density(alpha=0.4) + geom_vline(data=mu, aes(xintercept=grp.mean, color=Score), linetype="dashed"))
#ggplot(Depression_compare, aes(x=Score, y=Depression)) + geom_violin(trim=TRUE, scale = "width", aes(fill = Score)) + stat_summary(fun.y=median, geom="point", size=2, color="black") 
dev.off()

Depression_scores[is.na(Depression_scores[,1]),1] <- Depression_scores[is.na(Depression_scores[,1]),2]



oxford_imperial_data$Depression_tscore <- Depression_scores[,1]*attr(Depression_scores, 'scaled:scale')[[1]] + attr(Depression_scores, 'scaled:center')[[1]]

###############################################
save(oxford_imperial_data, file = paste0(PATH, "Oxford_DOLORisk_data.RData"))
#######################################################
missing_values <- data.frame(Variable = names(oxford_imperial_data), Missing_values = (sapply(oxford_imperial_data, function(x) sum(is.na(x))))/nrow(oxford_imperial_data) )

pdf(paste0(PATH, "missing_data_proportions_Oxford.pdf"))
ggplot(data=missing_values, aes(x=reorder(Variable, Missing_values), y=Missing_values)) + geom_bar( stat="identity", , color="blue", fill="white") + ggtitle("Proportions of missing values") + labs(x = "Variable", y = "Missing values proportion") + theme(axis.text.y = element_text(size= 10, face="bold"), axis.title.y = element_text(size=10, face = "bold"), axis.title.x = element_blank(), axis.text.x = element_text(size= 5.5, face = "bold", angle=90, vjust=0.5), title=element_text(size = 12, face="bold"))
dev.off()

missing <- aggr(oxford_imperial_data, combined=TRUE, sortVars=TRUE, sortCombs=TRUE, cex.axis=0.6, cex.lab=1, prop=TRUE, plot=FALSE)

pdf(paste0(PATH, "missing_combinations_Oxford.pdf"))
plot(missing, cex.axis = 0.7, cex.lab=1, only.miss=TRUE, combined=TRUE, oma = c(14,5,5,3))
dev.off()

write.csv(file = paste0(PATH, "oxford_imperial_data_summary.csv"), summary(oxford_imperial_data))

library(Hmisc)

numeric <- unlist(lapply(oxford_imperial_data, is.numeric))  

pdf(paste0(PATH, "oxford_imperial_data_histogram.pdf"))
hist.data.frame(oxford_imperial_data[,numeric])
dev.off()

table(oxford_imperial_data$Outcome, useNA="ifany")

table(oxford_imperial_data$Neuropathy, useNA="ifany")


######################################
### Dundee datasets GS and GODARTS ###
######################################
GoDARTS <- read.csv(paste0(PATH,"raw_data/DOLORisk_GoDARTS_Data_BASELINE.csv"), na.strings=c("","na", "NA", "Not calculated – Missing data", "ND", " Not calculated – Missing data"), stringsAsFactors=FALSE)

GoDARTS_data <- data.frame(matrix(ncol=0,nrow=nrow(GoDARTS)))

#EQ5D Q1
eq5d.scores.df <- data.frame(MO=GoDARTS$Q1a, SC=GoDARTS$Q1b, UA=GoDARTS$Q1c, PD=GoDARTS$Q1d, AD=GoDARTS$Q1e)

GoDARTS$EQ5D_Index <- NA 
GoDARTS$EQ5D_Index[complete.cases(eq5d.scores.df)] <- eq5d(na.omit(eq5d.scores.df), country = "UK", version = "5L", type = "CW")

GoDARTS_data$EQ5D_Index <- GoDARTS$EQ5D_Index
rownames(GoDARTS_data) <- GoDARTS$LinkageID

#PROMIS Anxiety CSV files PIN id to row nr
PROMIS_anxiety <- read.csv(paste0(PATH,"raw_data/PROMIS_anxiety_tscore.csv"), skip=4)
GoDARTS_data$Anxiety_tscore <- NA
GoDARTS_data$Anxiety_tscore[PROMIS_anxiety$PIN] <- PROMIS_anxiety$TScore

#PROMIS Depression 
PROMIS_depression <- read.csv(paste0(PATH,"raw_data/PROMIS_depression_tscore.csv"), skip=4)
GoDARTS_data$Depression_tscore <- NA
GoDARTS_data$Depression_tscore[PROMIS_depression$PIN] <- PROMIS_depression$TScore

#PROMIS Sleep
PROMIS_sleep <- read.csv(paste0(PATH,"raw_data/PROMIS_sleep_tscore.csv"), skip=4)
GoDARTS_data$Sleep_Disturbance_tscore <- NA
GoDARTS_data$Sleep_Disturbance_tscore[PROMIS_sleep$PIN] <- PROMIS_sleep$TScore

GoDARTS_data$Depression_metric <- "PROMIS_tscore"
GoDARTS_data$Anxiety_metric <- "PROMIS_tscore"

#HBA1c model for GS 
GoDARTS_data$HBA1C <- GoDARTS$HBA1C

#TIPI
reverse_order <- function(x){
return(7 - x + 1)
}


GoDARTS_data$TIPIExtraversion <- rowMeans(cbind(GoDARTS$Q7a, reverse_order(GoDARTS$Q7f))) 
GoDARTS_data$TIPIAgreeableness <- rowMeans(cbind(reverse_order(GoDARTS$Q7b), GoDARTS$Q7g)) 
GoDARTS_data$TIPIConscientiousness <- rowMeans(cbind(GoDARTS$Q7c, reverse_order(GoDARTS$Q7h))) 
GoDARTS_data$TIPIEmotionalStability <- rowMeans(cbind(reverse_order(GoDARTS$Q7d), GoDARTS$Q7i)) 
GoDARTS_data$TIPIOpenness <- rowMeans(cbind(GoDARTS$Q7e, reverse_order(GoDARTS$Q7j))) 

#Trauma
GoDARTS_data$Trauma <- NA
GoDARTS_data$Trauma[GoDARTS$Q5 > 1] <- 1 
GoDARTS_data$Trauma[GoDARTS$Q5 == 1 ] <- 0 
GoDARTS_data$Trauma <- as.factor(GoDARTS_data$Trauma)
levels(GoDARTS_data$Trauma)[levels(GoDARTS_data$Trauma)=="1"] <- TRUE 
levels(GoDARTS_data$Trauma)[levels(GoDARTS_data$Trauma)=="0"] <- FALSE 

#Hospital stay
GoDARTS_data$Hospital_stay <- NA
GoDARTS_data$Hospital_stay[GoDARTS$Q6 == 1 ] <- 1 
GoDARTS_data$Hospital_stay[GoDARTS$Q6 == 2 ] <- 0 
GoDARTS_data$Hospital_stay <- as.factor(GoDARTS_data$Hospital_stay)
levels(GoDARTS_data$Hospital_stay)[levels(GoDARTS_data$Hospital_stay)=="1"] <- TRUE 
levels(GoDARTS_data$Hospital_stay)[levels(GoDARTS_data$Hospital_stay)=="0"] <- FALSE 

#Smoking
GoDARTS_data$Ever_smoked_status <- NA
GoDARTS_data$Ever_smoked_status[GoDARTS$Q10 == 4 ] <- 0
GoDARTS_data$Ever_smoked_status[!is.na(GoDARTS$Q10) & GoDARTS$Q10!=4] <- 1
GoDARTS_data$Ever_smoked_status <- as.factor(GoDARTS_data$Ever_smoked_status)
levels(GoDARTS_data$Ever_smoked_status)[levels(GoDARTS_data$Ever_smoked_status)=="1"] <- TRUE 
levels(GoDARTS_data$Ever_smoked_status)[levels(GoDARTS_data$Ever_smoked_status)=="0"] <- FALSE 

#Alcohol
#Keep full alcohol consumption information
GoDARTS_data$Alcohol_consumption <- as.factor(GoDARTS$Q14)
levels(GoDARTS_data$Alcohol_consumption)[levels(GoDARTS_data$Alcohol_consumption)=="6"] <- "Never"
levels(GoDARTS_data$Alcohol_consumption)[levels(GoDARTS_data$Alcohol_consumption)=="5"] <- "Less than 1 day per month"
levels(GoDARTS_data$Alcohol_consumption)[levels(GoDARTS_data$Alcohol_consumption)=="4"] <- "1 to 3 days per month"
levels(GoDARTS_data$Alcohol_consumption)[levels(GoDARTS_data$Alcohol_consumption)=="3"] <- "1 or 2 days per week"
levels(GoDARTS_data$Alcohol_consumption)[levels(GoDARTS_data$Alcohol_consumption)=="2"] <- "3 or 4 days per week"
levels(GoDARTS_data$Alcohol_consumption)[levels(GoDARTS_data$Alcohol_consumption)=="1"] <- "Daily or almost daily"
GoDARTS_data$Alcohol_consumption <- ordered(GoDARTS_data$Alcohol_consumption, levels = c("Never", "Less than 1 day per month", "1 to 3 days per month", "1 or 2 days per week", "3 or 4 days per week", "Daily or almost daily"))

GoDARTS_data$Alcohol_consumption_likert <- GoDARTS_data$Alcohol_consumption
levels(GoDARTS_data$Alcohol_consumption_likert)[levels(GoDARTS_data$Alcohol_consumption_likert)=="Never"] <- 0
levels(GoDARTS_data$Alcohol_consumption_likert)[levels(GoDARTS_data$Alcohol_consumption_likert)=="Less than 1 day per month"] <- 1
levels(GoDARTS_data$Alcohol_consumption_likert)[levels(GoDARTS_data$Alcohol_consumption_likert)=="1 to 3 days per month"] <- 2
levels(GoDARTS_data$Alcohol_consumption_likert)[levels(GoDARTS_data$Alcohol_consumption_likert)=="1 or 2 days per week"] <- 3
levels(GoDARTS_data$Alcohol_consumption_likert)[levels(GoDARTS_data$Alcohol_consumption_likert)=="3 or 4 days per week"] <- 4
levels(GoDARTS_data$Alcohol_consumption_likert)[levels(GoDARTS_data$Alcohol_consumption_likert)=="Daily or almost daily"] <- 5
GoDARTS_data$Alcohol_consumption_likert <- as.numeric(as.character(GoDARTS_data$Alcohol_consumption_likert))

GoDARTS_data$Alcohol_status <- NA
GoDARTS_data$Alcohol_status[GoDARTS_data$Alcohol_consumption_likert == 0] <- FALSE
GoDARTS_data$Alcohol_status[GoDARTS_data$Alcohol_consumption_likert > 0] <- TRUE
GoDARTS_data$Alcohol_status <- as.factor(GoDARTS_data$Alcohol_status)

#PCS score
GoDARTS_data$PCS_score <- rowSums(GoDARTS[,grep("Q16", names(GoDARTS))] - 1)

#MNSI_score
MNSI_forward <- GoDARTS[,c(62:66, 68:71, 73:74)]
MNSI_forward[MNSI_forward==1] <- 1
MNSI_forward[MNSI_forward==2] <- 0
MNSI_reverse <- GoDARTS[,c(67, 72)]
MNSI_reverse[MNSI_reverse==1] <- 0
MNSI_reverse[MNSI_reverse==2] <- 1
MNSI <- cbind(MNSI_forward, MNSI_reverse) 

#Q24. Has your doctor ever told you that you have diabetic neuropathy?
#Q18. Do you ever have any burning pain in your legs and/or feet?

MNSI_pain_feet <- GoDARTS$Q18==1

GoDARTS_data$MNSI_score <- rowSums(MNSI)
GoDARTS_data$MNSI_status <- NA
GoDARTS_data$MNSI_status[GoDARTS_data$MNSI_score > 2 & !is.na(GoDARTS_data$MNSI_score)] <- 1
GoDARTS_data$MNSI_status[GoDARTS_data$MNSI_score < 3 & !is.na(GoDARTS_data$MNSI_score)] <- 0


#DN4
DN4 <- cbind(GoDARTS[,grep("Q34", names(GoDARTS))], GoDARTS[,grep("Q35", names(GoDARTS))])
DN4[DN4==1] <- 1
DN4[DN4==2] <- 0
GoDARTS_data$DN4_score <- rowSums(DN4)
GoDARTS_data$DN4_status <- NA
GoDARTS_data$DN4_status[GoDARTS_data$DN4_score > 2 & !is.na(GoDARTS_data$DN4_score)] <- 1
GoDARTS_data$DN4_status[GoDARTS_data$DN4_score < 3 & !is.na(GoDARTS_data$DN4_score)] <- 0

#Age
GoDARTS_data$Age <- GoDARTS$Age

#Gender
GoDARTS_data$Gender <- as.factor(GoDARTS$Gender)

#BMI
GoDARTS_data$BMI <- GoDARTS$BMI

#HBA1C
GoDARTS_data$HBA1C <- as.numeric(GoDARTS$HBA1C)

#Duration of diabetes
GoDARTS_data$Diabetes_Duration <- GoDARTS$Diabetes_Duration

#LDL
GoDARTS_data$LDL <- GoDARTS$LDL

#HDL
GoDARTS_data$HDL <- GoDARTS$HDL

#Creatinine
GoDARTS_data$Creatinine <- GoDARTS$Creatinine

#TRIGLYCERIDES
GoDARTS_data$TRIGLYCERIDES <- GoDARTS$TRIGLYCERIDES

#Cholesterol
GoDARTS_data$Cholesterol <- GoDARTS$Cholesterol

#CKD
GoDARTS_data$CKD[GoDARTS$CKD == 0] <- FALSE
GoDARTS_data$CKD[GoDARTS$CKD == 1] <- TRUE

#PAIN
GoDARTS_data$Pain <- NA
GoDARTS_data$Pain[(GoDARTS$Q30==1 & !is.na(GoDARTS$Q30)) | (GoDARTS$Q31==1 & !is.na(GoDARTS$Q31))] <- TRUE
GoDARTS_data$Pain[(GoDARTS$Q30==2 & !is.na(GoDARTS$Q30)) & (GoDARTS$Q31==2 & !is.na(GoDARTS$Q31))] <- FALSE

#Create outcomes
GoDARTS_data$Neuropathy <- NA
GoDARTS_data$Neuropathy[GoDARTS_data$MNSI_status == 1 & !is.na(GoDARTS_data$MNSI_status)] <- TRUE
GoDARTS_data$Neuropathy[GoDARTS_data$MNSI_status == 0 & !is.na(GoDARTS_data$MNSI_status)] <- FALSE

GoDARTS_data$Chronic <- GoDARTS$Q32 > 2

GoDARTS_data$DN4_leg_feet <- GoDARTS$Q33b %in% c(10,11)
GoDARTS_data$DN4_leg_feet[is.na(GoDARTS$Q33b)] <- NA

GoDARTS_data$DN4_feet <- GoDARTS$Q33b %in% c(11)
GoDARTS_data$DN4_feet[is.na(GoDARTS$Q33b)] <- NA

GoDARTS_data$Pain_in_leg_feet <- grepl("11", GoDARTS$Q33a) | grepl("10", GoDARTS$Q33a)
GoDARTS_data$Pain_in_leg_feet[is.na(GoDARTS$Q33a)] <- NA

GoDARTS_data$Pain_in_feet <- grepl("11", GoDARTS$Q33a)
GoDARTS_data$Pain_in_feet[is.na(GoDARTS$Q33a)] <- NA

GoDARTS_data$Widespread_pain <- grepl("12", GoDARTS$Q33a)
GoDARTS_data$Widespread_pain[is.na(GoDARTS$Q33a)] <- NA

#Exclusions from cases and controls
#Most bothersome pain where no pain before N=14
GoDARTS_data <- GoDARTS_data[GoDARTS$Q33b %in% GoDARTS$Q33a,]

#Widespread pain N=223
GoDARTS_data <- GoDARTS_data[GoDARTS_data$Widespread_pain==FALSE | is.na(GoDARTS_data$Widespread_pain),]

#No pain but chronic status N=38
GoDARTS_data <- GoDARTS_data[!(GoDARTS_data$Pain==FALSE & !is.na(GoDARTS_data$Chronic)),]

#Pain medication without pain N=59
#GoDARTS$Q31==1 & GoDARTS$Q30==2
##########

GoDARTS_data$Outcome_no_DN4 <- NA
GoDARTS_data$Outcome_no_DN4[!GoDARTS_data$Neuropathy & !is.na(GoDARTS_data$Neuropathy)] <- "No_neuropathy"
GoDARTS_data$Outcome_no_DN4[(GoDARTS_data$Neuropathy & !is.na(GoDARTS_data$Neuropathy) & ((GoDARTS_data$Pain==FALSE & !is.na(GoDARTS_data$Pain)) | (!GoDARTS_data$DN4_feet & !GoDARTS_data$Pain_in_feet & !is.na(GoDARTS_data$DN4_feet) & !is.na(GoDARTS_data$Pain_in_feet))))] <- "Painless_neuropathy"
GoDARTS_data$Outcome_no_DN4[(GoDARTS_data$Neuropathy & !is.na(GoDARTS_data$Neuropathy) & GoDARTS_data$Pain & !is.na(GoDARTS_data$Pain) & (GoDARTS_data$Pain_in_feet | GoDARTS_data$DN4_feet))] <- "Painful_neuropathy"

#write.csv(table(GoDARTS_data$Outcome_no_DN4, useNA = "ifany"), file = "GoDARTS_data$Outcome_no_DN4.csv")

GoDARTS_data$Outcome_legs_feet <- NA
GoDARTS_data$Outcome_legs_feet[!GoDARTS_data$Neuropathy & !is.na(GoDARTS_data$Neuropathy)] <- "No_neuropathy"
GoDARTS_data$Outcome_legs_feet[(GoDARTS_data$Neuropathy & (!GoDARTS_data$Pain | (!GoDARTS_data$DN4_leg_feet & !GoDARTS_data$Pain_in_leg_feet)))] <- "Painless_neuropathy"
GoDARTS_data$Outcome_legs_feet[(GoDARTS_data$Neuropathy & GoDARTS_data$Pain & (GoDARTS_data$Pain_in_leg_feet | GoDARTS_data$DN4_leg_feet))] <- "Painful_neuropathy"

#write.csv(table(GoDARTS_data$Outcome_legs_feet, useNA = "ifany"), file = "GoDARTS_data$Outcome_legs_feet.csv")

GoDARTS_data$Outcome_Chronic <- GoDARTS_data$Outcome_no_DN4
GoDARTS_data[GoDARTS_data$Outcome_Chronic=="Painful_neuropathy" & !is.na(GoDARTS_data$Outcome_Chronic) & (is.na(GoDARTS_data$Chronic) | GoDARTS_data$Chronic==FALSE),]$Outcome_Chronic <- NA

#write.csv(table(GoDARTS_data$Outcome_Chronic, useNA = "ifany"), file = "GoDARTS_data$Outcome_Chronic.csv")

GoDARTS_data$Outcome_DN4_exclude <- GoDARTS_data$Outcome_no_DN4
GoDARTS_data[GoDARTS_data$Outcome_no_DN4=="Painless_neuropathy" & !is.na(GoDARTS_data$Outcome_no_DN4) & GoDARTS_data$DN4_status==1  & !is.na(GoDARTS_data$DN4_status),]$Outcome_DN4_exclude <- NA 
#GoDARTS_data[GoDARTS_data$Outcome_no_DN4=="Painful_neuropathy" & !is.na(GoDARTS_data$Outcome_no_DN4) & GoDARTS_data$DN4_feet==TRUE & !is.na(GoDARTS_data$DN4_feet) & GoDARTS_data$DN4_status==0  & !is.na(GoDARTS_data$DN4_status),]$Outcome_DN4_exclude <- NA 
GoDARTS_data[GoDARTS_data$Outcome_no_DN4=="Painful_neuropathy" & !is.na(GoDARTS_data$Outcome_no_DN4) & GoDARTS_data$DN4_status==0 & !is.na(GoDARTS_data$DN4_status),]$Outcome_DN4_exclude <- NA 


#write.csv(table(GoDARTS_data$Outcome_DN4_exclude, useNA = "ifany"), file = "GoDARTS_data$Outcome_DN4_exclude.csv")

GoDARTS_data$Outcome_DN4_exclude_chronic <- GoDARTS_data$Outcome_DN4_exclude
GoDARTS_data[GoDARTS_data$Outcome_DN4_exclude_chronic=="Painful_neuropathy" & !is.na(GoDARTS_data$Outcome_DN4_exclude_chronic) & (is.na(GoDARTS_data$Chronic) | GoDARTS_data$Chronic==FALSE),]$Outcome_DN4_exclude_chronic <- NA

table(GoDARTS_data$Outcome_DN4_exclude_chronic, MNSI_pain_feet)

#write.csv(table(GoDARTS_data$Outcome_DN4_exclude_chronic, useNA = "ifany"), file = "GoDARTS_data$Outcome_DN4_exclude_chronic.csv")

#boxplot(HBA1C ~ Outcome_DN4_exclude_chronic, GoDARTS_data)

Pain_bothered_most <- GoDARTS$Q33b[GoDARTS_data$Pain & GoDARTS_data$DN4_status==1]

Pain_bothered_most[Pain_bothered_most==1] <- "Back pain pain"
Pain_bothered_most[Pain_bothered_most==2] <- "Neck shoulder pain"
Pain_bothered_most[Pain_bothered_most==3] <- "Facial dental pain"
Pain_bothered_most[Pain_bothered_most==4] <- "Headache"
Pain_bothered_most[Pain_bothered_most==5] <- "Stomach ache or abdominal pain"
Pain_bothered_most[Pain_bothered_most==6] <- "Pain in arms"
Pain_bothered_most[Pain_bothered_most==7] <- "Pain in hands"
Pain_bothered_most[Pain_bothered_most==8] <- "Chest pain"
Pain_bothered_most[Pain_bothered_most==8] <- "Pain in tips"
Pain_bothered_most[Pain_bothered_most==9] <- "Pain in tips"
Pain_bothered_most[Pain_bothered_most==10] <- "Pain in legs or knees"
Pain_bothered_most[Pain_bothered_most==11] <- "Pain in feet"
Pain_bothered_most[Pain_bothered_most==12] <- "Pain throughout the body"
Pain_bothered_most[Pain_bothered_most==13] <- "Other"

write.csv(table(Pain_bothered_most, useNA="ifany"), file = "pain_for_DN4_positive.csv")

write.csv(table(GoDARTS_data$Neuropathy, GoDARTS_data$Pain, useNA="ifany"), file = "cross_tab_neuropathy_pain.csv")

write.csv(table(GoDARTS_data$Neuropathy, (GoDARTS_data$Pain_in_feet | GoDARTS_data$DN4_feet), useNA="ifany"), file = "cross_tab_neuropathy_pain_in_feet.csv")
#######################################################
GoDARTS_data$Center <- "Dundee"
save(GoDARTS_data, file = paste0(PATH, "GoDARTS_data_DOLORisk_data.RData"))
#######################################################
missing_values <- data.frame(Variable = names(GoDARTS_data), Missing_values = (sapply(GoDARTS_data, function(x) sum(is.na(x))))/nrow(GoDARTS_data) )

pdf(paste0(PATH, "missing_data_proportions_GoDARTS.pdf"))
ggplot(data=missing_values, aes(x=reorder(Variable, Missing_values), y=Missing_values)) + geom_bar( stat="identity", , color="blue", fill="white") + ggtitle("Proportions of missing values") + labs(x = "Variable", y = "Missing values proportion") + theme(axis.text.y = element_text(size= 10, face="bold"), axis.title.y = element_text(size=10, face = "bold"), axis.title.x = element_blank(), axis.text.x = element_text(size= 5.5, face = "bold", angle=90, vjust=0.5), title=element_text(size = 12, face="bold"))
dev.off()

missing <- aggr(GoDARTS_data, combined=TRUE, sortVars=TRUE, sortCombs=TRUE, cex.axis=0.6, cex.lab=1, prop=TRUE, plot=FALSE)

pdf(paste0(PATH, "missing_combinations_GoDARTS.pdf"))
plot(missing, cex.axis = 0.7, cex.lab=1, only.miss=TRUE, combined=TRUE, oma = c(14,5,5,3))
dev.off()


#######################################################
### Create follow-up outcome
GoDARTS_follow_up <- read.csv(paste0(PATH,"raw_data/DOLORisk_GoDARTS_Data_FOLLOWUP.csv"), na.strings=c("","na", "NA", "Not calculated – Missing data", "ND", " Not calculated – Missing data"), stringsAsFactors=FALSE)

rownames(GoDARTS_follow_up) <- GoDARTS_follow_up$LinkageID

#MNSI_score Q6-Q18 follow up
GoDARTS_follow_up[,c(33:45)][GoDARTS_follow_up[,c(33:45)] == "Yes"] <- 1
GoDARTS_follow_up[,c(33:45)][GoDARTS_follow_up[,c(33:45)] == "No"] <- 2
MNSI_forward <- GoDARTS_follow_up[,c(33:37, 39:42, 44:45)]
MNSI_forward[MNSI_forward==1] <- 1
MNSI_forward[MNSI_forward==2] <- 0
MNSI_reverse <- GoDARTS_follow_up[,c(38, 43)]
MNSI_reverse[MNSI_reverse==1] <- 0
MNSI_reverse[MNSI_reverse==2] <- 1
MNSI <- cbind(MNSI_forward, MNSI_reverse) 
MNSI <- sapply(MNSI, as.numeric)

GoDARTS_follow_up$MNSI_score <- rowSums(MNSI)
GoDARTS_follow_up$MNSI_status <- NA
GoDARTS_follow_up$MNSI_status[GoDARTS_follow_up$MNSI_score > 2 & !is.na(GoDARTS_follow_up$MNSI_score)] <- 1
GoDARTS_follow_up$MNSI_status[GoDARTS_follow_up$MNSI_score < 3 & !is.na(GoDARTS_follow_up$MNSI_score)] <- 0

#PAIN Q19 - Q20
GoDARTS_follow_up$Pain <- GoDARTS_follow_up$Q19=="Yes" | GoDARTS_follow_up$Q20=="Yes"

#DN4 Q23 - Q24
DN4 <- cbind(GoDARTS_follow_up[,grep("Q23", names(GoDARTS_follow_up))], GoDARTS_follow_up[,grep("Q24", names(GoDARTS_follow_up))])
DN4[DN4=="Yes"] <- 1
DN4[DN4=="No"] <- 0
DN4 <- sapply(DN4, as.numeric)
GoDARTS_follow_up$DN4_score <- rowSums(DN4)
GoDARTS_follow_up$DN4_relevant <- GoDARTS_follow_up$Q19=="Yes" | GoDARTS_follow_up$Q20=="Yes"
GoDARTS_follow_up$DN4_status <- NA
GoDARTS_follow_up$DN4_status[GoDARTS_follow_up$DN4_score > 2 & !is.na(GoDARTS_follow_up$DN4_score)] <- 1
GoDARTS_follow_up$DN4_status[GoDARTS_follow_up$DN4_score < 3 & !is.na(GoDARTS_follow_up$DN4_score)] <- 0

GoDARTS_follow_up$Neuropathy <- NA
GoDARTS_follow_up$Neuropathy[GoDARTS_follow_up$MNSI_status == 1 & !is.na(GoDARTS_follow_up$MNSI_status)] <- TRUE
GoDARTS_follow_up$Neuropathy[GoDARTS_follow_up$MNSI_status == 0 & !is.na(GoDARTS_follow_up$MNSI_status)] <- FALSE

GoDARTS_follow_up$DN4 <- GoDARTS_follow_up$DN4_status 

#Q21
GoDARTS_follow_up$Chronic <- !GoDARTS_follow_up$Q21 %in% c("Less than 1 month", "1-3 months", NA)

GoDARTS_follow_up$DN4_leg_feet <- grepl("Legs", GoDARTS_follow_up$Q22a_b) | grepl("Feet", GoDARTS_follow_up$Q22a_b)
GoDARTS_follow_up$DN4_leg_feet[is.na(GoDARTS_follow_up$Q22a_b)] <- NA

GoDARTS_follow_up$DN4_feet <- grepl("Feet", GoDARTS_follow_up$Q22a_b)
GoDARTS_follow_up$DN4_feet[is.na(GoDARTS_follow_up$Q22a_b)] <- NA

GoDARTS_follow_up$Pain_in_leg_feet <- grepl("Legs", GoDARTS_follow_up$Q22a_a) | grepl("Feet", GoDARTS_follow_up$Q22a_a)
GoDARTS_follow_up$Pain_in_leg_feet[is.na(GoDARTS_follow_up$Q22a_a)] <- NA

GoDARTS_follow_up$Pain_in_feet <- grepl("Feet", GoDARTS_follow_up$Q22a_a)
GoDARTS_follow_up$Pain_in_feet[is.na(GoDARTS_follow_up$Q22a_a)] <- NA

####

#write.csv(table(GoDARTS_follow_up$Neuropathy, GoDARTS_follow_up$Pain, useNA="ifany"), file = "cross_tab_neuropathy_pain_follow_up.csv")

#write.csv(table(GoDARTS_follow_up$Neuropathy, (GoDARTS_follow_up$Pain_in_feet | GoDARTS_follow_up$DN4_feet), useNA="ifany"), file = "cross_tab_neuropathy_pain_in_feet_follow_up.csv")

GoDARTS_follow_up$Outcome_no_DN4 <- NA
GoDARTS_follow_up$Outcome_no_DN4[!GoDARTS_follow_up$Neuropathy & !is.na(GoDARTS_follow_up$Neuropathy)] <- "No_neuropathy"
GoDARTS_follow_up$Outcome_no_DN4[(GoDARTS_follow_up$Neuropathy & !is.na(GoDARTS_follow_up$Neuropathy) & ((GoDARTS_follow_up$Pain==FALSE & !is.na(GoDARTS_follow_up$Pain)) | (!GoDARTS_follow_up$DN4_feet & !GoDARTS_follow_up$Pain_in_feet & !is.na(GoDARTS_follow_up$DN4_feet) & !is.na(GoDARTS_follow_up$Pain_in_feet))))] <- "Painless_neuropathy"
GoDARTS_follow_up$Outcome_no_DN4[(GoDARTS_follow_up$Neuropathy & !is.na(GoDARTS_follow_up$Neuropathy) & GoDARTS_follow_up$Pain & !is.na(GoDARTS_follow_up$Pain) & (GoDARTS_follow_up$Pain_in_feet | GoDARTS_follow_up$DN4_feet))] <- "Painful_neuropathy"

#write.csv(table(GoDARTS_follow_up$Outcome_no_DN4, useNA = "ifany"), file = "GoDARTS_follow_up$Outcome_no_DN4.csv")

GoDARTS_follow_up$Outcome_legs_feet <- NA
GoDARTS_follow_up$Outcome_legs_feet[!GoDARTS_follow_up$Neuropathy & !is.na(GoDARTS_follow_up$Neuropathy)] <- "No_neuropathy"
GoDARTS_follow_up$Outcome_legs_feet[(GoDARTS_follow_up$Neuropathy & (!GoDARTS_follow_up$Pain | (!GoDARTS_follow_up$DN4_leg_feet & !GoDARTS_follow_up$Pain_in_leg_feet)))] <- "Painless_neuropathy"
GoDARTS_follow_up$Outcome_legs_feet[(GoDARTS_follow_up$Neuropathy & GoDARTS_follow_up$Pain & (GoDARTS_follow_up$Pain_in_leg_feet | GoDARTS_follow_up$DN4_leg_feet))] <- "Painful_neuropathy"

#write.csv(table(GoDARTS_follow_up$Outcome_legs_feet, useNA = "ifany"), file = "GoDARTS_follow_up$Outcome_legs_feet.csv")

GoDARTS_follow_up$Outcome_Chronic <- GoDARTS_follow_up$Outcome_no_DN4
GoDARTS_follow_up[GoDARTS_follow_up$Outcome_Chronic=="Painful_neuropathy" & !is.na(GoDARTS_follow_up$Outcome_Chronic) & (is.na(GoDARTS_follow_up$Chronic) | GoDARTS_follow_up$Chronic==FALSE),]$Outcome_Chronic <- NA

#write.csv(table(GoDARTS_follow_up$Outcome_Chronic, useNA = "ifany"), file = "GoDARTS_follow_up$Outcome_Chronic.csv")

GoDARTS_follow_up$Outcome_DN4_exclude <- GoDARTS_follow_up$Outcome_no_DN4
GoDARTS_follow_up[GoDARTS_follow_up$Outcome_no_DN4=="Painless_neuropathy" & !is.na(GoDARTS_follow_up$Outcome_no_DN4) & GoDARTS_follow_up$DN4_status==1  & !is.na(GoDARTS_follow_up$DN4_status),]$Outcome_DN4_exclude <- NA 
#GoDARTS_follow_up[GoDARTS_follow_up$Outcome_no_DN4=="Painful_neuropathy" & !is.na(GoDARTS_follow_up$Outcome_no_DN4) & GoDARTS_follow_up$DN4_feet==TRUE & !is.na(GoDARTS_follow_up$DN4_feet) & GoDARTS_follow_up$DN4_status==0  & !is.na(GoDARTS_follow_up$DN4_status),]$Outcome_DN4_exclude <- NA 
GoDARTS_follow_up[GoDARTS_follow_up$Outcome_no_DN4=="Painful_neuropathy" & !is.na(GoDARTS_follow_up$Outcome_no_DN4) & GoDARTS_follow_up$DN4_status==0 & !is.na(GoDARTS_follow_up$DN4_status),]$Outcome_DN4_exclude <- NA 

#write.csv(table(GoDARTS_follow_up$Outcome_DN4_exclude, useNA = "ifany"), file = "GoDARTS_follow_up$Outcome_DN4_exclude.csv")

GoDARTS_follow_up$Outcome_DN4_exclude_chronic <- GoDARTS_follow_up$Outcome_DN4_exclude
GoDARTS_follow_up[GoDARTS_follow_up$Outcome_DN4_exclude_chronic=="Painful_neuropathy" & !is.na(GoDARTS_follow_up$Outcome_DN4_exclude_chronic) & (is.na(GoDARTS_follow_up$Chronic) | GoDARTS_follow_up$Chronic==FALSE),]$Outcome_DN4_exclude_chronic <- NA

################################
save(GoDARTS_follow_up, file = paste0(PATH, "GoDARTS_followup_DOLORisk.RData"))

################################
load(paste0(PATH, "GoDARTS_data_DOLORisk_data.RData"))

GoDARTS_data$Outcome_no_DN4_2yrs <- NA
GoDARTS_data$Outcome_legs_feet_2yrs <- NA
GoDARTS_data$Outcome_Chronic_2yrs <- NA
GoDARTS_data$Outcome_DN4_exclude_2yrs <- NA
GoDARTS_data$Outcome_DN4_exclude_chronic_2yrs <- NA


GoDARTS_data$MNSI_status_2yrs <- NA
GoDARTS_data[rownames(GoDARTS_follow_up),]$MNSI_status_2yrs <- GoDARTS_follow_up$MNSI_status

GoDARTS_data[rownames(GoDARTS_follow_up),]$Outcome_no_DN4_2yrs <- as.character(GoDARTS_follow_up$Outcome_no_DN4)
GoDARTS_data[rownames(GoDARTS_follow_up),]$Outcome_legs_feet_2yrs <- as.character(GoDARTS_follow_up$Outcome_legs_feet)
GoDARTS_data[rownames(GoDARTS_follow_up),]$Outcome_Chronic_2yrs <- as.character(GoDARTS_follow_up$Outcome_Chronic)
GoDARTS_data[rownames(GoDARTS_follow_up),]$Outcome_DN4_exclude_2yrs <- as.character(GoDARTS_follow_up$Outcome_DN4_exclude)
GoDARTS_data[rownames(GoDARTS_follow_up),]$Outcome_DN4_exclude_chronic_2yrs <- as.character(GoDARTS_follow_up$Outcome_DN4_exclude_chronic)

write.csv(GoDARTS_data, file =paste0(PATH, "GoDARTS_data_baseline_follow_up.csv"))

save(GoDARTS_data, file = paste0(PATH, "GoDARTS_baseline_followp.RData"))
########################################################
GoDARTS_followup_change <-  GoDARTS_data[GoDARTS_data$Outcome_DN4_exclude_chronic_2yrs != GoDARTS_data$Outcome_DN4_exclude_chronic & !(is.na(GoDARTS_data$Outcome_DN4_exclude_chronic_2yrs) | is.na(GoDARTS_data$Outcome_DN4_exclude_chronic)), ]

write.csv(table(GoDARTS_followup_change$Outcome_DN4_exclude_chronic,GoDARTS_followup_change$Outcome_DN4_exclude_chronic_2yrs, useNA="ifany"),  file ="Contigency_table_2yrs_follow_up_GoDARTS_change_only.csv")

save(GoDARTS_followup_change, file = paste0(PATH, "GoDARTS_data_followup_change.RData"))

########################################################
GoDARTS_valid_mnsi_both <- GoDARTS_data[!is.na(GoDARTS_data$MNSI_status) & !is.na(GoDARTS_data$MNSI_status_2yrs),]

write.csv(table(GoDARTS_valid_mnsi_both$Outcome_DN4_exclude_chronic, useNA = "ifany"), file = "GoDARTS_valid_mnsi$Outcome_DN4_exclude_chronic.csv")
write.csv(table(GoDARTS_valid_mnsi_both$Outcome_DN4_exclude_chronic_2yrs, useNA = "ifany"), file = "GoDARTS_valid_mnsi$Outcome_DN4_exclude_chronic_2yrs.csv")

write.csv(table(GoDARTS_valid_mnsi_both$Outcome_DN4_exclude_chronic, GoDARTS_valid_mnsi_both$Outcome_DN4_exclude_chronic_2yrs, useNA = "ifany"), file = "contigency_table_baseline_2yrs_both_valid_MNSI.csv")

#######################################################

GoDARTS_valid_mnsi_followup_change <-  GoDARTS_valid_mnsi_both[GoDARTS_valid_mnsi_both$Outcome_DN4_exclude_chronic_2yrs != GoDARTS_valid_mnsi_both$Outcome_DN4_exclude_chronic & !(is.na(GoDARTS_valid_mnsi_both$Outcome_DN4_exclude_chronic_2yrs) | is.na(GoDARTS_valid_mnsi_both$Outcome_DN4_exclude_chronic)), ]

write.csv(table(GoDARTS_valid_mnsi_followup_change$Outcome_DN4_exclude_chronic, GoDARTS_valid_mnsi_followup_change$Outcome_DN4_exclude_chronic_2yrs, useNA = "ifany"), file = "contigency_table_baseline_2yrs_both_valid_MNSI.csv")



################################

