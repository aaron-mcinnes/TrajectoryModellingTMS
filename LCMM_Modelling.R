#TMS trajectories updated by AM 

rm(list = ls())
options(scipen = 999)

#packages
library(stringr)
library(plyr)
library(dplyr)
library(ggplot2)
library(zoo)
library(tidyr)
library(fitdistrplus)
library(rBayesianOptimization)
library(RColorBrewer)
library(pROC)

#load data and functions
curDir <- rstudioapi::getSourceEditorContext()$path
datLoc <- str_replace(curDir, 'latentClassModelling.R', 'TrajectoryData.csv')
#datLoc <- '/Users/aaronmcinnes/Downloads/TMS_Deidentified_2024_04_23.csv'
trajFunPath <- str_replace(curDir, 'LCMM_Modelling.R', 'trajFun.R')
source(trajFunPath)
dat <- as.data.frame(read.csv(str_replace(curDir, 'LCMM_Modelling.R', 'processedData.csv'))) #assumes data is saved in same folder as this script

#seed for random folds
seedNum = 1234
set.seed(seedNum)
PIDs <- unique(dat$PID)

#number of classes to test for each model polynomial
numClasses <- 5
goal = "raw" #blcor or raw

#set blcor or raw
if(goal == "blcor"){
  dat$PHQ9_SCORE <- dat$PHQ9_SCORE - dat$StartingPHQ9
}

#num polynomials (1-3 for linear, quadratic, cubic)
numPoly = 3
polyNames = c("Linear", "Quadratic", "Cubic")

#set colors for classes when plotting
palette_name <- "Set1"  
num_colors <- numClasses
custom_palette <- brewer.pal(n = num_colors, name = palette_name)
#number of columns for plots
cols2use <- 3
plotWidth <- 15 #for saving
plotHeight <- 12 #for saving

#split Pts into each fold
folds <- KFold(PIDs, nfolds = 5, stratified = FALSE, seed =seedNum) #split the ids into 5 folds

#initialise lists
outputPoly <- list()
smTblPoly <- list()
modelsPoly <- list()
outputFold <- list()
smTblFold <- list()
modelsFold <- list()
testOutputPoly <- list()
testOutputFold <- list()

#run each fold
for(foldNum in 1:length(folds)){
  print(paste0('Fold ', foldNum, ' of ', length(folds)))
  #create train and test sets
  train_PIDs <- PIDs[-folds[[foldNum]]]
  test_PIDs <- PIDs[folds[[foldNum]]]
  train <- dat[dat$PID %in% train_PIDs,]
  test <- dat[dat$PID %in% test_PIDs,]
  #iterate over polynomials
  for(polyNum in 1:numPoly){
    print(polyNames[polyNum])
    print('Training model...')
    outputPoly[[polyNum]] <- fitModelFun(train, numClasses, polyNum, NA) #data, number of classes, polynomial (1 = linear, 2 = quadratic, 3 = cubic)
    trainModels <- outputPoly[[polyNum]][[1]]
    modelEvalTrain <- evalModelFun(train, trainModels, outputPoly[[polyNum]][[2]]) 
    #test model
    print('Testing model...')
    testOutputPoly[[polyNum]] <- fitModelFun(test, numClasses, polyNum, trainModels)
    modelEvalTest <- evalModelFun(test, testOutputPoly[[polyNum]][[1]], testOutputPoly[[polyNum]][[2]]) 
    
    #add test predictions to train summary
    testR <- modelEvalTest$`Model fit`[,c('r_M', 'r_SD')]
    names(testR) <- c("r_M_test", "r_SD_test")
    modelEvalTrain$`Model fit` <- cbind(modelEvalTrain$`Model fit`, testR)
    modelEvalTrain$`Model fit`$fold <- foldNum
    
    #store output for poly
    smTblPoly[[polyNum]] <- modelEvalTrain
    modelsPoly[[polyNum]] <- trainModels
  }
  #store output for fold
  outputFold[[foldNum]] <- outputPoly
  testOutputFold[[foldNum]] <- testOutputPoly
  smTblFold[[foldNum]] <- smTblPoly
  modelsFold[[foldNum]] <- modelsPoly
}
#store to disk
save(outputFold, file = str_replace(curDir, 'LCMM_Modelling.R', paste0('/Models/outputFold_', goal)))
save(smTblFold, file = str_replace(curDir, 'LCMM_Modelling.R', paste0('/Models/smTblFold_', goal)))
save(modelsFold, file = str_replace(curDir, 'LCMM_Modelling.R', paste0('/Models/modelsFold', goal)))
save(testOutputFold, file = str_replace(curDir, 'LCMM_Modelling.R', paste0('/Models/testOutputFold', goal)))

## Full summary table (all folds and polys)
sumTbl <- list()
for(fold in 1:length(smTblFold)){
  foldTbl <- smTblFold[[fold]]
  foldTbl <- foldTbl
  
}

# Extract Model fit for each list
sumTblLists <- lapply(smTblFold, function(x) {
  lapply(x, function(y) {
    y$`Model fit`
  })
})
sumTbl <- list()
for(fold in 1:length(sumTblLists)){
  curList <- do.call(rbind, sumTblLists[[fold]])
  sumTbl <- rbind(sumTbl, curList)
}

#average across folds and calculate delta SABIC
Model_Summary <- as.data.frame(sumTbl %>% 
                                 group_by(fold, Polynomials) %>% 
                                 mutate(BIC_diff_class = lag(BIC)-BIC, SABIC_diff_class = lag(SABIC)-SABIC, AIC_diff_class = lag(AIC)-AIC))

Model_Summary <- as.data.frame(Model_Summary %>% 
                                 group_by(fold, Num_Classes) %>% 
                                 mutate(BIC_diff_poly = lag(BIC)-BIC, SABIC_diff_poly = lag(SABIC)-SABIC, AIC_diff_poly = lag(AIC)-AIC))
Model_Summary <- as.data.frame(Model_Summary %>%  
                                 group_by(fold, Num_Classes, Polynomials) %>% 
                                 mutate(Min_Class_Percent = min(c(`%class1`, `%class2`, `%class3`, `%class4`, `%class5`), na.rm = TRUE)))

Model_Summary$Polynomials <- factor(Model_Summary$Polynomials, levels = c("linear", "quadratic", "cubic"))
Model_Summary_final <- ddply(Model_Summary, c("Polynomials", "Num_Classes"), summarise, 
                             BIC=mean(BIC), BIC_diff_poly=mean(BIC_diff_poly),
                             BIC_diff_class=mean(BIC_diff_class),
                             SABIC=mean(SABIC), SABIC_diff_poly=mean(SABIC_diff_poly),
                             SABIC_diff_class=mean(SABIC_diff_class),
                             train_r=mean(r_M), test_r=mean(r_M_test),
                             Min_Class_Percent=mean(Min_Class_Percent))
Model_Summary_final[,3:ncol(Model_Summary_final)] <- round(Model_Summary_final[,3:ncol(Model_Summary_final)], 2)
Model_Summary_final <- Model_Summary_final[,names(Model_Summary_final) != "train_r"]
Model_Summary_final
write.csv(Model_Summary_final, file = paste0(dirname(curDir), '/ModelSummary_', goal, '.csv'))


#############################################
# Run selected model (LQC) on full dataset
names(dat)[names(dat) == 'Week'] <- 'Week_Num'
modelOut <- fitModelFun(dat, numClasses, 3, NA)
save(modelOut, file = paste0(dirname(curDir), '/FullSet_ModelSummary_', goal))

smTbl <- modelOut[[2]]
models <- modelOut[[1]]

modelEval <- evalModelFun(dat, models, smTbl)
modelEval
modelPlots <- classPlots(preds)
tiff(paste0(dirname(curDir), '/Models/Plots/FullSet_', goal, '.tiff'),width = 12, height = 7.5, units = 'in', res = 300)  
modelPlots
dev.off()  # Close the TIFF file
classGroupNums <- FullGroupNums


##for raw
if(goal == "raw"){
  predDatRaw <- preds$cubic_5class_4
  nClasses <- 4
} else if(goal == "blcor"){
  predDatRaw <- preds$cubic_5class_2
  nClasses <- 2
  
  aggregate(PHQ9_SCORE ~ class, predDatRaw[predDatRaw$Week_Num == 6,], FUN = mean)
  aggregate(PHQ9_SCORE ~ class, predDatRaw[predDatRaw$Week_Num == 6,], FUN = sd)
}
unqPred <- unique(dplyr::select(predDatRaw,PID, class, StartingPHQ9, endPHQ9, txPercChange, responderPerc50, Tx_Year))
unqPred %>% group_by(class) %>%
  summarise(count_respondPerc_1 = sum(responderPerc50 == 1, na.rm = TRUE),
            n= n())

if(goal == "raw"){
  #get slope of class 2 (minimal improvement) and 3 (non response)
  class2lm <- lm(PHQ9_SCORE ~ Week_Num, predDatRaw[predDatRaw$class == 2,])
  coef(class2lm)
  class3lm <- lm(PHQ9_SCORE ~ Week_Num, predDatRaw[predDatRaw$class == 3,])
  coef(class3lm)
  allclasslm <- lm(PHQ9_SCORE ~ Week_Num*class, predDatRaw)
  anova(allclasslm)
}

#plot
if(goal == 'raw'){
  classLabs <- c("Non-response (n=63)", "Gradual improvement (n=38)", "Rapid improvement (n=25)", "Minimal improvement (n=112)")
  predDatRaw$class_factor <- factor(predDatRaw$class, levels = c(3, 4, 1, 2), labels = classLabs)
  tiff(paste0(dirname(curDir), '/Models/Plots/LQC4class_', goal, '.tiff'),width = 4, height = 4, units = 'in', res = 300)  
  plot <- ggplot(predDatRaw, aes(x = Week_Num, y = PHQ9_SCORE, fill = factor(class_factor), color = factor(class_factor), group = factor(class_factor))) + 
    geom_smooth(method = "gam", formula = y ~ s(x, k = k_value)) +
    scale_color_manual(values=custom_palette)+
    scale_fill_manual(values=custom_palette)+
    coord_cartesian(ylim = c(0, 25))+
    labs(x = 'Week', y = 'PHQ9', face = "bold", title = "LQC 4 Class") +
    theme_pubr(base_family = "", base_size=12) +
    theme(legend.position = c(0.7, .01),
          legend.justification = c("right", "bottom"),
          legend.box.just = "right",
          legend.margin = margin(-2, 2, 2, 2),
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.background=element_blank(),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.key.size = unit(0.75, "lines"),
          axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10),
          plot.title = element_text(hjust = 0.5, vjust = -5, size = 12)
    )    
  plot
  dev.off()
  save(plot, file = paste0(dirname(curDir), '/Models/Plots/LQC4class_', goal))
} else if(goal == 'blcor'){
  classLabs <- c("Improvement (n=97)", "Non-response (n=141)")
  predDatRaw$class_factor <- factor(predDatRaw$class, levels = c(1, 2), labels = classLabs)
  tiff(paste0(dirname(curDir), '/Models/Plots/LQC2class_', goal, '.tiff'),width = 4, height = 4, units = 'in', res = 300)  
  plot <- ggplot(predDatRaw, aes(x = Week_Num, y = PHQ9_SCORE, fill = factor(class_factor), color = factor(class_factor), group = factor(class_factor))) + 
    geom_smooth(method = "gam", formula = y ~ s(x, k = k_value)) +
    scale_color_manual(values=custom_palette)+
    scale_fill_manual(values=custom_palette)+
    coord_cartesian(ylim = c(min(dat$PHQ9_SCORE), 5))+
    labs(x = 'Week', y = 'PHQ9', face = "bold", title = "LQC 2 Class, Baseline-Corrected") +
    theme_pubr(base_family = "", base_size=12) +
    theme(legend.position = c(0.5, .05),
          legend.justification = c("right", "bottom"),
          legend.box.just = "right",
          legend.margin = margin(-2, 2, 2, 2),
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.background=element_blank(),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.key.size = unit(0.75, "lines"),
          axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10),
          plot.title = element_text(hjust = 0.5, vjust = -5, size = 12)
    )    
  plot
  dev.off()
  save(plot, file = paste0(dirname(curDir), '/Mode