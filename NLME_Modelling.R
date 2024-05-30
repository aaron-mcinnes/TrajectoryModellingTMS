#NLME 

#load data
dat <- read.csv(str_replace(curDir, 'LCMM_Modelling.R', 'processedData.csv'))
dat$txChange <- dat$StartingPHQ9 - dat$endPHQ9
names(dat)[names(dat) == "Week"] <- "Week_Num"

#get classes assigned by LCMM
classDat <- read.csv(str_replace(curDir, 'LCMM_Modelling.R', 'classDat.csv'))
dat <- left_join(dat, dplyr::select(classDat, PID, class, class_factor), by = "PID")

#get initial values from lm coefficients
lm_model <- lm(log(StartingPHQ9) ~ Week_Num, data = dat)
lm_coef <- coef(lm_model)
#Estimate initial values for A, B, and C
A_start <- mean(dat$txChange)
B_start <- -1 / lm_coef[[2]]  # Negative reciprocal of the slope for decay rate
C_start <- mean(dat$endPHQ9)

#nlme params
tolerance = .5
scaling = 0.0000000000000001
iterations = 1000
pnlsIterations = 100
msIterations = 100
EMiterations = 100
nlsTolerance = 100
tolerance = 100

#Full population model
#Compare factors included as RE
model <- nlme::nlme(PHQ9_SCORE ~ txChange * exp(-Week_Num/B) + endPHQ9,
                    data = dat,
                    fixed = txChange + B + endPHQ9 ~ 1,
                    random = txChange + B + endPHQ9 ~ 1 | PID,
                    start = c(txChange = 1, B = 1, endPHQ9 = C_start),
                    control = nlme::nlmeControl(pnlsTol = tolerance, minScale = scaling, 
                                                maxIter = iterations, pnlsMaxIter = pnlsIterations, msMaxIter = msIterations, niterEM = EMiterations, 
                                                opt = "nlminb", nlsTols=nlsTolerance, tolerance = tolerance),
                    verbose = T)
summary(model)

model2 <- nlme::nlme(PHQ9_SCORE ~ txChange * exp(-Week_Num/B) + endPHQ9,
                    data = dat,
                    fixed = txChange + B + endPHQ9 ~ 1,
                    random = txChange + endPHQ9 ~ 1 | PID,
                    start = c(txChange = 1, B = 1, endPHQ9 = C_start),
                    control = nlme::nlmeControl(pnlsTol = tolerance, minScale = scaling, 
                                                maxIter = iterations, pnlsMaxIter = pnlsIterations, msMaxIter = msIterations, niterEM = EMiterations, 
                                                opt = "nlminb", nlsTols=nlsTolerance, tolerance = tolerance),
                    verbose = T)
summary(model2)
anova(model2, model)

# 5fold cross validation
#split Pts into each fold
PIDs <- unique(dat$PID)
folds <- KFold(PIDs, nfolds = 5, stratified = FALSE, seed =1234) #split the ids into 5 folds
model1list <- list()
model2list <- list()
modelStats1 <- data.frame(model = NA, fold = NA, AIC = NA, BIC = NA, logLik = NA, coefTxChange = NA, coefB = NA, coefEndPHQ = NA)
modelStats2 <- modelStats1
for(f in 1:5){
  #create train and test sets
  train_PIDs <- PIDs[-folds[[f]]]
  test_PIDs <- PIDs[folds[[f]]]
  train <- dat[dat$PID %in% train_PIDs,]
  test <- dat[dat$PID %in% test_PIDs,]
  
  #train model with random effects A,B,C
  model1 <- nlme::nlme(PHQ9_SCORE ~ txChange * exp(-Week_Num/B) + endPHQ9,
                       data = train,
                       fixed = txChange + B + endPHQ9 ~ 1,
                       random = txChange + B + endPHQ9 ~ 1 | PID,
                       start = model$coefficients$fixed,
                       control = nlme::nlmeControl(pnlsTol = tolerance, minScale = scaling, 
                                                   maxIter = iterations, pnlsMaxIter = pnlsIterations, msMaxIter = msIterations, niterEM = EMiterations, 
                                                   opt = "nlminb", nlsTols=nlsTolerance, tolerance = tolerance),
                       verbose = T)
  model1list[[f]] <- model1
  trainPred <- cbind(predict(model1), train$PHQ9_SCORE)
  trainCor <- cor.test(trainPred[,1], trainPred[,2])
  modelStats1[f,] <- c("model1", f, summary(model1)$AIC, summary(model1)$BIC, summary(model1)$logLik, c(model1$coefficients$fixed))
  
  
  #train model with random effects A,C
  model2 <- nlme::nlme(PHQ9_SCORE ~ txChange * exp(-Week_Num/B) + endPHQ9,
                       data = train,
                       fixed = txChange + B + endPHQ9 ~ 1,
                       random = txChange + endPHQ9 ~ 1 | PID,
                       start = model$coefficients$fixed,
                       control = nlme::nlmeControl(pnlsTol = tolerance, minScale = scaling, 
                                                   maxIter = iterations, pnlsMaxIter = pnlsIterations, msMaxIter = msIterations, niterEM = EMiterations, 
                                                   opt = "nlminb", nlsTols=nlsTolerance, tolerance = tolerance),
                       verbose = T)
  model2list[[f]] <- model2
  modelStats2[f,] <- c("model2", f, summary(model2)$AIC, summary(model2)$BIC, summary(model2)$logLik, c(model2$coefficients$fixed))
  
}
mStat <- rbind(modelStats1, modelStats2)
mStat$AIC <- as.numeric(mStat$AIC)
mStat$BIC <- as.numeric(mStat$BIC)
mStat$logLik <- as.numeric(mStat$logLik)
mStat$coefTxChange <- as.numeric(mStat$coefTxChange)
mStat$coefB <- as.numeric(mStat$coefB)
mStat$coefEndPHQ <- as.numeric(mStat$coefEndPHQ)
aggregate(. ~ model, mStat[-2], FUN = mean)

#exp decay formula
form1 <- function(A, B, C, t){
  D <- list()
  for(w in 1:length(t)){
    D[[w]] <- A*exp(-w/B)+C
  }
  return(do.call(rbind, D))
}

form2 <- function(D, B){
  C <- list()
  for(w in 1:length(D)){
    C[[w]] <- D[w]-D[1]*exp(-w/B)/
      1-exp(-w/B)
  }
  return(do.call(rbind,C))
}

predFutureScores <- function(D, B, predC){
  preds <- list()
  for(t in 1:6){
    A <- D[1] - predC[t]
    predictedHeldOut <- form1(A, B, predC[t], (t+1):7)
    obsHeldOut <- D[(t+1):7]
    preds[[t]] <- data.frame(maxWeek = t-1, pred = predictedHeldOut, obs = obsHeldOut)
  }
  preds[[7]] <- data.frame(maxWeek = 6, pred = form1(D[1] - predC[7], B, predC[7], 1:7), obs = D)
  return(preds)
}



modelList <- list()
testFolds <- list()
for(f in 1:5){
  #create train and test sets
  train_PIDs <- PIDs[-folds[[f]]]
  test_PIDs <- PIDs[folds[[f]]]
  train <- dat[dat$PID %in% train_PIDs,]
  test <- dat[dat$PID %in% test_PIDs,]
  
  #estimate B at the population level
  modelTrain <- nlme::nlme(PHQ9_SCORE ~ txChange * exp(-Week_Num/B) + endPHQ9,
                       data = train,
                       fixed = txChange + B + endPHQ9 ~ 1,
                       random = txChange + B + endPHQ9 ~ 1 | PID,
                       start = model$coefficients$fixed,
                       control = nlme::nlmeControl(pnlsTol = tolerance, minScale = scaling, 
                                                   maxIter = iterations, pnlsMaxIter = pnlsIterations, msMaxIter = msIterations, niterEM = EMiterations, 
                                                   opt = "nlminb", nlsTols=nlsTolerance, tolerance = tolerance),
                       verbose = T)
  modelList[[f]] <- modelTrain

  predObs <- list()
  predictedScores <- list()
  for(testP in 1:length(test_PIDs)){
    testDat <- test[test$PID == test_PIDs[testP],]
    A <- testDat$txChange[1]
    B <- as.numeric(modelTrain$coefficients$fixed[2])
    C <- testDat$endPHQ9[1]
    t <- 1:nrow(testDat)
    predScores <- form1(A, B, C, t)
    predictedC <- form2(testDat$PHQ9_SCORE, B)
    predPercChange <-  ((predictedC - testDat$PHQ9_SCORE[1])/testDat$PHQ9_SCORE[1])*100
    # predHeldOut <- predFutureScores(testDat$PHQ9_SCORE, B, predictedC)
    # predHeldOut <- do.call(rbind, predHeldOut)
    # predHeldOut$PID <- test_PIDs[testP]
    # predFuture[[testP]] <- predHeldOut
    predObs[[testP]] <- data.frame(fold = f, B = B, PID = testDat$PID, Week_Num = t-1, obsScore = testDat$PHQ9_SCORE, predictedScore = predScores, 
                                   obsC = C, predictedC = predictedC, obsPercChange = testDat$txPercChange, predictedPercChange = predPercChange,
                                   obsResponder = testDat$responderPerc50, class = testDat$class, class_factor = testDat$class_factor)
   
    #predict scores based on predicted C. Get r[CI] and assign to a trajectory
    for(wk in 0:6){
      A <- testDat$StartingPHQ9[1] - predictedC[wk+1] 
      B <- as.numeric(modelTrain$coefficients$fixed[2])
      C <- predictedC[wk+1]
      t = 1:nrow(testDat)
      pred <- as.data.frame(form1(A, B, C, t))
      pred$PID <- testDat$PID
      pred$fold <- f
      pred$maxWeek <- wk
      pred$obs <- testDat$PHQ9_SCORE
      predictedScores[[f]] <- pred
    }
    
  }
  testFolds[[f]] <- do.call(rbind, predObs)
  
  # for each week and fold, get r(observed C vs predicted), AUC (response50%)
  foldDat <- testFolds[[f]]
  for(week2use in 0:6){
    weekDat <- foldDat[foldDat$Week_Num == week2use,]
    weekDat$predResponder <- ifelse(weekDat$predictedPercChange < -50, 1, 0)
    #sens,spec
    pred <- as.factor(weekDat$predResponder)
    actual <- as.factor(weekDat$obsResponder)
    Confusion <- confusionMatrix(pred, actual, positive = "1")
    Sensitivity <- unname(Confusion$byClass["Sensitivity"])
    Specificity <- unname(Confusion$byClass["Specificity"])
    #auc
    responder_roc <- roc(weekDat$obsResponder, weekDat$predResponder)
    responder_auc <- auc(responder_roc)
    responder_auc_ci <- ci.auc(responder_roc)
    #add to DF
    testFolds[[f]]$Sens[testFolds[[f]]$Week_Num == week2use] <- Sensitivity
    testFolds[[f]]$Spec[testFolds[[f]]$Week_Num == week2use] <- Specificity
    testFolds[[f]]$AUC[testFolds[[f]]$Week_Num == week2use] <- responder_auc
    testFolds[[f]]$AUClower[testFolds[[f]]$Week_Num == week2use] <- responder_auc_ci[1]
    testFolds[[f]]$AUCupper[testFolds[[f]]$Week_Num == week2use] <- responder_auc_ci[3]
  }
  
}
testFolds <- do.call(rbind, testFolds)

predAUC <- dplyr::select(testFolds, Week_Num, Sens, Spec, AUC, AUClower, AUCupper)
predicting_future_summary <- aggregate(. ~ Week_Num, predAUC, FUN = mean)
#write.csv(predicting_future_summary, file = paste0(dirname(curDir), '/Predictions/PredictExpDecay_RespBinary', goal, '.csv'), row.names = F)
AUC_plot_expBinary <- ggplot(predicting_future_summary) +
  geom_line(aes(x=Week_Num, y=AUC), stat="identity", size = 1.5) +
  geom_errorbar( aes(x=Week_Num, ymin=AUClower, ymax=AUCupper), width=0.4, colour="black", alpha=0.9, size=1) +
  #ggtitle("Probability Rapid Response")+
  geom_hline(yintercept=0.5, linetype = 'dashed')+
  scale_y_continuous(breaks=c(0.5,0.7,0.9))+
  xlab("Maximum week given to model") + ylab("Area under the curve") +
  #ggtitle("LQ 4 Class") +
  theme_pubr(base_family = "", base_size=14) 
#theme_pubr(base_family = "", legend = "right") +
#theme(plot.title = element_text(hjust = 0.5))
AUC_plot_expBinary
save(AUC_plot_expBinary, file = paste0(dirname(curDir), '/Predictions/AUCcurve_respBinary_ExpDecay'))
