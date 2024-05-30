## functions to run Sarah's trajectory modelling 
#input data (df dat2use), num classes (num nClass), fitting (lin: 1, quad: 2, cubic: 3)

library(lcmm)
library(rBayesianOptimization)
library(MLmetrics)
library(patchwork)
library(ggpubr)

k_value = 6 #number of weeks recorded in data

fitNames <- c("linear", "quadratic", "cubic")

## formula for linear, qudrtatic, cubic
formulaLin <- "PHQ9_SCORE~ Week_Num"
formulaQuad <- "PHQ9_SCORE~ Week_Num + I(Week_Num^2)"
formulaCub <- "PHQ9_SCORE~ Week_Num + I(Week_Num^2) + I(Week_Num^3)"
#random/mixture term 
randLin <- "~Week_Num"
randQuad <- "~Week_Num + I(Week_Num^2)"
randCub <- "~Week_Num + I(Week_Num^2)+ I(Week_Num^3)" 

## store formulas to pass to hlme
formulas <- c(formulaLin, formulaQuad, formulaCub)
REs <- c(randLin, randQuad, randCub)

###Function that gets the correlation of predicted with observed, and returns NA if the SD is 0 
subject_cor_fun<- function(x,y){
  if (any(sapply(c(x,y), FUN = sd) == 0)) {
    return(NA)
  } else {
    val <- cor(x,y)
    return(val)
  }
}

####Function that gets the average of the correlations for each subject
cor_fun <- function(data){
  corrs <- data %>%
    group_by(PID) %>%
    summarise(c = subject_cor_fun(as.data.frame(PHQ9_SCORE), as.data.frame(Pred_wtAvg)))
  corrs_noError <- as.data.frame(corrs[corrs$c < 1.5,])
  return(corrs_noError)
}

#function to run models and get fit statistics
modelFitFun <- function(data, modelName){
  out <- data.frame(Model = modelName,
             MSE = MSE(data$pred_ss, data$obs),
             RMSE = RMSE(data$pred_ss, data$obs),
             R2 = R2_Score(data$pred_ss, data$obs),
             cor = cor(data$pred_ss, data$obs)
  )
  return(out)
}
  
fitModelFun <- function(dat2use, nClass, fit, model2test){
  trainModel <- list()
  smTbl = data.frame(matrix(NA, nrow = nClass, ncol = 8+nClass))
  names(smTbl)[1:9] <- c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "%class1")
  
  #Run model iteratively up to n Classes
  for(model in 1:nClass){
      if(model == 1){
        if(is.list(model2test) == F){
          trainModel[[model]] <- hlme(formula(formulas[fit]), subject = "PID",random=formula(REs[fit]), ng = model, data = dat2use)
        } else if(length(model2test) == nClass){
          trainModel[[model]] <- hlme(formula(formulas[fit]), subject = "PID",random=formula(REs[fit]), ng = model, data = dat2use, nwg=FALSE, B=model2test[[model]]$best, maxiter = 0)
        } else{
          print('training model not supplied')
        }
      } else{
        names(smTbl)[8+model] <- paste0("%class",model)
        if(is.list(model2test) == F){
          trainModel[[model]] <- hlme(formula(formulas[fit]), subject = "PID",random=formula(REs[fit]), nwg = T, data = dat2use, mixture = formula(REs[fit]), ng = model, B = trainModel[[1]])
        }else if(length(model2test) == nClass){
          trainModel[[model]] <- hlme(formula(formulas[fit]), subject = "PID",random=formula(REs[fit]), nwg = T, data = dat2use, mixture = formula(REs[fit]), ng = model, B = model2test[[model]]$best, maxiter = 0)        
        } else{
          print('training model not supplied')
        }
      }
    names(trainModel)[[model]] <- paste(fitNames[fit], paste0(nClass, "class"), model, sep = "_")
    curMod <- trainModel[[model]]
    smTbl[model, 1:(8+model)] <- summarytable(curMod, which = c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "%class"))
    trainModel[[model]]$fit <- fit
  }
  smTbl$Polynomials <- fitNames[fit]
  smTbl$Num_Classes <- 1:nClass
  smTbl <- cbind(str_sub(names(trainModel), 1, -3), smTbl)
  names(smTbl)[1] <- "model"
  
  out <- list(trainModel, smTbl)
  return(out)
}


######################
####Get the class membership info, and predicted values, and add it to the data file
####Calculate the weighted average of the predicted value
######################
evalModelFun <- function(data2use, modelList, smTbl){
  
  trainData <- list()
  trainPred <- list()
  trainPredObs <- list()
  corPredObs <- list()
  classWtAvg <- list()
  fitStat <- list()
  plot <- list()
  for(model in 1:length(modelList)){
    modelName <- names(modelList)[model]
    trainData[[model]] <- merge(data2use, modelList[[model]]$pprob, by = "PID") #merge data with class probability
    trainPred[[model]] <- ddply(modelList[[model]]$pred, c("PID"), transform, Week_Num=1:length(pred_ss)-1) #row numbers are weeks, subtract 1 so that weeks = 0:8
    trainPredObs[[model]] <- merge(trainData[[model]], trainPred[[model]], by = c("PID", "Week_Num")) #merge observed data, class probs, and predictions
    fitStat[[model]] <- modelFitFun(trainPredObs[[model]], modelName) #evaluates fit statistics and stores in list
    for(subc in 1:model){ #calculate the weighted average for each class
      ssCol <- which(names(trainPredObs[[model]]) == paste0("pred_ss", subc))
      probCol <- which(names(trainPredObs[[model]]) == paste0("prob", subc))
      classWtAvg[[subc]] <- trainPredObs[[model]][ssCol] * trainPredObs[[model]][probCol] #multiply prediction by probability and store vector in list
    }
    if(model > 1){ #sums vectors if more than one class
      sumWtAvg <- Reduce("+", classWtAvg)
      names(sumWtAvg) <- "Pred_wtAvg"
    } else{ #if one class, use the initial vector
      sumWtAvg <- as.data.frame(classWtAvg[[subc]])
      names(sumWtAvg) <- "Pred_wtAvg"
    }
    trainPredObs[[model]]$Pred_wtAvg <- sumWtAvg #calculate the wieghted average predicted scores
    names(trainPredObs)[[model]] <- modelName
    corPredObs[[model]] <- cor_fun(as.data.frame(trainPredObs[[model]])) #subject correlations
    smTbl[model, "r_M"] <- mean(corPredObs[[model]]$c, na.rm = T) #add correlation of predicted with observed to summary table
    smTbl[model, "r_SD"] <- sd(corPredObs[[model]]$c, na.rm = T) #add correlation of predicted with observed to summary table
    
    #plot
    plot[[model]] <- ggplot(trainPredObs[[model]], aes(x = Week_Num, #plot scores by class for each model 
                                                       y = PHQ9_SCORE,
                                                       group = factor(class),
                                                       color = factor(class),
                                                       fill = factor(class)))+
      geom_smooth(method = "loess", span = 1)+
      coord_cartesian(ylim = c(5,25))+
      theme_classic()+
      theme(legend.position = "none")
  }
  
  fitStats <- do.call(rbind, fitStat)
  smTbl.fit <- smTbl
  smTbl.error <- fitStats
  
  plots <<- wrap_plots(plot, ncol = 3)
  preds <<- trainPredObs
  out <- list(smTbl.fit, smTbl.error)
  names(out) <- c("Model fit", "Prediction error")
  return(out)
}

#get class membership
classNums <- function(preds){
  group_size = max(preds[[1]]$Week_Num) + 1 #number of weeks + baseline
  classNums <- list()
  for(class in 1:length(preds)){
    classDat <- preds[[class]]
    classesToFit <- 1:class
    groupNums <- classDat %>%
      group_by(class) %>%
      dplyr::summarise(count = n()/group_size)
    groupNums <- cbind(names(preds)[class], groupNums)
    groupNums$label <- paste0(groupNums$class, " (n = ", groupNums$count, ")")
    names(groupNums)[1] <- "model"
    #handle cases where no PIDs are assigned to a class
    classNoFits <- setdiff(classesToFit, groupNums$class)
    if(length(classNoFits > 0)){
      class0 <- data.frame(matrix(NA, nrow = length(classNoFits), ncol = 4))
      names(class0) <- names(groupNums)
      class0$model <- unique(groupNums$model)
      class0$class <- classNoFits
      class0$count <- 0
      class0$label <- paste(classNoFits, "(n = 0)")
      groupNums <- rbind(class0, groupNums)
    }

    classNums[[class]] <- groupNums
  }
  classNums <- do.call(rbind, classNums)
  classNums <-  classNums %>%
    arrange(model, class)
  return(classNums)
}


# plot scores by class and model
classPlots <- function(preds){
  plotList <- list()
  ## plot scores by class for each model
  FullGroupNums <<- classNums(preds)
  groupNums <<- FullGroupNums[FullGroupNums$count != 0,]
  
  if(goal == "blcor"){
    yLimits <- c(-25, 12)
    axisBreaks <- seq(-25, 15, 5)
  } else{
    yLimits = c(-3, 29)
    axisBreaks <- seq(0, 25, 5)
  }
  # Loop through the list of dataframes and create plots
  for (df_name in names(preds)) {
    df <- preds[[df_name]]
    classIdx <- word(groupNums$model[groupNums$model == df_name],3,sep = "_")[1]
    classLab <- ifelse(as.numeric(classIdx) == 1, "Class", "Classes")
    df$class <- factor(df$class, levels = unique(df$class), labels = groupNums$label[groupNums$model == df_name & groupNums$count != 0])
    pTitle <- str_to_title(paste(word(groupNums$model[groupNums$model == df_name],1,sep = "_")[1], 
                                 classIdx, 
                                 classLab))
    
    plot <- ggplot(df, aes(x = Week_Num, y = PHQ9_SCORE, fill = factor(class), color = factor(class), group = factor(class))) + 
      geom_smooth(method = "gam", formula = y ~ s(x, k = k_value)) +
      scale_color_manual(values=custom_palette)+
      scale_fill_manual(values=custom_palette)+
      scale_y_continuous(limits = yLimits, breaks = axisBreaks)+
      labs(x = 'Week', y = 'PHQ9', face = "bold", title = pTitle) +
      theme_pubr(base_family = "", base_size=12) +
      theme(legend.position = c(0.25, .01),
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
    
    plotList[[df_name]] <- plot
  }
  
  #wrap plots
  allPlots <- wrap_plots(plotList, ncol = cols2use)
  return(allPlots)
}

runModel <- function(data2use, numClasses, poly){
  models <- fitModelFun(data2use, numClasses, poly) #data, number of classes, polynomial (1 = linear, 2 = quadratic, 3 = cubic)
  smTbl <- models[[2]]
  models <- models[[1]]
  modelEval <- evalModelFun(data2use, models)
  
  ## plot scores by class for each model
  groupNums <- classNums(preds)
  allPlots <- classPlots(preds)
  
  #wrap plots
  allPlots <- wrap_plots(allPlots, ncol = 2) # You can adjust the number of columns as needed
  allPlots
  
  return(list(modelEval, allPlots))
}



writeToXLsheet <- function(data2write, xlFile, sheetWrite){
  # Check if the file exists; if not, create a new file
  if (!file.exists(xlFile)) {
    write.xlsx(data2write, file = xlFile, sheetName = sheetWrite)
  } else {
    # Load the existing file
    wb <- loadWorkbook(xlFile)
    
    # Check if the desired sheet name already exists
    sheet_names <- getSheetNames(xlFile)
    
    if (sheetWrite %in% sheet_names) {
      # Remove the existing sheet with the same name
      removeWorksheet(wb, sheetWrite)
    }
    
    # Add a new sheet with the desired name and write the data
    addWorksheet(wb, sheetName = sheetWrite)
    writeData(wb, sheet = sheetWrite, x = data2write)
    
    # Save the updated workbook
    saveWorkbook(wb, file = xlFile, overwrite = T)
  }
}
