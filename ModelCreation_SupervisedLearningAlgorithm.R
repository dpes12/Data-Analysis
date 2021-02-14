#Reading Dataset
dat <- read.csv("MalwareSamples10000.csv")
#Setting random seed 
set.seed(1)
#Required Packages
install.packages(c("tidyverse" , "GGally",
                   "caret" , "corrplot" ,"rpart" , "rpart.plot" ,
                   "car" , "dplyr" , "glmnet" , "e1071" , "ipred"))
library(tidyverse)
library(GGally)
library(caret)
library(corrplot)
library(car)
library(dplyr)
library(glmnet)
library(e1071)
library(rpart)  
library(rpart.plot)  
library(ipred) 
view(dat)

#excluding specimenID
dat <- na.omit(dat[,-1])

#Partitioning Data
TrainRow <- createDataPartition(dat$isMalware,
                                p = 0.80,
                                list = FALSE);
#creating training dataset
traindata <- dat[TrainRow,]
#creating test dataset
testdata <- dat[-TrainRow,]

#selecting Random 3 modelling algorithm
set.seed(1)
models.list1 <- c("Logistic Ridge Regression" , 
                  "Logistic LASSO Regression" , 
                  "Logistic Elastic-Net Regression")
models.list2 <- c("Classification Tree" , 
                  "Bagging Tree" , 
                  "Random Forest")
myModels <- c("Binary Logistic Regression" , 
              sample(models.list1 , size = 1) , 
              sample(models.list2 , size = 1))
myModels %>% data.frame

str(dat)
view(myModels)
view(traindata)
dim(traindata)


#Binary Logistic Regression------------------------------------------------------------------------------------------
#binary logistic regression on train data
mod.binarylogisticregression.lg <- glm(isMalware~., ##with sender domain suffix instead of all
                                       family = "binomial" , 
                                       data = traindata);
summary(mod.binarylogisticregression.lg)

#predicted probability on test data
pred.prob <- predict(mod.binarylogisticregression.lg , new = testdata , 
                     type = "response")
pred.class <- ifelse(pred.prob>0.5,"Yes" , "No")



#confusion matrix
cf.lg <- table(pred.class %>% as.factor %>% relevel(ref="Yes"),
               testdata$isMalware %>% as.factor %>% relevel(ref="Yes"));
#prop table on the basis of column
prop <- prop.table(cf.lg , 2); prop %>% round(digit=3)

#summary of confusion matrix
confusionMatrix(cf.lg);

#Performing RFE-----------------------------------------------------------------------------
#setting seed
set.seed(1)


#Creating dummy variables for test data
dummy.testdata <-model.matrix(~. , data = testdata )[, -1] %>% data.frame

#creating dummy variable for train data
dummy.traindata <- model.matrix(~. , data = traindata)[, -1] %>% data.frame


#turning off warning messages
options(warn = -1)
subset <- c(2:10)
view(dummy.traindata)
#setting parameters
cvalidation_par <- rfeControl(functions = lrFuncs , 
                              method = "repeatedcv" , 
                              repeats = 10 ,
                              verbose = FALSE
                              )

#performing rfe
lmProfile <- rfe(traindata$isMalware~.,
                 data = traindata[,1:10],
                 sizes=subset,
                 rfeControl= cvalidation_par)
lmProfile
summary(lmProfile$fit)
#predictive performance
pred.optmod <- predict(lmProfile$fit , newdata = dummy.testdata )
pred.optmod %>% round(digits = 3)

#final predictive probability--------------------------------


pred.prob <- predict(lmProfile$fit,new=dummy.testdata,type="response")
pred.class <- ifelse(pred.prob>0.5,"Yes","No")
#Confusion matrix with re-ordering of "Yes" and "No" responses
cf.lg <- table(pred.class %>% as.factor %>% relevel(ref="Yes"),
               testdata$isMalware %>% as.factor %>% relevel(ref="Yes"));
prop <- prop.table(cf.lg,2); prop %>% round(digit=3) #Proportions by columns


#Summary of confusion matrix
confusionMatrix(cf.lg);


#-----------------------------------------------------------------------------------------
#2nd model Logistic Ridge Regression

#Partitioning Data
TrainRow.lridge <- createDataPartition(dat$isMalware,
                                p = 0.80,
                                list = FALSE);
#creating training dataset
traindata.lridge <- dat[TrainRow.lridge,]
#creating test dataset
testdata.lridge <- dat[-TrainRow.lridge,]

lambdas <- 10^seq(-3,3,length=100)

#setting seed
set.seed(1)
mod.malware.lridge <- train(isMalware ~., #Formula
                         data = traindata.lridge, #Training data
                         method = "glmnet", #Penalised regression modelling
                         #Set to c("center", "scale") to standardise data
                         preProcess = NULL,
                         #Perform 10-fold CV, 5 times over.
                         trControl = trainControl("repeatedcv",
                                                  number = 10,
                                                  repeats = 5),
                         tuneGrid = expand.grid(alpha = 0, #Ridge regression
                                                lambda = lambdas)
)
#Optimal lambda value
mod.malware.lridge$bestTune

#coefficient
coef(mod.malware.lridge$finalModel, mod.malware.lridge$bestTune$lambda)


#predicted probability of death on the test data
pred.class.ridge <- predict(mod.malware.lridge,new=testdata.lridge)

#Confusion matrix with re-ordering of "Yes" and "No" responses
cf.ridge <- table(pred.class.ridge %>% as.factor %>% relevel(ref="Yes"),
                  testdata.lridge$isMalware %>% as.factor %>% relevel(ref="Yes"));

prop <- prop.table(cf.ridge,2); prop %>% round(digit=3) #Proportions

confusionMatrix(cf.ridge)

#----------------------------------------------------------------------------------------------------
#3rd Model Bagging Tree


## --------------------------------------------------------------------------------
set.seed(1)  #Set the random seed.

#Partitioning the Data
trainRowNum.dat <- createDataPartition(dat$isMalware, #The outcome variable
                                      #proportion of data to form the training set
                                      p=0.80,
                                      #Don't store the result in a list
                                      list=FALSE);

train.dat <- dat[trainRowNum.dat,]

test.dat <- dat[-trainRowNum.dat,]


## --------------------------------------------------------------------------------
set.seed(1)
btree.dat <- bagging(isMalware~.,
                    data=train.dat,
                    nbagg=100,  
                    coob=TRUE); 
btree.dat

#Summary of predictions on test set
test.pred.dat <- predict(btree.dat,newdata=test.dat,type="class"); 

test.cf.dat <- confusionMatrix(test.pred.dat %>% relevel(ref="Yes"),
                              test.dat$isMalware %>% relevel(ref="Yes"))
test.cf.dat


## --------------------------------------------------------------------------------
#Intialise the hyperparamter search grid
grid.dat <- expand.grid(nbagg=seq(25,150,25),  #A sequence of nbagg values
                       cp=seq(0,0.5,0.1),  #A sequence of cp values
                       minsplit=seq(5,20,5),  #A sequence of minsplits values
                       #Initialise columns to store the OOB misclassification rate
                       OOB.misclass=NA, 
                       #Initialise columns to store sensitivity, specificity and
                       #accuracy of bagging at each run.
                       test.sens=NA,
                       test.spec=NA,
                       test.acc=NA)  

for (I in 1:nrow(grid.dat))
{
  set.seed(1)
  
  #Perform bagging
  btree.dat <- bagging(isMalware~.,
                      data=train.dat,
                      nbagg=grid.dat$nbagg[I],  
                      coob=TRUE,
                      control=rpart.control(cp=grid.dat$cp[I],
                                            minsplit=grid.dat$minsplit[I]));
  
  #OOB misclassification rate
  grid.dat$OOB.misclass[I] <- btree.dat$err*100
  
  #Summary of predictions on test set
  test.pred.dat <- predict(btree.dat,newdata=test.dat,type="class");  #Class prediction
  
  #Confusion matrix
  test.cf.dat <- confusionMatrix(test.pred.dat %>% as.factor%>% relevel(ref="Yes"),
                                test.dat$isMalware %>% as.factor %>% relevel(ref="Yes"))
  
  prop.cf.dat <- test.cf.dat$table %>% prop.table(2)
  grid.dat$test.sens[I] <- prop.cf.dat[1,1]*100  #Sensitivity
  grid.dat$test.spec[I] <- prop.cf.dat[2,2]*100  #Specificity
  grid.dat$test.acc[I] <- test.cf.dat$overall[1]*100  #Accuracy
}

#Sort the results by the OOB misclassification rate and display them.
grid.dat[order(grid.dat$OOB.misclass,decreasing=FALSE)[1:10],] %>% round(2)

#Perform bagging with optimized hyperparameter
btree_optimised.dat <- bagging(isMalware~.,
                     data=train.dat,
                     nbagg=150,  
                     coob=TRUE,
                     control=rpart.control(cp=0,
                                           minsplit=20));
#summary of prediction
test.pred.dat <- predict(btree_optimised.dat,newdata=test.dat,type="class");#Class prediction
#confusion matrix
test.cf.dat <- confusionMatrix(test.pred.dat %>% as.factor%>% relevel(ref="Yes"),
                               test.dat$isMalware %>% as.factor %>% relevel(ref="Yes"))
test.cf.dat

#-------------------------------------------------------------------------------------------------
#PART 3 Real World testing
#Reading the dataset
emaildata<- read.csv("EmailSamples50000.csv")
set.seed(1)
dim(emaildata)
str(emaildata)
#Removing SpecimenID from the data
emaildata<-na.omit(emaildata[,-1])
dim(emaildata)
#----------------------------------------------------------------
#Binary Logistic Regression
set.seed(1)
options(warn = -1)
dummy.emaildata<-model.matrix(~. ,data=emaildata)[, -1]%>%data.frame
dim(dummy.emaildata)
#pred.optmod <- predict(lmProfile$fit , newdata = emaildata )
#pred.optmod %>% round(digits = 3)

#final predictive probability--------------------------------

pred.prob <- predict(lmProfile$fit,new=dummy.emaildata,type="response")
pred.class <- ifelse(pred.prob>0.5,"Yes" , "No")

#confusion matrix----------------------------------------------------
cf.lg <- table(pred.class %>% as.factor %>% relevel(ref="Yes"),
               emaildata$isMalware %>% as.factor %>% relevel(ref="Yes"));
prop <- prop.table(cf.lg,2); prop %>% round(digit=3)
confusionMatrix(cf.lg);

#Logistic Ridge Regression--------------------------------------------------------
set.seed(1)
#Predicting on Email data
pred.class.email <- predict(mod.malware.lridge,new=emaildata)
#Confusion Matrix
cf.emailridge <- table(pred.class.email %>% as.factor %>% relevel(ref="Yes"),
                  emaildata$isMalware %>% as.factor %>% relevel(ref="Yes"));

prop <- prop.table(cf.ridge,2); prop %>% round(digit=3) #Proportions

confusionMatrix(cf.emailridge)

#Bagging Tree--------------------------------------------------------------
set.seed(1)
#summary of prediction
email.pred <- predict(btree_optimised.dat,newdata=emaildata,type="class");#Class prediction
#confusion matrix
email.cf <- confusionMatrix(email.pred %>% as.factor%>% relevel(ref="Yes"),
                               emaildata$isMalware %>% as.factor %>% relevel(ref="Yes"))
email.cf


