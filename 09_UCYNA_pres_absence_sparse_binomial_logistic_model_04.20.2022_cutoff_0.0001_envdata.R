#Sparse binomial logistic regression on UCYN-A presence/ absence
#Looking specifically at environmental variables
#Ignore all 18S sequences, including host

#1. Load packages, format data-----
library(tidyverse)
#install.packages("caret")
library(caret)
#install.packages("glmnet")
library(glmnet)
#install.packages("mlbench")
library(mlbench)
install.packages('Rcpp')
library(Rcpp)
library(tidyr)

getwd()
setwd("Desktop/DataAnalysis/SPOT_16S_18S_eLSA/")

#Read in data, use pre-formatted, interpolated abiotic data written in  previous R file
data=read.table("int_CLR_data/int_env_data_UCYNA_Brad_pres_abs_cutoff_0.0001_04.18.2022.tsv", header=T, stringsAsFactors = F, sep="\t", row.names=1)
dim(data) #122 33
names(data) #The last numerical one is "NO2+NO3" and need to exclude sample id as well
#Some of these variables may not be numeric and we need to fix that
#Also exclude model data because the first variable is missing
head(data)

#Get rid of complicated columns
names(data)[grep("Model|Sample", colnames(data))]
data_subs=data[,-c(grep("Model|Sample", colnames(data)))]
dim(data_subs) #122 30

for(i in c(1:ncol(data_subs))){
  print(names(data_subs)[i])
  #print(typeof(data_subs[,i]))
  print(summary(data_subs[,i]))
  #print(length(which(is.na(data_subs[,i]))))
} #All of the numeric variables are "double"s

#There are 6 PO4s and 8 [NO2+NO3]s missing
for(i in c(1:ncol(data_subs))){
  if(length(which(is.na(data_subs[,i])))>0){
    print(names(data_subs)[i])
  }
}

#Get rid of those
names(data_subs)[grep("PO4", names(data_subs))[1]]
names(data_subs)[grep("NO2.NO3", names(data_subs))[1]]
data_subs=data_subs[,-c(grep("PO4", names(data_subs))[1],grep("NO2.NO3", names(data_subs))[1])]
dim(data_subs) #122 28

#Ok rows are samples --check
#Columns are variables (18S ASVs or env variables) --check #Not UCYN-A ASVs relabun, because just want presence/ absence 
#No NAs --check

#Scale 
#Scale just the non-binomial variables 
names(data_subs)
dim(data_subs[,c(1:grep("A1_pres", colnames(data_subs))-1)]) #122 18
data_scaled=as.data.frame(data.matrix(scale(data_subs[,c(1:grep("A1_pres", colnames(data_subs))-1)], scale=T, center=T))) #Change from as.matrix to data.matrix 
colnames(data_scaled)
typeof(data_scaled)
data_scaled[c(1:6), c(1:6)] #Looks good I guess

#Cbind it with the last 5 columns, which have UCYN-A presence/ absence and host presence/ absence
dim(data_subs[,grep("A1_pres", colnames(data_subs)):ncol(data_subs)]) #122 10
names(data_subs[,grep("A1_pres", colnames(data_subs)):ncol(data_subs)])
dim(data_scaled) #123 18
input_data=cbind(data_scaled, data_subs[,grep("A1_pres", colnames(data_subs)):ncol(data_subs)])
dim(input_data) #123 28

#2. Split into training and test set-----
round(0.8*nrow(input_data)) #want 98 rows in training data, 25 rows in test data #Approx
round(0.2*nrow(input_data))

#Pick which years should be in test data (one with low UCYN-A abundance, one with high (pre-post 2013, maybe?))
for(i in c(2008:2018)){
  print(i)
  #print(grep(i, rownames(input_data)))
  print(mean(input_data$A1_dummy[grep(i, rownames(input_data))]))
  print(mean(input_data$A2_dummy[grep(i, rownames(input_data))]))
} #Not a huge variation...

#Pick 2011 as "sparse" year and 2015 as "replete" year
length(grep("2011|2015", rownames(input_data)))
24/nrow(input_data) #0.197 #approximtely 20% 

#Set up training and test data
test.data=input_data[grep("2011|2015", rownames(input_data)),]
dim(test.data) #24 22
train.data=input_data[-grep("2011|2015", rownames(input_data)),]
dim(train.data) #99 22
nrow(test.data)+nrow(train.data)==nrow(input_data) #TRUE

#Define x and y #X= numeric predictor data, Y=binomial response data for Clade 1 and Clade 2
colnames(train.data)
c(grep("A1_dummy", colnames(train.data)):ncol(train.data))
X=train.data[,-c(grep("A1_pres", colnames(train.data)):ncol(train.data))] #gets rid of all presence/absence data for UCYN-A and Brad
dim(X) #98 18

Y1=train.data[,grep("A1_dummy", colnames(train.data))] #Modeling ASV1 presence/ absence
length(Y1) #9
head(Y1)
length(which(Y1==1)) #41 for 0.0001 cutoff
length(which(Y1==1))/length(Y1) #0.41

Y2=train.data[,grep("A2_dummy", colnames(train.data))]
length(Y2) #98
head(Y2)
length(which(Y2==1)) #35
length(which(Y2==1))/length(Y2) #0.357
 
#I. UCYN-A CLADE 1 (ASV1 and ASV6)------------

#4. Find the best lambda through cross validation-----
library(glmnet)
set.seed(42)
cv.lasso <- cv.glmnet(data.matrix(as.data.frame(X)), Y1, alpha = 1, family = "binomial") 
head(cv.lasso)
str(cv.lasso)
cv.lasso$lambda.min #0.026
cv.lasso$lambda.1se #0.107

#Plot
plot(cv.lasso) 

coef(cv.lasso, cv.lasso$lambda.min)
sort(coef(cv.lasso, cv.lasso$lambda.min)) #6 predictors are nonzero
coef(cv.lasso, cv.lasso$lambda.1se)
sort(coef(cv.lasso, cv.lasso$lambda.1se)) #Only CUTI_ix is nonzero (negative predictor) 

#5. Fit model on training data------
model_1se <- glmnet(data.matrix(X), Y1, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.1se)
coef(model_1se)[coef(model_1se)[,1]!=0,] #Just CUTI_ix #-0.5697040 

#6. Assess model accuracy-----
#Figure out what the model predicts
colnames(test.data)
X_test=test.data[,-c(grep("A1_pres", colnames(test.data)):ncol(test.data))]
dim(X_test) #24 18
head(X_test)
#Now plug it into the model to get out predicted UCYN-A1 presence/ absence
PROB_1se=model_1se %>% predict(newx=data.matrix(X_test), type="response") #specify type="response" #Without specifying, this line gives you lodgit scores
head(PROB_1se) 
summary(PROB_1se) 
pres_abs_predict=ifelse(PROB_1se>0.5, 1, 0)
summary(pres_abs_predict)
length(which(pres_abs_predict==1)) #8
length(pres_abs_predict) #24
length(which(pres_abs_predict==1))/length(pres_abs_predict) #0.333333 #Predicted to be present 33% of the time

pres_abs_obs=test.data$A1_dummy
length(pres_abs_obs) #24
length(which(pres_abs_obs==1)) #10
10/24 #Actually present 41.6% of the time
summary(pres_abs_predict)
mean(pres_abs_obs==pres_abs_predict) #0.542

#7. QC model-----
#Confusion matrix
confusionMatrix(data=as.factor(pres_abs_predict), reference=as.factor(pres_abs_obs))
#Sensitivity: 0.667
#Specificity: 0.333
#Kappa: 0

confusionMatrix(data = as.factor(pres_abs_predict), reference = as.factor(pres_abs_obs), mode = "prec_recall")
#Precision: 0.625
#Recall: 0.667
#F1: 0.645

#II. UCYN-A CLADE 2 (ASV5)-------

#4. Find the best lambda through cross validation-----
library(glmnet)
set.seed(42)
cv.lasso_2 <- cv.glmnet(data.matrix(X), Y2, alpha = 1, family = "binomial")
str(cv.lasso_2)
head(cv.lasso_2)
cv.lasso_2$lambda.min #0.0642
cv.lasso_2$lambda.1se #0.0642344

plot(cv.lasso_2) #Uh oh #There is not an upper/ lower bound on this lambda, it just goes straight down
coef(cv.lasso_2) #NONE 
sort(coef(cv.lasso_2)) #No predictors  #Same for cutoff of 0.01   

#5. Fit model on training data-----
model_1se_2 <- glmnet(data.matrix(X), Y2, alpha=1, family="binomial", lambda=cv.lasso_2$lambda.1se)
coef(model_1se_2) #Yeah, no predictors at all
coef(model_1se_2)[coef(model_1se_2)[,1]>0,]
coef(model_1se_2)[coef(model_1se_2)[,1]!=0,] %>% sort(decreasing=T)
length(coef(model_1se_2)[coef(model_1se_2)[,1]>0,] %>% sort(decreasing=T))

#6. Assess model accuracy-----
dim(X_test) #24 18
#Plug X_test into the model to get out predicted UCYN-A1 presence/ absence
PROB_1se_2=model_1se_2 %>% predict(newx=data.matrix(X_test), type="response")
head(PROB_1se_2)
summary(PROB_1se_2) #All the same probability of presence, lol, which is 0.3571

#Convert
pres_abs_predict_2=ifelse(PROB_1se_2>0.5, 1, 0)
head(pres_abs_predict_2)
summary(pres_abs_predict_2) #all 0 

pres_abs_obs_2=test.data$A2_dummy
length(pres_abs_obs_2) #24
length(which(pres_abs_obs_2==1)) #12 
mean(pres_abs_obs_2==pres_abs_predict_2) #0.5 

#7. QC model-----
confusionMatrix(data=as.factor(pres_abs_predict_2), reference=as.factor(pres_abs_obs_2))
#Sensitivity: 1.00
#Specificity: 0.00
#Kappa: 0.00

confusionMatrix(data = as.factor(pres_abs_predict_2), reference = as.factor(pres_abs_obs_2), mode = "prec_recall")
#Kappa:0
#Precision: 0.5 
#Recall: 1.000
#F1: 0.667

#III. BRAARUDOSPHAERA ASVs-----

#2. Format data-----

#a. Euks CLR'd data----
dim(euks_subs) #125 650
#Remove "SPOT.2009.09.24.5m.AE"
grep("SPOT.2009.09.24.5m.AE", colnames(euks_subs)) #Not actually there? Ok! 
head(colnames(euks_subs))
head(rownames(euks_subs))
grep("SPOT.2009.09.01.5m.AE", rownames(euks_subs)) #there it is 
euks_subs=euks_subs[-grep("SPOT.2009.09.01.5m.AE", rownames(euks_subs)),]
dim(euks_subs) #124 650

#Remove Brad ASV7 and Brad ASV4 from the data table 
head(rownames(euks_subs))
head(colnames(euks_subs)) #ASVs 
grep(Brad_ix[7], colnames(euks_subs))
c(grep(Brad_ix[7], colnames(euks_subs)), grep(Brad_ix[4], colnames(euks_subs)))
dim(euks_subs) #124 650
euks_subs=euks_subs[,-c(grep(Brad_ix[7], colnames(euks_subs)), grep(Brad_ix[4], colnames(euks_subs)))]
dim(euks_subs) #124 648

#b. Environmental data----
dim(env_subs) #125 21
head(colnames(env_subs))
tail(colnames(env_subs))
head(rownames(env_subs)) 
grep("SPOT.2009.09.01.5m.AE", rownames(env_subs)) #Remove it
env_subs=env_subs[-grep("SPOT.2009.09.01.5m.AE", rownames(env_subs)),]
dim(env_subs) #124 21

#Remove UCYN-A columns
tail(colnames(env_subs))
env_subs=env_subs[,-c(grep("ASV1", colnames(env_subs)):ncol(env_subs))]
dim(env_subs) #124 16
names(env_subs)

#Add in Brad ix presence/ absence with cutoff of 0.0001 (0.01% of 18S)-----
dim(euks_rel_subs) #28402 124
dim(env_subs) #124 16
euks_rel_subs_t=t(euks_rel_subs)
dim(euks_rel_subs_t) #124 28402 
head(colnames(euks_rel_subs_t)) #ASVs
head(rownames(euks_rel_subs_t)) #Dates

summary(as.numeric(euks_rel_subs_t[,grep(Brad_ix[7], colnames(euks_rel_subs_t))]))
length(which(as.numeric(euks_rel_subs_t[,grep(Brad_ix[7], colnames(euks_rel_subs_t))])>0.0001)) #48 
length(which(as.numeric(euks_rel_subs_t[,grep(Brad_ix[7], colnames(euks_rel_subs_t))])>0)) #48 also

tail(colnames(env_subs))
env_subs$Brad7=0
env_subs$Brad7[which(as.numeric(euks_rel_subs_t[,grep(Brad_ix[7], colnames(euks_rel_subs_t))])>0.0001)]=1

summary(as.numeric(euks_rel_subs_t[,grep(Brad_ix[4], colnames(euks_rel_subs_t))]))
length(which(as.numeric(euks_rel_subs_t[,grep(Brad_ix[4], colnames(euks_rel_subs_t))])>0.0001)) #37
length(which(as.numeric(euks_rel_subs_t[,grep(Brad_ix[4], colnames(euks_rel_subs_t))])>0)) #also 37

env_subs$Brad4=0
env_subs$Brad4[which(as.numeric(euks_rel_subs_t[,grep(Brad_ix[4], colnames(euks_rel_subs_t))])>0.0001)]=1

#c. Combined-----
dim(euks_subs) #124 650
head(colnames(euks_subs)) #ASVs
head(rownames(euks_subs)) #Dates

dim(env_subs) #124 18
head(colnames(env_subs)) #Env data
head(rownames(env_subs)) #Dates 

#Cbind 'em
euks_env=cbind(euks_subs, env_subs)
dim(euks_env) #124 666

#NA omit
euks_env_subs=na.omit(euks_env)
dim(euks_env_subs) #122 666

#Scale just the numerical data
tail(colnames(euks_env_subs))
data_scaled=as.data.frame(data.matrix(scale(euks_env_subs[,-c(grep("Brad", colnames(euks_env_subs)))], center=TRUE, scale=TRUE)))
dim(data_scaled) #122 664
typeof(data_scaled) #"list"
data_scaled[1:6, 1:6]

#Cbind it back with the last two columns
input_data=cbind(data_scaled, euks_env_subs[,grep("Brad", colnames(euks_env_subs))])
dim(input_data) #122 668

#3. Split into test and training set-----
0.8*nrow(euks_env_subs) #97.6
0.2*nrow(euks_env_subs) #24.4 

#Are 2011 and 2015 still replete/ sparse years?? 
year=c()
Brad7=c()
Brad4=c()
for(i in c(2008:2018)){
  year=c(year, i)
  Brad7=c(Brad7, mean(as.numeric(euks_rel_subs[grep(Brad_ix[7], rownames(euks_rel_subs)), grep(i, colnames(euks_rel_subs))])))
  Brad4=c(Brad4, mean(as.numeric(euks_rel_subs[grep(Brad_ix[4], rownames(euks_rel_subs)), grep(i, colnames(euks_rel_subs))])))
}
test=as.data.frame(matrix(ncol=3, nrow=length(c(2008:2018))))
dim(test) #11 3
names(test)=c("Year", "Brad7", "Brad4")
test$Year=year
test$Brad7=Brad7
test$Brad4=Brad4
view(test)
#Yeah, there was no Brad7 or Brad4 at all in 2011 (really?)
#But both are fairly abundant in 2015
mean(as.numeric(proks_subs[grep(UCYNA_ix[1], rownames(proks_subs)), grep("2011", colnames(proks_subs))])) #Also 0
mean(as.numeric(proks_subs[grep(UCYNA_ix[5], rownames(proks_subs)), grep("2011", colnames(proks_subs))])) #Really really low

#Split training and test data
head(rownames(input_data))
length(grep("2011|2015", rownames(input_data))) #24 #Perfect 
24/nrow(input_data) #About 20%
train.data=input_data[-grep("2011|2015", rownames(input_data)),]
dim(train.data) #98 666
test.data=input_data[grep("2011|2015", rownames(input_data)),]
dim(test.data) #24 666

#Define X and Y (binomial response for Brad ASV7 and Brad ASV4)
colnames(train.data)
X=train.data[,-grep("Brad", colnames(train.data))]
dim(X) #98 664

Y7=train.data$Brad7
length(Y7) #98
head(Y7)
length(which(Y7==1)) #38
length(which(Y7==1))/length(Y7)

Y4=train.data$Brad4
length(Y4)
head(Y4)
length(which(Y4==1)) #26
length(which(Y4==1))/length(Y4)

#B. BIGELOWII ASV7 (HOST OF UCYN-A1)------

#4. Find the best lambda through cross validation-----
library(glmnet)
cv.lasso <- cv.glmnet(data.matrix(X), Y7, alpha=1, family="binomial")
head(cv.lasso)
str(cv.lasso)
cv.lasso$lambda.min #0.1222 #0.001: 0.1161
cv.lasso$lambda.1se #0.1858 #0.001: 0.2127

plot(cv.lasso) #between 5 and 2 predictors #0.001: 2 predictors 

coef(cv.lasso, cv.lasso$lambda.min) #CUTI and a few 18S ASVs
coef(cv.lasso, cv.lasso$lambda.1se) #No upwelling indices 
sort(coef(cv.lasso, cv.lasso$lambda.1se)) #Two predictors 

#5. Fit model on training data-----
model_7_1se <- glmnet(data.matrix(X), Y7, alpha=1, family="binomial", lambda=cv.lasso$lambda.1se)
coef(model_7_1se)
coef(model_7_1se)[coef(model_7_1se)[,1]>0,] #Two predictors 
euks_tax$Taxon[grep("d69c8c78557e476eb5a3d0f3ad89615d", euks_tax$Feature.ID)] #Ciliophore
euks_tax$Taxon[grep("ead0d51affbef7dd17d17e483eb0f244", euks_tax$Feature.ID)] #Dinoflagellate 

#6. Assess model accuracy-----
names(test.data)
X_test=test.data[,-c(grep("Brad", colnames(test.data)))]
head(X_test)

PROB_1se_7=model_7_1se %>% predict(newx=data.matrix(X_test), type="response")
head(PROB_1se_7) #Lot of the probabilities are 0.37488

pres_abs_predict=ifelse(PROB_1se_7>0.5, 1, 0)
summary(pres_abs_predict) #All zero 

pres_abs_obs=test.data$Brad7
summary(pres_abs_obs)
length(pres_abs_obs) #24
length(which(pres_abs_obs==1)) #8
length(which(pres_abs_obs==1))/length(pres_abs_obs)
  
summary(pres_abs_obs)
mean(pres_abs_obs==pres_abs_predict) #0.666667 #Just predicting that it's not there 

#B. BIGELOWII ASV4 (HOST OF UCYN-A2)-----

#4. Find the best lambda through cross validation-----
cv.lasso <- cv.glmnet(data.matrix(X), Y4, alpha=1, family="binomial")
head(cv.lasso)
str(cv.lasso)
cv.lasso$lambda.min #0.14079 
cv.lasso$lambda.1se #0.20427 

plot(cv.lasso) #Between 7 and 4 predictors #0.001: Between 21 and 2 predictors 

coef(cv.lasso, lambda=cv.lasso$lambda.min) #No coefficients
coef(cv.lasso, lambda=cv.lasso$lambda.1se)
sort(coef(cv.lasso, lambda=cv.lasso$lambda.1se)) #No coefficients 

#5. Fit model on training data-----
model_4_1se <- glmnet(data.matrix(X), Y4, alpha=1, family="binomial", lambda=cv.lasso$lambda.1se)
coef(model_4_1se)
coef(model_4_1se)[coef(model_4_1se)[,1]>0,] #One predictor #30e5193b6a4eb94e700148d0edacf7ee
#Cutoff of 0.0001: two predictors, c3b3a5a8e14ea2c9cdb58028bc55bca0 and 3cda683a43307a23d529d82a6c614023
euks_tax$Taxon[grep("30e5193b6a4eb94e700148d0edacf7ee", euks_tax$Feature.ID)] #A syndiniales 
euks_tax$Taxon[grep("c3b3a5a8e14ea2c9cdb58028bc55bca0", euks_tax$Feature.ID)] #Chlorophyta
euks_tax$Taxon[grep("3cda683a43307a23d529d82a6c614023", euks_tax$Feature.ID)] #A syndiniales

#6. Assess model accuracy-----
X_test=test.data[,-c(grep("Brad", colnames(test.data)))]
head(X_test)

PROB_1se_4=model_4_1se %>% predict(newx=data.matrix(X_test), type="response")
head(PROB_1se_4) #Lot of the probabilities are 0.32653 

pres_abs_predict=ifelse(PROB_1se_4>0.5, 1, 0)
summary(pres_abs_predict) #All 0 :( 

pres_abs_obs=test.data$Brad4
length(pres_abs_obs)

mean(pres_abs_obs==pres_abs_predict) #0.7916667

