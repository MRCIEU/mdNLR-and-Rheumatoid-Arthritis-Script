library(survival)
setwd("~/mDNLR project/Rheumatoid Arthritis/mdNLR_RA_GSE42681")
data<-read.csv("pData_GSE42861.csv", header=TRUE)
data<-read.csv("cellcount_pDataGSE42861_nooutliers.csv",header=T)
# 2 samples with very high mdNLR were removed from the analysis leaving 687 samples
min(data$mdNLR)
max(data$mdNLR)
res.by <- by(data$mdNLR, data$group, mean)
res.by
data_nor<- data[data$group == 'Normal',]
mdNLR_nor<- data_nor$mdNLR
IQR(mdNLR_nor) # 1.045757
# Multiply by 1.5 to detect outliers (values above 1.568636 are outliers?)
1.045757*1.5

data_ra<- data[data$group == 'Rheumatoid arthritis',]
mdNLR_ra<- data_ra$mdNLR
IQR(mdNLR_ra) # 2.694616
# Multiply by 1.5 to detect outliers (values above 4.061022 are outliers?)
2.694616*1.5

# Testing if mdNLR is normally distributed 
ks.test(mdNLR_ra, mdNLR_nor)
shapiro.test(data$mdNLR)
hist(data$mdNLR, breaks=250) 
hist(mdNLR_nor, breaks=250)
hist(mdNLR_ra, breaks=250)
plot(density(mdNLR_nor));plot(density(mdNLR_ra))
shapiro.test(mdNLR_nor); shapiro.test(mdNLR_ra)

## Plot using a qqplot
qqnorm(mdNLR_nor);qqline(mdNLR_nor, col = 2)
qqnorm(mdNLR_ra);qqline(mdNLR_ra, col = 2)

hist(merged.data$mdNLR)

# Boxplot
# boxplot(mdNLR~group,data=data, main="mdNLR Data", xlab="group", ylab="Inflammation index")
library(ggplot2)

test<-data[,c(13,11)]
head(test)
library(reshape2)
test.m <- melt(test)
library(ggplot2)
tiff("boxplot_RACaCo_mdNLR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
P= ggplot(test.m, aes(x=variable, y=value, fill = group)) +
  geom_boxplot(outlier.colour = "black", outlier.size = NA) +
  scale_y_continuous(limits = quantile(test.m$value, c(0.1, 0.99)))+
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Case Control status") +
  ylab("mdNLR")+
  ggtitle("")
P+ labs(x = "Case Control status")
dev.off()



# Logistic regression analysis 
# Univariate model
model <- glm(group ~ mdNLR,family=binomial(link='logit'),data=data)
summary(model)
# Multivariate model with age, gender and smoking status
model1 <- glm(group ~ mdNLR+age+gender+smokingstatus,family=binomial(link='logit'),data=data)
summary(model1)
# Multivariate model with age, gender,smoking status and batch
model2 <- glm(group ~ mdNLR+age+gender+smokingstatus+as.factor(sentrixid),family=binomial(link='logit'),data=data)
summary(model2)
coef(model2)
confint(model2)
# Interpreting the results 
anova(model, test="Chisq")
anova(model1, test="Chisq")
anova(model2, test="Chisq")
library(pscl)
pR2(model)
pR2(model1)
pR2(model2)

# Summary table
#gender
table(data$group)
table(data$gender, data$group)
table(data_nor$gender)
table(data_ra$gender)
# age
res.by<-by(data$age,data$group, mean)
res.by
res.by<-by(data$age,data$group, range)
res.by

# Wilcox test
age_nor<-data$age[data$group=="Normal"]
age_ra<-data$age[data$group=="rheumatoid arthritis"]
wilcox.test(age_nor, age_ra)

mean(data_nor$age)
mean(data_ra$age)
# smoking status
table(data$smokingstatus, data$group)

table(data_nor$smokingstatus)
table(data_ra$smokingstatus)
###############################################
# Smoking NA samples removed from analysis
###############################################
setwd("~/mDNLR project/Rheumatoid Arthritis/mdNLR_RA_GSE42681")
data<-read.csv("cellcount_pDataGSE42861_nooutliers_nosmkna.csv",header=T)
# Logistic regression analysis 
# Univariate model
model <- glm(group ~ mdNLR+sva$n.sv,family=binomial(link='logit'),data=data)
summary(model)
# Multivariate model with age, gender and smoking status
model1 <- glm(group ~ mdNLR+age+gender+smokingstatus,family=binomial(link='logit'),data=data)
options(scipen = 6)# Open scientific calculator
summary(model1)
# Multivariate model with age, gender,smoking status and batch
model2 <- glm(group ~ mdNLR+age+gender+smokingstatus+as.factor(sentrixid),family=binomial(link='logit'),data=data)
summary(model2)
coef(model2)
confint(model2)
options("scipen"=100, "digits"=6)
options(scipen = 999) # Close scientific calculator
exp(cbind(Odds=coef(model1), confint(model1)))
exp(cbind(Odds=coef(model2), confint(model2)))


# Interpreting the results 
anova(model, test="Chisq")
anova(model1, test="Chisq")
anova(model2, test="Chisq")
# Results do not change
#################################################


#################################################
res.by <- by( data$mdNLR, data$group, mean)
res.by

mdNLR_nor<-data$mdNLR[data$group== "Normal"]
mdNLR_ra<-data$mdNLR[data$group=="Rheumatoid arthritis"]
options( scipen = 0)
wilcox.test(mdNLR_nor, mdNLR_ra)

res.by <- by( data$age, data$group, mean)
res.by

Age_RA<-data$age[data$group== "Rheumatoid arthritis"]

Age_control<-data$age[data$group== "Normal"]

wilcox.test(Age_RA, Age_control, alternative = "g")

res.by <- by( data$age, data$group, min)
res.by

res.by <- by( data$age, data$group, max)
res.by

table(data_ra$gender)
table(data_nor$gender)

Gender_RA <-
  matrix(c(101, 251, 96, 239),
         nrow = 2,
         dimnames =
           list(c("Males", "Females"),
                c("Rheumatoid arthritis", "Normal")))
chisq.test(Gender_RA) 

table(data_ra$smokingstatus)
table(data_nor$smokingstatus)


Smoking_RA <- matrix(c(89,108,101,35,2, 110,120,91,31,0), 5, 2,
                dimnames = list(income = c("current", "former", "never", "occasional","na"),
                                satisfaction = c("Control", "RA")))  
fisher.test(Smoking_RA)


mean(data_nor$mdNLR)
sd(data_nor$mdNLR)
range(data_nor$mdNLR)

mean(data_ra$mdNLR)
sd(data_ra$mdNLR)
range(data_ra$mdNLR)

mdNLR_nor<-data_nor$mdNLR
mdNLR_ra<-data_ra$mdNLR
t.test(mdNLR_nor,mdNLR_ra)
 
t.test(Age_control,Age_RA)




library(ggplot2)

df <- data.frame(f1=factor(rbinom(100, 1, 0.45), label=c("m","w")), 
                 f2=factor(rbinom(100, 1, 0.45), label=c("young","old")),
                 boxthis=rnorm(100))

df$f1f2 <- interaction(df$f1, df$f2)

ggplot(aes(y = boxthis, x = f1f2), data = df) + geom_boxplot()

# Boxplot of mdNLR based on Ethnicity
boxplot(mdNLR~ ethnicity.x ,data=SABRE_merge_base, main="mdNLR ethnicity", 
        xlab="mdNLR", ylab="Ethnicity")


# add the cell proportions to the pData
#head(pData)
#pData<-pData[1:18]
#pData_Houseman=cbind(pData,wbc2$counts)
#View(pData_Houseman)
# plot the proportions by group of interest
test<-data[,c(8:10,13)]
head(test)
colnames(test) <- c("Monocytes", "Granulocytes","Lymphocytes", "group")
library(reshape2)
test.m <- melt(test)
library(ggplot2)
tiff("boxplot_lymphocyteRA.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(test.m, aes(x=variable, y=value, fill = group)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 2) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
    xlab("Leucocyte subtype") +
  ylab("Estimated proportion in whole bloood")+
  ggtitle("")
dev.off()

test1<-data[,c(11,13)]
ggplot(test1, aes(mdNLR, colour = group)) +
  geom_density() +
xlim(0, 10)
###################################

# ROC surve for cases-control RA

##################################

# Known risk factors

fit1 <- glm(group ~ age+gender+smokingstatus,family=binomial(link='logit'),data=data)

# mdNLR
fit2 <- glm(group ~ mdNLR,family=binomial(link='logit'),data=data)

# mdNLR with known risk factors
fit3 <- glm(group ~ mdNLR+age+gender+smokingstatus+as.factor(sentrixid),family=binomial(link='logit'),data=data)


preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)

preds3=predict(fit3)
roc3=roc(fit3$y , fit3$fitted.values,ci=T)


roc.plot(roc1,col='black', CI = TRUE)
ggroc(roc1,alpha = 0.5, colour = "red", linetype = 1, size = 1)

tiff("ROC_modelcomparison.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggroc(list("Covariates only"=roc1, "mdNLR only"=roc2, "mdNLR+Covariates"=roc3), aes="colour", size=2)
dev.off()
##########################################################
# mdNLR and RA treatment response
# Differential Methylation as a Biomarker of Response to Etanercept in Patients With Rheumatoid Arthritis
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4914881/
##########################################################  
setwd("~/mDNLR project/Rheumatoid Arthritis/sri_2017")
data_Plant<- read.csv("RA_treatment_DP_060917_CRPNA_removed.csv",header=T)
library(readxl)
data_Plant<- read_excel("mdNLR.xlsx")# Includes all samples for Suppl. Table 1
colnames(data_Plant)


res.by <- by( data_Plant$mdNLR,data_Plant$Sample_Group, mean)
res.by
res.by <- by( data_Plant$mdNLR,data_Plant$Sample_Group, sd)
res.by

mdNLR_good<-data_Plant$mdNLR[data_Plant$Sample_Group=="good"]
mdNLR_poor<-data_Plant$mdNLR[data_Plant$Sample_Group== "poor"]
wilcox.test(mdNLR_good, mdNLR_poor)


res.by <- by( data_Plant$cg00901982,data_Plant$Sample_Group, mean)
res.by

res.by <- by( data_Plant$cg25938803,data_Plant$Sample_Group, mean)
res.by

res.by <- by( data_Plant$cg01591037,data_Plant$Sample_Group, mean)
res.by

res.by <- by( data_Plant$cg03621504,data_Plant$Sample_Group, mean)
res.by

res.by <- by( data_Plant$cg10456459,data_Plant$Sample_Group, mean)
res.by

# Sample Group and Myeloid CpG

# Wilcox test (cg00901982)
cg00901982_good<-data_Plant$cg00901982[data_Plant$Sample_Group=="good"]
cg00901982_poor<-data_Plant$cg00901982[data_Plant$Sample_Group=="poor"]
wilcox.test(cg00901982_good, cg00901982_poor)

library("ggpubr")
tiff("corrplot_cg00901982_mdNLR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggscatter(data_Plant, x = "mdNLR",  y = "cg00901982" , add = "reg.line", conf.int = TRUE,ylim = c(0,1), 
          cor.coef = TRUE, cor.method = "spearman",cor.coef.coord = c(6,1),cor.coef.size = 5,
          xlab = "mdNLR", ylab = "DNA methylation")
dev.off()


# Wilcox test (cg25938803)
cg25938803_good<-data_Plant$cg25938803[data_Plant$Sample_Group=="good"]
cg25938803_poor<-data_Plant$cg25938803[data_Plant$Sample_Group=="poor"]
wilcox.test(cg25938803_good,cg25938803_poor)

# Wilcox test (cg01591037)
cg01591037_good<-data_Plant$cg01591037[data_Plant$Sample_Group=="good"]
cg01591037_poor<-data_Plant$cg01591037[data_Plant$Sample_Group=="poor"]
wilcox.test(cg01591037_good,cg01591037_poor)


# Wilcox test (cg03621504)
cg03621504_good<-data_Plant$cg03621504[data_Plant$Sample_Group=="good"]
cg03621504_poor<-data_Plant$cg03621504[data_Plant$Sample_Group=="poor"]
wilcox.test(cg03621504_good,cg03621504_poor)


# Wilcox test (cg10456459)
cg10456459_good<-data_Plant$cg10456459[data_Plant$Sample_Group=="good"]
cg10456459_poor<-data_Plant$cg10456459[data_Plant$Sample_Group=="poor"]
wilcox.test(cg10456459_good,cg10456459_poor)



# Wilcox test (age)
age_good<-data_Plant$ageatbaseline[data_Plant$Sample_Group=="good"]
age_poor<-data_Plant$ageatbaseline[data_Plant$Sample_Group=="poor"]
wilcox.test(age_good, age_poor)

# Wilcox test (das score)
dascore_good<-data_Plant$dascore[data_Plant$Sample_Group=="good"]
dascore_poor<-data_Plant$dascore[data_Plant$Sample_Group=="poor"]
wilcox.test(dascore_good, dascore_poor, alternative = "g")

# Wilcox test (age)
haqScore_good<-data_Plant$haqScore[data_Plant$Sample_Group=="good"]
haqScore_poor<-data_Plant$haqScore[data_Plant$Sample_Group=="poor"]
wilcox.test(haqScore_good, haqScore_poor)

##########################################################################
# Use of DMARDs and Myeloid CpG (None were significant)

cg10456459_dmardy<-data_Plant$cg10456459[data_Plant$curDmard=="yes"]
cg10456459_dmardn<-data_Plant$cg10456459[data_Plant$curDmard=="no"]
wilcox.test(cg10456459_dmardy,cg10456459_dmardn)

cg03621504_dmardy<-data_Plant$cg03621504[data_Plant$curDmard=="yes"]
cg03621504_dmardn<-data_Plant$cg03621504[data_Plant$curDmard=="no"]
wilcox.test(cg03621504_dmardy,cg03621504_dmardn)

cg01591037_dmardy<-data_Plant$cg01591037[data_Plant$curDmard=="yes"]
cg01591037_dmardn<-data_Plant$cg01591037[data_Plant$curDmard=="no"]
wilcox.test(cg01591037_dmardy,cg01591037_dmardn)

cg25938803_dmardy<-data_Plant$cg25938803[data_Plant$curDmard=="yes"]
cg25938803_dmardn<-data_Plant$cg25938803[data_Plant$curDmard=="no"]
wilcox.test(cg25938803_dmardy,cg25938803_dmardn)

cg00901982_dmardy<-data_Plant$cg00901982[data_Plant$curDmard=="yes"]
cg00901982_dmardn<-data_Plant$cg00901982[data_Plant$curDmard=="no"]
wilcox.test(cg00901982_dmardy,cg00901982_dmardn)

##########################################################################



library("ggpubr")
tiff("corrplot_CRP_mdNLR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggscatter(data_Plant, x = "crp" , y = "mdNLR", add = "reg.line", conf.int = TRUE,ylim = c(0,10), 
          cor.coef = TRUE, cor.method = "spearman",cor.coef.coord = c(50,10),cor.coef.size = 5,
          xlab = "CRP (mg/L)", ylab = "mdNLR")
dev.off()

data<- read.csv("RA_treatment_DP_060917_ESRNA_removed.csv",header=T)
colnames(data_Plant)
tiff("corrplot_ESR_mdNLR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggscatter(data_Plant, x = "esr" , y = "mdNLR", add = "reg.line", conf.int = TRUE,ylim = c(0,10), 
          cor.coef = TRUE, cor.method = "spearman",cor.coef.coord = c(100,10),cor.coef.size = 5,
          xlab = "ESR", ylab = "mdNLR")
dev.off()


# Load the full dataset
data<- read.csv("RA_treatment_DP_060917.csv",header=T)

res.by <- by( data$age, data$Sample_Group, mean)
res.by

res.by <- by( data$age, data$Sample_Group, range)
res.by

table( data$Gender, data$Sample_Group)

table( data$smokestatbase, data$Sample_Group)

res.by <- by(data$DAS28, data$Sample_Group, mean, na.rm=TRUE)
res.by

res.by <- by(data$DAS28, data$Sample_Group, sd, na.rm=TRUE)
res.by

res.by <- by(data$haqScore, data$Sample_Group, mean, na.rm=TRUE)
res.by

res.by <- by(data$haqScore, data$Sample_Group, sd, na.rm=TRUE)
res.by

res.by <- by( data_Plant$mdNLR, data_Plant$curDmard, mean)
res.by

res.by <- by( data_Plant$DAS28, data_Plant$Sample_Group, mean, na.rm=TRUE)
res.by


res.by <- by( data_Plant$DAS28, data_Plant$Sample_Group, sd, na.rm=TRUE)
res.by



library(ggplot2)

test<-data[,c(4,21)]
head(test)
library(reshape2)
test.m <- melt(test)
library(ggplot2)
tiff("boxplot_RA_treatment_mdNLR.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(test.m, aes(x=variable, y=value, fill = Sample_Group)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 2) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Treatment response") +
  ylab("mdNLR")+
  ggtitle("")
dev.off()

res.by <- by( data$esr, data$Sample_Group, mean, na.rm=TRUE)
res.by
# ESR higher in poor responders

res.by <- by( data$crp, data$Sample_Group, mean, na.rm=TRUE)
res.by
# CRP higher in poor responders
test<-data[,c(4,20)]
head(test)
library(reshape2)
test.m <- melt(test)
library(ggplot2)
tiff("boxplot_RA_treatment_CRP.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(test.m, aes(x=variable, y=value, fill = Sample_Group)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 2) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Treatment response") +
  ylab("mdNLR")+
  ggtitle("")
dev.off()


# Logistic regression analysis 
# Univariate model
model <- glm(Sample_Group ~ mdNLR,family=binomial(link='logit'),data=data)
summary(model)
# Multivariate model with age, gender and smoking status
# Baseline disease activity, ,  were included in the model as covariates. 
model1 <- glm(Sample_Group ~ mdNLR+ageatbaseline+Gender+curDmard+dascore,family=binomial(link='logit'),data=data)
summary(model1)
# Multivariate model with age, gender,smoking status and batch
model2 <- glm(Sample_Group ~ mdNLR+ageatbaseline+Gender+curDmard+dascore,family=binomial(link='logit'),data=data)
summary(model2)

exp(cbind(Odds=coef(model2), confint(model2)))


model3 <- glm(Sample_Group ~ mdNLR+ageatbaseline+smokestatbase+Gender+curDmard+dascore+Batch,family=binomial(link='logit'),data=data)
summary(model3)

exp(cbind(Odds=coef(model3), confint(model3)))

# Used for the final analysis
model4 <- glm(Sample_Group ~ mdNLR+ageatbaseline+smokestatbase+Gender+curDmard+DAS28+Batch,family=binomial(link='logit'),data=data)
summary(model4)

exp(cbind(Odds=coef(model4), confint(model4)))


# Smoking NA removed

data<- read.csv("RA_treatment_DP_060917_SMKNA_removed.csv",header=T)
model3 <- glm(Sample_Group ~ mdNLR+smokestatbase+ageatbaseline+Gender+curDmard+dascore+Batch,family=binomial(link='logit'),data=data)
summary(model3)


model4 <- glm(Sample_Group ~ mdNLR+smokestatbase+ageatbaseline+Gender+curDmard+DAS28+Batch,family=binomial(link='logit'),data=data)
summary(model4)

exp(cbind(Odds=coef(model4), confint(model4)))


# # Load the full dataset
data<- read.csv("RA_treatment_DP_060917.csv",header=T)
model <- glm(curDmard ~ mdNLR,family=binomial(link='logit'),data=data)
summary(model)

model <- glm(curDmard ~ mdNLR+ageatbaseline+Gender+dascore+Sample_Group,family=binomial(link='logit'),data=data)
summary(model)


lm1<-lm(mdNLR~ Sample_Group, data=data)
summary(lm1)


lm1<-lm(mdNLR~ Sample_Group+ageatbaseline+Gender+curDmard+dascore, data=data)
summary(lm1)

### Association between DAS28 score and mdNLR

# +ageatbaseline+as.factor(Gender)+curDmard+dascore+Batch
lm1<-lm(mdNLR~ DAS28+Sample_Group+ageatbaseline+as.factor(Gender)+curDmard+dascore+Batch, data=data_Plant)
summary(lm1)

# mdNLR and dascore
lm1<-lm(mdNLR~ dascore+Sample_Group+ageatbaseline+as.factor(Gender)+curDmard+dascore+Batch, data=data_Plant)
summary(lm1)

# mdNLR and haqscore
lm1<-lm(mdNLR~ haqScore+Sample_Group+ageatbaseline+as.factor(Gender)+curDmard+dascore+Batch, data=data_Plant)
summary(lm1)

# mdNLR and DMARD use
lm1<-lm(mdNLR~ curDmard, data=data_Plant)
summary(lm1)

# mdNLR and esr
lm1<-lm(mdNLR~ esr, data=data_Plant)
summary(lm1)

# mdNLR and CRP
lm1<-lm(mdNLR~ crp+esr+dascore+DAS28+ageatbaseline+Gender+curDmard+Batch, data=data_Plant) # Associated
summary(lm1)

#####################################################################################
# Derived the mdNLR score by using IDOL function by Devin Koestler
# Done this on Blue Crystal
#####################################################################################
#####################################################################################
# Estimation of methylation-derived NLR (mdNLR) 
# Devin C. Koestler
# 02/16/2017
#####################################################################################

#####################################################################################
# FUNCTION:  EstimateMdNLR
#    This function estimates the methylation-derived NLR (mdNLR)
#    blood-derived DNA methylation data
#
# ARGUMENTS:
#    Y:              Data frame or matrix (J x N) of methylation beta values for the target 
#                    data set, i.e., blood-based DNAm data on N subjects across
#                    J CpG loci.  The rownames of Y should represent the CpG names; i.e, 
#                    "cg12345678" 
#    dir:            path to the directory where the ProjectionMatrix.RData file is located.  If not
#                    specified, defaults to the current working directory.  
#    verbose:        A logical indicating whether the number of CpGs that overlap between the supplied data set, Y, 
#                    and the number of mdNLR DMRs is to be reported.  Defaults to FALSE
#                    
# RETURNS:   A vector of estimated mdNLR for each of the N subjects in Y
#
# NOTES: You will need to load in the 
#####################################################################################

EstimateMdNLR = function(Y, dir = NULL, verbose = FALSE) {
  
  if(is.null(dir)) {
    load("ProjectionMatrix.RData")
  }
  else {
    file = paste(dir, "/ProjectionMatrix.RData", sep = "")
    load(file)
  }
  
  projectWBCnew = function(Y, coefWBC, contrastWBC=NULL, nonnegative=TRUE){ 
    if(is.null(contrastWBC)) Xmat = coefWBC
    else Xmat = coefWBC %*% t(contrastWBC) 
    
    nCol = dim(Xmat)[2]
    nSubj = dim(Y)[2]
    
    mixCoef = matrix(0, nSubj, nCol)
    rownames(mixCoef) = colnames(Y)
    colnames(mixCoef) = colnames(Xmat)
    
    if(nonnegative){
      library(quadprog)
      
      Amat = cbind(rep(-1,nCol), diag(nCol))
      b0vec = c(-1,rep(0,nCol))
      
      for(i in 1:nSubj){
        obs = which(!is.na(Y[,i])) 
        Dmat = t(Xmat[obs,])%*%Xmat[obs,]
        mixCoef[i,] = solve.QP(Dmat, t(Xmat[obs,])%*%Y[obs,i], Amat, b0vec, meq = 0)$sol
      }
    }
    
    else {
      for(i in 1:nSubj){
        obs = which(!is.na(Y[,i])) 
        Dmat = t(Xmat[obs,])%*%Xmat[obs,]
        mixCoef[i,] = solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
      }
    }
    return(mixCoef)
  }
  
  mdNLRCpGs = rownames(ProjectionMatrix)
  int = intersect(rownames(Y), mdNLRCpGs)
  if(verbose) {
    overlap = length(int)
    print(paste("Of the 228 mdNLR DMRs, ", overlap, " out of 228 were contained in the supplied data set", 
                sep = ""))
  }
  targetData = Y[int,]
  projData = ProjectionMatrix[int,]
  pred = projectWBCnew(targetData, projData)
  mdNLR = pred[,"gran_means"]/pred[,"lymph_means"]
  mdNLR
}

mdNLR = EstimateMdNLR(betas)

###################################################

# Correlation between mdNLR and mdNLR IDOL

###################################################
# prepare the mdNLR IDOL data and merge
head(mdNLR)
mdNLR$GSM_ID<-row.names(mdNLR)

# Merge the 2 datasets data and mdNLR
head(data)

colnames(mdNLR) <- c("mdNLR_IDOL", "GSM_ID")

merged.data <- merge(data, mdNLR, by.x="gsmid", by.y="GSM_ID")

write.csv(merged.data, "data_RA_GSE42861_mdNLR_IDOL.csv")

merged.data<-read.csv("data_RA_GSE42861_mdNLR_IDOL.csv")

library("ggpubr")
tiff("corrplot_mdNLR_mdNLR_IDOL.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggscatter(merged.data, x = "mdNLR" , y = "mdNLR_IDOL",color = "group",palette =c("blue", "red"), add = "none", conf.int = TRUE,xlim = c(0,15), ylim = c(0,15), 
          cor.coef = TRUE, cor.method = "pearson",cor.coef.coord = c(14,15),cor.coef.size = 5,
          xlab = "mdNLR", ylab = "mdNLR_IDOL")
dev.off()


tiff("corrplot_mdNLR_DAS28_IDOL.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggscatter(data_Plant, x = "mdNLR" , y = "DAS28", add = "none", conf.int = TRUE,xlim = c(0,15), ylim = c(0,15), 
          cor.coef = TRUE, cor.method = "pearson",cor.coef.coord = c(14,15),cor.coef.size = 5,
          xlab = "mdNLR", ylab = "DAS28")
dev.off()


# QC for the samples
pvals<-detectionP(rgSet, type = "m+u")
annotation <- meffil.get.features("450k")

pvalue_over_0.05 <- pvals > 0.05
count_over_0.05 <- rowSums(sign(pvalue_over_0.05))
Probes_to_exclude_Pvalue <- rownames(pvals)[which(count_over_0.05 > ncol(pvals)*0.05)]
XY <- as.character(annotation$name[which(annotation$chromosome %in% c("chrX", "chrY"))])
SNPs.and.controls <- as.character(annotation$name[-grep("cg|ch", annotation$name)])
annotation<- annotation[-which(annotation$name %in% c(XY,SNPs.and.controls,Probes_to_exclude_Pvalue)),]

preprocessQuantile(rgSet, fixOutliers = TRUE, removeBadSamples = TRUE,
                   badSampleCutoff = 10.5, quantileNormalize = TRUE,
                   stratified = TRUE, mergeManifest = FALSE, sex = NULL,
                   verbose = TRUE)



# Extract all the information required from the rgset data
pvals<-detectionP(rgset)
save(pvals, file="pvals_rgset_200617.RData")
### Cell Type Composition
library("FlowSorted.Blood.450k")
cellcount<-estimateCellCounts(rgset,compositeCellType="Blood")
rm(FlowSorted.Blood.450k)
##### Normalize the data###
mSet<-preprocessIllumina(rgset,bg.correct=T,normalize="controls")
rm(rgset)
gc()
mset<-preprocessSWAN(rgset, mSet=mSet,verbose=T)
############ Exclude samples and CpG sites that performed poorly 
##### Samples with >5% missing CpG sites (detection p-value>0.01) will be excluded
failed<-pvals>0.01
subQuality_sample<-colSums(failed)
perc_sample<-subQuality_sample/nrow(failed)*100
mset<-mset[,which(perc_sample<=5)]

##### CpG sites with > 20% samples with missing values (detection p-value>0.01)  will be excluded
subQuality_cpg<-rowSums(failed)
perc_cpg<-subQuality_cpg/ncol(failed)*100
mset<-mset[which(perc_cpg<=20),]

############ Remove CH
# remove cross reactive probes
cross_reactive <- read.csv("Crossreactive_probes.csv",header=F)$V1
as.character(cross_reactive)
mset<-mset[!featureNames(mset)%in%cross_reactive,]
# remove SNP probes
SNP <- read.csv("SNP_EUR_race_5percent.csv",header=F)$V1
as.character(SNP)
mset_SNP<-mset[!featureNames(mset)%in%SNP,]

########### Exclude Y chromosome
info <- getAnnotation(mset_SNP) 
rowWant  <- !is.element(info$chr, "chrY")
mset_fil<- mset_SNP[rowWant,]
mset_fil = updateObject(mset_fil)
save(mset_fil, file="mset_filtered_902_020417.RData")



cellcounts<-estimateCellCounts(rgSet, compositeCellType = "Blood",
                   processMethod = "auto", probeSelect = "auto",
                   cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"),
                   referencePlatform = c("IlluminaHumanMethylation450k"),
                   returnAll = FALSE, meanPlot = FALSE, verbose = TRUE)

write.csv(cellcounts, "cellcounts_RAGSE42861.csv")

mset_SWAN<-preprocessSWAN(rgSet, mSet = NULL, verbose = FALSE)

save(mset_SWAN, file="mset_SWAN_RA170917.RData")



############ Exclude samples and CpG sites that performed poorly 
##### Samples with >5% missing CpG sites (detection p-value>0.01) will be excluded
pvals<-detectionP(rgSet)
failed<-pvals>0.01
subQuality_sample<-colSums(failed)
perc_sample<-subQuality_sample/nrow(failed)*100
mset<-mset_SWAN[,which(perc_sample<=5)]

##### CpG sites with > 20% samples with missing values (detection p-value>0.01)  will be excluded
subQuality_cpg<-rowSums(failed)
perc_cpg<-subQuality_cpg/ncol(failed)*100
mset<-mset[which(perc_cpg<=20),]

############ Remove CH
# remove cross reactive probes
cross_reactive <- read.csv("Crossreactive_probes.csv",header=F)$V1
as.character(cross_reactive)
mset<-mset[!featureNames(mset)%in%cross_reactive,]
# remove SNP probes
SNP <- read.csv("SNP_EUR_race_5percent.csv",header=F)$V1
as.character(SNP)
mset_SNP<-mset[!featureNames(mset)%in%SNP,]

########### Exclude Y chromosome
info <- getAnnotation(mset_SNP) 
rowWant  <- !is.element(info$chr, c("chrY", "chrX"))
mset_fil<- mset_SNP[rowWant,]
mset_fil = updateObject(mset_fil)
save(mset_fil, file="mset_filtered_RA_170917.RData")



pData<-read.csv("pDataGSE42861.csv", header=T)
betas<-betas(mset_fil)
meth<-beta2m(betas)
meth[!is.finite(meth)] <- 0

cov<-data.frame(Sample_Group=pData$group,Smoke_Status=pData$smokingstatus,Age_Blood=pData$age,Gender=pData$gender,
                Sentrix_ID=pData$sentrixid,Sentrix_Position=pData$sentrixpos, CD8T=pData$CD8T, CD4T=pData$CD4T, NK=pData$NK, Bcell=pData$Bcell, Mono=pData$Mono, Gran=pData$Gran)


tiff("pcrplot_RA_GSE42861.tif", res=300, compression = "lzw", height=5, width=5, units="in")

pcrplot(meth, cov, npc=10)

dev.off()

# load sva from the C drive
df<-as.data.frame(sva$sv)
write.csv(df, "df_sva_RA.csv")
# In excel we merged the data of SVs and pData used for carrying out SVA analysis "pDataGSE42861.csv" and removed 2 samples with na in smk status
# load the merged data and carry out the analysis

###############################################
# Smoking NA samples removed from analysis
###############################################
setwd("~/mDNLR project/Rheumatoid Arthritis/mdNLR_RA_GSE42681")
data<-read.csv("cellcount_pDataGSE42861_nooutliers_nosmkna_SVA.csv",header=T)
# Logistic regression analysis 
# Univariate model
model <- glm(group ~ mdNLR,family=binomial(link='logit'),data=data)
summary(model)
#with sva (remains significant)
model1 <- glm(group ~ mdNLR+as.numeric(SV1)+as.numeric(SV2)+as.numeric(SV3)+as.numeric(SV4)+as.numeric(SV5)+as.numeric(SV6)+as.numeric(SV7)+as.numeric(SV8)+as.numeric(SV9)+as.numeric(SV10),family=binomial(link='logit'),data=data)
summary(model1)

# Multivariate model with age, gender and smoking status
#model2 <- glm(group ~ mdNLR+age+gender+smokingstatus,family=binomial(link='logit'),data=data)
#options(scipen = 6)# Open scientific calculator
#summary(model1)
# Multivariate model with age, gender,smoking status and SVs
model3 <- glm(group ~ mdNLR+age+gender+smokingstatus+as.numeric(SV1)+as.numeric(SV2)+as.numeric(SV3)+as.numeric(SV4)+as.numeric(SV5)+as.numeric(SV6)+as.numeric(SV7)+as.numeric(SV8)+as.numeric(SV9)+as.numeric(SV10),family=binomial(link='logit'),data=data)
summary(model3)
options("scipen"=100, "digits"=6)
options(scipen = 999) # Close scientific calculator
exp(cbind(Odds=coef(model1), confint(model1)))
exp(cbind(Odds=coef(model3), confint(model3)))

# Adding CpGs instead of mdNLR and observe betas for the model and OR
#
model4 <- glm(group ~ cg00901982+age+gender+smokingstatus+as.numeric(SV1)+as.numeric(SV2)+as.numeric(SV3)+as.numeric(SV4)+as.numeric(SV5)+as.numeric(SV6)+as.numeric(SV7)+as.numeric(SV8)+as.numeric(SV9)+as.numeric(SV10),family=binomial(link='logit'),data=merge)
summary(model4)
exp(cbind(Odds=coef(model4), confint(model4)))

# Interpreting the results 
anova(model, test="Chisq")
anova(model1, test="Chisq")
anova(model2, test="Chisq")
# Results do not change


####################################

# AUC curve

####################################
library(pROC)

fit1<- glm(group ~ mdNLR,family=binomial(link='logit'),data=data)
summary(fit1)

fit2<- glm(group ~ age+gender+smokingstatus,family=binomial(link='logit'),data=data)
summary(fit2)

fit3 <- glm(group ~ mdNLR+age+gender+smokingstatus,family=binomial(link='logit'),data=data)
summary(fit3)


preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)
auc(roc1)
ci(roc1)
preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)
auc(roc2)
ci(roc2)

preds3=predict(fit3)
roc3=roc(fit3$y , fit3$fitted.values,ci=T)
auc(roc3)
ci(roc3)

roc.test(roc1, roc2)


#roc.plot(roc1,col='black', CI = TRUE)
#ggroc(roc1,alpha = 0.5, colour = "red", linetype = 1, size = 1)

tiff("ROC_modelcomparison_SVA.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggroc(list("Covariates only"=roc2, "mdNLR only"=roc1, "mdNLR+Covariates"=roc3), aes="colour", size=2)
dev.off()
#########################################################

# Selecting myeloid specific CpGs from https://www.ncbi.nlm.nih.gov/pubmed/28184256 

########################################################

mset_myeloid<- mset_fil[  featureNames(mset_fil) %in% CpG_myeloid$CpG_myeloid, ]
dim(mset_myeloid)
betas_myeloid<-betas(mset_myeloid)
save(betas_myeloid, file="betas_myeloid.RData")
write.csv(betas_myeloid, "betas_myeloid.csv")

# Transpose the myeloid CPGs and merge with the data with 689 samples
setwd("~/mDNLR project/Rheumatoid Arthritis/mdNLR_RA_GSE42681")
load("~/mDNLR project/Rheumatoid Arthritis/mdNLR_RA_GSE42681/betas_myeloid.RData")
myeloid<-data.frame(betas_myeloid)
tmyeloid<-as.data.frame(t(myeloid))
row.names(tmyeloid)
colnames(tmyeloid)
tmyeloid <- cbind(sample_name = rownames(tmyeloid), tmyeloid)
# Load pdata with 689 samples
pData<-read.csv("pDataGSE42861.csv", header=T)
data<-read.csv("cellcount_pDataGSE42861_nooutliers_nosmkna_SVA.csv",header=T)
# Merge the 2 datasets
merge<-merge(data, tmyeloid, by.x = "Sample_Name", by.y = "sample_name")

save(merge, file="alldata_gse42861_myeloidCpG.RData")
colnames(merge)
cg00901982_nor<-merge$cg00901982[merge$group=="Normal"]
cg00901982_ra<-merge$cg00901982[merge$group=="Rheumatoid arthritis"]
wilcox.test(cg00901982_nor, cg00901982_ra)# p-value < 2.2e-16
res.by<-by(merge$cg00901982,merge$group, mean)
res.by
res.by<-by(merge$cg00901982,merge$group, range)
res.by



cg25938803_nor<-merge$cg25938803[merge$group=="Normal"]
cg25938803_ra<-merge$cg25938803[merge$group=="Rheumatoid arthritis"]
wilcox.test(cg25938803_nor, cg25938803_ra)# p-value < 2.2e-16
res.by<-by(merge$cg25938803,merge$group, mean)
res.by
res.by<-by(merge$cg25938803,merge$group, range)
res.by



cg01591037_nor<-merge$cg01591037[merge$group=="Normal"]
cg01591037_ra<-merge$cg01591037[merge$group=="Rheumatoid arthritis"]
wilcox.test(cg01591037_nor, cg01591037_ra)# p-value < 2.2e-16
res.by<-by(merge$cg01591037,merge$group, mean)
res.by
res.by<-by(merge$cg01591037,merge$group, range)
res.by



cg03621504_nor<-merge$cg03621504[merge$group=="Normal"]
cg03621504_ra<-merge$cg03621504[merge$group=="Rheumatoid arthritis"]
wilcox.test(cg03621504_nor, cg03621504_ra)# p-value < 2.2e-16
res.by<-by(merge$cg03621504,merge$group, mean)
res.by
res.by<-by(merge$cg03621504,merge$group, range)
res.by


cg10456459_nor<-merge$cg10456459[merge$group=="Normal"]
cg10456459_ra<-merge$cg10456459[merge$group=="Rheumatoid arthritis"]
wilcox.test(cg10456459_nor, cg10456459_ra)# p-value < 2.2e-16
res.by<-by(merge$cg10456459,merge$group, mean)
res.by
res.by<-by(merge$cg10456459,merge$group, range)
res.by

IQR(cg10456459_nor)
IQR(cg10456459_ra)


test<-merge[,c(14,29:33)]
head(test)
library(reshape2)
test.m <- melt(test)
library(ggplot2)
tiff("boxplot_myeloid_CpGs.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggplot(test.m, aes(x=variable, y=value, fill = group)) +
  geom_boxplot(outlier.colour = "black", outlier.size = NA) +
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("mdNLR associted CpGs") +
  ylab("Methylation")+
  ggtitle("")
dev.off()
################################################

# Correlation between mdNLR and 5 myeloid CPGs

###############################################

library("ggpubr")
tiff("corrplot_mdNLR_cg00901982.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggscatter(merge, x = "mdNLR" , y = "cg00901982", add = "none", conf.int = TRUE,xlim = c(0,100), ylim = c(0,1), 
          cor.coef = TRUE, cor.method = "pearson",cor.coef.coord = c(14,15),cor.coef.size = 5,
          xlab = "mdNLR", ylab = "cg00901982")
dev.off()


####################################

# AUC curve with CpGs included
# Adding CpGs do not add more AUC it remains around 80%

####################################
library(pROC)

fit1<- glm(as.factor(group) ~ mdNLR,family=binomial(link='logit'),data=merge)
summary(fit1)

fit2<- glm(as.factor(group) ~ age+as.factor(gender)+as.factor(smokingstatus),family=binomial(link='logit'),data=merge)
summary(fit2)

fit3 <- glm(as.factor(group) ~ mdNLR+age+as.factor(gender)+as.factor(smokingstatus),family=binomial(link='logit'),data=merge)
summary(fit3)

# Include mdNLR associated CpGs (All 5 CpGs)
fit4 <- glm(group ~ age+as.factor(gender)+as.factor(smokingstatus)+cg00901982+cg25938803+cg01591037+cg03621504+cg10456459,family=binomial(link='logit'),data=merge)
summary(fit4)

# Include mdNLR associated CpGs (All 1 CpG)

fit5 <- glm(group ~ age+as.factor(gender)+as.factor(smokingstatus)+cg10456459,family=binomial(link='logit'),data=merge)
summary(fit5)
exp(cbind(Odds=coef(fit5), confint(fit5)))



preds=predict(fit1)
roc1=roc(fit1$y , fit1$fitted.values,ci=T)
auc(roc1)
ci(roc1)

preds2=predict(fit2)
roc2=roc(fit2$y , fit2$fitted.values,ci=T)
auc(roc2)
ci(roc2)

preds3=predict(fit3)
roc3=roc(fit3$y , fit3$fitted.values,ci=T)
auc(roc3)
ci(roc3)

preds4=predict(fit4)
roc4=roc(fit4$y , fit4$fitted.values,ci=T)
auc(roc4)
ci(roc4)

roc.test(roc1, roc2)


#roc.plot(roc1,col='black', CI = TRUE)
#ggroc(roc1,alpha = 0.5, colour = "red", linetype = 1, size = 1)

tiff("ROC_modelcomparison_SVA.tif", res=300, compression = "lzw", height=5, width=10, units="in")
ggroc(list("Covariates only"=roc2, "mdNLR only"=roc1, "mdNLR+Covariates"=roc3), aes="colour", size=2)
dev.off()


# Boxplot of individual lymphocyte components
library(ggplot2)
test<-data[,c(13,4,5,6,7)]
head(test)
library(reshape2)
test.m <- melt(test)
library(ggplot2)
tiff("boxplot_RACaCo_CD8T_CD4T_NK_Bcell.tif", res=300, compression = "lzw", height=5, width=10, units="in")
P= ggplot(test.m, aes(x=variable, y=value, fill = group)) +
  geom_boxplot(outlier.colour = "black", outlier.size = NA) +
  scale_y_continuous(limits = quantile(test.m$value, c(0.1, 0.99)))+
  scale_fill_manual(values = c("coral2", "cornflowerblue"))+
  xlab("Lymphocyte subtypes") +
  ylab("Proportion of lymphocytes")+
  ggtitle("")
P+ labs(x="Lymphocyte subtypes")
dev.off()


