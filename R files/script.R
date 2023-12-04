
# R SCRIPT FOR "Integrating joint latent class mixed models and Bayesian network for uncovering clinical subgroups of COVID-19 patients" ## 

# Loading R libraries
library(bnlearn)
library("Rgraphviz")
library("graph")
library("gRain")
library("DoE.base")
library("tableone")
library("lcmm")

# Sourcing R functions for "Bayesian Networks: robustness and sensitivity issues" 
source("functions.R")

# Load data
load("sample_data.RData")

cov1<-c("Age","Gender","Cancer","Diabetes","Hypertension","Cardiovascular_disease","Chronic_lung_disease")
dd<-unique(data[,c("ID",cov1)])

cov2<-c("Gender","Age", "Cancer" ,"Diabetes","Hypertension" , "Cardiovascular_disease" , "Chronic_lung_disease" ,  
"Pneumonia","Diarrhea","Respiratory_symptoms","Fever","Cough")
dataBN<-unique(data[,c("ID",cov2)])

for(i in cov1[-1]){
data[,i]<-as.numeric(data[,i])-1
}


##########################################################################
##### Joint latent class model for serum creatinine and time to death
##########################################################################

m1Cr <- Jointlcmm(fixed=logcreatinine~time_l,random=~time_l,subject='ID',
        link="5-equi-splines",survival = Surv(Time_30,status01_30)~1+Gender+Age60+Cancer+
        Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, ng=1, data=data)

m2Cr <- Jointlcmm(fixed=logcreatinine~time_l,random=~time_l,subject='ID',mixture=~time_l,
        link="5-equi-splines",survival = Surv(Time_30,status01_30)~1+Gender+Age60+Cancer+
        Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, ng=2, data=data, B=m1Cr)

m3Cr <- Jointlcmm(fixed=logcreatinine~time_l,random=~time_l,subject='ID',mixture=~time_l,
        link="5-equi-splines",survival = Surv(Time_30,status01_30)~1+Gender+Age60+Cancer+
        Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, ng=3, data=data, B=m1Cr)

m4Cr <- Jointlcmm(fixed=logcreatinine~time_l,random=~time_l,subject='ID',mixture=~time_l,
        link="5-equi-splines",survival = Surv(Time_30,status01_30)~1+Gender+Age60+Cancer+
        Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, ng=4, data=data, B=m1Cr)

m5Cr <- Jointlcmm(fixed=logcreatinine~time_l,random=~time_l,subject='ID',mixture=~time_l,
        link="5-equi-splines",survival = Surv(Time_30,status01_30)~1+Gender+Age60+Cancer+
        Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, ng=5, data=data, B=m1Cr)

m6Cr <- Jointlcmm(fixed=logcreatinine~time_l,random=~time_l,subject='ID',mixture=~time_l,
        link="5-equi-splines",survival = Surv(Time_30,status01_30)~1+Gender+Age60+Cancer+
        Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, ng=6, data=data, B=m1Cr)

summarytable(m1Cr,m2Cr,m3Cr,m4Cr,m5Cr,m6Cr,which = c("G", "AIC", "BIC", "entropy", "ICL", "%class"))

####The best model is the model with the lowest ICL
mCr <- m2Cr 

####Description Posterior classification from the joint latent class model 
postprob(mCr)

wCr <- mCr$pprob
dd1<-merge(dd,wCr,by="ID")

print(CreateTableOne(vars=cov1,
factorVars=c("Gender","Cancer","Diabetes","Hypertension","Cardiovascular_disease","Chronic_lung_disease"),
data=dd1,strata="class"),nonnormal=T)

####Class-specific predicted trajectory with 95% confidence interval (shaded area)

datnew<-data.frame(data[rep(1,61), c("Age60",cov1[-1])],time_l=seq(0,30,0.5))

predY <- predictY(mCr,newdata=datnew,var.time="time_l",draws=TRUE)

plot(predY,col=c("darkred","darkgreen"),bty="l",lty=1,legend.loc="topleft",ylab="log(Serum Ceatinine)",xlab="Time(Days)",shades=TRUE)

####Class-specific predicted survival curve for a 60 years woman without comorbiditie

plot(mCr,"survival",col=c("darkred","darkgreen"),bty="l",lty=1,legend.loc="bottomleft",xlab="Time(Days)",main="Class-specific event-free probability (Age=60)")


##########################################################################
##### Joint latent class model for CRP and time to death
##########################################################################

m1Cp <- Jointlcmm(fixed=logcrp~time_l,random=~time_l,subject='ID',
        link="5-equi-splines",survival = Surv(Time_30,status01_30)~1+Gender+Age60+Cancer+
        Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, ng=1, data=data)

m2Cp <- Jointlcmm(fixed=logcrp~time_l,random=~time_l,subject='ID',mixture=~time_l,
        link="5-equi-splines",survival = Surv(Time_30,status01_30)~1+Gender+Age60+Cancer+
        Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, ng=2, data=data, B=m1Cp)

m3Cp <- Jointlcmm(fixed=logcrp~time_l,random=~time_l,subject='ID',mixture=~time_l,
        link="5-equi-splines",survival = Surv(Time_30,status01_30)~1+Gender+Age60+Cancer+
        Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, ng=3, data=data, B=m1Cp)

m4Cp <- Jointlcmm(fixed=logcrp~time_l,random=~time_l,subject='ID',mixture=~time_l,
        link="5-equi-splines",survival = Surv(Time_30,status01_30)~1+Gender+Age60+Cancer+
        Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, ng=4, data=data, B=m1Cp)

m5Cp <- Jointlcmm(fixed=logcrp~time_l,random=~time_l,subject='ID',mixture=~time_l,
        link="5-equi-splines",survival = Surv(Time_30,status01_30)~1+Gender+Age60+Cancer+
        Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, ng=5, data=data, B=m1Cp)

m6Cp <- Jointlcmm(fixed=logcrp~time_l,random=~time_l,subject='ID',mixture=~time_l,
        link="5-equi-splines",survival = Surv(Time_30,status01_30)~1+Gender+Age60+Cancer+
        Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, ng=6, data=data, B=m1Cp)

summarytable(m1Cp,m2Cp,m3Cp,m4Cp,m5Cp,m6Cp,which = c("G", "AIC", "BIC", "entropy", "ICL", "%class"))

####The best model is the model with the lowest ICL
mCp <- m2Cp

####Description Posterior classification from the joint latent class model 
postprob(mCp)

wCp <- mCp$pprob
dd1<-merge(dd,wCp,by="ID")

print(CreateTableOne(vars=cov1,
factorVars=c("Gender","Cancer","Diabetes","Hypertension","Cardiovascular_disease","Chronic_lung_disease"),
data=dd1,strata="class"),nonnormal=T)

####Class-specific predicted trajectory with 95% confidence interval (shaded area)

datnew<-data.frame(data[rep(1,61), c("Age60",cov1[-1])],time_l=seq(0,30,0.5))

predY <- predictY(mCp,newdata=datnew,var.time="time_l",draws=TRUE)

plot(predY,col=c("darkred","blue","darkgreen"),bty="l",lty=1,legend.loc="topleft",ylab="log(C-Reactive Protein)",xlab="Time(Days)",shades=TRUE)

####Class-specific predicted survival curve for a 60 years woman without comorbiditie

plot(mCp,"survival",col=c("darkred","blue","darkgreen"),bty="l",lty=1,legend.loc="bottomleft",xlab="Time(Days)",main="Class-specific event-free probability (Age=60)")


##########################################################################
##### Bayesian Network among baseline patients’ characteristics 
##### and membership in latent classes of JLCMM
##########################################################################

colnames(wCr)[2]<-"class_creatinine_surv"
colnames(wCp)[2]<-"class_CRP_surv"

dim(dataBN)
dim(na.omit(dataBN))
dataBN1<-merge(dataBN,wCr,by="ID")
dataBN1<-merge(dataBN1,wCp,by="ID")

dataBN1<-dataBN1[,c(cov2,"class_creatinine_surv","class_CRP_surv")]
dataBN1$class_creatinine_surv<-factor(dataBN1$class_creatinine_surv,labels = c("Class 1","Class 2"))
dataBN1$class_CRP_surv<-factor(dataBN1$class_CRP_surv,labels = c("Class 1","Class 2"))

dataBN1<-bnlearn::discretize(dataBN1, breaks = 4)

blacklist<-tiers2blacklist(tiers=list("Gender","Age",
c("Cancer","Diabetes","Hypertension","Cardiovascular_disease","Chronic_lung_disease"),
c("Pneumonia","Diarrhea","Respiratory_symptoms","Fever","Cough"),
c("class_creatinine_surv","class_CRP_surv")))
blacklist<-rbind(blacklist,c("Gender","Age"))

bn1<-hc(dataBN1,blacklist=blacklist,score="aic")
graphviz.plot(bn1,shape="ellipse")
lll<-bn.fit(bn1,dataBN1)

png(file="Figure4.png",width=6, height=7, units="in", res=300)
plotbn(dataBN1,bn1,v_class=c("class_creatinine_surv","class_CRP_surv"))
dev.off()

#### Hypothetical scenarios. Conditional probabilities given class creatinine = 3 
png(file="Figure5a.png",width=6, height=7, units="in", res=300)
whatIFplot(dataBN1,bn1,evidenceName=c("class_creatinine_surv"), states = list("Class 3"),cex_p=1.1)
dev.off()

#### Hypothetical scenarios. Conditional probabilities given class CRP = 3 
png(file="Figure5b.png",width=6, height=7, units="in", res=300)
whatIFplot(dataBN1,bn1,evidenceName=c("class_CRP_surv"), states = list("Class 3"),cex_p=1.1)
dev.off()

##########################################################################
##### Multivariate joint latent class models
##########################################################################

m_cr1 <- lcmm(logcreatinine~time_l,random=~time_l,subject="ID",link="5-equi-splines",ng=1,data=data)
m_cr2 <- lcmm(logcreatinine~time_l,random=~time_l,mixture=~time_l,subject="ID",link="5-equi-splines",ng=2,data=data, B=m_cr1)
m_cr3 <- lcmm(logcreatinine~time_l,random=~time_l,mixture=~time_l,subject="ID",link="5-equi-splines",ng=3,data=data, B=m_cr1)
m_cr4 <- lcmm(logcreatinine~time_l,random=~time_l,mixture=~time_l,subject="ID",link="5-equi-splines",ng=4,data=data, B=m_cr1)
m_cr5 <- lcmm(logcreatinine~time_l,random=~time_l,mixture=~time_l,subject="ID",link="5-equi-splines",ng=5,data=data, B=m_cr1)
m_cr6 <- lcmm(logcreatinine~time_l,random=~time_l,mixture=~time_l,subject="ID",link="5-equi-splines",ng=6,data=data, B=m_cr1)

m_cp1 <- lcmm(logcrp~time_l,random=~time_l,subject="ID",link="5-equi-splines",ng=1,data=data)
m_cp2 <- lcmm(logcrp~time_l,random=~time_l,mixture=~time_l,subject="ID",link="5-equi-splines",ng=2,data=data, B=m_cp1)
m_cp3 <- lcmm(logcrp~time_l,random=~time_l,mixture=~time_l,subject="ID",link="5-equi-splines",ng=3,data=data, B=m_cp1)
m_cp4 <- lcmm(logcrp~time_l,random=~time_l,mixture=~time_l,subject="ID",link="5-equi-splines",ng=4,data=data, B=m_cp1)
m_cp5 <- lcmm(logcrp~time_l,random=~time_l,mixture=~time_l,subject="ID",link="5-equi-splines",ng=5,data=data, B=m_cp1)
m_cp6 <- lcmm(logcrp~time_l,random=~time_l,mixture=~time_l,subject="ID",link="5-equi-splines",ng=6,data=data, B=m_cp1)

Jmj1 <- mpjlcmm(longitudinal=list(m_cr1,m_cp1),survival=Surv(Time_30,status01_30)~Age60+Gender+Cancer+Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, subject="ID",ng=1,data=data)

Jmj2 <- mpjlcmm(longitudinal=list(m_cr2,m_cp2),survival=Surv(Time_30,status01_30)~Age60+Gender+Cancer+Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, subject="ID",ng=2,data=data)

Jmj3 <- mpjlcmm(longitudinal=list(m_cr3,m_cp3),survival=Surv(Time_30,status01_30)~Age60+Gender+Cancer+Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, subject="ID",ng=3,data=data)

Jmj4 <- mpjlcmm(longitudinal=list(m_cr4,m_cp4),survival=Surv(Time_30,status01_30)~Age60+Gender+Cancer+Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, subject="ID",ng=4,data=data)

Jmj5 <- mpjlcmm(longitudinal=list(m_cr5,m_cp5),survival=Surv(Time_30,status01_30)~Age60+Gender+Cancer+Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, subject="ID",ng=5,data=data)

Jmj6 <- mpjlcmm(longitudinal=list(m_cr5,m_cp5),survival=Surv(Time_30,status01_30)~Age60+Gender+Cancer+Diabetes+Hypertension+Cardiovascular_disease+Chronic_lung_disease, subject="ID",ng=5,data=data)

summarytable(Jmj1,Jmj2,Jmj3,Jmj4,Jmj5,Jmj6,which = c("G", "loglik", "AIC", "BIC", "entropy","%class","ICL" ))


####The best model is the model with the lowest ICL
Jmj <- Jmj2

####Description Posterior classification from the joint latent class model 
postprob(Jmj)

w <- Jmj$pprob
dd1<-merge(dd,w,by="ID")

print(CreateTableOne(vars=cov1,
factorVars=c("Gender","Cancer","Diabetes","Hypertension","Cardiovascular_disease","Chronic_lung_disease"),
data=dd1,strata="class"),nonnormal=T)

####Class-specific predicted trajectory with 95% confidence interval (shaded area)
mod_biv <- update(Jmj)
mult_post_cr <- mod_biv[[1]]
mult_post_cp <- mod_biv[[2]]

datnew<-data.frame(data[rep(1,61), c("Age60",cov1[-1])],time_l=seq(0,30,0.5))

predYcr <- predictY(mult_post_cr,newdata=datnew,var.time="time_l",draws=TRUE)
plot(predYcr,col=c("darkred","darkgreen"),bty="l",lty=1,legend.loc="topleft",ylab="log(Serum Ceatinine)",xlab="Time(Days)",shades=TRUE)

predYcp <- predictY(mult_post_cp,newdata=datnew,var.time="time_l",draws=TRUE)
plot(predYcp,col=c("darkred","darkgreen"),bty="l",lty=1,legend.loc="topleft",ylab="log(C-Reactive Protein)",xlab="Time(Days)",shades=TRUE)

####Class-specific predicted survival curve for a 60 years woman without comorbiditie

plot(Jmj,"survival",col=c("darkred","darkgreen"),bty="l",lty=1,legend.loc="bottomleft",xlab="Time(Days)",main="Class-specific event-free probability (Age=60)")



###############################################################################
##### Bayesian Network among baseline patient characteristics and membership
#####in latent classes of the multivariate joint latent class mode
###############################################################################

dataBN2<-merge(dataBN,w,by="ID")

dataBN2<-dataBN2[,c(2:14)]
dataBN2$class<-factor(dataBN2$class,labels = c("Class 1","Class 2"))

dataBN2<-bnlearn::discretize(dataBN2, breaks = 4)

blacklist<-tiers2blacklist(tiers=list("Gender","Age",
c("Cancer","Diabetes","Hypertension","Cardiovascular_disease","Chronic_lung_disease"),
c("Pneumonia","Diarrhea","Respiratory_symptoms","Fever","Cough"),
"class"))
blacklist<-rbind(blacklist,c("Gender","Age"))

bn2<-hc(dataBN2,blacklist=blacklist,score="aic")
graphviz.plot(bn2,shape="ellipse")
lll<-bn.fit(bn2,dataBN2)

png(file="Figure7.png",width=6, height=7, units="in", res=300)
plotbn(dataBN2,bn2,v_class=c("class"))
dev.off()

#### Hypothetical scenarios by fixing Age and Cardiovascular disease.
png(file="Figure8a.png",width=6, height=7, units="in", res=300)
whatIFplot(dataBN2,bn2,evidenceName=c("Age","Cardiovascular_disease"), states = list("(72,80]","No"),cex_p=1.1)
dev.off()

png(file="Figure8b.png",width=6, height=7, units="in", res=300)
whatIFplot(dataBN2,bn2,evidenceName=c("Age","Cardiovascular_disease"), states = list("(72,80]","Yes"),cex_p=1.1)
dev.off()




