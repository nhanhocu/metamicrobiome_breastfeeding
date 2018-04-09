#' ---
#' title: "Effects of exclusive breastfeeding on infant gut microbiota: a meta-analysis across studies and populations"
#' author: "Nhan Thi Ho"
#' output:
#'  html_document:
#'   toc: yes
#'   toc_float: true
#' ---

#+ comment="",cache=FALSE,message=FALSE,warning=FALSE,echo=FALSE
rm(list=ls()) # clear all
library(knitr)

opts_chunk$set(comment="",cache=FALSE,message=FALSE,warning=FALSE,echo=FALSE)
# load all required packages
library(caret)
library(ggplot2)
library(digest)
library(plyr)
library(gridExtra)
library(gmodels)
library(chron)
library(lubridate)
library(date)
library(dplyr)
library(data.table)
library(dplyr)
library(dtplyr)
library(tidyr)
library(knitr)
library(lme4)
library(sjmisc)
library(sjPlot)
library(lmerTest)
library(gridExtra)
library(reshape2)
library(mgcv)
library(itsadug)
library(zoo)
library(geepack)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(gdata)
library(randomForest)
# metaanalysis
library(meta)


dir<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/metamicrobiome_breastfeeding/"  #to be replaced by your directory
source(paste(dir,"miscfun.microbiome.R",sep=""))

# load data
load(paste(dir,"data/sam.rm.rda",sep="")) 
load(paste(dir,"data/taxlist.filter.rda",sep="")) 
load(paste(dir,"data/sumstud.rda",sep=""))
load(paste(dir,"data/bangladesh.rda",sep=""))

#' # Summary of included studies

#' All studies
kable(n6.bf.all)

#' All studies for stratified meta-analysis by birth mode
kable(n6.bm.all)
#' # Microbiome age
#' ## Analysis for microbiome age based on shared genera
#' Microbiome age was predicted based on the Random Forest model using L6 (genus) taxa relative abundance with list of taxa shared by all 7 included studies.
#'

load(paste(dir,"data/rmdat.shareg7s.rda",sep="")) 
load(paste(dir,"data/rmshare.conbfg7s.rda",sep="")) 

#' ## Plot of microbiome age by study all ages
#' With Generalized additive mixed model (GAMM) fit and 95%CI.
#'

rmdat.ha$personid<-as.factor(rmdat.ha$sampleid)
rmdat.hav$personid<-as.factor(rmdat.hav$sampleid)
rmdat.ca$personid<-as.factor(rmdat.ca$sampleid)
rmdat.all<-rbind.fill(rmdat.ha,rmdat.rm,rmdat.usbmk,rmdat.uw,rmdat.unc,rmdat.hav,rmdat.ca)
rmdat.all$study<-as.factor(rmdat.all$study)
b2<-gamm(age.predicted~s(age.sample,by=study),family=gaussian,
         data=rmdat.all,random=list(personid=~1))
pred <- predict(b2$gam, newdata = rmdat.all,se.fit=TRUE)
datfit<-cbind(rmdat.all, fit=pred$fit,ul=(pred$fit+(1.96*pred$se.fit)),ll=(pred$fit-(1.96*pred$se.fit)))
rmdat.all.6<-rmdat.all[rmdat.all$age.sample<=6,]

#' ### Sample age < 6 months only
#+ fig.width=10, fig.height=5
rmdat.all.6<-rmdat.all[rmdat.all$age.sample<=6,]
datfit.6<-datfit[datfit$age.sample<=6,]
p<-ggplot()+ geom_point(data = rmdat.all.6, aes(x = age.sample, y = age.predicted, group = personid, colour=bf))+
  geom_line(data = rmdat.all.6, aes(x = age.sample, y = age.predicted, group = personid, colour=bf),size=0.3)+
  geom_line(data = datfit.6,aes(x = age.sample, y = fit),size = 1)+
  #geom_line(data = datfit.6,aes(x = age.sample, y = ul),size = 0.1)+
  #geom_line(data = datfit.6,aes(x = age.sample, y = ll),size = 0.1)+
  geom_ribbon(data = datfit.6,aes(x=age.sample, ymax=ul, ymin=ll), alpha=.5)+
  xlab("Chronological age (month)") +ylab("Microbiome age (month)")+
  labs(color='')+
  theme(legend.position = "bottom",
        #plot.background = element_blank(),
        #panel.background = element_blank()
        axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
        strip.background =element_rect(fill="white"))+
  facet_grid(.~ pop)
p


#' ### Standardized microbiome age and GAMM fit comparison between bf group within each study
rmdat.all.nona<-rmdat.all[!is.na(rmdat.all$bf),]
rmdat.all.nona$pop.bf<-paste(rmdat.all.nona$pop,rmdat.all.nona$bf,sep="_")
rmdat.all.nona$pop.bf<-as.factor(rmdat.all.nona$pop.bf)
b2<-gamm(age.predicteds~s(age.sample,by=pop.bf),family=gaussian,
         data=rmdat.all.nona,random=list(personid=~1))
pred <- predict(b2$gam, newdata = rmdat.all.nona,se.fit=TRUE)
datfit<-cbind(rmdat.all.nona, fit=pred$fit,ul=(pred$fit+(1.96*pred$se.fit)),ll=(pred$fit-(1.96*pred$se.fit)))
p.rms.bf<-ggplot()+ geom_point(data = subset(rmdat.all.nona,age.sample<=6), aes(x = age.sample, y = age.predicteds, group = personid, colour=bf),size=1)+
  geom_line(data = subset(rmdat.all.nona,age.sample<=6), aes(x = age.sample, y = age.predicteds, group = personid, colour=bf),size=0.1)+
  geom_line(data = subset(datfit,age.sample<=6),aes(x = age.sample, y = fit, colour=bf),size = 1)+
  geom_ribbon(data = subset(datfit,age.sample<=6),aes(x=age.sample, ymax=ul, ymin=ll, fill=bf), alpha=.4)+guides(fill=FALSE)+
  xlab("Chronological age (month)") +ylab("Standardized microbiome age")+
  labs(color='')+
  theme(legend.position = "bottom",
        #plot.background = element_blank(),
        #panel.background = element_blank()
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="white"))+
  facet_grid(.~ pop)
p.rms.bf



#' ## Meta-analysis for samples in <= 6 months old infants
#' Meta-analysis models based on adjusted estimate (adjusted for age of infant at sample collection) and standard error from linear mixed effect models.
#'
#' ### Change of RM in non-exclusive breastfed (nebf) vs. exclusive breastfed (exbf)
#+ fig.width=10, fig.height=5
rmshare.sum.6<-as.data.frame(rmshare.sum.6s)
rm.nebf<-metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=rmshare.sum.6,sm="RD", backtransf=FALSE)
forest(rm.nebf,smlab="Standardized \n microbiome age difference",sortvar=rmshare.sum.6$pop)
rm.nebf
kable(cbind(study=rm.nebf$studlab,pval=rm.nebf$pval))


#' ### Sensitivity analysis
#' #### No Haiti data
#+ fig.width=10, fig.height=5
rmshare.sum.6.noha<-rmshare.sum.6[!rownames(rmshare.sum.6) %in% "ha",]
rm.nebf<-metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=rmshare.sum.6.noha,sm="RD", backtransf=FALSE)
forest(rm.nebf,smlab="Standardized \n microbiome age difference",sortvar=rmshare.sum.6.noha$pop)
rm.nebf
kable(cbind(study=rm.nebf$studlab,pval=rm.nebf$pval))
#' #### No USA(NC) data
#+ fig.width=10, fig.height=5
rmshare.sum.6.nounc<-rmshare.sum.6[!rownames(rmshare.sum.6) %in% "unc",]
rm.nebf<-metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=rmshare.sum.6.nounc,sm="RD", backtransf=FALSE)
forest(rm.nebf,smlab="Standardized \n microbiome age difference",sortvar=rmshare.sum.6.nounc$pop)
rm.nebf
kable(cbind(study=rm.nebf$studlab,pval=rm.nebf$pval))
#' #### No USA(CA_MA_MO) data
#+ fig.width=10, fig.height=5
rmshare.sum.6.nohav<-rmshare.sum.6[!rownames(rmshare.sum.6) %in% "hav",]
rm.nebf<-metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=rmshare.sum.6.nohav,sm="RD", backtransf=FALSE)
forest(rm.nebf,smlab="Standardized \n microbiome age difference",sortvar=rmshare.sum.6.nohav$pop)
rm.nebf
kable(cbind(study=rm.nebf$studlab,pval=rm.nebf$pval))

#' ### Stratify by birth mode (not standardized)
load(paste(dir,"data/rmshare.vagcs.rda",sep="")) 
#' #### Vaginal
#+ fig.width=10, fig.height=5
rm.nebf<-metagen(estimate.nebf, se.nebf, studlab=study,data=rmshare.vag,sm="RD", backtransf=FALSE)
forest(rm.nebf,smlab="Microbiome age difference")
rm.nebf
kable(cbind(study=rm.nebf$studlab,pval=rm.nebf$pval))
#' #### C-section
#+ fig.width=10, fig.height=5
rm.nebf<-metagen(estimate.nebf, se.nebf, studlab=study,data=rmshare.cs,sm="RD", backtransf=FALSE)
forest(rm.nebf,smlab="Microbiome age difference")
rm.nebf
kable(cbind(study=rm.nebf$studlab,pval=rm.nebf$pval))


#' ### Trend of microbiome age in exclusive breastfed (exbf), non-exclusive breastfed (nebf) and no bf
#+ fig.width=10, fig.height=5
rmshare.conbf<-rmshare.conbfs
rmshare.conbf$pop<-c("Bangladesh","USA(CA_FL)","USA(CA_MA_MO)","Canada","USA(NC)")
rmshare.conbf<-rmshare.conbf[order(rmshare.conbf$pop),]
rm.conbf<-metagen(estimate.conbf, se.conbf, studlab=study,data=rmshare.conbf,sm="RD", backtransf=FALSE)
forest(rm.conbf,smlab="Standardized \n microbiome age difference",sortvar=rmshare.conbf$pop)
rm.conbf
kable(cbind(study=rm.conbf$studlab,pval=rm.conbf$pval))


#' ## Exploratory analysis for effect of duration of exbf, formula and solid intro
#' With GAMM fit and 95%CI.
#'
#' ### Subramanian (Bangladesh) data
rmdat.rm<-rmdat.rm %>% group_by(personid) %>% arrange(personid,age.sample) %>%
  mutate(month.food6=cut(month.food, breaks=c(-Inf, 6, Inf), labels=c("<=6 months",">6 months")),
         month.food5=cut(month.food, breaks=c(-Inf, 5, Inf), labels=c("<=5 months",">5 months")),
         month.food4=cut(month.food, breaks=c(-Inf, 4, Inf), labels=c("<=4 months",">4 months")),
         month.exbf3=cut(month.exbf, breaks=c(-Inf, 3, Inf), labels=c("<=3 months",">3 months")),
         month.exbf2=cut(month.exbf, breaks=c(-Inf, 2, Inf), labels=c("<=2 months",">2 months")),
         month.exbf1=cut(month.exbf, breaks=c(-Inf, 1, Inf), labels=c("<=1 months",">1 months")))
#' Number of infants by duration of bf in the test set
table(rmdat.rm$month.exbf2[duplicated(rmdat.rm$child.id)==FALSE])
#' Number of samples by duration of bf in the test set
table(rmdat.rm$month.exbf2)
#' Number of samples by duration of bf in the test set (>6 months only)
table(rmdat.rm$month.exbf2[rmdat.rm$age.sample>6])


#' #### Duration (month) of exbf
rmdat.rm<-as.data.frame(rmdat.rm)
b2<-gamm(age.predicted~s(age.sample,by=month.exbf2) +month.exbf2,family=gaussian,
         data=rmdat.rm,random=list(personid=~1))
pred <- predict(b2$gam, newdata = rmdat.rm,se.fit=TRUE)
datfit<-cbind(rmdat.rm, fit=pred$fit,ul=(pred$fit+(1.96*pred$se.fit)),ll=(pred$fit-(1.96*pred$se.fit)))
pexbf2<-ggplot()+ geom_point(data = rmdat.rm, aes(x = age.sample, y = age.predicted, group = personid, colour=month.exbf2))+
  geom_line(data = rmdat.rm, aes(x = age.sample, y = age.predicted, group = personid, colour=month.exbf2),size=0.3)+
  geom_line(data = datfit,aes(x = age.sample, y = fit, colour=month.exbf2),size = 1.5)+
  geom_ribbon(data = datfit,aes(x=age.sample, ymax=ul, ymin=ll,group=month.exbf2), alpha=.5)+
  #theme(legend.position = "bottom")+
  theme(legend.title = element_text(colour="black", size=10))+
  labs(color='Duration EBF')+
  scale_x_continuous(breaks=seq(from=0,to=24,by=3),
                     labels=seq(from=0,to=24,by=3))+
  #theme(legend.text = element_text(colour="black", size = 10))+
  theme(legend.position = c(0.15,0.85),legend.title=element_text(size=8),legend.text=element_text(size=6))+
  theme(legend.key.size = unit(0.5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("Chronological age (month)") +ylab("Microbiome age (month)")
pexbf2
rmdat.rm.6plus15<-rmdat.rm[rmdat.rm$age.sample>6 & rmdat.rm$age.sample<15,]
b2<-gamm(age.predicted~s(age.sample,by=month.exbf2) +month.exbf2,family=gaussian,
         data=rmdat.rm.6plus15,random=list(personid=~1))
#' Test for age > 6 months and <15months
#' GAM part
#'
kable(summary(b2$gam)$p.table)
kable(summary(b2$gam)$s.table)
#' LME part
#'
kable(summary(b2$lme)$tTable)



#' ## Performance of the Random Forest model to estimate microbiome age
load(paste(dir,"data/SrfFit.rml6.shareg7.train.test.rda",sep="")) 
#' ### Evaluated on the training and test set of Subramanian (Bangladesh) data.
Straining$age.predicted <- predict(SrfFit.rml6.share, newdata = Straining)
actual<-Straining$Sage
predict<-Straining$age.predicted
R2 <- 1 - (sum((actual-predict )^2)/sum((actual-mean(actual))^2))
R2<-round(R2,2)
ptrain<-ggplot() +geom_point(data=Straining,aes(x=Sage, y=age.predicted))+
  theme(legend.text = element_text(colour="black", size = 10))+
  annotate("text", x=15, y=5,label=paste("R squared =",R2,sep=" ")) +
  labs(title="Training set")+
  theme(legend.key.size = unit(0.5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("Chronological age (month)") +ylab("Microbiome age (month)")

actual<-rmdat.rm$age.sample
predict<-rmdat.rm$age.predicted
R2 <- 1 - (sum((actual-predict )^2)/sum((actual-mean(actual))^2))
R2<-round(R2,2)
ptest<-ggplot() +geom_point(data=rmdat.rm,aes(x=age.sample, y=age.predicted))+
  theme(legend.text = element_text(colour="black", size = 10))+
  annotate("text", x=15, y=5,label=paste("R squared =",R2,sep=" ")) +
  labs(title="Test set")+
  theme(legend.key.size = unit(0.5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("Chronological age (month)") +ylab("Microbiome age (month)")
grid.arrange(ptrain, ptest,nrow=1)

#' ### List of shared taxa and their relative importance
library(randomForest)
taxim<-as.data.frame(importance(SrfFit.rml6.share$finalModel))
taxim$genera<-rownames(taxim)
taxim$importance<-taxim[,"IncNodePurity"]
taxim<-taxim[order(taxim[,"IncNodePurity"],decreasing = TRUE),]
taxim<-taxim[,c("genera","importance")]
rownames(taxim)<-NULL
kable(taxim)





#' # Alpha diversity indexes
#' ## Plots of alpha diversity by study by age
#' With GAMM fit and 95%CI.
#'
#+ fig.width=7, fig.height=10
#standardize
load(paste(dir,"data/alphamean7s.pooled.rda",sep="")) 
load(paste(dir,"data/alphaall7s.rda",sep="")) 
load(paste(dir,"data/alphameta.allindex7s.rda",sep="")) 

alphap$pop<-as.factor(alphap$pop)
b2<-gamm(shannon~s(age.sample,by=pop),family=gaussian,
         data=alphap,random=list(personid=~1))
pred <- predict(b2$gam, newdata = alphap,se.fit=TRUE)
datfit<-cbind(alphap, fit=pred$fit,ul=(pred$fit+(1.96*pred$se.fit)),ll=(pred$fit-(1.96*pred$se.fit)))
#' Samples <= 6 months only
p.sha6<-ggplot()+ geom_point(data = subset(alphap,age.sample<=6), aes(x = age.sample, y = shannon, group = personid, colour=bf))+
  geom_line(data = subset(alphap,age.sample<=6), aes(x = age.sample, y = shannon, group = personid, colour=bf), size=0.3)+
  geom_line(data = subset(datfit,age.sample<=6),aes(x = age.sample, y = fit),size = 1)+
  geom_ribbon(data = subset(datfit,age.sample<=6),aes(x=age.sample, ymax=ul, ymin=ll), alpha=.5)+
  theme(legend.position = "bottom")+
  theme(legend.title = element_text(colour="black", size=10))+
  labs(color='')+
  scale_x_continuous(breaks=seq(from=0,to=24,by=2),
                     labels=seq(from=0,to=24,by=2))+
  theme(legend.text = element_text(colour="black", size = 10))+
  theme(legend.key.size = unit(0.5, "cm"))+
  theme(legend.position = "bottom",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="white"))+
  xlab("Infant age (months)") +ylab("Standardized Shannon index")+
  facet_grid(. ~ pop)
p.sha6

#' ### GAMM fit comparison between bf group within each study
alphap.nona<-alphap[!is.na(alphap$bf),]
alphap.nona$pop.bf<-paste(alphap.nona$pop,alphap.nona$bf,sep="_")
alphap.nona$pop.bf<-as.factor(alphap.nona$pop.bf)
b2<-gamm(shannon~s(age.sample,by=pop.bf),family=gaussian,
         data=alphap.nona,random=list(personid=~1))
pred <- predict(b2$gam, newdata = alphap.nona,se.fit=TRUE)
datfit<-cbind(alphap.nona, fit=pred$fit,ul=(pred$fit+(1.96*pred$se.fit)),ll=(pred$fit-(1.96*pred$se.fit)))
#' Samples <=6 months only
p.sha.rm.bf6<-ggplot()+ geom_point(data = subset(alphap.nona,age.sample<=6), aes(x = age.sample, y = shannon, group = personid, colour=bf), size=1)+ #, colour=bf
  geom_line(data = subset(alphap.nona,age.sample<=6), aes(x = age.sample, y = shannon, group = personid, colour=bf),size=0.1)+ #, colour=bf
  geom_line(data = subset(datfit,age.sample<=6),aes(x = age.sample, y = fit, colour=bf),size = 1)+
  geom_ribbon(data = subset(datfit,age.sample<=6),aes(x=age.sample, ymax=ul, ymin=ll, fill=bf), alpha=.4)+ guides(fill=FALSE)+
  theme(legend.position = "bottom")+
  theme(legend.title = element_text(colour="black", size=10))+
  labs(color='')+
  scale_x_continuous(breaks=seq(from=0,to=24,by=2),
                     labels=seq(from=0,to=24,by=2))+
  theme(legend.text = element_text(colour="black", size = 10))+
  theme(legend.key.size = unit(0.5, "cm"))+
  theme(legend.position = "bottom",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="white"))+
  xlab("Infant age (months)") +ylab("Standardized Shannon index")+
  facet_grid(. ~ pop)
p.sha.rm.bf6


#' ## Exploration by duration of exclusive bf 
#' Subramanian (Bangladesh) data only.
#'
samfile<-merge(samde, he50[,c("child.id","gender","month.exbf","month.food")],by="child.id")
samfile$age.sample<-samfile$age.months
samfile$bf<-factor(samfile$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
samfile$personid<-as.factor(paste("rm", tolower(samfile$child.id),sep="."))
samfile$sampleid<- paste("rm",tolower(samfile$fecal.sample.id),sep=".")
samfile$author<-"Subramanian et al"
samfile$year<-"2014"
samfile$pop<-"Bangladesh"


#' ### Exclusive bf duration
#' Shannon index.
#'
alpha.m.rm<-alpha.m.rm %>% group_by(personid) %>% arrange(personid,age.sample) %>%
  mutate(month.food6=cut(month.food, breaks=c(-Inf, 6, Inf), labels=c("<=6 months",">6 months")),
         month.food5=cut(month.food, breaks=c(-Inf, 5, Inf), labels=c("<=5 months",">5 months")),
         month.food4=cut(month.food, breaks=c(-Inf, 4, Inf), labels=c("<=4 months",">4 months")),
         month.foodr=as.factor(sort(round(month.food,0))),
         month.exbf3=cut(month.exbf, breaks=c(-Inf, 3, Inf), labels=c("<=3 months",">3 months")),
         month.exbf2=cut(month.exbf, breaks=c(-Inf, 2, Inf), labels=c("<=2 months",">2 months")),
         month.exbf1=cut(month.exbf, breaks=c(-Inf, 1, Inf), labels=c("<=1 months",">1 months")),
         month.exbfr=as.factor(sort(round(month.exbf,0))))
ggplot(data=alpha.m.rm,aes(x=age.sample, y=shannon, colour=month.exbf2, group=personid)) +geom_point()+geom_line() +
  geom_smooth(data=alpha.m.rm,aes(x=age.sample, y=shannon,group=month.exbf2))


#' ## Meta-analysis
#' For samples <=6 months old only.
#'
#' ### Change in non-exbf vs. exbf
#' #### Chao1
#+ fig.width=12, fig.height=5
chao1.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="chao1"),sm="RD", backtransf=FALSE)
forest(chao1.nebf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="chao1")$pop)
chao1.nebf
kable(cbind(study=chao1.nebf$studlab,pval=chao1.nebf$pval))
#' #### Observed_species
#+ fig.width=12, fig.height=5
os.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="observed_species"),sm="RD", backtransf=FALSE)
forest(os.nebf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="observed_species")$pop)
os.nebf
kable(cbind(study=os.nebf$studlab,pval=os.nebf$pval))
#' #### Pd_whole_tree
#+ fig.width=10, fig.height=5
pwt.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="pd_whole_tree"),sm="RD", backtransf=FALSE)
forest(pwt.nebf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="pd_whole_tree")$pop)
pwt.nebf
kable(cbind(study=pwt.nebf$studlab,pval=pwt.nebf$pval))
#' #### Shannon
#+ fig.width=10, fig.height=5
shannon.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="shannon"),sm="RD", backtransf=FALSE)
forest(shannon.nebf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="shannon")$pop)
shannon.nebf
kable(cbind(study=shannon.nebf$studlab,pval=shannon.nebf$pval))

#' ### Sensitivity analysis Change in non-exbf vs. exbf
#' Show the results of Shannon indexes only.
#'
#' #### No Haiti data
#+ fig.width=10, fig.height=5
shannon.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="shannon"&(pop!="Haiti")),sm="RD", backtransf=FALSE)
forest(shannon.nebf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="shannon"&(pop!="Haiti"))$pop)
shannon.nebf
kable(cbind(study=shannon.nebf$studlab,pval=shannon.nebf$pval))
#' #### No UNC data
#+ fig.width=10, fig.height=5
shannon.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="shannon"&(pop!="USA(NC)")),sm="RD", backtransf=FALSE)
forest(shannon.nebf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="shannon"&(pop!="USA(NC)"))$pop)
shannon.nebf
kable(cbind(study=shannon.nebf$studlab,pval=shannon.nebf$pval))
#' #### No USA(CA_MA_MO) data
#+ fig.width=10, fig.height=5
shannon.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="shannon"&(pop!="USA(CA_MA_MO)")),sm="RD", backtransf=FALSE)
forest(shannon.nebf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="shannon"&(pop!="USA(CA_MA_MO)"))$pop)
shannon.nebf
kable(cbind(study=shannon.nebf$studlab,pval=shannon.nebf$pval))

#' ### Stratify by birth mode for Change in non-exbf vs. exbf (not standardized)
#' Show results of Shannon index only.
#'
load(paste(dir,"data/alphasum.vagcs.rda",sep="")) 
#' #### Vaginal
#+ fig.width=10, fig.height=5
shannon.nebf <- metagen(estimate.nebf, se.nebf, studlab=study,data=subset(alphasum.vag, index=="shannon"),sm="RD", backtransf=FALSE)
forest(shannon.nebf,smlab="Standardized \n diversity difference")
shannon.nebf
kable(cbind(study=shannon.nebf$studlab,pval=shannon.nebf$pval))
#' #### C-section
#+ fig.width=10, fig.height=5
shannon.nebf <- metagen(estimate.nebf, se.nebf, studlab=study,data=subset(alphasum.cs, index=="shannon"),sm="RD", backtransf=FALSE)
forest(shannon.nebf,smlab="Standardized \n diversity difference")
shannon.nebf
kable(cbind(study=shannon.nebf$studlab,pval=shannon.nebf$pval))


#' #### Put meta-analysis results (random models) of all indexes together
a.nebf<-ggplot(data=a.metatab.r,aes(x=estimate.nebf,y=index))+
  geom_point(shape=17, colour="red")+
  geom_errorbarh(aes(xmin=ll.nebf,xmax=ul.nebf),height=0.0, colour="red")+
  scale_x_continuous(breaks=seq(from=0,to=1,by=0.2),
                     labels=seq(from=0,to=1,by=0.2))+
  geom_vline(xintercept=0,linetype="dashed")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("Non-EBF vs. EBF")+
  xlab("Pooled Standardized diversity difference")+ylab("Alpha diversity index")
a.nebf
kable(a.metatab.r[,grep(".nebf",colnames(a.metatab.r))])

#' ### Trend effect in exbf, non-exbf, no bf
#' #### Chao1
#+ fig.width=10, fig.height=5
chao1.conbf <- metagen(estimate.conbf, se.conbf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="chao1"),sm="RD", backtransf=FALSE)
forest(chao1.conbf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="chao1")$pop)
chao1.conbf
kable(cbind(study=chao1.conbf$studlab,pval=chao1.conbf$pval))
#' #### Observed_species
#+ fig.width=10, fig.height=5
os.conbf <- metagen(estimate.conbf, se.conbf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="observed_species"),sm="RD", backtransf=FALSE)
forest(os.conbf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="observed_species")$pop)
os.conbf
kable(cbind(study=os.conbf$studlab,pval=os.conbf$pval))
#' #### Pd_whole_tree
#+ fig.width=10, fig.height=5
pwt.conbf <- metagen(estimate.conbf, se.conbf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="pd_whole_tree"),sm="RD", backtransf=FALSE)
forest(pwt.conbf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="pd_whole_tree")$pop)
pwt.conbf
kable(cbind(study=pwt.conbf$studlab,pval=pwt.conbf$pval))
#' #### Shannon
#+ fig.width=10, fig.height=5
shannon.conbf <- metagen(estimate.conbf, se.conbf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="shannon"&(!pop %in% c("Haiti","South Africa"))),sm="RD", backtransf=FALSE)
forest(shannon.conbf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="shannon"&(!pop %in% c("Haiti","South Africa")))$pop)
shannon.conbf
kable(cbind(study=shannon.conbf$studlab,pval=shannon.conbf$pval))

#' #### Put meta-analysis results (random models) of all indexes together
a.conbf<-ggplot(data=a.metatab.r,aes(x=estimate.conbf,y=index))+
  geom_point(shape=17, colour="red")+
  geom_errorbarh(aes(xmin=ll.conbf,xmax=ul.conbf),height=0.0, colour="red")+
  geom_vline(xintercept=0,linetype="dashed")+
  scale_x_continuous(breaks=seq(from=0,to=1,by=0.2),
                     labels=seq(from=0,to=1,by=0.2))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("Trend across EBF, non-EBF, non-BF")+
  xlab("Pooled standardized diversity difference")+ylab("Alpha diversity index")
a.conbf
kable(a.metatab.r[,grep(".conbf",colnames(a.metatab.r))])


#' # Meta-analysis of taxa relative abundance
#' Results of 7 studies.
#'
#' For samples <= 6 months old only in all studies (note for USA(NC) study: GAMLSS BEZI with random subject effect could not run on very small sample size=> did not include subject random effect).
#' Results of random meta-analysis models for taxa available in at least >50% of studies
#' based on adjusted estimates and standard errors from GAMLSS models with zero-inflated beta family
#' adjusted for infant age at sample collection.
#'

#' ## Meta-analysis of Change in non-exbf vs. exbf adjusted for age
load(paste(dir,"data/metatab.zi7.rda",sep="")) 
#' GAMLSS models with zero-inflated beta family
#'
#' ### Phylum (l2)
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="heatmap",fill.value="log(OR)")
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="forest",fill.value="log(OR)")
kable(metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' Nice plot
#+ fig.width=10,fig.height=5
metadat<-metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="main",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=8,forest.axis.text.x=8)

#' ### Order (l4)
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="heatmap",fill.value="log(OR)")
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="forest",fill.value="log(OR)")
kable(metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' Nice plot
#+ fig.width=10, fig.height=5
metadat<-metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=8,forest.axis.text.x=8)

#' ### Family (l5)
#+ fig.width=8, fig.height=10
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="heatmap",fill.value="log(OR)")
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="forest",fill.value="log(OR)")
kable(metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' Nice plot
#+ fig.width=10, fig.height=10
metadat<-metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=8,forest.axis.text.x=8)

#' ### Genus (l6)
#+ fig.width=8, fig.height=12
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="heatmap",fill.value="log(OR)")
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="forest",fill.value="log(OR)")
kable(metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' Nice plot
#+ fig.width=10, fig.height=12
metadat<-metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=8,forest.axis.text.x=8)


#' ## Sensitivity meta-analysis of Change in non-exbf vs. exbf adjusted for age
#' ### No UNC data
load(paste(dir,"data/metatab.zi.nounc7.rda",sep="")) 
#' #### l2
kable(metatab.show(metatab=metatab.zi.nounc$random,taxacom.pooled.tab=taxacom.zi.nounc,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' #### l4
kable(metatab.show(metatab=metatab.zi.nounc$random,taxacom.pooled.tab=taxacom.zi.nounc,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' #### l5
kable(metatab.show(metatab=metatab.zi.nounc$random,taxacom.pooled.tab=taxacom.zi.nounc,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' #### l6
kable(metatab.show(metatab=metatab.zi.nounc$random,taxacom.pooled.tab=taxacom.zi.nounc,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' ### No Haiti data
load(paste(dir,"data/metatab.zi.noha7.rda",sep="")) 
#' #### l2
kable(metatab.show(metatab=metatab.zi.noha$random,taxacom.pooled.tab=taxacom.zi.noha,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' #### l4
kable(metatab.show(metatab=metatab.zi.noha$random,taxacom.pooled.tab=taxacom.zi.noha,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' #### l5
kable(metatab.show(metatab=metatab.zi.noha$random,taxacom.pooled.tab=taxacom.zi.noha,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' #### l6
kable(metatab.show(metatab=metatab.zi.noha$random,taxacom.pooled.tab=taxacom.zi.noha,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' ### No USA(CA_MA_MO) data
load(paste(dir,"data/metatab.zi.nohav7.rda",sep="")) 
#' #### l2
kable(metatab.show(metatab=metatab.zi.nohav$random,taxacom.pooled.tab=taxacom.zi.nohav,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' #### l4
kable(metatab.show(metatab=metatab.zi.nohav$random,taxacom.pooled.tab=taxacom.zi.nohav,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' #### l5
kable(metatab.show(metatab=metatab.zi.nohav$random,taxacom.pooled.tab=taxacom.zi.nohav,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' #### l6
kable(metatab.show(metatab=metatab.zi.nohav$random,taxacom.pooled.tab=taxacom.zi.nohav,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))


#' ## Stratified analysis by birth mode for change in non-exbf vs. exbf adjusted for age
load(paste(dir,"data/metatab.zi.vagcs.rda",sep="")) 
#' ### Vaginal
#' #### l2
#+ fig.width=10, fig.height=5
kable(metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="main",phyla.col="rainbow",leg.key.size=0.5,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=8,forest.axis.text.x=8)

#' #### l4
#+ fig.width=10, fig.height=5
kable(metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",phyla.col="rainbow",leg.key.size=1,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=8,forest.axis.text.x=8)

#' #### l5
#+ fig.width=10, fig.height=10
kable(metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",phyla.col="rainbow",leg.key.size=1,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=8,forest.axis.text.x=8)

#' #### l6
#+ fig.width=10, fig.height=12
kable(metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",phyla.col="rainbow",leg.key.size=1,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=8,forest.axis.text.x=8)


#' ### C-section
#' #### l2
#+ fig.width=10, fig.height=5
kable(metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="main",phyla.col="rainbow",leg.key.size=0.4,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=8,forest.axis.text.x=8)

#' #### l4
#+ fig.width=10, fig.height=5
kable(metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",phyla.col="rainbow",leg.key.size=1,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=8,forest.axis.text.x=8)

#' #### l5
#+ fig.width=10, fig.height=10
kable(metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",phyla.col="rainbow",leg.key.size=1,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=8,forest.axis.text.x=8)

#' #### l6
#+ fig.width=10, fig.height=12
kable(metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",phyla.col="rainbow",leg.key.size=1,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=8,forest.axis.text.x=8)



#' ## Meta-analysis of Trend in exbf, non-exbf and no bf adjusted for age
#' ### Phylum (l2)
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l2",showvar=".conbf",p.cutoff.type="p", p.cutoff=1,plot="heatmap",fill.value="log(OR)")
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l2",showvar=".conbf",p.cutoff.type="p", p.cutoff=1,plot="forest",fill.value="log(OR)")
#' Significant (pooled p<0.05) only
kable(metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l2",showvar=".conbf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' ### Order (l4)
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l4",showvar=".conbf",p.cutoff.type="p", p.cutoff=1,plot="heatmap",fill.value="log(OR)")
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l4",showvar=".conbf",p.cutoff.type="p", p.cutoff=1,plot="forest",fill.value="log(OR)")
#' Significant only
kable(metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l4",showvar=".conbf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' ### Family (l5)
#+ fig.width=10, fig.height=10
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l5",showvar=".conbf",p.cutoff.type="p", p.cutoff=1,plot="heatmap",fill.value="log(OR)")
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l5",showvar=".conbf",p.cutoff.type="p", p.cutoff=1,plot="forest",fill.value="log(OR)")
#' Significant only
kable(metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l5",showvar=".conbf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' ### Genus (l6)
#+ fig.width=10, fig.height=10
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l6",showvar=".conbf",p.cutoff.type="p", p.cutoff=1,plot="heatmap",fill.value="log(OR)")
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l6",showvar=".conbf",p.cutoff.type="p", p.cutoff=1,plot="forest",fill.value="log(OR)")
#' Significant only
kable(metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l6",showvar=".conbf",p.cutoff.type="p", p.cutoff=0.05,display="table"))




#' # Meta-analysis of KEGG pathway relative abundance
#' ## Change in non-exbf vs. exbf
load(paste(dir,"data/pathmetatab7.rda",sep="")) 
#' ### Level 2 KEGG pathway
#+ fig.width=8, fig.height=10
metatab.show(metatab=pathmetatab.zi$random,taxacom.pooled.tab=pathcom.zi,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="heatmap",fill.value="log(OR)")
metatab.show(metatab=pathmetatab.zi$random,taxacom.pooled.tab=pathcom.zi,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="forest",fill.value="log(OR)")
kable(metatab.show(metatab=pathmetatab.zi$random,taxacom.pooled.tab=pathcom.zi,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' Nice plot
#+ fig.width=10, fig.height=8
metadat<-metatab.show(metatab=pathmetatab.zi$random,taxacom.pooled.tab=pathcom.zi,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=8,forest.axis.text.x=8,heat.text.x.angle=0)

#' ### Level 3 KEGG pathway
kable(metatab.show(metatab=pathmetatab.zi$random,taxacom.pooled.tab=pathcom.zi,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' Nice plot all pathways
#+ fig.width=10, fig.height=15
metadat<-metatab.show(metatab=pathmetatab.zi$random,taxacom.pooled.tab=pathcom.zi,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=4,forest.axis.text.x=4)

#' Nice plot significant pathways only (pooled p<0.05)
#+ fig.width=11, fig.height=5
metadat<-metatab.show(metatab=pathmetatab.zi$random,taxacom.pooled.tab=pathcom.zi,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=5,forest.axis.text.x=6)

#' Nice plot multiple testing adjusted significant pathways only (adjusted pooled p<0.1)
#+ fig.width=11, fig.height=5
metadat<-metatab.show(metatab=pathmetatab.zi$random,taxacom.pooled.tab=pathcom.zi,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p.adjust", p.cutoff=0.1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=6.5,forest.axis.text.x=6)



#' ### Sensitivity analysis
load(paste(dir,"data/pathmetatab.zi.sen.rda",sep="")) 
#' #### No Haiti data
#' Level 2
#+ fig.width=10, fig.height=7
kable(metatab.show(metatab=pathmetatab.zi.noha$random,taxacom.pooled.tab=pathcom.zi.noha,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=pathmetatab.zi.noha$random,taxacom.pooled.tab=pathcom.zi.noha,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=7,forest.axis.text.x=7)
#' Level 3
#+ fig.width=10, fig.height=5
kable(metatab.show(metatab=pathmetatab.zi.noha$random,taxacom.pooled.tab=pathcom.zi.noha,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=pathmetatab.zi.noha$random,taxacom.pooled.tab=pathcom.zi.noha,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=6,forest.axis.text.x=6)
#' #### No USA(CA_MA_MO) data
#' Level 2
#+ fig.width=10, fig.height=7
kable(metatab.show(metatab=pathmetatab.zi.nohav$random,taxacom.pooled.tab=pathcom.zi.nohav,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=pathmetatab.zi.nohav$random,taxacom.pooled.tab=pathcom.zi.nohav,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=7,forest.axis.text.x=7)
#' Level 3
#+ fig.width=10, fig.height=5
kable(metatab.show(metatab=pathmetatab.zi.nohav$random,taxacom.pooled.tab=pathcom.zi.nohav,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=pathmetatab.zi.nohav$random,taxacom.pooled.tab=pathcom.zi.nohav,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=6,forest.axis.text.x=6)
#' ### No UNC data
#' Level 2
#+ fig.width=10, fig.height=7
kable(metatab.show(metatab=pathmetatab.zi.nounc$random,taxacom.pooled.tab=pathcom.zi.nounc,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=pathmetatab.zi.nounc$random,taxacom.pooled.tab=pathcom.zi.nounc,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=7,forest.axis.text.x=7)
#' Level 3
#+ fig.width=10, fig.height=5
kable(metatab.show(metatab=pathmetatab.zi.nounc$random,taxacom.pooled.tab=pathcom.zi.nounc,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=pathmetatab.zi.nounc$random,taxacom.pooled.tab=pathcom.zi.nounc,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=6,forest.axis.text.x=6)

#' ## stratify by birth mode
load(paste(dir,"data/pathmetatab.zi.vagcs.rda",sep="")) 
#' ### Vaginal birth
#' Level 2
#+ fig.width=10, fig.height=7
kable(metatab.show(metatab=pathmetatab.zi.vag$random,taxacom.pooled.tab=pathcom.zi.vag,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=pathmetatab.zi.vag$random,taxacom.pooled.tab=pathcom.zi.vag,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=7,forest.axis.text.x=7)
#' Level 3
#' All
#+ fig.width=10, fig.height=15
metadat<-metatab.show(metatab=pathmetatab.zi.vag$random,taxacom.pooled.tab=pathcom.zi.vag,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=4,forest.axis.text.x=4)
#' Significant only (pooled p<0.05)
#+ fig.width=10, fig.height=6
kable(metatab.show(metatab=pathmetatab.zi.vag$random,taxacom.pooled.tab=pathcom.zi.vag,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=pathmetatab.zi.vag$random,taxacom.pooled.tab=pathcom.zi.vag,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=5,forest.axis.text.x=6)

#' Multiple testing adjusted Significant only (adjusted pooled p<0.1)
#+ fig.width=10, fig.height=6
metadat<-metatab.show(metatab=pathmetatab.zi.vag$random,taxacom.pooled.tab=pathcom.zi.vag,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p.adjust", p.cutoff=0.1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=6.5,forest.axis.text.x=6)

#' ### C-section birth
#' Level 2
#+ fig.width=10, fig.height=7
kable(metatab.show(metatab=pathmetatab.zi.cs$random,taxacom.pooled.tab=pathcom.zi.cs,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=pathmetatab.zi.cs$random,taxacom.pooled.tab=pathcom.zi.cs,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=7,forest.axis.text.x=7)
#' Level 3
#' All
#+ fig.width=10, fig.height=15
metadat<-metatab.show(metatab=pathmetatab.zi.cs$random,taxacom.pooled.tab=pathcom.zi.cs,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=4,forest.axis.text.x=4)
#' Significant only (pooled p<0.05)
#+ fig.width=10, fig.height=5
kable(metatab.show(metatab=pathmetatab.zi.cs$random,taxacom.pooled.tab=pathcom.zi.cs,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=pathmetatab.zi.cs$random,taxacom.pooled.tab=pathcom.zi.cs,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=5,forest.axis.text.x=6)

#' Multiple testing adjusted Significant only (adjusted pooled p<0.1)
#+ fig.width=10, fig.height=5
metadat<-metatab.show(metatab=pathmetatab.zi.cs$random,taxacom.pooled.tab=pathcom.zi.cs,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p.adjust", p.cutoff=0.1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=6.5,forest.axis.text.x=5)



#' # Taxa relative abundance Bangladesh data only
#' All analyses are adjusted for age of infants or breastfeeding status at sample collection and accounting for repeated/longitudinal sample collection.
#'
#' ## Mean taxa relative abundance by duration of exbf in samples > 6 months
load(paste(dir,"data/taxacom.612plus.food5.exbf2f.rda",sep="")) 
#' For Subramanian (Bangladesh) data only.
#'
#' ### Duration of exbf
#' Change in taxa relative abundance in duration of exclusive bf >2months vs. <=2 months
#'
#' #### Phylum (L2)
p.exbf.l2<-taxa.mean.plot(tabmean=taxa.meansdn.exbf2.rm,taxlist=taxlist.rm,tax.lev="l2", comvar="month.exbf2", groupvar="age.sample",mean.filter=0.005)
p.exbf.l2$p
#dev.off()
#' GAMLSS
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.exbf2.zi.rm,tax.select=p.exbf.l2$taxuse.rm, tax.lev="l2",p.adjust.method="fdr"))

#' #### Order (l4)
p.exbf.l4<-taxa.mean.plot(tabmean=taxa.meansdn.exbf2.rm,taxlist=taxlist.rm,tax.lev="l4", comvar="month.exbf2", groupvar="age.sample",mean.filter=0.005)
p.exbf.l4$p
#' GAMLSS
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.exbf2.zi.rm,tax.select=p.exbf.l4$taxuse.rm, tax.lev="l4",p.adjust.method="fdr"))

#' #### Family (L5)
#+ fig.width=10, fig.height=7
p.exbf.l5<-taxa.mean.plot(tabmean=taxa.meansdn.exbf2.rm,taxlist=taxlist.rm,tax.lev="l5", comvar="month.exbf2", groupvar="age.sample",mean.filter=0.005,show.taxname="short",legend.position="right")
p.exbf.l5$p
# better plot view
for (i in 1: length(taxa.meansdn.exbf2.rm)){
  taxa.meansdn.exbf2.rm[[i]]$month.exbf2l<-mapvalues(taxa.meansdn.exbf2.rm[[i]]$month.exbf2,from=c("<=2 months",">2 months"),to=c("Duration EBF <=2 months","Duration EBF >2 months"))
}
p.exbf.l5.noleg<-taxa.mean.plot(tabmean=taxa.meansdn.exbf2.rm,taxlist=taxlist.rm,tax.lev="l5", comvar="month.exbf2l", groupvar="age.sample",mean.filter=0.005,legend.position="none")
grid.arrange(pexbf2,p.exbf.l5.noleg$p,nrow=1)
#' GAMLSS
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.exbf2.zi.rm,tax.select=p.exbf.l5$taxuse.rm, tax.lev="l5",p.adjust.method="fdr"))

#' #### Genus (L6)
#+ fig.width=15, fig.height=10
p.exbf.l6<-taxa.mean.plot(tabmean=taxa.meansdn.exbf2.rm,taxlist=taxlist.rm,tax.lev="l6", comvar="month.exbf2", groupvar="age.sample",mean.filter=0.005)
p.exbf.l6$p
#' GAMLSS
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.exbf2.zi.rm,tax.select=p.exbf.l6$taxuse.rm, tax.lev="l6",p.adjust.method="fdr"))


#' # Modification effect of breastfeeding status on samples of patients with vs. without diarrhea at sample collection
#' Subramania (Bangladesh) data only.
#'
#' ## Taxa relative abundance in samples of patients with vs. without diarrhea at sample collection
#' Only LME was used as GAMLSS has issues with small sample size (when stratifying).
#' LME as showed above has lower power than GAMLSS in general.
#'
#'
load(paste(dir,"data/taxacom.dia.abf.rda",sep="")) 
load(paste(dir,"data/taxacom.dia.exbf2.zinolong.rda",sep="")) 
load(paste(dir,"data/taxacom.dia.bf.zinolong.rda",sep="")) 
#' ### After 6 month
#' #### Stratified by duration of exclusive bf
#' ##### Family
#+ fig.width=15, fig.height=5
p.dia.exbf2.6plus.l5<-taxa.mean.plot(tabmean=taxa.meansdn.dia.exbf2.6plus.rm,taxlist=taxlist.rm,tax.lev="l5", comvar="diarrhea", groupvar="month.exbf2",mean.filter=0.005,legend.position="right")
p.dia.exbf2.6plus.l5$p
#more detail legend
teste<-taxa.meansdn.dia.exbf2.6plus.rm
for (i in 1:length(names(teste))){
  teste[[i]]$month.exbf2l<-mapvalues(teste[[i]]$month.exbf2,from=c("<=2 months",">2 months"),to=c("Duration EBF <=2 months","Duration EBF >2 months"))
}
p.dia.exbf2.6plus.l5noleg<-taxa.mean.plot(tabmean=teste,taxlist=taxlist.rm,tax.lev="l5", comvar="diarrhea", groupvar="month.exbf2l",mean.filter=0.005,legend.position="none",ylab="Relative abundance (6 months - 2 years)")


#' ###### GAMLSS
#' In infants with duration of exclusive bf <=2 months
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.exbf2.zi.rm, tax.lev="l5",tax.select=p.dia.exbf2.6plus.l5$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1))
#' In infants with duration of exclusive bf >2 months
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.exbf2plus.zi.rm, tax.lev="l5",tax.select=p.dia.exbf2.6plus.l5$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1))



#' ###### Genus
#+ fig.width=15, fig.height=5
p.dia.exbf2.6plus.l6<-taxa.mean.plot(tabmean=taxa.meansdn.dia.exbf2.6plus.rm,taxlist=taxlist.rm,tax.lev="l6", comvar="diarrhea", groupvar="month.exbf2",mean.filter=0.005,legend.position = "right")
p.dia.exbf2.6plus.l6$p
#' ####### GAMLSS
#' In infants with duration of exclusive bf <=2 months
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.exbf2.zi.rm, tax.lev="l6",tax.select=p.dia.exbf2.6plus.l6$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1))
#' In infants with duration of exclusive bf >2 months
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.exbf2plus.zi.rm, tax.lev="l6",tax.select=p.dia.exbf2.6plus.l6$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1))


#' #### Stratified by bf status
#'
#+ fig.width=15, fig.height=5
taxa.meansdn.dia.bf.6plus.rm.s<-list()
for (i in 1:5){
  taxa.meansdn.dia.bf.6plus.rm.s[[i]]<-taxa.meansdn.dia.bf.6plus.rm[[i]][taxa.meansdn.dia.bf.6plus.rm[[i]]$bf!="ExclusiveBF",]
  taxa.meansdn.dia.bf.6plus.rm.s[[i]]$bf<-drop.levels(taxa.meansdn.dia.bf.6plus.rm.s[[i]]$bf,reorder=FALSE)
}
names(taxa.meansdn.dia.bf.6plus.rm.s)<-names(taxa.meansdn.dia.bf.6plus.rm)

#' ##### Family
p.dia.bf.6plus.l5<-taxa.mean.plot(tabmean=taxa.meansdn.dia.bf.6plus.rm.s,tax.select=p.dia.exbf2.6plus.l5$taxuse.rm,taxlist=taxlist.rm,tax.lev="l5", comvar="diarrhea", groupvar="bf",mean.filter=0.005,show.taxname = "short")
p.dia.bf.6plus.l5$p

# reverse levels for better combined plot view
testr=taxa.meansdn.dia.bf.6plus.rm.s
for (i in 1: length(names(testr))){
  testr[[i]]$bfl<-mapvalues(testr[[i]]$bf,from=c("No_BF","Non_exclusiveBF"),to=c("No BF when diarrhea","BF when diarrhea"))
  testr[[i]]$bfl<-factor(testr[[i]]$bfl,levels=c("No BF when diarrhea","BF when diarrhea"))
}
p.dia.bf.6plus.l5.nolegend<-taxa.mean.plot(tabmean=testr,tax.select=p.dia.exbf2.6plus.l5$taxuse.rm,taxlist=taxlist.rm,tax.lev="l5", comvar="diarrhea", groupvar="bfl",mean.filter=0.005,show.taxname = "short",legend.position = "none",ylab="Relative abundance (6 months - 2 years)")

#' ###### GAMLSS
#' In non-exclusive bf infants
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.nexbf.zi.rm, tax.lev="l5",tax.select=p.dia.bf.6plus.l5$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1))
#' In no bf infants
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.nobf.zi.rm, tax.lev="l5",tax.select=p.dia.bf.6plus.l5$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1))

#' ##### Genus
p.dia.bf.6plus.l6<-taxa.mean.plot(tabmean=taxa.meansdn.dia.bf.6plus.rm.s,tax.select=p.dia.exbf2.6plus.l6$taxuse.rm,taxlist=taxlist.rm,tax.lev="l6", comvar="diarrhea", groupvar="bf",mean.filter=0.005)
p.dia.bf.6plus.l6$p
#' ###### GAMLSS
#' In bf infants
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.nexbf.zi.rm, tax.lev="l6",tax.select=p.dia.bf.6plus.l6$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1))
#' In non-bf infants
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.nobf.zi.rm, tax.lev="l6",tax.select=p.dia.bf.6plus.l6$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1))



#' ## Microbiome age in samples of infants with vs. without diarrhea at sample collection stratified by duration of exclusive breastfeeding
#' With GAMM fit and 95%CI.
#'
rmdat.rm$dia.exbf2<-paste(rmdat.rm$diarrhea,rmdat.rm$month.exbf2,sep=".")
rmdat.rm<-as.data.frame(rmdat.rm)
rmdat.rm$dia.exbf2<-as.factor(rmdat.rm$dia.exbf2)
rmdat.rm$diarrhea<-as.factor(rmdat.rm$diarrhea)
rmdat.rm$month.exbf2l<-mapvalues(rmdat.rm$month.exbf2,from=c("<=2 months",">2 months"),to=c("Duration EBF <=2 months","Duration EBF >2 months"))
b2<-gamm(age.predicted~s(age.sample,by=dia.exbf2) +dia.exbf2,family=gaussian,
         data=rmdat.rm,random=list(personid=~1))
pred <- predict(b2$gam, newdata = rmdat.rm,se.fit=TRUE)
datfit<-cbind(rmdat.rm, fit=pred$fit,ul=(pred$fit+(1.96*pred$se.fit)),ll=(pred$fit-(1.96*pred$se.fit)))

pd.rm.exbf2<-ggplot()+ geom_point(data = rmdat.rm, aes(x = age.sample, y = age.predicted, group = personid, colour=diarrhea))+
  geom_line(data = rmdat.rm, aes(x = age.sample, y = age.predicted, group = personid, colour=diarrhea),size=0.1)+
  geom_line(data = datfit,aes(x = age.sample, y = fit, colour=diarrhea),size = 1.5)+
  geom_ribbon(data = datfit,aes(x=age.sample, ymax=ul, ymin=ll,group=diarrhea, fill=diarrhea), alpha=.4)+guides(fill=FALSE)+
  scale_x_continuous(breaks=seq(from=0,to=24,by=3),
                     labels=seq(from=0,to=24,by=3))+
  theme(legend.position = c(0.15,0.95),legend.title=element_text(size=8),legend.text=element_text(size=8))+
  theme(legend.key.size = unit(0.3, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="white"))+
  xlab("Chronological age (months)") +ylab("Microbiota age (months)") +
  facet_wrap(~month.exbf2l, ncol = 1)
pd.rm.exbf2


#' GAM part
#'
kable(summary(b2$gam)$p.table)
kable(summary(b2$gam)$s.table)
#' LME part
#'
kable(summary(b2$lme)$tTable)
#' Test for heterogeneity (interaction)
#'
fit<-glmer(age.predicted~age.sample +month.exbf2*diarrhea+(1|personid),data=rmdat.rm)
tab<-summary(fit)$coefficients
tabz<-tab[,"Estimate"]/tab[,"Std. Error"]
tabp<-2*pnorm(-abs(tabz))
tab<-cbind(tab,p.val=tabp)
kable(tab)

#' ## Alpha diversity in samples of infants with vs. without diarrhea at sample collection stratified by duration of exclusive breastfeeding
#' With GAMM fit and 95%CI.
#'
#' ### Shannon
#use non-standardized alpha
load(paste(dir,"data/alphamean7.pooled.rda",sep="")) 
alpha.m.rm<-alpha.m.rm %>% group_by(personid) %>% arrange(personid,age.sample) %>%
  mutate(month.food6=cut(month.food, breaks=c(-Inf, 6, Inf), labels=c("<=6 months",">6 months")),
         month.food5=cut(month.food, breaks=c(-Inf, 5, Inf), labels=c("<=5 months",">5 months")),
         month.food4=cut(month.food, breaks=c(-Inf, 4, Inf), labels=c("<=4 months",">4 months")),
         month.foodr=as.factor(sort(round(month.food,0))),
         month.exbf3=cut(month.exbf, breaks=c(-Inf, 3, Inf), labels=c("<=3 months",">3 months")),
         month.exbf2=cut(month.exbf, breaks=c(-Inf, 2, Inf), labels=c("<=2 months",">2 months")),
         month.exbf1=cut(month.exbf, breaks=c(-Inf, 1, Inf), labels=c("<=1 months",">1 months")),
         month.exbfr=as.factor(sort(round(month.exbf,0))))
alpha.m.rm$dia.exbf2<-paste(alpha.m.rm$diarrhea,alpha.m.rm$month.exbf2,sep=".")
alpha.m.rm<-as.data.frame(alpha.m.rm)
alpha.m.rm$dia.exbf2<-as.factor(alpha.m.rm$dia.exbf2)
alpha.m.rm$diarrhea<-as.factor(alpha.m.rm$diarrhea)
alpha.m.rm$month.exbf2l<-mapvalues(alpha.m.rm$month.exbf2,from=c("<=2 months",">2 months"),to=c("Duration EBF <=2 months","Duration EBF >2 months"))
b2<-gamm(shannon~s(age.sample,by=dia.exbf2) +dia.exbf2,family=gaussian,
         data=alpha.m.rm,random=list(personid=~1))
pred <- predict(b2$gam, newdata = alpha.m.rm,se.fit=TRUE)
kable(summary(b2$gam)$p.table)
kable(summary(b2$gam)$s.table)
kable(summary(b2$lme)$tTable)
datfit<-cbind(alpha.m.rm, fit=pred$fit,ul=(pred$fit+(1.96*pred$se.fit)),ll=(pred$fit-(1.96*pred$se.fit)))
pd.s.exbf2<-ggplot()+ geom_point(data = alpha.m.rm, aes(x = age.sample, y = shannon, group = personid, colour=diarrhea))+
  geom_line(data = alpha.m.rm, aes(x = age.sample, y = shannon, group = personid, colour=diarrhea),size=0.1)+
  geom_line(data = datfit,aes(x = age.sample, y = fit, colour=diarrhea),size = 1.5)+
  geom_ribbon(data = datfit,aes(x=age.sample, ymax=ul, ymin=ll,group=diarrhea, fill=diarrhea), alpha=.4)+guides(fill=FALSE)+
  theme(legend.position = c(0.15,0.95),legend.title=element_text(size=8),legend.text=element_text(size=8))+
  theme(legend.key.size = unit(0.3, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="white"))+
  scale_x_continuous(breaks=seq(from=0,to=24,by=3),
                     labels=seq(from=0,to=24,by=3))+
  xlab("Chronological age (months)") +ylab("Shannon index") +
  facet_wrap(~month.exbf2l, ncol = 1)
pd.s.exbf2

#combine 4 plots
#+ fig.height=5,fig.width=10
grid.arrange(pd.rm.exbf2,pd.s.exbf2,p.dia.exbf2.6plus.l5noleg$p,p.dia.bf.6plus.l5.nolegend$p,nrow=1)


#' GAM part
#'
kable(summary(b2$gam)$p.table)
kable(summary(b2$gam)$s.table)
#' LME part
#'
kable(summary(b2$lme)$tTable)
#' Test for heterogeneity (interaction)
#'
fit<-glmer(shannon~age.sample +month.exbf2*diarrhea+(1|personid),data=alpha.m.rm)
tab<-summary(fit)$coefficients
tabz<-tab[,"Estimate"]/tab[,"Std. Error"]
tabp<-2*pnorm(-abs(tabz))
tab<-cbind(tab,p.val=tabp)
kable(tab)


#' ### Observed_species
b2<-gamm(observed_species~s(age.sample,by=dia.exbf2) +dia.exbf2,family=gaussian,
         data=alpha.m.rm,random=list(personid=~1))
pred <- predict(b2$gam, newdata = alpha.m.rm,se.fit=TRUE)
kable(summary(b2$gam)$p.table)
kable(summary(b2$gam)$s.table)
kable(summary(b2$lme)$tTable)
datfit<-cbind(alpha.m.rm, fit=pred$fit,ul=(pred$fit+(1.96*pred$se.fit)),ll=(pred$fit-(1.96*pred$se.fit)))
pd.os<-ggplot()+ geom_point(data = alpha.m.rm, aes(x = age.sample, y = observed_species, group = personid, colour=diarrhea))+
  geom_line(data = alpha.m.rm, aes(x = age.sample, y = observed_species, group = personid, colour=diarrhea))+
  geom_line(data = datfit,aes(x = age.sample, y = fit, colour=diarrhea),size = 1.5)+
  geom_ribbon(data = datfit,aes(x=age.sample, ymax=ul, ymin=ll,group=diarrhea), alpha=.5)+
  theme(legend.key.size = unit(0.5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ # strip.background = element_blank() element_rect(fill="white")
  scale_x_continuous(breaks=seq(from=0,to=24,by=3),
                     labels=seq(from=0,to=24,by=3))+
  theme(legend.position = "none")+
  xlab("Chronological age") +ylab("Observed_species") +facet_grid(.~month.exbf2)
pd.os
#' GAM part
#'
kable(summary(b2$gam)$p.table)
kable(summary(b2$gam)$s.table)
#' LME part
#'
kable(summary(b2$lme)$tTable)
#' Test for heterogeneity (interaction)
#'
fit<-glmer(observed_species~age.sample +month.exbf2*diarrhea+(1|personid),data=alpha.m.rm)
tab<-summary(fit)$coefficients
tabz<-tab[,"Estimate"]/tab[,"Std. Error"]
tabp<-2*pnorm(-abs(tabz))
tab<-cbind(tab,p.val=tabp)
kable(tab)

#' ### Pd_whole_tree
b2<-gamm(pd_whole_tree~s(age.sample,by=dia.exbf2) +dia.exbf2,family=gaussian,
         data=alpha.m.rm,random=list(personid=~1))
pred <- predict(b2$gam, newdata = alpha.m.rm,se.fit=TRUE)
kable(summary(b2$gam)$p.table)
kable(summary(b2$gam)$s.table)
kable(summary(b2$lme)$tTable)
datfit<-cbind(alpha.m.rm, fit=pred$fit,ul=(pred$fit+(1.96*pred$se.fit)),ll=(pred$fit-(1.96*pred$se.fit)))
pd.wt<-ggplot()+ geom_point(data = alpha.m.rm, aes(x = age.sample, y = pd_whole_tree, group = personid, colour=diarrhea))+
  geom_line(data = alpha.m.rm, aes(x = age.sample, y = pd_whole_tree, group = personid, colour=diarrhea))+
  geom_line(data = datfit,aes(x = age.sample, y = fit, colour=diarrhea),size = 1.5)+
  geom_ribbon(data = datfit,aes(x=age.sample, ymax=ul, ymin=ll,group=diarrhea), alpha=.5)+
  theme(legend.key.size = unit(0.5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ # strip.background = element_blank() element_rect(fill="white")
  scale_x_continuous(breaks=seq(from=0,to=24,by=3),
                     labels=seq(from=0,to=24,by=3))+
  theme(legend.position = "none")+
  xlab("Chronological age") +ylab("Pd_whole_tree") +facet_grid(.~month.exbf2)
pd.wt
#' GAM part
#'
kable(summary(b2$gam)$p.table)
kable(summary(b2$gam)$s.table)
#' LME part
#'
kable(summary(b2$lme)$tTable)
#' Test for heterogeneity (interaction)
#'
fit<-glmer(pd_whole_tree~age.sample +month.exbf2*diarrhea+(1|personid),data=alpha.m.rm)
tab<-summary(fit)$coefficients
tabz<-tab[,"Estimate"]/tab[,"Std. Error"]
tabp<-2*pnorm(-abs(tabz))
tab<-cbind(tab,p.val=tabp)
kable(tab)
ggplot(data=alpha.m.rm,aes(x=age.sample, y=pd_whole_tree, colour=diarrhea, group=personid)) +geom_point()+geom_line() +
  geom_smooth(data=alpha.m.rm,aes(x=age.sample, y=pd_whole_tree,group=diarrhea))+facet_grid(.~month.exbf2)
#' ### Chao1
b2<-gamm(chao1~s(age.sample,by=dia.exbf2) +dia.exbf2,family=gaussian,
         data=alpha.m.rm,random=list(personid=~1))
pred <- predict(b2$gam, newdata = alpha.m.rm,se.fit=TRUE)
kable(summary(b2$gam)$p.table)
kable(summary(b2$gam)$s.table)
kable(summary(b2$lme)$tTable)
datfit<-cbind(alpha.m.rm, fit=pred$fit,ul=(pred$fit+(1.96*pred$se.fit)),ll=(pred$fit-(1.96*pred$se.fit)))
pd.chao<-ggplot()+ geom_point(data = alpha.m.rm, aes(x = age.sample, y = chao1, group = personid, colour=diarrhea))+
  geom_line(data = alpha.m.rm, aes(x = age.sample, y = chao1, group = personid, colour=diarrhea))+
  geom_line(data = datfit,aes(x = age.sample, y = fit, colour=diarrhea),size = 1.5)+
  geom_ribbon(data = datfit,aes(x=age.sample, ymax=ul, ymin=ll,group=diarrhea), alpha=.5)+
  theme(legend.key.size = unit(0.5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ # strip.background = element_blank() element_rect(fill="white")
  scale_x_continuous(breaks=seq(from=0,to=24,by=3),
                     labels=seq(from=0,to=24,by=3))+
  theme(legend.position = "bottom")+
  xlab("Chronological age") +ylab("Chao1") +facet_grid(.~month.exbf2)
pd.chao
#' GAM part
#'
kable(summary(b2$gam)$p.table)
kable(summary(b2$gam)$s.table)
#' LME part
#'
kable(summary(b2$lme)$tTable)
#' Test for heterogeneity (interaction)
#'
fit<-glmer(chao1~age.sample +month.exbf2*diarrhea+(1|personid),data=alpha.m.rm)
tab<-summary(fit)$coefficients
tabz<-tab[,"Estimate"]/tab[,"Std. Error"]
tabp<-2*pnorm(-abs(tabz))
tab<-cbind(tab,p.val=tabp)
kable(tab)

#' Combined graph
#+ fig.with=5,fig.height=10
grid.arrange(pd.os,pd.wt,pd.chao,nrow=3)

#' R session information
sessionInfo()
