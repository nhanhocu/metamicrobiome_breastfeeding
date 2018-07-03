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
library(wesanderson)
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
#' # Microbiota age
#' ## Analysis for microbiota age based on shared genera
#' microbiota age was predicted based on the Random Forest model using L6 (genus) taxa relative abundance with list of taxa shared by all 7 included studies.
#'

load(paste(dir,"data/rmdat.shareg7s.rda",sep="")) 
load(paste(dir,"data/rmshare.conbfg7s.rda",sep="")) 

#' ## Plot of microbiota age by study all ages
#' ### With Generalized additive mixed model (GAMM) fit and 95%CI.
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
  xlab("Chronological age (month)") +ylab("microbiota age (month)")+
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


#' ### Standardized microbiota age and GAMM fit comparison between bf group within each study
rmdat.all$bfe<-mapvalues(rmdat.all$bf,from=c("ExclusiveBF","Non_exclusiveBF","No_BF"),to=c("EBF","Non-EBF","Non-BF"))
rmdat.all.nona<-rmdat.all[!is.na(rmdat.all$bf),]
rmdat.all.nona$pop.bf<-paste(rmdat.all.nona$pop,rmdat.all.nona$bf,sep="_")
rmdat.all.nona$pop.bf<-as.factor(rmdat.all.nona$pop.bf)
b2<-gamm(age.predicteds~s(age.sample,by=pop.bf),family=gaussian,
         data=rmdat.all.nona,random=list(personid=~1))
pred <- predict(b2$gam, newdata = rmdat.all.nona,se.fit=TRUE)
datfit<-cbind(rmdat.all.nona, fit=pred$fit,ul=(pred$fit+(1.96*pred$se.fit)),ll=(pred$fit-(1.96*pred$se.fit)))
p.rms.bf<-ggplot()+ geom_point(data = subset(rmdat.all.nona,age.sample<=6), aes(x = age.sample, y = age.predicteds, group = personid, colour=bfe),size = 1)+
  scale_color_manual(values=c("#08519c","#fe9929","#bd0026")) +
  geom_line(data = subset(rmdat.all.nona,age.sample<=6), aes(x = age.sample, y = age.predicteds, group = personid, colour=bfe),size=0.1)+
  geom_line(data = subset(datfit,age.sample<=6),aes(x = age.sample, y = fit, colour=bfe),size = 1)+
  geom_ribbon(data = subset(datfit,age.sample<=6),aes(x=age.sample, ymax=ul, ymin=ll, fill=bfe), alpha=.6)+guides(fill=FALSE)+
  scale_fill_manual(values=c("#08519c","#fe9929","#bd0026")) + 
  xlab("Chronological age (month)") +ylab("Standardized microbiota age")+
  labs(color='')+
  theme(legend.position = "bottom",
        legend.text = element_text(colour="black", size = 10,face="bold"),
        legend.background = element_rect(fill=alpha('white', 0)),
        legend.key = element_rect(fill=alpha('white', 0)),
        #plot.background = element_blank(),
        #panel.background = element_blank()
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="white"),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=10, face="bold"))+
  facet_grid(.~ pop)
p.rms.bf



#' Grey color 
p.rms.bf.g<-ggplot()+ geom_point(data = subset(rmdat.all.nona,age.sample<=6), aes(x = age.sample, y = age.predicteds, group = personid, colour=bfe),size=1)+
  scale_color_grey(start=0, end=0.7) +
  geom_line(data = subset(rmdat.all.nona,age.sample<=6), aes(x = age.sample, y = age.predicteds, group = personid, colour=bfe),size=0.1)+
  geom_line(data = subset(datfit,age.sample<=6),aes(x = age.sample, y = fit, colour=bfe),size = 1)+
  geom_ribbon(data = subset(datfit,age.sample<=6),aes(x=age.sample, ymax=ul, ymin=ll, fill=bfe), alpha=.6)+guides(fill=FALSE)+
  scale_fill_grey(start=0, end=0.7) + 
  xlab("Chronological age (month)") +ylab("Standardized microbiota age")+
  labs(color='')+
  theme(legend.position = "bottom",
        legend.text = element_text(colour="black", size = 10,face="bold"),
        legend.background = element_rect(fill=alpha('white', 0)),
        legend.key = element_rect(fill=alpha('white', 0)),
        #plot.background = element_blank(),
        #panel.background = element_blank()
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="white"),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=10, face="bold"))+
  facet_grid(.~ pop)
p.rms.bf.g


#' ## Meta-analysis for samples in <= 6 months old infants
#' Meta-analysis models based on adjusted estimate (adjusted for age of infant at sample collection) and standard error from linear mixed effect models.
#'
#' ### Change of RM in non-exclusive breastfed (nebf) vs. exclusive breastfed (exbf)
#+ fig.width=10, fig.height=5
#fix UW year
rmshare.sum.6<-as.data.frame(rmshare.sum.6s)
rmshare.sum.6["uw","year"]<-"2018"
rm.nebf<-metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=rmshare.sum.6,sm="RD", backtransf=FALSE)
forest(rm.nebf,smlab="Standardized \n microbiota age difference",lwd=2,sortvar=rmshare.sum.6$pop)
rm.nebf
kable(cbind(study=rm.nebf$studlab,pval=rm.nebf$pval))


#' ### Sensitivity analysis
#' #### No Haiti data
#+ fig.width=10, fig.height=5
rmshare.sum.6.noha<-rmshare.sum.6[!rownames(rmshare.sum.6) %in% "ha",]
rm.nebf<-metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=rmshare.sum.6.noha,sm="RD", backtransf=FALSE)
forest(rm.nebf,smlab="Standardized \n microbiota age difference",lwd=2,sortvar=rmshare.sum.6.noha$pop)
rm.nebf
kable(cbind(study=rm.nebf$studlab,pval=rm.nebf$pval))
#' #### No USA(NC) data
#+ fig.width=10, fig.height=5
rmshare.sum.6.nounc<-rmshare.sum.6[!rownames(rmshare.sum.6) %in% "unc",]
rm.nebf<-metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=rmshare.sum.6.nounc,sm="RD", backtransf=FALSE)
forest(rm.nebf,smlab="Standardized \n microbiota age difference",lwd=2,sortvar=rmshare.sum.6.nounc$pop)
rm.nebf
kable(cbind(study=rm.nebf$studlab,pval=rm.nebf$pval))
#' #### No USA(CA_MA_MO) data
#+ fig.width=10, fig.height=5
rmshare.sum.6.nohav<-rmshare.sum.6[!rownames(rmshare.sum.6) %in% "hav",]
rm.nebf<-metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=rmshare.sum.6.nohav,sm="RD", backtransf=FALSE)
forest(rm.nebf,smlab="Standardized \n microbiota age difference",lwd=2,sortvar=rmshare.sum.6.nohav$pop)
rm.nebf
kable(cbind(study=rm.nebf$studlab,pval=rm.nebf$pval))

#' ### Stratify by birth mode (not standardized)
load(paste(dir,"data/rmshare.vagcs.rda",sep="")) 
#' #### Vaginal
#+ fig.width=10, fig.height=5
rm.nebf<-metagen(estimate.nebf, se.nebf, studlab=study,data=rmshare.vag,sm="RD", backtransf=FALSE)
forest(rm.nebf,lwd=2,smlab="Microbiota age difference")
rm.nebf
kable(cbind(study=rm.nebf$studlab,pval=rm.nebf$pval))
#' #### C-section
#+ fig.width=10, fig.height=5
rm.nebf<-metagen(estimate.nebf, se.nebf, studlab=study,data=rmshare.cs,sm="RD", backtransf=FALSE)
forest(rm.nebf,lwd=2,smlab="Microbiota age difference")
rm.nebf
kable(cbind(study=rm.nebf$studlab,pval=rm.nebf$pval))


#' ### Trend of microbiota age in exclusive breastfed (exbf), non-exclusive breastfed (nebf) and no bf
#+ fig.width=10, fig.height=5
rmshare.conbf<-rmshare.conbfs
rmshare.conbf$pop<-c("Bangladesh","USA(CA_FL)","USA(CA_MA_MO)","Canada","USA(NC)")
rmshare.conbf<-rmshare.conbf[order(rmshare.conbf$pop),]
rm.conbf<-metagen(estimate.conbf, se.conbf, studlab=study,data=rmshare.conbf,sm="RD", backtransf=FALSE)
forest(rm.conbf,smlab="Standardized \n microbiota age difference",lwd=2,sortvar=rmshare.conbf$pop)
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
  scale_color_manual("Duration EBF",values=c("<=2 months"="#bd0026",">2 months"="#08519c"),labels=c("<=2 months"="<=2 months",">2 months"=">2 months"))+
  #scale_color_grey(start=0, end=0.5) +
  #geom_line(data = rmdat.rm, aes(x = age.sample, y = age.predicted, group = personid, colour=month.exbf2),size=0.3)+
  geom_line(data = datfit,aes(x = age.sample, y = fit, colour=month.exbf2),size = 1.5)+
  geom_ribbon(data = datfit,aes(x=age.sample, ymax=ul, ymin=ll,group=month.exbf2, fill=month.exbf2), alpha=.6)+guides(fill=FALSE)+
  scale_fill_manual("Duration EBF",values=c("<=2 months"="#bd0026",">2 months"="#08519c"),labels=c("<=2 months"="<=2 months",">2 months"=">2 months"))+
  #scale_fill_grey(start=0, end=0.5) + 
  #theme(legend.position = "bottom")+
  #theme(legend.title = element_text(colour="black", size=10))+
  labs(color  = "Duration EBF", shape = "Duration EBF")+
  scale_x_continuous(breaks=seq(from=0,to=24,by=3),
                     labels=seq(from=0,to=24,by=3))+
  #theme(legend.text = element_text(colour="black", size = 10))+
  theme(legend.position = c(0.15,0.85),legend.title=element_text(size=10,face="bold"),legend.text=element_text(size=10,face="bold"))+
  theme(legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(fill=alpha('white', 0)),
        legend.key = element_rect(fill=alpha('white', 0)),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face="bold", size=10),
        axis.title.y = element_text(face="bold", size=10),
        axis.text.x = element_text(face="bold", size=10, color = "black"),
        axis.text.y = element_text(face="bold", size=10, color = "black"))+ 
  xlab("Chronological age (month)") +ylab("Microbiota age (month)")
pexbf2



#' Grey
pexbf2.g<-ggplot()+ geom_point(data = rmdat.rm, aes(x = age.sample, y = age.predicted, group = personid, colour=month.exbf2))+
  scale_color_grey(start=0, end=0.5) +
  #geom_line(data = rmdat.rm, aes(x = age.sample, y = age.predicted, group = personid, colour=month.exbf2),size=0.3)+
  geom_line(data = datfit,aes(x = age.sample, y = fit, colour=month.exbf2),size = 1.5)+
  geom_ribbon(data = datfit,aes(x=age.sample, ymax=ul, ymin=ll,group=month.exbf2, fill=month.exbf2), alpha=.6)+guides(fill=FALSE)+
  scale_fill_grey(start=0, end=0.5) + 
  #theme(legend.position = "bottom")+
  #theme(legend.title = element_text(colour="black", size=10))+
  labs(color='Duration EBF')+
  scale_x_continuous(breaks=seq(from=0,to=24,by=3),
                     labels=seq(from=0,to=24,by=3))+
  #theme(legend.text = element_text(colour="black", size = 10))+
  theme(legend.position = c(0.15,0.85),legend.title=element_text(size=10,face="bold"),legend.text=element_text(size=10,face="bold"))+
  theme(legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(fill=alpha('white', 0)),
        legend.key = element_rect(fill=alpha('white', 0)),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        axis.text.x = element_text(face="bold", size=10, color = "black"),
        axis.text.y = element_text(face="bold", size=10, color = "black"))+ 
  xlab("Chronological age (month)") +ylab("Microbiota age (month)")
pexbf2.g

#' Test for age > 6 months and <15months
#' GAM part
#'
rmdat.rm.6plus15<-rmdat.rm[rmdat.rm$age.sample>6 & rmdat.rm$age.sample<15,]
b2<-gamm(age.predicted~s(age.sample,by=month.exbf2) +month.exbf2,family=gaussian,
         data=rmdat.rm.6plus15,random=list(personid=~1))
kable(summary(b2$gam)$p.table)
kable(summary(b2$gam)$s.table)
#' LME part
#'
kable(summary(b2$lme)$tTable)



#' ## Performance of the Random Forest model to estimate microbiota age
load(paste(dir,"data/SrfFit.rml6.shareg7.train.test.rda",sep="")) 
#' ### Evaluated on the training and test set of Subramanian (Bangladesh) data.
Straining$age.predicted <- predict(SrfFit.rml6.share, newdata = Straining)
actual<-Straining$Sage
predict<-Straining$age.predicted
R2 <- 1 - (sum((actual-predict )^2)/sum((actual-mean(actual))^2))
R2<-round(R2,2)
ptrain<-ggplot() +geom_point(data=Straining,aes(x=Sage, y=age.predicted))+
  theme(legend.text = element_text(colour="black", size = 10))+
  annotate("text", x=15, y=5,label=paste("R squared =",R2,sep=" "), colour="black",size=5) +
  labs(title="Training set")+
  theme(legend.key.size = unit(0.5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=10,face="bold"),
        strip.text.x = element_text(size=10, face="bold"))+
  xlab("Chronological age (month)") +ylab("Microbiota age (month)")

actual<-rmdat.rm$age.sample
predict<-rmdat.rm$age.predicted
R2 <- 1 - (sum((actual-predict )^2)/sum((actual-mean(actual))^2))
R2<-round(R2,2)
ptest<-ggplot() +geom_point(data=rmdat.rm,aes(x=age.sample, y=age.predicted))+
  theme(legend.text = element_text(colour="black", size = 10))+
  annotate("text", x=15, y=5,label=paste("R squared =",R2,sep=" "), colour="black",size=5) +
  labs(title="Test set")+
  theme(legend.key.size = unit(0.5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=10,face="bold"),
        strip.text.x = element_text(size=10, face="bold"))+
  xlab("Chronological age (month)") +ylab("Microbiota age (month)")
grid.arrange(ptrain, ptest,nrow=1)


#' ### List of shared taxa and their relative importance
library(randomForest)
taxim<-as.data.frame(importance(SrfFit.rml6.share$finalModel))
taxim$genera<-rownames(taxim)
taxim$importance<-taxim[,"IncNodePurity"]
taxim<-taxim[order(taxim[,"IncNodePurity"],decreasing = TRUE),]
taxim<-taxim[,c("genera","importance")]
taxim$importance.percent<-(taxim$importance/sum(taxim$importance))*100
taxim[,c("importance","importance.percent")]<-round(taxim[,c("importance","importance.percent")],2)
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
alphap$bfe<-mapvalues(alphap$bf,from=c("ExclusiveBF","Non_exclusiveBF","No_BF"),to=c("EBF","Non-EBF","Non-BF"))
alphap.nona<-alphap[!is.na(alphap$bf),]
alphap.nona$pop.bf<-paste(alphap.nona$pop,alphap.nona$bf,sep="_")
alphap.nona$pop.bf<-as.factor(alphap.nona$pop.bf)
b2<-gamm(shannon~s(age.sample,by=pop.bf),family=gaussian,
         data=alphap.nona,random=list(personid=~1))
pred <- predict(b2$gam, newdata = alphap.nona,se.fit=TRUE)
datfit<-cbind(alphap.nona, fit=pred$fit,ul=(pred$fit+(1.96*pred$se.fit)),ll=(pred$fit-(1.96*pred$se.fit)))
#' Samples <=6 months only
p.sha.rm.bf6<-ggplot()+ geom_point(data = subset(alphap.nona,age.sample<=6), aes(x = age.sample, y = shannon, group = personid, colour=bfe), size=1)+ #, colour=bfe
  scale_color_manual(values=c("#08519c","#fe9929","#bd0026")) +
  geom_line(data = subset(alphap.nona,age.sample<=6), aes(x = age.sample, y = shannon, group = personid, colour=bfe),size=0.1)+ #, colour=bfe
  geom_line(data = subset(datfit,age.sample<=6),aes(x = age.sample, y = fit, colour=bfe),size = 1)+
  geom_ribbon(data = subset(datfit,age.sample<=6),aes(x=age.sample, ymax=ul, ymin=ll, fill=bfe), alpha=.6)+ guides(fill=FALSE)+
  scale_fill_manual(values=c("#08519c","#fe9929","#bd0026")) + 
  theme(legend.position = "bottom")+
  theme(legend.title = element_text(colour="black", size=10))+
  labs(color='')+
  scale_x_continuous(breaks=seq(from=0,to=24,by=2),
                     labels=seq(from=0,to=24,by=2))+
  theme(legend.text = element_text(colour="black", size = 10))+
  theme(legend.key.size = unit(0.5, "cm"))+
  theme(legend.position = "bottom",
        legend.text = element_text(colour="black", size = 10,face="bold"),
        legend.background = element_rect(fill=alpha('white', 0)),
        legend.key = element_rect(fill=alpha('white', 0)),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="white"),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=10, face="bold"))+
  xlab("Chronological age (month)") +ylab("Standardized Shannon index")+
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
#fix year UW
alphaall[alphaall$author=="Wood et al","year"]<-"2018"
chao1.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="chao1"),sm="RD", backtransf=FALSE)
forest(chao1.nebf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="chao1")$pop, lwd=2)
chao1.nebf
kable(cbind(study=chao1.nebf$studlab,pval=chao1.nebf$pval))
#' #### Observed_species
#+ fig.width=12, fig.height=5
os.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="observed_species"),sm="RD", backtransf=FALSE)
forest(os.nebf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="observed_species")$pop, lwd=2)
os.nebf
kable(cbind(study=os.nebf$studlab,pval=os.nebf$pval))
#' #### Pd_whole_tree
#+ fig.width=10, fig.height=5
pwt.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="pd_whole_tree"),sm="RD", backtransf=FALSE)
forest(pwt.nebf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="pd_whole_tree")$pop, lwd=2)
pwt.nebf
kable(cbind(study=pwt.nebf$studlab,pval=pwt.nebf$pval))
#' #### Shannon
#+ fig.width=10, fig.height=5
shannon.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="shannon"),sm="RD", backtransf=FALSE)
forest(shannon.nebf,smlab="Standardized \n diversity difference",lwd=2,sortvar=subset(alphaall,index=="shannon")$pop)
shannon.nebf
kable(cbind(study=shannon.nebf$studlab,pval=shannon.nebf$pval))

#' ### Sensitivity analysis Change in non-exbf vs. exbf
#' Show the results of Shannon indexes only.
#'
#' #### No Haiti data
#+ fig.width=10, fig.height=5
shannon.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="shannon"&(pop!="Haiti")),sm="RD", backtransf=FALSE)
forest(shannon.nebf,smlab="Standardized \n diversity difference",lwd=2,sortvar=subset(alphaall,index=="shannon"&(pop!="Haiti"))$pop)
shannon.nebf
kable(cbind(study=shannon.nebf$studlab,pval=shannon.nebf$pval))
#' #### No UNC data
#+ fig.width=10, fig.height=5
shannon.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="shannon"&(pop!="USA(NC)")),sm="RD", backtransf=FALSE)
forest(shannon.nebf,smlab="Standardized \n diversity difference",lwd=2,sortvar=subset(alphaall,index=="shannon"&(pop!="USA(NC)"))$pop)
shannon.nebf
kable(cbind(study=shannon.nebf$studlab,pval=shannon.nebf$pval))
#' #### No USA(CA_MA_MO) data
#+ fig.width=10, fig.height=5
shannon.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="shannon"&(pop!="USA(CA_MA_MO)")),sm="RD", backtransf=FALSE)
forest(shannon.nebf,smlab="Standardized \n diversity difference",lwd=2,sortvar=subset(alphaall,index=="shannon"&(pop!="USA(CA_MA_MO)"))$pop)
shannon.nebf
kable(cbind(study=shannon.nebf$studlab,pval=shannon.nebf$pval))

#' ### Stratify by birth mode for Change in non-exbf vs. exbf (not standardized)
#' Show results of Shannon index only.
#'
load(paste(dir,"data/alphasum.vagcs.rda",sep="")) 
#' #### Vaginal
#+ fig.width=10, fig.height=5
shannon.nebf <- metagen(estimate.nebf, se.nebf, studlab=study,data=subset(alphasum.vag, index=="shannon"),sm="RD", backtransf=FALSE)
forest(shannon.nebf,smlab="Standardized \n diversity difference",lwd=2)
shannon.nebf
kable(cbind(study=shannon.nebf$studlab,pval=shannon.nebf$pval))
#' #### C-section
#+ fig.width=10, fig.height=5
shannon.nebf <- metagen(estimate.nebf, se.nebf, studlab=study,data=subset(alphasum.cs, index=="shannon"),sm="RD", backtransf=FALSE)
forest(shannon.nebf,smlab="Standardized \n diversity difference",lwd=2)
shannon.nebf
kable(cbind(study=shannon.nebf$studlab,pval=shannon.nebf$pval))


#' #### Put meta-analysis results (random models) of all indexes together
a.nebf<-ggplot(data=a.metatab.r,aes(x=estimate.nebf,y=index))+
  geom_point(shape=17, colour="black",size=2)+
  geom_errorbarh(aes(xmin=ll.nebf,xmax=ul.nebf),height=0.0, colour="black",size=1)+
  scale_x_continuous(breaks=seq(from=0,to=1,by=0.2),
                     labels=seq(from=0,to=1,by=0.2))+
  geom_vline(xintercept=0,linetype="dashed")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=10, face="bold"),
        plot.title = element_text(size=14,face="bold"))+
  ggtitle("Non-EBF vs. EBF")+
  xlab("Pooled Standardized diversity difference")+ylab("Alpha diversity index")
a.nebf
kable(a.metatab.r[,grep(".nebf",colnames(a.metatab.r))])

#' ### Trend effect in exbf, non-exbf, no bf
#' #### Chao1
#+ fig.width=10, fig.height=5
chao1.conbf <- metagen(estimate.conbf, se.conbf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="chao1"),sm="RD", backtransf=FALSE)
forest(chao1.conbf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="chao1")$pop,lwd=2)
chao1.conbf
kable(cbind(study=chao1.conbf$studlab,pval=chao1.conbf$pval))
#' #### Observed_species
#+ fig.width=10, fig.height=5
os.conbf <- metagen(estimate.conbf, se.conbf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="observed_species"),sm="RD", backtransf=FALSE)
forest(os.conbf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="observed_species")$pop,lwd=2)
os.conbf
kable(cbind(study=os.conbf$studlab,pval=os.conbf$pval))
#' #### Pd_whole_tree
#+ fig.width=10, fig.height=5
pwt.conbf <- metagen(estimate.conbf, se.conbf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="pd_whole_tree"),sm="RD", backtransf=FALSE)
forest(pwt.conbf,smlab="Standardized \n diversity difference",sortvar=subset(alphaall,index=="pd_whole_tree")$pop,lwd=2)
pwt.conbf
kable(cbind(study=pwt.conbf$studlab,pval=pwt.conbf$pval))
#' #### Shannon
#+ fig.width=10, fig.height=5
shannon.conbf <- metagen(estimate.conbf, se.conbf, studlab=paste(author,year,"(",pop,")"),data=subset(alphaall,index=="shannon"&(!pop %in% c("Haiti","South Africa"))),sm="RD", backtransf=FALSE)
forest(shannon.conbf,smlab="Standardized \n diversity difference",lwd=2,sortvar=subset(alphaall,index=="shannon"&(!pop %in% c("Haiti","South Africa")))$pop)
shannon.conbf
kable(cbind(study=shannon.conbf$studlab,pval=shannon.conbf$pval))

#' #### Put meta-analysis results (random models) of all indexes together
a.conbf<-ggplot(data=a.metatab.r,aes(x=estimate.conbf,y=index))+
  geom_point(shape=17, colour="black",size=2)+
  geom_errorbarh(aes(xmin=ll.conbf,xmax=ul.conbf),height=0.0, colour="black",size=1)+
  geom_vline(xintercept=0,linetype="dashed")+
  scale_x_continuous(breaks=seq(from=0,to=1,by=0.2),
                     labels=seq(from=0,to=1,by=0.2))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=10, face="bold"),
        plot.title = element_text(size=14,face="bold"))+
  ggtitle("Trend across EBF, non-EBF, non-BF")+
  xlab("Pooled standardized diversity difference")+ylab("Alpha diversity index")
a.conbf
kable(a.metatab.r[,grep(".conbf",colnames(a.metatab.r))])

#grid.arrange(a.nebf,a.conbf,nrow=1)



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
meta.niceplot(metadat=metadat,sumtype="taxa",level="main",leg.key.size=1,leg.text.size=10,heat.text.x.size=4,forest.axis.text.y=6,forest.axis.text.x=4,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,1))

#' ### Order (l4)
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="heatmap",fill.value="log(OR)")
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="forest",fill.value="log(OR)")
kable(metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' Nice plot
#+ fig.width=10, fig.height=5
metadat<-metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",leg.key.size=1,leg.text.size=10,heat.text.x.size=4,forest.axis.text.y=6,forest.axis.text.x=4,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,1))


#' ### Family (l5)
#+ fig.width=8, fig.height=10
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="heatmap",fill.value="log(OR)")
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="forest",fill.value="log(OR)")
kable(metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' Nice plot
#+ fig.width=10, fig.height=10
metadat<-metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",leg.key.size=1,leg.text.size=10,heat.text.x.size=4,forest.axis.text.y=6,forest.axis.text.x=4,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,1))


#' ### Genus (l6)
#+ fig.width=8, fig.height=12
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="heatmap",fill.value="log(OR)")
metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,plot="forest",fill.value="log(OR)")
kable(metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' Nice plot
#+ fig.width=10, fig.height=12
metadat<-metatab.show(metatab=metatab.zi$random,taxacom.pooled.tab=taxacom.zi,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",leg.key.size=1,leg.text.size=10,heat.text.x.size=4,forest.axis.text.y=6,forest.axis.text.x=4,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,0.8))


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
#+ fig.width=15, fig.height=5
kable(metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
metadat$taxsig.all<-metadat$taxsig.all[metadat$taxsig.all$taxa!=".thermi.",]
metadat$taxsig<-metadat$taxsig[metadat$taxsig$taxa!=".thermi.",]
meta.niceplot(metadat=metadat,sumtype="taxa",level="main",phyla.col="rainbow",leg.key.size=1,leg.text.size=8,heat.text.x.size=8,forest.axis.text.y=8,forest.axis.text.x=8,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,1))

#' #### l4
#+ fig.width=15, fig.height=5
kable(metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",phyla.col="rainbow",leg.key.size=1,leg.text.size=10,heat.text.x.size=8,forest.axis.text.y=8,forest.axis.text.x=8,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,1))

#' #### l5
#+ fig.width=15, fig.height=10
kable(metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",phyla.col="rainbow",leg.key.size=1,leg.text.size=10,heat.text.x.size=8,forest.axis.text.y=8,forest.axis.text.x=8,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,1))

#' #### l6
#+ fig.width=15, fig.height=12
kable(metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.vag$random,taxacom.pooled.tab=taxacom.zi.vag,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",phyla.col="rainbow",leg.key.size=1,leg.text.size=10,heat.text.x.size=8,forest.axis.text.y=8,forest.axis.text.x=8,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,1))


#' ### C-section
#' #### l2
#+ fig.width=15, fig.height=5
kable(metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#tiff("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/fig/meta.relabund.nebf.l2.cs.rev2.tiff",width =10,height =1.5,units= "in",res= 300)
#pdf("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/fig/meta.relabund.nebf.l2.cs.rev2.pdf",width =10,height =1.5)
metadat<-metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="main",phyla.col="rainbow",leg.key.size=1,leg.text.size=10,heat.text.x.size=8,forest.axis.text.y=8,forest.axis.text.x=8,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,1))
#dev.off()

#' #### l4
#+ fig.width=15, fig.height=5
kable(metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",phyla.col="rainbow",leg.key.size=1,leg.text.size=10,heat.text.x.size=8,forest.axis.text.y=8,forest.axis.text.x=8,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,1))


#' #### l5
#+ fig.width=15, fig.height=10
kable(metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",phyla.col="rainbow",leg.key.size=1,leg.text.size=10,heat.text.x.size=8,forest.axis.text.y=8,forest.axis.text.x=8,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,1))


#' #### l6
#+ fig.width=15, fig.height=12
kable(metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.cs$random,taxacom.pooled.tab=taxacom.zi.cs,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",phyla.col="rainbow",leg.key.size=1,leg.text.size=10,heat.text.x.size=8,forest.axis.text.y=8,forest.axis.text.x=8,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,1))



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
#+ fig.width=15, fig.height=8
metadat<-metatab.show(metatab=pathmetatab.zi$random,taxacom.pooled.tab=pathcom.zi,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=10,heat.text.x.size=6,forest.axis.text.y=6,forest.axis.text.x=6,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,1))

#' ### Level 3 KEGG pathway
kable(metatab.show(metatab=pathmetatab.zi$random,taxacom.pooled.tab=pathcom.zi,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#' Nice plot all pathways
#+ fig.width=10, fig.height=15
metadat<-metatab.show(metatab=pathmetatab.zi$random,taxacom.pooled.tab=pathcom.zi,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=4,forest.axis.text.x=4)

#' Nice plot significant pathways only (pooled p<0.05)
#+ fig.width=15, fig.height=5
metadat<-metatab.show(metatab=pathmetatab.zi$random,taxacom.pooled.tab=pathcom.zi,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=6,forest.axis.text.x=6)

#' Nice plot multiple testing adjusted significant pathways only (adjusted pooled p<0.1)
#+ fig.width=15, fig.height=5
metadat<-metatab.show(metatab=pathmetatab.zi$random,taxacom.pooled.tab=pathcom.zi,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p.adjust", p.cutoff=0.1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=10,heat.text.x.size=6,forest.axis.text.y=6,forest.axis.text.x=6,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,1))


#' ### Sensitivity analysis
load(paste(dir,"data/pathmetatab.zi.sen.rda",sep="")) 
#' #### No Haiti data
#' Level 2
#+ fig.width=10, fig.height=7
kable(metatab.show(metatab=pathmetatab.zi.noha$random,taxacom.pooled.tab=pathcom.zi.noha,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=pathmetatab.zi.noha$random,taxacom.pooled.tab=pathcom.zi.noha,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=6,forest.axis.text.x=7)
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
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=6,forest.axis.text.x=7)
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
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=6,forest.axis.text.x=7)
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
#+ fig.width=15, fig.height=6
metadat<-metatab.show(metatab=pathmetatab.zi.vag$random,taxacom.pooled.tab=pathcom.zi.vag,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p.adjust", p.cutoff=0.1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=10,heat.text.x.size=6,forest.axis.text.y=6,forest.axis.text.x=6,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,0.8))


#' ### C-section birth
#' Level 2
#+ fig.width=10, fig.height=7
kable(metatab.show(metatab=pathmetatab.zi.cs$random,taxacom.pooled.tab=pathcom.zi.cs,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=pathmetatab.zi.cs$random,taxacom.pooled.tab=pathcom.zi.cs,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=7,forest.axis.text.x=7)
#' Level 3
#' All
#+ fig.width=12, fig.height=15
metadat<-metatab.show(metatab=pathmetatab.zi.cs$random,taxacom.pooled.tab=pathcom.zi.cs,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=4,forest.axis.text.y=4,forest.axis.text.x=4)
#' Significant only (pooled p<0.05)
#+ fig.width=12, fig.height=5
kable(metatab.show(metatab=pathmetatab.zi.cs$random,taxacom.pooled.tab=pathcom.zi.cs,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=pathmetatab.zi.cs$random,taxacom.pooled.tab=pathcom.zi.cs,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=8,heat.text.x.size=6,forest.axis.text.y=5,forest.axis.text.x=6)

#' Multiple testing adjusted Significant only (adjusted pooled p<0.1)
#+ fig.width=12, fig.height=5
metadat<-metatab.show(metatab=pathmetatab.zi.cs$random,taxacom.pooled.tab=pathcom.zi.cs,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p.adjust", p.cutoff=0.1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=10,heat.text.x.size=6,forest.axis.text.y=6,forest.axis.text.x=6,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,0.8))



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
kable(taxcomtab.show(taxcomtab=taxacom.6plus.exbf2.zi.rm,tax.select=p.exbf.l2$taxuse.rm, tax.lev="l2",p.adjust.method="fdr",readjust.p=TRUE))

#' #### Order (l4)
p.exbf.l4<-taxa.mean.plot(tabmean=taxa.meansdn.exbf2.rm,taxlist=taxlist.rm,tax.lev="l4", comvar="month.exbf2", groupvar="age.sample",mean.filter=0.005)
p.exbf.l4$p
#' GAMLSS
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.exbf2.zi.rm,tax.select=p.exbf.l4$taxuse.rm, tax.lev="l4",p.adjust.method="fdr",readjust.p=TRUE))

#' #### Family (L5)
#+ fig.width=10, fig.height=7
p.exbf.l5<-taxa.mean.plot(tabmean=taxa.meansdn.exbf2.rm,taxlist=taxlist.rm,tax.lev="l5", comvar="month.exbf2", groupvar="age.sample",mean.filter=0.005,show.taxname="short",legend.position="right")
# better plot view
for (i in 1: length(taxa.meansdn.exbf2.rm)){
  taxa.meansdn.exbf2.rm[[i]]$month.exbf2l<-mapvalues(taxa.meansdn.exbf2.rm[[i]]$month.exbf2,from=c("<=2 months",">2 months"),to=c("Duration EBF <=2 months","Duration EBF >2 months"))
}
p.exbf.l5.noleg<-taxa.mean.plot(tabmean=taxa.meansdn.exbf2.rm,taxlist=taxlist.rm,tax.lev="l5", comvar="month.exbf2l", groupvar="age.sample",mean.filter=0.005,legend.position="none")
p.exbf.l5<-taxa.mean.plot(tabmean=taxa.meansdn.exbf2.rm,taxlist=taxlist.rm,tax.lev="l5", comvar="month.exbf2l", groupvar="age.sample",mean.filter=0.005,legend.position="right",show.taxname = "short")
p.exbf.l5$p

#tiff("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/fig/relabundrm.l5.age.rm.exbf2.rev.tiff",width =10,height =4,units= "in",res= 300)
#pdf("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/fig/relabundrm.l5.age.rm.exbf2.rev.pdf",width =10,height =4)
grid.arrange(pexbf2,p.exbf.l5.noleg$p,nrow=1)
#dev.off()

#' Grey
#+ fig.width=10, fig.height=7
#tiff("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/fig/relabundrm.l5.age.rm.exbf2.grey.tiff",width =10,height =4,units= "in",res= 300)
#pdf("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/fig/relabundrm.l5.age.rm.exbf2.grey.pdf",width =10,height =4)
grid.arrange(pexbf2.g,p.exbf.l5.noleg$p,nrow=1)
#dev.off()



#' GAMLSS
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.exbf2.zi.rm,tax.select=p.exbf.l5$taxuse.rm, tax.lev="l5",p.adjust.method="fdr",readjust.p=TRUE))

#' #### Genus (L6)
#+ fig.width=15, fig.height=10
p.exbf.l6<-taxa.mean.plot(tabmean=taxa.meansdn.exbf2.rm,taxlist=taxlist.rm,tax.lev="l6", comvar="month.exbf2", groupvar="age.sample",mean.filter=0.005)
p.exbf.l6$p
#' GAMLSS
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.exbf2.zi.rm,tax.select=p.exbf.l6$taxuse.rm, tax.lev="l6",p.adjust.method="fdr",readjust.p=TRUE))


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
#more detail legend
teste<-taxa.meansdn.dia.exbf2.6plus.rm
for (i in 1:length(names(teste))){
  teste[[i]]$month.exbf2l<-mapvalues(teste[[i]]$month.exbf2,from=c("<=2 months",">2 months"),to=c("Duration EBF <=2 months","Duration EBF >2 months"))
  teste[[i]]$Diarrhea<-teste[[i]]$diarrhea
}
p.dia.exbf2.6plus.l5noleg<-taxa.mean.plot(tabmean=teste,taxlist=taxlist.rm,tax.lev="l5", comvar="Diarrhea", groupvar="month.exbf2l",mean.filter=0.005,legend.position="none",ylab="Relative abundance (6 months - 2 years)")
p.dia.exbf2.6plus.l5<-taxa.mean.plot(tabmean=teste,taxlist=taxlist.rm,tax.lev="l5", comvar="Diarrhea", groupvar="month.exbf2l",mean.filter=0.005,legend.position="right", show.taxname = "short",ylab="Relative abundance (6 months - 2 years)")
p.dia.exbf2.6plus.l5$p

#' ###### GAMLSS
#' In infants with duration of exclusive bf <=2 months
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.exbf2.zi.rm, tax.lev="l5",tax.select=p.dia.exbf2.6plus.l5$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1,readjust.p=TRUE))
#' In infants with duration of exclusive bf >2 months
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.exbf2plus.zi.rm, tax.lev="l5",tax.select=p.dia.exbf2.6plus.l5$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1,readjust.p=TRUE))


#' ###### Genus
#+ fig.width=15, fig.height=5
p.dia.exbf2.6plus.l6<-taxa.mean.plot(tabmean=taxa.meansdn.dia.exbf2.6plus.rm,taxlist=taxlist.rm,tax.lev="l6", comvar="diarrhea", groupvar="month.exbf2",mean.filter=0.005,legend.position = "right")
p.dia.exbf2.6plus.l6$p
#' ####### GAMLSS
#' In infants with duration of exclusive bf <=2 months
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.exbf2.zi.rm, tax.lev="l6",tax.select=p.dia.exbf2.6plus.l6$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1,readjust.p=TRUE))
#' In infants with duration of exclusive bf >2 months
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.exbf2plus.zi.rm, tax.lev="l6",tax.select=p.dia.exbf2.6plus.l6$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1,readjust.p=TRUE))


#' #### Stratified by bf status
#'
#+ fig.width=15, fig.height=5
taxa.meansdn.dia.bf.6plus.rm.s<-list()
for (i in 1:5){
  taxa.meansdn.dia.bf.6plus.rm.s[[i]]<-taxa.meansdn.dia.bf.6plus.rm[[i]][taxa.meansdn.dia.bf.6plus.rm[[i]]$bf!="ExclusiveBF",]
  taxa.meansdn.dia.bf.6plus.rm.s[[i]]$bf<-drop.levels(taxa.meansdn.dia.bf.6plus.rm.s[[i]]$bf,reorder=FALSE)
  taxa.meansdn.dia.bf.6plus.rm.s[[i]]$Diarrhea<-taxa.meansdn.dia.bf.6plus.rm.s[[i]]$diarrhea
}
names(taxa.meansdn.dia.bf.6plus.rm.s)<-names(taxa.meansdn.dia.bf.6plus.rm)

#' ##### Family
p.dia.bf.6plus.l5<-taxa.mean.plot(tabmean=taxa.meansdn.dia.bf.6plus.rm.s,tax.select=p.dia.exbf2.6plus.l5$taxuse.rm,taxlist=taxlist.rm,tax.lev="l5", comvar="Diarrhea", groupvar="bf",mean.filter=0.005,show.taxname = "short")
p.dia.bf.6plus.l5$p

# reverse levels for better combined plot view
testr=taxa.meansdn.dia.bf.6plus.rm.s
for (i in 1: length(names(testr))){
  testr[[i]]$bfl<-mapvalues(testr[[i]]$bf,from=c("No_BF","Non_exclusiveBF"),to=c("No BF when diarrhea","BF when diarrhea"))
  testr[[i]]$bfl<-factor(testr[[i]]$bfl,levels=c("No BF when diarrhea","BF when diarrhea"))
}
p.dia.bf.6plus.l5.nolegend<-taxa.mean.plot(tabmean=testr,tax.select=p.dia.exbf2.6plus.l5$taxuse.rm,taxlist=taxlist.rm,tax.lev="l5", comvar="Diarrhea", groupvar="bfl",mean.filter=0.005,show.taxname = "short",legend.position = "none",ylab="Relative abundance (6 months - 2 years)")

#' ###### GAMLSS
#' In non-exclusive bf infants
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.nexbf.zi.rm, tax.lev="l5",tax.select=p.dia.bf.6plus.l5$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1,readjust.p=TRUE))
#' In no bf infants
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.nobf.zi.rm, tax.lev="l5",tax.select=p.dia.bf.6plus.l5$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1,readjust.p=TRUE))

#' ##### Genus
#+ fig.width=10, fig.height=5
p.dia.bf.6plus.l6<-taxa.mean.plot(tabmean=taxa.meansdn.dia.bf.6plus.rm.s,tax.select=p.dia.exbf2.6plus.l6$taxuse.rm,taxlist=taxlist.rm,tax.lev="l6", comvar="diarrhea", groupvar="bf",mean.filter=0.005)
p.dia.bf.6plus.l6$p
#' ###### GAMLSS
#' In bf infants
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.nexbf.zi.rm, tax.lev="l6",tax.select=p.dia.bf.6plus.l6$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1,readjust.p=TRUE))
#' In non-bf infants
#+ results = 'asis'
kable(taxcomtab.show(taxcomtab=taxacom.6plus.dia.nobf.zi.rm, tax.lev="l6",tax.select=p.dia.bf.6plus.l6$taxuse.rm,p.adjust.method="fdr",p.cutoff=0.1,readjust.p=TRUE))



#' ## microbiota age in samples of infants with vs. without diarrhea at sample collection stratified by duration of exclusive breastfeeding
#' With GAMM fit and 95%CI.
#'
rmdat.rm$dia.exbf2<-paste(rmdat.rm$diarrhea,rmdat.rm$month.exbf2,sep=".")
rmdat.rm<-as.data.frame(rmdat.rm)
rmdat.rm$dia.exbf2<-as.factor(rmdat.rm$dia.exbf2)
rmdat.rm$Diarrhea<-as.factor(rmdat.rm$diarrhea)
rmdat.rm$month.exbf2l<-mapvalues(rmdat.rm$month.exbf2,from=c("<=2 months",">2 months"),to=c("Duration EBF <=2 months","Duration EBF >2 months"))
b2<-gamm(age.predicted~s(age.sample,by=dia.exbf2) +dia.exbf2,family=gaussian,
         data=rmdat.rm,random=list(personid=~1))
pred <- predict(b2$gam, newdata = rmdat.rm,se.fit=TRUE)
datfit<-cbind(rmdat.rm, fit=pred$fit,ul=(pred$fit+(1.96*pred$se.fit)),ll=(pred$fit-(1.96*pred$se.fit)))

pd.rm.exbf2<-ggplot()+ geom_point(data = rmdat.rm, aes(x = age.sample, y = age.predicted, group = personid, colour=Diarrhea))+
  scale_color_grey(start=0.5, end=0) +
  #geom_line(data = rmdat.rm, aes(x = age.sample, y = age.predicted, group = personid, colour=Diarrhea),size=0.1)+
  geom_line(data = datfit,aes(x = age.sample, y = fit, colour=Diarrhea),size = 1.5)+
  geom_ribbon(data = datfit,aes(x=age.sample, ymax=ul, ymin=ll,group=Diarrhea, fill=Diarrhea), alpha=.6)+guides(fill=FALSE)+
  scale_fill_grey(start=0.5, end=0) +
  scale_x_continuous(breaks=seq(from=0,to=24,by=3),
                     labels=seq(from=0,to=24,by=3))+
  theme(legend.position = c(0.2,0.95),legend.title=element_text(size=10,face="bold"),legend.text=element_text(size=10,face="bold"))+
  theme(legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(fill=alpha('white', 0)),
        legend.key = element_rect(fill=alpha('white', 0)),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="white"),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=10,face="bold"),
        strip.text.x = element_text(size=10, face="bold"))+
  xlab("Chronological age (month)") +ylab("Microbiota age (months)") +
  facet_wrap(~month.exbf2l, ncol = 1)
pd.rm.exbf2


#' LME fit
#' <=2 months
fit<-lmer(age.predicted~diarrhea+age.sample+age.sample*age.sample +diarrhea*age.sample*age.sample +(1|personid),data=subset(rmdat.rm,month.exbf2=="<=2 months"))

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
alpha.m.rm$Diarrhea<-as.factor(alpha.m.rm$diarrhea)
alpha.m.rm$month.exbf2l<-mapvalues(alpha.m.rm$month.exbf2,from=c("<=2 months",">2 months"),to=c("Duration EBF <=2 months","Duration EBF >2 months"))
b2<-gamm(shannon~s(age.sample,by=dia.exbf2) +dia.exbf2,family=gaussian,
         data=alpha.m.rm,random=list(personid=~1))
pred <- predict(b2$gam, newdata = alpha.m.rm,se.fit=TRUE)
kable(summary(b2$gam)$p.table)
kable(summary(b2$gam)$s.table)
kable(summary(b2$lme)$tTable)
datfit<-cbind(alpha.m.rm, fit=pred$fit,ul=(pred$fit+(1.96*pred$se.fit)),ll=(pred$fit-(1.96*pred$se.fit)))
pd.s.exbf2<-ggplot()+ geom_point(data = alpha.m.rm, aes(x = age.sample, y = shannon, group = personid, colour=Diarrhea))+
  scale_color_grey(start=0.5, end=0) +
  #geom_line(data = alpha.m.rm, aes(x = age.sample, y = shannon, group = personid, colour=Diarrhea),size=0.1)+
  geom_line(data = datfit,aes(x = age.sample, y = fit, colour=Diarrhea),size = 1.5)+
  geom_ribbon(data = datfit,aes(x=age.sample, ymax=ul, ymin=ll,group=Diarrhea, fill=Diarrhea), alpha=.6)+guides(fill=FALSE)+
  scale_fill_grey(start=0.5, end=0) +
  theme(legend.position = c(0.2,0.95),legend.title=element_text(size=10,face="bold"),legend.text=element_text(size=10,face="bold"))+
  theme(legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(fill=alpha('white', 0)),
        legend.key = element_rect(fill=alpha('white', 0)),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="white"),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=10,face="bold"),
        strip.text.x = element_text(size=10, face="bold"))+
  scale_x_continuous(breaks=seq(from=0,to=24,by=3),
                     labels=seq(from=0,to=24,by=3))+
  xlab("Chronological age (month)") +ylab("Shannon index") +
  facet_wrap(~month.exbf2l, ncol = 1)
pd.s.exbf2

#combine 4 plots
#+ fig.height=5,fig.width=10
#tiff("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/fig/rm.alpha.relabund.6plus.dia.exbf2.bf.rev.tiff",width =10,height =5,units= "in",res= 300)
#pdf("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/fig/rm.alpha.relabund.6plus.dia.exbf2.bf.rev.pdf",width =10,height =5)
grid.arrange(pd.rm.exbf2,pd.s.exbf2,p.dia.exbf2.6plus.l5noleg$p,p.dia.bf.6plus.l5.nolegend$p,nrow=1)
#dev.off()

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
pd.os<-ggplot()+ geom_point(data = alpha.m.rm, aes(x = age.sample, y = observed_species, group = personid, colour=Diarrhea))+
  scale_color_grey(start=0.5, end=0) +
  #geom_line(data = alpha.m.rm, aes(x = age.sample, y = shannon, group = personid, colour=Diarrhea),size=0.1)+
  geom_line(data = datfit,aes(x = age.sample, y = fit, colour=Diarrhea),size = 1.5)+
  geom_ribbon(data = datfit,aes(x=age.sample, ymax=ul, ymin=ll,group=Diarrhea, fill=Diarrhea), alpha=.6)+guides(fill=FALSE)+
  scale_fill_grey(start=0.5, end=0) +
  theme(legend.position = c(0.2,0.95),legend.title=element_text(size=10,face="bold"),legend.text=element_text(size=10,face="bold"))+
  theme(legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(fill=alpha('white', 0)),
        legend.key = element_rect(fill=alpha('white', 0)),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="white"),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=10,face="bold"),
        strip.text.x = element_text(size=10, face="bold"))+
  scale_x_continuous(breaks=seq(from=0,to=24,by=3),
                     labels=seq(from=0,to=24,by=3))+
  #theme(legend.position = "none")+
  xlab("Chronological age (month)") +ylab("Observed_species") +facet_wrap(~month.exbf2l, ncol = 1)
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
pd.wt<-ggplot()+ geom_point(data = alpha.m.rm, aes(x = age.sample, y = pd_whole_tree, group = personid, colour=Diarrhea))+
  scale_color_grey(start=0.5, end=0) +
  #geom_line(data = alpha.m.rm, aes(x = age.sample, y = shannon, group = personid, colour=Diarrhea),size=0.1)+
  geom_line(data = datfit,aes(x = age.sample, y = fit, colour=Diarrhea),size = 1.5)+
  geom_ribbon(data = datfit,aes(x=age.sample, ymax=ul, ymin=ll,group=Diarrhea, fill=Diarrhea), alpha=.6)+guides(fill=FALSE)+
  scale_fill_grey(start=0.5, end=0) +
  theme(legend.position = c(0.2,0.95),legend.title=element_text(size=10,face="bold"),legend.text=element_text(size=10,face="bold"))+
  theme(legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(fill=alpha('white', 0)),
        legend.key = element_rect(fill=alpha('white', 0)),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="white"),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=10,face="bold"),
        strip.text.x = element_text(size=10, face="bold"))+
  scale_x_continuous(breaks=seq(from=0,to=24,by=3),
                     labels=seq(from=0,to=24,by=3))+
  #theme(legend.position = "none")+
  xlab("Chronological age (month)") +ylab("Pd_whole_tree") +facet_wrap(~month.exbf2l, ncol = 1)
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

#' ### Chao1
b2<-gamm(chao1~s(age.sample,by=dia.exbf2) +dia.exbf2,family=gaussian,
         data=alpha.m.rm,random=list(personid=~1))
pred <- predict(b2$gam, newdata = alpha.m.rm,se.fit=TRUE)
kable(summary(b2$gam)$p.table)
kable(summary(b2$gam)$s.table)
kable(summary(b2$lme)$tTable)
datfit<-cbind(alpha.m.rm, fit=pred$fit,ul=(pred$fit+(1.96*pred$se.fit)),ll=(pred$fit-(1.96*pred$se.fit)))
pd.chao<-ggplot()+ geom_point(data = alpha.m.rm, aes(x = age.sample, y = chao1, group = personid, colour=Diarrhea))+
  scale_color_grey(start=0.5, end=0) +
  #geom_line(data = alpha.m.rm, aes(x = age.sample, y = shannon, group = personid, colour=Diarrhea),size=0.1)+
  geom_line(data = datfit,aes(x = age.sample, y = fit, colour=Diarrhea),size = 1.5)+
  geom_ribbon(data = datfit,aes(x=age.sample, ymax=ul, ymin=ll,group=Diarrhea, fill=Diarrhea), alpha=.6)+guides(fill=FALSE)+
  scale_fill_grey(start=0.5, end=0) +
  theme(legend.position = c(0.2,0.95),legend.title=element_text(size=10,face="bold"),legend.text=element_text(size=10,face="bold"))+
  theme(legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(fill=alpha('white', 0)),
        legend.key = element_rect(fill=alpha('white', 0)),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="white"),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=10,face="bold"),
        strip.text.x = element_text(size=10, face="bold"))+
  scale_x_continuous(breaks=seq(from=0,to=24,by=3),
                     labels=seq(from=0,to=24,by=3))+
  #theme(legend.position = "none")+
  xlab("Chronological age (month)") +ylab("Chao1") +facet_wrap(~month.exbf2l, ncol = 1)
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

# Combined graph
# fig.with=15,fig.height=10
#grid.arrange(pd.os,pd.wt,pd.chao,ncol=3)



#' # Meta-analysis non-EBF vs. EBF adjusting for age and gender 
#' ## For comparison: adjusting for age only 
#' ### Bangladesh data: Non-EBF vs. EBF adjusting for age
#' #### Alpha diversity
load(paste(dir,"data/alphaall7s.rda",sep="")) 
rma.bf<-cbind(alphaall[alphaall$pop=="Bangladesh",1:4],index=alphaall[alphaall$pop=="Bangladesh","index"])
rma.bf$ll<-rma.bf$estimate.nebf-1.96*rma.bf$se.nebf
rma.bf$ul<-rma.bf$estimate.nebf+1.96*rma.bf$se.nebf
rma.bf[,c("estimate.nebf","ll","ul")]<-round(rma.bf[,c("estimate.nebf","ll","ul")],2)
rma.bf$pval.nebf<-round(rma.bf$pval.nebf,4)
kable(rma.bf[,c("estimate.nebf","ll","ul","pval.nebf","index")])
#' #### Microbiota age
fitsum<-summary(lmer(age.predicteds~bf+age.sample+(1|personid),data=subset(rmdat.rm,age.sample<=6)))
fitdat<-as.data.frame(fitsum$coefficients[-1,])
fitdat$ll<-fitdat$Estimate-1.96*fitdat$`Std. Error`
fitdat$ul<-fitdat$Estimate+1.96*fitdat$`Std. Error`
fitdat[,c("Estimate","ll","ul")]<-round(fitdat[,c("Estimate","ll","ul")],2)
fitdat[,"Pr(>|t|)"]<-round(fitdat[,"Pr(>|t|)"],4)
fitdat[,"pop"]<-"Bangladesh"
kable(fitdat)
#' #### Taxa relative abundance 
load(paste(dir,"data/taxacom.basic.rm.rda",sep="")) 
#' Phylum
kable(taxcomtab.show(taxcomtab=taxacom6.zi.rm,tax.select="none", tax.lev="l2",p.adjust.method="fdr"))
#' Order
kable(taxcomtab.show(taxcomtab=taxacom6.zi.rm,tax.select="none", tax.lev="l4",p.adjust.method="fdr"))
#' Family
kable(taxcomtab.show(taxcomtab=taxacom6.zi.rm,tax.select="none", tax.lev="l5",p.adjust.method="fdr"))
#' Genus
kable(taxcomtab.show(taxcomtab=taxacom6.zi.rm,tax.select="none", tax.lev="l6",p.adjust.method="fdr"))
#' #### KEGG pathway
load(paste(dir,"data/pathcom.rm6.rda",sep="")) 
#' Level 2
kable(taxcomtab.show(taxcomtab=pathcom.rm6.rel.gamlss,tax.select="none", tax.lev="l2",p.adjust.method="fdr"))
#' Level 3
kable(taxcomtab.show(taxcomtab=pathcom.rm6.rel.gamlss,tax.select="none", tax.lev="l3",p.adjust.method="fdr"))


#' ### Meta-analysis on 4 studies with gender info but adjusting for age only (no gender)
#' #### Alpha diversity
a4.bf<-cbind(alphaall[alphaall$pop %in% c("Bangladesh","Haiti","USA(CA_FL)","USA(NC)"),1:4],alphaall[alphaall$pop %in% c("Bangladesh","Haiti","USA(CA_FL)","USA(NC)"),c("index","author","pop","year")])
#' Shannon
#+ fig.width=10, fig.height=5
shannon.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(a4.bf,index=="shannon"),sm="RD", backtransf=FALSE)
forest(shannon.nebf,smlab="Standardized \n diversity difference",sortvar=subset(a4.bf,index=="shannon")$pop,lwd=2)
shannon.nebf
kable(cbind(study=shannon.nebf$studlab,pval=shannon.nebf$pval))
#' Other indexes
chao1.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(a4.bf,index=="chao1"),sm="RD", backtransf=FALSE)
observed_species.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(a4.bf,index=="observed_species"),sm="RD", backtransf=FALSE)
pd_whole_tree.nebf <- metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=subset(a4.bf,index=="pd_whole_tree"),sm="RD", backtransf=FALSE)
#show random meta-analysis model results of all indexes
atab<-as.data.frame(cbind(estimate=c(shannon.nebf$TE.random,chao1.nebf$TE.random,observed_species.nebf$TE.random,pd_whole_tree.nebf$TE.random),
                          ll=c(shannon.nebf$lower.random,chao1.nebf$lower.random,observed_species.nebf$lower.random,pd_whole_tree.nebf$lower.random),
                          ul=c(shannon.nebf$upper.random,chao1.nebf$upper.random,observed_species.nebf$upper.random,pd_whole_tree.nebf$upper.random),
                          index=c("shannon","chao1","observed_species","pd_whole_tree")))
atab[,1:3]<-lapply(atab[,1:3],as.character)
atab[,1:3]<-lapply(atab[,1:3],as.numeric)
a4.a.nebf<-ggplot(data=atab,aes(x=estimate,y=index))+
  geom_point(shape=17, colour="black",size=2)+
  geom_errorbarh(aes(xmin=ll,xmax=ul),height=0.0, colour="black",size=1)+
  geom_vline(xintercept=0,linetype="dashed")+
  scale_x_continuous(breaks=seq(from=-0.5,to=0.5,by=0.1),
                     labels=seq(from=-0.5,to=0.5,by=0.1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=10, face="bold"),
        plot.title = element_text(size=14,face="bold"))+
  ggtitle("Non-EBF vs. EBF adjusting for age")+
  xlab("Pooled standardized diversity difference")+ylab("Alpha diversity index")
a4.a.nebf


#' #### Microbiota age
#+ fig.width=10, fig.height=5
rmshare.4<-subset(as.data.frame(rmshare.sum.6s),pop %in% c("Bangladesh","Haiti","USA(CA_FL)","USA(NC)"))
rm.nebf.4<-metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=rmshare.4,sm="RD", backtransf=FALSE)
forest(rm.nebf.4,smlab="Standardized \n microbiota age difference",sortvar=rmshare.4$pop,lwd=2)
rm.nebf.4
kable(cbind(study=rm.nebf.4$studlab,pval=rm.nebf.4$pval))


#' #### Taxa relative abundance 
load(paste(dir,"data/metatab.zi.4com.rda",sep="")) 
#' Phylum (l2) 
#+ fig.width=10,fig.height=5
kable(metatab.show(metatab=metatab.zi.4com$random,taxacom.pooled.tab=taxacom.zi.4com,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.4com$random,taxacom.pooled.tab=taxacom.zi.4com,tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
metadat$taxsig<-subset(metadat$taxsig, taxa %in% c("verrucomicrobia","proteobacteria","fusobacteria","firmicutes","bacteroidetes","actinobacteria",".thermi."))
metadat$taxsig$taxa<-drop.levels(metadat$taxsig$taxa)
metadat$taxsig.all<-subset(metadat$taxsig.all, taxa %in% c("verrucomicrobia","proteobacteria","fusobacteria","firmicutes","bacteroidetes","actinobacteria",".thermi."))
metadat$taxsig.all$taxa<-drop.levels(metadat$taxsig.all$taxa)
meta.niceplot(metadat=metadat,sumtype="taxa",level="main",phyla.col="rainbow",heat.forest.width.ratio = c(1,1),point.ratio=c(3,1.5), leg.key.size=0.8,leg.text.size=10,heat.text.x.size=6,forest.axis.text.y=8,forest.axis.text.x=6)


#' Order (l4)
#+ fig.width=10, fig.height=5
kable(metatab.show(metatab=metatab.zi.4com$random,taxacom.pooled.tab=taxacom.zi.4com,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.4com$random,taxacom.pooled.tab=taxacom.zi.4com,tax.lev="l4",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",phyla.col="rainbow",heat.forest.width.ratio = c(1,1),point.ratio=c(3,1.5),leg.key.size=1,leg.text.size=10,heat.text.x.size=6,forest.axis.text.y=9,forest.axis.text.x=6)
#' Family (l5)
#+ fig.width=10, fig.height=10
kable(metatab.show(metatab=metatab.zi.4com$random,taxacom.pooled.tab=taxacom.zi.4com,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.4com$random,taxacom.pooled.tab=taxacom.zi.4com,tax.lev="l5",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",phyla.col="rainbow",heat.forest.width.ratio = c(1,1),point.ratio=c(3,1.5),leg.key.size=1,leg.text.size=10,heat.text.x.size=6,forest.axis.text.y=9,forest.axis.text.x=6)
#' Genus (l6)
#+ fig.width=10, fig.height=12
kable(metatab.show(metatab=metatab.zi.4com$random,taxacom.pooled.tab=taxacom.zi.4com,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metatab.zi.4com$random,taxacom.pooled.tab=taxacom.zi.4com,tax.lev="l6",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",phyla.col="rainbow",heat.forest.width.ratio = c(1,1),point.ratio=c(3,1.5),leg.key.size=1,leg.text.size=10,heat.text.x.size=6,forest.axis.text.y=9,forest.axis.text.x=6)

#' #### KEGG pathway
load(paste(dir,"data/pathmetatab.4.rda",sep="")) 
#' Level 2
#+ fig.width=15, fig.height=8
kable(metatab.show(metatab=pathmetatab.zi.4$random,taxacom.pooled.tab=pathcom.zi.4,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#nice plot
metadat<-metatab.show(metatab=pathmetatab.zi.4$random,taxacom.pooled.tab=pathcom.zi.4,sumvar="path",tax.lev="l2",showvar=".nebf",p.cutoff.type="p", p.cutoff=1,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",leg.key.size=1,leg.text.size=10,heat.text.x.size=8,forest.axis.text.y=7,forest.axis.text.x=8,heat.forest.width.ratio=c(1,1),point.ratio=c(4,1.5),line.ratio=c(2,0.8))

#' Level 3 KEGG pathway
#+ fig.width=15, fig.height=10
kable(metatab.show(metatab=pathmetatab.zi.4$random,taxacom.pooled.tab=pathcom.zi.4,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=pathmetatab.zi.4$random,taxacom.pooled.tab=pathcom.zi.4,sumvar="path",tax.lev="l3",showvar=".nebf",p.cutoff.type="p", p.cutoff=0.05,display="data",fill.value="log(OR)")
meta.niceplot(metadat=metadat,sumtype="path",heat.forest.width.ratio = c(1,1),forest.col="by.pvalue",leg.key.size=0.8,leg.text.size=10,heat.text.x.size=8,forest.axis.text.y=7,forest.axis.text.x=8)



#' ## Adjusting for age and gender
rm(list=ls()) # clear all
library(devtools) # for installing R package 'metamicrobiomeR' from Github
#install and load package metamicrobiomeR
install_github("nhanhocu/metamicrobiomeR")
library(metamicrobiomeR) # overwrite the functions used above 
#Load other needed packages 
library(knitr)
library(plyr)
library(dplyr)
library(gdata)
library(gridExtra)
library(ggplot2)
library(lme4) 
library(lmerTest)
library(mgcv) 
library(meta) 

#' ### Bangladesh data
#' #### Alpha diversity
data(sam.rm)
patht<-system.file("extdata/QIIME_outputs/Bangladesh/alpha_div_collated", package = "metamicrobiomeR", mustWork = TRUE)
alpha.rm<-read.multi(patht=patht,patternt=".txt",assignt="no",study="Bangladesh")
names(alpha.rm)<-sub(patht,"",names(alpha.rm))
samfile<-merge(samde, he50[,c("child.id","gender","month.exbf","month.food")],by="child.id")
samfile$age.sample<-samfile$age.months
samfile$bf<-factor(samfile$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
samfile$personid<-samfile$child.id
samfile$sampleid<-tolower(samfile$fecal.sample.id)
#comparison of standardized alpha diversity indexes between genders adjusting for breastfeeding and infant age at sample collection in infants <=6 months of age 
alphacom6.rm.bf.sexsg<-alpha.compare(datlist=alpha.rm,depth=3,mapfile=samfile,mapsampleid="fecal.sample.id",comvar="bf",adjustvar=c("age.sample","gender"),longitudinal="yes",age.limit=6,standardize=TRUE)
a.bf.a.g<-alphacom6.rm.bf.sexsg$alphasum[,1:5]
a.bf.a.g[,2:4]<-round(a.bf.a.g[,2:4],2)
a.bf.a.g[,5]<-round(a.bf.a.g[,5],4)
kable(a.bf.a.g)

#' #### Microbiota age
data(miage)
samhe<-merge(samde,he50[,c("child.id","gender","month.exbf","month.food")],by="child.id")
rmdat.rm<-merge(samhe,miage$microbiomeage.bangladesh$healthy,by.y="sampleid",by.x="fecal.sample.id")
rmdat.rm$gender<-as.factor(rmdat.rm$gender)
rmdat.rm$bf<-as.factor(rmdat.rm$bf)
rmdat.rm$personid<-paste("rm",as.factor(tolower(rmdat.rm$personid)),sep=".")
rmdat.rm$sampleid<-paste("rm",tolower(rmdat.rm$fecal.sample.id),sep=".")
rmdat.rm$author<-"Subramanian et al"
rmdat.rm$pop<-"Bangladesh"
rmdat.rm$year<-"2014"
# standardize age.predicted to have mean of zero and standard deviation of 1 
rmdat.rm$age.predicteds<-(rmdat.rm$age.predicted-mean(rmdat.rm$age.predicted,na.rm=T))/sd(rmdat.rm$age.predicted)
# Comparison in infants <=6 months of age 
fitsum<-summary(lmer(age.predicteds~bf+gender+age.sample+(1|personid),data=subset(rmdat.rm,age.sample<=6)))
fitdat<-as.data.frame(fitsum$coefficients[-1,])
fitdat[,"varname"]<-rownames(fitdat)
fitdat[,"pop"]<-"Bangladesh"
kable(fitdat)


#' #### Taxa relative abundance 
data(taxacom.rm.sex.adjustbfage)
#' phylum
kable(taxcomtab.show(taxcomtab=taxacom6.zi.rm.sex.adjustbfage,tax.select="none", showvar="Non_exclusiveBF", tax.lev="l2",p.adjust.method="fdr"))
#' order
kable(taxcomtab.show(taxcomtab=taxacom6.zi.rm.sex.adjustbfage,tax.select="none", showvar="Non_exclusiveBF", tax.lev="l4",p.adjust.method="fdr"))
#' family
kable(taxcomtab.show(taxcomtab=taxacom6.zi.rm.sex.adjustbfage,tax.select="none", showvar="Non_exclusiveBF", tax.lev="l5",p.adjust.method="fdr"))
#' genus
kable(taxcomtab.show(taxcomtab=taxacom6.zi.rm.sex.adjustbfage,tax.select="none", showvar="Non_exclusiveBF", tax.lev="l6",p.adjust.method="fdr"))

#' #### KEGG pathway
data(pathcom.rm6.rel.gamlss.sexg)
#' Level 2
kable(taxcomtab.show(taxcomtab=pathcom.rm6.rel.gamlss.sexg$l2, sumvar="path",tax.lev="l2",tax.select="none",showvar="Non_exclusiveBF", p.adjust.method="fdr",p.cutoff=0.05))
#' Level 3
kable(taxcomtab.show(taxcomtab=pathcom.rm6.rel.gamlss.sexg$l3, sumvar="path",tax.lev="l3",tax.select="none",showvar="Non_exclusiveBF", p.adjust.method="fdr",p.cutoff=0.05))


#' ### Meta-analysis of 4 studies with gender info
#' #### Alpha diversity
# load saved results of 4 studies 
data(alphacom6.sex4.scaledg)
# put data from 4 studies together for meta-analysis 
asum.ba<-alphacom6.rm.sexsg$alphasum
asum.ba$pop<-"Bangladesh"
asum.ha<-alphacom6.ha.sexsg$alphasum
asum.ha$pop<-"Haiti"
asum.cafl<-alphacom6.usbmk.sexsg$alphasum
asum.cafl$pop<-"USA(CA_FL)"
asum.unc<-alphacom6.unc.sexsg$alphasum
asum.unc$pop<-"USA(UNC)"
asum4<-rbind.fill(asum.ba,asum.ha,asum.cafl,asum.unc)
#Shannon index 
#+ fig.width=10, fig.height=5
shannon.sex <- metagen(Estimate.bfNon_exclusiveBF, `Std. Error.bfNon_exclusiveBF`, studlab=pop,data=subset(asum4,id=="shannon"),sm="RD", backtransf=FALSE)
forest(shannon.sex,smlab="Standardized \n diversity difference",sortvar=subset(asum4,id=="shannon")$pop, lwd=2)
shannon.sex
kable(cbind(study=shannon.sex$studlab,pval=shannon.sex$pval))
#' Other indexes
chao1.sex <- metagen(Estimate.bfNon_exclusiveBF, `Std. Error.bfNon_exclusiveBF`, studlab=pop,data=subset(asum4,id=="chao1"),sm="RD", backtransf=FALSE)
observed_species.sex <- metagen(Estimate.bfNon_exclusiveBF, `Std. Error.bfNon_exclusiveBF`, studlab=pop,data=subset(asum4,id=="observed_species"),sm="RD", backtransf=FALSE)
pd_whole_tree.sex <- metagen(Estimate.bfNon_exclusiveBF, `Std. Error.bfNon_exclusiveBF`, studlab=pop,data=subset(asum4,id=="pd_whole_tree"),sm="RD", backtransf=FALSE)
#show random meta-analysis model results of all indexes
atab<-as.data.frame(cbind(estimate=c(shannon.sex$TE.random,chao1.sex$TE.random,observed_species.sex$TE.random,pd_whole_tree.sex$TE.random),
                          ll=c(shannon.sex$lower.random,chao1.sex$lower.random,observed_species.sex$lower.random,pd_whole_tree.sex$lower.random),
                          ul=c(shannon.sex$upper.random,chao1.sex$upper.random,observed_species.sex$upper.random,pd_whole_tree.sex$upper.random),
                          index=c("shannon","chao1","observed_species","pd_whole_tree")))
atab[,1:3]<-lapply(atab[,1:3],as.character)
atab[,1:3]<-lapply(atab[,1:3],as.numeric)
a4.a.g<-ggplot(data=atab,aes(x=estimate,y=index))+
  geom_point(shape=17, colour="black",size=2)+
  geom_errorbarh(aes(xmin=ll,xmax=ul),height=0.0, colour="black",size=1)+
  geom_vline(xintercept=0,linetype="dashed")+
  scale_x_continuous(breaks=seq(from=-0.5,to=0.5,by=0.1),
                     labels=seq(from=-0.5,to=0.5,by=0.1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=10, face="bold"),
        plot.title = element_text(size=14,face="bold"))+
  ggtitle("Non-EBF vs. EBF adjusting for age, sex")+
  xlab("Pooled standardized diversity difference")+ylab("Alpha diversity index")
a4.a.g



#' #### Microbiota age
#+ fig.width=10, fig.height=5
data(rm4.sexs)
rmtest<-rm4.sexs
rmtest$study<-mapvalues(rmtest$pop,from=c("Bangladesh","Haiti","USA(CA_FL)","USA(NC)"),to=c("Subramanian et al 2014 (Bangladesh)","Bender et al 2016 (Haiti)","Pannaraj et al 2017 (USA(CA_FL))","Thompson et al 2015 (USA(NC))"))
rm.sex<-metagen(Estimate.bfNon_exclusiveBF, `Std. Error.bfNon_exclusiveBF`, studlab=study,data=rmtest,sm="RD", backtransf=FALSE)
forest(rm.sex,smlab="Standardized \n microbiome age difference",lwd=2)
rm.sex
kable(cbind(study=rm.sex$studlab,pval=rm.sex$pval))

#' #### Taxa relative abundance 
#load saved results 
data(metab.sex)
#' phylum
#+ fig.width=12, fig.height=10
kable(metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l2",showvar="bfNon_exclusiveBF",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#nice plot
metadat<-metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l2",showvar="bfNon_exclusiveBF",p.cutoff.type="p", p.cutoff=1,display="data")
meta.niceplot(metadat=metadat,sumtype="taxa",level="main",p="p",p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="no",
              heat.forest.width.ratio = c(1,1),forest.col="by.pvalue",point.ratio = c(4,2),line.ratio = c(2,1),
              leg.key.size=0.8,leg.text.size=10,heat.text.x.size=10,heat.text.x.angle=0,forest.axis.text.y=10,forest.axis.text.x=10)


#' order
#+ fig.width=15, fig.height=10
kable(metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l4",showvar="bfNon_exclusiveBF",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l4",showvar="bfNon_exclusiveBF",p.cutoff.type="p", p.cutoff=1,display="data")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",p="p",p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="no",
              heat.forest.width.ratio = c(1,1),forest.col="by.pvalue",leg.key.size=1,leg.text.size=8,heat.text.x.size=7,forest.axis.text.y=8,forest.axis.text.x=7,
              point.ratio = c(4,2),line.ratio = c(2,1))
#' family
#+ fig.width=15, fig.height=15
kable(metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l5",showvar="bfNon_exclusiveBF",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#nice plot
metadat<-metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l5",showvar="bfNon_exclusiveBF",p.cutoff.type="p", p.cutoff=1,display="data")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",p="p",p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="no",
              heat.forest.width.ratio = c(1,1),forest.col="by.pvalue",leg.key.size=1,leg.text.size=8,heat.text.x.size=7,forest.axis.text.y=8,forest.axis.text.x=7,
              point.ratio = c(4,2),line.ratio = c(2,1))
#' genus
#+ fig.width=15, fig.height=15
kable(metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l6",showvar="bfNon_exclusiveBF",p.cutoff.type="p", p.cutoff=0.05,display="table"))
metadat<-metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l6",showvar="bfNon_exclusiveBF",p.cutoff.type="p", p.cutoff=1,display="data")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",p="p",p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="no",
              heat.forest.width.ratio = c(1,1),forest.col="by.pvalue",leg.key.size=1,leg.text.size=8,heat.text.x.size=7,forest.axis.text.y=8,forest.axis.text.x=7,
              point.ratio = c(4,2),line.ratio = c(2,1))

#' #### KEGG pathway
data(pathmetatab.zi.sexg) 
#' Level 2
#+ fig.width=15, fig.height=15
kable(metatab.show(metatab=pathmetatab.zi.sex.l2$random,com.pooled.tab=pathcom.zi.sexg$l2,sumvar="path",showvar="bfNon_exclusiveBF",p.cutoff.type="p", p.cutoff=0.05,display="table"))
# Nice plot all pathways
metadat<-metatab.show(metatab=pathmetatab.zi.sex.l2$random,com.pooled.tab=pathcom.zi.sexg$l2,sumvar="path",showvar="bfNon_exclusiveBF",p.cutoff.type="p", p.cutoff=1,display="data")
metadat$taxsig.all$pop<-factor(metadat$taxsig.all$pop,levels=c("Bangladesh","Haiti","USA(CA_FL)","USA(NC)","Pooled"))
meta.niceplot(metadat=metadat,sumtype="path",p="p",p.adjust="p.adjust",p.sig.heat="no",
              est.break = c(-Inf, -0.5,-0.1,-0.05,0,0.05,0.1,0.5, Inf),est.break.label = c("<-0.5)", "[-0.5,-0.1)","[-0.1,-0.05)","[-0.05,0)","[0,0.05)","[0.05,0.1)","[0.1,0.5)",">=0.5"),
              heat.forest.width.ratio = c(1,1.5),forest.col="by.pvalue",leg.key.size=1,leg.text.size=8,heat.text.x.size=7,forest.axis.text.y=8,forest.axis.text.x=7,
              point.ratio = c(4,2),line.ratio = c(2,1))
#' Level 3
#+ fig.width=15, fig.height=15
kable(metatab.show(metatab=pathmetatab.zi.sex.l3$random,com.pooled.tab=pathcom.zi.sexg$l3,sumvar="path",showvar="bfNon_exclusiveBF",p.cutoff.type="p", p.cutoff=0.05,display="table"))
# Nice plot for pathways with pooled p-values<=0.05
metadat<-metatab.show(metatab=pathmetatab.zi.sex.l3$random,com.pooled.tab=pathcom.zi.sexg$l3,sumvar="path",showvar="bfNon_exclusiveBF",p.cutoff.type="p", p.cutoff=0.05,display="data")
metadat$taxsig.all$pop<-factor(metadat$taxsig.all$pop,levels=c("Bangladesh","Haiti","USA(CA_FL)","USA(NC)","Pooled"))
meta.niceplot(metadat=metadat,sumtype="path",p="p",p.adjust="p.adjust",p.sig.heat="no",
              est.break = c(-Inf, -0.5,-0.1,-0.05,0,0.05,0.1,0.5, Inf),est.break.label = c("<-0.5)", "[-0.5,-0.1)","[-0.1,-0.05)","[-0.05,0)","[0,0.05)","[0.05,0.1)","[0.1,0.5)",">=0.5"),
              heat.forest.width.ratio = c(1,1.5),forest.col="by.pvalue",leg.key.size=0.8,leg.text.size=10,heat.text.x.size=10,forest.axis.text.y=7,forest.axis.text.x=10,
              point.ratio = c(4,2),line.ratio = c(2,1))



#' ## CLR transformation: Performance of the Random Forest model to estimate microbiota age
dir<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/metamicrobiome_breastfeeding/"  #to be replaced by your directory
load(paste(dir,"data/sam.rm.rda",sep="")) 
load(paste(dir,"data/taxlist.filter.rda",sep="")) 
load(paste(dir,"data/sumstud.rda",sep=""))
load(paste(dir,"data/bangladesh.rda",sep=""))
load(paste(dir,"data/SrfFit.rml6.share.clr.rda",sep="")) 
Straining$age.predicted <- predict(SrfFit.rml6.share.clr, newdata = Straining)
actual<-Straining$Sage
predict<-Straining$age.predicted
R2 <- 1 - (sum((actual-predict )^2)/sum((actual-mean(actual))^2))
R2<-round(R2,2)
ptrain<-ggplot() +geom_point(data=Straining,aes(x=Sage, y=age.predicted))+
  theme(legend.text = element_text(colour="black", size = 10))+
  annotate("text", x=15, y=5,label=paste("R squared =",R2,sep=" "),size=5) +
  labs(title="Training set")+
  theme(legend.key.size = unit(0.5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=10,face="bold"),
        strip.text.x = element_text(size=10, face="bold"))+
  xlab("Chronological age (month)") +ylab("Microbiota age (month)")

testage <- predict(SrfFit.rml6.share.clr, newdata = Stesting)
testdat1<-cbind(sampleid=rownames(Stesting),age.sample=Stesting[,"Sage"],age.predicted=testage)
testdat2<-merge(testdat1,dclr.rml6s, by.x="sampleid",by.y="x.sampleid")
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.character)
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.numeric)
samhe<-merge(samde,he50[,c("child.id","gender","month.exbf","month.food")],by="child.id")
testdat3<-merge(samhe,testdat2,by.y="sampleid",by.x="fecal.sample.id")
rmdat.rm<-testdat3

actual<-rmdat.rm$age.sample
predict<-rmdat.rm$age.predicted
R2 <- 1 - (sum((actual-predict )^2)/sum((actual-mean(actual))^2))
R2<-round(R2,2)
ptest<-ggplot() +geom_point(data=rmdat.rm,aes(x=age.sample, y=age.predicted))+
  theme(legend.text = element_text(colour="black", size = 10))+
  annotate("text", x=15, y=5,label=paste("R squared =",R2,sep=" "), size=5) +
  labs(title="Test set")+
  theme(legend.key.size = unit(0.5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y =element_text(size=10, colour = "black",face="bold"),
        axis.text.x =element_text(size=10,face="bold",colour="black"),
        axis.title=element_text(size=10,face="bold"),
        strip.text.x = element_text(size=10, face="bold"))+
  xlab("Chronological age (month)") +ylab("Microbiota age (month)")
grid.arrange(ptrain, ptest,nrow=1)



#' R session information
sessionInfo()
