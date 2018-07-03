# R codes to generate intermediate outputs/results for the projects: 
# 1 "Effects of exclusive breastfeeding on infant gut microbiota: a meta-analysis across studies and populations"
# 2 "MetamicrobiomeR: An R package for analysis of microbiome relative abundance data using zero inflated beta GAMLSS and meta-analysis across microbiome studies using random effect models"
# 3 Some other explorations
# Author: Nhan T Ho
# Note: a simpler/revised version of the codes/workflow for the analysis and meta-analysis using our approaches are available at: https://github.com/nhanhocu/metamicrobiomeR (supplementary_file)


#rm(list=ls()) # clear all
sapply(c("phyloseq","dada2", "msa", "caret","ggplot2", "plyr", "dplyr", "reshape2","ade4", "ggrepel","GUniFrac","vegan"), require, character.only = TRUE)
#
library(ape) #read.tree
library(digest)
library(plyr)
library(gridExtra)
library(gmodels)
library(chron)
library(lubridate)
library(date)
library(xlsx)
library(data.table)
library(dplyr)
library(dtplyr)
library(tidyr)
library(Hmisc)
library(sas7bdat)
library(XLConnect)
library(knitr)
library(lme4)
library(sjmisc)
library(sjPlot)
library(lmerTest)
library(reshape2)
library(mgcv)
library(itsadug)
library(zoo)
library(geepack)
library(gdata)
# metaanalysis
library(meta)
library(metafor)
library(rmeta)
# random forest
library(rfPermute)
library(rfUtilities)
library(randomForest)

dir<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/"
source(paste(dir,"miscfun.microbiome.R",sep="")) 



# Load Subramanian's paper Bangladesh data
# get metadata
require(XLConnect)
wb <- loadWorkbook("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/nature13421-s2.xlsx")
lst = readWorksheet(wb, sheet = getSheets(wb))
he50<-lst[[1]]
colnames(he50)<-tolower(he50[1,])
colnames(he50)<-gsub(" ",".",colnames(he50))
he50<-he50[-1,]
he50<-he50[!is.na(he50$birth.cohort),]
colnames(he50)[colnames(he50) %in% "months.of.exclusive.breastfeeding"]<-"month.exbf"
colnames(he50)[colnames(he50) %in% "age.at.first.introduction.of.solid.food.(months)"]<-"month.food"
colnames(he50)[colnames(he50) %in% "age.at.first.fecal.sample.collection.(days)"]<-"day.firstsample"
colnames(he50)[colnames(he50) %in% "age.at.last.fecal.sample.collection.(days)"]<-"day.lastsample"
colnames(he50)[colnames(he50) %in% "number.of.fecal.samples.collected"]<-"n.sample"
colnames(he50)[colnames(he50) %in% "sampling.interval.(days).mean.?.sd"]<-"sampling.interval.msd"
colnames(he50)[colnames(he50) %in% "number.of.diarrhoeal.episodes./.yr"]<-"n.diarrhea.yr"
colnames(he50)[colnames(he50) %in% "%.days.with.diarrhoea.during.sampling.period"]<-"percent.time.diarrhea"
colnames(he50)[colnames(he50) %in% "fraction.of.samples.collected.where.antibiotics.had.been.consumed.within.prior.7.days"]<-"fraction.antibiotic"
colnames(he50)[colnames(he50) %in% "training-validation.set.....subject.allocation"]<-"subject.allocation"
he50[,c("month.exbf","month.food")]<-lapply(he50[,c("month.exbf","month.food")],as.character)
he50[,c("month.exbf","month.food")]<-lapply(he50[,c("month.exbf","month.food")],as.numeric)

samde<-lst[[2]]
samde[1,][!is.na(samde[2,])]<-samde[2,][!is.na(samde[2,])]
colnames(samde)<-tolower(samde[1,])
samde<-samde[-c(1,2),]
colnames(samde)[colnames(samde) %in% "medications (antibiotics and other) 4"]<- "medications"
colnames(samde)<-gsub(" ",".",colnames(samde))
colnames(samde)<-gsub(",","",colnames(samde))
colnames(samde)<-gsub("-",".",colnames(samde))
samde$bf<-NA
samde$bf[samde$breast.milk=="Yes" & samde$formula1=="No" &samde$solid.foods2=="No"]<-"ExclusiveBF"
samde$bf[samde$breast.milk=="Yes" & (samde$formula1=="Yes" |samde$solid.foods2=="Yes")]<-"Non_exclusiveBF"
samde$bf[samde$breast.milk=="No"]<-"No_BF"
samde$sam.month<-sub('.*\\.', '', samde$fecal.sample.id)
samde[,c("age.days","age.months","whz","haz","waz")]<-lapply(samde[,c("age.days","age.months","whz","haz","waz")], as.character)
samde[,c("age.days","age.months","whz","haz","waz")]<-lapply(samde[,c("age.days","age.months","whz","haz","waz")], as.numeric)
samde<-samde[!is.na(samde$fecal.sample.id),]
colnames(samde)[colnames(samde) %in% "diarrhoea.at.the.time.of.sample.collection3"]<-"diarrhea"
colnames(samde)[colnames(samde) %in% "antibiotics.within.7.days.prior.to.sample.collection"]<-"antibiotics"
samde6<-samde[samde$age.months<6,]
samde3<-samde[samde$age.months<3,]
samfile<-merge(samde, he50,by="child.id")

#first formula
fm.intro<-samde%>% group_by(child.id) %>% filter(formula1 == "Yes") %>%
  arrange(child.id,age.months) %>% mutate(fm.intro=head(age.months,1)) %>% filter(age.months==fm.intro)
samdef<-merge(samde,fm.intro[duplicated(fm.intro$child.id)==FALSE,c("child.id","fm.intro")],by=c("child.id"))
before.fm<-samdef%>% group_by(child.id) %>% filter(age.months < fm.intro) %>%
  arrange(child.id,age.months) %>% mutate(lastm.before.fm=tail(age.months,1)) %>% filter(age.months==lastm.before.fm)
#exbf.last<-samde%>% group_by(child.id) %>% filter(bf=="ExclusiveBF") %>%
#  arrange(child.id,age.months) %>% mutate(exbf.lastm=tail(age.months,1)) %>% filter(age.months==exbf.lastm)
#nofm.last<-samde%>% group_by(child.id) %>% filter(bf=="ExclusiveBF") %>%
#  arrange(child.id,age.months) %>% mutate(nofm.lastm=tail(age.months,1)) %>% filter(age.months==nofm.lastm)
ba.fm<-rbind.fill(fm.intro[(fm.intro$child.id %in% before.fm$child.id) & (fm.intro$fecal.sample.id %in% fm.intro$fecal.sample.id[duplicated(fm.intro$child.id)]==FALSE),],
                    before.fm[(before.fm$child.id %in% fm.intro$child.id) & (before.fm$fecal.sample.id %in% before.fm$fecal.sample.id[duplicated(before.fm$child.id)==FALSE]),])
ba.fm<-ba.fm %>% group_by(child.id) %>% arrange(child.id,age.months) %>%
  mutate(lastm.before.fm=mean(lastm.before.fm,na.rm=T),fm.timediff=fm.intro-lastm.before.fm, fm.intro3=cut(fm.intro, breaks=c(-Inf, 3, Inf), labels=c("<=3 months",">3months")))
ba.fm$formula<-factor(ba.fm$formula1, levels=c("No","Yes"))

#first solid
sl.intro<-samde%>% group_by(child.id) %>% filter(solid.foods2 == "Yes") %>%
  arrange(child.id,age.months) %>% mutate(sl.intro=head(age.months,1)) %>% filter(age.months==sl.intro)
samdes<-merge(samde,sl.intro[duplicated(sl.intro$child.id)==FALSE, c("child.id","sl.intro")], by="child.id")
before.sl<-samdes%>% group_by(child.id) %>% filter(age.months<sl.intro) %>%
  arrange(child.id,age.months) %>% mutate(lastm.before.sl=tail(age.months,1)) %>% filter(age.months==lastm.before.sl)
ba.sl<-rbind.fill(sl.intro[(sl.intro$child.id %in% before.sl$child.id) & (sl.intro$fecal.sample.id %in% sl.intro$fecal.sample.id[duplicated(sl.intro$child.id)]==FALSE),],
                  before.sl[(before.sl$child.id %in% sl.intro$child.id) & (before.sl$fecal.sample.id %in% before.sl$fecal.sample.id[duplicated(before.sl$child.id)==FALSE]),])
ba.sl<-ba.sl %>% group_by(child.id) %>% arrange(child.id,age.months) %>%
  mutate(lastm.before.sl=mean(lastm.before.sl,na.rm=T),sl.timediff=sl.intro-lastm.before.sl, sl.intro6=cut(sl.intro, breaks=c(-Inf, 6, Inf), labels=c("<=6 months",">6months")))
ba.sl$solid<-factor(ba.sl$solid.foods2,levels=c("No","Yes"))
#save(he50,samde,samfile,ba.fm,ba.sl,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/sam.rm.rda")
#print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/sam.rm.rda"))




# Bacterial taxa summary data
# Bangladesh (Relative maturity data)
taxsum<-"rel" # choose "rel" or "abs"
if (taxsum=="rel"){
  patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/tax_mapping7"
}
if (taxsum=="abs"){
  patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/tax_mapping_abs"
}
setwd(patht)
filenames <- list.files(path=patht, pattern=".txt")
tmp <- lapply(filenames, function(x) read.delim(file=x))
names(tmp)<-tolower(gsub(".txt","", filenames))
for (i in 1:length(names(tmp))){
  colnames(tmp[[i]])<-tolower(colnames(tmp[[i]]))
  tmp[[i]][,"study"]<-"RM_Bangladesh"
  assign(names(tmp[i]),tmp[[names(tmp[i])]])
}
dat.rm<-tmp
rm(i,tmp)
rmmap<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/Subramanian_et_al_mapping_file.txt")
colnames(rmmap)<-tolower(colnames(rmmap))
dat.rm7<-dat.rm
#save(dat.rm7,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/dat.rm7.rda")

# Haiti data
if (taxsum=="rel"){
  patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Haiti_data_010417/tax_mapping7"
}
if (taxsum=="abs"){
  patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Haiti_data_010417/tax_mapping_abs"
}
setwd(patht)
filenames <- list.files(path=patht, pattern=".txt")
tmp <- lapply(filenames, function(x) read.delim(file=x))
names(tmp)<-tolower(gsub(".txt","", filenames))
for (i in 1:length(names(tmp))){
  colnames(tmp[[i]])<-tolower(colnames(tmp[[i]]))
  tmp[[i]][,"study"]<-"UCLA_Haiti"
  assign(names(tmp[i]),tmp[[names(tmp[i])]])
}
dat.haiti<-tmp
rm(i,tmp)
haitimap<-read.delim('C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Haiti_data_010417/Haiti_Mapping.complete.txt')
colnames(haitimap)<-tolower(colnames(haitimap))
dat.haiti7<-dat.haiti

# USA(NC) Frontier data
if (taxsum=="rel"){
  patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Frontier/tax_mapping7all"
}
if (taxsum=="abs"){
  patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Frontier/tax_mapping_abs"
}
setwd(patht)
filenames <- list.files(path=patht, pattern=".txt")
tmp <- lapply(filenames, function(x) read.delim(file=x))
names(tmp)<-paste("unc",tolower(gsub(".txt","", filenames)),sep="_")
for (i in 1:length(names(tmp))){
  colnames(tmp[[i]])<-tolower(colnames(tmp[[i]]))
  tmp[[i]][,"study"]<-"UNC_US"
  assign(names(tmp[i]),tmp[[names(tmp[i])]])
}
dat.unc.all<-tmp
rm(i,tmp)
uncmap<- read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Frontier/mapping_fileall.txt")
colnames(uncmap)<-tolower(colnames(uncmap))
dat.unc<-dat.unc.all

#solid intro
sl.intro.unc<-uncmap %>% group_by(patient) %>% filter(anysolids == "1") %>%
  arrange(patient, agemo) %>% mutate(sl.intro=head(agemo,1)) %>% filter(agemo==sl.intro)
uncmaps<-merge(uncmap,sl.intro.unc[duplicated(sl.intro.unc$patient)==FALSE,c("patient","sl.intro")],by="patient")
before.sl.unc<-uncmaps%>% group_by(patient) %>% filter(agemo<sl.intro) %>%
  arrange(patient, agemo) %>% mutate(lastm.before.sl=tail(agemo,1)) %>% filter(agemo==lastm.before.sl)
ba.sl.unc<-rbind.fill(sl.intro.unc[(sl.intro.unc$patient %in% before.sl.unc$patient) & (sl.intro.unc$x.sampleid %in% sl.intro.unc$x.sampleid[duplicated(sl.intro.unc$patient)]==FALSE),],
                     before.sl.unc[(before.sl.unc$patient %in% sl.intro.unc$patient) & (before.sl.unc$x.sampleid %in% before.sl.unc$x.sampleid[duplicated(before.sl.unc$patient)==FALSE]),])
ba.sl.unc<-ba.sl.unc %>% group_by(patient) %>% arrange(patient,agemo) %>%
  mutate(lastm.before.sl=mean(lastm.before.sl,na.rm=T),sl.timediff=sl.intro-lastm.before.sl)
ba.sl.unc$solid<-mapvalues(ba.sl.unc$anysolids,from=c("0","1"),to=c("No","Yes"))
ba.sl.unc$solid<-factor(ba.sl.unc$solid,levels=c("No","Yes"))

# South Africa (UW data)
if (taxsum=="rel"){
  patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UW/tax_mapping7"
}
if (taxsum=="abs"){
  patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UW/tax_mapping_abs"
}
setwd(patht)
filenames <- list.files(path=patht, pattern=".txt")
tmp <- lapply(filenames, function(x) read.delim(file=x))
names(tmp)<-tolower(gsub(".txt","", filenames))
for (i in 1:length(names(tmp))){
  colnames(tmp[[i]])<-tolower(colnames(tmp[[i]]))
  tmp[[i]][,"study"]<-"UW_SouthAfrica"
  assign(names(tmp[i]),tmp[[names(tmp[i])]])
}
dat.uw<-tmp
rm(i,tmp)
uwmap<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UW/HJ_feeding_mapping_BPB_4Nhan.txt")
colnames(uwmap)<-tolower(colnames(uwmap))

# USA(CA_FL) (USBMK data)
if (taxsum=="rel"){
  patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/USBMK/tax_mapping7"
}
if (taxsum=="abs"){
  patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/USBMK/tax_mapping_abs"
}
setwd(patht)
filenames <- list.files(path=patht, pattern=".txt")
tmp <- lapply(filenames, function(x) read.delim(file=x))
names(tmp)<-tolower(gsub(".txt","", filenames))
for (i in 1:length(names(tmp))){
  colnames(tmp[[i]])<-tolower(colnames(tmp[[i]]))
  tmp[[i]][,"study"]<-"UCLA_USBMK"
  assign(names(tmp[i]),tmp[[names(tmp[i])]])
}
dat.usbmk<-tmp
rm(i,tmp)
usbmkmap<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/USBMK/USBMK_Mapping_122414.filtered.with_pairing_status.txt")
colnames(usbmkmap)<-tolower(colnames(usbmkmap))
#solid intro
mkbsmap<-usbmkmap[usbmkmap$mombb=="Baby"&usbmkmap$sampletype=="STL",]
mkbsmap$bbage.char<-as.character(mkbsmap$bbage)
mkbsmap$agesolidintro<-as.character(mkbsmap$agesolidintro)
sl.intro.mk<-mkbsmap %>% group_by(pairid) %>% filter(solidsintroduced == "Y") %>%
  arrange(pairid, bbage) %>% mutate(sl.intro=head(bbage,1)) %>% filter(bbage==sl.intro)
mkbsmaps<-merge(mkbsmap,sl.intro.mk[duplicated(sl.intro.mk$pairid)==FALSE,c("pairid","sl.intro")],by="pairid")
before.sl.mk<-mkbsmaps%>% group_by(pairid) %>% filter(bbage<sl.intro) %>%
  arrange(pairid, bbage) %>% mutate(lastm.before.sl=tail(bbage,1)) %>% filter(bbage==lastm.before.sl)
ba.sl.mk<-rbind.fill(sl.intro.mk[(sl.intro.mk$pairid %in% before.sl.mk$pairid) & (sl.intro.mk$x.sampleid %in% sl.intro.mk$x.sampleid[duplicated(sl.intro.mk$pairid)]==FALSE),],
                  before.sl.mk[(before.sl.mk$pairid %in% sl.intro.mk$pairid) & (before.sl.mk$x.sampleid %in% before.sl.mk$x.sampleid[duplicated(before.sl.mk$pairid)==FALSE]),])
ba.sl.mk<-ba.sl.mk %>% group_by(pairid) %>% arrange(pairid,bbage) %>%
  mutate(lastm.before.sl=mean(lastm.before.sl,na.rm=T)/30,sl.intro=sl.intro/30,sl.timediff=sl.intro-lastm.before.sl)
ba.sl.mk$solid<-mapvalues(ba.sl.mk$solidsintroduced,from=c("N","Y"),to=c("No","Yes"))
ba.sl.mk$solid<-factor(ba.sl.mk$solid,levels=c("No","Yes"))
ba.sl.mk$solid[is.na(ba.sl.mk$solid)]<-"No"
# formula intro => no reliable info to derive
#save(ba.sl.mk,usbmkmap,dat.usbmk,uwmap,dat.uw,ba.sl.unc,uncmap,dat.unc,dat.unc.old,haitimap,dat.haiti,rmmap,dat.rm,
#     file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxamapdat7.uncall.rda")


# Canada (Ualberta)
patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UAlberta"
ca<-read.multi(patht=patht,patternt=".csv",assignt="yes",study="Canada")
load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/SrfFit.rml6.shareg.train.test.rda")
shareg<-colnames(Straining)[-1]
galb<-colnames(ca$rel_abund_bjog)[grep("g__",colnames(ca$rel_abund_bjog))]
sharegall<-shareg[shareg %in% galb]
rel.alb<-ca$rel_abund_bjog

# USA(CA_MA_MO) (Harvard)
patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Harvard/ForNhan/SentToNhan"
setwd(patht)
filenames <- list.files(path=patht, pattern="_L")
tmp <- lapply(filenames, function(x) read.delim(file=x))
names(tmp)<-tolower(gsub(".txt","", filenames))
for (i in 1:length(names(tmp))){
  colnames(tmp[[i]])<-tolower(colnames(tmp[[i]]))
  tmp[[i]][,"study"]<-"Harvard (CA,MA,MO)"
  assign(names(tmp[i]),tmp[[names(tmp[i])]])
}
dat.hav<-tmp
rm(i,tmp)
harvardmap<-read.csv(file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Harvard/ForNhan/SentToNhan/final_map.csv")
colnames(harvardmap)<-tolower(colnames(harvardmap))

#7studies
#save(ba.sl.mk,usbmkmap,dat.usbmk,uwmap,dat.uw,ba.sl.unc,uncmap,dat.unc,haitimap,dat.haiti,rmmap,dat.rm,harvardmap,dat.hav,rel.alb,
#     file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxamapdat.7studies.rda")
#print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxamapdat.7studies.rda"))



# Random Forest models for microbiome age 

# analysis of RM based on Bangladesh data only
#rm.f<-taxa.filter(taxtab=dat.rm, percent.filter=0.05, relabund.filter=0.00005,max.lev=7)
#rml6<-rm.f$l6[grep("g__",rm.f$l6)]
#rml6<-rml6[grep("k__bacteria",rml6)]
#dat.rm.l6<-dat.rm$subramanian_et_al_mapping_file_l6[,c(colnames(rmmap),rml6)]
#d.rml6<-dat.rm.l6

# Concatenate  use g__ for l6 only
haitil6<-colnames(haiti_mapping.complete_l6)[grep("g__",colnames(haiti_mapping.complete_l6))]
rml6<-colnames(subramanian_et_al_mapping_file_l6)[grep("g__",colnames(subramanian_et_al_mapping_file_l6))]
uncl6<-colnames(dat.unc$unc_mapping_fileall_l6)[grep("g__",colnames(dat.unc$unc_mapping_fileall_l6))]
uwl6<-colnames(hj_feeding_mapping_bpb_4nhan_l6)[grep("g__",colnames(hj_feeding_mapping_bpb_4nhan_l6))]
usbmkl6<-colnames(usbmk_mapping_122414.filtered.with_pairing_status_l6)[grep("g__",colnames(usbmk_mapping_122414.filtered.with_pairing_status_l6))]
havardl6<-colnames(final_map_l6)[grep("g__",colnames(final_map_l6))]
canadal6<-colnames(ca$rel_abund_bjog)[grep("g__",colnames(ca$rel_abund_bjog))]
#taxa in all 7 studies
taxshare<-uwl6[(uwl6 %in% uncl6) & (uwl6 %in% rml6) & (uwl6 %in% usbmkl6) & (uwl6 %in% haitil6) &(uwl6 %in% havardl6)&(uwl6 %in% canadal6)]


#CLR transformation (for reviewer's comments)
# replacement of zero-values
library(compositions)
library(zCompositions)
dclr.rml6<-subramanian_et_al_mapping_file_l6[,c("study",tolower(colnames(rmmap)),rml6)] #
clrr<-function(x){as.data.frame(clr(x))}
test<-dclr.rml6[,taxshare]#rml6
#test0<-lrEM(test,label=0,dl=rep(1,ncol(test)),ini.cov="multRepl") # not work as no complete column (same for lrDA)
test0<-multLN(test,label=0,dl=rep(1,ncol(test)))
clrdat<-as.data.frame(clr(test0))
#clrdat<-NULL
#for (i in 1:nrow(test0)){
#  clrdat<-rbind(clrdat, clrr(test0[i,]))
#}
dclr.rml6[,taxshare]<-clrdat #rml6
dclr.rml6s<-dclr.rml6[,c("study",tolower(colnames(rmmap)),taxshare)]
SinTrain<-as.character(dclr.rml6s$x.sampleid[dclr.rml6s$ena.libraryname=="BANG_HLTHY" &dclr.rml6s$health_analysis_groups=="Healthy Singletons"])
#testing data Bangladesh study
SinTest<-as.character(dclr.rml6s$x.sampleid[dclr.rml6s$ena.libraryname!="BANG_HLTHY" & (dclr.rml6s$health_analysis_groups=="Healthy Singletons" | dclr.rml6s$health_analysis_groups=="Healthy Twins Triplets") |dclr.rml6s$health_analysis_groups=="Severe Acute Malnutrition Study"])
Sage<-dclr.rml6s$age_in_months
names(Sage)<-dclr.rml6s$x.sampleid
SdataMatrix <- cbind(Sage, dclr.rml6s[,taxshare])
Straining<-SdataMatrix[SinTrain,]
Stesting<-SdataMatrix[SinTest,]
#randomForest
library(caret)
set.seed(123)
SrfFit.rml6.share.clr <- train(Sage ~ ., data = Straining, method = "rf",preProc = "center", proximity = TRUE)
#save(dclr.rml6s,SrfFit.rml6.share.clr,Straining,Stesting,taxshare,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/SrfFit.rml6.share.clr.rda")
#predict on Bangladesh test data
testage <- predict(SrfFit.rml6.share.clr, newdata = Stesting)
testdat1<-cbind(sampleid=rownames(Stesting),age.sample=Stesting[,"Sage"],age.predicted=testage)
testdat2<-merge(testdat1,dclr.rml6s, by.x="sampleid",by.y="x.sampleid")
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.character)
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.numeric)
ggplot() +geom_point(data=testdat2,aes(x=age.sample, y=age.predicted, colour=health_analysis_groups))
samhe<-merge(samde,he50[,c("child.id","gender","month.exbf","month.food")],by="child.id")
testdat3<-merge(samhe,testdat2,by.y="sampleid",by.x="fecal.sample.id")
ggplot() +geom_point(data=testdat3,aes(x=age.sample, y=age.predicted, colour=bf))
rmdat.rm<-testdat3
rmdat.rm$bf<-factor(rmdat.rm$bf, levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
rmdat.rm$personid<-paste("rm",as.factor(tolower(rmdat.rm$personid)),sep=".")
rmdat.rm$sampleid<-paste("rm",tolower(rmdat.rm$fecal.sample.id),sep=".")
rmdat.rm$author<-"Subramanian et al"
rmdat.rm$pop<-"Bangladesh"
rmdat.rm$year<-"2014"



# fit RF model in rm data with taxshare and predict in other data
d.rml6<-subramanian_et_al_mapping_file_l6[,c("study",tolower(colnames(rmmap)),taxshare)]
d.hal6<-haiti_mapping.complete_l6[haiti_mapping.complete_l6$sampletype=="STL" &haiti_mapping.complete_l6$mombb=="B",c("study",tolower(colnames(haitimap)),taxshare)]
d.uwl6<-hj_feeding_mapping_bpb_4nhan_l6[,c("study",tolower(colnames(uwmap)),taxshare)]
d.uncl6<-unc_mapping_fileall_l6[,c("study",tolower(colnames(uncmap)),taxshare)]
d.usbmkl6<-usbmk_mapping_122414.filtered.with_pairing_status_l6[usbmk_mapping_122414.filtered.with_pairing_status_l6$sampletype=="STL" &usbmk_mapping_122414.filtered.with_pairing_status_l6$mombb=="Baby",c("study",tolower(colnames(usbmkmap)),taxshare)]
d.harvardl6<-final_map_l6[,c("x.sampleid","study","breastfeeding","age_in_days","excl_bf_til_mo", "csec","neonate_abx","labor_abx",taxshare)]
d.canadal6<-ca$rel_abund_bjog[,c("sampleid","study","bf","age.sample",taxshare)]

SinTrain<-as.character(d.rml6$x.sampleid[d.rml6$ena.libraryname=="BANG_HLTHY" &d.rml6$health_analysis_groups=="Healthy Singletons"])
#testing data in
SinTest<-as.character(d.rml6$x.sampleid[d.rml6$ena.libraryname!="BANG_HLTHY" & (d.rml6$health_analysis_groups=="Healthy Singletons" | d.rml6$health_analysis_groups=="Healthy Twins Triplets") |d.rml6$health_analysis_groups=="Severe Acute Malnutrition Study"])

Sage<-d.rml6$age_in_months
names(Sage)<-d.rml6$x.sampleid
SdataMatrix <- cbind(Sage, d.rml6[,taxshare])
Straining<-SdataMatrix[SinTrain,]
Stesting<-SdataMatrix[SinTest,]
#randomForest
set.seed(123)
SrfFit.rml6.share <- train(Sage ~ ., data = Straining, method = "rf",preProc = "center", proximity = TRUE)
# 7 studies
#save(SrfFit.rml6.share,Straining,Stesting,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/SrfFit.rml6.shareg7.train.test.rda")
#save(SrfFit.rml6.share,Straining,Stesting,taxshare,d.rml6,d.hal6,d.uwl6,d.uncl6,d.usbmkl6,d.harvardl6,d.canadal6,
#     file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/SrfFit.rml6.shareg7.train.test.dat.rda")
#print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/SrfFit.rml6.shareg7.train.test.dat.rda"))

#predict on Bangladesh data
testage <- predict(SrfFit.rml6.share, newdata = Stesting)
testdat1<-cbind(sampleid=rownames(Stesting),age.sample=Stesting[,"Sage"],age.predicted=testage)
testdat2<-merge(testdat1,d.rml6, by.x="sampleid",by.y="x.sampleid")
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.character)
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.numeric)
ggplot() +geom_point(data=testdat2,aes(x=age.sample, y=age.predicted, colour=health_analysis_groups))
samhe<-merge(samde,he50[,c("child.id","gender","month.exbf","month.food")],by="child.id")
testdat3<-merge(samhe,testdat2,by.y="sampleid",by.x="fecal.sample.id")
ggplot() +geom_point(data=testdat3,aes(x=age.sample, y=age.predicted, colour=bf))
rmdat.rm<-testdat3
rmdat.rm$bf<-factor(rmdat.rm$bf, levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
rmdat.rm$personid<-paste("rm",as.factor(tolower(rmdat.rm$personid)),sep=".")
rmdat.rm$sampleid<-paste("rm",tolower(rmdat.rm$fecal.sample.id),sep=".")
rmdat.rm$author<-"Subramanian et al"
rmdat.rm$pop<-"Bangladesh"
rmdat.rm$year<-"2014"
# standardize age.predicted
rmdat.rm$age.predicteds<-(rmdat.rm$age.predicted-mean(rmdat.rm$age.predicted,na.rm=T))/sd(rmdat.rm$age.predicted)
#rmdat.rm.l6all<-rmdat.rm
#fitsum<-summary(lmer(age.predicted~bf+age.sample+(1|personid),data=rmdat.rm))
#fitsum<-summary(lmer(age.predicted~bf+age.sample+(1|personid),data=subset(rmdat.rm,age.sample<=6)))
#standardize age.predicted
fitsum<-summary(lmer(age.predicteds~bf+age.sample+(1|personid),data=subset(rmdat.rm,age.sample<=6)))
rmshare.rm<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[4,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.rm<-matrix(rmshare.rm,nrow=1,ncol=8)
colnames(rmshare.rm)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#conbf
#rmshare.conbf.rm<-summary(lmer(age.predicted~as.numeric(bf)+age.sample+(1|personid),data=subset(rmdat.rm,age.sample<=6)))$coefficients[2,c("Estimate","Std. Error")]
rmshare.conbf.rm<-summary(lmer(age.predicteds~as.numeric(bf)+age.sample+(1|personid),data=subset(rmdat.rm,age.sample<=6)))$coefficients[2,c("Estimate","Std. Error")]

# formula intro
rmdat.fm.rm<-merge(ba.fm[,c("child.id","fecal.sample.id","fm.intro","lastm.before.fm","fm.timediff","formula")],rmdat.rm,by=c("child.id","fecal.sample.id"))
rmdat.fm.rm$personid<-as.factor(rmdat.fm.rm$personid)
ggplot() +geom_point(data=rmdat.fm.rm,aes(x=fm.intro, y=age.predicted, colour=formula))
fitsum<-summary(lmer(age.predicted~formula+fm.intro+fm.timediff+(1|personid),data=rmdat.fm.rm))
geesum<-summary(geeglm(age.predicted~formula+fm.intro+fm.timediff, data=rmdat.fm.rm, id=personid, family=gaussian))
geefull<-geeglm(age.predicted~formula+fm.intro+fm.timediff, data=rmdat.fm.rm, id=personid, family=gaussian)
geenofm<-geeglm(age.predicted~fm.intro+fm.timediff, data=rmdat.fm.rm, id=personid, family=gaussian)
anova(geefull,geenofm)

# solid intro
rmdat.sl.rm<-merge(ba.sl[,c("child.id","fecal.sample.id","sl.intro","lastm.before.sl","sl.timediff","solid")],rmdat.rm,by=c("child.id","fecal.sample.id"))
rmdat.sl.rm$personid<-as.factor(rmdat.sl.rm$personid)
ggplot() +geom_point(data=rmdat.sl.rm,aes(x=sl.intro, y=age.predicted, colour=solid))
fitsum<-summary(lmer(age.predicted~solid+sl.intro+sl.timediff+(1|personid),data=rmdat.sl.rm))
geesum<-summary(geeglm(age.predicted~solid+sl.intro+sl.timediff, data=rmdat.sl.rm, id=personid, family=gaussian))
sl.rm.rm<-geesum$coefficients[2,c(1,2)]
sl.rm.rm[,"study"]<-"Subramanian et al 2014 (Bangladesh)"
#save(SrfFit.rml6.all,Straining,Stesting,rmdat.rm,rmdat.fm.rm,rmdat.sl.rm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmdat.l6f.rda")

# gender
#fitsum<-summary(lmer(age.predicted~gender+bf+age.sample+(1|personid),data=subset(rmdat.rm,age.sample<=6)))
#standardize
fitsum<-summary(lmer(age.predicteds~gender+bf+age.sample+(1|personid),data=subset(rmdat.rm,age.sample<=6)))
fitdat<-as.data.frame(fitsum$coefficients[-1,])
fitdat[,"varname"]<-rownames(fitdat)
fitdat[,"pop"]<-"Bangladesh"
rm.ba.sex<-reshape(fitdat, idvar="pop", timevar="varname", direction="wide")


# predict on  Haiti data
test.ha<-cbind(Sage=d.hal6[,"bbage"]/30,d.hal6[,taxshare])
rownames(test.ha)<-d.hal6$x.sampleid
testage <- predict(SrfFit.rml6.share, newdata = test.ha)
testdat1<-cbind(sampleid=rownames(test.ha),age.sample=test.ha[,"Sage"],age.predicted=testage)
testdat2<-merge(testdat1,d.hal6, by.x="sampleid",by="x.sampleid")
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.character)
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.numeric)
ggplot() +geom_point(data=testdat2,aes(x=age.sample, y=age.predicted, colour=hivstatus))
fit<-glm(age.predicted~hivstatus+age.sample,data=testdat2,family=gaussian)
summary(fit)
rmdat.ha<-testdat2
rmdat.ha$bf<-mapvalues(rmdat.ha$exclusivebf,from=c("Yes","No"),to=c('ExclusiveBF',"Non_exclusiveBF"))
rmdat.ha$bf<-factor(rmdat.ha$bf,levels=c('ExclusiveBF',"Non_exclusiveBF"))
rmdat.ha$sampleid<-paste("ha",tolower(rmdat.ha$sampleid),sep=".")
rmdat.ha$author<-"Bender et al"
rmdat.ha$pop<-"Haiti"
rmdat.ha$year<-"2016"
rmdat.ha$gender<-factor(rmdat.ha$babygender, levels=c("Female","Male"))
rmdat.ha$bm<-mapvalues(rmdat.ha$delivery,from=c("C-section","Vaginal"), to=c("C_section","Vaginal"))
rmdat.ha$bm<-factor(rmdat.ha$bm,levels=c("C_section","Vaginal"))
rmdat.ha$age.predicteds<-(rmdat.ha$age.predicted-mean(rmdat.ha$age.predicted,na.rm=T))/sd(rmdat.ha$age.predicted)
#fitsum<-summary(glm(age.predicted~bf+age.sample,data=rmdat.ha))
#fitsum<-summary(glm(age.predicted~bf+age.sample,data=subset(rmdat.ha,age.sample<=6)))
#standardize
fitsum<-summary(glm(age.predicteds~bf+age.sample,data=subset(rmdat.ha,age.sample<=6)))
rmshare.ha<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.ha<-matrix(rmshare.ha,nrow=1,ncol=8)
colnames(rmshare.ha)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#stratified on birthmode
#vaginal
fitsum<-summary(glm(age.predicted~bf+age.sample,data=subset(rmdat.ha,age.sample<=6&delivery=="Vaginal")))
rmshare.ha.vag<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.ha.vag<-matrix(rmshare.ha.vag,nrow=1,ncol=8)
colnames(rmshare.ha.vag)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#c-section
fitsum<-summary(glm(age.predicted~bf+age.sample,data=subset(rmdat.ha,age.sample<=6&delivery=="C-section")))
rmshare.ha.cs<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.ha.cs<-matrix(rmshare.ha.cs,nrow=1,ncol=8)
colnames(rmshare.ha.cs)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#save(rmshare.ha.vag,rmshare.ha.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmshare.ha.vagcs.rda")

#stratify the effects of birthmode by bf
#exbf
fitsum<-summary(glm(age.predicted~delivery+age.sample,data=subset(rmdat.ha,age.sample<=6&bf=="ExclusiveBF")))
rmshare.ha.bm.exbf<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.ha.bm.exbf<-matrix(rmshare.ha.bm.exbf,nrow=1,ncol=8)
colnames(rmshare.ha.bm.exbf)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#nexbf
fitsum<-summary(glm(age.predicted~delivery+age.sample,data=subset(rmdat.ha,age.sample<=6&bf=="Non_exclusiveBF")))
rmshare.ha.bm.nexbf<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.ha.bm.nexbf<-matrix(rmshare.ha.bm.nexbf,nrow=1,ncol=8)
colnames(rmshare.ha.bm.nexbf)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#save(rmshare.ha.bm.exbf,rmshare.ha.bm.nexbf,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmshare.ha.bmbf.rda")

#birth mode adjusted for bf and age
fitsum<-summary(glm(age.predicted~bm+bf+age.sample,data=subset(rmdat.ha,age.sample<=6)))
fitdat<-as.data.frame(fitsum$coefficients[-1,])
fitdat[,"varname"]<-rownames(fitdat)
fitdat[,"pop"]<-"Haiti"
rm.ha.bm<-reshape(fitdat, idvar="pop", timevar="varname", direction="wide")

#gender adjusted for bf and age
#fitsum<-summary(glm(age.predicted~gender+bf+age.sample,data=subset(rmdat.ha,age.sample<=6)))
#standardize
fitsum<-summary(glm(age.predicteds~gender+bf+age.sample,data=subset(rmdat.ha,age.sample<=6)))
fitdat<-as.data.frame(fitsum$coefficients[-1,])
fitdat[,"varname"]<-rownames(fitdat)
fitdat[,"pop"]<-"Haiti"
rm.ha.sex<-reshape(fitdat, idvar="pop", timevar="varname", direction="wide")


# USA(CA_FL) test data usbmk
test.usbmk<-cbind(Sage=d.usbmkl6[,"bbage"]/30,d.usbmkl6[,taxshare])
rownames(test.usbmk)<-d.usbmkl6$x.sampleid
testage <- predict(SrfFit.rml6.share, newdata = test.usbmk)
testdat1<-cbind(sampleid=rownames(test.usbmk),age.sample=test.usbmk[,"Sage"],age.predicted=testage)
testdat2<-merge(testdat1,d.usbmkl6, by.x="sampleid",by="x.sampleid")
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.character)
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.numeric)
rmdat.usbmk<-testdat2
rmdat.usbmk$bf<-mapvalues(rmdat.usbmk$bfmixfm,from=c("BF","FM","Mix","SOLID"),to=c("ExclusiveBF","No_BF","Non_exclusiveBF","No_BF"))
rmdat.usbmk$bf<-factor(rmdat.usbmk$bf,levels=c('ExclusiveBF',"Non_exclusiveBF","No_BF"))
rmdat.usbmk$personid<-paste("usbmk",as.factor(tolower(rmdat.usbmk$pairid)),sep=".")
rmdat.usbmk$sampleid<-paste("usbmk",tolower(rmdat.usbmk$sampleid),sep=".")
rmdat.usbmk$author<-"Pannaraj et al"
rmdat.usbmk$pop<-"USA(CA_FL)"
rmdat.usbmk$year<-"2017"
rmdat.usbmk$bm<-mapvalues(rmdat.usbmk$delivery,from=c("C-section","Vaginal"), to=c("C_section","Vaginal"))
rmdat.usbmk$bm<-factor(rmdat.usbmk$bm,levels=c("C_section","Vaginal"))
levels(rmdat.usbmk$bm)[levels(rmdat.usbmk$bm)==""]<-NA
rmdat.usbmk$gender<-mapvalues(rmdat.usbmk$babygender,from=c("F","M"),to=c("Female","Male"))
levels(rmdat.usbmk$gender)[levels(rmdat.usbmk$gender)==""]<-NA
rmdat.usbmk$age.predicteds<-(rmdat.usbmk$age.predicted-mean(rmdat.usbmk$age.predicted,na.rm=T))/sd(rmdat.usbmk$age.predicted,na.rm=T)
ggplot() +geom_point(data=rmdat.usbmk,aes(x=age.sample, y=age.predicted, colour=bf))
#fitsum<-summary(lmer(age.predicted~bf+age.sample+(1|personid),data=rmdat.usbmk))
#fitsum<-summary(lmer(age.predicted~bf+age.sample+(1|personid),data=subset(rmdat.usbmk,age.sample<=6)))
#standardize
fitsum<-summary(lmer(age.predicteds~bf+age.sample+(1|personid),data=subset(rmdat.usbmk,age.sample<=6)))
rmshare.usbmk<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[4,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.usbmk<-matrix(rmshare.usbmk,nrow=1,ncol=8)
colnames(rmshare.usbmk)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#conbf
#rmshare.conbf.usbmk<-summary(lmer(age.predicted~as.numeric(bf)+age.sample+(1|personid),data=subset(rmdat.usbmk,age.sample<=6)))$coefficients[2,c("Estimate","Std. Error")]
#standardize
rmshare.conbf.usbmk<-summary(lmer(age.predicteds~as.numeric(bf)+age.sample+(1|personid),data=subset(rmdat.usbmk,age.sample<=6)))$coefficients[2,c("Estimate","Std. Error")]
#rmshare.conbf<-as.data.frame(rbind(rmshare.conbf.rm,rmshare.conbf.usbmk))
#colnames(rmshare.conbf)<-c("estimate.conbf","se.conbf")
#rmshare.conbf[,"study"]<-c("Subramanian et al 2014 (Bangladesh)","Pannaraj et al 2017 (USA(CA_FL))")
#save(rmshare.conbf,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmshare.conbfg.rda")

#gender adjusted for bf and age
#fitsum<-summary(lmer(age.predicted~gender+bf+age.sample+(1|personid),data=subset(rmdat.usbmk,age.sample<=6)))
#standardize
fitsum<-summary(lmer(age.predicteds~gender+bf+age.sample+(1|personid),data=subset(rmdat.usbmk,age.sample<=6)))
fitdat<-as.data.frame(fitsum$coefficients[-1,])
fitdat[,"varname"]<-rownames(fitdat)
fitdat[,"pop"]<-"USA(CA_FL)"
rm.usbmk.sex<-reshape(fitdat, idvar="pop", timevar="varname", direction="wide")

#birthmode adjusted for bf and age
fitsum<-summary(lmer(age.predicted~bm+bf+age.sample+(1|personid),data=subset(rmdat.usbmk,age.sample<=6)))
fitdat<-as.data.frame(fitsum$coefficients[-1,])
fitdat[,"varname"]<-rownames(fitdat)
fitdat[,"pop"]<-"USA(CA_FL)"
rm.usbmk.bm<-reshape(fitdat, idvar="pop", timevar="varname", direction="wide")

#stratify by birthmode
#vaginal
fitsum<-summary(lmer(age.predicted~bf+age.sample+(1|personid),data=subset(rmdat.usbmk,age.sample<=6&delivery=="Vaginal")))
rmshare.usbmk.vag<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[4,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.usbmk.vag<-matrix(rmshare.usbmk.vag,nrow=1,ncol=8)
colnames(rmshare.usbmk.vag)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#C-section
fitsum<-summary(lmer(age.predicted~bf+age.sample+(1|personid),data=subset(rmdat.usbmk,age.sample<=6&delivery=="C-section")))
rmshare.usbmk.cs<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[4,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.usbmk.cs<-matrix(rmshare.usbmk.cs,nrow=1,ncol=8)
colnames(rmshare.usbmk.cs)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#save(rmshare.usbmk.vag,rmshare.usbmk.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmshare.usbmk.vagcs.rda")

# stratify the effect of birthmode by bf
#exbf
fitsum<-summary(lmer(age.predicted~bm+age.sample+(1|personid),data=subset(rmdat.usbmk,age.sample<=6&bf=="ExclusiveBF")))
rmshare.usbmk.bm.exbf<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.usbmk.bm.exbf<-matrix(rmshare.usbmk.bm.exbf,nrow=1,ncol=8)
colnames(rmshare.usbmk.bm.exbf)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#nexbf
fitsum<-summary(lmer(age.predicted~bm+age.sample+(1|personid),data=subset(rmdat.usbmk,age.sample<=6&bf=="Non_exclusiveBF")))
rmshare.usbmk.bm.nexbf<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.usbmk.bm.nexbf<-matrix(rmshare.usbmk.bm.nexbf,nrow=1,ncol=8)
colnames(rmshare.usbmk.bm.nexbf)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#nobf
fitsum<-summary(lmer(age.predicted~bm+age.sample+(1|personid),data=subset(rmdat.usbmk,age.sample<=6&bf=="No_BF")))
rmshare.usbmk.bm.nobf<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.usbmk.bm.nobf<-matrix(rmshare.usbmk.bm.nobf,nrow=1,ncol=8)
colnames(rmshare.usbmk.bm.nobf)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#nexbf.nobf
fitsum<-summary(lmer(age.predicted~bm+age.sample+(1|personid),data=subset(rmdat.usbmk,age.sample<=6&bf %in% c("Non_exclusiveBF","No_BF"))))
rmshare.usbmk.bm.nexbf.nobf<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.usbmk.bm.nexbf.nobf<-matrix(rmshare.usbmk.bm.nexbf.nobf,nrow=1,ncol=8)
colnames(rmshare.usbmk.bm.nexbf.nobf)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#save(rmshare.usbmk.bm.exbf,rmshare.usbmk.bm.nexbf,rmshare.usbmk.bm.nobf,rmshare.usbmk.bm.nexbf.nobf,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmshare.usbmk.bmbf.rda")

# solid intro
ba.sl.mk$sampleid<-paste("usbmk",tolower(ba.sl.mk$x.sampleid),sep=".")
rmdat.sl.usbmk<-merge(ba.sl.mk[,c("pairid","sampleid","sl.intro","lastm.before.sl","sl.timediff","solid")],rmdat.usbmk,by=c("pairid","sampleid"))
rmdat.sl.usbmk$personid<-as.factor(rmdat.sl.usbmk$pairid)
ggplot() +geom_point(data=rmdat.sl.usbmk,aes(x=sl.intro, y=age.predicted, colour=solid))
fitsum<-summary(lmer(age.predicted~solid+sl.intro+sl.timediff+(1|personid),data=rmdat.sl.usbmk))
geesum<-summary(geeglm(age.predicted~solid+sl.intro+sl.timediff, data=rmdat.sl.usbmk, id=personid, family=gaussian))
sl.rm.usbmk<-geesum$coefficients[2,c(1,2)]
sl.rm.usbmk[,"study"]<-"Pannaraj et al 2017 (USA(CA_FL))"

# South Africa test data UW
d.uwl6$agemonth<-as.numeric(as.character(substr(d.uwl6$timepoint,2,3)))/4
test.uw<-cbind(Sage=d.uwl6[,"agemonth"],d.uwl6[,taxshare])
rownames(test.uw)<-d.uwl6$x.sampleid
testage <- predict(SrfFit.rml6.share, newdata = test.uw)
testdat1<-cbind(sampleid=rownames(test.uw),age.sample=test.uw[,"Sage"],age.predicted=testage)
testdat2<-merge(testdat1,d.uwl6, by.x="sampleid",by="x.sampleid")
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.character)
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.numeric)
ggplot() +geom_point(data=testdat2,aes(x=age.sample, y=age.predicted, colour=ebf))
rmdat.uw<-testdat2
rmdat.uw$bf<-mapvalues(rmdat.uw$ebf,from=c("1","0"), to=c('ExclusiveBF',"Non_exclusiveBF"))
rmdat.uw$bf<-factor(rmdat.uw$bf, levels=c('ExclusiveBF',"Non_exclusiveBF"))
rmdat.uw$personid<-paste("uw",as.factor(tolower(rmdat.uw$personid)),sep=".")
rmdat.uw$sampleid<-paste("uw",tolower(rmdat.uw$sampleid),sep=".")
rmdat.uw$author<-"Wood et al"
rmdat.uw$pop<-"South Africa"
rmdat.uw$year<-"2017"
rmdat.uw$age.predicteds<-(rmdat.uw$age.predicted-mean(rmdat.uw$age.predicted,na.rm=T))/sd(rmdat.uw$age.predicted,na.rm=T)
#fitsum<-summary(lmer(age.predicted~bf+age.sample+(1|personid),data=rmdat.uw))
#fitsum<-summary(lmer(age.predicted~bf+age.sample+(1|personid),data=subset(rmdat.uw,age.sample<=6)))
#standardize
fitsum<-summary(lmer(age.predicteds~bf+age.sample+(1|personid),data=subset(rmdat.uw,age.sample<=6)))
rmshare.uw<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.uw<-matrix(rmshare.uw,nrow=1,ncol=8)
colnames(rmshare.uw)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")


# USA(NC) test dat unc
d.uncl6$agemo<-as.numeric(as.character(d.uncl6$agemo))
test.unc<-cbind(Sage=d.uncl6$agemo,d.uncl6[,taxshare])
rownames(test.unc)<-d.uncl6$x.sampleid
testage <- predict(SrfFit.rml6.share, newdata = test.unc)
testdat1<-cbind(sampleid=rownames(test.unc),age.sample=test.unc[,"Sage"],age.predicted=testage)
testdat2<-merge(testdat1,d.uncl6, by.x="sampleid",by="x.sampleid")
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.character)
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.numeric)
ggplot() +geom_point(data=testdat2,aes(x=age.sample, y=age.predicted, colour=bf))
rmdat.unc<-testdat2
rmdat.unc$bf<-factor(rmdat.unc$bf, levels=c('ExclusiveBF',"Non_exclusiveBF","No_BF"))
rmdat.unc$personid<-paste("unc",as.factor(tolower(rmdat.unc$patient)),sep=".")
rmdat.unc$sampleid<-paste("unc",tolower(rmdat.unc$sampleid),sep=".")
rmdat.unc$author<-"Thompson et al"
rmdat.unc$pop<-"USA(NC)"
rmdat.unc$year<-"2015"
rmdat.unc$gender<-mapvalues(rmdat.unc$sex, from=c("0","1"),to=c("Male","Female"))
rmdat.unc$gender<-factor(rmdat.unc$gender,levels=c("Female","Male"))
rmdat.unc$age.predicteds<-(rmdat.unc$age.predicted-mean(rmdat.unc$age.predicted,na.rm=T))/sd(rmdat.unc$age.predicted,na.rm=T)
#fitsum<-summary(lmer(age.predicted~bf+age.sample+(1|personid),data=rmdat.unc))
#rmshare.unc<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[4,c("Estimate","Std. Error","t value","Pr(>|t|)")])
#fitsum<-summary(lmer(age.predicted~bf+age.sample+(1|personid),data=subset(rmdat.unc,age.sample<=6)))
#standardize
fitsum<-summary(lmer(age.predicteds~bf+age.sample+(1|personid),data=subset(rmdat.unc,age.sample<=6)))
rmshare.unc<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[4,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.unc<-matrix(rmshare.unc,nrow=1,ncol=8)
colnames(rmshare.unc)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#conbf
#rmshare.conbf.unc<-summary(lmer(age.predicted~as.numeric(bf)+age.sample+(1|personid),data=subset(rmdat.unc,age.sample<=6)))$coefficients[2,c("Estimate","Std. Error")]
#standardize
rmshare.conbf.unc<-summary(lmer(age.predicteds~as.numeric(bf)+age.sample+(1|personid),data=subset(rmdat.unc,age.sample<=6)))$coefficients[2,c("Estimate","Std. Error")]


#gender
#fitsum<-summary(lmer(age.predicted~gender+bf+age.sample+(1|personid),data=subset(rmdat.unc,age.sample<=6)))
#standardize
fitsum<-summary(lmer(age.predicteds~gender+bf+age.sample+(1|personid),data=subset(rmdat.unc,age.sample<=6)))
fitdat<-as.data.frame(fitsum$coefficients[-1,])
fitdat[,"varname"]<-rownames(fitdat)
fitdat[,"pop"]<-"USA(NC)"
rm.unc.sex<-reshape(fitdat, idvar="pop", timevar="varname", direction="wide")
#combine 4 studies
#rm4.sex<-rbind.fill(rm.ba.sex,rm.ha.sex,rm.usbmk.sex,rm.unc.sex)
#standardize
rm4.sexs<-rbind.fill(rm.ba.sex,rm.ha.sex,rm.usbmk.sex,rm.unc.sex)
#save(rm4.sex,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rm4.sex.rda")
#save(rm4.sexs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rm4.sexs.rda")

# solid intro
ba.sl.unc$sampleid<-paste("unc",tolower(ba.sl.unc$x.sampleid),sep=".")
rmdat.sl.unc<-merge(ba.sl.unc[,c("patient","sampleid","sl.intro","lastm.before.sl","sl.timediff","solid")],rmdat.unc,by=c("patient","sampleid"))
rmdat.sl.unc$personid<-as.factor(rmdat.sl.unc$patient)
ggplot() +geom_point(data=rmdat.sl.unc,aes(x=sl.intro, y=age.predicted, colour=solid))
fitsum<-summary(lmer(age.predicted~solid+sl.intro+sl.timediff+(1|personid),data=rmdat.sl.unc))
geesum<-summary(geeglm(age.predicted~solid+sl.intro+sl.timediff, data=rmdat.sl.unc, id=personid, family=gaussian))
sl.rm.unc<-geesum$coefficients[2,c(1,2)]
sl.rm.unc[,"study"]<-"Thompson et al 2015 (USA(NC))"


# USA(CA_MA_MO) Test Harvard
test.hav<-cbind(Sage=d.harvardl6$age_in_days/30,d.harvardl6[,taxshare])
rownames(test.hav)<-d.harvardl6$x.sampleid
testage <- predict(SrfFit.rml6.share, newdata = test.hav)
testdat1<-cbind(sampleid=rownames(test.hav),age.sample=test.hav[,"Sage"],age.predicted=testage)
testdat2<-merge(testdat1,d.harvardl6, by.x="sampleid",by="x.sampleid")
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.character)
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.numeric)
ggplot() +geom_point(data=testdat2,aes(x=age.sample, y=age.predicted, colour=breastfeeding))
rmdat.hav<-testdat2
rmdat.hav$bf<-mapvalues(rmdat.hav$breastfeeding, from=c("Excl_BF","Any_BF","No_BF"),to=c('ExclusiveBF',"Non_exclusiveBF","No_BF"))
rmdat.hav$bf<-factor(rmdat.hav$bf, levels=c('ExclusiveBF',"Non_exclusiveBF","No_BF"))
rmdat.hav$sampleid<-paste("hav",tolower(rmdat.hav$sampleid),sep=".")
rmdat.hav$author<-"Sordillo et al"
rmdat.hav$pop<-"USA(CA_MA_MO)"
rmdat.hav$year<-"2017"
rmdat.hav$bf2<-mapvalues(rmdat.hav$bf,from=c('ExclusiveBF',"Non_exclusiveBF","No_BF"),to=c('exbf',"nexbf.nobf","nexbf.nobf"))
rmdat.hav$bm<-mapvalues(rmdat.hav$csec,from=c("1","0"),to=c("C_section", "Vaginal"))
rmdat.hav$bm<-factor(rmdat.hav$bm,levels=c("C_section", "Vaginal"))
rmdat.hav$age.predicteds<-(rmdat.hav$age.predicted-mean(rmdat.hav$age.predicted,na.rm=T))/sd(rmdat.hav$age.predicted)
#fitsum<-summary(lm(age.predicted~bf+age.sample,data=subset(rmdat.hav,age.sample<=6)))
#standardize
fitsum<-summary(lm(age.predicteds~bf+age.sample,data=subset(rmdat.hav,age.sample<=6)))
rmshare.hav<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[4,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.hav<-matrix(rmshare.hav,nrow=1,ncol=8)
colnames(rmshare.hav)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#conbf
#rmshare.conbf.hav<-summary(lm(age.predicted~as.numeric(bf)+age.sample,data=subset(rmdat.hav,age.sample<=6)))$coefficients[2,c("Estimate","Std. Error")]
#standardize
rmshare.conbf.hav<-summary(lm(age.predicteds~as.numeric(bf)+age.sample,data=subset(rmdat.hav,age.sample<=6)))$coefficients[2,c("Estimate","Std. Error")]


#bm adjusted for bf and age
fitsum<-summary(lm(age.predicted~bm+bf+age.sample,data=subset(rmdat.hav,age.sample<=6)))
fitdat<-as.data.frame(fitsum$coefficients[-1,])
fitdat[,"varname"]<-rownames(fitdat)
fitdat[,"pop"]<-"USA(CA_MA_MO)"
rm.hav.bm<-reshape(fitdat, idvar="pop", timevar="varname", direction="wide")


#stratify by birthmode
#vaginal
fitsum<-summary(lm(age.predicted~bf+age.sample,data=subset(rmdat.hav,age.sample<=6&csec=="0")))
rmshare.hav.vag<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[4,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.hav.vag<-matrix(rmshare.hav.vag,nrow=1,ncol=8)
colnames(rmshare.hav.vag)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#C-section
fitsum<-summary(lm(age.predicted~bf+age.sample,data=subset(rmdat.hav,age.sample<=6&csec=="1")))
rmshare.hav.cs<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[4,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.hav.cs<-matrix(rmshare.hav.cs,nrow=1,ncol=8)
colnames(rmshare.hav.cs)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#save(rmshare.hav.vag,rmshare.hav.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmshare.hav.vagcs.rda")

# stratify the effects of birth mode by bf
#test for heterogeneity
fnoi<-lm(age.predicted~bf2+age.sample+bm,data=subset(rmdat.hav,age.sample<=6))
fi<-lm(age.predicted~bf2+age.sample+bm +bf2*bm,data=subset(rmdat.hav,age.sample<=6))
anova(fnoi,fi)
#exbf
fitsum<-summary(lm(age.predicted~bm+age.sample,data=subset(rmdat.hav,age.sample<=6&bf=="ExclusiveBF")))
rmshare.hav.bm.exbf<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.hav.bm.exbf<-matrix(rmshare.hav.bm.exbf,nrow=1,ncol=8)
colnames(rmshare.hav.bm.exbf)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#nexbf
fitsum<-summary(lm(age.predicted~bm+age.sample,data=subset(rmdat.hav,age.sample<=6&bf=="Non_exclusiveBF")))
rmshare.hav.bm.nexbf<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.hav.bm.nexbf<-matrix(rmshare.hav.bm.nexbf,nrow=1,ncol=8)
colnames(rmshare.hav.bm.nexbf)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#nobf
fitsum<-summary(lm(age.predicted~bm+age.sample,data=subset(rmdat.hav,age.sample<=6&bf=="No_BF")))
rmshare.hav.bm.nobf<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.hav.bm.nobf<-matrix(rmshare.hav.bm.nobf,nrow=1,ncol=8)
colnames(rmshare.hav.bm.nobf)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#nexbf.nobf
fitsum<-summary(lm(age.predicted~bm+age.sample,data=subset(rmdat.hav,age.sample<=6&bf %in% c("Non_exclusiveBF","No_BF"))))
rmshare.hav.bm.nexbf.nobf<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.hav.bm.nexbf.nobf<-matrix(rmshare.hav.bm.nexbf.nobf,nrow=1,ncol=8)
colnames(rmshare.hav.bm.nexbf.nobf)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#save(rmshare.hav.bm.exbf,rmshare.hav.bm.nexbf,rmshare.hav.bm.nobf,rmshare.hav.bm.nexbf.nobf,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmshare.hav.bmbf.rda")

# Test canada
test.ca<-cbind(Sage=d.canadal6$age.sample,d.canadal6[,taxshare])
rownames(test.ca)<-d.canadal6$sampleid
testage <- predict(SrfFit.rml6.share, newdata = test.ca)
testdat1<-cbind(sampleid=rownames(test.ca),age.sample=test.ca[,"Sage"],age.predicted=testage)
testdat2<-merge(testdat1[,c("sampleid","age.predicted")],d.canadal6, by="sampleid")
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.character)
testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.numeric)
ggplot() +geom_point(data=testdat2,aes(x=age.sample, y=age.predicted, colour=bf))
rmdat.ca<-testdat2
rmdat.ca$bf<-factor(rmdat.ca$bf, levels=c('ExclusiveBF',"Non_exclusiveBF","No_BF"))
rmdat.ca$sampleid<-paste("ca",tolower(rmdat.ca$sampleid),sep=".")
rmdat.ca$author<-"Azad et al"
rmdat.ca$pop<-"Canada"
rmdat.ca$year<-"2015"
rmdat.ca$age.predicteds<-(rmdat.ca$age.predicted-mean(rmdat.ca$age.predicted,na.rm=T))/sd(rmdat.ca$age.predicted,na.rm=T)
#fitsum<-summary(lm(age.predicted~bf+age.sample,data=subset(rmdat.ca,age.sample<=6)))
#standardize
fitsum<-summary(lm(age.predicteds~bf+age.sample,data=subset(rmdat.ca,age.sample<=6)))
rmshare.ca<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[4,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.ca<-matrix(rmshare.ca,nrow=1,ncol=8)
colnames(rmshare.ca)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#conbf
#rmshare.conbf.ca<-summary(lm(age.predicted~as.numeric(bf)+age.sample,data=subset(rmdat.ca,age.sample<=6)))$coefficients[2,c("Estimate","Std. Error")]
#standardize
rmshare.conbf.ca<-summary(lm(age.predicteds~as.numeric(bf)+age.sample,data=subset(rmdat.ca,age.sample<=6)))$coefficients[2,c("Estimate","Std. Error")]


#stratify by birthmode
bmca<-read.csv("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UAlberta/Birth_Mode_BJOG.csv")
colnames(bmca)<-tolower(colnames(bmca))
bmca$sampleid<-paste("ca",tolower(bmca$sampleid),sep=".")
bmca$bm<-mapvalues(bmca$birth.mode,from=c("CS-Emergency ","CS-Scheduled","Vaginal"),to=c("C_section","C_section","Vaginal"))
bmca$bm<-factor(bmca$bm,levels=c("C_section","Vaginal"))
rmdat.ca<-merge(rmdat.ca,bmca,by="sampleid")
rmdat.ca$bf2<-mapvalues(rmdat.ca$bf,from=c('ExclusiveBF',"Non_exclusiveBF","No_BF"),to=c('exbf',"nexbf.nobf","nexbf.nobf"))
#vaginal
fitsum<-summary(lm(age.predicted~bf+age.sample,data=subset(rmdat.ca,age.sample<=6&birth.mode=="Vaginal")))
rmshare.ca.vag<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[4,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.ca.vag<-matrix(rmshare.ca.vag,nrow=1,ncol=8)
colnames(rmshare.ca.vag)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#C-section
fitsum<-summary(lm(age.predicted~bf+age.sample,data=subset(rmdat.ca,age.sample<=6&birth.mode!="Vaginal")))
rmshare.ca.cs<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[4,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.ca.cs<-matrix(rmshare.ca.cs,nrow=1,ncol=8)
colnames(rmshare.ca.cs)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#test for interaction
fi<-lm(age.predicted~bf2+age.sample+bm+bm*bf2,data=subset(rmdat.ca,age.sample<=6))
fnoi<-lm(age.predicted~bf2+age.sample+bm,data=subset(rmdat.ca,age.sample<=6))
anova(fi,fnoi)
#save(rmshare.ca.vag,rmshare.ca.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmshare.ca.vagcs.rda")

#bm adjusted for bf and age
fitsum<-summary(lm(age.predicted~bm+bf+age.sample,data=subset(rmdat.ca,age.sample<=6)))
fitdat<-as.data.frame(fitsum$coefficients[-1,])
fitdat[,"varname"]<-rownames(fitdat)
fitdat[,"pop"]<-"Canada"
rm.ca.bm<-reshape(fitdat, idvar="pop", timevar="varname", direction="wide")
#combine 4 studies
rm4.bm.adjustbfage<-rbind.fill(rm.ca.bm,rm.ha.bm,rm.usbmk.bm,rm.hav.bm)
#save(rm4.bm.adjustbfage,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rm4.bm.adjustbfage.rda")

#stratify by bf for the effects of birthmode
# exbf
fitsum<-summary(lm(age.predicted~bm+age.sample,data=subset(rmdat.ca,age.sample<=6&bf=="ExclusiveBF")))
rmshare.ca.bm.exbf<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.ca.bm.exbf<-matrix(rmshare.ca.bm.exbf,nrow=1,ncol=8)
colnames(rmshare.ca.bm.exbf)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#nexbf
fitsum<-summary(lm(age.predicted~bm+age.sample,data=subset(rmdat.ca,age.sample<=6&bf=="Non_exclusiveBF")))
rmshare.ca.bm.nexbf<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.ca.bm.nexbf<-matrix(rmshare.ca.bm.nexbf,nrow=1,ncol=8)
colnames(rmshare.ca.bm.nexbf)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#nobf
fitsum<-summary(lm(age.predicted~bm+age.sample,data=subset(rmdat.ca,age.sample<=6&bf=="No_BF")))
rmshare.ca.bm.nobf<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.ca.bm.nobf<-matrix(rmshare.ca.bm.nobf,nrow=1,ncol=8)
colnames(rmshare.ca.bm.nobf)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#nexbf.nobf
fitsum<-summary(lm(age.predicted~bm+age.sample,data=subset(rmdat.ca,age.sample<=6&bf %in% c("Non_exclusiveBF","No_BF"))))
rmshare.ca.bm.nexbf.nobf<-c(fitsum$coefficients[2,c("Estimate","Std. Error","t value","Pr(>|t|)")],fitsum$coefficients[3,c("Estimate","Std. Error","t value","Pr(>|t|)")])
rmshare.ca.bm.nexbf.nobf<-matrix(rmshare.ca.bm.nexbf.nobf,nrow=1,ncol=8)
colnames(rmshare.ca.bm.nexbf.nobf)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")
#save(rmshare.ca.bm.exbf,rmshare.ca.bm.nexbf,rmshare.ca.bm.nobf,rmshare.ca.bm.nexbf.nobf,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmshare.ca.bmbf.rda")


#combine all studies for conbf
rmshare.conbf<-as.data.frame(rbind(rmshare.conbf.rm,rmshare.conbf.usbmk,rmshare.conbf.hav,rmshare.conbf.ca,rmshare.conbf.unc))
colnames(rmshare.conbf)<-c("estimate.conbf","se.conbf")
rmshare.conbf[,"study"]<-c("Subramanian et al 2014 (Bangladesh)","Pannaraj et al 2017 (USA(CA_FL))","Sordillo et al 2017 (USA(CA_MA_MO))","Azad et al 2015 (Canada)","Thompson et al 2015 (USA(NC))")
#save(rmshare.conbf,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmshare.conbfg7.rda")
#standardize
rmshare.conbfs<-rmshare.conbf
#save(rmshare.conbfs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmshare.conbfg7s.rda")
#

# combine all studies for vag and cs
rmshare.vag<- as.data.frame(rbind(rmshare.ca.vag,rmshare.ha.vag,rmshare.usbmk.vag,rmshare.hav.vag))
rmshare.vag[,"study"]<-c("Azad et al 2015 (Canada)","Bender et al 2016 (Haiti)","Pannaraj et al 2017 (USA(CA_FL))","Sordillo et al 2017 (USA(CA_MA_MO))")
rmshare.cs<- as.data.frame(rbind(rmshare.ca.cs,rmshare.ha.cs,rmshare.usbmk.cs,rmshare.hav.cs))
rmshare.cs[,"study"]<-c("Azad et al 2015 (Canada)","Bender et al 2016 (Haiti)","Pannaraj et al 2017 (USA(CA_FL))","Sordillo et al 2017 (USA(CA_MA_MO))")
#save(rmshare.vag,rmshare.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmshare.vagcs.rda")

#combine all general data
#rmshare.sum<-rbind(rmshare.unc,rmshare.uw,rmshare.usbmk,rmshare.ha,rmshare.rm)
#rownames(rmshare.sum)<-c("unc","uw","usbmk","ha","rm")
rmshare.sum.6<-as.data.frame(rbind(rmshare.unc,rmshare.uw,rmshare.usbmk,rmshare.ha,rmshare.rm,rmshare.hav,rmshare.ca))
rownames(rmshare.sum.6)<-c("unc","uw","usbmk","ha","rm","hav","ca")
rmshare.sum.6$author<-c("Thompson et al","Wood et al","Pannaraj et al","Bender et al","Subramanian et al","Sordillo et al","Azad et al")
rmshare.sum.6$pop<-c("USA(NC)","South Africa","USA(CA_FL)","Haiti","Bangladesh","USA(CA_MA_MO)","Canada")
rmshare.sum.6$year<-c("2015","2017","2017","2016","2014","2017","2015")

#save(sl.rm.unc,sl.rm.usbmk,sl.rm.rm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/sl.rm.share.rda")
#g__ only
#save(rmshare.sum.6,rmdat.unc,rmdat.uw,rmdat.usbmk,rmdat.ha,rmdat.rm,rmdat.rm.l6all,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmdat.shareg.rda")
#save(sl.rm.unc,sl.rm.usbmk,sl.rm.rm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/sl.rm.shareg.rda")
# 7 studies
#save(rmshare.sum.6,rmdat.unc,rmdat.uw,rmdat.usbmk,rmdat.ha,rmdat.rm,rmdat.hav,rmdat.ca,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmdat.shareg7.rda")
#standardize
rmshare.sum.6s<-rmshare.sum.6
#save(rmshare.sum.6s,rmdat.unc,rmdat.uw,rmdat.usbmk,rmdat.ha,rmdat.rm,rmdat.hav,rmdat.ca,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/rmdat.shareg7s.rda")




# Bacterial Taxa composition
library(glmmADMB)
library(pscl)
library(MASS)
library(boot)
library(betareg) # does not work with zero inflated proportion
library(gamlss)
library(gamlss.dist)
library(zoib)
library(ZIBR)
#Haiti
taxsum="rel" #can choose "rel" or "abs"
promed="lm" # can choose "lm" or "gamlss"
taxlev<-paste("l",2:6,sep="")
taxtab.ha<-list()
taxtab.vag.ha<-list()
taxtab.cs.ha<-list()
for (j in 1:length(dat.haiti)){
  #get baby stool samples only
  taxtab.ha[[j]]<-dat.haiti[[j]][dat.haiti[[j]]$mombb=="B"&dat.haiti[[j]]$sampletype=="STL",]
  taxtab.ha[[j]]$age.sample<-taxtab.ha[[j]]$bbage/30
  taxtab.ha[[j]]$personid<-as.factor(as.character(taxtab.ha[[j]]$pairid))
  taxtab.ha[[j]]$bf<-mapvalues(as.factor(taxtab.ha[[j]]$exclusivebf),from=c("Yes","No"),to=c("ExclusiveBF","Non_exclusiveBF"))
  taxtab.ha[[j]]$bf<-factor(taxtab.ha[[j]]$bf,levels=c("ExclusiveBF","Non_exclusiveBF"))
  taxtab.ha[[j]]$bm<-mapvalues(taxtab.ha[[j]]$delivery,from=c("C-section","Vaginal"),to=c("C_section","Vaginal"))
  taxtab.ha[[j]]$bm<-factor(taxtab.ha[[j]]$bm,levels=c("C_section","Vaginal"))
  taxtab.ha[[j]]$gender<-factor(taxtab.ha[[j]]$babygender, levels=c("Female","Male"))
  taxtab.vag.ha[[j]]<-subset(taxtab.ha[[j]],delivery=="Vaginal")
  taxtab.cs.ha[[j]]<-subset(taxtab.ha[[j]],delivery=="C-section")
}
#rel
taxacom.ha<-taxa.compare(taxtab=taxtab.ha,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxacom.zi.ha<-taxa.compare(taxtab=taxtab.ha,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxa.meansdn.ha<-taxa.meansdn(taxtab=taxtab.ha)
#save(taxacom.ha,taxacom.zi.ha,taxa.meansdn.ha,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.basic.ha.rda")
taxa.meansdn.bfbm.ha<-taxa.meansdn(taxtab=taxtab.ha,sumvar="bf",groupvar="bm")

#bm adjusted for bf and age
taxacom.ha.bm.adjustbfage<-taxa.compare.gen(taxtab=taxtab.ha[[5]],propmed.rel="lm",comvar="bm",adjustvar=c("age.sample","bf"),longitudinal="no")
taxacom.zi.ha.bm.adjustbfage<-taxa.compare.gen(taxtab=taxtab.ha[[5]],propmed.rel="gamlss",comvar="bm",adjustvar=c("age.sample","bf"),longitudinal="no")
#save(taxacom.ha.bm.adjustbfage,taxacom.zi.ha.bm.adjustbfage,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.ha.bm.adjustbfage.rda")

#sex adjusted for bf and age
#bf customized function
taxacom.ha.sex.adjustbfageo<-taxa.compare(taxtab=taxtab.ha,taxsum="rel",propmed.rel="lm",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="no",p.adjust.method="fdr")
taxacom.zi.ha.sex.adjustbfageo<-taxa.compare(taxtab=taxtab.ha,taxsum="rel",propmed.rel="gamlss",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="no",p.adjust.method="fdr")
#save(taxacom.ha.sex.adjustbfageo,taxacom.zi.ha.sex.adjustbfageo,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.ha.sex.adjustbfageo.rda")
# general function
taxacom.ha.sex.adjustbfage<-taxa.compare.gen(taxtab=taxtab.ha[[5]],propmed.rel="lm",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="no")
taxacom.zi.ha.sex.adjustbfage<-taxa.compare.gen(taxtab=taxtab.ha[[5]],propmed.rel="gamlss",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="no")
#save(taxacom.ha.sex.adjustbfage,taxacom.zi.ha.sex.adjustbfage,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.ha.sex.adjustbfage.rda")

# delivery antibiotic
taxacom.ha.m<-taxa.compare(taxtab=taxtab.ha,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar=c("age.sample","delivery","momantibiotics"),longitudinal="no",p.adjust.method="fdr")
taxacom.zi.ha.m<-taxa.compare(taxtab=taxtab.ha,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar=c("age.sample","delivery","momantibiotics"),longitudinal="no",p.adjust.method="fdr")
#save(taxacom.ha.m,taxacom.zi.ha.m,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.multi.ha.rda")

#delivery adjusted
taxacom.ha.bm<-taxa.compare(taxtab=taxtab.ha,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="no",p.adjust.method="fdr")
taxacom.zi.ha.bm<-taxa.compare(taxtab=taxtab.ha,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="no",p.adjust.method="fdr")
#save(taxacom.ha.bm,taxacom.zi.ha.bm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.bm.ha.rda")

#stratify delivery
#vaginal
taxacom.ha.vag<-taxa.compare(taxtab=taxtab.vag.ha,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr")
taxacom.zi.ha.vag<-taxa.compare(taxtab=taxtab.vag.ha,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr")
#C-section =>to be checked because of only a few samples with no comparison group
taxacom.ha.cs<-taxa.compare(taxtab=taxtab.cs.ha,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr")
taxacom.zi.ha.cs<-taxa.compare(taxtab=taxtab.cs.ha,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr")
#save(taxacom.ha.vag,taxacom.zi.ha.vag,taxacom.ha.cs,taxacom.zi.ha.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.vagcs.ha.rda")

#get taxlist after filter
taxlist.ha<-taxa.filter(taxtab=taxtab.ha,percent.filter = 0.05, relabund.filter = 0.00005)


# Bangladesh data
taxtab.rm<-list()
taxtab6.rm<-list()
taxtab6plus.rm<-list()
taxtab6plus.noexbf.rm<-list()
taxtab612.rm<-list()
taxtab6.exbf.rm<-list()
taxtab6.nexbf.rm<-list()
taxtab6.exbf2.rm<-list()
taxtab6.exbf2plus.rm<-list()
taxtab612.nexbf.rm<-list()
taxtab612.nobf.rm<-list()
taxtab6plus.nexbf.rm<-list()
taxtab6plus.nobf.rm<-list()
taxtab6plus.exbf2.rm<-list()
taxtab6plus.exbf2plus.rm<-list()
#before after fm, sl
taxtab.ba.fm.rm<-list()
taxtab.ba.sl.rm<-list()
for (j in 1:length(dat.rm)){
  taxtab.rm[[j]]<-merge(merge(dat.rm[[j]][dat.rm[[j]]$x.sampleid %in% samde$fecal.sample.id,],samde, by.x="x.sampleid",by.y="fecal.sample.id"), he50[,c("child.id","gender","month.exbf","month.food")],by.x="personid", by.y="child.id")
  taxtab.rm[[j]]$bf<-factor(taxtab.rm[[j]]$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
  taxtab.rm[[j]]$age.sample<-taxtab.rm[[j]]$age_in_months
  taxtab.rm[[j]]$gender<-factor(taxtab.rm[[j]]$gender, levels=c("Female","Male"))
  taxtab.rm[[j]]$personid<-as.factor(as.character(taxtab.rm[[j]]$personid))
  taxtab.rm[[j]]<-taxtab.rm[[j]] %>% group_by(personid) %>% arrange(personid,age.sample)  %>%
    mutate(month.food6=cut(month.food, breaks=c(-Inf, 6, Inf), labels=c("<=6 months",">6 months")),
           month.food5=cut(month.food, breaks=c(-Inf, 5, Inf), labels=c("<=5 months",">5 months")),
           month.food4=cut(month.food, breaks=c(-Inf, 4, Inf), labels=c("<=4 months",">4 months")),
           month.foodr=as.factor(as.character(round(month.food,0))),
           month.exbf3=cut(month.exbf, breaks=c(-Inf, 3, Inf), labels=c("<=3 months",">3 months")),
           month.exbf2=cut(month.exbf, breaks=c(-Inf, 2, Inf), labels=c("<=2 months",">2 months")),
           month.exbf1=cut(month.exbf, breaks=c(-Inf, 1, Inf), labels=c("<=1 months",">1 months")),
           month.exbfr=as.factor(as.character(round(month.exbf,0))))
  taxtab6.rm[[j]]<-taxtab.rm[[j]][taxtab.rm[[j]]$age.sample<=6,]
  taxtab6.exbf.rm[[j]]<-taxtab.rm[[j]][taxtab.rm[[j]]$age.sample<=6 &taxtab.rm[[j]]$bf=="ExclusiveBF",]
  taxtab6.nexbf.rm[[j]]<-taxtab.rm[[j]][taxtab.rm[[j]]$age.sample<=6 &taxtab.rm[[j]]$bf=="Non_exclusiveBF",]
  taxtab6.exbf2.rm[[j]]<-taxtab.rm[[j]][taxtab.rm[[j]]$age.sample<=6 &taxtab.rm[[j]]$month.exbf2=="<=2 months",]
  taxtab6.exbf2plus.rm[[j]]<-taxtab.rm[[j]][taxtab.rm[[j]]$age.sample<=6 &taxtab.rm[[j]]$month.exbf2==">2 months",]
  taxtab6plus.rm[[j]]<-taxtab.rm[[j]][taxtab.rm[[j]]$age.sample>6,]
  taxtab6plus.noexbf.rm[[j]]<-taxtab.rm[[j]][taxtab.rm[[j]]$age.sample>6 & taxtab.rm[[j]]$bf!="ExclusiveBF",]
  taxtab6plus.nexbf.rm[[j]]<-taxtab.rm[[j]][taxtab.rm[[j]]$age.sample>6 &taxtab.rm[[j]]$bf=="Non_exclusiveBF",]
  taxtab6plus.nobf.rm[[j]]<-taxtab.rm[[j]][taxtab.rm[[j]]$age.sample>6 &taxtab.rm[[j]]$bf=="No_BF",]
  taxtab6plus.exbf2.rm[[j]]<-taxtab.rm[[j]][taxtab.rm[[j]]$age.sample>6 &taxtab.rm[[j]]$month.exbf2=="<=2 months",]
  taxtab6plus.exbf2plus.rm[[j]]<-taxtab.rm[[j]][taxtab.rm[[j]]$age.sample>6 &taxtab.rm[[j]]$month.exbf2==">2 months",]
  taxtab612.rm[[j]]<-taxtab.rm[[j]][taxtab.rm[[j]]$age.sample>6 &taxtab.rm[[j]]$age.sample<=12,]
  taxtab612.nexbf.rm[[j]]<-taxtab.rm[[j]][taxtab.rm[[j]]$age.sample>6 &taxtab.rm[[j]]$age.sample<=12 &taxtab.rm[[j]]$bf=="Non_exclusiveBF",]
  taxtab612.nobf.rm[[j]]<-taxtab.rm[[j]][taxtab.rm[[j]]$age.sample>6 &taxtab.rm[[j]]$age.sample<=12 &taxtab.rm[[j]]$bf=="No_BF",]
  taxtab.ba.fm.rm[[j]]<-merge(ba.fm[,c("fecal.sample.id","fm.intro","lastm.before.fm","fm.timediff","formula")],taxtab.rm[[j]], by.x="fecal.sample.id",by.y="x.sampleid")
  taxtab.ba.sl.rm[[j]]<-merge(ba.sl[,c("fecal.sample.id","sl.intro","lastm.before.sl","sl.timediff","solid")],taxtab.rm[[j]], by.x="fecal.sample.id",by.y="x.sampleid")
}
#save(taxtab.rm,taxtab6.rm,taxtab6plus.rm,taxtab6plus.noexbf.rm,taxtab612.rm,
#     taxtab6.exbf.rm,taxtab6.nexbf.rm,taxtab6.exbf2.rm,taxtab6.exbf2plus.rm,
#     taxtab612.nexbf.rm,taxtab612.nobf.rm,taxtab6plus.nexbf.rm,taxtab6plus.nobf.rm,
#     taxtab6plus.exbf2.rm,taxtab6plus.exbf2plus.rm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxtab.rm7.rda")

#get taxlist after filter
taxlist.rm<-taxa.filter(taxtab=taxtab.rm,percent.filter = 0.05, relabund.filter = 0.00005)

taxacom.rm<-taxa.compare(taxtab=taxtab.rm,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom6.rm<-taxa.compare(taxtab=taxtab6.rm,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.zi.rm<-taxa.compare(taxtab=taxtab.rm,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom6.zi.rm<-taxa.compare(taxtab=taxtab6.rm,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxa.meansdn.rm<-taxa.meansdn(taxtab=taxtab.rm)
#save(taxacom.rm,taxacom6.rm,taxacom.zi.rm,taxacom6.zi.rm,taxa.meansdn.rm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.basic.rm.rda")
#generic function
taxacom6.rmg<-taxa.compare.gen(taxtab=taxtab6.rm[[5]],taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom6.zi.rmg<-taxa.compare.gen(taxtab=taxtab6.rm[[5]],taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
#save(taxacom6.rmg,taxacom6.zi.rmg,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom6.rmg.rda")

#gender adjusted for bf and age
#bf customized function
taxtab6plus.rm.sex.adjustbfageo<-taxa.compare(taxtab=taxtab6plus.rm,propmed.rel="lm",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="yes",p.adjust.method="fdr")
taxtab6plus.zi.rm.sex.adjustbfageo<-taxa.compare(taxtab=taxtab6plus.rm,propmed.rel="gamlss",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="yes",p.adjust.method="fdr")
taxacom6.rm.sex.adjustbfageo<-taxa.compare(taxtab=taxtab6.rm,propmed.rel="lm",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="yes",p.adjust.method="fdr")
taxacom6.zi.rm.sex.adjustbfageo<-taxa.compare(taxtab=taxtab6.rm,propmed.rel="gamlss",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="yes",p.adjust.method="fdr")
#save(taxtab6plus.rm.sex.adjustbfageo,taxtab6plus.zi.rm.sex.adjustbfageo,taxacom6.rm.sex.adjustbfageo,taxacom6.zi.rm.sex.adjustbfageo,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.rm.sex.adjustbfageo.rda")
#generic function
taxtab6plus.rm.sex.adjustbfage<-taxa.compare.gen(taxtab=taxtab6plus.rm[[5]],propmed.rel="lm",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="yes")
taxtab6plus.zi.rm.sex.adjustbfage<-taxa.compare.gen(taxtab=taxtab6plus.rm[[5]],propmed.rel="gamlss",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="yes")
taxacom6.rm.sex.adjustbfage<-taxa.compare.gen(taxtab=taxtab6.rm[[5]],propmed.rel="lm",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="yes")
taxacom6.zi.rm.sex.adjustbfage<-taxa.compare.gen(taxtab=taxtab6.rm[[5]],propmed.rel="gamlss",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="yes")
taxa.meansdn.rm.sexage<-taxa.meansdn(taxtab=taxtab.rm,sumvar="gender",groupvar="age.sample")
#save(taxtab6plus.rm.sex.adjustbfage,taxtab6plus.zi.rm.sex.adjustbfage,taxacom6.rm.sex.adjustbfage,taxacom6.zi.rm.sex.adjustbfage,taxa.meansdn.rm.sexage,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.rm.sex.adjustbfage.rda")

# compare afer 6 month and 6-12 month for month.food5
taxacom.6plus.sl5.rm<-taxa.compare(taxtab=taxtab6plus.rm,taxsum="rel",propmed.rel="lm",comvar="month.food5",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.612.sl5.rm<-taxa.compare(taxtab=taxtab612.rm,taxsum="rel",propmed.rel="lm",comvar="month.food5",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6plus.sl5.zi.rm<-taxa.compare(taxtab=taxtab6plus.rm,taxsum="rel",propmed.rel="gamlss",comvar="month.food5",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.612.sl5.zi.rm<-taxa.compare(taxtab=taxtab612.rm,taxsum="rel",propmed.rel="gamlss",comvar="month.food5",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxa.meansdn.sl5.rm<-taxa.meansdn(taxtab=taxtab.rm,sumvar="month.food5")
taxa.meansdn.food.rm<-taxa.meansdn(taxtab=taxtab.rm,sumvar="month.foodr")
#generic function
taxacom.6plus.sl5.rmg<-taxa.compare.gen(taxtab=taxtab6plus.rm[[5]],taxsum="rel",propmed.rel="lm",comvar="month.food5",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6plus.sl5.zi.rmg<-taxa.compare.gen(taxtab=taxtab6plus.rm[[5]],taxsum="rel",propmed.rel="gamlss",comvar="month.food5",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
#save(taxacom.6plus.sl5.rmg,taxacom.6plus.sl5.zi.rmg,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.6plus.sl5.rmg.rda")

# compare after 6 month for month.exbf2
taxacom.6plus.exbf2.rm<-taxa.compare(taxtab=taxtab6plus.rm,taxsum="rel",propmed.rel="lm",comvar="month.exbf2",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.612.exbf2.rm<-taxa.compare(taxtab=taxtab612.rm,taxsum="rel",propmed.rel="lm",comvar="month.exbf2",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6plus.exbf2.zi.rm<-taxa.compare(taxtab=taxtab6plus.rm,taxsum="rel",propmed.rel="gamlss",comvar="month.exbf2",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.612.exbf2.zi.rm<-taxa.compare(taxtab=taxtab612.rm,taxsum="rel",propmed.rel="gamlss",comvar="month.exbf2",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxa.meansdn.exbf2.rm<-taxa.meansdn(taxtab=taxtab.rm,sumvar="month.exbf2")
taxa.meansdn.exbf3.rm<-taxa.meansdn(taxtab=taxtab.rm,sumvar="month.exbf3")
taxa.meansdn.exbf.rm<-taxa.meansdn(taxtab=taxtab.rm,sumvar="month.exbfr")
# compare after 6 month for no-BF vs. non-exbf
taxacom.6plus.noexbf.rm<-taxa.compare(taxtab=taxtab6plus.noexbf.rm,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6plus.noexbf.zi.rm<-taxa.compare(taxtab=taxtab6plus.noexbf.rm,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")

#save(taxacom.6plus.sl5.rm,taxacom.612.sl5.rm,taxacom.6plus.sl5.zi.rm,taxacom.612.sl5.zi.rm,taxa.meansdn.sl5.rm,
#     taxacom.6plus.exbf2.rm,taxacom.612.exbf2.rm,taxacom.6plus.exbf2.zi.rm,taxacom.612.exbf2.zi.rm,taxa.meansdn.exbf2.rm,taxa.meansdn.exbf3.rm,
#     file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.612plus.food5.exbf2f.rda")

#save(taxacom.6plus.noexbf.rm,taxacom.6plus.noexbf.zi.rm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.6plus.noexbf.rmf.rda")

#compare diarrhea, antibiotics
taxacom.6.dia.rm<-taxa.compare(taxtab=taxtab6.rm,taxsum="rel",propmed.rel="lm",comvar="diarrhea",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6.dia.exbf.rm<-taxa.compare(taxtab=taxtab6.exbf.rm,taxsum="rel",propmed.rel="lm",comvar="diarrhea",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6.dia.nexbf.rm<-taxa.compare(taxtab=taxtab6.nexbf.rm,taxsum="rel",propmed.rel="lm",comvar="diarrhea",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6.dia.zi.rm<-taxa.compare(taxtab=taxtab6.rm,taxsum="rel",propmed.rel="gamlss",comvar="diarrhea",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
#error with gamlss => no longitudinal
taxacom.6.dia.exbf.zi.rm<-taxa.compare(taxtab=taxtab6.exbf.rm,taxsum="rel",propmed.rel="gamlss",comvar="diarrhea",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
# error with gamlss => no longitudinal
taxacom.6.dia.nexbf.zi.rm<-taxa.compare(taxtab=taxtab6.nexbf.rm,taxsum="rel",propmed.rel="gamlss",comvar="diarrhea",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxacom.6plus.dia.rm<-taxa.compare(taxtab=taxtab6plus.rm,taxsum="rel",propmed.rel="lm",comvar="diarrhea",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6plus.dia.nobf.rm<-taxa.compare(taxtab=taxtab6plus.nobf.rm,taxsum="rel",propmed.rel="lm",comvar="diarrhea",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6plus.dia.nexbf.rm<-taxa.compare(taxtab=taxtab6plus.nexbf.rm,taxsum="rel",propmed.rel="lm",comvar="diarrhea",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6plus.dia.zi.rm<-taxa.compare(taxtab=taxtab6plus.rm,taxsum="rel",propmed.rel="gamlss",comvar="diarrhea",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
# error with gamlss => no longitudinal
taxacom.6plus.dia.nobf.zi.rm<-taxa.compare(taxtab=taxtab6plus.nobf.rm,taxsum="rel",propmed.rel="gamlss",comvar="diarrhea",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
# error with gamlss => no longitudinal
taxacom.6plus.dia.nexbf.zi.rm<-taxa.compare(taxtab=taxtab6plus.nexbf.rm,taxsum="rel",propmed.rel="gamlss",comvar="diarrhea",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
# save results with gamlss no longitudinal
#save(taxacom.6.dia.exbf.zi.rm,taxacom.6.dia.nexbf.zi.rm,taxacom.6plus.dia.nobf.zi.rm,taxacom.6plus.dia.nexbf.zi.rm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.dia.bf.zinolong.rda")
#generic function
taxacom.6plus.dia.nobf.zi.rmg<-taxa.compare.gen(taxtab=taxtab6plus.nobf.rm[[5]],taxsum="rel",propmed.rel="gamlss",comvar="diarrhea",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxacom.6plus.dia.nobf.rmg<-taxa.compare.gen(taxtab=taxtab6plus.nobf.rm[[5]],taxsum="rel",propmed.rel="lm",comvar="diarrhea",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxacom.6plus.dia.nexbf.zi.rmg<-taxa.compare.gen(taxtab=taxtab6plus.nexbf.rm[[5]],taxsum="rel",propmed.rel="gamlss",comvar="diarrhea",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxacom.6plus.dia.nexbf.rmg<-taxa.compare.gen(taxtab=taxtab6plus.nexbf.rm[[5]],taxsum="rel",propmed.rel="lm",comvar="diarrhea",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
#save(taxacom.6plus.dia.nobf.zi.rmg,taxacom.6plus.dia.nobf.rmg,taxacom.6plus.dia.nexbf.zi.rmg,taxacom.6plus.dia.nexbf.rmg,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.6plus.dia.bf.zi.rmg.rda")

taxa.meansdn.dia.rm<-taxa.meansdn(taxtab=taxtab.rm,sumvar="diarrhea")
taxa.meansdn.dia.exbf2.6.rm<-taxa.meansdn(taxtab=taxtab6.rm,sumvar="diarrhea", groupvar="month.exbf2")
taxa.meansdn.dia.bf.6.rm<-taxa.meansdn(taxtab=taxtab6.rm,sumvar="diarrhea", groupvar="bf")
taxa.meansdn.dia.exbf2.6plus.rm<-taxa.meansdn(taxtab=taxtab6plus.rm,sumvar="diarrhea", groupvar="month.exbf2")
taxa.meansdn.dia.bf.6plus.rm<-taxa.meansdn(taxtab=taxtab6plus.rm,sumvar="diarrhea", groupvar="bf")
taxa.meansdn.ab.rm<-taxa.meansdn(taxtab=taxtab.rm,sumvar="antibiotics")
taxa.meansdn.ab.exbf2.6.rm<-taxa.meansdn(taxtab=taxtab6.rm,sumvar="antibiotics", groupvar="month.exbf2")
taxa.meansdn.ab.bf.6.rm<-taxa.meansdn(taxtab=taxtab6.rm,sumvar="antibiotics", groupvar="bf")
taxa.meansdn.ab.exbf2.6plus.rm<-taxa.meansdn(taxtab=taxtab6plus.rm,sumvar="antibiotics", groupvar="month.exbf2")
taxa.meansdn.ab.bf.6plus.rm<-taxa.meansdn(taxtab=taxtab6plus.rm,sumvar="antibiotics", groupvar="bf")

#save(taxacom.6.dia.rm,taxacom.6.dia.exbf.rm,taxacom.6.dia.nexbf.rm,taxacom.6.dia.zi.rm,
#     taxacom.6plus.dia.rm,taxacom.6plus.dia.nobf.rm,taxacom.6plus.dia.nexbf.rm,taxacom.6plus.dia.zi.rm,
#     taxa.meansdn.dia.rm,taxa.meansdn.dia.exbf2.6.rm,taxa.meansdn.dia.bf.6.rm,taxa.meansdn.dia.exbf2.6plus.rm,taxa.meansdn.dia.bf.6plus.rm,
#     taxa.meansdn.ab.rm,taxa.meansdn.ab.exbf2.6.rm,taxa.meansdn.ab.bf.6.rm,taxa.meansdn.ab.exbf2.6plus.rm,taxa.meansdn.ab.bf.6plus.rm,
#     file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.dia.ab.rda")
save(taxacom.6.dia.rm,taxacom.6.dia.exbf.rm,taxacom.6.dia.nexbf.rm,taxacom.6.dia.zi.rm,
     taxacom.6plus.dia.rm,taxacom.6plus.dia.nobf.rm,taxacom.6plus.dia.nexbf.rm,taxacom.6plus.dia.zi.rm,
     taxa.meansdn.dia.rm,taxa.meansdn.dia.exbf2.6.rm,taxa.meansdn.dia.bf.6.rm,taxa.meansdn.dia.exbf2.6plus.rm,taxa.meansdn.dia.bf.6plus.rm,
     taxa.meansdn.ab.rm,taxa.meansdn.ab.exbf2.6.rm,taxa.meansdn.ab.bf.6.rm,taxa.meansdn.ab.exbf2.6plus.rm,taxa.meansdn.ab.bf.6plus.rm,
     file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.dia.abf.rda")

# dia, antibiotic for exbf2
#diarrhea
taxacom.6.dia.exbf2.rm<-taxa.compare(taxtab=taxtab6.exbf2.rm,taxsum="rel",propmed.rel="lm",comvar="diarrhea",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6.dia.exbf2plus.rm<-taxa.compare(taxtab=taxtab6.exbf2plus.rm,taxsum="rel",propmed.rel="lm",comvar="diarrhea",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6plus.dia.exbf2.rm<-taxa.compare(taxtab=taxtab6plus.exbf2.rm,taxsum="rel",propmed.rel="lm",comvar="diarrhea",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6plus.dia.exbf2plus.rm<-taxa.compare(taxtab=taxtab6plus.exbf2plus.rm,taxsum="rel",propmed.rel="lm",comvar="diarrhea",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
#save(taxacom.6.dia.exbf2.rm,taxacom.6.dia.exbf2plus.rm,taxacom.6plus.dia.exbf2.rm,taxacom.6plus.dia.exbf2plus.rm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.dia.exbf2.rda")
# dirty fix error gamlss no longitudinal
taxacom.6.dia.exbf2.zi.rm<-taxa.compare(taxtab=taxtab6.exbf2.rm,taxsum="rel",propmed.rel="gamlss",comvar="diarrhea",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxacom.6.dia.exbf2plus.zi.rm<-taxa.compare(taxtab=taxtab6.exbf2plus.rm,taxsum="rel",propmed.rel="gamlss",comvar="diarrhea",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxacom.6plus.dia.exbf2.zi.rm<-taxa.compare(taxtab=taxtab6plus.exbf2.rm,taxsum="rel",propmed.rel="gamlss",comvar="diarrhea",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxacom.6plus.dia.exbf2plus.zi.rm<-taxa.compare(taxtab=taxtab6plus.exbf2plus.rm,taxsum="rel",propmed.rel="gamlss",comvar="diarrhea",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
#save(taxacom.6.dia.exbf2.zi.rm,taxacom.6.dia.exbf2plus.zi.rm,taxacom.6plus.dia.exbf2.zi.rm,taxacom.6plus.dia.exbf2plus.zi.rm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.dia.exbf2.zinolong.rda")
#generic function
taxacom.6plus.dia.exbf2.rmg<-taxa.compare.gen(taxtab=taxtab6plus.exbf2.rm[[5]],taxsum="rel",propmed.rel="lm",comvar="diarrhea",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6plus.dia.exbf2plus.rmg<-taxa.compare.gen(taxtab=taxtab6plus.exbf2plus.rm[[5]],taxsum="rel",propmed.rel="lm",comvar="diarrhea",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6plus.dia.exbf2.zi.rmg<-taxa.compare.gen(taxtab=taxtab6plus.exbf2.rm[[5]],taxsum="rel",propmed.rel="gamlss",comvar="diarrhea",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxacom.6plus.dia.exbf2plus.zi.rmg<-taxa.compare.gen(taxtab=taxtab6plus.exbf2plus.rm[[5]],taxsum="rel",propmed.rel="gamlss",comvar="diarrhea",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
#save(taxacom.6plus.dia.exbf2.rmg,taxacom.6plus.dia.exbf2plus.rmg,taxacom.6plus.dia.exbf2.zi.rmg,taxacom.6plus.dia.exbf2plus.zi.rmg,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.dia.exbf2.zi.rmg.rda")

#antibiotic
taxacom.6.ab.exbf2.rm<-taxa.compare(taxtab=taxtab6.exbf2.rm,taxsum="rel",propmed.rel="lm",comvar="antibiotics",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6.ab.exbf2plus.rm<-taxa.compare(taxtab=taxtab6.exbf2plus.rm,taxsum="rel",propmed.rel="lm",comvar="antibiotics",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6plus.ab.exbf2.rm<-taxa.compare(taxtab=taxtab6plus.exbf2.rm,taxsum="rel",propmed.rel="lm",comvar="antibiotics",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.6plus.ab.exbf2plus.rm<-taxa.compare(taxtab=taxtab6plus.exbf2plus.rm,taxsum="rel",propmed.rel="lm",comvar="antibiotics",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
#save(taxacom.6.ab.exbf2.rm,taxacom.6.ab.exbf2plus.rm,taxacom.6plus.ab.exbf2.rm,taxacom.6plus.ab.exbf2plus.rm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.ab.exbf2.rda")
# dirty fix error gamlss no longitudinal
taxacom.6.ab.exbf2.zi.rm<-taxa.compare(taxtab=taxtab6.exbf2.rm,taxsum="rel",propmed.rel="gamlss",comvar="antibiotics",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxacom.6.ab.exbf2plus.zi.rm<-taxa.compare(taxtab=taxtab6.exbf2plus.rm,taxsum="rel",propmed.rel="gamlss",comvar="antibiotics",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxacom.6plus.ab.exbf2.zi.rm<-taxa.compare(taxtab=taxtab6plus.exbf2.rm,taxsum="rel",propmed.rel="gamlss",comvar="antibiotics",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxacom.6plus.ab.exbf2plus.zi.rm<-taxa.compare(taxtab=taxtab6plus.exbf2plus.rm,taxsum="rel",propmed.rel="gamlss",comvar="antibiotics",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
#save(taxacom.6.ab.exbf2.zi.rm,taxacom.6.ab.exbf2plus.zi.rm,taxacom.6plus.ab.exbf2.zi.rm,taxacom.6plus.ab.exbf2plus.zi.rm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.ab.exbf2.zinolong.rda")


#after vs. before solid intro
taxacom.ba.sl.rm<-taxa.compare(taxtab=taxtab.ba.sl.rm,taxsum="rel",propmed.rel="lm",comvar="solid",adjustvar="sl.intro",longitudinal="yes",p.adjust.method="fdr")
taxacom.ba.sl.zi.rm<-taxa.compare(taxtab=taxtab.ba.sl.rm,taxsum="rel",propmed.rel="gamlss",comvar="solid",adjustvar="sl.intro",longitudinal="yes",p.adjust.method="fdr")
#save(taxacom.ba.sl.usbmk,taxacom.ba.sl.rm,taxacom.ba.sl.zi.rm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.ba.sl.rda")




#USA(CA_FL) USBMK
taxtab.usbmk<-list()
taxtab6.usbmk<-list()
taxtab6.vag.usbmk<-list()
taxtab6.cs.usbmk<-list()
taxtab.ba.sl.usbmk<-list()
for (j in 1:length(dat.usbmk)){
  taxtab.usbmk[[j]]<-dat.usbmk[[j]][dat.usbmk[[j]]$mombb=="Baby"&dat.usbmk[[j]]$sampletype=="STL",]
  taxtab.usbmk[[j]]$age.sample<-taxtab.usbmk[[j]]$bbage/30
  taxtab.usbmk[[j]]$bf<-factor(taxtab.usbmk[[j]]$bfmixfm)
  levels(taxtab.usbmk[[j]]$bf)[levels(taxtab.usbmk[[j]]$bf)==""]<-NA
  taxtab.usbmk[[j]]$bf<-mapvalues(taxtab.usbmk[[j]]$bf,from=c("BF","FM","Mix","SOLID"),to=c("ExclusiveBF","No_BF","Non_exclusiveBF","No_BF"))
  taxtab.usbmk[[j]]$bf<-factor(taxtab.usbmk[[j]]$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
  taxtab.usbmk[[j]]$gender<-mapvalues(taxtab.usbmk[[j]]$babygender,from=c("F","M"),to=c("Female","Male"))
  taxtab.usbmk[[j]]$gender<-factor(taxtab.usbmk[[j]]$gender,levels=c("Female","Male"))
  taxtab.usbmk[[j]]$bm<-mapvalues(taxtab.usbmk[[j]]$delivery,from=c("C-section","Vaginal"),to=c("C_section","Vaginal"))
  taxtab.usbmk[[j]]$bm<-factor(taxtab.usbmk[[j]]$bm,levels=c("C_section","Vaginal"))
  taxtab.usbmk[[j]]$personid<-as.factor(taxtab.usbmk[[j]]$pairid)
  taxtab6.usbmk[[j]]<-taxtab.usbmk[[j]][taxtab.usbmk[[j]]$age.sample<=6,]
  taxtab6.vag.usbmk[[j]]<-subset(taxtab6.usbmk[[j]],delivery=="Vaginal")
  taxtab6.cs.usbmk[[j]]<-subset(taxtab6.usbmk[[j]],delivery=="C-section")
  #taxtab.ba.sl.usbmk[[j]]<-merge(ba.sl.mk[,c("x.sampleid","sl.intro","lastm.before.sl","sl.timediff","solid")],taxtab.usbmk[[j]],by="x.sampleid")
}
#get taxlist after filter
taxlist.usbmk<-taxa.filter(taxtab=taxtab.usbmk,percent.filter = 0.05, relabund.filter = 0.00005)

taxacom.usbmk<-taxa.compare(taxtab=taxtab.usbmk,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom6.usbmk<-taxa.compare(taxtab=taxtab6.usbmk,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.zi.usbmk<-taxa.compare(taxtab=taxtab.usbmk,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom6.zi.usbmk<-taxa.compare(taxtab=taxtab6.usbmk,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxa.meansdn.usbmk<-taxa.meansdn(taxtab=taxtab.usbmk)
#save(taxacom.usbmk,taxacom6.usbmk,taxacom.zi.usbmk,taxacom6.zi.usbmk,taxa.meansdn.usbmk,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.basic.usbmk.rda")
taxa.meansdn.bfbm.usbmk<-taxa.meansdn(taxtab=taxtab.usbmk,sumvar='bf',groupvar="bm")

#gender adjusted for bf and age
#bf customized function
taxacom6.usbmk.sex.adjustbfageo<-taxa.compare(taxtab=taxtab6.usbmk,propmed.rel="lm",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="yes")
#error if longitudinal =yes => no instead
taxacom6.zi.usbmk.sex.adjustbfageo<-taxa.compare(taxtab=taxtab6.usbmk,propmed.rel="gamlss",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="no")
#save(taxacom6.usbmk.sex.adjustbfageo,taxacom6.zi.usbmk.sex.adjustbfageo,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom6.zi.usbmk.sex.adjustbfageo.rda")
#general function
taxacom6.usbmk.sex.adjustbfage<-taxa.compare.gen(taxtab=taxtab6.usbmk[[5]],propmed.rel="lm",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="yes")
#error if longitudinal =yes => no instead
taxacom6.zi.usbmk.sex.adjustbfage<-taxa.compare.gen(taxtab=taxtab6.usbmk[[5]],propmed.rel="gamlss",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="no")
taxa.meansdn.usbmk.sexage<-taxa.meansdn(taxtab=taxtab.usbmk,sumvar="gender",groupvar="age.sample")
#save(taxacom6.usbmk.sex.adjustbfage,taxacom6.zi.usbmk.sex.adjustbfage,taxa.meansdn.usbmk.sexage,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom6.zi.usbmk.sex.adjustbfage.rda")

# ba sl
taxacom.ba.sl.usbmk<-taxa.compare(taxtab=taxtab.ba.sl.usbmk,taxsum="rel",propmed.rel="lm",comvar="solid",adjustvar="sl.intro",longitudinal="yes",p.adjust.method="fdr")
save(taxacom.ba.sl.usbmk,taxacom.ba.sl.rm,taxacom.ba.sl.zi.rm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.ba.sl.rda")

# delivery antibiotics
taxacom6.usbmk.m<-taxa.compare(taxtab=taxtab6.usbmk,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar=c("age.sample","delivery","abxatdelivery"),longitudinal="yes",p.adjust.method="fdr")
taxacom6.zi.usbmk.m<-taxa.compare(taxtab=taxtab6.usbmk,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar=c("age.sample","delivery","abxatdelivery"),longitudinal="yes",p.adjust.method="fdr")
#save(taxacom6.usbmk.m,taxacom6.zi.usbmk.m,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.multi.usbmk.rda")

#delivery adjusted
taxacom6.usbmk.bm<-taxa.compare(taxtab=taxtab6.usbmk,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="yes",p.adjust.method="fdr")
taxacom6.zi.usbmk.bm<-taxa.compare(taxtab=taxtab6.usbmk,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="yes",p.adjust.method="fdr")
#save(taxacom6.usbmk.bm,taxacom6.zi.usbmk.bm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.bm.usbmk.rda")

#stratifiction delivery
#vaginal
taxacom6.usbmk.vag<-taxa.compare(taxtab=taxtab6.vag.usbmk,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar=c("age.sample"),longitudinal="yes",p.adjust.method="fdr")
taxacom6.zi.usbmk.vag<-taxa.compare(taxtab=taxtab6.vag.usbmk,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar=c("age.sample"),longitudinal="yes",p.adjust.method="fdr")
#C-section
taxacom6.usbmk.cs<-taxa.compare(taxtab=taxtab6.cs.usbmk,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar=c("age.sample"),longitudinal="yes",p.adjust.method="fdr")
#taxacom6.zi.usbmk.cs<-taxa.compare(taxtab=taxtab6.cs.usbmk,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar=c("age.sample"),longitudinal="yes",p.adjust.method="fdr")
# Error with longitudinal, use no instead
taxacom6.zi.usbmk.cs<-taxa.compare(taxtab=taxtab6.cs.usbmk,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr")
#save(taxacom6.usbmk.vag,taxacom6.zi.usbmk.vag,taxacom6.usbmk.cs,taxacom6.zi.usbmk.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.vagcs.usbmk.rda")



# South Africa UW
taxtab.uw<-list()
for (j in 1:5){
  taxtab.uw[[j]]<-dat.uw[[j]]
  taxtab.uw[[j]]$age.sample<-as.numeric(as.character(substr(taxtab.uw[[j]]$timepoint,2,3)))/4
  taxtab.uw[[j]]$bf<-mapvalues(taxtab.uw[[j]]$ebf,from=c("0","1"),to=c("Non_exclusiveBF","ExclusiveBF"))
  taxtab.uw[[j]]$bf<-factor(taxtab.uw[[j]]$bf,levels=c("ExclusiveBF","Non_exclusiveBF"))
  taxtab.uw[[j]]$personid<-as.factor(taxtab.uw[[j]]$personid)
}
#get taxlist after filter
taxlist.uw<-taxa.filter(taxtab=taxtab.uw,percent.filter = 0.05, relabund.filter = 0.00005)

taxacom.uw<-taxa.compare(taxtab=taxtab.uw,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.zi.uw<-taxa.compare(taxtab=taxtab.uw,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxa.meansdn.uw<-taxa.meansdn(taxtab=taxtab.uw)
#save(taxacom.uw,taxacom.zi.uw,taxa.meansdn.uw,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.basic.uw.rda")


# USA(NC) UNC data
taxtab.unc<-list()
taxtab6.unc<-list()
for (j in 1:length(dat.unc)){
  taxtab.unc[[j]]<-dat.unc[[j]]
  taxtab.unc[[j]]$age.sample<-taxtab.unc[[j]]$agemo
  taxtab.unc[[j]]$bf<-factor(taxtab.unc[[j]]$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
  taxtab.unc[[j]]$gender<-mapvalues(taxtab.unc[[j]]$sex,from=c("0","1"),to=c("Male","Female"))
  taxtab.unc[[j]]$gender<-factor(taxtab.unc[[j]]$gender,levels=c("Female","Male"))
  taxtab.unc[[j]]$personid<-as.factor(taxtab.unc[[j]]$patient)
  taxtab6.unc[[j]]<-taxtab.unc[[j]][taxtab.unc[[j]]$age.sample<=6,] #&taxtab.unc[[j]]$bf!="No_BF"
  taxtab6.unc[[j]]$bf<-drop.levels(taxtab6.unc[[j]]$bf,reorder=FALSE)
}
#get taxlist after filter
taxlist.unc<-taxa.filter(taxtab=taxtab.unc,percent.filter = 0.05, relabund.filter = 0.00005)

taxacom.unc<-taxa.compare(taxtab=taxtab.unc,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom6.unc<-taxa.compare(taxtab=taxtab6.unc,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
taxacom.zi.unc<-taxa.compare(taxtab=taxtab.unc,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
#error due to small dataset=> use no longitudinal
taxacom6.zi.unc<-taxa.compare(taxtab=taxtab6.unc,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxa.meansdn.unc<-taxa.meansdn(taxtab=taxtab.unc)
#save(taxacom.unc,taxacom6.unc,taxacom.zi.unc,taxa.meansdn.unc, file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.basic.unc2.rda")
#include
#save(taxacom.unc,taxacom6.unc,taxacom.zi.unc,taxacom6.zi.unc,taxa.meansdn.unc, file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.basic.unc2.rda")

# gender adjusted for bf and age
#bf customized function
taxacom6.unc.sex.adjustedbfageo<-taxa.compare(taxtab=taxtab6.unc,propmed.rel="lm",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="yes")
taxacom6.zi.unc.sex.adjustedbfageo<-taxa.compare(taxtab=taxtab6.unc,propmed.rel="gamlss",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="no")
#save(taxacom6.unc.sex.adjustedbfageo,taxacom6.zi.unc.sex.adjustedbfageo,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom6.unc.sex.adjustedbfageo.rda")
#general function
taxacom6.unc.sex.adjustedbfage<-taxa.compare.gen(taxtab=taxtab6.unc[[5]],propmed.rel="lm",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="yes")
taxacom6.zi.unc.sex.adjustedbfage<-taxa.compare.gen(taxtab=taxtab6.unc[[5]],propmed.rel="gamlss",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="no")
taxa.meansdn.unc.sexage<-taxa.meansdn(taxtab=taxtab.unc,sumvar="gender",groupvar="age.sample")
#save(taxacom6.unc.sex.adjustedbfage,taxacom6.zi.unc.sex.adjustedbfage,taxa.meansdn.unc.sexage,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom6.unc.sex.adjustedbfage.rda")


#save(estisum.rm,estisum.ha,estisum.uw,estisum.usbmk,estisum.unc,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/estisum.rel.rda")
#save(taxacom.rm,taxacom.ha,taxacom.uw,taxacom.usbmk,taxacom.unc,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/taxacom.rel.rda")
#save(taxacom.zi.rm,taxacom.zi.ha,taxacom.zi.unc,taxacom.zi.uw,taxacom.zi.usbmk, file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/taxacom.rel.zi.rda")
#save(taxa.meansdn.unc,taxa.meansdn.uw,taxa.meansdn.usbmk,taxa.meansdn.rm,taxa.meansdn.ha, file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/taxa.meansdn.rda")

#save(taxacom.rm,taxacom.ha,taxacom.uw,taxacom.usbmk,taxacom.unc,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.rel.rda")
#save(taxacom6.rm,taxacom.ha,taxacom.uw,taxacom6.usbmk,taxacom6.unc,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom6.rel.rda")
#save(taxacom.zi.rm,taxacom.zi.ha,taxacom.zi.unc,taxacom.zi.uw,taxacom.zi.usbmk, file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.rel.zi.rda")
#save(taxacom6.zi.rm,taxacom.zi.ha,taxacom.zi.unc,taxacom.zi.uw,taxacom6.zi.usbmk, file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom6.rel.zi.rda")
#save(taxa.meansdn.unc,taxa.meansdn.uw,taxa.meansdn.usbmk,taxa.meansdn.rm,taxa.meansdn.ha, file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxa.meansdn.rda")
save(taxacom.rm,taxacom.ha,taxacom.uw,taxacom.usbmk,taxacom.unc,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.relfunc.rda")
save(taxacom6.rm,taxacom.ha,taxacom.uw,taxacom6.usbmk,taxacom6.unc,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom6.relfunc.rda")
save(taxacom.zi.rm,taxacom.zi.ha,taxacom.zi.unc,taxacom.zi.uw,taxacom.zi.usbmk, file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.rel.zifunc.rda")
save(taxacom6.zi.rm,taxacom.zi.ha,taxacom.zi.unc,taxacom6.zi.unc,taxacom.zi.uw,taxacom6.zi.usbmk, file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom6.rel.zifunc.rda")
save(taxa.meansdn.unc,taxa.meansdn.uw,taxa.meansdn.usbmk,taxa.meansdn.rm,taxa.meansdn.ha, file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxa.meansdnunc.rda")
#taxlist
save(taxlist.ha,taxlist.rm,taxlist.usbmk,taxlist.uw,taxlist.unc,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxlist.filter.rda")


# Canada (Ualberta)
patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UAlberta"
ca<-read.multi(patht=patht,patternt=".csv",assignt="yes",study="Canada")
bm<-read.csv("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UAlberta/Birth_Mode_BJOG.csv")
colnames(bm)<-tolower(colnames(bm))
comvar<-"bf"; adjustvar=c("age.sample");percent.filter=0.05;relabund.filter=0.00005
propmed.rel<-"lm"
taxdat<-merge(as.data.frame(ca$rel_abund_bjog),bm,by="sampleid")
taxdat<-subset(taxdat,birth.mode!='Vaginal')
taxdat$bf<-factor(taxdat$bf, levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
taxdat[,comvar]<-drop.levels(taxdat[,comvar],reorder=FALSE) #drop missing/unused level and keep level order
taxdat[,"comvarnum"]<-as.numeric(taxdat[,comvar]) # for trend testing
comvarnum<-"comvarnum"
levcom<-levels(taxdat[,comvar])
nlevcom<-nlevels(as.factor(as.character(taxdat[,comvar]))) # to remove empty level
# get assigned taxa only
taxlist<-colnames(taxdat)[grep("k__",colnames(taxdat))]
#filter using percent.filter
taxtest<-apply(taxdat[,taxlist],2,function(x){length(x[!is.na(x)&x>0])})
taxget<-taxtest[taxtest>=percent.filter*(nrow(taxdat))]
#taxname<-names(taxget)
#filter using relabund.filter
taxtestm<-apply(taxdat[,taxlist],2,mean,na.rm=T)
taxgetm<-taxtestm[taxtestm>relabund.filter]
taxname<-names(taxget)[names(taxget) %in% names(taxgetm)]
estisum<-matrix(NA,nrow=length(taxname),ncol=20)
colnames(estisum)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","pval.adjust.nebf",
                          "estimate.nbf","se.nbf","teststat.nbf","pval.nbf","pval.adjust.nbf",
                          "estimate.age","se.age","teststat.age","pval.age","pval.adjust.age",
                          "estimate.conbf","se.conbf","teststat.conbf","pval.conbf","pval.adjust.conbf")
rownames(estisum)<-taxname
for (i in 1:length(taxname)){
  #linear regression: not optimal test but work
  if (propmed.rel=="lm"){
    fitsum.ji<-try(summary(glm(as.formula(paste(taxname[i],paste(c(comvar,adjustvar),collapse="+"),sep="~")), data=taxdat,family="gaussian")))
    if (class(fitsum.ji) == "try-error") {
      cat("Error in model fit, NA introduced.\n")
      estisum[i,]<-rep(NA,ncol(estisum))
    }
    if (class(fitsum.ji) != "try-error") {
      if (nlevcom==3){
        estisum[i,c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.nbf","se.nbf","teststat.nbf","pval.nbf","estimate.age","se.age","teststat.age","pval.age")]<-c(fitsum.ji$coefficients[2,],fitsum.ji$coefficients[3,],fitsum.ji$coefficients[4,])
        #treat bf as continuous to test for trend
        fitsum.conbf.ji<-summary(glm(as.formula(paste(taxname[i],paste(c(comvarnum,adjustvar),collapse="+"),sep="~")), data=taxdat,family="gaussian"))
        estisum[i,c("estimate.conbf","se.conbf","teststat.conbf","pval.conbf")]<- fitsum.conbf.ji$coefficients[2,]
      }
      if (nlevcom==2){
        estisum[i,c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")]<-c(fitsum.ji$coefficients[2,],fitsum.ji$coefficients[3,])
      }
    }
  }
  #Generalized Additive Models for Location Scale and Shape: Betazeroinflated family, mu link logit
  if (propmed.rel=="gamlss"){
    testdat<-taxdat[,c(taxname[i],comvar,comvarnum,adjustvar)]
    testdat[,taxname[i]][testdat[,taxname[i]]==1]<-0.9999 # dirty fix for 1 value of relative abundance in UW data (to be checked)
    testdat<-testdat[!is.na(testdat[,comvar]),] #dirty fix for missing values of bf in usbmk data (to be checked)
    fitsum.ji<-try(summary(gamlss(as.formula(paste(taxname[i],paste(c(comvar,adjustvar),collapse="+"),sep="~")), family = BEZI, data = testdat, trace = FALSE),save=TRUE))
    if (class(fitsum.ji) == "try-error") {
      cat("Error in model fit, NA introduced.\n")
      estisum[i,]<-rep(NA,ncol(estisum))
    }
    if (class(fitsum.ji) != "try-error") {
      if (nlevcom==3){
        estisum[i,c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.nbf","se.nbf","teststat.nbf","pval.nbf","estimate.age","se.age","teststat.age","pval.age")]<-c(fitsum.ji$coef.table[2,],fitsum.ji$coef.table[3,],fitsum.ji$coef.table[4,])
        #treat bf as continuous to test for trend
        fitsum.conbf.ji<-summary(gamlss(as.formula(paste(taxname[i],paste(c(comvarnum,adjustvar),collapse="+"),sep="~")), family = BEZI, data = testdat, trace = FALSE),save=TRUE)
        estisum[i,c("estimate.conbf","se.conbf","teststat.conbf","pval.conbf")]<- fitsum.conbf.ji$coef.table[2,]
      }
      if (nlevcom==2){
        estisum[i,c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.age","se.age","teststat.age","pval.age")]<-c(fitsum.ji$coef.table[2,],fitsum.ji$coef.table[3,])
      }
    }
  }
}
estisum[,c("pval.adjust.nebf","pval.adjust.nbf","pval.adjust.age","pval.adjust.conbf")]<-apply(estisum[,c("pval.nebf","pval.nbf","pval.age","pval.conbf")],2, p.adjust, method = "fdr")
estisum<-estisum[order(estisum[,"pval.nebf"]),]
#taxacom.alb<-estisum
#taxacom.zi.alb<-estisum
#taxacom.zi.alb.bm<-estisum
#taxacom.alb.bm<-estisum
#taxacom.alb.vag<-estisum
#taxacom.zi.alb.vag<-estisum
#taxacom.zi.alb.cs<-estisum
taxacom.alb.cs<-estisum
#save(taxacom.alb,taxacom.zi.alb,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.alb.rda")
#save(taxacom.zi.alb.bm,taxacom.alb.bm,taxacom.alb.vag,taxacom.zi.alb.vag,taxacom.zi.alb.cs,taxacom.alb.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.alb.bmgavcs.rda")


# USA(CA_MA_MO) Harvard
taxtab.hav<-list()
taxtab.vag.hav<-list()
taxtab.cs.hav<-list()
for (j in 1:length(dat.hav)){
  taxtab.hav[[j]]<-dat.hav[[j]]
  taxtab.hav[[j]]$age.sample<-taxtab.hav[[j]]$age_in_days/30
  taxtab.hav[[j]]$bf<-mapvalues(as.factor(taxtab.hav[[j]]$breastfeeding),from=c("Excl_BF","Any_BF","No_BF"),to=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
  taxtab.hav[[j]]$bf<-factor(taxtab.hav[[j]]$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
  taxtab.hav[[j]]$bm<-mapvalues(taxtab.hav[[j]]$csec,from=c("0","1"),to=c("Vaginal","C-section"))
  taxtab.vag.hav[[j]]<-subset(taxtab.hav[[j]],csec=="0")
  taxtab.cs.hav[[j]]<-subset(taxtab.hav[[j]],csec=="1")
}
#rel
taxacom.hav<-taxa.compare(taxtab=taxtab.hav,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxacom.zi.hav<-taxa.compare(taxtab=taxtab.hav,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr")
taxa.meansdn.hav<-taxa.meansdn(taxtab=taxtab.hav)
#save(taxacom.hav,taxacom.zi.hav,taxa.meansdn.hav,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.basic.hav.rda")
taxa.meansdn.bfbm.hav<-taxa.meansdn(taxtab=taxtab.hav,sumvar="bf", groupvar="bm")
#save all meansdn bfbm 4 studies
#save(taxa.meansdn.bfbm.hav,taxa.meansdn.bfbm.ca,taxa.meansdn.bfbm.usbmk,taxa.meansdn.bfbm.ha,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxa.meansdn.bfbm4.rda")

#delivery adjusted
taxacom.hav.bm<-taxa.compare(taxtab=taxtab.hav,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar=c("age.sample","csec"),longitudinal="no",p.adjust.method="fdr")
taxacom.zi.hav.bm<-taxa.compare(taxtab=taxtab.hav,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar=c("age.sample","csec"),longitudinal="no",p.adjust.method="fdr")
#stratify delivery
#vaginal
taxacom.hav.vag<-taxa.compare(taxtab=taxtab.vag.hav,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr")
taxacom.zi.hav.vag<-taxa.compare(taxtab=taxtab.vag.hav,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr")
#c-section
taxacom.hav.cs<-taxa.compare(taxtab=taxtab.cs.hav,taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr")
taxacom.zi.hav.cs<-taxa.compare(taxtab=taxtab.cs.hav,taxsum="rel",propmed.rel="gamlss",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr")
#save(taxacom.hav.bm,taxacom.zi.hav.bm,taxacom.hav.vag,taxacom.zi.hav.vag,taxacom.hav.cs,taxacom.zi.hav.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.bmvagcs.hav.rda")




# alpha diversity
#RM data
patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/alpha_div_0601/alpha_div_collated"
alpha.rm<-read.multi(patht=patht,patternt=".txt",assignt="no",study="RM_Bangladesh")
samfile<-merge(samde, he50[,c("child.id","gender","month.exbf","month.food")],by="child.id")
samfile$age.sample<-samfile$age.months
samfile$bf<-factor(samfile$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
samfile$personid<-samfile$child.id
samfile$sampleid<-tolower(samfile$fecal.sample.id)
alphacom.rm<-alpha.compare(datlist=alpha.rm,depth=3,mapfile=samfile,mapsampleid="fecal.sample.id",comvar="bf",adjustvar="age.sample",longitudinal="yes")
alphacom6.rm<-alpha.compare(datlist=alpha.rm,depth=3,mapfile=samfile,mapsampleid="fecal.sample.id",comvar="bf",adjustvar="age.sample",longitudinal="yes",age.limit=6)
alphacom.rm.max<-alpha.compare(datlist=alpha.rm,depth="max",mapfile=samfile,mapsampleid="fecal.sample.id",comvar="bf",adjustvar="age.sample",longitudinal="yes")
#standardize
alphacom.rms<-alpha.compare(datlist=alpha.rm,depth=3,mapfile=samfile,mapsampleid="fecal.sample.id",comvar="bf",adjustvar="age.sample",longitudinal="yes",age.limit=100,standardize = TRUE)
alphacom6.rms<-alpha.compare(datlist=alpha.rm,depth=3,mapfile=samfile,mapsampleid="fecal.sample.id",comvar="bf",adjustvar="age.sample",longitudinal="yes",age.limit=6,standardize = TRUE)

#gender adjusted for bf and age
alphacom6.rm.sex<-alpha.compare(datlist=alpha.rm,depth=3,mapfile=samfile,mapsampleid="fecal.sample.id",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",age.limit=6)
alphacom6.rm.sexs<-alpha.compare(datlist=alpha.rm,depth=3,mapfile=samfile,mapsampleid="fecal.sample.id",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",age.limit=6,standardize=TRUE)
#generic function
alphacom6.rm.sexg<-alpha.compare.gen(datlist=alpha.rm,depth=3,mapfile=samfile,mapsampleid="fecal.sample.id",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",age.limit=6)
alphacom6.rm.sexsg<-alpha.compare.gen(datlist=alpha.rm,depth=3,mapfile=samfile,mapsampleid="fecal.sample.id",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",age.limit=6,standardize=TRUE)


# Haiti data
patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Haiti_data_010417/alpha_div_0517/alpha_div_collated"
alpha.ha<-read.multi(patht=patht,patternt=".txt",assignt="no",study="UCLA_Haiti")
haitimap$bf<-mapvalues(as.factor(haitimap$exclusivebf),from=c("Yes","No"),to=c("ExclusiveBF","Non_exclusiveBF"))
haitimap$bf<-factor(haitimap$bf,levels=c("ExclusiveBF","Non_exclusiveBF"))
haitimap$age.sample<-haitimap$bbage/30
haitimap$gender<-factor(haitimap$babygender,levels=c("Female","Male"))
haitimapfile<-haitimap[haitimap$mombb=="B"&haitimap$sampletype=="STL",]
#include only baby stool samples
alphacom.ha2<-alpha.compare(datlist=alpha.ha,depth=3,mapfile=haitimapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="no")
alphacom6.ha2<-alpha.compare(datlist=alpha.ha,depth=3,mapfile=haitimapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="no",age.limit=6)
alphacom.ha.max2<-alpha.compare(datlist=alpha.ha,depth="max",mapfile=haitimapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="no")
#save(alphacom.ha2,alphacom6.ha2,alphacom.ha.max2,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.ha2.rda")
#standardize
alphacom.has<-alpha.compare(datlist=alpha.ha,depth=3,mapfile=haitimapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="no",standardize = TRUE)
alphacom6.has<-alpha.compare(datlist=alpha.ha,depth=3,mapfile=haitimapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="no",age.limit=6,standardize = TRUE)

#gender adjusted for bf and age
alphacom6.ha.sex<-alpha.compare(datlist=alpha.ha,depth=3,mapfile=haitimapfile,mapsampleid="x.sampleid",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="no",age.limit=6)
alphacom6.ha.sexs<-alpha.compare(datlist=alpha.ha,depth=3,mapfile=haitimapfile,mapsampleid="x.sampleid",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="no",age.limit=6,standardize = TRUE)
#generic function
alphacom6.ha.sexg<-alpha.compare.gen(datlist=alpha.ha,depth=3,mapfile=haitimapfile,mapsampleid="x.sampleid",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="no",age.limit=6)
alphacom6.ha.sexsg<-alpha.compare.gen(datlist=alpha.ha,depth=3,mapfile=haitimapfile,mapsampleid="x.sampleid",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="no",age.limit=6,standardize = TRUE)

#delivery adjusted
alphacom.ha.bm<-alpha.compare(datlist=alpha.ha,depth=3,mapfile=haitimapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="no")
alphacom6.ha.bm<-alpha.compare(datlist=alpha.ha,depth=3,mapfile=haitimapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="no",age.limit=6)
alphacom.ha.max.bm<-alpha.compare(datlist=alpha.ha,depth="max",mapfile=haitimapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="no")
#stratification delivery
#vaginal
alphacom.ha.vag<-alpha.compare(datlist=alpha.ha,depth=3,mapfile=subset(haitimapfile,delivery=="Vaginal"),mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="no")
alphacom6.ha.vag<-alpha.compare(datlist=alpha.ha,depth=3,mapfile=subset(haitimapfile,delivery=="Vaginal"),mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="no",age.limit=6)
alphacom.ha.max.vag<-alpha.compare(datlist=alpha.ha,depth="max",mapfile=subset(haitimapfile,delivery=="Vaginal"),mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="no")
#save(alphacom.ha.bm,alphacom6.ha.bm,alphacom.ha.max.bm,alphacom.ha.vag,alphacom6.ha.vag,alphacom.ha.max.vag,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.ha.bmvag.rda")
#C-section => not doable due to alpha file does not contain the values for "2015b.stl" sample => need to be checked
alphacom.ha.cs<-alpha.compare(datlist=alpha.ha,depth=3,mapfile=subset(haitimapfile,delivery=="C-section"),mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="no")
alphacom6.ha.cs<-alpha.compare(datlist=alpha.ha,depth=3,mapfile=subset(haitimapfile,delivery=="C-section"),mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="no",age.limit=6)
alphacom.ha.max.cs<-alpha.compare(datlist=alpha.ha,depth="max",mapfile=subset(haitimapfile,delivery=="C-section"),mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="no")


# USA(CA_FL) USBMK
patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/USBMK/alpha_div_final/alpha_div_collated"
alpha.usbmk<-read.multi(patht=patht,patternt=".txt",assignt="no",study="UCLA_USBMK")
usbmkmapfile<-usbmkmap[usbmkmap$mombb=="Baby"&usbmkmap$sampletype=="STL",]
usbmkmapfile$age.sample<-usbmkmapfile$bbage/30
usbmkmapfile$bf<-factor(usbmkmapfile$bfmixfm)
levels(usbmkmapfile$bf)[levels(usbmkmapfile$bf)==""]<-NA
usbmkmapfile$bf<-mapvalues(usbmkmapfile$bf,from=c("BF","FM","Mix","SOLID"),to=c("ExclusiveBF","No_BF","Non_exclusiveBF","No_BF"))
usbmkmapfile$bf<-factor(usbmkmapfile$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
usbmkmapfile$personid<-usbmkmapfile$pairid
usbmkmapfile$gender<-mapvalues(usbmkmapfile$babygender,from=c("","F","M"),to=c(NA,"Female","Male"))
usbmkmapfile$gender<-factor(usbmkmapfile$gender,levels=c("Female","Male"))
alphacom.usbmk<-alpha.compare(datlist=alpha.usbmk,depth=3,mapfile=usbmkmapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes")
alphacom6.usbmk<-alpha.compare(datlist=alpha.usbmk,depth=3,mapfile=usbmkmapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes",age.limit=6)
alphacom.usbmk.max<-alpha.compare(datlist=alpha.usbmk,depth="max",mapfile=usbmkmapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes")
#standardize
alphacom.usbmks<-alpha.compare(datlist=alpha.usbmk,depth=3,mapfile=usbmkmapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes",standardize = TRUE)
alphacom6.usbmks<-alpha.compare(datlist=alpha.usbmk,depth=3,mapfile=usbmkmapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes",age.limit=6,standardize = TRUE)

#gender adjusted for bf and age
alphacom6.usbmk.sex<-alpha.compare(datlist=alpha.usbmk,depth=3,mapfile=usbmkmapfile,mapsampleid="x.sampleid",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",age.limit=6)
alphacom6.usbmk.sexs<-alpha.compare(datlist=alpha.usbmk,depth=3,mapfile=usbmkmapfile,mapsampleid="x.sampleid",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",age.limit=6,standardize = TRUE)
#generic functions
alphacom6.usbmk.sexg<-alpha.compare.gen(datlist=alpha.usbmk,depth=3,mapfile=usbmkmapfile,mapsampleid="x.sampleid",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",age.limit=6)
alphacom6.usbmk.sexsg<-alpha.compare.gen(datlist=alpha.usbmk,depth=3,mapfile=usbmkmapfile,mapsampleid="x.sampleid",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",age.limit=6,standardize = TRUE)

#delivery adjusted
alphacom.usbmk.bm<-alpha.compare(datlist=alpha.usbmk,depth=3,mapfile=usbmkmapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="yes")
alphacom6.usbmk.bm<-alpha.compare(datlist=alpha.usbmk,depth=3,mapfile=usbmkmapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="yes",age.limit=6)
alphacom.usbmk.max.bm<-alpha.compare(datlist=alpha.usbmk,depth="max",mapfile=usbmkmapfile,mapsampleid="x.sampleid",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="yes")
#stratification
#vaginal
alphacom.usbmk.vag<-alpha.compare(datlist=alpha.usbmk,depth=3,mapfile=subset(usbmkmapfile,delivery=="Vaginal"),mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes")
alphacom6.usbmk.vag<-alpha.compare(datlist=alpha.usbmk,depth=3,mapfile=subset(usbmkmapfile,delivery=="Vaginal"),mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes",age.limit=6)
alphacom.usbmk.max.vag<-alpha.compare(datlist=alpha.usbmk,depth="max",mapfile=subset(usbmkmapfile,delivery=="Vaginal"),mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes")
#C-section
alphacom.usbmk.cs<-alpha.compare(datlist=alpha.usbmk,depth=3,mapfile=subset(usbmkmapfile,delivery=="C-section"),mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes")
alphacom6.usbmk.cs<-alpha.compare(datlist=alpha.usbmk,depth=3,mapfile=subset(usbmkmapfile,delivery=="C-section"),mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes",age.limit=6)
alphacom.usbmk.max.cs<-alpha.compare(datlist=alpha.usbmk,depth="max",mapfile=subset(usbmkmapfile,delivery=="C-section"),mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes")
#save(alphacom.usbmk.bm,alphacom6.usbmk.bm,alphacom.usbmk.max.bm,
#     alphacom.usbmk.vag,alphacom6.usbmk.vag,alphacom.usbmk.max.vag,
#     alphacom.usbmk.cs,alphacom6.usbmk.cs,alphacom.usbmk.max.cs,
#     file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.usbmk.bmvagcs.rda")


# USA(NC) UNC Frontier
patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Frontier/alpha_div_all/alpha_div_collated"
alpha.unc<-read.multi(patht=patht,patternt=".txt",assignt="no",study="UNC_US")
uncmap$age.sample<-uncmap$agemo
uncmap$bf<-factor(uncmap$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
uncmap$personid<-uncmap$patient
uncmap$gender<-mapvalues(uncmap$sex,from=c("0","1"),to=c("Male","Female"))
uncmap$gender<-factor(uncmap$gender,levels=c("Female","Male"))
alphacom.unc<-alpha.compare(datlist=alpha.unc,depth=3,mapfile=uncmap,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes")
alphacom6.unc<-alpha.compare(datlist=alpha.unc,depth=3,mapfile=uncmap,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes",age.limit=6)
alphacom.unc.max<-alpha.compare(datlist=alpha.unc,depth="max",mapfile=uncmap,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes")
#standardize
alphacom.uncs<-alpha.compare(datlist=alpha.unc,depth=3,mapfile=uncmap,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes",standardize = TRUE)
alphacom6.uncs<-alpha.compare(datlist=alpha.unc,depth=3,mapfile=uncmap,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes",age.limit=6,standardize = TRUE)

#gender adjusted for bf and age
alphacom6.unc.sex<-alpha.compare(datlist=alpha.unc,depth=3,mapfile=uncmap,mapsampleid="x.sampleid",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",age.limit=6)
alphacom6.unc.sexs<-alpha.compare(datlist=alpha.unc,depth=3,mapfile=uncmap,mapsampleid="x.sampleid",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",age.limit=6,standardize = TRUE)
#generic functions
alphacom6.unc.sexg<-alpha.compare.gen(datlist=alpha.unc,depth=3,mapfile=uncmap,mapsampleid="x.sampleid",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",age.limit=6)
alphacom6.unc.sexsg<-alpha.compare.gen(datlist=alpha.unc,depth=3,mapfile=uncmap,mapsampleid="x.sampleid",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",age.limit=6,standardize = TRUE)
#save 4 studies
#save(alphacom6.unc.sex,alphacom6.usbmk.sex,alphacom6.ha.sex,alphacom6.rm.sex,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom6.sex4.rda")
#save(alphacom6.unc.sexs,alphacom6.usbmk.sexs,alphacom6.ha.sexs,alphacom6.rm.sexs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom6.sex4.scaled.rda")
#generic function
#save(alphacom6.unc.sexg,alphacom6.usbmk.sexg,alphacom6.ha.sexg,alphacom6.rm.sexg,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom6.sex4g.rda")
#save(alphacom6.unc.sexsg,alphacom6.usbmk.sexsg,alphacom6.ha.sexsg,alphacom6.rm.sexsg,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom6.sex4.scaledg.rda")


# South Africa UW
patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UW/alpha_div/alpha_div_collated"
alpha.uw<-read.multi(patht=patht,patternt=".txt",assignt="no",study="UW_SouthAfrica")
uwmap$age.sample<-as.numeric(as.character(substr(uwmap$timepoint,2,3)))/4
uwmap$bf<-mapvalues(uwmap$ebf,from=c("0","1"),to=c("Non_exclusiveBF","ExclusiveBF"))
uwmap$bf<-factor(uwmap$bf,levels=c("ExclusiveBF","Non_exclusiveBF"))
alphacom.uw<-alpha.compare(datlist=alpha.uw,depth=3,mapfile=uwmap,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes")
alphacom6.uw<-alpha.compare(datlist=alpha.uw,depth=3,mapfile=uwmap,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes",age.limit=6)
alphacom.uw.max<-alpha.compare(datlist=alpha.uw,depth="max",mapfile=uwmap,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes")
#standardize
alphacom.uws<-alpha.compare(datlist=alpha.uw,depth=3,mapfile=uwmap,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes",standardize = TRUE)
alphacom6.uws<-alpha.compare(datlist=alpha.uw,depth=3,mapfile=uwmap,mapsampleid="x.sampleid",comvar="bf",adjustvar="age.sample",longitudinal="yes",age.limit=6,standardize = TRUE)

#mean sd n
alphamean.uw<-alphacom.uw$alphamean

#save(alphacom.rm,alphacom.rm.max,alphacom.ha,alphacom.ha.max,alphacom.usbmk,alphacom.usbmk.max,alphacom.unc,alphacom.unc.max,alphacom.uw,alphacom.uw.max,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/alphacom.rda")
#save(alphacom6.rm,alphacom6.ha,alphacom6.usbmk,alphacom6.unc,alphacom6.uw,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom6.rda")
#fix factor level order messing up isse
#save(alphacom.rm,alphacom.rm.max,alphacom.ha,alphacom.ha.max,alphacom.usbmk,alphacom.usbmk.max,alphacom.unc,alphacom.unc.max,alphacom.uw,alphacom.uw.max,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/alphacomfunc.rda")
#save(alphacom6.rm,alphacom6.ha,alphacom6.usbmk,alphacom6.unc,alphacom6.uw,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom6func.rda")

# Canada (Ualberta)
patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UAlberta"
ca<-read.multi(patht=patht,patternt=".csv",assignt="yes",study="Canada")
bm<-read.csv("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UAlberta/Birth_Mode_BJOG.csv")
colnames(bm)<-tolower(colnames(bm))
#alphamap<-merge(ca$alphamean_bjog,bm,by="sampleid")
#standardize
alphameans<-mutate_at(ca$alphamean_bjog,.vars=c("chao1","obeserved_species","pd_whole_tree","shannon"),.funs=function(x){(x-mean(x,na.r=T))/sd(x,na.rm=T)})
alphamap<-merge(alphameans,bm,by="sampleid")
alphamap$bf<-factor(alphamap$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
#alphamap<-subset(alphamap,birth.mode!="Vaginal")
comvar="bf";adjustvar=c("age.sample")
alphamap[,"comvarnum"]<-as.numeric(alphamap[,comvar]) # for trend testing
comvarnum<-"comvarnum"
nlevcom<-nlevels(as.factor(as.character(alphamap[,"bf"])))
aindex<-c("chao1","obeserved_species","pd_whole_tree","shannon")
alphasum<-matrix(NA,ncol=16,nrow=length(aindex))
rownames(alphasum)<-aindex
colnames(alphasum)<-c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.nbf","se.nbf","teststat.nbf","pval.nbf","estimate.age","se.age","teststat.age","pval.age","estimate.conbf","se.conbf","teststat.conbf","pval.conbf")
for (j in 1:length(aindex)){
  sumfit.j<-summary(glm(as.formula(paste(aindex[j],paste(c(comvar,adjustvar),collapse="+"),sep="~")), data=alphamap,family="gaussian"))
  alphasum[j,c("estimate.nebf","se.nebf","teststat.nebf","pval.nebf","estimate.nbf","se.nbf","teststat.nbf","pval.nbf","estimate.age","se.age","teststat.age","pval.age")]<-c(sumfit.j$coefficients[2,],sumfit.j$coefficients[3,],sumfit.j$coefficients[4,])
  #treat bf as continuous to test for trend adjusted for age
  sumfit.conbf.j<-summary(glm(as.formula(paste(aindex[j],paste(c(comvarnum,adjustvar),collapse="+"),sep="~")), data=alphamap,family="gaussian"))
  alphasum[j,c("estimate.conbf","se.conbf","teststat.conbf","pval.conbf")]<-sumfit.conbf.j$coefficients[2,]
}
#alphacom.alb<-list(alphamean=alphamap,alphasum=alphasum)
#alphasum.alb<-alphasum
#save(alphacom.alb,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.alb.rda")
#standardize
alphacom.albs<-list(alphamean=alphamap,alphasum=alphasum)

#birthmode adjusted
#alphasum.alb.bm<-alphasum
#stratification birthmode
#vaginal
#alphasum.alb.vag<-alphasum
#C-section
#alphasum.alb.cs<-alphasum
#save(alphasum.alb.bm,alphasum.alb.vag,alphasum.alb.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.alb.bmvagcs.rda")


# USA(CA_MA_MO) Harvard
patht<-"C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Harvard/ForNhan/SentToNhan/alpha_div_collated"
alpha.hav<-read.multi(patht=patht,patternt=".txt",assignt="no",study="US(CA_MA_MO)")
harvardmap$bf<-mapvalues(harvardmap$breastfeeding,from=c("Excl_BF","Any_BF","No_BF"),to=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
harvardmap$bf<-factor(harvardmap$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
harvardmap$age.sample<-harvardmap$age_in_days/30
harvardmapfile<-harvardmap
alphacom.hav<-alpha.compare(datlist=alpha.hav,depth=3,mapfile=harvardmap,mapsampleid="ANONID",comvar="bf",adjustvar="age.sample",longitudinal="no")
alphacom6.hav<-alpha.compare(datlist=alpha.hav,depth=3,mapfile=harvardmap,mapsampleid="ANONID",comvar="bf",adjustvar="age.sample",longitudinal="no",age.limit=6)
alphacom.hav.max<-alpha.compare(datlist=alpha.hav,depth="max",mapfile=harvardmap,mapsampleid="ANONID",comvar="bf",adjustvar="age.sample",longitudinal="no")
#save(alphacom.hav,alphacom6.hav,alphacom.hav.max,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.hav.rda")
#standardize
alphacom.havs<-alpha.compare(datlist=alpha.hav,depth=3,mapfile=harvardmap,mapsampleid="sampleid",comvar="bf",adjustvar="age.sample",longitudinal="no",standardize = TRUE)
alphacom6.havs<-alpha.compare(datlist=alpha.hav,depth=3,mapfile=harvardmap,mapsampleid="sampleid",comvar="bf",adjustvar="age.sample",longitudinal="no",age.limit=6,standardize = TRUE)

#save all standardize estimates from all 7 studies
#save(alphacom.havs,alphacom6.havs,alphacom.albs,alphacom.uws,alphacom6.uws,alphacom.uncs,alphacom6.uncs,alphacom.usbmks,alphacom6.usbmks,
#     alphacom.has,alphacom6.has,alphacom.rms,alphacom6.rms,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.scaled7.rda")

#delivery adjusted
alphacom.hav.bm<-alpha.compare(datlist=alpha.hav,depth=3,mapfile=harvardmapfile,mapsampleid="ANONID",comvar="bf",adjustvar=c("age.sample","csec"),longitudinal="no")
alphacom6.hav.bm<-alpha.compare(datlist=alpha.hav,depth=3,mapfile=harvardmapfile,mapsampleid="ANONID",comvar="bf",adjustvar=c("age.sample","csec"),longitudinal="no",age.limit=6)
alphacom.hav.max.bm<-alpha.compare(datlist=alpha.hav,depth="max",mapfile=harvardmapfile,mapsampleid="ANONID",comvar="bf",adjustvar=c("age.sample","csec"),longitudinal="no")
#stratification delivery
#vaginal
alphacom.hav.vag<-alpha.compare(datlist=alpha.hav,depth=3,mapfile=subset(harvardmapfile,csec=="0"),mapsampleid="ANONID",comvar="bf",adjustvar="age.sample",longitudinal="no")
alphacom6.hav.vag<-alpha.compare(datlist=alpha.hav,depth=3,mapfile=subset(harvardmapfile,csec=="0"),mapsampleid="ANONID",comvar="bf",adjustvar="age.sample",longitudinal="no",age.limit=6)
alphacom.hav.max.vag<-alpha.compare(datlist=alpha.hav,depth="max",mapfile=subset(harvardmapfile,csec=="0"),mapsampleid="ANONID",comvar="bf",adjustvar="age.sample",longitudinal="no")
#C-section
alphacom.hav.cs<-alpha.compare(datlist=alpha.hav,depth=3,mapfile=subset(harvardmapfile,csec=="1"),mapsampleid="ANONID",comvar="bf",adjustvar="age.sample",longitudinal="no")
alphacom6.hav.cs<-alpha.compare(datlist=alpha.hav,depth=3,mapfile=subset(harvardmapfile,csec=="1"),mapsampleid="ANONID",comvar="bf",adjustvar="age.sample",longitudinal="no",age.limit=6)
alphacom.hav.max.cs<-alpha.compare(datlist=alpha.hav,depth="max",mapfile=subset(harvardmapfile,csec=="1"),mapsampleid="ANONID",comvar="bf",adjustvar="age.sample",longitudinal="no")
#save(alphacom.hav.bm,alphacom6.hav.bm,alphacom.hav.max.bm,alphacom.hav.vag,alphacom6.hav.vag,alphacom.hav.max.vag,alphacom.hav.cs,alphacom6.hav.cs,alphacom.hav.max.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.hav.bmvag.rda")


# combine for vag cs
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.alb.bmvagcs.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.hav.bmvag.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.usbmk.bmvagcs.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.ha.bmvag.rda"))
# vaginal
alphasum.alb.vag<-as.data.frame(alphasum.alb.vag)
alphasum.alb.vag[,"index"]<-rownames(alphasum.alb.vag)
alphasum.alb.vag[,"study"]<-" Azad et al 2015 (Canada)"
alphasum.ha.vag<-as.data.frame(alphacom.ha.vag$alphasum)
alphasum.ha.vag[,"index"]<-rownames(alphasum.ha.vag)
alphasum.ha.vag[,"study"]<-"Bender et al 2016 (Haiti)"
alphasum.usbmk.vag<-as.data.frame(alphacom.usbmk.vag$alphasum)
alphasum.usbmk.vag[,"index"]<-rownames(alphasum.usbmk.vag)
alphasum.usbmk.vag[,"study"]<-"Pannaraj et al 2017 (USA(CA_FL))"
alphasum.hav.vag<-as.data.frame(alphacom.hav.vag$alphasum)
alphasum.hav.vag[,"index"]<-rownames(alphasum.hav.vag)
alphasum.hav.vag[,"study"]<-"Sordillo et al 2017 (USA(CA_MA_MO))"
alphasum.vag<-rbind.fill(alphasum.alb.vag,alphasum.ha.vag,alphasum.usbmk.vag,alphasum.hav.vag)
# C-section
alphasum.alb.cs<-as.data.frame(alphasum.alb.cs)
alphasum.alb.cs[,"index"]<-rownames(alphasum.alb.cs)
alphasum.alb.cs[,"study"]<-" Azad et al 2015 (Canada)"
alphasum.usbmk.cs<-as.data.frame(alphacom.usbmk.cs$alphasum)
alphasum.usbmk.cs[,"index"]<-rownames(alphasum.usbmk.cs)
alphasum.usbmk.cs[,"study"]<-"Pannaraj et al 2017 (USA(CA_FL))"
alphasum.hav.cs<-as.data.frame(alphacom.hav.cs$alphasum)
alphasum.hav.cs[,"index"]<-rownames(alphasum.hav.cs)
alphasum.hav.cs[,"study"]<-"Sordillo et al 2017 (USA(CA_MA_MO))"
alphasum.cs<-rbind.fill(alphasum.alb.cs,alphasum.usbmk.cs,alphasum.hav.cs)
#save(alphasum.vag,alphasum.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphasum.vagcs.rda")

# combine for main alpha see in meta-analysis section







# metaanalysis
library(meta)
#library(metafor)
#library(rmeta)

# alpha diversity 
#print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/alphacom.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom6func.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.ha2.rda")) #baby stool sample only
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.alb.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.hav.rda"))
#standardize
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.scaled7.rda"))
#Haiti
#as.ha<-as.data.frame(alphacom6.ha2$alphasum)
as.ha<-as.data.frame(alphacom6.has$alphasum)
as.ha$index<-rownames(as.ha)
as.ha$author<-"Bender et al"
as.ha$pop<-"Haiti"
as.ha$year<-"2016"
#USBMK
#as.usbmk<-as.data.frame(alphacom6.usbmk$alphasum)
as.usbmk<-as.data.frame(alphacom6.usbmks$alphasum)
as.usbmk$index<-rownames(as.usbmk)
as.usbmk$author<-"Pannaraj et al"
as.usbmk$pop<-"USA(CA_FL)"
as.usbmk$year<-"2017"
#RM
#as.rm<-as.data.frame(alphacom6.rm$alphasum)
as.rm<-as.data.frame(alphacom6.rms$alphasum)
as.rm$index<-rownames(as.rm)
as.rm$author<-"Subramanian et al"
as.rm$pop<-"Bangladesh"
as.rm$year<-"2014"
#UNC
#as.unc<-as.data.frame(alphacom6.unc$alphasum)
as.unc<-as.data.frame(alphacom6.uncs$alphasum)
as.unc$index<-rownames(as.unc)
as.unc$author<-"Thompson et al"
as.unc$pop<-"USA(NC)"
as.unc$year<-"2015"
#UW
#as.uw<-as.data.frame(alphacom6.uw$alphasum)
as.uw<-as.data.frame(alphacom6.uws$alphasum)
as.uw$index<-rownames(as.uw)
as.uw$author<-"Wood et al"
as.uw$pop<-"South Africa"
as.uw$year<-"2017"
#Ualberta => fix the names
#as.alb<-as.data.frame(alphacom.alb$alphasum)
as.alb<-as.data.frame(alphacom.albs$alphasum)
rownames(as.alb)[rownames(as.alb)=="obeserved_species"]<-"observed_species"
as.alb$index<-rownames(as.alb)
as.alb$author<-"Azad et al"
as.alb$pop<-"Canada"
as.alb$year<-"2015"
#Harvard
#as.hav<-as.data.frame(alphacom6.hav$alphasum)
as.hav<-as.data.frame(alphacom6.havs$alphasum)
as.hav$index<-rownames(as.hav)
as.hav$author<-"Sordillo et al"
as.hav$pop<-"USA(CA_MA_MO)"
as.hav$year<-"2017"

#alphaall<-rbind.fill(as.ha,as.usbmk,as.uw,as.unc,as.rm,as.alb,as.hav)
#save(alphaall, file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphaall7.rda")
#standardize
alphaall<-rbind.fill(as.rm,as.alb,as.ha,as.uw,as.usbmk,as.hav,as.unc)
#save(alphaall, file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphaall7s.rda")


# run all meta-analysis, store results for summary heatmap, forest plot
adat<-alphaall
adat$study<-paste(adat$author,adat$year,"(",adat$pop,")",sep=" ")
aname<-unique(adat[,"index"])
metatab.f<-matrix(NA,ncol=28, nrow=length(aname))
colnames(metatab.f)<-c("estimate.nebf",'se.nebf','ll.nebf','ul.nebf','z.nebf','p.nebf','p.adjust.nebf',
                            "estimate.nbf",'se.nbf','ll.nbf','ul.nbf','z.nbf','p.nbf','p.adjust.nbf',
                            "estimate.age",'se.age','ll.age','ul.age','z.age','p.age','p.adjust.age',
                            "estimate.conbf",'se.conbf','ll.conbf','ul.conbf','z.conbf','p.conbf','p.adjust.conbf')
rownames(metatab.f)<-aname
metatab.r<-matrix(NA,ncol=28, nrow=length(aname))
colnames(metatab.r)<-c("estimate.nebf",'se.nebf','ll.nebf','ul.nebf','z.nebf','p.nebf','p.adjust.nebf',
                            "estimate.nbf",'se.nbf','ll.nbf','ul.nbf','z.nbf','p.nbf','p.adjust.nbf',
                            "estimate.age",'se.age','ll.age','ul.age','z.age','p.age','p.adjust.age',
                            "estimate.conbf",'se.conbf','ll.conbf','ul.conbf','z.conbf','p.conbf','p.adjust.conbf')
rownames(metatab.r)<-aname
for (j in 1:length(aname)){
  testdat<-subset(adat,index %in% aname[j])
  fit.nebf<-metagen(estimate.nebf, se.nebf, studlab=paste(author,year,"(",pop,")"),data=testdat,sm="RD", backtransf=FALSE)
  metatab.f[j,c("estimate.nebf",'se.nebf','ll.nebf','ul.nebf','z.nebf','p.nebf')]<-c(summary(fit.nebf)$fixed$TE,summary(fit.nebf)$fixed$seTE,summary(fit.nebf)$fixed$lower,summary(fit.nebf)$fixed$upper,summary(fit.nebf)$fixed$z,summary(fit.nebf)$fixed$p)
  metatab.r[j,c("estimate.nebf",'se.nebf','ll.nebf','ul.nebf','z.nebf','p.nebf')]<-c(summary(fit.nebf)$random$TE,summary(fit.nebf)$random$seTE,summary(fit.nebf)$random$lower,summary(fit.nebf)$random$upper,summary(fit.nebf)$random$z,summary(fit.nebf)$random$p)
  fit.nbf<-metagen(estimate.nbf, se.nbf, studlab=paste(author,year,"(",pop,")"),data=testdat,sm="RD", backtransf=FALSE)
  metatab.f[j,c("estimate.nbf",'se.nbf','ll.nbf','ul.nbf','z.nbf','p.nbf')]<-c(summary(fit.nbf)$fixed$TE,summary(fit.nbf)$fixed$seTE,summary(fit.nbf)$fixed$lower,summary(fit.nbf)$fixed$upper,summary(fit.nbf)$fixed$z,summary(fit.nbf)$fixed$p)
  metatab.r[j,c("estimate.nbf",'se.nbf','ll.nbf','ul.nbf','z.nbf','p.nbf')]<-c(summary(fit.nbf)$random$TE,summary(fit.nbf)$random$seTE,summary(fit.nbf)$random$lower,summary(fit.nbf)$random$upper,summary(fit.nbf)$random$z,summary(fit.nbf)$random$p)
  fit.age<-metagen(estimate.age, se.age, studlab=paste(author,year,"(",pop,")"),data=testdat,sm="RD", backtransf=FALSE)
  metatab.f[j,c("estimate.age",'se.age','ll.age','ul.age','z.age','p.age')]<-c(summary(fit.age)$fixed$TE,summary(fit.age)$fixed$seTE,summary(fit.age)$fixed$lower,summary(fit.age)$fixed$upper,summary(fit.age)$fixed$z,summary(fit.age)$fixed$p)
  metatab.r[j,c("estimate.age",'se.age','ll.age','ul.age','z.age','p.age')]<-c(summary(fit.age)$random$TE,summary(fit.age)$random$seTE,summary(fit.age)$random$lower,summary(fit.age)$random$upper,summary(fit.age)$random$z,summary(fit.age)$random$p)
  fit.conbf<-metagen(estimate.conbf, se.conbf, studlab=paste(author,year,"(",pop,")"),data=testdat,sm="RD", backtransf=FALSE)
  metatab.f[j,c("estimate.conbf",'se.conbf','ll.conbf','ul.conbf','z.conbf','p.conbf')]<-c(summary(fit.conbf)$fixed$TE,summary(fit.conbf)$fixed$seTE,summary(fit.conbf)$fixed$lower,summary(fit.conbf)$fixed$upper,summary(fit.conbf)$fixed$z,summary(fit.conbf)$fixed$p)
  metatab.r[j,c("estimate.conbf",'se.conbf','ll.conbf','ul.conbf','z.conbf','p.conbf')]<-c(summary(fit.conbf)$random$TE,summary(fit.conbf)$random$seTE,summary(fit.conbf)$random$lower,summary(fit.conbf)$random$upper,summary(fit.conbf)$random$z,summary(fit.conbf)$random$p)
}
metatab.f[,c("p.adjust.nebf","p.adjust.nbf","p.adjust.age","p.adjust.conbf")]<-apply(metatab.f[,c("p.nebf","p.nbf","p.age","p.conbf")],2, p.adjust, method = "fdr")
metatab.f<-metatab.f[order(metatab.f[,"p.nebf"]),]
metatab.r[,c("p.adjust.nebf","p.adjust.nbf","p.adjust.age","p.adjust.conbf")]<-apply(metatab.r[,c("p.nebf","p.nbf","p.age","p.conbf")],2, p.adjust, method = "fdr")
metatab.r<-metatab.r[order(metatab.r[,"p.nebf"]),]
alpha.metatab.fr<-list(metatab.f,metatab.r)
names(alpha.metatab.fr)<-c("fixed","random")
a.metatab.r<-as.data.frame(metatab.r)
a.metatab.r$study<-"Pooled estimate"
a.metatab.r$index<-rownames(a.metatab.r)
adat.p<-rbind.fill(adat,a.metatab.r)
#save(adat,a.metatab.r,adat.p,alphap,alphaall, file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alpha.use.unc.rda")
#7 studies
#save(adat,a.metatab.r,adat.p,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphameta.allindex7.rda")
#standardize
#save(adat,a.metatab.r,adat.p,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphameta.allindex7s.rda")


# Pooled alpha data for pooled/meta-analysis
# alpha
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/alphacomfunc.rda"))
#print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom6func.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.ha2.rda")) #baby stool sample only
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.alb.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.hav.rda"))

#standardize
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphacom.scaled7.rda"))
# Haiti map for baby stool samples only
haitimap$personid<-as.factor(paste("ha",haitimap$pairid,sep="."))
haitimap$age.sample<-as.numeric(as.character(haitimap$bbage/30))
haitimap$bf<-mapvalues(haitimap$exclusivebf,from=c("Yes","No"),to=c("ExclusiveBF","Non_exclusiveBF"))
haitimap$bf<-factor(haitimap$bf,levels=c("ExclusiveBF","Non_exclusiveBF"))
haitimap$gender<-factor(haitimap$babygender,levels=c("Female","Male"))
haitimap$author<-"Bender et al"
haitimap$year<-"2016"
haitimap$pop<-"Haiti"
haitimap$sampleid<-paste("ha",tolower(haitimap$x.sampleid),sep=".")
hamap.bs<-haitimap[haitimap$mombb=="B"&haitimap$sampletype=="STL",c("sampleid","personid","age.sample","bf","author","year","pop","gender")]
#alphamean.ha<-alphacom.ha$alphamean
#alphamean.ha<-alphacom.ha2$alphamean
#standardize
alphamean.ha<-alphacom.has$alphamean.standardized
alphamean.ha$sampleid<-paste("ha",alphamean.ha$sampleid,sep=".")
alpha.m.ha<-merge(hamap.bs,alphamean.ha,by="sampleid")
#usbmk
usbmkmap$personid<-as.factor(paste("usbmk",usbmkmap$pairid,sep="."))
usbmkmap$age.sample<-as.numeric(as.character(usbmkmap$bbage/30))
usbmkmap$gender<-mapvalues(usbmkmap$babygender, from=c("F","M"), to=c("Female","Male"))
usbmkmap$gender<-factor(usbmkmap$gender,levels=c("Female","Male"))
usbmkmap$bf<-factor(usbmkmap$bfmixfm)
levels(usbmkmap$bf)[levels(usbmkmap$bf)==""]<-NA
usbmkmap$bf<-mapvalues(usbmkmap$bf,from=c("BF","FM","Mix","SOLID"),to=c("ExclusiveBF","No_BF","Non_exclusiveBF","No_BF"))
usbmkmap$bf<-factor(usbmkmap$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
usbmkmap$sampleid<-paste("usbmk",tolower(usbmkmap$x.sampleid),sep=".")
usbmkmap$author<-"Pannaraj et al"
usbmkmap$year<-"2017"
usbmkmap$pop<-"USA(CA_FL)"
usbmkmap.bs<-usbmkmap[usbmkmap$mombb=="Baby"&usbmkmap$sampletype=="STL",c("sampleid","personid","age.sample","bf","author","year","pop","gender")]
#alphamean.usbmk<-alphacom.usbmk$alphamean
#standardize
alphamean.usbmk<-alphacom.usbmks$alphamean.standardized
alphamean.usbmk$sampleid<-paste("usbmk",alphamean.usbmk$sampleid,sep=".")
alpha.m.usbmk<-merge(usbmkmap.bs, alphamean.usbmk, by="sampleid")
#unc
uncmap$personid<-as.factor(paste("unc",tolower(uncmap$patient),sep="."))
uncmap$sampleid<-paste("unc",tolower(uncmap$x.sampleid),sep=".")
uncmap$bf<-factor(uncmap$bf, levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
uncmap$age.sample<-as.numeric(as.character(uncmap$agemo))
uncmap$gender<-mapvalues(as.factor(uncmap$sex), from=c("0","1"),to=c("Male","Female"))
uncmap$gender<-factor(uncmap$gender,levels=c("Female","Male"))
uncmap$author<-"Thompson et al"
uncmap$year<-"2015"
uncmap$pop<-"USA(NC)"
#alphamean.unc<-alphacom.unc$alphamean
#standardize
alphamean.unc<-alphacom6.uncs$alphamean.standardized
alphamean.unc$sampleid<-paste("unc",tolower(alphamean.unc$sampleid),sep=".")
alpha.m.unc<-merge(uncmap[,c("sampleid","personid","age.sample","bf","author","year","pop","gender")],alphamean.unc,by='sampleid')
#uw
uwmap$personid<-as.factor(paste("uw",tolower(uwmap$personid),sep="."))
uwmap$sampleid<-paste("uw",tolower(uwmap$x.sampleid),sep=".")
uwmap$age.sample<-as.numeric(as.character(substr(uwmap$timepoint,2,3)))/4
uwmap$bf<-mapvalues(uwmap$ebf,from=c("0","1"),to=c("Non_exclusiveBF","ExclusiveBF"))
uwmap$bf<-factor(uwmap$bf,levels=c("ExclusiveBF","Non_exclusiveBF"))
uwmap$author<-"Wood et al"
uwmap$year<-"2017"
uwmap$pop<-"South Africa"
#alphamean.uw<-alphacom.uw$alphamean
#standardize
alphamean.uw<-alphacom.uws$alphamean.standardized
alphamean.uw$sampleid<-paste("uw",tolower(alphamean.uw$sampleid),sep=".")
alpha.m.uw<-merge(uwmap[,c("sampleid","personid","age.sample","bf","author","year","pop")],alphamean.uw,by="sampleid")
#Subramanian
samfile<-merge(samde, he50[,c("child.id","gender","month.exbf","month.food")],by="child.id")
samfile$age.sample<-samfile$age.months
samfile$bf<-factor(samfile$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
samfile$personid<-as.factor(paste("rm", tolower(samfile$child.id),sep="."))
samfile$sampleid<- paste("rm",tolower(samfile$fecal.sample.id),sep=".")
samfile$author<-"Subramanian et al"
samfile$year<-"2014"
samfile$pop<-"Bangladesh"
#alphamean.rm<-alphacom.rm$alphamean
#standardize
alphamean.rm<-alphacom6.rms$alphamean.standardized
alphamean.rm$sampleid<-paste("rm",tolower(alphamean.rm$sampleid),sep=".")
#alpha.m.rm<-merge(samfile[,c("sampleid","personid","age.sample","bf","author","year","pop","gender","month.exbf","month.food")],alphamean.rm,by='sampleid')
alpha.m.rm<-merge(samfile,alphamean.rm,by='sampleid')

#Harvard
harvardmap$personid<-as.factor(harvardmap$anonid)
harvardmap$sampleid<-tolower(harvardmap$anonid)
harvardmap$age.sample<-harvardmap$age_in_days/30
harvardmap$bf<-mapvalues(harvardmap$breastfeeding,from=c("Excl_BF","Any_BF","No_BF"),to=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
harvardmap$bf<-factor(harvardmap$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
harvardmap$author<-"Sordillo et al"
harvardmap$year<-"2017"
harvardmap$pop<-"USA(CA_MA_MO)"
#alphamean.hav<-alphacom.hav$alphamean
#standardize
alphamean.hav<-alphacom6.havs$alphamean.standardized
alpha.m.hav<-merge(harvardmap,alphamean.hav,by="sampleid")

#Ualberta
#alphamean.alb<-alphacom.alb$alphamean
#standardize
alphamean.alb<-alphacom.albs$alphamean
colnames(alphamean.alb)[colnames(alphamean.alb)=="obeserved_species"]<-"observed_species"
alphamean.alb$bf<-factor(alphamean.alb$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
alphamean.alb$author<-"Azad et al"
alphamean.alb$year<-"2015"
alphamean.alb$pop<-"Canada"
alphamean.alb$sampleid<-paste("ca",tolower(alphamean.alb$sampleid),sep=".")
alphamean.alb$personid<-paste("ca",tolower(alphamean.alb$personid),sep=".")
alpha.m.alb<-alphamean.alb

#pooled
#alphap<-rbind.fill(alpha.m.rm,alpha.m.uw,alpha.m.unc,alpha.m.usbmk,alpha.m.ha,alpha.m.hav,alpha.m.alb)
#save(alphap,alpha.m.alb,alpha.m.hav,alpha.m.uw,alpha.m.unc,alpha.m.usbmk,alpha.m.rm,alpha.m.ha,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphamean7.pooled.rda")
#standardize
alphap<-rbind.fill(alpha.m.rm,alpha.m.alb,alpha.m.ha,alpha.m.uw,alpha.m.usbmk,alpha.m.hav,alpha.m.unc)
#save(alphap,alpha.m.alb,alpha.m.hav,alpha.m.uw,alpha.m.unc,alpha.m.usbmk,alpha.m.rm,alpha.m.ha,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/alphamean7s.pooled.rda")




# Taxa relative abundance
# Heatmap of estimate (RR or RD) for the taxa significant in at least one study
# Across studies and pooled (meta-analysis)
# Each heatmap for each level L2-L6
# put data together for each level
#print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom6.rel.zi.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom6.rel.zifunc.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.alb.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.basic.hav.rda"))

#RM
taxacom.zi.rm<-taxacom6.zi.rm
for (i in 1: length(names(taxacom.zi.rm))){
  taxacom.zi.rm[[i]]<-as.data.frame(taxacom.zi.rm[[i]])
  taxacom.zi.rm[[i]][,'taxa']<-rownames(taxacom.zi.rm[[i]])
  taxacom.zi.rm[[i]][,'author']<-"Subramanian et al"
  taxacom.zi.rm[[i]][,'year']<-"2014"
  taxacom.zi.rm[[i]][,'pop']<-"Bangladesh"
}
#Haiti
for (i in 1: length(names(taxacom.zi.ha))){
  taxacom.zi.ha[[i]]<-as.data.frame(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'taxa']<-rownames(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'author']<-"Bender et al"
  taxacom.zi.ha[[i]][,'year']<-"2016"
  taxacom.zi.ha[[i]][,'pop']<-"Haiti"
}
#UW
for (i in 1: length(names(taxacom.zi.uw))){
  taxacom.zi.uw[[i]]<-as.data.frame(taxacom.zi.uw[[i]])
  taxacom.zi.uw[[i]][,'taxa']<-rownames(taxacom.zi.uw[[i]])
  taxacom.zi.uw[[i]][,'author']<-"Wood et al"
  taxacom.zi.uw[[i]][,'year']<-"2017"
  taxacom.zi.uw[[i]][,'pop']<-"South Africa"
}
# UNC
#include nobf in data<=6 month
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.basic.unc2.rda"))
taxacom.zi.unc<-taxacom6.zi.unc
for (i in 1:5 ){ #length(names(taxacom.zi.unc))
  taxacom.zi.unc[[i]]<-as.data.frame(taxacom.zi.unc[[i]])
  taxacom.zi.unc[[i]][,'taxa']<-rownames(taxacom.zi.unc[[i]])
  taxacom.zi.unc[[i]][,'author']<-"Thompson et al"
  taxacom.zi.unc[[i]][,'year']<-"2015"
  taxacom.zi.unc[[i]][,'pop']<-"USA(NC)"
}
#USBMK
taxacom.zi.usbmk<-taxacom6.zi.usbmk
for (i in 1: length(names(taxacom.zi.usbmk))){
  taxacom.zi.usbmk[[i]]<-as.data.frame(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'taxa']<-rownames(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'author']<-"Pannaraj et al"
  taxacom.zi.usbmk[[i]][,'year']<-"2017"
  taxacom.zi.usbmk[[i]][,'pop']<-"USA(CA_FL)"
}
#Harvard
taxacom.zi.hav<-taxacom.zi.hav
for (i in 1: length(names(taxacom.zi.hav))){
  taxacom.zi.hav[[i]]<-as.data.frame(taxacom.zi.hav[[i]])
  taxacom.zi.hav[[i]][,'taxa']<-rownames(taxacom.zi.hav[[i]])
  taxacom.zi.hav[[i]][,'author']<-"Sordillo et al"
  taxacom.zi.hav[[i]][,'year']<-"2017"
  taxacom.zi.hav[[i]][,'pop']<-"USA(CA_MA_MO)"
}
#Alberta
taxacom.zi.alb1<-taxacom.zi.alb
al2<-rownames(taxacom.zi.alb1)[-grep("c__",rownames(taxacom.zi.alb1))]
al3<-rownames(taxacom.zi.alb1)[-grep("o__",rownames(taxacom.zi.alb1))]
al4<-rownames(taxacom.zi.alb1)[-grep("f__",rownames(taxacom.zi.alb1))]
al5<-rownames(taxacom.zi.alb1)[-grep("g__",rownames(taxacom.zi.alb1))]
taxacom.zi.alb2<-list(l2=taxacom.zi.alb1[al2,],l3=taxacom.zi.alb1[al3,],l4=taxacom.zi.alb1[al4,],l5=taxacom.zi.alb1[al5,],l6=taxacom.zi.alb1)
for (i in 1: length(names(taxacom.zi.alb2))){
  taxacom.zi.alb2[[i]]<-as.data.frame(taxacom.zi.alb2[[i]])
  taxacom.zi.alb2[[i]][,'taxa']<-rownames(taxacom.zi.alb2[[i]])
  taxacom.zi.alb2[[i]][,'author']<-"Azad et al"
  taxacom.zi.alb2[[i]][,'year']<-"2015"
  taxacom.zi.alb2[[i]][,'pop']<-"Canada"
}

taxacom.zi.l2<-rbind.fill(taxacom.zi.rm$l2,taxacom.zi.ha$l2,taxacom.zi.unc$l2,taxacom.zi.uw$l2,taxacom.zi.usbmk$l2,taxacom.zi.hav$l2,taxacom.zi.alb2$l2)
taxacom.zi.l3<-rbind.fill(taxacom.zi.rm$l3,taxacom.zi.ha$l3,taxacom.zi.unc$l3,taxacom.zi.uw$l3,taxacom.zi.usbmk$l3,taxacom.zi.hav$l3,taxacom.zi.alb2$l3)
taxacom.zi.l4<-rbind.fill(taxacom.zi.rm$l4,taxacom.zi.ha$l4,taxacom.zi.unc$l4,taxacom.zi.uw$l4,taxacom.zi.usbmk$l4,taxacom.zi.hav$l4,taxacom.zi.alb2$l4)
taxacom.zi.l5<-rbind.fill(taxacom.zi.rm$l5,taxacom.zi.ha$l5,taxacom.zi.unc$l5,taxacom.zi.uw$l5,taxacom.zi.usbmk$l5,taxacom.zi.hav$l5,taxacom.zi.alb2$l5)
taxacom.zi.l6<-rbind.fill(taxacom.zi.rm$l6,taxacom.zi.ha$l6,taxacom.zi.unc$l6,taxacom.zi.uw$l6,taxacom.zi.usbmk$l6,taxacom.zi.hav$l6,taxacom.zi.alb2$l6)
taxacom.zi<-list(taxacom.zi.l2,taxacom.zi.l3,taxacom.zi.l4,taxacom.zi.l5,taxacom.zi.l6)
names(taxacom.zi)<-paste("l",2:6,sep="")
# metaanalysis for gamlss
metatab.zi<-meta.taxa(taxcomdat=taxacom.zi, sm="RR",p.adjust.method="fdr",percent.meta=0.5)
#save(metatab.zi,taxacom.zi,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/metatab.zi7.rda")

# Sensitivity analysis: meta-analysis without unc data
taxacom.zi.l2.nounc<-rbind.fill(taxacom.zi.rm$l2,taxacom.zi.ha$l2,taxacom.zi.uw$l2,taxacom.zi.usbmk$l2,taxacom.zi.hav$l2,taxacom.zi.alb2$l2)
taxacom.zi.l3.nounc<-rbind.fill(taxacom.zi.rm$l3,taxacom.zi.ha$l3,taxacom.zi.uw$l3,taxacom.zi.usbmk$l3,taxacom.zi.hav$l3,taxacom.zi.alb2$l3)
taxacom.zi.l4.nounc<-rbind.fill(taxacom.zi.rm$l4,taxacom.zi.ha$l4,taxacom.zi.uw$l4,taxacom.zi.usbmk$l4,taxacom.zi.hav$l4,taxacom.zi.alb2$l4)
taxacom.zi.l5.nounc<-rbind.fill(taxacom.zi.rm$l5,taxacom.zi.ha$l5,taxacom.zi.uw$l5,taxacom.zi.usbmk$l5,taxacom.zi.hav$l5,taxacom.zi.alb2$l5)
taxacom.zi.l6.nounc<-rbind.fill(taxacom.zi.rm$l6,taxacom.zi.ha$l6,taxacom.zi.uw$l6,taxacom.zi.usbmk$l6,taxacom.zi.hav$l6,taxacom.zi.alb2$l6)
taxacom.zi.nounc<-list(taxacom.zi.l2.nounc,taxacom.zi.l3.nounc,taxacom.zi.l4.nounc,taxacom.zi.l5.nounc,taxacom.zi.l6.nounc)
names(taxacom.zi.nounc)<-paste("l",2:6,sep="")
# metaanalysis for gamlss
metatab.zi.nounc<-meta.taxa(taxcomdat=taxacom.zi.nounc, sm="RR",p.adjust.method="fdr",percent.meta=0.5)
#save(taxacom.zi.nounc,metatab.zi.nounc,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/metatab.zi.nounc.rda")
#7 studies
#save(taxacom.zi.nounc,metatab.zi.nounc,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/metatab.zi.nounc7.rda")

# Sensitivity analysis: meta-analysis without ha data
taxacom.zi.l2.noha<-rbind.fill(taxacom.zi.rm$l2,taxacom.zi.unc$l2,taxacom.zi.uw$l2,taxacom.zi.usbmk$l2,taxacom.zi.hav$l2,taxacom.zi.alb2$l2)
taxacom.zi.l3.noha<-rbind.fill(taxacom.zi.rm$l3,taxacom.zi.unc$l3,taxacom.zi.uw$l3,taxacom.zi.usbmk$l3,taxacom.zi.hav$l3,taxacom.zi.alb2$l3)
taxacom.zi.l4.noha<-rbind.fill(taxacom.zi.rm$l4,taxacom.zi.unc$l4,taxacom.zi.uw$l4,taxacom.zi.usbmk$l4,taxacom.zi.hav$l4,taxacom.zi.alb2$l4)
taxacom.zi.l5.noha<-rbind.fill(taxacom.zi.rm$l5,taxacom.zi.unc$l5,taxacom.zi.uw$l5,taxacom.zi.usbmk$l5,taxacom.zi.hav$l5,taxacom.zi.alb2$l5)
taxacom.zi.l6.noha<-rbind.fill(taxacom.zi.rm$l6,taxacom.zi.unc$l6,taxacom.zi.uw$l6,taxacom.zi.usbmk$l6,taxacom.zi.hav$l6,taxacom.zi.alb2$l6)
taxacom.zi.noha<-list(taxacom.zi.l2.noha,taxacom.zi.l3.noha,taxacom.zi.l4.noha,taxacom.zi.l5.noha,taxacom.zi.l6.noha)
names(taxacom.zi.noha)<-paste("l",2:6,sep="")
# metaanalysis for gamlss
metatab.zi.noha<-meta.taxa(taxcomdat=taxacom.zi.noha, sm="RR",p.adjust.method="fdr",percent.meta=0.5)
#save(taxacom.zi.noha,metatab.zi.noha,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/metatab.zi.noha.func.rda")
#7 studies
#save(taxacom.zi.noha,metatab.zi.noha,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/metatab.zi.noha7.rda")

# Sensitivity analysis: no harvard (trial)
taxacom.zi.l2.nohav<-rbind.fill(taxacom.zi.rm$l2,taxacom.zi.ha$l2,taxacom.zi.unc$l2,taxacom.zi.uw$l2,taxacom.zi.usbmk$l2,taxacom.zi.alb2$l2)
taxacom.zi.l3.nohav<-rbind.fill(taxacom.zi.rm$l3,taxacom.zi.ha$l3,taxacom.zi.unc$l3,taxacom.zi.uw$l3,taxacom.zi.usbmk$l3,taxacom.zi.alb2$l3)
taxacom.zi.l4.nohav<-rbind.fill(taxacom.zi.rm$l4,taxacom.zi.ha$l4,taxacom.zi.unc$l4,taxacom.zi.uw$l4,taxacom.zi.usbmk$l4,taxacom.zi.alb2$l4)
taxacom.zi.l5.nohav<-rbind.fill(taxacom.zi.rm$l5,taxacom.zi.ha$l5,taxacom.zi.unc$l5,taxacom.zi.uw$l5,taxacom.zi.usbmk$l5,taxacom.zi.alb2$l5)
taxacom.zi.l6.nohav<-rbind.fill(taxacom.zi.rm$l6,taxacom.zi.ha$l6,taxacom.zi.unc$l6,taxacom.zi.uw$l6,taxacom.zi.usbmk$l6,taxacom.zi.alb2$l6)
taxacom.zi.nohav<-list(taxacom.zi.l2.nohav,taxacom.zi.l3.nohav,taxacom.zi.l4.nohav,taxacom.zi.l5.nohav,taxacom.zi.l6.nohav)
names(taxacom.zi.nohav)<-paste("l",2:6,sep="")
# metaanalysis for gamlss
metatab.zi.nohav<-meta.taxa(taxcomdat=taxacom.zi.nohav, sm="RR",p.adjust.method="fdr",percent.meta=0.5)
#save(taxacom.zi.nohav,metatab.zi.nohav,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/metatab.zi.nohav7.rda")


# meta-analysis for 4 studies with gender info but not adjusting for gender (only age): for comparison
load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.basic.rm.rda")
load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.basic.ha.rda")
load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.basic.usbmk.rda")
load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.basic.unc2.rda")
taxacom.zi.l2<-rbind.fill(taxacom.zi.rm$l2,taxacom.zi.ha$l2,taxacom.zi.unc$l2,taxacom.zi.usbmk$l2)
taxacom.zi.l3<-rbind.fill(taxacom.zi.rm$l3,taxacom.zi.ha$l3,taxacom.zi.unc$l3,taxacom.zi.usbmk$l3)
taxacom.zi.l4<-rbind.fill(taxacom.zi.rm$l4,taxacom.zi.ha$l4,taxacom.zi.unc$l4,taxacom.zi.usbmk$l4)
taxacom.zi.l5<-rbind.fill(taxacom.zi.rm$l5,taxacom.zi.ha$l5,taxacom.zi.unc$l5,taxacom.zi.usbmk$l5)
taxacom.zi.l6<-rbind.fill(taxacom.zi.rm$l6,taxacom.zi.ha$l6,taxacom.zi.unc$l6,taxacom.zi.usbmk$l6)
taxacom.zi.4com<-list(taxacom.zi.l2,taxacom.zi.l3,taxacom.zi.l4,taxacom.zi.l5,taxacom.zi.l6)
names(taxacom.zi.4com)<-paste("l",2:6,sep="")
# metaanalysis for gamlss
metatab.zi.4com<-meta.taxa(taxcomdat=taxacom.zi.4com, sm="RR",p.adjust.method="fdr",percent.meta=0.5)
#save(metatab.zi.4com,taxacom.zi.4com,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/metatab.zi.4com.rda")



# gender adjusted for bf and age
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.rm.sex.adjustbfage.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.ha.sex.adjustbfage.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom6.zi.usbmk.sex.adjustbfage.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom6.unc.sex.adjustedbfage.rda"))
taxacom6.zi.rm.sex.adjustbfage$study<-"Subramanian et al 2014 (Bangladesh)"
taxacom6.zi.rm.sex.adjustbfage$pop<-"Bangladesh"
taxacom.zi.ha.sex.adjustbfage$study<-"Bender et al 2016 (Haiti)"
taxacom.zi.ha.sex.adjustbfage$pop<-"Haiti"
taxacom6.zi.usbmk.sex.adjustbfage$study<-"Pannaraj et al 2017 (USA(CA_FL))"
taxacom6.zi.usbmk.sex.adjustbfage$pop<-"USA(CA_FL)"
taxacom6.zi.unc.sex.adjustedbfage$study<-"Thompson et al 2015 (USA(NC))"
taxacom6.zi.unc.sex.adjustedbfage$pop<-"USA(NC)"
tabsex4<-rbind.fill(taxacom6.zi.rm.sex.adjustbfage,taxacom.zi.ha.sex.adjustbfage,taxacom6.zi.usbmk.sex.adjustbfage,taxacom6.zi.unc.sex.adjustedbfage)
# use general function
metab.sex<-meta.taxa.gen(taxcomdat=tabsex4,summary.measure="RR",pool.var="id",studylab="study",backtransform=FALSE,percent.meta=0.5,p.adjust.method="fdr")
#save(metab.sex,tabsex4,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/metab.sex.rda")


# use bf customized function
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.rm.sex.adjustbfageo.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.ha.sex.adjustbfageo.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom6.zi.usbmk.sex.adjustbfageo.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom6.unc.sex.adjustedbfageo.rda"))
#RM
taxacom.zi.rm<-taxacom6.zi.rm.sex.adjustbfageo
for (i in 1: length(names(taxacom.zi.rm))){
  taxacom.zi.rm[[i]]<-as.data.frame(taxacom.zi.rm[[i]])
  taxacom.zi.rm[[i]][,'taxa']<-rownames(taxacom.zi.rm[[i]])
  taxacom.zi.rm[[i]][,'author']<-"Subramanian et al"
  taxacom.zi.rm[[i]][,'year']<-"2014"
  taxacom.zi.rm[[i]][,'pop']<-"Bangladesh"
}
#Haiti
taxacom.zi.ha<-taxacom.zi.ha.sex.adjustbfageo
for (i in 1: length(names(taxacom.zi.ha))){
  taxacom.zi.ha[[i]]<-as.data.frame(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'taxa']<-rownames(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'author']<-"Bender et al"
  taxacom.zi.ha[[i]][,'year']<-"2016"
  taxacom.zi.ha[[i]][,'pop']<-"Haiti"
}
# UNC
taxacom.zi.unc<-taxacom6.zi.unc.sex.adjustedbfageo
for (i in 1:5 ){ #length(names(taxacom.zi.unc))
  taxacom.zi.unc[[i]]<-as.data.frame(taxacom.zi.unc[[i]])
  taxacom.zi.unc[[i]][,'taxa']<-rownames(taxacom.zi.unc[[i]])
  taxacom.zi.unc[[i]][,'author']<-"Thompson et al"
  taxacom.zi.unc[[i]][,'year']<-"2015"
  taxacom.zi.unc[[i]][,'pop']<-"USA(NC)"
}
#USBMK
taxacom.zi.usbmk<-taxacom6.zi.usbmk.sex.adjustbfageo
for (i in 1: length(names(taxacom.zi.usbmk))){
  taxacom.zi.usbmk[[i]]<-as.data.frame(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'taxa']<-rownames(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'author']<-"Pannaraj et al"
  taxacom.zi.usbmk[[i]][,'year']<-"2017"
  taxacom.zi.usbmk[[i]][,'pop']<-"USA(CA_FL)"
}
taxacom.zi.l2<-rbind.fill(taxacom.zi.rm$l2,taxacom.zi.ha$l2,taxacom.zi.unc$l2,taxacom.zi.usbmk$l2)
taxacom.zi.l3<-rbind.fill(taxacom.zi.rm$l3,taxacom.zi.ha$l3,taxacom.zi.unc$l3,taxacom.zi.usbmk$l3)
taxacom.zi.l4<-rbind.fill(taxacom.zi.rm$l4,taxacom.zi.ha$l4,taxacom.zi.unc$l4,taxacom.zi.usbmk$l4)
taxacom.zi.l5<-rbind.fill(taxacom.zi.rm$l5,taxacom.zi.ha$l5,taxacom.zi.unc$l5,taxacom.zi.usbmk$l5)
taxacom.zi.l6<-rbind.fill(taxacom.zi.rm$l6,taxacom.zi.ha$l6,taxacom.zi.unc$l6,taxacom.zi.usbmk$l6)
taxacom.zi.sex<-list(taxacom.zi.l2,taxacom.zi.l3,taxacom.zi.l4,taxacom.zi.l5,taxacom.zi.l6)
names(taxacom.zi.sex)<-paste("l",2:6,sep="")
# metaanalysis for gamlss
metatab.zi.sex.adjustbfageo<-meta.taxa(taxcomdat=taxacom.zi.sex, sm="RR",p.adjust.method="fdr",percent.meta=0.5)
#save(metatab.zi.sex.adjustbfageo,taxacom.zi.sex,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/metatab.zi.sex.adjustbfageo.rda")


# stratification by birth mode
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.bmvagcs.hav.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.alb.bmgavcs.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.vagcs.usbmk.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/taxacom.vagcs.ha.rda"))
# vaginal
#Haiti
taxacom.zi.ha<-taxacom.zi.ha.vag
for (i in 1: length(names(taxacom.zi.ha))){
  taxacom.zi.ha[[i]]<-as.data.frame(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'taxa']<-rownames(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'author']<-"Bender et al"
  taxacom.zi.ha[[i]][,'year']<-"2016"
  taxacom.zi.ha[[i]][,'pop']<-"Haiti"
}
#USBMK
taxacom.zi.usbmk<-taxacom6.zi.usbmk.vag
for (i in 1: length(names(taxacom.zi.usbmk))){
  taxacom.zi.usbmk[[i]]<-as.data.frame(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'taxa']<-rownames(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'author']<-"Pannaraj et al"
  taxacom.zi.usbmk[[i]][,'year']<-"2017"
  taxacom.zi.usbmk[[i]][,'pop']<-"USA(CA_FL)"
}
#Harvard
taxacom.zi.hav<-taxacom.zi.hav.vag
for (i in 1: length(names(taxacom.zi.hav))){
  taxacom.zi.hav[[i]]<-as.data.frame(taxacom.zi.hav[[i]])
  taxacom.zi.hav[[i]][,'taxa']<-rownames(taxacom.zi.hav[[i]])
  taxacom.zi.hav[[i]][,'author']<-"Sordillo et al"
  taxacom.zi.hav[[i]][,'year']<-"2017"
  taxacom.zi.hav[[i]][,'pop']<-"USA(CA_MA_MO)"
}
#Alberta
taxacom.zi.alb1<-taxacom.zi.alb.vag
al2<-rownames(taxacom.zi.alb1)[-grep("c__",rownames(taxacom.zi.alb1))]
al3<-rownames(taxacom.zi.alb1)[-grep("o__",rownames(taxacom.zi.alb1))]
al4<-rownames(taxacom.zi.alb1)[-grep("f__",rownames(taxacom.zi.alb1))]
al5<-rownames(taxacom.zi.alb1)[-grep("g__",rownames(taxacom.zi.alb1))]
taxacom.zi.alb2<-list(l2=taxacom.zi.alb1[al2,],l3=taxacom.zi.alb1[al3,],l4=taxacom.zi.alb1[al4,],l5=taxacom.zi.alb1[al5,],l6=taxacom.zi.alb1)
for (i in 1: length(names(taxacom.zi.alb2))){
  taxacom.zi.alb2[[i]]<-as.data.frame(taxacom.zi.alb2[[i]])
  taxacom.zi.alb2[[i]][,'taxa']<-rownames(taxacom.zi.alb2[[i]])
  taxacom.zi.alb2[[i]][,'author']<-"Azad et al"
  taxacom.zi.alb2[[i]][,'year']<-"2015"
  taxacom.zi.alb2[[i]][,'pop']<-"Canada"
}

taxacom.zi.l2<-rbind.fill(taxacom.zi.ha$l2,taxacom.zi.usbmk$l2,taxacom.zi.hav$l2,taxacom.zi.alb2$l2)
taxacom.zi.l3<-rbind.fill(taxacom.zi.ha$l3,taxacom.zi.usbmk$l3,taxacom.zi.hav$l3,taxacom.zi.alb2$l3)
taxacom.zi.l4<-rbind.fill(taxacom.zi.ha$l4,taxacom.zi.usbmk$l4,taxacom.zi.hav$l4,taxacom.zi.alb2$l4)
taxacom.zi.l5<-rbind.fill(taxacom.zi.ha$l5,taxacom.zi.usbmk$l5,taxacom.zi.hav$l5,taxacom.zi.alb2$l5)
taxacom.zi.l6<-rbind.fill(taxacom.zi.ha$l6,taxacom.zi.usbmk$l6,taxacom.zi.hav$l6,taxacom.zi.alb2$l6)
taxacom.zi.vag<-list(taxacom.zi.l2,taxacom.zi.l3,taxacom.zi.l4,taxacom.zi.l5,taxacom.zi.l6)
names(taxacom.zi.vag)<-paste("l",2:6,sep="")
# metaanalysis for gamlss
metatab.zi.vag<-meta.taxa(taxcomdat=taxacom.zi.vag, sm="RR",p.adjust.method="fdr",percent.meta=0.5)

# C-section
#Haiti
taxacom.zi.ha<-taxacom.zi.ha.cs
for (i in 1: length(names(taxacom.zi.ha))){
  taxacom.zi.ha[[i]]<-as.data.frame(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'taxa']<-rownames(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'author']<-"Bender et al"
  taxacom.zi.ha[[i]][,'year']<-"2016"
  taxacom.zi.ha[[i]][,'pop']<-"Haiti"
}
#USBMK
taxacom.zi.usbmk<-taxacom6.zi.usbmk.cs
for (i in 1: length(names(taxacom.zi.usbmk))){
  taxacom.zi.usbmk[[i]]<-as.data.frame(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'taxa']<-rownames(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'author']<-"Pannaraj et al"
  taxacom.zi.usbmk[[i]][,'year']<-"2017"
  taxacom.zi.usbmk[[i]][,'pop']<-"USA(CA_FL)"
}
#Harvard
taxacom.zi.hav<-taxacom.zi.hav.cs
for (i in 1: length(names(taxacom.zi.hav))){
  taxacom.zi.hav[[i]]<-as.data.frame(taxacom.zi.hav[[i]])
  taxacom.zi.hav[[i]][,'taxa']<-rownames(taxacom.zi.hav[[i]])
  taxacom.zi.hav[[i]][,'author']<-"Sordillo et al"
  taxacom.zi.hav[[i]][,'year']<-"2017"
  taxacom.zi.hav[[i]][,'pop']<-"USA(CA_MA_MO)"
}
#Alberta
taxacom.zi.alb1<-taxacom.zi.alb.cs
al2<-rownames(taxacom.zi.alb1)[-grep("c__",rownames(taxacom.zi.alb1))]
al3<-rownames(taxacom.zi.alb1)[-grep("o__",rownames(taxacom.zi.alb1))]
al4<-rownames(taxacom.zi.alb1)[-grep("f__",rownames(taxacom.zi.alb1))]
al5<-rownames(taxacom.zi.alb1)[-grep("g__",rownames(taxacom.zi.alb1))]
taxacom.zi.alb2<-list(l2=taxacom.zi.alb1[al2,],l3=taxacom.zi.alb1[al3,],l4=taxacom.zi.alb1[al4,],l5=taxacom.zi.alb1[al5,],l6=taxacom.zi.alb1)
for (i in 1: length(names(taxacom.zi.alb2))){
  taxacom.zi.alb2[[i]]<-as.data.frame(taxacom.zi.alb2[[i]])
  taxacom.zi.alb2[[i]][,'taxa']<-rownames(taxacom.zi.alb2[[i]])
  taxacom.zi.alb2[[i]][,'author']<-"Azad et al"
  taxacom.zi.alb2[[i]][,'year']<-"2015"
  taxacom.zi.alb2[[i]][,'pop']<-"Canada"
}

taxacom.zi.l2<-rbind.fill(taxacom.zi.ha$l2,taxacom.zi.usbmk$l2,taxacom.zi.hav$l2,taxacom.zi.alb2$l2)
taxacom.zi.l3<-rbind.fill(taxacom.zi.ha$l3,taxacom.zi.usbmk$l3,taxacom.zi.hav$l3,taxacom.zi.alb2$l3)
taxacom.zi.l4<-rbind.fill(taxacom.zi.ha$l4,taxacom.zi.usbmk$l4,taxacom.zi.hav$l4,taxacom.zi.alb2$l4)
taxacom.zi.l5<-rbind.fill(taxacom.zi.ha$l5,taxacom.zi.usbmk$l5,taxacom.zi.hav$l5,taxacom.zi.alb2$l5)
taxacom.zi.l6<-rbind.fill(taxacom.zi.ha$l6,taxacom.zi.usbmk$l6,taxacom.zi.hav$l6,taxacom.zi.alb2$l6)
taxacom.zi.cs<-list(taxacom.zi.l2,taxacom.zi.l3,taxacom.zi.l4,taxacom.zi.l5,taxacom.zi.l6)
names(taxacom.zi.cs)<-paste("l",2:6,sep="")
# metaanalysis for gamlss
metatab.zi.cs<-meta.taxa(taxcomdat=taxacom.zi.cs, sm="RR",p.adjust.method="fdr",percent.meta=0.5)
#save(metatab.zi.cs,taxacom.zi.cs,metatab.zi.vag,taxacom.zi.vag,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/metatab.zi.vagcs.rda")



# Kegg pathway
# Bangladesh
load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/sam.rm.rda")
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/open_uclust6/metagenome_at_level3_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke3<-as.data.frame(t(kegg))
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/open_uclust6/metagenome_at_level2_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke2<-as.data.frame(t(kegg))
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/open_uclust6/metagenome_at_level1_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke1<-as.data.frame(t(kegg))
kegg.rm<-list(ke1,ke2,ke3)
covar.rm<-merge(samde, he50[,c("child.id","gender","zygosity","day.firstsample","day.lastsample","n.sample","sampling.interval.msd","month.exbf","month.food",
                               "n.diarrhea.yr","percent.time.diarrhea","fraction.antibiotic","subject.allocation")], by="child.id")
covar.rm<-covar.rm %>% rename(sampleid=fecal.sample.id, personid=child.id ,age.sample=age.months)
covar.rm$bf<-factor(covar.rm$bf, levels=c('ExclusiveBF','Non_exclusiveBF','No_BF'))

covar.rm$personid<-as.factor(covar.rm$personid)
covar.rm<-covar.rm %>% group_by(personid) %>% arrange(personid,age.sample)  %>%
  mutate(month.food6=cut(month.food, breaks=c(-Inf, 6, Inf), labels=c("<=6 months",">6 months")),
         month.food5=cut(month.food, breaks=c(-Inf, 5, Inf), labels=c("<=5 months",">5 months")),
         month.food4=cut(month.food, breaks=c(-Inf, 4, Inf), labels=c("<=4 months",">4 months")),
         month.foodr=as.factor(as.character(round(month.food,0))),
         month.exbf3=cut(month.exbf, breaks=c(-Inf, 3, Inf), labels=c("<=3 months",">3 months")),
         month.exbf2=cut(month.exbf, breaks=c(-Inf, 2, Inf), labels=c("<=2 months",">2 months")),
         month.exbf1=cut(month.exbf, breaks=c(-Inf, 1, Inf), labels=c("<=1 months",">1 months")),
         month.exbfr=as.factor(as.character(round(month.exbf,0))))
# all age
pathcom.rm.log<-pathway.compare(pathtab=kegg.rm,mapfile=covar.rm,sampleid="sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
pathcom.rm.rel.lm<-pathway.compare(pathtab=kegg.rm,mapfile=covar.rm,sampleid="sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
pathcom.rm.rel.gamlss<-pathway.compare(pathtab=kegg.rm,mapfile=covar.rm,sampleid="sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
# age <=6months
pathcom.rm6.log<-pathway.compare(pathtab=kegg.rm,mapfile=covar.rm,sampleid="sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.rm6.rel.lm<-pathway.compare(pathtab=kegg.rm,mapfile=covar.rm,sampleid="sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.rm6.rel.gamlss<-pathway.compare(pathtab=kegg.rm,mapfile=covar.rm,sampleid="sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)

#save(pathcom.rm.log,pathcom.rm.rel.lm,pathcom.rm.rel.gamlss,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.rm.rda")
#save(pathcom.rm6.log,pathcom.rm6.rel.lm,pathcom.rm6.rel.gamlss,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.rm6.rda")
#save(pathdatrel,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathdatrel.rm.rda")

# gender adjusted for bf and age
pathcom.rm6.rel.gamlss.sex<-pathway.compare(pathtab=kegg.rm,mapfile=covar.rm,sampleid="sampleid",pathsum="rel",stat.med="gamlss",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.rm6.rel.gamlss.sex,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.rm6.rel.gamlss.sex.rda")
#generic function
pathcom.rm6.rel.gamlss.sexg<-pathway.compare.gen(pathtab=kegg.rm,mapfile=covar.rm,sampleid="sampleid",pathsum="rel",stat.med="gamlss",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.rm6.rel.gamlss.sexg,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.rm6.rel.gamlss.sexg.rda")


# Haiti
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Haiti_data_010417/open_uclust/metagenome_at_level3_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke3<-as.data.frame(t(kegg))
rownames(ke3)<-gsub("X","",rownames(ke3))
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Haiti_data_010417/open_uclust/metagenome_at_level2_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke2<-as.data.frame(t(kegg))
rownames(ke2)<-gsub("X","",rownames(ke2))
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Haiti_data_010417/open_uclust/metagenome_at_level1_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke1<-as.data.frame(t(kegg))
rownames(ke1)<-gsub("X","",rownames(ke1))
kegg.ha<-list(ke1,ke2,ke3)
haitimap<-read.delim('C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Haiti_data_010417/Haiti_Mapping.complete.txt')
colnames(haitimap)<-tolower(colnames(haitimap))
covar.ha<-subset(haitimap, mombb=="B"&sampletype=="STL")
covar.ha$bf<-mapvalues(covar.ha$exclusivebf,from=c("No","Yes"), to=c("Non_exclusiveBF","ExclusiveBF"))
covar.ha$bf<-factor(covar.ha$bf,levels=c("ExclusiveBF","Non_exclusiveBF"))
covar.ha$age.sample<-covar.ha$bbage/30
covar.ha$gender<-factor(covar.ha$babygender,levels=c("Female","Male"))

pathcom.ha6.log<-pathway.compare(pathtab=kegg.ha,mapfile=covar.ha,sampleid="x.sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.ha6.rel.lm<-pathway.compare(pathtab=kegg.ha,mapfile=covar.ha,sampleid="x.sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.ha6.rel.gamlss<-pathway.compare(pathtab=kegg.ha,mapfile=covar.ha,sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.ha6.log,pathcom.ha6.rel.lm,pathcom.ha6.rel.gamlss,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.ha6.rda")


#gender adjusted for bf and age
pathcom.ha6.rel.gamlss.sex<-pathway.compare(pathtab=kegg.ha,mapfile=covar.ha,sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.ha6.rel.gamlss.sex,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.ha6.rel.gamlss.sex.rda")
#generic function
pathcom.ha6.rel.gamlss.sexg<-pathway.compare.gen(pathtab=kegg.ha,mapfile=covar.ha,sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.ha6.rel.gamlss.sexg,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.ha6.rel.gamlss.sexg.rda")


# delivery antibiotics
pathcom.ha6.log.m<-pathway.compare(pathtab=kegg.ha,mapfile=covar.ha,sampleid="x.sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar=c("age.sample","delivery","momantibiotics"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.ha6.rel.lm.m<-pathway.compare(pathtab=kegg.ha,mapfile=covar.ha,sampleid="x.sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar=c("age.sample","delivery","momantibiotics"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.ha6.rel.gamlss.m<-pathway.compare(pathtab=kegg.ha,mapfile=covar.ha,sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar=c("age.sample","delivery","momantibiotics"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.ha6.log.m,pathcom.ha6.rel.lm.m,pathcom.ha6.rel.gamlss.m,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.ha6.m.rda")

#delivery adjusted
pathcom.ha6.log.bm<-pathway.compare(pathtab=kegg.ha,mapfile=covar.ha,sampleid="x.sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.ha6.rel.lm.bm<-pathway.compare(pathtab=kegg.ha,mapfile=covar.ha,sampleid="x.sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.ha6.rel.gamlss.bm<-pathway.compare(pathtab=kegg.ha,mapfile=covar.ha,sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.ha6.log.bm,pathcom.ha6.rel.lm.bm,pathcom.ha6.rel.gamlss.bm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.ha6.bm.rda")

#stratification delivery
#vaginal
pathcom.ha6.log.vag<-pathway.compare(pathtab=kegg.ha,mapfile=subset(covar.ha,delivery=="Vaginal"),sampleid="x.sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.ha6.rel.lm.vag<-pathway.compare(pathtab=kegg.ha,mapfile=subset(covar.ha,delivery=="Vaginal"),sampleid="x.sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.ha6.rel.gamlss.vag<-pathway.compare(pathtab=kegg.ha,mapfile=subset(covar.ha,delivery=="Vaginal"),sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
# C-section
pathcom.ha6.log.cs<-pathway.compare(pathtab=kegg.ha,mapfile=subset(covar.ha,delivery=="C-section"),sampleid="x.sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.ha6.rel.lm.cs<-pathway.compare(pathtab=kegg.ha,mapfile=subset(covar.ha,delivery=="C-section"),sampleid="x.sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.ha6.rel.gamlss.cs<-pathway.compare(pathtab=kegg.ha,mapfile=subset(covar.ha,delivery=="C-section"),sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.ha6.log.vag,pathcom.ha6.rel.lm.vag,pathcom.ha6.rel.gamlss.vag,pathcom.ha6.log.cs,pathcom.ha6.rel.lm.cs,pathcom.ha6.rel.gamlss.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.ha6.vagcs.rda")


# USA(CA_FL) USBMK
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/USBMK/metagenome_at_level3_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke3<-as.data.frame(t(kegg))
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/USBMK/metagenome_at_level2_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke2<-as.data.frame(t(kegg))
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/USBMK/metagenome_at_level1_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke1<-as.data.frame(t(kegg))
kegg.usbmk<-list(ke1,ke2,ke3)
usbmkmap<-read.delim('C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/USBMK/USBMK_Mapping_122414.filtered.with_pairing_status.txt')
colnames(usbmkmap)<-tolower(colnames(usbmkmap))
covar.usbmk<-subset(usbmkmap, mombb=="Baby"&sampletype=="STL")
covar.usbmk$age.sample<-covar.usbmk$bbage/30
covar.usbmk$bf<-factor(covar.usbmk$bfmixfm)
levels(covar.usbmk$bf)[levels(covar.usbmk$bf)==""]<-NA
covar.usbmk$bf<-mapvalues(covar.usbmk$bf,from=c("BF","FM","Mix","SOLID"),to=c("ExclusiveBF","No_BF","Non_exclusiveBF","No_BF"))
covar.usbmk$bf<-factor(covar.usbmk$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
covar.usbmk$personid<-covar.usbmk$pairid
covar.usbmk$birth.mode<-covar.usbmk$delivery
levels(covar.usbmk$birth.mode)[levels(covar.usbmk$birth.mode)==""]<-NA
covar.usbmk$gender<-mapvalues(covar.usbmk$babygender,from=c("F","M"),to=c("Female","Male"))
levels(covar.usbmk$gender)[levels(covar.usbmk$gender)==""]<-NA
covar.usbmk$gender<-factor(covar.usbmk$gender,levels=c("Female","Male"))

pathcom.usbmk6.log<-pathway.compare(pathtab=kegg.usbmk,mapfile=covar.usbmk,sampleid="x.sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.usbmk6.rel.lm<-pathway.compare(pathtab=kegg.usbmk,mapfile=covar.usbmk,sampleid="x.sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.usbmk6.rel.gamlss<-pathway.compare(pathtab=kegg.usbmk,mapfile=covar.usbmk,sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.usbmk6.log,pathcom.usbmk6.rel.lm,pathcom.usbmk6.rel.gamlss,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.usbmk6.rda")

#gender adjusted for bf age
#error if longitudinal => no longitudinal
pathcom.usbmk6.rel.gamlss.sex<-pathway.compare(pathtab=kegg.usbmk,mapfile=covar.usbmk,sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.usbmk6.rel.gamlss.sex,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.usbmk6.rel.gamlss.sex.rda")
#generic function
pathcom.usbmk6.rel.gamlss.sexg<-pathway.compare.gen(pathtab=kegg.usbmk,mapfile=covar.usbmk,sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.usbmk6.rel.gamlss.sexg,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.usbmk6.rel.gamlss.sexg.rda")


# delivery antibiotics adjusted
pathcom.usbmk6.log.m<-pathway.compare(pathtab=kegg.usbmk,mapfile=covar.usbmk,sampleid="x.sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar=c("age.sample","delivery","abxatdelivery"),longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.usbmk6.rel.lm.m<-pathway.compare(pathtab=kegg.usbmk,mapfile=covar.usbmk,sampleid="x.sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar=c("age.sample","delivery","abxatdelivery"),longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.usbmk6.rel.gamlss.m<-pathway.compare(pathtab=kegg.usbmk,mapfile=covar.usbmk,sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar=c("age.sample","delivery","abxatdelivery"),longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.usbmk6.log.m,pathcom.usbmk6.rel.lm.m,pathcom.usbmk6.rel.gamlss.m,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.usbmk6.m.rda")

# delivery adjusted (model error if code delivery "" as NA)
pathcom.usbmk6.log.bm<-pathway.compare(pathtab=kegg.usbmk,mapfile=covar.usbmk,sampleid="x.sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.usbmk6.rel.lm.bm<-pathway.compare(pathtab=kegg.usbmk,mapfile=covar.usbmk,sampleid="x.sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.usbmk6.rel.gamlss.bm<-pathway.compare(pathtab=kegg.usbmk,mapfile=covar.usbmk,sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar=c("age.sample","delivery"),longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.usbmk6.log.bm,pathcom.usbmk6.rel.lm.bm,pathcom.usbmk6.rel.gamlss.bm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.usbmk6.bm.rda")

# stratification birth.mode
# vaginal
pathcom.usbmk6.log.vag<-pathway.compare(pathtab=kegg.usbmk,mapfile=subset(covar.usbmk,delivery=="Vaginal"),sampleid="x.sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar=c("age.sample"),longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.usbmk6.rel.lm.vag<-pathway.compare(pathtab=kegg.usbmk,mapfile=subset(covar.usbmk,delivery=="Vaginal"),sampleid="x.sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar=c("age.sample"),longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.usbmk6.rel.gamlss.vag<-pathway.compare(pathtab=kegg.usbmk,mapfile=subset(covar.usbmk,delivery=="Vaginal"),sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar=c("age.sample"),longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#c-section
pathcom.usbmk6.log.cs<-pathway.compare(pathtab=kegg.usbmk,mapfile=subset(covar.usbmk,delivery=="C-section"),sampleid="x.sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar=c("age.sample"),longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.usbmk6.rel.lm.cs<-pathway.compare(pathtab=kegg.usbmk,mapfile=subset(covar.usbmk,delivery=="C-section"),sampleid="x.sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar=c("age.sample"),longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.usbmk6.rel.gamlss.cs<-pathway.compare(pathtab=kegg.usbmk,mapfile=subset(covar.usbmk,delivery=="C-section"),sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#cs model error with longitudinal so use non longitudinal instead.
#save(pathcom.usbmk6.log.vag,pathcom.usbmk6.rel.lm.vag,pathcom.usbmk6.rel.gamlss.vag,pathcom.usbmk6.log.cs,pathcom.usbmk6.rel.lm.cs,pathcom.usbmk6.rel.gamlss.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.usbmk6.vagcs.rda")

#USA(NC)
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Frontier/open_uclust2/metagenome_at_level3_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke3<-as.data.frame(t(kegg))
rownames(ke3)<-gsub("X","",rownames(ke3))
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Frontier/open_uclust2/metagenome_at_level2_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke2<-as.data.frame(t(kegg))
rownames(ke2)<-gsub("X","",rownames(ke2))
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Frontier/open_uclust2/metagenome_at_level1_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke1<-as.data.frame(t(kegg))
rownames(ke1)<-gsub("X","",rownames(ke1))
kegg.unc<-list(ke1,ke2,ke3)
uncmap<-read.delim('C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Frontier/mapping_fileall.txt')
colnames(uncmap)<-tolower(colnames(uncmap))
covar.unc<-uncmap
covar.unc$bf<-factor(covar.unc$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
covar.unc<-covar.unc %>% rename(personid=patient,age.sample=agemo)
covar.unc$personid<-as.factor(covar.unc$personid)
covar.unc$gender<-mapvalues(covar.unc$sex,from=c("1","0"),to=c("Female","Male"))
covar.unc$gender<-factor(covar.unc$gender,levels=c("Female","Male"))

pathcom.unc6.log<-pathway.compare(pathtab=kegg.unc,mapfile=covar.unc,sampleid="x.sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.unc6.rel.lm<-pathway.compare(pathtab=kegg.unc,mapfile=covar.unc,sampleid="x.sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.unc.rel.gamlss<-pathway.compare(pathtab=kegg.unc,mapfile=covar.unc,sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=100)
# error gamlss if <6month longitudinal => no longitudinal
pathcom.unc6.rel.gamlss<-pathway.compare(pathtab=kegg.unc,mapfile=covar.unc,sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.unc6.log,pathcom.unc6.rel.lm,pathcom.unc.rel.gamlss,pathcom.unc6.rel.gamlss,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.unc6.rda")

#gender adjusted for bf and age
pathcom.unc6.rel.gamlss.sex<-pathway.compare(pathtab=kegg.unc,mapfile=covar.unc,sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.unc6.rel.gamlss.sex,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.unc6.rel.gamlss.sex.rda")
#generic function
pathcom.unc6.rel.gamlss.sexg<-pathway.compare.gen(pathtab=kegg.unc,mapfile=covar.unc,sampleid="x.sampleid",pathsum="rel",stat.med="gamlss",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#save(pathcom.unc6.rel.gamlss.sexg,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.unc6.rel.gamlss.sexg.rda")


# South Africa UW
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UW/newrawdata/open_uclust/metagenome_at_level3_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke3<-as.data.frame(t(kegg))
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UW/newrawdata/open_uclust/metagenome_at_level2_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke2<-as.data.frame(t(kegg))
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UW/newrawdata/open_uclust/metagenome_at_level1_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke1<-as.data.frame(t(kegg))
kegg.uw<-list(ke1,ke2,ke3)
covar.uw<-read.delim('C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UW/HJ_feeding_mapping_BPB_4Nhan.txt')
colnames(covar.uw)<-tolower(colnames(covar.uw))
covar.uw$age.sample<-as.numeric(as.character(substr(covar.uw$timepoint,2,3)))/4
covar.uw$bf<-mapvalues(covar.uw$ebf,from=c("0","1"),to=c("Non_exclusiveBF","ExclusiveBF"))
covar.uw$bf<-factor(covar.uw$bf,levels=c("ExclusiveBF","Non_exclusiveBF"))
covar.uw$personid<-as.factor(covar.uw$personid)
covar.uw$sampleid<-paste("X",sub("\\.S.*$", "", covar.uw$x.sampleid),sep="")
pathcom.uw6.log<-pathway.compare(pathtab=kegg.uw,mapfile=covar.uw,sampleid="sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.uw6.rel.lm<-pathway.compare(pathtab=kegg.uw,mapfile=covar.uw,sampleid="sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
pathcom.uw6.rel.gamlss<-pathway.compare(pathtab=kegg.uw,mapfile=covar.uw,sampleid="sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
# error if longitudinal="yes" => no longitudinal
#save(pathcom.uw6.log,pathcom.uw6.rel.lm,pathcom.uw6.rel.gamlss,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.uw6.rda")

# Canada (UAlberta)
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UAlberta/predicted_metagenomes_BJOG.L3_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke3<-as.data.frame(t(kegg))
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UAlberta/predicted_metagenomes_BJOG.L2_edit.txt",header=TRUE)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke2<-as.data.frame(t(kegg))
kegg.alb<-list(ke2,ke3)
bm<-read.csv("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UAlberta/Birth_Mode_BJOG.csv")
ap<-read.csv("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/UAlberta/alphamean_BJOG.csv")
covar.alb<-merge(bm,ap[,c("sampleid","bf","age.sample")],by="sampleid")
covar.alb$sampleid<-paste("X",covar.alb$sampleid,sep="")
colnames(covar.alb)<-tolower(colnames(covar.alb))
covar.alb$age.sample<-as.numeric(as.character(covar.alb$age.sample))
covar.alb$bf<-factor(covar.alb$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
pathcom.alb6.log<-pathway.compare(pathtab=kegg.alb,mapfile=covar.alb,sampleid="sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.alb6.log)<-c("l2","l3")
pathcom.alb6.rel.lm<-pathway.compare(pathtab=kegg.alb,mapfile=covar.alb,sampleid="sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.alb6.rel.lm)<-c("l2","l3")
pathcom.alb6.rel.gamlss<-pathway.compare(pathtab=kegg.alb,mapfile=covar.alb,sampleid="sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.alb6.rel.gamlss)<-c("l2","l3")
#save(pathcom.alb6.log,pathcom.alb6.rel.lm,pathcom.alb6.rel.gamlss,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.alb6.rda")

# stratify by birth mode
# vaginal
pathcom.alb6.log.vag<-pathway.compare(pathtab=kegg.alb,mapfile=subset(covar.alb,birth.mode=="Vaginal"),sampleid="sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.alb6.log.vag)<-c("l2","l3")
pathcom.alb6.rel.lm.vag<-pathway.compare(pathtab=kegg.alb,mapfile=subset(covar.alb,birth.mode=="Vaginal"),sampleid="sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.alb6.rel.lm.vag)<-c("l2","l3")
pathcom.alb6.rel.gamlss.vag<-pathway.compare(pathtab=kegg.alb,mapfile=subset(covar.alb,birth.mode=="Vaginal"),sampleid="sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.alb6.rel.gamlss.vag)<-c("l2","l3")
# c-section
pathcom.alb6.log.cs<-pathway.compare(pathtab=kegg.alb,mapfile=subset(covar.alb,birth.mode!="Vaginal"),sampleid="sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.alb6.log.cs)<-c("l2","l3")
pathcom.alb6.rel.lm.cs<-pathway.compare(pathtab=kegg.alb,mapfile=subset(covar.alb,birth.mode!="Vaginal"),sampleid="sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.alb6.rel.lm.cs)<-c("l2","l3")
pathcom.alb6.rel.gamlss.cs<-pathway.compare(pathtab=kegg.alb,mapfile=subset(covar.alb,birth.mode!="Vaginal"),sampleid="sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.alb6.rel.gamlss.cs)<-c("l2","l3")
#save(pathcom.alb6.log.vag,pathcom.alb6.rel.lm.vag,pathcom.alb6.rel.gamlss.vag,pathcom.alb6.log.cs,pathcom.alb6.rel.lm.cs,pathcom.alb6.rel.gamlss.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.alb6.vagcs.rda")

# adjust for birth mode
pathcom.alb6.log.bm<-pathway.compare(pathtab=kegg.alb,mapfile=covar.alb,sampleid="sampleid",pathsum="log",stat.med="lm",comvar="bf",adjustvar=c("age.sample","birth.mode"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.alb6.log.bm)<-c("l2","l3")
pathcom.alb6.rel.lm.bm<-pathway.compare(pathtab=kegg.alb,mapfile=covar.alb,sampleid="sampleid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar=c("age.sample","birth.mode"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.alb6.rel.lm.bm)<-c("l2","l3")
pathcom.alb6.rel.gamlss.bm<-pathway.compare(pathtab=kegg.alb,mapfile=covar.alb,sampleid="sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar=c("age.sample","birth.mode"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.alb6.rel.gamlss.bm)<-c("l2","l3")
#save(pathcom.alb6.log.bm,pathcom.alb6.rel.lm.bm,pathcom.alb6.rel.gamlss.bm,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.alb6.bm.rda")


# USA(CA_MA_MO) (Harvard)
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Harvard/ForNhan/SentToNhan/PICRUSt_3_edit.txt",header=TRUE)
l3rm<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/open_uclust6/metagenome_at_level3_edit.txt",header=TRUE)
kegg<-merge(kegg,l3rm[,c("OTU.ID","KEGG_Pathways")],by="OTU.ID", all.x=T)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke3<-as.data.frame(t(kegg))
kegg<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Harvard/ForNhan/SentToNhan/PICRUSt_2_edit.txt",header=TRUE)
l2rm<-read.delim("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/relativematurity/open_uclust6/metagenome_at_level2_edit.txt",header=TRUE)
kegg<-merge(kegg,l2rm[,c("OTU.ID","KEGG_Pathways")],by="OTU.ID", all.x=T)
rownames(kegg)<-kegg[,"KEGG_Pathways"]
kegg<-kegg[,colnames(kegg)[!colnames(kegg) %in% c("OTU.ID","KEGG_Pathways")]]
ke2<-as.data.frame(t(kegg))
kegg.hav<-list(ke2,ke3)
covar.hav<-read.csv("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/data/Harvard/ForNhan/SentToNhan/final_map.csv")
colnames(covar.hav)<-tolower(colnames(covar.hav))
covar.hav$age.sample<-as.numeric(as.character(covar.hav$age_in_days))/30
covar.hav$bf<-mapvalues(covar.hav$breastfeeding,from=c("Excl_BF","Any_BF","No_BF"),to=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
covar.hav$bf<-factor(covar.hav$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
pathcom.hav6.log<-pathway.compare(pathtab=kegg.hav,mapfile=covar.hav,sampleid="anonid",pathsum="log",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.hav6.log)<-c("l2","l3")
pathcom.hav6.rel.lm<-pathway.compare(pathtab=kegg.hav,mapfile=covar.hav,sampleid="anonid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.hav6.rel.lm)<-c("l2","l3")
pathcom.hav6.rel.gamlss<-pathway.compare(pathtab=kegg.hav,mapfile=covar.hav,sampleid="anonid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.hav6.rel.gamlss)<-c("l2","l3")
#save(pathcom.hav6.log,pathcom.hav6.rel.lm,pathcom.hav6.rel.gamlss,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.hav6.rda")

#delivery adjusted
pathcom.hav6.log.bm<-pathway.compare(pathtab=kegg.hav,mapfile=covar.hav,sampleid="anonid",pathsum="log",stat.med="lm",comvar="bf",adjustvar=c("age.sample","csec"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.hav6.log.bm)<-c("l2","l3")
pathcom.hav6.rel.lm.bm<-pathway.compare(pathtab=kegg.hav,mapfile=covar.hav,sampleid="anonid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar=c("age.sample","csec"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.hav6.rel.lm.bm)<-c("l2","l3")
pathcom.hav6.rel.gamlss.bm<-pathway.compare(pathtab=kegg.hav,mapfile=covar.hav,sampleid="anonid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar=c("age.sample","csec"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.hav6.rel.gamlss.bm)<-c("l2","l3")

# stratification delivery
#vaginal
pathcom.hav6.log.vag<-pathway.compare(pathtab=kegg.hav,mapfile=subset(covar.hav,csec=="0"),sampleid="anonid",pathsum="log",stat.med="lm",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.hav6.log.vag)<-c("l2","l3")
pathcom.hav6.rel.lm.vag<-pathway.compare(pathtab=kegg.hav,mapfile=subset(covar.hav,csec=="0"),sampleid="anonid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.hav6.rel.lm.vag)<-c("l2","l3")
pathcom.hav6.rel.gamlss.vag<-pathway.compare(pathtab=kegg.hav,mapfile=subset(covar.hav,csec=="0"),sampleid="anonid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.hav6.rel.gamlss.vag)<-c("l2","l3")
#csection
pathcom.hav6.log.cs<-pathway.compare(pathtab=kegg.hav,mapfile=subset(covar.hav,csec=="1"),sampleid="anonid",pathsum="log",stat.med="lm",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.hav6.log.cs)<-c("l2","l3")
pathcom.hav6.rel.lm.cs<-pathway.compare(pathtab=kegg.hav,mapfile=subset(covar.hav,csec=="1"),sampleid="anonid",pathsum="rel",stat.med="lm",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.hav6.rel.lm.cs)<-c("l2","l3")
pathcom.hav6.rel.gamlss.cs<-pathway.compare(pathtab=kegg.hav,mapfile=subset(covar.hav,csec=="1"),sampleid="anonid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar=c("age.sample"),longitudinal="no",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)
names(pathcom.hav6.rel.gamlss.cs)<-c("l2","l3")
#save(pathcom.hav6.log.bm,pathcom.hav6.rel.lm.bm,pathcom.hav6.rel.gamlss.bm,pathcom.hav6.log.vag,pathcom.hav6.rel.lm.vag,pathcom.hav6.rel.gamlss.vag,
#     pathcom.hav6.log.cs,pathcom.hav6.rel.lm.cs,pathcom.hav6.rel.gamlss.cs,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.hav6.bmvagcs.rda")




# Meta analysis of pathways
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.unc6.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.usbmk6.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.ha6.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.rm6.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.uw6.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.alb6.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.hav6.rda"))

# GAMLSS
#RM
taxacom.zi.rm<-pathcom.rm6.rel.gamlss
for (i in 1: length(names(taxacom.zi.rm))){
  taxacom.zi.rm[[i]]<-as.data.frame(taxacom.zi.rm[[i]])
  taxacom.zi.rm[[i]][,'path']<-rownames(taxacom.zi.rm[[i]])
  taxacom.zi.rm[[i]][,'author']<-"Subramanian et al"
  taxacom.zi.rm[[i]][,'year']<-"2014"
  taxacom.zi.rm[[i]][,'pop']<-"Bangladesh"
}
#Haiti
taxacom.zi.ha<-pathcom.ha6.rel.gamlss
for (i in 1: length(names(taxacom.zi.ha))){
  taxacom.zi.ha[[i]]<-as.data.frame(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'path']<-rownames(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'author']<-"Bender et al"
  taxacom.zi.ha[[i]][,'year']<-"2016"
  taxacom.zi.ha[[i]][,'pop']<-"Haiti"
}
# UNC
#fix names
names(pathcom.unc.rel.gamlss)<-c("l1","l2","l3")
taxacom.zi.unc<-pathcom.unc.rel.gamlss
for (i in 1:length(names(taxacom.zi.unc))){ #
  taxacom.zi.unc[[i]]<-as.data.frame(taxacom.zi.unc[[i]])
  taxacom.zi.unc[[i]][,'path']<-rownames(taxacom.zi.unc[[i]])
  taxacom.zi.unc[[i]][,'author']<-"Thompson et al"
  taxacom.zi.unc[[i]][,'year']<-"2015"
  taxacom.zi.unc[[i]][,'pop']<-"USA(NC)"
}
#USBMK
taxacom.zi.usbmk<-pathcom.usbmk6.rel.gamlss
for (i in 1: length(names(taxacom.zi.usbmk))){
  taxacom.zi.usbmk[[i]]<-as.data.frame(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'path']<-rownames(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'author']<-"Pannaraj et al"
  taxacom.zi.usbmk[[i]][,'year']<-"2017"
  taxacom.zi.usbmk[[i]][,'pop']<-"USA(CA_FL)"
}
#UW
taxacom.zi.uw<-pathcom.uw6.rel.gamlss
for (i in 1: length(names(taxacom.zi.uw))){
  taxacom.zi.uw[[i]]<-as.data.frame(taxacom.zi.uw[[i]])
  taxacom.zi.uw[[i]][,'path']<-rownames(taxacom.zi.uw[[i]])
  taxacom.zi.uw[[i]][,'author']<-"Wood et al"
  taxacom.zi.uw[[i]][,'year']<-"2017"
  taxacom.zi.uw[[i]][,'pop']<-"South Africa"
}
#Ualberta
taxacom.zi.alb<-pathcom.alb6.rel.gamlss
for (i in 1: length(names(taxacom.zi.alb))){
  taxacom.zi.alb[[i]]<-as.data.frame(taxacom.zi.alb[[i]])
  taxacom.zi.alb[[i]][,'path']<-rownames(taxacom.zi.alb[[i]])
  taxacom.zi.alb[[i]][,'author']<-"Azad et al"
  taxacom.zi.alb[[i]][,'year']<-"2015"
  taxacom.zi.alb[[i]][,'pop']<-"Canada"
}
#harvard
taxacom.zi.hav<-pathcom.hav6.rel.gamlss
for (i in 1: length(names(taxacom.zi.hav))){
  taxacom.zi.hav[[i]]<-as.data.frame(taxacom.zi.hav[[i]])
  taxacom.zi.hav[[i]][,'path']<-rownames(taxacom.zi.hav[[i]])
  taxacom.zi.hav[[i]][,'author']<-"Sordillo et al"
  taxacom.zi.hav[[i]][,'year']<-"2017"
  taxacom.zi.hav[[i]][,'pop']<-"USA(CA_MA_MO)"
}

#taxacom.zi.l1<-rbind.fill(taxacom.zi.rm$l1,taxacom.zi.ha$l1,taxacom.zi.unc$l1,taxacom.zi.usbmk$l1,taxacom.zi.uw$l1)
taxacom.zi.l2<-rbind.fill(taxacom.zi.rm$l2,taxacom.zi.ha$l2,taxacom.zi.unc$l2,taxacom.zi.usbmk$l2,taxacom.zi.uw$l2,taxacom.zi.alb$l2,taxacom.zi.hav$l2)
taxacom.zi.l3<-rbind.fill(taxacom.zi.rm$l3,taxacom.zi.ha$l3,taxacom.zi.unc$l3,taxacom.zi.usbmk$l3,taxacom.zi.uw$l3,taxacom.zi.alb$l3,taxacom.zi.hav$l3)
pathcom.zi<-list(taxacom.zi.l2,taxacom.zi.l3)
names(pathcom.zi)<-paste("l",2:3,sep="")
# metaanalysis for gamlss
pathmetatab.zi<-meta.taxa(taxcomdat=pathcom.zi, sm="RR",p.adjust.method="fdr",percent.meta=0.5,pool.var="path")
#save(pathcom.zi,pathmetatab.zi,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathmetatab7.rda")

# Sensitivity analysis
# No UNC data
taxacom.zi.l2.nounc<-rbind.fill(taxacom.zi.rm$l2,taxacom.zi.ha$l2,taxacom.zi.usbmk$l2,taxacom.zi.uw$l2,taxacom.zi.alb$l2,taxacom.zi.hav$l2)
taxacom.zi.l3.nounc<-rbind.fill(taxacom.zi.rm$l3,taxacom.zi.ha$l3,taxacom.zi.usbmk$l3,taxacom.zi.uw$l3,taxacom.zi.alb$l3,taxacom.zi.hav$l3)
pathcom.zi.nounc<-list(taxacom.zi.l2.nounc,taxacom.zi.l3.nounc)
names(pathcom.zi.nounc)<-paste("l",2:3,sep="")
# metaanalysis for gamlss
pathmetatab.zi.nounc<-meta.taxa(taxcomdat=pathcom.zi.nounc, sm="RR",p.adjust.method="fdr",percent.meta=0.5,pool.var="path")
#No Haiti data
taxacom.zi.l2.noha<-rbind.fill(taxacom.zi.rm$l2,taxacom.zi.unc$l2,taxacom.zi.usbmk$l2,taxacom.zi.uw$l2,taxacom.zi.alb$l2,taxacom.zi.hav$l2)
taxacom.zi.l3.noha<-rbind.fill(taxacom.zi.rm$l3,taxacom.zi.unc$l3,taxacom.zi.usbmk$l3,taxacom.zi.uw$l3,taxacom.zi.alb$l3,taxacom.zi.hav$l3)
pathcom.zi.noha<-list(taxacom.zi.l2.noha,taxacom.zi.l3.noha)
names(pathcom.zi.noha)<-paste("l",2:3,sep="")
# metaanalysis for gamlss
pathmetatab.zi.noha<-meta.taxa(taxcomdat=pathcom.zi.noha, sm="RR",p.adjust.method="fdr",percent.meta=0.5,pool.var="path")
# No harvard data
taxacom.zi.l2.nohav<-rbind.fill(taxacom.zi.rm$l2,taxacom.zi.ha$l2,taxacom.zi.unc$l2,taxacom.zi.usbmk$l2,taxacom.zi.uw$l2,taxacom.zi.alb$l2)
taxacom.zi.l3.nohav<-rbind.fill(taxacom.zi.rm$l3,taxacom.zi.ha$l3,taxacom.zi.unc$l3,taxacom.zi.usbmk$l3,taxacom.zi.uw$l3,taxacom.zi.alb$l3)
pathcom.zi.nohav<-list(taxacom.zi.l2.nohav,taxacom.zi.l3.nohav)
names(pathcom.zi.nohav)<-paste("l",2:3,sep="")
# metaanalysis for gamlss
pathmetatab.zi.nohav<-meta.taxa(taxcomdat=pathcom.zi.nohav, sm="RR",p.adjust.method="fdr",percent.meta=0.5,pool.var="path")
#save(pathcom.zi.nounc,pathmetatab.zi.nounc,pathcom.zi.noha,pathmetatab.zi.noha,pathcom.zi.nohav,pathmetatab.zi.nohav,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathmetatab.zi.sen.rda")


# Meta-analysis for 4 studies with gender info but adjusting only for age (not gender): for comparison
taxacom.zi.l2<-rbind.fill(taxacom.zi.rm$l2,taxacom.zi.ha$l2,taxacom.zi.unc$l2,taxacom.zi.usbmk$l2)
taxacom.zi.l3<-rbind.fill(taxacom.zi.rm$l3,taxacom.zi.ha$l3,taxacom.zi.unc$l3,taxacom.zi.usbmk$l3)
pathcom.zi.4<-list(taxacom.zi.l2,taxacom.zi.l3)
names(pathcom.zi.4)<-paste("l",2:3,sep="")
# metaanalysis for gamlss
pathmetatab.zi.4<-meta.taxa(taxcomdat=pathcom.zi.4, sm="RR",p.adjust.method="fdr",percent.meta=0.5,pool.var="path")
#save(pathcom.zi.4,pathmetatab.zi.4,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathmetatab.4.rda")



#gender adjusted for bf and age (generic function)
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.unc6.rel.gamlss.sexg.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.ha6.rel.gamlss.sexg.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.rm6.rel.gamlss.sexg.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.usbmk6.rel.gamlss.sexg.rda"))
#names(pathcom.rm6.rel.gamlss.sex)<-names(pathcom.ha6.rel.gamlss.sex)<-names(pathcom.unc6.rel.gamlss.sex)<-names(pathcom.usbmk6.rel.gamlss.sex)<-c("l1","l2","l3")
#RM
taxacom.zi.rm<-pathcom.rm6.rel.gamlss.sexg
for (i in 1: length(names(taxacom.zi.rm))){
  taxacom.zi.rm[[i]]<-as.data.frame(taxacom.zi.rm[[i]])
  taxacom.zi.rm[[i]][,'path']<-rownames(taxacom.zi.rm[[i]])
  taxacom.zi.rm[[i]][,'study']<-"Subramanian et al 2014 (Bangladesh)"
  taxacom.zi.rm[[i]][,'pop']<-"Bangladesh"
}
#Haiti
taxacom.zi.ha<-pathcom.ha6.rel.gamlss.sexg
for (i in 1: length(names(taxacom.zi.ha))){
  taxacom.zi.ha[[i]]<-as.data.frame(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'path']<-rownames(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'study']<-"Bender et al 2016 (Haiti)"
  taxacom.zi.ha[[i]][,'pop']<-"Haiti"
}
# UNC
taxacom.zi.unc<-pathcom.unc6.rel.gamlss.sexg
for (i in 1:length(names(taxacom.zi.unc))){ #
  taxacom.zi.unc[[i]]<-as.data.frame(taxacom.zi.unc[[i]])
  taxacom.zi.unc[[i]][,'path']<-rownames(taxacom.zi.unc[[i]])
  taxacom.zi.unc[[i]][,'study']<-"Thompson et al 2015 (USA(NC))"
  taxacom.zi.unc[[i]][,'pop']<-"USA(NC)"
}
#USBMK
taxacom.zi.usbmk<-pathcom.usbmk6.rel.gamlss.sexg
for (i in 1: length(names(taxacom.zi.usbmk))){
  taxacom.zi.usbmk[[i]]<-as.data.frame(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'path']<-rownames(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'study']<-"Pannaraj et al 2017 (USA(CA_FL))"
  taxacom.zi.usbmk[[i]][,'pop']<-"USA(CA_FL)"
}
taxacom.zi.l2<-rbind.fill(taxacom.zi.rm$l2,taxacom.zi.ha$l2,taxacom.zi.unc$l2,taxacom.zi.usbmk$l2)
taxacom.zi.l3<-rbind.fill(taxacom.zi.rm$l3,taxacom.zi.ha$l3,taxacom.zi.unc$l3,taxacom.zi.usbmk$l3)
pathcom.zi.sexg<-list(taxacom.zi.l2,taxacom.zi.l3)
names(pathcom.zi.sexg)<-paste("l",2:3,sep="")
# metaanalysis for gamlss
pathmetatab.zi.sex<-meta.taxa(taxcomdat=pathcom.zi.sex, sm="RR",p.adjust.method="fdr",percent.meta=0.5,pool.var="path")
#save(pathmetatab.zi.sex,pathcom.zi.sex,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathmetatab.zi.sex.rda")
#generic function
pathmetatab.zi.sex.l2<-meta.taxa.gen(taxcomdat=pathcom.zi.sexg$l2, sm="RR",studylab = "pop", p.adjust.method="fdr",percent.meta=0.5,pool.var="id")
pathmetatab.zi.sex.l3<-meta.taxa.gen(taxcomdat=pathcom.zi.sexg$l3, sm="RR",studylab = "pop", p.adjust.method="fdr",percent.meta=0.5,pool.var="id")
#save(pathmetatab.zi.sex.l2,pathmetatab.zi.sex.l3,pathcom.zi.sexg,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathmetatab.zi.sexg.rda")


# stratification by birth mode
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.ha6.vagcs.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.hav6.bmvagcs.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.alb6.vagcs.rda"))
print(load("C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathcom.usbmk6.vagcs.rda"))

# vaginal
# GAMLSS
#Haiti
taxacom.zi.ha<-pathcom.ha6.rel.gamlss.vag
for (i in 1: length(names(taxacom.zi.ha))){
  taxacom.zi.ha[[i]]<-as.data.frame(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'path']<-rownames(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'author']<-"Bender et al"
  taxacom.zi.ha[[i]][,'year']<-"2016"
  taxacom.zi.ha[[i]][,'pop']<-"Haiti"
}
#USBMK
taxacom.zi.usbmk<-pathcom.usbmk6.rel.gamlss.vag
for (i in 1: length(names(taxacom.zi.usbmk))){
  taxacom.zi.usbmk[[i]]<-as.data.frame(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'path']<-rownames(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'author']<-"Pannaraj et al"
  taxacom.zi.usbmk[[i]][,'year']<-"2017"
  taxacom.zi.usbmk[[i]][,'pop']<-"USA(CA_FL)"
}
#Ualberta
taxacom.zi.alb<-pathcom.alb6.rel.gamlss.vag
for (i in 1: length(names(taxacom.zi.alb))){
  taxacom.zi.alb[[i]]<-as.data.frame(taxacom.zi.alb[[i]])
  taxacom.zi.alb[[i]][,'path']<-rownames(taxacom.zi.alb[[i]])
  taxacom.zi.alb[[i]][,'author']<-"Azad et al"
  taxacom.zi.alb[[i]][,'year']<-"2015"
  taxacom.zi.alb[[i]][,'pop']<-"Canada"
}
#harvard
taxacom.zi.hav<-pathcom.hav6.rel.gamlss.vag
for (i in 1: length(names(taxacom.zi.hav))){
  taxacom.zi.hav[[i]]<-as.data.frame(taxacom.zi.hav[[i]])
  taxacom.zi.hav[[i]][,'path']<-rownames(taxacom.zi.hav[[i]])
  taxacom.zi.hav[[i]][,'author']<-"Sordillo et al"
  taxacom.zi.hav[[i]][,'year']<-"2017"
  taxacom.zi.hav[[i]][,'pop']<-"USA(CA_MA_MO)"
}
taxacom.zi.l2.vag<-rbind.fill(taxacom.zi.ha$l2,taxacom.zi.usbmk$l2,taxacom.zi.alb$l2,taxacom.zi.hav$l2)
taxacom.zi.l3.vag<-rbind.fill(taxacom.zi.ha$l3,taxacom.zi.usbmk$l3,taxacom.zi.alb$l3,taxacom.zi.hav$l3)
pathcom.zi.vag<-list(taxacom.zi.l2.vag,taxacom.zi.l3.vag)
names(pathcom.zi.vag)<-paste("l",2:3,sep="")
# metaanalysis for gamlss
pathmetatab.zi.vag<-meta.taxa(taxcomdat=pathcom.zi.vag, sm="RR",p.adjust.method="fdr",percent.meta=0.5,pool.var="path")

# C-section
# GAMLSS
#Haiti: sample too small for C-section => not included
taxacom.zi.ha<-pathcom.ha6.rel.gamlss.cs
for (i in 1: length(names(taxacom.zi.ha))){
  taxacom.zi.ha[[i]]<-as.data.frame(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'path']<-rownames(taxacom.zi.ha[[i]])
  taxacom.zi.ha[[i]][,'author']<-"Bender et al"
  taxacom.zi.ha[[i]][,'year']<-"2016"
  taxacom.zi.ha[[i]][,'pop']<-"Haiti"
}

#USBMK
taxacom.zi.usbmk<-pathcom.usbmk6.rel.gamlss.cs
for (i in 1: length(names(taxacom.zi.usbmk))){
  taxacom.zi.usbmk[[i]]<-as.data.frame(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'path']<-rownames(taxacom.zi.usbmk[[i]])
  taxacom.zi.usbmk[[i]][,'author']<-"Pannaraj et al"
  taxacom.zi.usbmk[[i]][,'year']<-"2017"
  taxacom.zi.usbmk[[i]][,'pop']<-"USA(CA_FL)"
}
#Ualberta
taxacom.zi.alb<-pathcom.alb6.rel.gamlss.cs
for (i in 1: length(names(taxacom.zi.alb))){
  taxacom.zi.alb[[i]]<-as.data.frame(taxacom.zi.alb[[i]])
  taxacom.zi.alb[[i]][,'path']<-rownames(taxacom.zi.alb[[i]])
  taxacom.zi.alb[[i]][,'author']<-"Azad et al"
  taxacom.zi.alb[[i]][,'year']<-"2015"
  taxacom.zi.alb[[i]][,'pop']<-"Canada"
}
#harvard
taxacom.zi.hav<-pathcom.hav6.rel.gamlss.cs
for (i in 1: length(names(taxacom.zi.hav))){
  taxacom.zi.hav[[i]]<-as.data.frame(taxacom.zi.hav[[i]])
  taxacom.zi.hav[[i]][,'path']<-rownames(taxacom.zi.hav[[i]])
  taxacom.zi.hav[[i]][,'author']<-"Sordillo et al"
  taxacom.zi.hav[[i]][,'year']<-"2017"
  taxacom.zi.hav[[i]][,'pop']<-"USA(CA_MA_MO)"
}
taxacom.zi.l2.cs<-rbind.fill(taxacom.zi.ha$l2,taxacom.zi.usbmk$l2,taxacom.zi.alb$l2,taxacom.zi.hav$l2)
taxacom.zi.l3.cs<-rbind.fill(taxacom.zi.ha$l3,taxacom.zi.usbmk$l3,taxacom.zi.alb$l3,taxacom.zi.hav$l3)
pathcom.zi.cs<-list(taxacom.zi.l2.cs,taxacom.zi.l3.cs)
names(pathcom.zi.cs)<-paste("l",2:3,sep="")
# metaanalysis for gamlss
pathmetatab.zi.cs<-meta.taxa(taxcomdat=pathcom.zi.cs, sm="RR",p.adjust.method="fdr",percent.meta=0.5,pool.var="path")
#save(pathcom.zi.cs,pathmetatab.zi.cs,pathcom.zi.vag,pathmetatab.zi.vag,file="C:/Users/nth2111/My files/Dr Kuhn/Microbiome/Rprac/analysis/data/pathmetatab.zi.vagcs.rda")




