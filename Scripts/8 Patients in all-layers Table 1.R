library(Seurat)
library(dplyr)
library(tidyr)
library(foreign)
library(survival)
library(survminer)
library(table1)
library(ggpubr)

##Table 1
###Table1 for patients in Metabolome
pdata.trauma<-Pamper.pdata %>% filter(PAMPID %in% PAMPer.Metabolome$PAID)
pdata.trauma$TREATMENT<-PAMPer.imp$TREATMENT[match(pdata.trauma$PAMPID,PAMPer.imp$PAMPER.ID.NUMBER)] %>% as.character()
pdata.trauma$outcome<-as.character(pdata.trauma$outcome)
pdata.trauma$traumatic_brain_injury<-pdata.trauma$traumatic_brain_injury %>% droplevels()
pdata.trauma$race<-as.character(pdata.trauma$race)
pdata.trauma$race[is.na(pdata.trauma$race)]<-"Non-white"
pdata.trauma$ed_coagulopathy[which(pdata.trauma$ed_coagulopathy==1)]<-"Yes"
pdata.trauma$ed_coagulopathy[which(pdata.trauma$ed_coagulopathy==2)]<-"No"
pdata.trauma$plasma_id_1<-as.character(pdata.trauma$plasma_id_1)
pdata.trauma$plasma_id_2<-as.character(pdata.trauma$plasma_id_2)
pdata.trauma$plasma_unit<-"2"
pdata.trauma$plasma_unit[which(pdata.trauma$plasma_id_1!="                " & pdata.trauma$plasma_id_2=="                 ")]<-"1"
pdata.trauma$plasma_unit[which(pdata.trauma$plasma_id_1=="                " & pdata.trauma$plasma_id_2=="                 ")]<-"0"
pdata.trauma$stnum[which(pdata.trauma$plasma_unit=="0" & pdata.trauma$TBI_Arms=="NonTBI_FFP")]

pdata.trauma$TREATMENT<-factor(pdata.trauma$TREATMENT, 
levels=c("standard care","prehospital plasma","p-value"))

pdata.trauma$TBI<-Layer5_merge$TRAUMATIC.BRAIN.INJURY[match(pdata.trauma$PAID,Layer5_merge$PAID)]

pdata.trauma$TBI_Arms<-paste(pdata.trauma$TBI,pdata.trauma$TREATMENT,sep = "_") %>% factor(levels = c("no_standard care","no_prehospital plasma","yes_standard care","yes_prehospital plasma","p-value"),labels = c("NonTBI_Std","NonTBI_FFP","TBI_Std","TBI_FFP","p-value"))


table(pdata.trauma$TBI_Arms)
pdata.trauma$outcome<-factor(pdata.trauma$outcome, 
                             levels=c("Resolving","Non-resolving","Early-Nonsurvivors","p-value"))


rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- pdata.trauma[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- kruskal.test(y ~ pdata.trauma$outcome)$p.value
    } else {
      p <- chisq.test(table(y, droplevels(pdata.trauma$outcome)))$p.value
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "MEDIAN (IQR)"=sprintf("%s (&plusmn; %s)", MEDIAN, IQR)))
}

my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%0.0f %%)", FREQ, PCT))))
}



table1(~ age+gender+race+iss+ais_head+initial_GCS+plasma_unit+PH_sbp_70+vitals_hr+injury_type+PH_time
       +PH_CPR+PH_intubation+PH_blood+PH_crystalloid+PH_prbc
       +transfusion_24h+prbc_24h+plasma_24h+platelets_24h+crystalloid_24h+vaso_24h+INR
       +ed_coagulopathy+ALI+NI+MOF+mech_vent_days+icu_los+hospital_los| outcome,
       data=pdata.trauma, droplevels=F, render=rndr, render.strat=rndr.strat, overall=F,
       render.continuous=my.render.cont, render.categorical="FREQ (PCTnoNA%)",topclass="Rtable1-zebra",render.missing=NULL)

###Table 1 for Proteome data
table(pdata.trauma.pro$severe_head)
table(pdata.trauma$severe_head)

pdata.trauma.pro<-Pamper.pdata[which(Pamper.pdata$PAID %in% PAMPer.Proteome$PAID),]
pdata.trauma.pro$TREATMENT<-factor(pdata.trauma.pro$plasma_on_helicopter,levels = c(1,2,"p-value"),labels = c("Plasma","Standard","p-value"))
pdata.trauma.pro$TBI<-PAMPer.Proteome$TBI[match(pdata.trauma.pro$PAID,PAMPer.Proteome$PAID)] %>% factor(levels = c("NO","YES","p-value"),labels = c("NonTBI","TBI","p-value"))
pdata.trauma.pro$TBI_Arms<-paste(pdata.trauma.pro$TBI,pdata.trauma.pro$TREATMENT,sep = "_")
table(pdata.trauma.pro$TBI_Arms)
pdata.trauma.pro$TBI_Arms<-factor(pdata.trauma.pro$TBI_Arms,levels = c("NonTBI_Standard","NonTBI_Plasma","TBI_Standard","TBI_Plasma","p-value"))

pdata.trauma.pro$plasma_id_1<-as.character(pdata.trauma.pro$plasma_id_1)
pdata.trauma.pro$plasma_id_2<-as.character(pdata.trauma.pro$plasma_id_2)
pdata.trauma.pro$plasma_unit<-"2"
pdata.trauma.pro$plasma_unit[which(pdata.trauma.pro$plasma_id_1!="                " & pdata.trauma.pro$plasma_id_2=="                 ")]<-"1"
pdata.trauma.pro$plasma_unit[which(pdata.trauma.pro$plasma_id_1=="                " & pdata.trauma.pro$plasma_id_2=="                 ")]<-"0"
pdata.trauma.pro$outcome<-pdata.trauma.pro$outcome %>% drop.levels() %>% factor(levels=c("Resolving","Non-resolving","p-value"))


#pdata.trauma.pro$outcome<-factor(pdata.trauma.pro$outcome, 
#levels=c("Resolving","Non-resolving","Early-Nonsurvivors","p-value"))
#pdata.trauma.pro$ed_coagulopathy[which(pdata.trauma.pro$ed_coagulopathy==1)]<-"Yes"
#pdata.trauma.pro$ed_coagulopathy[which(pdata.trauma.pro$ed_coagulopathy==2)]<-"No"


rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- pdata.trauma.pro[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- kruskal.test(y ~ pdata.trauma.pro$outcome)$p.value
    } else {
      p <- chisq.test(table(y, droplevels(pdata.trauma.pro$outcome)))$p.value
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "MEDIAN (IQR)"=sprintf("%s (&plusmn; %s)", MEDIAN, IQR)))
}

my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%0.0f %%)", FREQ, PCT))))
}



table1(~ age+gender+race+iss+ais_head+plasma_unit+initial_GCS+PH_sbp_70+vitals_hr+injury_type+TREATMENT+PH_time
       +PH_CPR+PH_intubation+PH_blood+PH_crystalloid+PH_prbc
       +transfusion_24h+prbc_24h+plasma_24h+platelets_24h+crystalloid_24h+vaso_24h+INR
       +ed_coagulopathy+ALI+NI+MOF+mech_vent_days+icu_los+hospital_los| outcome,
       data=pdata.trauma.pro, droplevels=F, render=rndr, render.strat=rndr.strat, overall=F,
       render.continuous=my.render.cont, render.categorical="FREQ (PCTnoNA%)",topclass="Rtable1-zebra",render.missing=NULL)

###outcome2 for metabolome

pdata.trauma$outcome2<-pdata.trauma$outcome %>% as.character()
pdata.trauma$outcome2[which(pdata.trauma$outcome2== "Non-resolving" & pdata.trauma$alive_at_30==2)]<-"Late-Nonsurvivors"

pdata.trauma$outcome2<-factor(pdata.trauma$outcome2, 
                             levels=c("Resolving","Non-resolving","Late-Nonsurvivors","Early-Nonsurvivors","p-value"))

rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- pdata.trauma[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- kruskal.test(y ~ pdata.trauma$outcome2)$p.value
    } else {
      p <- chisq.test(table(y, droplevels(pdata.trauma$outcome2)))$p.value
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "MEDIAN (IQR)"=sprintf("%s (&plusmn; %s)", MEDIAN, IQR)))
}

my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%0.0f %%)", FREQ, PCT))))
}



table1(~ age+gender+race+iss+ais_head+traumatic_brain_injury+initial_GCS+PH_sbp_70+vitals_hr+injury_type+TREATMENT+PH_time
       +PH_CPR+PH_intubation+PH_blood+PH_crystalloid+PH_prbc
       +transfusion_24h+prbc_24h+plasma_24h+platelets_24h+crystalloid_24h+vaso_24h+INR
       +ed_coagulopathy+ALI+NI+MOF+mech_vent_days+icu_los+hospital_los| outcome2,
       data=pdata.trauma, droplevels=F, render=rndr, render.strat=rndr.strat, overall=F,
       render.continuous=my.render.cont, render.categorical="FREQ (PCTnoNA%)",topclass="Rtable1-zebra",render.missing=NULL)

###outcome2 for whole pamper
Pamper.pdata$outcome2<-Pamper.pdata$outcome %>% as.character()
Pamper.pdata$outcome2[which(Pamper.pdata$outcome2== "Non-resolving" & Pamper.pdata$alive_at_30==2)]<-"Late-Nonsurvivors"

Pamper.pdata$outcome2<-factor(Pamper.pdata$outcome2, 
                              levels=c("Resolving","Non-resolving","Late-Nonsurvivors","Early-Nonsurvivors","p-value"))
Pamper.pdata$TREATMENT<-factor(Pamper.pdata$plasma_on_helicopter,levels = c(1,2),labels = c("Standard","Plasma"))
Pamper.pdata$coagulopathy<-factor(Pamper.pdata$ed_coagulopathy,levels = c(1,2),labels = c("Yes","No"))

table(Pamper.pdata$TREATMENT,Pamper.pdata$outcome2)

rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- Pamper.pdata[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- kruskal.test(y ~ Pamper.pdata$outcome2)$p.value
    } else {
      p <- chisq.test(table(y, droplevels(Pamper.pdata$outcome2)))$p.value
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "MEDIAN (IQR)"=sprintf("%s (&plusmn; %s)", MEDIAN, IQR)))
}

my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%0.0f %%)", FREQ, PCT))))
}



table1(~ age+gender+race+iss+ais_head+initial_GCS+PH_sbp_70+vitals_hr+injury_type+TREATMENT+PH_time
       +PH_CPR+PH_intubation+PH_blood+PH_crystalloid+PH_prbc
       +transfusion_24h+prbc_24h+plasma_24h+platelets_24h+crystalloid_24h+vaso_24h+INR
       +coagulopathy+ALI+NI+MOF+mech_vent_days+icu_los+hospital_los| outcome2,
       data=Pamper.pdata, droplevels=F, render=rndr, render.strat=rndr.strat, overall=F,
       render.continuous=my.render.cont, render.categorical="FREQ (PCTnoNA%)",topclass="Rtable1-zebra",render.missing=NULL)

###survival analysis for treatment effect(Metabolome)
###All patients
####K-P CURVE
surv_object <- Surv(time = pdata.trauma$t_30d_mort_h,event = pdata.trauma$t_30d_censor)
fit1 <- survfit(surv_object ~ TREATMENT, data = pdata.trauma)
ggsurvplot(fit1, data = pdata.trauma,  risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),pval = TRUE,
           xscale=24,xlim=c(0,720),break.time.by=24,censor=T)


###COX REGRESSION
fit.coxph <- coxph(surv_object ~age+gender+iss+TREATMENT, 
                   data = pdata.trauma)
ggforest(fit.coxph,data = pdata.trauma)

###All patients exclude early-nonsurvivors
####K-P CURVE
pdata.trauma.latephase<-pdata.trauma[which(pdata.trauma$outcome != "Early-Nonsurvivors"),] %>% drop.levels()
surv_object <- Surv(time = pdata.trauma.latephase$trans_icu_los,event = pdata.trauma.latephase$trans_icu_los_censor)
fit1 <- survfit(surv_object ~ TREATMENT, data = pdata.trauma.latephase)
dim(pdata.trauma.latephase)
length(pdata.trauma.latephase$TREATMENT)

ggsurvplot(fit1, data = pdata.trauma.latephase,  risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),pval = TRUE,
           xscale=1,xlim=c(0,30),break.time.by=4,censor=T)

###COX REGRESSION
fit.coxph <- coxph(surv_object ~age+gender+iss+TREATMENT, 
                   data = pdata.trauma.latephase)
ggforest(fit.coxph,data = pdata.trauma)


###Comparasion of different prehospital variables among outcome group
Pamper.pdata$outcome<-Pamper.pdata$outcome %>% factor(levels = c("Resolving","Non-resolving","Early-Nonsurvivors"))
my_comparisons <- list( c("Resolving", "Non-resolving"), c("Resolving", "Early-Nonsurvivors"), c("Non-resolving", "Early-Nonsurvivors") )
p1<-ggboxplot(Pamper.pdata[which(Pamper.pdata$Metabolome==1),], x = "outcome", y = "PH_time",
              color = "outcome")+stat_compare_means(label.y = 200)+scale_y_continuous(breaks=seq(0, 240, 60))+theme(aspect.ratio = 1)

p2<-ggboxplot(Pamper.pdata[which(Pamper.pdata$Proteome==1),], x = "outcome", y = "PH_time",
              color = "outcome")+stat_compare_means(label.y = 200)+scale_y_continuous(breaks=seq(0, 240, 60))+theme(aspect.ratio = 1)


p3<-ggviolin(Pamper.pdata[which(Pamper.pdata$Metabolome==1),], x = "outcome", y = "PH_prbc",
             color = "outcome",add = "median")+stat_compare_means(label.y = 2)+theme(aspect.ratio = 1)+scale_y_continuous(breaks=seq(0, 1, 5))+ylim(0,5)

p4<-ggviolin(Pamper.pdata[which(Pamper.pdata$Proteome==1),], x = "outcome", y = "PH_prbc",
             color = "outcome",add = "median")+stat_compare_means(label.y = 2)+theme(aspect.ratio = 1)+scale_y_continuous(breaks=seq(0, 1, 5))+ylim(0,5)

p5<-ggviolin(Pamper.pdata[which(Pamper.pdata$Metabolome==1),], x = "outcome", y = "PH_crystalloid",
             color = "outcome",add = "median")+stat_compare_means(label.y = 2)+theme(aspect.ratio = 1)+scale_y_continuous(breaks=seq(0, 1000, 6000))+ylim(0,6000)

p6<-ggviolin(Pamper.pdata[which(Pamper.pdata$Proteome==1),], x = "outcome", y = "PH_crystalloid",
             color = "outcome",add = "median")+stat_compare_means(label.y = 2)+theme(aspect.ratio = 1)+scale_y_continuous(breaks=seq(0, 1000, 6000))+ylim(0,6000)

ggarrange(plotlist = list(p1,p3,p5,p2,p4,p6), 
          ncol = 3, nrow = 2,common.legend = T,align = "v")

