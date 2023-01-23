rm(list=ls())
library(rio)
library(tidyverse)
library('ez')
library(agricolae)
library(car)
library(lme4)
library(nlme)
require(multcomp)
library(emmeans)
#test<-Table.1.split[[1]]
#define HSD function with lme--------
analysis_HSD<-function(test){
  variable<-'value'
  ## factor and value
  ozone <- as.factor(test$ozone)
  ring <- as.factor(test$chamber)
  DOF<-as.factor(test$DOF)
  replicate<- test$replicate
  clone <- as.factor(test$clone)
  value<-test[,which(colnames(test)==variable)][[1]]
  DF.test<-data.frame(ozone,clone,DOF,ring,replicate,value)
  DF.test$Tukey_x<-with(DF.test,interaction(ozone,clone,DOF))
  #确定排序大小
  DF.test %>% dplyr::select(value,Tukey_x) %>% group_by(Tukey_x) %>% summarise_all(c('mean')) %>% arrange(desc(value))->order
  DF.test$Tukey_x <- factor(DF.test$Tukey_x, levels=order$Tukey_x, ordered=F)
  mod1<-lme(log(value)~Tukey_x,random=~1|ring/replicate,data=DF.test)
  HSD.result<-emmeans(mod1, list(pairwise ~ Tukey_x), adjust = "tukey")
  HSD.letter<-as.data.frame(HSD.result$`pairwise differences of Tukey_x`)
  names(HSD.letter)<-c("name","estimate","SE","df","t.ratio","p.value")
  HSD.letter<-rcompanion::cldList(p.value ~ name,
                                  data = HSD.letter,
                                  threshold  = 0.05)
  ozone<-data.frame(matrix(unlist(strsplit(HSD.letter$Group,'\\.')),ncol=3,byrow=T))[,1]
  clone<-data.frame(matrix(unlist(strsplit(HSD.letter$Group,'\\.')),ncol=3,byrow=T))[,2]
  DOF<-data.frame(matrix(unlist(strsplit(HSD.letter$Group,'\\.')),ncol=3,byrow=T))[,3]
  result<-data.frame(DOF,clone,ozone,name=test$name[1],HSD=HSD.letter$Letter)
  result$clone[which(result$clone=='17')]="107"
  result$ozone[which(result$ozone=='NF4')]="NF40"
  return(result)
}

#test<-Table.4.split[[12]]
#define the function for split-plot analysis
make_two_statistics<- function(test) { 
  variable<-'value'
  ## factor and value
  ozone <- as.factor(test$ozone)
  ring <- as.factor(test$chamber)
  replicate<- test$replicate
  treatment <- as.factor(test$clone)
  value<-test[,which(colnames(test)==variable)][[1]]
  DF.test<-data.frame(ozone,ring,replicate,treatment,value)
  ### lme
  mod1<-lme(log(value)~ozone*treatment,random=~1|ring/replicate,data=DF.test)
  summary(mod1)
  #Result.Tukey<-emmeans(mod1, list(pairwise ~ ozone*treatment), adjust = "tukey")
  #Testing normality of residuals
  shapiro.test(resid(mod1))
  # anova
  anov<-anova(mod1)
  O3.p.value<-anov$`p-value`[2]
  Tre.p.value<-anov$`p-value`[3]
  Inter.p.value<-anov$`p-value`[4]
  #clone
  test_107<-subset(DF.test,DF.test$treatment=='107')
  mod1<-lme(log(value)~ozone,random=~1|ring/replicate,data=test_107)
  anov<-anova(mod1)
  ozone_107<-anov$`p-value`[2]
  test_546<-subset(DF.test,DF.test$treatment=='546')
  mod1<-lme(log(value)~ozone,random=~1|ring/replicate,data=test_546)
  anov<-anova(mod1)
  ozone_546<-anov$`p-value`[2]
  #ozone
  test_NF<-subset(DF.test,DF.test$ozone=='NF')
  mod1<-lme(log(value)~treatment,random=~1|ring/replicate,data=test_NF)
  anov<-anova(mod1)
  clone_NF<-anov$`p-value`[2]
  test_EO3<-subset(DF.test,DF.test$ozone=='NF40')
  mod1<-lme(log(value)~treatment,random=~1|ring/replicate,data=test_EO3)
  anov<-anova(mod1)
  clone_EO3<-anov$`p-value`[2]
  #输出结果
  statis.result<-data.frame(variable=test$name[1],O3.p.value,Tre.p.value,Inter.p.value,ozone_107,ozone_546,clone_NF,clone_EO3)
  return(statis.result)
}
#define the function for split-plot analysis
#test<-Table.1.split[[1]]
make_three_statistics<- function(test) { 
  variable<-'value'
  ## factor and value
  ozone <- as.factor(test$ozone)
  ring <- as.factor(test$chamber)
  time<-test$DOF
  treatment <- as.factor(test$clone)
  pot<-as.factor(test$pot)
  value<-test[,which(colnames(test)==variable)][[1]]
  replicate<-test$replicate
  DF.test<-data.frame(ozone,ring,treatment,time,replicate,value,pot)
  ### lme
  mod1<-lme(log(value)~ozone*treatment*time,random=~1|ring/replicate,data=DF.test)
  summary(mod1)
  #Testing normality of residuals
  shapiro.test(resid(mod1))
  # anova
  anov<-anova(mod1)
  O3.p.value<-anov$`p-value`[2]
  Clo.p.value<-anov$`p-value`[3]
  Tim.p.value<-anov$`p-value`[4]
  O3.Clo.p.value<-anov$`p-value`[5]
  O3.Tim.p.value<-anov$`p-value`[6]
  Clo.Tim.p.value<-anov$`p-value`[7]
  O3.Clo.Tim.p.value<-anov$`p-value`[8]
  #O3×clone
  test_1<-subset(DF.test,DF.test$time=='63D')
  mod1<-lme(log(value)~ozone*treatment,random=~1|ring/replicate,data=test_1)
  anov<-anova(mod1)
  O3.clone_Aug<-anov$`p-value`[4]
  test_2<-subset(DF.test,DF.test$time=='94D')
  mod1<-lme(log(value)~ozone*treatment,random=~1|ring/replicate,data=test_2)
  anov<-anova(mod1)
  O3.clone_Sep<-anov$`p-value`[4]
  #O3×DOF
  test_1<-subset(DF.test,DF.test$treatment=='107')
  mod1<-lme(log(value)~ozone*time,random=~1|ring/replicate,data=test_1)
  anov<-anova(mod1)
  O3.DOF_107<-anov$`p-value`[4]
  test_2<-subset(DF.test,DF.test$treatment=='546')
  mod1<-lme(log(value)~ozone*time,random=~1|ring/replicate,data=test_2)
  anov<-anova(mod1)
  O3.DOF_546<-anov$`p-value`[4]
  #clone×DOF
  test_1<-subset(DF.test,DF.test$ozone=='NF')
  mod1<-lme(log(value)~treatment*time,random=~1|ring/replicate,data=test_1)
  anov<-anova(mod1)
  clone.DOF_NF<-anov$`p-value`[4]
  test_2<-subset(DF.test,DF.test$ozone=='NF40')
  mod1<-lme(log(value)~treatment*time,random=~1|ring/replicate,data=test_2)
  anov<-anova(mod1)
  clone.DOF_NF40<-anov$`p-value`[4]
  #O3
  DF.test %>% filter(treatment=='546'&time=='63D')->test_1
  mod1<-lme(log(value)~ozone,random=~1|ring/replicate,data=test_1)
  anov<-anova(mod1)
  O3_546_Aug<-anov$`p-value`[2]
  DF.test %>% filter(treatment=='546'&time=='94D')->test_1
  mod1<-lme(log(value)~ozone,random=~1|ring/replicate,data=test_1)
  anov<-anova(mod1)
  O3_546_Sep<-anov$`p-value`[2]
  DF.test %>% filter(treatment=='107'&time=='63D')->test_1
  mod1<-lme(log(value)~ozone,random=~1|ring/replicate,data=test_1)
  anov<-anova(mod1)
  O3_107_Aug<-anov$`p-value`[2]
  DF.test %>% filter(treatment=='107'&time=='94D')->test_1
  mod1<-lme(log(value)~ozone,random=~1|ring/replicate,data=test_1)
  anov<-anova(mod1)
  O3_107_Sep<-anov$`p-value`[2]
  #clone
  DF.test %>% filter(ozone=='NF'&time=='63D')->test_1
  mod1<-lme(log(value)~treatment,random=~1|ring/replicate,data=test_1)
  anov<-anova(mod1)
  clone_NF_Aug<-anov$`p-value`[2]
  DF.test %>% filter(ozone=='NF'&time=='94D')->test_1
  mod1<-lme(log(value)~treatment,random=~1|ring/replicate,data=test_1)
  anov<-anova(mod1)
  clone_NF_Sep<-anov$`p-value`[2]
  DF.test %>% filter(ozone=='NF40'&time=='63D')->test_1
  mod1<-lme(log(value)~treatment,random=~1|ring/replicate,data=test_1)
  anov<-anova(mod1)
  clone_NF40_Aug<-anov$`p-value`[2]
  DF.test %>% filter(ozone=='NF40'&time=='94D')->test_1
  mod1<-lme(log(value)~treatment,random=~1|ring/replicate,data=test_1)
  anov<-anova(mod1)
  clone_NF40_Sep<-anov$`p-value`[2]
  #TIME
  DF.test %>% filter(ozone=='NF'&treatment=='107')->test_1
  mod1<-lme(log(value)~time,random=~1|ring/replicate,data=test_1)
  anov<-anova(mod1)
  Time_NF_107<-anov$`p-value`[2]
  DF.test %>% filter(ozone=='NF40'&treatment=='107')->test_1
  mod1<-lme(log(value)~time,random=~1|ring/replicate,data=test_1)
  anov<-anova(mod1)
  Time_NF40_107<-anov$`p-value`[2]
  DF.test %>% filter(ozone=='NF'&treatment=='546')->test_1
  mod1<-lme(log(value)~time,random=~1|ring/replicate,data=test_1)
  anov<-anova(mod1)
  Time_NF_546<-anov$`p-value`[2]
  DF.test %>% filter(ozone=='NF40'&treatment=='546')->test_1
  mod1<-lme(log(value)~time,random=~1|ring/replicate,data=test_1)
  anov<-anova(mod1)
  Time_NF40_546<-anov$`p-value`[2]
  statis.result<-data.frame(variable=test$name[1],O3.p.value,Clo.p.value,Tim.p.value,O3.Clo.p.value,O3.Tim.p.value,Clo.Tim.p.value,O3.Clo.Tim.p.value,O3.clone_Aug,O3.clone_Sep,O3.DOF_107,O3.DOF_546,clone.DOF_NF,clone.DOF_NF40,O3_107_Aug,O3_546_Aug,O3_107_Sep,O3_546_Sep,clone_NF_Aug,clone_NF_Sep,clone_NF40_Aug,clone_NF40_Sep,Time_NF_107,Time_NF40_107,Time_NF_546,Time_NF40_546)
  return(statis.result)
}

main.folder<-"D:/OneDrive - 南京信息工程大学/2019年/00实验内容/0实验数据再处理/臭氧对叶片解剖结构的影响/Manuscript/Data submited to Journal/measured data"

#Table1-------
Table.1<-import(file.path(main.folder,"Table_1.xlsx"))
#HSD
Table.1->out.put.Table.1
Table.1 %>% pivot_longer(cols = 7:14)->Table.1.long
Table.1.split<-split(Table.1.long,Table.1.long$name)
Table.1.HSD<-lapply(Table.1.split,analysis_HSD)
Table.1.HSD <- plyr::ldply(Table.1.HSD,data.frame)
Table.1.HSD %>% dplyr::select(-.id)->Table.1.HSD
#three ways anova
Table.1.anova<-lapply(Table.1.split,make_three_statistics)
Table.1.anova <- plyr::ldply(Table.1.anova,data.frame)
Table.1.anova %>% dplyr::select(-.id)->Table.1.anova
#avg
Table.1.long %>% ungroup()%>% dplyr::select(-chamber,-pot,-replicate)%>% group_by(name,DOF,clone,ozone) %>% summarise_all(c('mean','sd')) %>% mutate(clone=as.character(clone)) ->Table.1.avg
Table.1.HSD %>% left_join(Table.1.avg)->Table.1
Table.1<-list(out.put.Table.1,Table.1,Table.1.anova)
export(Table.1,file.path(main.folder,"Output data","Table_1.xlsx"),overwrite=T)

#Table2-------
Table.2<-import(file.path(main.folder,"Table_2.xlsx"))
#HSD
Table.2->out.put.Table.2
Table.2 %>% pivot_longer(cols = 7:14)->Table.2.long
Table.2.split<-split(Table.2.long,Table.2.long$name)
Table.2.HSD<-lapply(Table.2.split,analysis_HSD)
Table.2.HSD <- plyr::ldply(Table.2.HSD,data.frame)
Table.2.HSD %>% dplyr::select(-.id)->Table.2.HSD
#two ways anova
Table.2.anova<-lapply(Table.2.split,make_two_statistics)
Table.2.anova <- plyr::ldply(Table.2.anova,data.frame)
Table.2.anova %>% dplyr::select(-.id)->Table.2.anova
#avg
Table.2.long %>% ungroup()%>% dplyr::select(-chamber,-pot,-replicate)%>% group_by(name,DOF   ,clone,ozone) %>% summarise_all(c('mean','sd')) %>% mutate(clone=as.character(clone)) ->Table.2.avg
Table.2.HSD %>% left_join(Table.2.avg)->Table.2
Table.2<-list(out.put.Table.2,Table.2,Table.2.anova)
export(Table.2,file.path(main.folder,"Output data","Table_2.xlsx"),overwrite=T)

#Table3---------
Table.3<-import(file.path(main.folder,"Table_3.xlsx"))
#HSD
Table.3->out.put.Table.3
Table.3 %>% pivot_longer(cols = 7:21)->Table.3.long
Table.3.split<-split(Table.3.long,Table.3.long$name)
Table.3.HSD<-lapply(Table.3.split,analysis_HSD)
Table.3.HSD <- plyr::ldply(Table.3.HSD,data.frame)
Table.3.HSD %>% dplyr::select(-.id)->Table.3.HSD
#anova
Table.3.anova<-lapply(Table.3.split,make_two_statistics)
Table.3.anova <- plyr::ldply(Table.3.anova,data.frame)
Table.3.anova %>% dplyr::select(-.id)->Table.3.anova
#avg
Table.3.long %>% dplyr::select(-chamber,-pot,-replicate) %>%  group_by(name,DOF,clone,ozone) %>% summarise_all(c('mean','sd')) %>% mutate(clone=as.character(clone)) ->Table.3.avg
Table.3.HSD %>% left_join(Table.3.avg)->Table.3
Table.3<-list(out.put.Table.3,Table.3,Table.3.anova)
export(Table.3,file.path(main.folder,"Output data","Table_3.xlsx"),overwrite=T)

#Table4---------
Table.4<-import(file.path(main.folder,"Table_4.xlsx"))
#HSD
Table.4->out.put.Table.4
Table.4 %>% pivot_longer(cols = 7:22)->Table.4.long
Table.4.split<-split(Table.4.long,Table.4.long$name)
Table.4.HSD<-lapply(Table.4.split,analysis_HSD)
Table.4.HSD <- plyr::ldply(Table.4.HSD,data.frame)
Table.4.HSD %>% dplyr::select(-.id)->Table.4.HSD
#anova
Table.4.anova<-lapply(Table.4.split,make_two_statistics)
Table.4.anova <- plyr::ldply(Table.4.anova,data.frame)
Table.4.anova %>% dplyr::select(-.id)->Table.4.anova
#avg
Table.4.long %>% dplyr::select(-chamber,-pot,-replicate) %>%  group_by(name,DOF,clone,ozone) %>% summarise_all(c('mean','sd')) %>% mutate(clone=as.character(clone)) ->Table.4.avg
Table.4.HSD %>% left_join(Table.4.avg)->Table.4
Table.4<-list(out.put.Table.4,Table.4,Table.4.anova)
export(Table.4,file.path(main.folder,"Output data","Table_4.xlsx"),overwrite=T)


#Figure 6 and 7 gm.pa----------
figure.stat<-import(file.path(main.folder,"Figure6 and 7.xlsx"))
#HSD
figure.stat->out.put.figure.stat
figure.stat %>% pivot_longer(cols = 7:15)->figure.long
figure.long.split<-split(figure.long,figure.long$name)
figure.HSD<-lapply(figure.long.split,analysis_HSD)
figure.HSD <- plyr::ldply(figure.HSD,data.frame)
figure.HSD %>% dplyr::select(-.id)->figure.HSD
#anova
figure.anova<-lapply(figure.long.split,make_two_statistics)
figure.anova <- plyr::ldply(figure.anova,data.frame)
figure.anova %>% dplyr::select(-.id)->figure.anova

#avg
figure.long %>% ungroup()%>%  dplyr::select(-chamber,-pot,-replicate) %>%  group_by(name,DOF,clone,ozone) %>% summarise_all(c('mean','sd')) %>% mutate(clone=as.character(clone)) ->figure.avg
figure.HSD %>% left_join(figure.avg)->figure
figure<-list(out.put.figure.stat,figure,figure.anova)
export(figure,file.path(main.folder,"Output data","figure.xlsx"),overwrite=T)



