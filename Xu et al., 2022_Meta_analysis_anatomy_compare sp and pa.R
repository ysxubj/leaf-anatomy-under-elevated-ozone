rm(list=ls())
library(rio)
library(tidyverse)
library(metafor)
library(moments)
library(normtest)
library(ggh4x)
library(rstatix)
library(ggprism)

windowsFonts(Times=windowsFont("Times New Roman"))
#set main folder path---------
main.folder<-"D:/OneDrive - 南京信息工程大学/2019年/00实验内容/0实验数据再处理/臭氧对叶片解剖结构的影响/Manuscript/Data submited to Journal/meta data"
Output.folder<-"D:/OneDrive - 南京信息工程大学/2019年/00实验内容/0实验数据再处理/臭氧对叶片解剖结构的影响/Manuscript/Data submited to Journal/meta data"


#Pvalue fun-----
Pvalue.fun<-function(value){
  if(value<0.001) {value<-c('***')
  }else if(value<0.01){value<-c('**')
  }else if(value<0.05){value<-c('*')
  }else{value<-c("ns")}
  return(value)
}

#mian prog-------
DF<-import(file.path(main.folder,'Effect of ozone on leaf anatomy.xlsx'))
DF %>% filter(Accept==1,No.paired.data!=1,Tissue %in% c("spongy","palisade"))->DF.metaanalysis

#Intercellular space of mesophyll----------
variable<-"Intercellular space of mesophyll"
DF.metaanalysis %>% filter(paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence.withouttissue,RR.adj=yi*w.adj)->DF.RR
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue)-1,data=DF.RR,weighted = T, method="REML")
print(res)

#outliers test
DF.RR$standresidual<-rstandard(res)$z
DF.RR<-cbind(DF.RR,influence(res)[[1]])

#filter standardized residual value >3
DF.RR %>% filter(abs(standresidual)>3)->DF.RR.rmove
No.Rmove<-length(DF.RR.rmove$Accept)
print(DF.RR.rmove)
DF.RR %>% filter(abs(standresidual)<3)->DF.RR.removeouter

palisade<-DF.RR.removeouter %>% filter(Tissue=='palisade')
spongy<-DF.RR.removeouter %>% filter(Tissue=='spongy')
No.use<-c(length(palisade$Accept),length(spongy$Accept))

No.species<-c(length(unique(paste0(palisade$species,palisade$cultivar.genotype))),length(unique(paste0(spongy$species,spongy$cultivar.genotype))))

#normality test
Sle.P.value<-skewness.norm.test(DF.RR.removeouter$standresidual)
print(Sle.P.value)
skewness.value<-Sle.P.value$statistic[[1]]
skewness.Pvalue<-Sle.P.value$p.value[[1]]
Kutr.P.value<-kurtosis.norm.test(DF.RR.removeouter$standresidual)
print(Kutr.P.value)
kurtosis.value<-Kutr.P.value$statistic[[1]]
kurtosis.Pvalue<-Kutr.P.value$p.value[[1]]
hist(DF.RR$standresidual, prob = TRUE)
lines(density(DF.RR$standresidual), col = 2, lwd = 3)


#RR average calculation
res<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue)-1,data=DF.RR.removeouter,weighted = T)
print(res)

res.compare<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue),data=DF.RR.removeouter,weighted = T, method="REML")
print(res.compare)
Pvalue.compare<-res.compare$QMp


ci.lb<-res$ci.lb
ci.ub<-res$ci.ub
mean.RR<-res$beta[,1]
P.value<-res$pval

#output
output.data<-data.frame(variable,Tissue=c("Palisade","Spongy"),mean.RR,ci.ub,ci.lb,P.value,Pvalue.compare,skewness.value,skewness.Pvalue,kurtosis.value,kurtosis.Pvalue,No.use,No.species)

Intercel.output<-output.data

#cell wall thickness-----
variable<-"cell wall thickness"
DF.metaanalysis %>% filter(paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence.withouttissue,RR.adj=yi*w.adj)->DF.RR
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue)-1,data=DF.RR,weighted = T, method="REML")
print(res)

#outliers test
DF.RR$standresidual<-rstandard(res)$z
DF.RR<-cbind(DF.RR,influence(res)[[1]])

#filter standardized residual value >3
DF.RR %>% filter(abs(standresidual)>3)->DF.RR.rmove
No.Rmove<-length(DF.RR.rmove$Accept)
print(DF.RR.rmove)
DF.RR %>% filter(abs(standresidual)<3)->DF.RR.removeouter

palisade<-DF.RR.removeouter %>% filter(Tissue=='palisade')
spongy<-DF.RR.removeouter %>% filter(Tissue=='spongy')
No.use<-c(length(palisade$Accept),length(spongy$Accept))

No.species<-c(length(unique(paste0(palisade$species,palisade$cultivar.genotype))),length(unique(paste0(spongy$species,spongy$cultivar.genotype))))

#normality test
Sle.P.value<-skewness.norm.test(DF.RR.removeouter$standresidual)
print(Sle.P.value)
skewness.value<-Sle.P.value$statistic[[1]]
skewness.Pvalue<-Sle.P.value$p.value[[1]]
Kutr.P.value<-kurtosis.norm.test(DF.RR.removeouter$standresidual)
print(Kutr.P.value)
kurtosis.value<-Kutr.P.value$statistic[[1]]
kurtosis.Pvalue<-Kutr.P.value$p.value[[1]]
hist(DF.RR$standresidual, prob = TRUE)
lines(density(DF.RR$standresidual), col = 2, lwd = 3)


#RR average calculation
res<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue)-1,data=DF.RR.removeouter,weighted = T)
print(res)

res.compare<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue),data=DF.RR.removeouter,weighted = T, method="REML")
print(res.compare)
Pvalue.compare<-res.compare$QMp


ci.lb<-res$ci.lb
ci.ub<-res$ci.ub
mean.RR<-res$beta[,1]
P.value<-res$pval

#output
output.data<-data.frame(variable,Tissue=c("Palisade","Spongy"),mean.RR,ci.ub,ci.lb,P.value,Pvalue.compare,skewness.value,skewness.Pvalue,kurtosis.value,kurtosis.Pvalue,No.use,No.species)

cellwall.output<-output.data

#chloroplast size------

variable<-"chloroplast size"
DF.metaanalysis %>% filter(paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence.withouttissue,RR.adj=yi*w.adj)->DF.RR
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue)-1,data=DF.RR,weighted = T, method="REML")
print(res)

#outliers test
DF.RR$standresidual<-rstandard(res)$z
DF.RR<-cbind(DF.RR,influence(res)[[1]])

#filter standardized residual value >3
DF.RR %>% filter(abs(standresidual)>3)->DF.RR.rmove
No.Rmove<-length(DF.RR.rmove$Accept)
print(DF.RR.rmove)
DF.RR %>% filter(abs(standresidual)<3)->DF.RR.removeouter

palisade<-DF.RR.removeouter %>% filter(Tissue=='palisade')
spongy<-DF.RR.removeouter %>% filter(Tissue=='spongy')
No.use<-c(length(palisade$Accept),length(spongy$Accept))

No.species<-c(length(unique(paste0(palisade$species,palisade$cultivar.genotype))),length(unique(paste0(spongy$species,spongy$cultivar.genotype))))

#normality test
Sle.P.value<-skewness.norm.test(DF.RR.removeouter$standresidual)
print(Sle.P.value)
skewness.value<-Sle.P.value$statistic[[1]]
skewness.Pvalue<-Sle.P.value$p.value[[1]]
Kutr.P.value<-kurtosis.norm.test(DF.RR.removeouter$standresidual)
print(Kutr.P.value)
kurtosis.value<-Kutr.P.value$statistic[[1]]
kurtosis.Pvalue<-Kutr.P.value$p.value[[1]]
hist(DF.RR$standresidual, prob = TRUE)
lines(density(DF.RR$standresidual), col = 2, lwd = 3)


#RR average calculation
res<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue)-1,data=DF.RR.removeouter,weighted = T)
print(res)

res.compare<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue),data=DF.RR.removeouter,weighted = T, method="REML")
print(res.compare)
Pvalue.compare<-res.compare$QMp


ci.lb<-res$ci.lb
ci.ub<-res$ci.ub
mean.RR<-res$beta[,1]
P.value<-res$pval

#output
output.data<-data.frame(variable,Tissue=c("Palisade","Spongy"),mean.RR,ci.ub,ci.lb,P.value,Pvalue.compare,skewness.value,skewness.Pvalue,kurtosis.value,kurtosis.Pvalue,No.use,No.species)

chloroplastsize.output<-output.data

#chloroplast width-------

variable<-"chloroplast width"
DF.metaanalysis %>% filter(paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence.withouttissue,RR.adj=yi*w.adj)->DF.RR
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue)-1,data=DF.RR,weighted = T, method="REML")
print(res)

#outliers test
DF.RR$standresidual<-rstandard(res)$z
DF.RR<-cbind(DF.RR,influence(res)[[1]])

#filter standardized residual value >3
DF.RR %>% filter(abs(standresidual)>3)->DF.RR.rmove
No.Rmove<-length(DF.RR.rmove$Accept)
print(DF.RR.rmove)
DF.RR %>% filter(abs(standresidual)<3)->DF.RR.removeouter

palisade<-DF.RR.removeouter %>% filter(Tissue=='palisade')
spongy<-DF.RR.removeouter %>% filter(Tissue=='spongy')
No.use<-c(length(palisade$Accept),length(spongy$Accept))

No.species<-c(length(unique(paste0(palisade$species,palisade$cultivar.genotype))),length(unique(paste0(spongy$species,spongy$cultivar.genotype))))

#normality test
Sle.P.value<-skewness.norm.test(DF.RR.removeouter$standresidual)
print(Sle.P.value)
skewness.value<-Sle.P.value$statistic[[1]]
skewness.Pvalue<-Sle.P.value$p.value[[1]]
Kutr.P.value<-kurtosis.norm.test(DF.RR.removeouter$standresidual)
print(Kutr.P.value)
kurtosis.value<-Kutr.P.value$statistic[[1]]
kurtosis.Pvalue<-Kutr.P.value$p.value[[1]]
hist(DF.RR$standresidual, prob = TRUE)
lines(density(DF.RR$standresidual), col = 2, lwd = 3)


#RR average calculation
res<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue)-1,data=DF.RR.removeouter,weighted = T)
print(res)

res.compare<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue),data=DF.RR.removeouter,weighted = T, method="REML")
print(res.compare)
Pvalue.compare<-res.compare$QMp


ci.lb<-res$ci.lb
ci.ub<-res$ci.ub
mean.RR<-res$beta[,1]
P.value<-res$pval

#output
output.data<-data.frame(variable,Tissue=c("Palisade","Spongy"),mean.RR,ci.ub,ci.lb,P.value,Pvalue.compare,skewness.value,skewness.Pvalue,kurtosis.value,kurtosis.Pvalue,No.use,No.species)

chloroplastwidth.output<-output.data

#chloroplast length----

variable<-"chloroplast length"
DF.metaanalysis %>% filter(paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence.withouttissue,RR.adj=yi*w.adj)->DF.RR
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue)-1,data=DF.RR,weighted = T, method="REML")
print(res)

#outliers test
DF.RR$standresidual<-rstandard(res)$z
DF.RR<-cbind(DF.RR,influence(res)[[1]])

#filter standardized residual value >3
DF.RR %>% filter(abs(standresidual)>3)->DF.RR.rmove
No.Rmove<-length(DF.RR.rmove$Accept)
print(DF.RR.rmove)
DF.RR %>% filter(abs(standresidual)<3)->DF.RR.removeouter

palisade<-DF.RR.removeouter %>% filter(Tissue=='palisade')
spongy<-DF.RR.removeouter %>% filter(Tissue=='spongy')
No.use<-c(length(palisade$Accept),length(spongy$Accept))

No.species<-c(length(unique(paste0(palisade$species,palisade$cultivar.genotype))),length(unique(paste0(spongy$species,spongy$cultivar.genotype))))

#normality test
Sle.P.value<-skewness.norm.test(DF.RR.removeouter$standresidual)
print(Sle.P.value)
skewness.value<-Sle.P.value$statistic[[1]]
skewness.Pvalue<-Sle.P.value$p.value[[1]]
Kutr.P.value<-kurtosis.norm.test(DF.RR.removeouter$standresidual)
print(Kutr.P.value)
kurtosis.value<-Kutr.P.value$statistic[[1]]
kurtosis.Pvalue<-Kutr.P.value$p.value[[1]]
hist(DF.RR$standresidual, prob = TRUE)
lines(density(DF.RR$standresidual), col = 2, lwd = 3)


#RR average calculation
res<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue)-1,data=DF.RR.removeouter,weighted = T)
print(res)

res.compare<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue),data=DF.RR.removeouter,weighted = T, method="REML")
print(res.compare)
Pvalue.compare<-res.compare$QMp


ci.lb<-res$ci.lb
ci.ub<-res$ci.ub
mean.RR<-res$beta[,1]
P.value<-res$pval

#output
output.data<-data.frame(variable,Tissue=c("Palisade","Spongy"),mean.RR,ci.ub,ci.lb,P.value,Pvalue.compare,skewness.value,skewness.Pvalue,kurtosis.value,kurtosis.Pvalue,No.use,No.species)

chloroplastlength.output<-output.data

#Starch/chloroplast-----
variable<-"Starch/chloroplast"
DF.metaanalysis %>% filter(paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence.withouttissue,RR.adj=yi*w.adj)->DF.RR
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue)-1,data=DF.RR,weighted = T, method="REML")
print(res)

#outliers test
DF.RR$standresidual<-rstandard(res)$z
DF.RR<-cbind(DF.RR,influence(res)[[1]])

#filter standardized residual value >3
DF.RR %>% filter(abs(standresidual)>3)->DF.RR.rmove
No.Rmove<-length(DF.RR.rmove$Accept)
print(DF.RR.rmove)
DF.RR %>% filter(abs(standresidual)<3)->DF.RR.removeouter

palisade<-DF.RR.removeouter %>% filter(Tissue=='palisade')
spongy<-DF.RR.removeouter %>% filter(Tissue=='spongy')
No.use<-c(length(palisade$Accept),length(spongy$Accept))

No.species<-c(length(unique(paste0(palisade$species,palisade$cultivar.genotype))),length(unique(paste0(spongy$species,spongy$cultivar.genotype))))

#normality test
Sle.P.value<-skewness.norm.test(DF.RR.removeouter$standresidual)
print(Sle.P.value)
skewness.value<-Sle.P.value$statistic[[1]]
skewness.Pvalue<-Sle.P.value$p.value[[1]]
Kutr.P.value<-kurtosis.norm.test(DF.RR.removeouter$standresidual)
print(Kutr.P.value)
kurtosis.value<-Kutr.P.value$statistic[[1]]
kurtosis.Pvalue<-Kutr.P.value$p.value[[1]]
hist(DF.RR$standresidual, prob = TRUE)
lines(density(DF.RR$standresidual), col = 2, lwd = 3)


#RR average calculation
res<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue)-1,data=DF.RR.removeouter,weighted = T)
print(res)

res.compare<- rma(yi,vi,weights = w.adj,mods= ~ factor(Tissue),data=DF.RR.removeouter,weighted = T, method="REML")
print(res.compare)
Pvalue.compare<-res.compare$QMp


ci.lb<-res$ci.lb
ci.ub<-res$ci.ub
mean.RR<-res$beta[,1]
P.value<-res$pval

#output
output.data<-data.frame(variable,Tissue=c("Palisade","Spongy"),mean.RR,ci.ub,ci.lb,P.value,Pvalue.compare,skewness.value,skewness.Pvalue,kurtosis.value,kurtosis.Pvalue,No.use,No.species)

starch.chloroplast.output<-output.data


#Plot leaf structural traits----------
anatomy.plot<-rbind(Intercel.output,cellwall.output,chloroplastsize.output,chloroplastwidth.output,chloroplastlength.output,starch.chloroplast.output)


#plot response
anatomy.plot$mean.RR<-(exp(anatomy.plot$mean.RR)-1)*100
anatomy.plot$ci.ub<-(exp(anatomy.plot$ci.ub)-1)*100
anatomy.plot$ci.lb<-(exp(anatomy.plot$ci.lb)-1)*100
anatomy.plot %>% arrange(mean.RR)->anatomy.plot
anatomy.plot$variable<-factor(anatomy.plot$variable,levels=unique(anatomy.plot$variable))
anatomy.plot$Pvalue.index<-as.vector(lapply(anatomy.plot$P.value,Pvalue.fun))
anatomy.plot$P.value.colour<-ifelse (anatomy.plot$P.value>0.05,"a","b")
anatomy.plot <- add_significance(anatomy.plot, 'Pvalue.compare')

anatomy.plot %>% filter(Tissue=='Spongy')->anatomy.plot.xy
anatomy.plot.xy$xmin<-as.numeric(anatomy.plot.xy$variable)-0.2
anatomy.plot.xy$xmax<-as.numeric(anatomy.plot.xy$variable)+0.2
anatomy.plot.xy$y.position<-anatomy.plot.xy$ci.ub+28


yaixs<-c(expression(italic(Area)[chl]),
         expression(italic(L)[chl]),
         expression(italic(T)[chl]),
         expression(italic(Area)[starch/chloroplast]),
         expression(italic(f)[ias]),
         expression(italic(T)[cw]))


dodge <- position_dodge(width=0.65)
ggplot(anatomy.plot,aes(x=(variable),y=mean.RR))+
  geom_hline(yintercept = 0,linetype="dotted")+
  geom_point(aes(colour=Tissue),size=2,position = dodge)+
  geom_point(aes(colour=Tissue),show.legend=FALSE,size=2,position = dodge)+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub,colour=Tissue),width=0.2,size=0.6,position = dodge)+
  geom_text(aes(x =(variable),y = ci.ub+10,color=Tissue),label = mapply(function(No,species) as.expression(bquote(.(No)~"("*.(species)*")")), anatomy.plot$No.use, anatomy.plot$No.species),size = 3,fontface = "bold",position = dodge,family="Times",show.legend=FALSE)+
  scale_x_discrete(labels = yaixs)+
  scale_colour_manual(values=c("#DC4401", "#0564C9"))+
  expand_limits(y=c(-80,80))+
  scale_y_continuous(breaks=seq(-80,80,20))+
  coord_flip()+
  add_pvalue(anatomy.plot.xy,
             xmin = "xmin",
             xmax = "xmax",
             label = "{Pvalue.compare.signif}",
             tip.length = c(0.02,0.02),
             coord.flip = TRUE,
             bracket.size=0.5
  )+
  theme_classic()+
  theme(
    text=element_text(family="Times"),
    plot.title =element_text(hjust = .5,face = "bold",size = 1.5),
    plot.subtitle =element_text(hjust = .6),
    #图例
    legend.background = element_blank(),
    legend.title=element_blank(),
    #legend.justification=c(1,1), 
    #legend.position=c(1.03,1.05),
    legend.key.size = unit(.6, "cm"),
    legend.key.width = unit(0.6,"cm"),
    legend.text = element_text(hjust = 0,colour="black", face = "bold"),
    #图底板
    panel.background = element_rect(fill = 'white'),
    panel.border=element_rect(fill = NA,colour = "black",size = 0.6),
    panel.spacing.y = unit(.5,"cm"),
    panel.grid = element_blank(),
    #坐标轴
    axis.ticks.length = unit(0.4,"lines"), 
    axis.ticks = element_line(color='black'),
    axis.line = element_line(colour = "black"),
    #坐标轴文字
    axis.title.y=element_text(colour='black', size=10,face = "bold",vjust = 0),
    axis.title.x=element_text(colour='black', size=10,face = "bold",vjust = -1),
    #轴刻度
    axis.text.y=element_text(colour='black',size=12),
    axis.text.x=element_text(colour = "black",size = 12,angle = 0,hjust = 0.5,vjust = 0),
    #分面标题
    strip.text = element_text(colour = "black",size = 9,face = "bold"),
    strip.background = element_rect(fill=NA,size=0.6,color="black"))+
  labs(y=as.expression("Change under elevated"~O[3]~"(%)"))+
  labs(x=c())
ggsave(file.path(Output.folder,'compared sp and pa.png'),width = 14, height = 13,units = "cm",dpi=300)










