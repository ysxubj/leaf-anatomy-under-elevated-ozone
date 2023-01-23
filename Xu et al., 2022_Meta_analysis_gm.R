rm(list=ls())
library(rio)
library(tidyverse)
library(metafor)
library(moments)
library(normtest)
library(ggh4x)
windowsFonts(Times=windowsFont("Times New Roman"))
#set main folder path---------
main.folder<-"D:/OneDrive - 南京信息工程大学/2019年/00实验内容/0实验数据再处理/臭氧对叶片解剖结构的影响/Manuscript/Data submited to Journal/meta data"
Output.folder<-"D:/OneDrive - 南京信息工程大学/2019年/00实验内容/0实验数据再处理/臭氧对叶片解剖结构的影响/Manuscript/Data submited to Journal/meta data"

#define function---------
meta.fun<-function(variable){
  #gm---------
  DF %>% filter(Accept==1&paramter==variable)->DF.analysis
  # calculate mean differences and corresponding sampling variances (RR and v)
  DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=N, m2i=AA, sd2i=A.SD, n2i=N, data=DF.analysis)
  #independence and weights
  DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
  hist(DF.RR$yi, prob = TRUE)
  lines(density(DF.RR$yi), col = 2, lwd = 3)
  
  
  #fit mixed-effects models
  res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
  #print(res)
  #outliers test
  DF.RR$standresidual<-rstandard(res)$z
  DF.RR<-cbind(DF.RR,influence(res)[[1]])
  
  
  #leave1out test 
  remove.test<-leave1out(res)
  
  #filter standardized residual value >3
  DF.RR %>% filter(abs(standresidual)>3)->DF.RR.rmove
  No.Rmove<-length(DF.RR.rmove$Accept)
  print(DF.RR.rmove)
  DF.RR %>% filter(abs(standresidual)<3)->DF.RR.removeouter
  No.use<-length(DF.RR.removeouter$Accept)
  No.species<-length(unique(paste0(DF.RR.removeouter$species,DF.RR.removeouter$cultivar.genotype)))
  
  
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
  
  #publication bias
  res.tf <- trimfill(res)
  funnel(res.tf, legend=TRUE, cex=1)#绘图
  filled <- data.frame(yi = res.tf$yi, vi = res.tf$vi, fill = res.tf$fill)
  Filled.res<-rma(yi, vi, data=filled)
  mean.RR.filled<-Filled.res$beta[1]
  ci.lb.filled<-Filled.res$ci.lb
  ci.ub.filled<-Filled.res$ci.ub
  
  #RR average calculation
  res<- rma(yi,vi,weights = w.adj,data=DF.RR.removeouter,weighted = T)
  print(res)
  ci.lb<-res$ci.lb
  ci.ub<-res$ci.ub
  mean.RR<-res$beta[1]
  P.value<-res$pval
  
  #output
  output.data<-data.frame(variable,mean.RR,ci.ub,ci.lb,P.value,No.use,No.Rmove,No.species,mean.RR.filled,ci.ub.filled,ci.lb.filled,skewness.value,skewness.Pvalue,kurtosis.value,kurtosis.Pvalue)
  
  return(output.data)
}

#Pvalue fun
Pvalue.fun<-function(value){
  if(value<0.001) {value<-c('***')
  }else if(value<0.01){value<-c('**')
  }else if(value<0.05){value<-c('*')
  }else{value<-c("ns")}
  return(value)
}

#mian prog-------
DF<-import(file.path(main.folder,'Effect of ozone on mesophyll conductance.xlsx'))
DF %>% filter(Accept==1)->DF

gm.output<-meta.fun('gm')
Asat.output<-meta.fun('Asat')
Vcmax.output<-meta.fun('Vcmax')
Jmax.output<-meta.fun('Jmax')
gs.output<-meta.fun('gs')

Out.put<-rbind(gm.output,Asat.output,Vcmax.output,Jmax.output,gs.output)
export(Out.put,file.path(Output.folder,"gasexchang.meta.xlsx"),overwrite=T)


#plot response of gas change-----------
Out.put$mean.RR<-(exp(Out.put$mean.RR)-1)*100
Out.put$ci.ub<-(exp(Out.put$ci.ub)-1)*100
Out.put$ci.lb<-(exp(Out.put$ci.lb)-1)*100
Out.put %>% arrange(mean.RR)->Out.put
Out.put$variable<-factor(Out.put$variable,levels=Out.put$variable)
Out.put$Pvalue.index<-as.vector(lapply(Out.put$P.value,Pvalue.fun))
yaixs<-c(expression(italic(g)[m]),
         expression(italic(A)[sat]),
         expression(italic(V)[cmax]),
         expression(italic(g)[s]),
         expression(italic(J)[max]))




dodge <- position_dodge(width=0.95)
ggplot(Out.put ,aes(x=as.numeric(variable),y=mean.RR))+
  geom_hline(yintercept = 0,linetype="dotted")+
  geom_point(aes(),shape=16,size=3,position = dodge,color="#DC4401")+
  geom_point(aes(),shape=16,size=3,position = dodge,color="#DC4401",show.legend=FALSE)+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub),width=.2,size=0.6,position = dodge,color="#DC4401")+
  expand_limits(y=c(-80,10))+
  scale_y_continuous(breaks=seq(-80,10,15))+
  scale_x_continuous(name = c(),breaks = seq(1,5,1),labels = yaixs,sec.axis = dup_axis(name = c(),labels = mapply(function(value,Pvalue) as.expression(bquote(.(value)*"%"^.(Pvalue))), round(Out.put$mean.RR,1), Out.put$Pvalue.index)))+
  coord_flip()+
  geom_text(aes(x =as.numeric(variable),y = ci.ub+8),label = mapply(function(No,species) as.expression(bquote(.(No)~"("*.(species)*")")), Out.put$No.use, Out.put$No.species),size = 3,color = "#DC4401",fontface = "bold",position = dodge,family="Times")+
  labs(y=as.expression("Change under elevated"~O[3]~"(%)"))+
  labs(x=c())+
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
    axis.ticks.length.x = unit(0.4,"lines"), 
    axis.ticks.length.y.left = unit(0.4,"lines"),
    axis.ticks.length.y.right = unit(0.5,"lines"),
    axis.ticks.y.right = element_line(color='white'),
    axis.ticks = element_line(color='black'),
    axis.line = element_line(colour = "black"),
    #坐标轴文字
    axis.title.y =element_text(colour='black', size=10,face = "bold",vjust = 0),
    axis.title.x=element_text(colour='black', size=10,face = "bold",vjust = 0),
    #轴刻度
    axis.text.y=element_text(colour='black',size=10,vjust = 0.3),
    axis.text.x=element_text(colour = "black",size = 9,angle = 0,hjust = 0.5,vjust = 0),
    #分面标题
    strip.text = element_text(colour = "black",size = 9,face = "bold"),
    strip.background = element_rect(fill=NA,size=0.6,color="black"))
ggsave(file.path(Output.folder,'gasexchang.response.pdf'),width = 10, height = 9,units = "cm",dpi=300)



#plot relationship between Changes of gm and Asat-------------
DF %>% mutate(change=(EO3-AA)/AA*100)->DF.change
DF.change %>% filter(paramter=='Asat')%>% mutate(Asat.change=change) %>% select(-paramter,-Units,-AA,-A.SD,-A.SE,-EO3,-E.SE,-E.SD,-change)->DF.Asat
DF.change%>% filter(paramter=='gm') %>% mutate(gm.change=change) %>% select(-paramter,-Units,-AA,-A.SD,-A.SE,-EO3,-E.SE,-E.SD,-change)->DF.gm
DF.change %>% filter(paramter=='Vcmax')%>% mutate(Vcmax.change=change) %>% select(-paramter,-Units,-AA,-A.SD,-A.SE,-EO3,-E.SE,-E.SD,-change)->DF.Vcmax
DF.change %>% filter(paramter=='Jmax')%>% mutate(Jmax.change=change) %>% select(-paramter,-Units,-AA,-A.SD,-A.SE,-EO3,-E.SE,-E.SD,-change)->DF.Jmax
DF.change %>% filter(paramter=='gs')%>% mutate(gs.change=change) %>% select(-paramter,-Units,-AA,-A.SD,-A.SE,-EO3,-E.SE,-E.SD,-change)->DF.gs
DF.Asat %>% left_join(DF.gm) %>% left_join(DF.Vcmax) %>% left_join(DF.Jmax) %>% left_join(DF.gs)->DF.change.relationship

DF.change.relationship %>% pivot_longer(cols = 24:27)->DF.change.long
ggplot(DF.change.long,aes(y=Asat.change,x = value,color=name))+
  facet_wrap(~name, scales="free_x") +
  geom_point()+
  ggpubr::stat_cor(label.y=max(DF.change.long$Asat.change)*2)+
  ggpubr::stat_regline_equation(label.y=max(DF.change.long$Asat.change))+ 
  geom_smooth(aes(y=Asat.change,x = value,color=name), method=lm, se=T)


#gs---------
plotdata<-DF.change.long %>% filter(name=='gs.change')

fit<-lm(Asat.change~value,data=plotdata)
inter<-round(coef(fit)[[1]],digits = 2)
slope<-round(coef(fit)[[2]],digits = 2)
Pvalue<-round(summary(fit)$coefficients[,4][[2]], 3)
if(Pvalue<0.001) {Pvalue<-c('<0.001')
}else if(Pvalue<0.01){Pvalue<-c('<0.01')
}else if(Pvalue<0.05){Pvalue<-c('<0.05')
}else{Pvalue<-paste("=",round(Pvalue,2))
}
Rsqura<-round(summary(fit)[[8]],2)
#方程
eqn<-as.character(as.expression(
  substitute(""~italic(P)*Pvalue*","~italic(R)^2*"="~Rsqura,list(Pvalue=Pvalue,Rsqura=Rsqura))))

ggplot(plotdata,aes(x=value,y=Asat.change))+
  geom_point(size=4,shape=16,stroke=1,color="#0564C9")+
  geom_point(size=4,shape=16,stroke=1,color="#0564C9",show.legend=FALSE)+
  geom_smooth(method=lm,formula = y~x,color="#3d3d3d")+
  geom_text(aes(x = -60, y =20),label='(a)',size = 4,color = "black",family="Times",font=2,face="bold")+
  labs(x=expression(atop("","Change in "~italic(g)[s]~"("~"%"~")")),family="Times")+
  labs(y=expression(atop("","Change in "~italic(A)[sat]~"("~"%"~")")),family="Times")+
  expand_limits(y=c(-80,20))+
  expand_limits(x=c(-60,40))+
  scale_x_continuous(breaks=seq(-60,40,20))+
  scale_y_continuous(breaks=c(seq(-80,20,20)))+
  #R2
  annotate("text",label=eqn,parse=TRUE,x=-5,y=-80,parse=F, hjust=0,size=5,color="#3d3d3d",family="Times")+
  theme_classic()+
  theme(
    text=element_text(family="Times"),
    plot.title =element_text(hjust = .5,face = "bold",size = 1.5),
    plot.subtitle =element_text(hjust = .6),
    #图例
    legend.background = element_blank(),
    #legend.title=element_blank(),
    legend.key = element_blank(),
    #legend.justification=c(1,1), 
    #legend.position=c(1,0.25),
    legend.key.size = unit(.5, "cm"),
    legend.key.width = unit(0.5,"cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.text = element_text(hjust = 0,colour="black", face = "bold"),
    #图底板折被子
    panel.background = element_rect(fill = 'white'),
    panel.border=element_rect(fill = NA,colour = "black",size = 0.5),
    panel.spacing.y = unit(.5,"cm"),
    panel.grid = element_blank(),
    #坐标轴
    axis.ticks.length = unit(0.4,"lines"), 
    axis.ticks = element_line(color='black'),
    axis.line = element_line(colour = "black",size = 0),
    #坐标轴文字
    axis.title.y=element_text(colour='black', size=15,face = "bold",vjust = 1.5),
    axis.title.x=element_text(colour='black', size=15,face = "bold",vjust = 1.5),
    #轴刻度
    axis.text.y=element_text(colour='black',size=12),
    axis.text.x=element_text(colour = "black",size = 12,angle = 0,hjust = 0.5,vjust = 0),
    #分面标题
    strip.text = element_text(colour = "black",size = 9,face = "bold"),
    strip.background = element_rect(fill=NA,size=0.5,color="black"))
ggsave(file.path(Output.folder,"Asat-gs.pdf"),width = 12, height = 12,units = "cm",dpi=400)

#gm---------
plotdata<-DF.change.long %>% filter(name=='gm.change')

fit<-lm(Asat.change~value,data=plotdata)
inter<-round(coef(fit)[[1]],digits = 2)
slope<-round(coef(fit)[[2]],digits = 2)
Pvalue<-round(summary(fit)$coefficients[,4][[2]], 3)
if(Pvalue<0.001) {Pvalue<-c('<0.001')
}else if(Pvalue<0.01){Pvalue<-c('<0.01')
}else if(Pvalue<0.05){Pvalue<-c('<0.05')
}else{Pvalue<-paste("=",round(Pvalue,2))
}
Rsqura<-round(summary(fit)[[8]],2)
#方程
eqn<-as.character(as.expression(
  substitute(""~italic(P)*Pvalue*","~italic(R)^2*"="~Rsqura,list(Pvalue=Pvalue,Rsqura=Rsqura))))

ggplot(plotdata,aes(x=value,y=Asat.change))+
  geom_point(size=4,shape=16,stroke=1,color="#0564C9")+
  geom_point(size=4,shape=16,stroke=1,color="#0564C9",show.legend=FALSE)+
  geom_smooth(method=lm,formula = y~x,color="#3d3d3d")+
  geom_text(aes(x = -90, y =20),label='(b)',size = 4,color = "black",family="Times",face="bold")+
  labs(x=expression(atop("","Change in "~italic(g)[m]~"("~"%"~")")),family="Times")+
  labs(y=expression(atop("","Change in "~italic(A)[sat]~"("~"%"~")")),family="Times")+
  expand_limits(y=c(-80,20))+
  expand_limits(x=c(-90,40))+
  scale_x_continuous(breaks=seq(-90,40,20))+
  scale_y_continuous(breaks=c(seq(-80,20,20)))+
  #R2
  annotate("text",label=eqn,parse=TRUE,x=-15,y=-80,parse=F, hjust=0,size=5,color="#3d3d3d",family="Times")+
  theme_classic()+
  theme(
    text=element_text(family="Times"),
    plot.title =element_text(hjust = .5,face = "bold",size = 1.5),
    plot.subtitle =element_text(hjust = .6),
    #图例
    legend.background = element_blank(),
    #legend.title=element_blank(),
    legend.key = element_blank(),
    #legend.justification=c(1,1), 
    #legend.position=c(1,0.25),
    legend.key.size = unit(.5, "cm"),
    legend.key.width = unit(0.5,"cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.text = element_text(hjust = 0,colour="black", face = "bold"),
    #图底板折被子
    panel.background = element_rect(fill = 'white'),
    panel.border=element_rect(fill = NA,colour = "black",size = 0.5),
    panel.spacing.y = unit(.5,"cm"),
    panel.grid = element_blank(),
    #坐标轴
    axis.ticks.length = unit(0.4,"lines"), 
    axis.ticks = element_line(color='black'),
    axis.line = element_line(colour = "black",size = 0),
    #坐标轴文字
    axis.title.y=element_text(colour='black', size=15,face = "bold",vjust = 1.5),
    axis.title.x=element_text(colour='black', size=15,face = "bold",vjust = 1.5),
    #轴刻度
    axis.text.y=element_text(colour='black',size=12),
    axis.text.x=element_text(colour = "black",size = 12,angle = 0,hjust = 0.5,vjust = 0),
    #分面标题
    strip.text = element_text(colour = "black",size = 9,face = "bold"),
    strip.background = element_rect(fill=NA,size=0.5,color="black"))
ggsave(file.path(Output.folder,"Asat-gm.pdf"),width = 12, height = 12,units = "cm",dpi=400)


#Vcmax---------
plotdata<-DF.change.long %>% filter(name=='Vcmax.change')

fit<-lm(Asat.change~value,data=plotdata)
inter<-round(coef(fit)[[1]],digits = 2)
slope<-round(coef(fit)[[2]],digits = 2)
Pvalue<-round(summary(fit)$coefficients[,4][[2]], 3)
if(Pvalue<0.001) {Pvalue<-c('<0.001')
}else if(Pvalue<0.01){Pvalue<-c('<0.01')
}else if(Pvalue<0.05){Pvalue<-c('<0.05')
}else{Pvalue<-paste("=",round(Pvalue,2))
}
Rsqura<-round(summary(fit)[[8]],2)
#方程
eqn<-as.character(as.expression(
  substitute(""~italic(P)*Pvalue*","~italic(R)^2*"="~Rsqura,list(Pvalue=Pvalue,Rsqura=Rsqura))))

ggplot(plotdata,aes(x=value,y=Asat.change))+
  geom_point(size=4,shape=16,stroke=1,color="#0564C9")+
  geom_point(size=4,shape=16,stroke=1,color="#0564C9",show.legend=FALSE)+
  geom_smooth(method=lm,formula = y~x,color="#3d3d3d")+
  geom_text(aes(x = -80, y =20),label='(c)',size = 4,color = "black",family="Times",face="bold")+
  labs(x=expression(atop("","Change in "~italic(V)[cmax]~"("~"%"~")")),family="Times")+
  labs(y=expression(atop("","Change in "~italic(A)[sat]~"("~"%"~")")),family="Times")+
  expand_limits(y=c(-80,20))+
  expand_limits(x=c(-80,20))+
  scale_x_continuous(breaks=seq(-80,20,20))+
  scale_y_continuous(breaks=c(seq(-80,20,20)))+
  #R2
  annotate("text",label=eqn,parse=TRUE,x=-25,y=-80,parse=F, hjust=0,size=5,color="#3d3d3d",family="Times")+
  theme_classic()+
  theme(
    text=element_text(family="Times"),
    plot.title =element_text(hjust = .5,face = "bold",size = 1.5),
    plot.subtitle =element_text(hjust = .6),
    #图例
    legend.background = element_blank(),
    #legend.title=element_blank(),
    legend.key = element_blank(),
    #legend.justification=c(1,1), 
    #legend.position=c(1,0.25),
    legend.key.size = unit(.5, "cm"),
    legend.key.width = unit(0.5,"cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.text = element_text(hjust = 0,colour="black", face = "bold"),
    #图底板折被子
    panel.background = element_rect(fill = 'white'),
    panel.border=element_rect(fill = NA,colour = "black",size = 0.5),
    panel.spacing.y = unit(.5,"cm"),
    panel.grid = element_blank(),
    #坐标轴
    axis.ticks.length = unit(0.4,"lines"), 
    axis.ticks = element_line(color='black'),
    axis.line = element_line(colour = "black",size = 0),
    #坐标轴文字
    axis.title.y=element_text(colour='black', size=15,face = "bold",vjust = 1.5),
    axis.title.x=element_text(colour='black', size=15,face = "bold",vjust = 1.5),
    #轴刻度
    axis.text.y=element_text(colour='black',size=12),
    axis.text.x=element_text(colour = "black",size = 12,angle = 0,hjust = 0.5,vjust = 0),
    #分面标题
    strip.text = element_text(colour = "black",size = 9,face = "bold"),
    strip.background = element_rect(fill=NA,size=0.5,color="black"))
ggsave(file.path(Output.folder,"Asat-Vcmax.pdf"),width = 12, height = 12,units = "cm",dpi=400)


#Jmax---------
plotdata<-DF.change.long %>% filter(name=='Jmax.change')

fit<-lm(Asat.change~value,data=plotdata)
inter<-round(coef(fit)[[1]],digits = 2)
slope<-round(coef(fit)[[2]],digits = 2)
Pvalue<-round(summary(fit)$coefficients[,4][[2]], 3)
if(Pvalue<0.001) {Pvalue<-c('<0.001')
}else if(Pvalue<0.01){Pvalue<-c('<0.01')
}else if(Pvalue<0.05){Pvalue<-c('<0.05')
}else{Pvalue<-paste("=",round(Pvalue,2))
}
Rsqura<-round(summary(fit)[[8]],2)
#方程
eqn<-as.character(as.expression(
  substitute(""~italic(P)*Pvalue*","~italic(R)^2*"="~Rsqura,list(Pvalue=Pvalue,Rsqura=Rsqura))))
ggplot(plotdata,aes(x=value,y=Asat.change))+
  geom_point(size=4,shape=16,stroke=1,color="#0564C9")+
  geom_point(size=4,shape=16,stroke=1,color="#0564C9",show.legend=FALSE)+
  geom_smooth(method=lm,formula = y~x,color="#3d3d3d")+
  geom_text(aes(x = -70, y =20),label='(d)',size = 4,color = "black",family="Times",face="bold")+
  labs(x=expression(atop("","Change in "~italic(J)[max]~"("~"%"~")")),family="Times")+
  labs(y=expression(atop("","Change in "~italic(A)[sat]~"("~"%"~")")),family="Times")+
  expand_limits(y=c(-80,20))+
  expand_limits(x=c(-70,20))+
  scale_x_continuous(breaks=seq(-70,20,20))+
  scale_y_continuous(breaks=c(seq(-80,20,20)))+
  #R2
  annotate("text",label=eqn,parse=TRUE,x=-18,y=-85,parse=F, hjust=0,size=5,color="#3d3d3d",family="Times")+
  theme_classic()+
  theme(
    text=element_text(family="Times"),
    plot.title =element_text(hjust = .5,face = "bold",size = 1.5),
    plot.subtitle =element_text(hjust = .6),
    #图例
    legend.background = element_blank(),
    #legend.title=element_blank(),
    legend.key = element_blank(),
    #legend.justification=c(1,1), 
    #legend.position=c(1,0.25),
    legend.key.size = unit(.5, "cm"),
    legend.key.width = unit(0.5,"cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.text = element_text(hjust = 0,colour="black", face = "bold"),
    #图底板折被子
    panel.background = element_rect(fill = 'white'),
    panel.border=element_rect(fill = NA,colour = "black",size = 0.5),
    panel.spacing.y = unit(.5,"cm"),
    panel.grid = element_blank(),
    #坐标轴
    axis.ticks.length = unit(0.4,"lines"), 
    axis.ticks = element_line(color='black'),
    axis.line = element_line(colour = "black",size = 0),
    #坐标轴文字
    axis.title.y=element_text(colour='black', size=15,face = "bold",vjust = 1.5),
    axis.title.x=element_text(colour='black', size=15,face = "bold",vjust = 1.5),
    #轴刻度
    axis.text.y=element_text(colour='black',size=12),
    axis.text.x=element_text(colour = "black",size = 12,angle = 0,hjust = 0.5,vjust = 0),
    #分面标题
    strip.text = element_text(colour = "black",size = 9,face = "bold"),
    strip.background = element_rect(fill=NA,size=0.5,color="black"))
ggsave(file.path(Output.folder,"Asat-Jmax.pdf"),width = 12, height = 12,units = "cm",dpi=400)





