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
DF %>% filter(Accept==1)->DF


#abaxial epidemis thickness---------
variable<-"abaxial epidemis thickness"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)


#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)
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

abaxial.output<-output.data


#adaxial epidemis thickness---------
variable<-"adaxial epidemis thickness"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)


#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)
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

adaxial.output<-output.data



#leaf thickness---------
variable<-"Leaf thickness"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)


#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR

#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)
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

leafthickness.output<-output.data


#Spongy thickness---------
variable<-"Spongy thickness"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)


#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)
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

spongythickness.output<-output.data


#Palisade thickness---------
variable<-"Palisade thickness"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)
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

palisadethickness.output<-output.data


#Intercellular space of mesophyll---------
variable<-"Intercellular space of mesophyll"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)

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

Intercellular.output<-output.data



#cell wall thickness---------
variable<-"cell wall thickness"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)

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

cellwallthickness.output<-output.data


#chloroplast size---------
variable<-"chloroplast size"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)

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

chlsize.output<-output.data




#chloroplast width---------
variable<-"chloroplast width"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)

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

chlwidth.output<-output.data



#chloroplast length---------
variable<-"chloroplast length"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)

#outliers test
DF.RR$standresidual<-rstandard(res)$z
DF.RR<-cbind(DF.RR,influence(res)[[1]])

#leave1out test 
remove.test<-leave1out(res)

#filter standardized residual value >3
DF.RR %>% filter(abs(standresidual)>2)->DF.RR.rmove
No.Rmove<-length(DF.RR.rmove$Accept)
print(DF.RR.rmove)
DF.RR %>% filter(abs(standresidual)<2)->DF.RR.removeouter
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

chllength.output<-output.data

#starch size---------
variable<-"Starch size"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)

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

starchsize.output<-output.data



#starch length---------
variable<-"starch length"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)

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

starchlength.output<-output.data



#starch width---------
variable<-"starch width"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)

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

starchwidth.output<-output.data



#starch/chloroplast---------
variable<-"Starch/chloroplast"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)

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

starch.chl.output<-output.data



#chlorplast/cell---------
variable<-"chloroplast/cell"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)

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

chl.cell.output<-output.data


#chlorplast/mesophyll---------
variable<-"chloroplast/mesophyll"
DF %>% filter(Accept==1&paramter==variable)->DF.analysis
# calculate mean differences and corresponding sampling variances (RR and v)
DF.RR <- escalc(measure="ROM",m1i=EO3, sd1i=E.SD, n1i=E.N, m2i=AA, sd2i=A.SD, n2i=A.N, data=DF.analysis)
hist(DF.RR$yi, prob = TRUE)
lines(density(DF.RR$yi), col = 2, lwd = 3)

#independence and weights
DF.RR %>%filter(!is.na(vi)) %>%  mutate(w=1/vi,w.adj=w/Independence,RR.adj=yi*w.adj)->DF.RR
#fit mixed-effects models
res<- rma(yi,vi,weights = w.adj,data=DF.RR,weighted = T)
print(res)

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

chl.mesophy.output<-output.data

#Plot leaf anatomy----------
anatomy.plot<-rbind(abaxial.output,adaxial.output,leafthickness.output,spongythickness.output,palisadethickness.output,Intercellular.output)
export(anatomy.plot,file.path(Output.folder,"leafanatomy.meta.xlsx"),overwrite=T)

#plot response
anatomy.plot$mean.RR<-(exp(anatomy.plot$mean.RR)-1)*100
anatomy.plot$ci.ub<-(exp(anatomy.plot$ci.ub)-1)*100
anatomy.plot$ci.lb<-(exp(anatomy.plot$ci.lb)-1)*100
anatomy.plot %>% arrange(mean.RR)->anatomy.plot
anatomy.plot$variable<-factor(anatomy.plot$variable,levels=anatomy.plot$variable)
anatomy.plot$Pvalue.index<-as.vector(lapply(anatomy.plot$P.value,Pvalue.fun))
anatomy.plot$P.value.colour<-ifelse (anatomy.plot$P.value>0.05,"a","b")


yaixs<-c(expression("Leaf thickness"),
         expression(italic(T)[palisade]),
         expression(italic(f)[ias]),
         expression(italic(T)[spongy]),
         expression(italic(T)['adaxial epidermis']),
         expression(italic(T)['adaxial epidermis']))

dodge <- position_dodge(width=0.95)
ggplot(anatomy.plot ,aes(x=as.numeric(variable),y=mean.RR))+
  geom_hline(yintercept = 0,linetype="dotted")+
  geom_point(aes(color=P.value.colour),shape=16,size=3,position = dodge)+
  geom_point(aes(color=P.value.colour),shape=16,size=3,position = dodge,show.legend=FALSE)+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub,color=P.value.colour),width=.2,size=0.6,position = dodge)+
  expand_limits(y=c(-20,20))+
  scale_y_continuous(breaks=seq(-20,20,5))+
  scale_x_continuous(name = c(),breaks = seq(1,6,1),labels = yaixs,sec.axis = dup_axis(name = c(),labels = mapply(function(value,Pvalue) as.expression(bquote(.(value)*"%"^.(Pvalue))), round(anatomy.plot$mean.RR,1), anatomy.plot$Pvalue.index)))+
  coord_flip()+
  geom_text(aes(x = 6.5, y =-19),label='(a)',size = 4,color = "black",family="Times",face="bold")+
  
  geom_text(aes(x =as.numeric(variable),y = ci.ub+5,color=P.value.colour),label = mapply(function(No,species) as.expression(bquote(.(No)~"("*.(species)*")")), anatomy.plot$No.use, anatomy.plot$No.species),size = 3,fontface = "bold",position = dodge,family="Times")+
  labs(y=as.expression("Change under elevated"~O[3]~"(%)"))+
  labs(x=c())+
  scale_color_manual(values=c("#0564C9","#DC4401"))+
  theme_classic()+
  theme(
    text=element_text(family="Times"),
    plot.title =element_text(hjust = .5,face = "bold",size = 1.5),
    plot.subtitle =element_text(hjust = .6),
    #图例
    legend.background = element_blank(),
    legend.title=element_blank(),
    #legend.justification=c(1,1), 
    legend.position='none',
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
ggsave(file.path(Output.folder,'leafanatomy.response.pdf'),width = 10, height = 9,units = "cm",dpi=300)

#Plot leaf structural traits----------
anatomy.plot<-rbind(cellwallthickness.output,chlsize.output,chlwidth.output,chllength.output,chl.mesophy.output,starchsize.output,starchlength.output,starchwidth.output,starch.chl.output)
export(anatomy.plot,file.path(Output.folder,"leafstructual.meta.xlsx"),overwrite=T)

#plot response
anatomy.plot$mean.RR<-(exp(anatomy.plot$mean.RR)-1)*100
anatomy.plot$ci.ub<-(exp(anatomy.plot$ci.ub)-1)*100
anatomy.plot$ci.lb<-(exp(anatomy.plot$ci.lb)-1)*100
anatomy.plot %>% arrange(mean.RR)->anatomy.plot
anatomy.plot$variable<-factor(anatomy.plot$variable,levels=anatomy.plot$variable)
anatomy.plot$Pvalue.index<-as.vector(lapply(anatomy.plot$P.value,Pvalue.fun))
anatomy.plot$P.value.colour<-ifelse (anatomy.plot$P.value>0.05,"a","b")

yaixs<-c(expression(italic(Area)[starch]),
         expression(italic(Area)[chl]),
         expression(italic(T)[starch]),
         expression(italic(Area)[chloroplast/mesophyll]),
         expression(italic(T)[chl]),
         expression(italic(L)[starch]),
         expression(italic(Area)[starch/chloroplast]),
         expression(italic(L)[chl]),
         expression(italic(T)[cw]))

dodge <- position_dodge(width=0.95)
ggplot(anatomy.plot ,aes(x=as.numeric(variable),y=mean.RR))+
  geom_hline(yintercept = 0,linetype="dotted")+
  geom_point(aes(color=P.value.colour),shape=16,size=3,position = dodge)+
  geom_point(aes(color=P.value.colour),shape=16,size=3,position = dodge,show.legend=FALSE)+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub,color=P.value.colour),width=.2,size=0.6,position = dodge)+
  expand_limits(y=c(-80,40))+
  scale_y_continuous(breaks=seq(-80,40,20))+
  scale_x_continuous(name = c(),breaks = seq(1,9,1),labels = yaixs,sec.axis = dup_axis(name = c(),labels = mapply(function(value,Pvalue) as.expression(bquote(.(value)*"%"^.(Pvalue))), round(anatomy.plot$mean.RR,1), anatomy.plot$Pvalue.index)))+
  coord_flip()+
  geom_text(aes(x = 9.5, y =-75),label='(b)',size = 4,color = "black",family="Times",face="bold")+
  geom_text(aes(x =as.numeric(variable),y = ci.ub+13,color=P.value.colour),label = mapply(function(No,species) as.expression(bquote(.(No)~"("*.(species)*")")), anatomy.plot$No.use, anatomy.plot$No.species),size = 3,fontface = "bold",position = dodge,family="Times")+
  labs(y=as.expression("Change under elevated"~O[3]~"(%)"))+
  labs(x=c())+
  scale_color_manual(values=c("#0564C9","#DC4401"))+
  theme_classic()+
  theme(
    text=element_text(family="Times"),
    plot.title =element_text(hjust = .5,face = "bold",size = 1.5),
    plot.subtitle =element_text(hjust = .6),
    #图例
    legend.background = element_blank(),
    legend.title=element_blank(),
    #legend.justification=c(1,1), 
    legend.position='none',
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
ggsave(file.path(Output.folder,'leafstructural.response.pdf'),width = 11, height = 9,units = "cm",dpi=300)

