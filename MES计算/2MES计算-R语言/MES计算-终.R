###1、载入相关包
rm(list=ls())
library(rmgarch)
library(rugarch)
library(tidyverse)
library(readxl)
library(ggplot2)
library(zoo)
library(xts)
library(SparseM)
library(plotrix)

#2、DCC-GARCH模型估计
garch11.spec = ugarchspec(mean.model = list(armaOrder = c(1,0)), 
                          variance.model = list(garchOrder = c(0,2),
                                                model = "gjrGARCH"), 
                          distribution.model = "std") 

dcc.garch11.spec = dccspec(uspec = multispec( replicate(2, garch11.spec) ), 
                           dccOrder = c(1,1), 
                           distribution = "mvnorm")

fct_MES=function(ret,c,ht_m,ht_i,rho){
  em=as.numeric(unlist(ret[,1]))/ht_m                        #  market first column
  xi=(as.numeric(unlist(ret[,2]))/ht_i-rho*em)/sqrt(1-rho^2) #  asset second column
  bwd=nrow(ret)^(-0.2)                  #  Scaillet's bwd p21
  K1=sum(em*pnorm((c/ht_m-em)/bwd))/sum(pnorm(c/ht_m-em)/bwd)
  K2=sum(xi*pnorm((c/ht_m-em)/bwd))/sum(pnorm(c/ht_m-em)/bwd)
  MES = (ht_i*rho*K1) + (ht_i*sqrt(1-rho^2)*K2)
  return(-MES)
}
#3、导入数据并计算1343家企业的MES
data1<-read_excel("算MES数据.xlsx")
data<-data1
data
MESData<-as.data.frame(data)
for (i in 3:ncol(data)){

	ret=xts(data[,c(2,i)], order.by=data$time)
	dcc.fit = dccfit(dcc.garch11.spec, data = ret)

	ht_m=sqrt(dcc.fit@mfit$H[1,1,])
	ht_i=sqrt(dcc.fit@mfit$H[2,2,])
	rho=rcor(dcc.fit)[1,2,]
	alpha = 0.05 
	k = 0.08 
	Asset_VaR = ht_i*quantile(as.numeric(unlist(ret[,2]))/ht_i,alpha)
	Market_VaR = ht_m*quantile(as.numeric(unlist(ret[,1]))/ht_m,alpha)
	Beta = rho*ht_i/ht_m
	c = quantile(as.numeric(unlist(ret[,1])),alpha)
	MES=fct_MES(ret,c,ht_m,ht_i,rho)
	head(MES,10)
	aa=as.data.frame(MES)
	MESData[,i]=aa
}

head(data,10)
head(MESData,10)
write.csv(MESData,"MESData.csv")

