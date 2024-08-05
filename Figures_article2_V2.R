
######################################################################################################
######################################################################################################
# Forecast probability of spatiotemporal records : France 
######################################################################################################
######################################################################################################


#############################################
###  packages
#############################################
# cleaning the workbench
rm(list=ls(all=TRUE))
set.seed(2018)
library(MASS)
library(latex2exp)
library(mev)
library(evd)
library(RColorBrewer)

#############################################
###  Data
#############################################
paris_information <- readRDS("paris_information.rds")
paris_observations <- readRDS("paris_observations.rds")
south_information <- readRDS("south_information.rds")
south_observations <- readRDS("south_observations.rds")
britany_information <- readRDS("britany_information.rds")
britany_observations <- readRDS("britany_observations.rds")
Used.Data.Information <- rbind(britany_information,south_information,paris_information) # bind the 3 information data.frames into 1 


#############################################
###  Figure 2 Map 
#############################################
colorRegion <- c("blue","blue","blue","blue","blue","red","red","red","red","red","green3","green3","green3","green3","green3")
Used.Data.Information <- cbind(Used.Data.Information,colorRegion)
# names on the stations of interest 
s0.Information <- filter(Used.Data.Information, name %in% c("Brest","Paris-Montsouris","Salon de Provence"))
head(s0.Information)
colors0 <- c("deepskyblue","orange","chartreuse")
colors_cluster <- c("blue","green4","red3")
name_cluster <- c("Britany","Paris region","South of France")
france_map <- map_data("france")
mapcluster <- ggplot() +
  geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey40', size = 0.6) +
  coord_quickmap() +
  geom_point(data = Used.Data.Information, aes(x = longitude, y = latitude, fill = factor(colorRegion)), color = "grey1", pch = 21, size = 5) +
  geom_point(data =s0.Information[1,], aes(x = longitude, y = latitude),fill = colors0[1], color = "black", pch = 21, size = 4)+
  geom_point(data =s0.Information[2,], aes(x = longitude, y = latitude),fill = colors0[2], color = "black", pch = 21, size = 4)+
  geom_point(data =s0.Information[3,], aes(x = longitude, y = latitude),fill = colors0[3], color = "black", pch = 21, size = 4)+
  scale_fill_manual(values = colors_cluster, name = " ", labels = name_cluster) +
  theme_minimal(base_size = 20) +
  labs(x = "Longitude", y = "Latitude", title = "Regions of study")+
  theme(legend.position = c(0.27, 0.2))+
  theme(axis.title.y = element_text(size = rel(0.9), angle = 90))+
  theme(axis.title.x = element_text(size = rel(0.9), angle = 00))+
  theme(legend.text=element_text(size=23),legend.title = element_blank())+
  theme(legend.background = element_rect(fill="grey92", 
                                         size=0.7, linetype="solid"))
print(mapcluster)


#############################################
###  Figure 3 - Observation time series
#############################################

years<-c(1960:2022)
all_observations <- cbind(britany_observations,paris_observations,south_observations)
head(all_observations)

BZH<-britany_observations
iBZH<-3   # Brest
iSUD<-1  # Salon de provence
iIDF<-1  # Mont-Souris
SUD<-south_observations
IDF<-paris_observations
#

colBZH<-brewer.pal(5, name ="Blues")
colSUD<-brewer.pal(5, name ="YlOrRd")
colIDF<-brewer.pal(5, name ="Greens")

records<-function(M,id, years){
  y.max<-NA
  I.max<-NA
  for(i in 1:(nrow(M)-1)){
    Max<-max(M[1:i,],na.rm = T)
    if(!is.na(M[i+1,id])){
      if(M[i+1,id] >=Max){
        y.max<-c(y.max,years[i+1])
        I.max<-c(I.max,M[i+1,id])
      }
      
    }
  }
  out<-cbind(y.max,I.max)
  out<-out[-1,]
  return(out)
}

recordBZH<-records(BZH,iBZH,years)
recordSUD<-records(SUD,iSUD,years)
recordIDF<-records(IDF,iIDF,years)

y.lim<-range(all_observations, na.rm=T)

def.par <- par(no.readonly = TRUE)
nf<-layout(matrix(1:3, nrow=3))
layout.show(nf)

yy<-c(1960,1970,1980,1990,2000,2010)


def.par <- par(no.readonly = TRUE)
nf<-layout(matrix(1:3, nrow=3))
layout.show(nf)

#############################################
par(mar=c(0, 5, 4, 4))
plot(years,britany_observations[,1], cex.lab=1.7,xlab="Years", ylab="Tmax (°C)", type="n",ylim=y.lim, axes=F)
axis(side = 3,cex.axis=2);axis(side=2,cex.axis=1.5)
for(i in 1:length(yy)) abline(v=yy[i], col="grey", lty=1)
box()
for(i in 1:5){
  lines(years,britany_observations[,i],col=colBZH[i], lwd=2, lty=2)
}
legend("topleft",legend = "Brittany",bg = "white", text.col = colBZH[5], cex=2)
lines(years,britany_observations[,iBZH],col=colBZH[5], lwd=3)
points(recordBZH,col=colBZH[5], cex=3,pch=16)
#############################################
par(mar=c(0, 5, 0, 4))
plot(years,paris_observations[,iIDF], cex.lab=1.7,xlab="", ylab="Tmax (°C)", type="n",ylim=y.lim, axes=F)
axis(side = 4,cex.axis=1.5)
for(i in 1:length(yy)) abline(v=yy[i], col="grey", lty=1)
box()
for(i in 1:5){
  lines(years,paris_observations[,i],col=colIDF[i], lwd=2, lty=2)
}
legend("topleft",legend = "Paris region",bg = "white", text.col = colIDF[5], cex=2)
lines(years,paris_observations[,iIDF],col=colIDF[5], lwd=3)
points(recordIDF,col=colIDF[5], cex=3,pch=16)
#############################################

par(mar=c(4, 5, 0, 4))
plot(years,south_observations[,1],cex.lab=1.7, xlab="", ylab="Tmax (°C)", type="n",ylim=y.lim, axes=F)
axis(side = 1,cex.axis=2);axis(side=2,cex.axis=1.5)
for(i in 1:length(yy)) abline(v=yy[i], col="grey", lty=1)
box()
for(i in 1:5){
  lines(years,south_observations[,i],col=colSUD[i], lwd=2, lty=2)
}
legend("topleft",legend = "South of France",bg = "white", text.col = colSUD[5], cex=2)
lines(years,south_observations[,iSUD],col=colSUD[5], lwd=3)
points(recordSUD,col=colSUD[5], cex=3,pch=16)


#############################################
###  Record probability
#############################################

source("/Users/pgonzale/article2_functions.R")

end.t <- dim(britany_observations)[1]
par(mfrow=c(2,2),mar=c(5,5,4,2))
h1 <- 25
h2 <- 25
GmaxHatY_britany <- GmaxY.estimation.function.step1(britany_observations, bandwidth = h1)$Array.GmaxHatY
E.GmaxY_britany <- EGmaxY.estimation.function.step2(GmaxHatY_britany, bandwidth = h2)$Array.EGmaxHatY
GmaxHatY_paris<- GmaxY.estimation.function.step1(paris_observations, bandwidth = h1)$Array.GmaxHatY
E.GmaxY_paris <- EGmaxY.estimation.function.step2(GmaxHatY_paris, bandwidth = h2)$Array.EGmaxHatY
GmaxHatY_south <- GmaxY.estimation.function.step1(south_observations, bandwidth = h1)$Array.GmaxHatY
E.GmaxY_south <- EGmaxY.estimation.function.step2(GmaxHatY_south, bandwidth = h2)$Array.EGmaxHatY

# Britany
#################################
s0 <- 3 # Brest
ProbaFinal_britany <- Get.P.spatiotemp(Data = britany_observations,
                                       EGmaxY = E.GmaxY_britany,
                                       site.interest=s0 ,
                                       add.iid=TRUE,
                                       add.indep=TRUE,
                                       add.stat=TRUE,
                                       h1 = h1,
                                       par1=1,par2=1,
                                       use.timesteps = 0)
Precord.so.Model_britany <- ProbaFinal_britany$P.record
Precord.so.iid_britany<- ProbaFinal_britany$P.record.iid
Precord.so.indep_britany <- ProbaFinal_britany$P.record.indep
Precord.so.stat_britany <- ProbaFinal_britany$P.record.stat
Pickands.model_britany <- ProbaFinal_britany$A
Pickands.stat_britany <- ProbaFinal_britany$Astat

# Paris region
#################################
s0 <- 1 # Paris-montsouris
ProbaFinal_paris <- Get.P.spatiotemp(Data = paris_observations,
                                     EGmaxY = E.GmaxY_paris,
                                     site.interest=s0 ,
                                     add.iid=TRUE,
                                     add.indep=TRUE,
                                     add.stat=TRUE,
                                     h1 = h1,
                                     par1=1,par2=1,
                                     use.timesteps = 0)
Precord.so.Model_paris <- ProbaFinal_paris$P.record
Precord.so.iid_paris<- ProbaFinal_paris$P.record.iid
Precord.so.indep_paris <- ProbaFinal_paris$P.record.indep
Precord.so.stat_paris <- ProbaFinal_paris$P.record.stat
Pickands.model_paris <- ProbaFinal_paris$A
Pickands.stat_paris <- ProbaFinal_paris$Astat

# South of France
#################################
s0 <- 1 # Salon de Provence
ProbaFinal_south <- Get.P.spatiotemp(Data = south_observations,
                                     EGmaxY = E.GmaxY_south,
                                     site.interest=s0 ,
                                     add.iid=TRUE,
                                     add.indep=TRUE,
                                     add.stat=TRUE,
                                     h1 = h1,
                                     par1=1,par2=1,
                                     use.timesteps = 0)
Precord.so.Model_south <- ProbaFinal_south$P.record
Precord.so.iid_south<- ProbaFinal_south$P.record.iid
Precord.so.indep_south <- ProbaFinal_south$P.record.indep
Precord.so.stat_south <- ProbaFinal_south$P.record.stat
Pickands.model_south <- ProbaFinal_south$A
Pickands.stat_south <- ProbaFinal_south$Astat

#############################################
###  Record probability ratio
#############################################
S <- 5 
nb.timesteps <- dim(britany_observations)[1]
timesstesp <- c(1:(nb.timesteps-1))

startYear <- 1960
nb.years <- length(Precord.so.Model_britany)
plotaxis.vec <- c(1,seq(10,50,10),nb.years)
V.alph <-5
lamtt <- timesstesp*V.alph
pidV <- 1/(1+(timesstesp*lamtt))

ylirang <- c(Precord.so.Model_britany/pidV,Precord.so.Model_paris/pidV,Precord.so.Model_south/pidV)
par(mfrow=c(1,1),mar=c(5,5.5,4,2)) 
plot(c(1:62),Precord.so.Model_south/pidV,col="red",type="l",lty=1,ylim=range(ylirang),xaxt = "n",cex.main=1.5,
     xlab= "Year",ylab=TeX(r'($\widehat{P}_{\tau}(s_0)/P_{\alpha}$)'), lwd=2.5, main= "Signal on record probability ",
     cex.lab=1.8,cex.axis=2)
abline(h=1/5,lwd=4,col='grey60', lty=2)
abline(h=1/V.alph,lwd=4,col='grey60', lty=2)
lines(c(1:62),Precord.so.Model_south/pidV,lwd=5,col="red")
lines(c(1:62),Precord.so.Model_paris/pidV,lwd=5,col="green4")
lines(c(1:62),Precord.so.Model_britany/pidV,lwd=5,col="blue")
axis(1, at= plotaxis.vec, labels= plotaxis.vec + (startYear),cex.axis=1.2)
legend('topleft',legend=c("Brest", "Paris-Montsouris" , "Salon the Provence"),
       col=c("blue",'green4',"red"),lty= 1,lwd=5, bg="transparent",cex=1.5)


#############################################
###  Forecast
#############################################

# Britany
#############################################
HW_2toT_britany <- stats::HoltWinters(ts(Precord.so.Model_britany, frequency = 1), gamma = FALSE)
HW_2toT.fit_britany  <- HW_2toT_britany$fitted[,1]
HW_Tnext.fit_britany <- predict(HW_2toT_britany, n.ahead = 1, prediction.interval = TRUE)
HW_Pupcoming_britany <- as.numeric(HW_Tnext.fit_britany[,1])

# Paris
#############################################
HW_2toT_paris <- stats::HoltWinters(ts(Precord.so.Model_paris, frequency = 1), gamma = FALSE)
HW_2toT.fit_paris  <- HW_2toT_paris$fitted[,1]
HW_Tnext.fit_paris <- predict(HW_2toT_paris, n.ahead = 1, prediction.interval = TRUE)
HW_Pupcoming_paris <- as.numeric(HW_Tnext.fit_paris[,1])

# South
#############################################
HW_2toT_south <- stats::HoltWinters(ts(Precord.so.Model_south, frequency = 1), gamma = FALSE)
HW_2toT.fit_south  <- HW_2toT_south$fitted[,1]
HW_Tnext.fit_south <- predict(HW_2toT_south, n.ahead = 1, prediction.interval = TRUE)
HW_Pupcoming_south <- as.numeric(HW_Tnext.fit_south[,1])

#############################################
###  Plot 4 - Record probability ratio up to T+1
#############################################
startYear <- 1960
nb.years <- length(Precord.so.Model_britany)

V.alph <-1
lamtt <- timesstesp*V.alph
pidV <- 1/(1+(timesstesp*lamtt))
plotaxis.vec <- c(1,seq(10,50,10),(nb.years)+1)

timesstespTnext <- c(1:(dim(paris_observations)[1]))
lamttTnext <- timesstespTnext*V.alph
pidVTnext <- 1/(1+(timesstespTnext*lamttTnext))

PrecordTnext_britany <- c(Precord.so.Model_britany,HW_Pupcoming_britany)
PrecordTnext_paris <- c(Precord.so.Model_paris,HW_Pupcoming_paris)
PrecordTnext_south <- c(Precord.so.Model_south,HW_Pupcoming_south)

start.t <- 2
end.t <- dim(britany_observations)[1]
ylirangTnext <- c(PrecordTnext_britany/pidVTnext,PrecordTnext_paris/pidVTnext,PrecordTnext_south/pidVTnext)
par(mfrow=c(1,1),mar=c(5,5.5,4,2)) 
plot(c(1:end.t),PrecordTnext_britany/pidVTnext,col="red",type="l",lty=1,ylim=range(ylirangTnext),xaxt = "n",cex.main=1.5,
     xlab= "Year",ylab=TeX(r'($\widehat{P}_{\tau}(s_0)/P_{\alpha}$)'), lwd=2.5, main= "Signal on record probability ",
     cex.lab=1.8,cex.axis=2)
abline(h=1/5,lwd=4,col='grey60', lty=2)
lines(c(1:end.t),PrecordTnext_south/pidVTnext,lwd=5,col="red")
lines(c(1:end.t),PrecordTnext_paris/pidVTnext,lwd=5,col="green4")
lines(c(1:end.t),PrecordTnext_britany/pidVTnext,lwd=5,col="blue")
abline(v=(nb.years)+1,lwd=2,col='grey60', lty=2)
points(end.t,(PrecordTnext_britany/pidVTnext)[end.t],pch=21,col="black",bg="blue",cex=2,lwd=2)
points(end.t,(PrecordTnext_paris/pidVTnext)[end.t],pch=21,col="black",bg="green4",cex=2,lwd=2)
points(end.t,(PrecordTnext_south/pidVTnext)[end.t],pch=21,col="black",bg="red",cex=2,lwd=2)

axis(1, at= plotaxis.vec, labels= plotaxis.vec + (startYear),cex.axis=1.2)
legend('topleft',legend=c("Brest", "Paris-Montsouris" , "Salon the Provence"),
       col=c("blue",'green4',"red"),lty= 1,lwd=5, bg="transparent",cex=1.5)

#############################################
###  Figure 1 - Schema Record probability
#############################################
S <- 5
nb.timesteps <- 2022-1959 
timesstesp <- c(1:(nb.timesteps-1))
aalpha <- 0.3
S <- 5
# example
######################################################################################################
xi <- -0.2
support <- 40
mu.0 <- 26
mu.vec<- mu.0 + ((c(timesstesp,63)^2)*0.002) + ((c(timesstesp,63))*0.001)
sigma.vec <- -xi*(support-mu.vec) 
sigmaStar.vec <- (cumsum(sigma.vec^(1/xi)))^xi
lam.vec <- (sigmaStar.vec[-63]/sigma.vec[-1])^(1/xi)
lamV_ex1 <- (S * (lam.vec^(aalpha))^(1/aalpha))
Pmodel_ex1 <- 1/(1+(timesstesp * lamV_ex1))
# Pmodel_ex1[-1 ]-Pmodel_ex1[-62]
Piid <-1/(1+(S*timesstesp))
PidV <- 1/(1+(S*timesstesp^2))
######################################################################################################
par(mfrow=c(1,2),mar=c(5,5.5,4,2)) 
plot(timesstesp, Piid, col="black",ylim = range(Pmodel_ex1[10:62],PidV [10:62],Piid[10:62] ),type="l",lty=2,xaxt = "n",cex.main=1.8,
     ylab=TeX(r'($P_{\tau}(s_0)$)'), xlab= "Year",lwd=3,main= TeX(r'(a) record probability$)'),cex.lab=1.8,cex.axis=1.5)
abline(v=42,lty=2,col="grey",lwd=3)
lines(timesstesp, PidV, col="blue", lty=2,lwd=3)
lines(timesstesp, Pmodel_ex1, col="red",lty=1, lwd=3.5)
plotaxis.vec <- c(1,seq(10,50,10),62)
axis(1, at= plotaxis.vec, labels= plotaxis.vec + (1960),cex.axis=1.2)
legend('topright',legend=c(TeX(r'($P_{\tau}(s_0)$)'), TeX(r'($P_{iid}$)'), TeX(r'($P_{\alpha}$)')),
       col=c("red",'black',"blue"),lty= c(1,2,2),lwd=5,cex=1.5)

plot(timesstesp,Pmodel_ex1/PidV, col="red",type="l",lty=1,xaxt = "n",cex.main=1.5,
     ylab=TeX(r'($P_{\tau}(s_0)/P_{\alpha}$)'), xlab= "Year",lwd=3,main= "b) detection of trend",cex.lab=1.8,cex.axis=1.5)
abline(h=1,col='blue',lty=2,lwd=3)
abline(v=42,lty=2,col="grey",lwd=3)
lines(timesstesp,Pmodel_ex1/PidV, col="red",type="1",lty=2,lwd=3.5)
axis(1, at= plotaxis.vec, labels= plotaxis.vec + (1960),cex.axis=1.2)
legend('topleft',legend=TeX(r'($P_{\tau}(s_0)/P_{\alpha}$)'), lty=1,lwd=5,cex=1.5, col='red')

######################################################################################################
Simulated_obs_and_Tnext <- DataStructure.func(Tnext=nb.timesteps,
                                              phi.mat = matrix(replicate(S,mu.vec),ncol=S), 
                                              S=S,
                                              xi=xi,
                                              dependence.param=aalpha)
SimulatedData <- Simulated_obs_and_Tnext$SimulatedObservations
obsTnext <- Simulated_obs_and_Tnext$SimulatedNext
PlotData.separate(SimulatedData,siteid = 0, par1 = 2, par2 = 3) 

############################################################################

############################################################################
############################################################################
# Figure S1 - QQ-plot stations
############################################################################
############################################################################

h1 <- 50
h2 <- 50
GmaxHatY_britany <- GmaxY.estimation.function.step1(britany_observations, bandwidth = 50)$Array.GmaxHatY
E.GmaxY_britany <- EGmaxY.estimation.function.step2(GmaxHatY_britany, bandwidth = 50)$Array.EGmaxHatY
s0 <- 3 
Plot.QQplotExp1.ProbMargs_so.VersionSQRTTMAX(GmaxHatY_britany,E.GmaxY_britany ,britany_information$name,site.interest=s0,site.others=0,par1=2,par2=3, add.print=FALSE)


h1 <- 40
h2 <- 40
GmaxHatY_paris<- GmaxY.estimation.function.step1(paris_observations, bandwidth = 40)$Array.GmaxHatY
E.GmaxY_paris <- EGmaxY.estimation.function.step2(GmaxHatY_paris, bandwidth = 40)$Array.EGmaxHatY
s0 <- 1
Plot.QQplotExp1.ProbMargs_so.VersionSQRTTMAX(GmaxHatY_paris,E.GmaxY_paris ,paris_information$name,site.interest=s0,site.others=0,par1=2,par=3, add.print=FALSE)



h1 <- 25
h2 <- 25
GmaxHatY_south <- GmaxY.estimation.function.step1(south_observations, bandwidth = 25)$Array.GmaxHatY
E.GmaxY_south <- EGmaxY.estimation.function.step2(GmaxHatY_south, bandwidth = 25)$Array.EGmaxHatY
s0 <- 1
Plot.QQplotExp1.ProbMargs_so.VersionSQRTTMAX(GmaxHatY_south,E.GmaxY_south,south_information$name,site.interest=s0,site.others=0,par1=2,par2=3, add.print=FALSE)




