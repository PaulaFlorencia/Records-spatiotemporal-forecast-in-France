######################################################################################################
######################################################################################################
# Forecast probability of spatiotemporal records : France 
######################################################################################################
######################################################################################################

# Data #
# Read files : "region"_observations.rds and "region"_information.rds

# Plot : Map #
# stations in France, color code for each region

# Plot : time series #
# yearly maxima of daily maxima temperature in all 15 stations (5 per region) from 1960 to 2022

# Estimation of record probability 
# 

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
###  Plot 1 :  Map 
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
###  Plot 2 - Observation time series
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

library(latex2exp)
G_estimator <- function(t, xt, h) {
  # kernel nonparametric estimator of the CDF G
  function(x, t0) {
    Kvect <- K((t - t0) / h)
    sum(Kvect * (xt <= x)) / sum(Kvect)
  }
}

# Equation (11) sans (\tau-1) dans l'article
V_madogram_func <- function(FiYi.mat,wi.vec,xi.vec){ 
  # Input : 
  #  - FiYi.mat = matrix of values F_t^s(Y_t(s)) with t=row index, s=column index
  #  - wi.vec and xi.vec = vectors of weight values and the x values that defined them.
  FiYiwi <- t(t(FiYi.mat)^(1/wi.vec)) # Lines of FiYi.mat elevated to the power wi.vec
  Mean.FiYiwi<- rowMeans(FiYiwi)  # moyenne spatiale
  Max.FiYiwi<- apply(FiYiwi,1,max) # max spatial
  vw.ti <- mean(Max.FiYiwi - Mean.FiYiwi) # moyenne temporelle des max-min 
  cw.ti <- mean(wi.vec/(1+wi.vec)) 
  Aw.ti <- (vw.ti + cw.ti)/ (1 - vw.ti - cw.ti)
  sum.xi <- sum(xi.vec^(-1)) 
  lambdaV <- sum.xi * Aw.ti
  return(lambdaV)
}

# Kernel
K <- function(x) { # kernel 
  3/4 * (1 - x^2) * (abs(x) <= 1)
}

# meme Kernel
dEpan <- function(x){
  ## Function of Epanechnikov density distribution
  k <- (3/4)* (1-x^2)
  k[-1>x] <- 0
  k[x>1] <- 0
  return (k)
}

# source("/Users/pgonzale/article2_functions.R")


obs_mat <- readRDS("/Users/pgonzale/paris_observations.rds") 
#obs_mat <- readRDS("/Users/pgonzale/south_observations.rds") 
#obs_mat <- readRDS("/Users/pgonzale/britany_observations.rds") 
s0 <- 1 # site d'interet
if (colnames(obs_mat)[1]=="Ile de Groix"){s0 <- 3} # pourquoi cette ligne ?

nbtimesteps <- T <- dim(obs_mat)[1]
nbtausteps <- nbtimesteps-1
timesteps <- seq(1,nbtimesteps,by=1)
nbsites <- dim(obs_mat)[2]
h2 <- rep(25,nbsites) # j?ai consid?r? que c??tait h? . Pas compris l?utilit? de le r?p?ter
h1 <- rep(25,nbsites) # j?ai consid?r? que c??tait h

Kij <- outer(timesteps,timesteps,function(zz,z) dEpan((zz - z) / h1[1]))
W <- Kij / rowSums(Kij) # matrice contenant les w_j(t) pour les estimations des F(Y)
Kij <- outer(timesteps,timesteps,function(zz,z) dEpan((zz - z) / h2[1]))
Wp <- Kij / rowSums(Kij) # matrice contenant les w?_j(t) pour les estimations de E(F(Y))
# Rem : t=num ligne, j=num colonne , w?_j(t)=K((t-j)/h?)/(sum en l des K((t-l)/h?)

# Pr?alable : construction de la liste des matrices contenant les Fchap_{t,s}(Y_j(s_0)) 
#             o? t = num?ro de ligne (entre 1 et T-1) et j=num?ro de colonne
FYlist <- list()
for (s in 1:nbsites){
  Ys <- obs_mat[,s]
  Ys0 <- obs_mat[,1] 
  matI  <- outer(Ys,Ys0,function(y,yy){(y<=yy)+0})
  # FYlist[[s]] <- W %*% matI[-T,] # Julien
  FYlist[[s]] <- (W %*% matI)[-T,] # Soulivanh
}

# Calcul des Echap_{t,s,tau} = estimations des E[F_{t,s}(Y_{\tau}(s_0)]
# On va les stocker dans nbsites matrices index?es par t et tau, 
# puis on r?organise et on en d?duit les uchap_{t,tau}(s_0,s) ds une liste en tau cette fois
matEFYlist <- list()
for (s in 1:nbsites){
  matEFYlist[[s]] <- FYlist[[s]] %*% t(Wp)
  # la tau-?me colonne de cette matrice contient les Echap_{t,s,tau} pour 1 = t = T
  # car la tau-?me colonne de Wp contient les poids ? orbitant ? autour du temps tau
}
mat_utau.list <- list()
for (tau in 2:nbtimesteps){
  mat_utau <- matrix(NA,tau-1,nbsites)
  for (s in 1:nbsites){
    # le vecteur suivant contient les Echapeau_{t,s,tau} pour 1 = t = (tau-1) 
    Echap <- (matEFYlist[[s]])[1:(tau-1),tau] 
    mat_utau[,s] <- Echap/(1-Echap) 
  } 
  mat_utau.list[[tau]] <- mat_utau  # matrice ? tau-1 lignes et S colonnes
} 


# Pr?alable : d?termination de la matrice des F_t,s(Y_t(s)) (t=row,s=col)
FtYt <- matrix(NA,nbtimesteps,nbsites)
for (s in 1:nbsites){
  Ys <- obs_mat[,s]
  matI  <- outer(Ys,Ys,function(y,yy){(y<=yy)+0})
  FtYt_temp <- W %*% matI # matrice des F_t,s(Y_t?(s)) (t=row,t?=col)
  # FtYt <- diag(FtYt_temp) # on ne garde que les coefficients diagonaux # Julien
  FtYt[, s] <- diag(FtYt_temp) # Soulivanh
}

LamV.tau_vec <- rep(NA,(nbtausteps))
LamV.tau_vec_indep <- rep(NA,(nbtausteps))
for (tau in 2:nbtimesteps){
  mat_utau <- mat_utau.list[[tau]]
  sum_Vtau <- 0
  sum_Vtau_indep <- 0
  for (t in 1:(tau-1)){ # ?tape 5.
    x.vec <- mat_utau[t,] # vecteur des u_t,tau(s0,s) (indice = s , t et tau fixes)
    
    #w.vec <- x.vec/sum(x.vec) # w correspondant     ##  !!!
    w.vec <- (1/x.vec)/sum((1/x.vec)) ##  !!! j'ai changé l'originale pour ceci 
    
    # FtYt dans la formule ci-dessous est toujours la m?me, quelles que soient
    # les valeurs des indices tau et t de boucles : il n?y a que les u_t,tau qui changent
    Vtau <- V_madogram_func(FtYt,w.vec,x.vec) # calcul du t-?me terme en V
    sum_Vtau <- sum_Vtau + Vtau 
    
    Vtau_indep <- sum(x.vec^(-1)) 
    sum_Vtau_indep <- sum_Vtau_indep + Vtau_indep 
  }
  LamV.tau_vec[tau-1] <- sum_Vtau 
  LamV.tau_vec_indep[tau-1] <- sum_Vtau_indep 
}
P.s0_vec_indep <- 1/(1+LamV.tau_vec_indep) 
P.s0_vec <- 1/(1+LamV.tau_vec)

# Forecast
#############################################
HW_2toT <- stats::HoltWinters(ts(P.s0_vec, frequency = 1), gamma = FALSE)
HW_2toT.fit  <- HW_2toT$fitted[,1]
HW_Tnext.fit <- predict(HW_2toT, n.ahead = 1, prediction.interval = TRUE)
HW_Pupcoming <- as.numeric(HW_Tnext.fit[,1])
PrecordTnext <- c(P.s0_vec,HW_Pupcoming)

HW_2toT_indep <- stats::HoltWinters(ts(P.s0_vec_indep, frequency = 1), gamma = FALSE)
HW_2toT.fit_indep  <- HW_2toT_indep$fitted[,1]
HW_Tnext.fit_indep <- predict(HW_2toT_indep, n.ahead = 1, prediction.interval = TRUE)
HW_Pupcoming_indep <- as.numeric(HW_Tnext.fit_indep[,1])
PrecordTnext_indep <- c(P.s0_vec_indep,HW_Pupcoming_indep)

#############################################
###  Record probability : plots
#############################################

PrecordTnext_paris <- PrecordTnext

PrecordTnext_paris_indep <- PrecordTnext_indep

end.t <- length(PrecordTnext)

plotaxis.vec <- c(1,seq(10,50,10),length(PrecordTnext))
startYear <- 1960
#############################################
###  Plot 3  - Record probability up to T+1
#############################################

#PrecordTnext_britany <- c(Precord.so.Model_britany,HW_Pupcoming_britany)
#PrecordTnext_paris <- c(Precord.so.Model_paris,HW_Pupcoming_paris)
#PrecordTnext_south <- c(Precord.so.Model_south,HW_Pupcoming_south)

end.t <- dim(britany_observations)[1]
par(mfrow=c(1,1),mar=c(5,5,4,2))
plot(1:length(PrecordTnext_paris),PrecordTnext_paris,log="y", type='l',col='green4',lwd=2, main="Estimated and forecasted record probabilities",ylab=TeX(r'($\widehat{P}_{\tau}(s_o)$)'), xlab="Year",
     ylim = range( 1/(1+(nbsites*c(1:nbtausteps))) , PrecordTnext_paris ),
     xaxt = "n",cex.lab=1.6,cex.axis=1.4, cex.main=1.7)
lines(1:length(PrecordTnext_paris),1/(1+(nbsites*c(1:length(PrecordTnext_paris)))),col='black',lwd=5,log="y")
#lines(1:length(PrecordTnext_paris),PrecordTnext_south,lwd=5,col="red",log="y")
#lines(1:length(PrecordTnext_paris),PrecordTnext_britany,lwd=5,col="blue",log="y")
lines(1:length(PrecordTnext_paris),PrecordTnext_paris ,lwd=5,col="green4",log="y")
abline(v=length(PrecordTnext),lwd=2,col='grey60', lty=2)
#points(end.t,PrecordTnext_britany[length(PrecordTnext_paris)],pch=21,col="black",bg="blue",cex=2,lwd=2)
points(end.t,PrecordTnext_paris[length(PrecordTnext_paris)],pch=21,col="black",bg="green4",cex=2,lwd=2)
#points(end.t,PrecordTnext_south[length(PrecordTnext_paris)],pch=21,col="black",bg="red",cex=2,lwd=2)

axis(1, at= plotaxis.vec, labels= plotaxis.vec + (startYear),cex.axis=1.4)

#legend('top',legend=c("Brest", "Paris-Montsouris" , "Salon de Provence", "i.i.d."),
       col=c("blue",'green4',"red","black"),lty= 1,lwd=5, bg="transparent",cex=1.5)

#############################################
###  Plot 4  - Probability ratio up to T+1
#############################################

par(mfrow=c(1,1),mar=c(5,5.5,4,2))
plot(1:length(PrecordTnext_paris),PrecordTnext_paris/(1/(1+(nbsites*c(1:length(PrecordTnext_paris))))),log="y", type='l',col='green4',lwd=2, main="
Signal on record probability",ylab=TeX(r'($\widehat{P}_{\tau}(s_o)/P_{\tau}^{(i.i.d.)}$)'), xlab="Year",
     ylim = range( 1,10 ),
     xaxt = "n",cex.lab=1.6,cex.axis=1.4,cex.main=1.7)
#lines(1:length(PrecordTnext_paris),PrecordTnext_south/(1/(1+(nbsites*c(1:length(PrecordTnext_paris))))),lwd=5,col="red",log="y")
#lines(1:length(PrecordTnext_paris),PrecordTnext_britany/(1/(1+(nbsites*c(1:length(PrecordTnext_paris))))),lwd=5,col="blue",log="y")
lines(1:length(PrecordTnext_paris),PrecordTnext_paris/(1/(1+(nbsites*c(1:length(PrecordTnext_paris))))) ,lwd=5,col="green4")
abline(v=length(PrecordTnext),lwd=2,col='grey60', lty=2)
#points(end.t,PrecordTnext_britany[length(PrecordTnext_paris)]/(1/(1+(nbsites*end.t))),pch=21,col="black",bg="blue",cex=2,lwd=2)
points(end.t,PrecordTnext_paris[length(PrecordTnext_paris)]/(1/(1+(nbsites*end.t))),pch=21,col="black",bg="green4",cex=2,lwd=2)
#points(end.t,PrecordTnext_south[length(PrecordTnext_paris)]/(1/(1+(nbsites*end.t))),pch=21,col="black",bg="red",cex=2,lwd=2)

axis(1, at= plotaxis.vec, labels= plotaxis.vec + (startYear),cex.axis=1.4)

#legend('topleft',legend=c("Brest", "Paris-Montsouris" , "Salon de Provence"),
       col=c("blue",'green4',"red"),lty= 1,lwd=5, bg="transparent",cex=1.45)

#############################################
###  Plot 5  - Probability ratio up to T+1
#############################################

par(mfrow=c(1,1),mar=c(5,5.5,4,2))
plot(1:length(PrecordTnext),PrecordTnext_paris/PrecordTnext_paris_indep, type='l',col='green4',lwd=5, main="Effect of spatial dependence on record probabilities",ylab=TeX(r'($\widehat{P}_{\tau}(s_o)/\widehat{P}_{\tau}^{(ind.)}$)'), xlab="Year",
     ylim = range(1,3.5),
     xaxt = "n",cex.lab=1.6,cex.axis=1.4,cex.main=1.7)
#lines(1:length(PrecordTnext),1/(1+(nbsites*c(1:length(PrecordTnext)))),col='black',lwd=5)
#lines(1:length(PrecordTnext),PrecordTnext_south/PrecordTnext_south_indep,lwd=5,col="red")
#lines(1:length(PrecordTnext),PrecordTnext_britany/PrecordTnext_britany_indep,lwd=5,col="blue")
#lines(1:length(PrecordTnext),PrecordTnext_paris/PrecordTnext_paris_indep ,lwd=5,col="green4")
abline(v=length(PrecordTnext),lwd=2,col='grey60', lty=2)
#points(end.t,PrecordTnext_britany[length(PrecordTnext_paris)]/PrecordTnext_britany_indep[length(PrecordTnext_paris)],pch=21,col="black",bg="blue",cex=2,lwd=2)
points(length(PrecordTnext),PrecordTnext_paris[length(PrecordTnext)]/PrecordTnext_paris_indep[length(PrecordTnext)],pch=21,col="black",bg="green4",cex=2,lwd=2)
#points(end.t,PrecordTnext_south[length(PrecordTnext_paris)]/PrecordTnext_south_indep[length(PrecordTnext_paris)],pch=21,col="black",bg="red",cex=2,lwd=2)
axis(1, at= plotaxis.vec, labels= plotaxis.vec + (startYear),cex.axis=1.4)
#legend('bottom',legend=c("Brest", "Paris-Montsouris" , "Salon de Provence"),
       col=c("blue",'green4',"red"),lty= 1,lwd=5, bg="transparent",cex=1.5)

#############################################
###  Plot 6  - Beta QQ-plot
#############################################
par(mfrow=c(3,5),mar=c(5,3,4,1))
#par(mfrow=c(2,3),mar=c(5,5,4,2))

utt.vec <- unlist(lapply(mat_utau.list,"[", ,1))
U.vec <- runif(length(utt.vec), min = 0, max = 1)

nf<-layout(matrix(1:15, nrow=3,byrow=TRUE))
layout.show(nf)

titleTex <- c(TeX(r'( Beta QQ-plot $(s_0,s_0)$)'),TeX(r'( Beta QQ-plot $(s_0,s_1)$)'),TeX(r'( Beta QQ-plot $(s_0,s_2)$)'),
              TeX(r'( Beta QQ-plot $(s_0,s_3)$)'),TeX(r'( Beta QQ-plot $(s_0,s_4)$)'))
for(ss in 1:nbsites){
  FYorder.s <- FYlist[[ss]]
  FYorder.s_vec <- FYorder.s[upper.tri(FYorder.s , diag = FALSE)] # 1891
  
  utt.vec <- unlist(lapply(mat_utau.list,"[", ,ss))  
  
  # U.vec <- runif(length(utt.vec), min = 0, max = 1)
  
  Uutt.vec <- U.vec^(1/utt.vec)
  
  if (ss==1){
    #par(mar=c(3,4.7,3.5,1)) # britany
    #par(mar=c(3,4.7,0,1))
    par(mar=c(3,3,3.5,1))
  }
  if (ss>1){
    par(mar=c(3,3,3.5,1))
    #par(mar=c(3,3,0,1))
  }
  
  plot(sort(FYorder.s_vec),sort(Uutt.vec),ylim=c(0,1),xlim=c(0,1), col="green",
       ylab=TeX(r'($\widehat{F}_{t,s}(Y_{\tau}(s_o))$)'), xlab="Empirical Beta quantiles",
       main=titleTex[ss],lwd = 0.5)
  abline(a = 0, b = 1, col = "black", lwd = 2)
}

#############################################
###  Plot 7  - Schema example
#############################################


nbsites <- 5
nbyears <- 2022-1959
timesteps <- c(1:c(nbyears-1))
aalpha <- .3
xi <- -0.2
theta <- 40
mu.0 <- 26
mu.vec<- mu.0 + ((c(timesteps,63)^2)*0.002) + ((c(timesteps,63))*0.001)
mut.mat <- matrix(replicate(nbsites,mu.vec),ncol=nbsites)
sigma.mat <- -xi*(theta -mut.mat) 
sigma.vec <- -xi*(theta -mu.vec) 

Vfunc.vec <- function(u.vec, alph){
  sum(u.vec^(-1/alph))^(alph)
}
s0 <- 1
pRecord <- rep(NA,62)

for (tt in 1: 62){
  cat(" tau - 1: ",tt,"\n")
  sigmaTplusUn_s0 <- sigma.mat[(tt+1),s0]
  sigmat_s <- as.matrix(sigma.mat[1:tt,])
  if (tt==1){sigmat_s <- t(sigmat_s)}
  
  u_tTplusUn_s0s <- (sigmaTplusUn_s0/sigmat_s)^(1/xi)
  sum.u <- 0
  for (ti in 1:dim(u_tTplusUn_s0s)[1]){
    cat(" ti: ",ti,"\n")
    uv <- Vfunc.vec(as.numeric(u_tTplusUn_s0s[ti,]),aalpha)
    #uv <- V.A.func.vec(as.numeric(u_tTplusUn_s0s[ti,]),aalpha)
    cat(" Vu: ",uv,"\n")
    sum.u <- sum.u + uv 
  }
  cat(" SumVu: ",sum.u,"\n")
  pRecord[tt] <- 1/(1+sum.u)
}

Pmodel_ex1 <- pRecord

#############################################

nbyears <- 63
nbsites <- 5
timesteps <- seq(1,nbyears,1)
xi <- -0.2
aalpha <- 0.3
mu0 <- 26
mu.vec <- rep(26,nbyears)
theta <- 40
mut.mat <- matrix(replicate(nbsites,mu.vec),ncol=nbsites)
sigma.mat <- -xi*(theta -mut.mat) 
sigma.vec <- -xi*(theta -mu.vec) 

for (tt in 1: 62){
  cat(" tau - 1: ",tt,"\n")
  sigmaTplusUn_s0 <- sigma.mat[(tt+1),s0]
  sigmat_s <- as.matrix(sigma.mat[1:tt,])
  if (tt==1){sigmat_s <- t(sigmat_s)}
  
  u_tTplusUn_s0s <- (sigmaTplusUn_s0/sigmat_s)^(1/xi)
  sum.u <- 0
  for (ti in 1:dim(u_tTplusUn_s0s)[1]){
    cat(" ti: ",ti,"\n")
    uv <- Vfunc.vec(as.numeric(u_tTplusUn_s0s[ti,]),aalpha)
    #uv <- V.A.func.vec(as.numeric(u_tTplusUn_s0s[ti,]),aalpha)
    cat(" Vu: ",uv,"\n")
    sum.u <- sum.u + uv 
  }
  cat(" SumVu: ",sum.u,"\n")
  pRecord[tt] <- 1/(1+sum.u)
}

PidV <- pRecord

#############################################
nbsites <- 5
nbyears <- 2022-1959
timesteps <- c(1:c(nbyears-1))
aalpha <- .999
xi <- -0.2
theta <- 40
mu.0 <- 26
mu.vec<- rep(mu.0,63)
mut.mat <- matrix(replicate(nbsites,mu.vec),ncol=nbsites)
sigma.mat <- -xi*(theta -mut.mat) 
sigma.vec <- -xi*(theta -mu.vec) 

for (tt in 1: 62){
  cat(" tau - 1: ",tt,"\n")
  sigmaTplusUn_s0 <- sigma.mat[(tt+1),s0]
  sigmat_s <- as.matrix(sigma.mat[1:tt,])
  if (tt==1){sigmat_s <- t(sigmat_s)}
  
  u_tTplusUn_s0s <- (sigmaTplusUn_s0/sigmat_s)^(1/xi)
  sum.u <- 0
  for (ti in 1:dim(u_tTplusUn_s0s)[1]){
    cat(" ti: ",ti,"\n")
    uv <- Vfunc.vec(as.numeric(u_tTplusUn_s0s[ti,]),aalpha)
    #uv <- V.A.func.vec(as.numeric(u_tTplusUn_s0s[ti,]),aalpha)
    cat(" Vu: ",uv,"\n")
    sum.u <- sum.u + uv 
  }
  cat(" SumVu: ",sum.u,"\n")
  pRecord[tt] <- 1/(1+sum.u)
}

Piid <- pRecord

par(mfrow=c(1,2),mar=c(5,5.5,4,2)) 
plot(timesteps, Piid, col="black",ylim = range(Pmodel_ex1[2:62],PidV [2:62],Piid[2:62] ),type="l",lty=2,xaxt = "n",cex.main=1.8,
     ylab=TeX(r'($P_{\tau}(s_0)$)'), xlab= TeX("Forecast year, $\\tau$"),lwd=3,main= TeX(r'(a) record probability$)'),cex.lab=1.8,cex.axis=1.5)
# abline(v=42,lty=2,col="grey",lwd=3)
lines(timesteps, PidV, col="blue", lty=2,lwd=3)
#lines(timesstesp, Pind, col="green", lty=2,lwd=3)
lines(timesteps, Pmodel_ex1, col="red",lty=1, lwd=3.5)
plotaxis.vec <- c(1,seq(10,50,10),62)
axis(1, at= plotaxis.vec, labels= plotaxis.vec + 1,cex.axis=1.2)
legend('topright',
       legend=c(TeX(r'($P_{\tau}(s_0)$)'), 
                TeX(r'($P^{(i.i.d.)}_{\tau}$)'),
                TeX(r'($P^{(alpha)}_{\tau}$)')),
       col=c("red",'black',"blue"),lty= c(1,2,2,2),lwd=5,cex=1.5)

plot(timesteps,Pmodel_ex1/Piid, 
     ylim = range(PidV/Piid, 10), col="red",type="l",lty=1,xaxt = "n",cex.main=1.8,
     ylab=TeX(r'($P_{\tau}(s_0)/P^{(i.i.d.)}_{\tau}$)'), xlab= TeX("Forecast year, $\\tau$"),lwd=3,
     main= TeX(r'(b) deviation from i.i.d.$)'),cex.lab=1.8,cex.axis=1.5)
abline(h=1,col='black',lty=2,lwd=3)
# abline(v=42,lty=2,col="grey",lwd=3)
lines(timesstesp, PidV/Piid, col="blue",lty=2,lwd=3.5)
#lines(timesstesp, Pind/Piid, col="green",lty=2,lwd=3.5)

axis(1, at= plotaxis.vec, labels= plotaxis.vec + 1,cex.axis=1.2)
legend('topleft',
       legend=c(TeX(r'($P_{\tau}(s_0)/P^{(i.i.d.)}_{\tau}$)'),
                TeX(r'($P^{(alpha)}_{\tau}/P^{(i.i.d.)}_{\tau}$)')
       ),
       
       lty=c(1, 2, 2) ,lwd=5,cex=1.5, col=c('red', 'blue')
)

#############################################
###  Multiple simulations 
#############################################

# Theoretical
#############################################
nbsites <- 5
nbyears <- 2022-1959
timesteps <- c(1:c(nbyears-1))
aalpha <- .3
xi <- -0.2
theta <- 40
mu.0 <- 26
mu.vec<- mu.0 + ((c(timesteps,63)^2)*0.0002) + ((c(timesteps,63))*0.001)
mut.mat <- matrix(replicate(nbsites,mu.vec),ncol=nbsites)
sigma.mat <- -xi*(theta -mut.mat) 
sigma.vec <- -xi*(theta -mu.vec) 

Vfunc.vec <- function(u.vec, alph){
  sum(u.vec^(-1/alph))^(alph)
}
s0 <- 1
pRecord <- rep(NA,62)

for (tt in 1: 62){
  cat(" tau - 1: ",tt,"\n")
  sigmaTplusUn_s0 <- sigma.mat[(tt+1),s0]
  sigmat_s <- as.matrix(sigma.mat[1:tt,])
  if (tt==1){sigmat_s <- t(sigmat_s)}
  
  u_tTplusUn_s0s <- (sigmaTplusUn_s0/sigmat_s)^(1/xi)
  sum.u <- 0
  for (ti in 1:dim(u_tTplusUn_s0s)[1]){
    cat(" ti: ",ti,"\n")
    uv <- Vfunc.vec(as.numeric(u_tTplusUn_s0s[ti,]),aalpha)
    #uv <- V.A.func.vec(as.numeric(u_tTplusUn_s0s[ti,]),aalpha)
    cat(" Vu: ",uv,"\n")
    sum.u <- sum.u + uv 
  }
  cat(" SumVu: ",sum.u,"\n")
  pRecord[tt] <- 1/(1+sum.u)
}

Pmodel_ex1 <- pRecord

# Estimator (N times)
#############################################

# function to generate data
DataStructure.func <- function(Tnext = 86, S = 8, dependence.type = "log",
                               dependence.param = 0.3, xi = -0.2, support = 50,
                               phi.mat= matrix(replicate(8,seq(20, 35, length.out = 86)),ncol=8)){
  # objective : simulate max-stable dependent trajectories with same xi and support
  # input: Distribution parameters, dependance stucture, number of sites and time steps
  # output : list of 2 elements: 1) matrix containing trajectories between t=1 and T = Tnext-1
  #                              2) vector of values at Tnext
  
  #############################################################################
  
  # code: 
  
  # Generate stationary maxstable data with defined dependance structure
  dat_maxsable <- evd::rmvevd(Tnext, dep= dependence.param, model = c(dependence.type), d = S, mar=c(0, 1,xi)) # matrix
  
  # empty matrix
  dat_gev <- matrix(NA, nrow=dim(dat_maxsable)[1], ncol=dim(dat_maxsable)[2]) 
  
  # Add trend to our times series
  if (xi < 0){
    sig.mat <- -xi*(support-phi.mat)
    for ( i in 1:dim(dat_maxsable)[2]){
      dat_gev[,i] <- dat_maxsable[,i]*sig.mat[,i] + phi.mat[,i]}
  }
  if (xi > 0){
    mu.mat <- support + phi.mat/xi
    for ( i in 1:dim(dat_maxsable)[2]){
      dat_gev[,i] <- dat_maxsable[,i]*phi.mat[,i] + mu.mat[,i]}
  }
  
  return(list("SimulatedObservations"=dat_gev[1:Tnext-1,],"SimulatedNext"=dat_gev[Tnext,]))
}

N <- 50
N_simulPrecord <- matrix(NA,ncol=N,nrow=62)

for (i in 1:N){
  Simulated_obs_and_Tnext <- DataStructure.func(Tnext=nbyears,phi.mat =mut.mat, S=nbsites,xi=xi,dependence.param=aalpha,support=theta)
  SimulatedData <- Simulated_obs_and_Tnext$SimulatedObservations
  obsTnext <- Simulated_obs_and_Tnext$SimulatedNext
  obs_mat <- rbind(SimulatedData,obsTnext)
  #
  #Simulated_data <- matrix(NA,ncol=5,nrow=63)
  #for(ti in 1: 63){
  # Simulated_data[ti,] <- rmvevd(1, dep = aalpha, model = "log", d = 5, mar = c(mu.vec[ti],sigma.vec[ti],xi))
  #}
  #obs_mat <-Simulated_data
  ###
  
  nbtimesteps <- T <- dim(obs_mat)[1]
  nbtausteps <- nbtimesteps-1
  timesteps <- seq(1,nbtimesteps,by=1)
  nbsites <- dim(obs_mat)[2]
  h2 <- rep(25,nbsites) # j?ai consid?r? que c??tait h? . Pas compris l?utilit? de le r?p?ter
  h1 <- rep(25,nbsites) # j?ai consid?r? que c??tait h
  
  Kij <- outer(timesteps,timesteps,function(zz,z) dEpan((zz - z) / h1[1]))
  W <- Kij / rowSums(Kij) # matrice contenant les w_j(t) pour les estimations des F(Y)
  Kij <- outer(timesteps,timesteps,function(zz,z) dEpan((zz - z) / h2[1]))
  Wp <- Kij / rowSums(Kij) # matrice contenant les w?_j(t) pour les estimations de E(F(Y))
  # Rem : t=num ligne, j=num colonne , w?_j(t)=K((t-j)/h?)/(sum en l des K((t-l)/h?)
  
  # Pr?alable : construction de la liste des matrices contenant les Fchap_{t,s}(Y_j(s_0)) 
  #             o? t = num?ro de ligne (entre 1 et T-1) et j=num?ro de colonne
  FYlist <- list()
  for (s in 1:nbsites){
    Ys <- obs_mat[,s]
    Ys0 <- obs_mat[,1] 
    matI  <- outer(Ys,Ys0,function(y,yy){(y<=yy)+0})
    # FYlist[[s]] <- W %*% matI[-T,] # Julien
    FYlist[[s]] <- (W %*% matI)[-T,] # Soulivanh
  }
  
  # Calcul des Echap_{t,s,tau} = estimations des E[F_{t,s}(Y_{\tau}(s_0)]
  # On va les stocker dans nbsites matrices index?es par t et tau, 
  # puis on r?organise et on en d?duit les uchap_{t,tau}(s_0,s) ds une liste en tau cette fois
  matEFYlist <- list()
  for (s in 1:nbsites){
    matEFYlist[[s]] <- FYlist[[s]] %*% t(Wp)
    # la tau-?me colonne de cette matrice contient les Echap_{t,s,tau} pour 1 = t = T
    # car la tau-?me colonne de Wp contient les poids ? orbitant ? autour du temps tau
  }
  mat_utau.list <- list()
  for (tau in 2:nbtimesteps){
    mat_utau <- matrix(NA,tau-1,nbsites)
    for (s in 1:nbsites){
      # le vecteur suivant contient les Echapeau_{t,s,tau} pour 1 = t = (tau-1) 
      Echap <- (matEFYlist[[s]])[1:(tau-1),tau] 
      mat_utau[,s] <- Echap/(1-Echap) 
    } 
    mat_utau.list[[tau]] <- mat_utau  # matrice ? tau-1 lignes et S colonnes
  } 
  
  
  # Pr?alable : d?termination de la matrice des F_t,s(Y_t(s)) (t=row,s=col)
  FtYt <- matrix(NA,nbtimesteps,nbsites)
  for (s in 1:nbsites){
    Ys <- obs_mat[,s]
    matI  <- outer(Ys,Ys,function(y,yy){(y<=yy)+0})
    FtYt_temp <- W %*% matI # matrice des F_t,s(Y_t?(s)) (t=row,t?=col)
    # FtYt <- diag(FtYt_temp) # on ne garde que les coefficients diagonaux # Julien
    FtYt[, s] <- diag(FtYt_temp) # Soulivanh
  }
  
  LamV.tau_vec <- rep(NA,(nbtausteps))
  for (tau in 2:nbtimesteps){
    mat_utau <- mat_utau.list[[tau]]
    sum_Vtau <- 0
    for (t in 1:(tau-1)){ # ?tape 5.
      x.vec <- mat_utau[t,] # vecteur des u_t,tau(s0,s) (indice = s , t et tau fixes)
      #w.vec <- x.vec/sum(x.vec) # w correspondant                                      ##  !!!!!!! ATTENTION !!!!!!
      w.vec <- (1/x.vec)/sum((1/x.vec)) ##  !!!!!!! ATTENTION !!!!!! j'ai changé l'originale pour ceci 
      # FtYt dans la formule ci-dessous est toujours la m?me, quelles que soient
      # les valeurs des indices tau et t de boucles : il n?y a que les u_t,tau qui changent
      Vtau <- V_madogram_func(FtYt,w.vec,x.vec) # calcul du t-?me terme en V
      sum_Vtau <- sum_Vtau + Vtau 
    }
    LamV.tau_vec[tau-1] <- sum_Vtau 
  }
  P.s0_vec <- 1/(1+LamV.tau_vec) # 
  N_simulPrecord[,i]<-P.s0_vec
  
}

lowQ <- rep(NA,length(P.s0_vec))
upQ <- rep(NA,length(P.s0_vec))
for(tti in 1:length(P.s0_vec)){
  lowQ[tti] <- quantile(N_simulPrecord[tti,], c(0.025, 0.975))[[1]]
  upQ[tti] <- quantile(N_simulPrecord[tti,], c(0.025, 0.975))[[2]]
}


#############################################
###  Plot 8  - Multiple simulations 
#############################################


coll1 <- "red"
coll2 <- "pink1"

#coll1 <- "black"
#coll2 <- "grey"

#coll1 <- "blue"
#coll2 <- "lightblue"


par(mfrow=c(1,2),mar=c(5,5,1.5,2),oma=c(0,0,0,0))

# y-axis
plot(1:length(P.s0_vec),N_simulPrecord[,1], type="l",col=adjustcolor(coll2, alpha = 0.2),
     ylim=c(0,max(N_simulPrecord,P.s0_vec)),lwd=2,main="Record probability",ylab=TeX(r'($\widehat{P}_{\tau}(s_o)$)'),xlab="Year")
for(i in 1:N){
  lines(1:length(P.s0_vec),N_simulPrecord[,i],col = adjustcolor(coll2, alpha = 0.2))
}
lines(1:length(pRecord),lowQ ,lwd=1,col="grey50",lty=2)
lines(1:length(pRecord),upQ ,lwd=1,col="grey50",lty=2)
lines(1:length(pRecord),Pmodel_ex1 ,lwd=1,col="red",lty=2)
lines(1:length(pRecord),Piid ,lwd=2,col="black",lty=2)

lines(1:length(pRecord),pRecord,lwd=3,col=coll1)
#lines(1:length(P.s0_vec),1/(1+(5*c(1:61))),lwd=2,col="black")
legend('topright',legend=c(TeX(r'($P_{\tau}(s_o)$)'),"i.i.d."),
       col=c("red","black"),lty= c(1,2),lwd=5, bg="transparent",cex=1.2)

#legend('topright',legend=TeX(r'($P_{\tau}^{(i.i.d.)}$)'),
  #     col="black",lty= 1,lwd=5, bg="transparent",cex=1.2)
#legend('topright',legend=TeX(r'($P_{\tau}^{(\alpha)}$)'),
   #    col="blue",lty= 1,lwd=5, bg="transparent",cex=1.2)

# Log y-axis
plot(1:length(P.s0_vec),N_simulPrecord[,1], type="l",col=adjustcolor(coll2, alpha = 0.2),ylim= c(min(N_simulPrecord,Pmodel_ex1),max(N_simulPrecord,Pmodel_ex1)),
     lwd=2,main="Log-scale Record probability",ylab=TeX(r'($\widehat{P}_{\tau}(s_o)$)'),xlab="Year", log= "y")
for(i in 1:N){
  lines(1:length(P.s0_vec),N_simulPrecord[,i],col = adjustcolor(coll2, alpha = 0.2))
}
lines(1:length(pRecord),lowQ ,lwd=1,col="grey50",lty=2)
lines(1:length(pRecord),upQ ,lwd=1,col="grey50",lty=2)
lines(1:length(pRecord),Pmodel_ex1 ,lwd=1,col="red",lty=2)
#lines(1:length(pRecord),Piid ,lwd=1,col="grey50",lty=2)
lines(1:length(pRecord),pRecord,lwd=3,col=coll1)
legend('topright',legend=TeX(r'($P_{\tau}(s_o)$)'),
       col="red",lty= 1,lwd=5, bg="transparent",cex=1.2)

