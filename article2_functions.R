############################
# FUNCTIONS ARTICLE 2
###########################

library(knitr)
######################################
# traitement base des donées 
library(maps)
library(ggplot2)
library(dplyr)

######################################
library(cluster)
library(RColorBrewer)
######################################

library(ggpubr)
######################################

library(evd)
library(latex2exp)
library(stats)


K <- function(x) { # kernel 
  3/4 * (1 - x^2) * (abs(x) <= 1)
}

Gmaxtrue_untilt<- function(theta, sigma_t, ksi){ # true Gmax (prod of true G) "glissant" function
  function(x,t) {
    prod(pgev(x, theta, sigma_t[1:t], ksi))
    #pgev(x, theta, sigma_t, ksi)
  }
}

V <- function(sigmaMat,alph){ # logisitic function 
  sumv <- 0
  for (i in 1:(dim(sigmaMat)[2])){
    v <- 1/(sigmaMat[,i]^alph)
    sumv <- sumv + v }
  Vtotal <- sumv^(1/alph)
  return(Vtotal)
}

dEpan <- function(x){
  ## Function of Epanechnikov density distribution
  k <- (3/4)* (1-x^2)
  k[-1>x] <- 0
  k[x>1] <- 0
  return (k)
}

pgev <- function(x, theta, sigma, ksi) {# GEV CDF reparametrized with (xi,theta,sigma)
  exp( -pmax(0, (ksi/sigma * (x - theta)))^(-1/ksi) )
}

G_estimator <- function(t, xt, h) { # kernel nonparametric estimator of the CDF G
  function(x, t0) {
    Kvect <- K((t - t0) / h)
    sum(Kvect * (xt <= x)) / sum(Kvect)
  }
}

Gmax_glissant_estimator <- function(t, Ghat) { # nonparametric estimator of the CDF Gmax until time t-1, where t is our time of interest
  function(x) {  
    # sum()
    tmoins1 <- c(1:(t-1))
    prod(sapply(tmoins1, Ghat, x = x))
    # sapply(t, Ghat, x = x)
  }
}

RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

GenerateStationarySmaple <- function(Ttoday, alpha, S, xi ,mu0, sigma0){
  dat <- evd::rmvevd(Ttoday, dep= alpha, model = c("log"), d = S, mar=c(mu0, sigma0,xi))
  return(dat)
}

AddTrendToData <- function(dat.mat, sigma.mat, support){
  datgev.mat <- matrix(NA, nrow=dim(dat.mat)[1], ncol=dim(dat.mat)[2])
  mu.mat <- support + sigma.mat/xi
  for ( i in 1:dim(dat.mat)[2]){
    datgev.mat[,i] <- dat.mat[,i]*sigma.mat[,i] + mu.mat[,i]
  }
  return(datgev.mat)
}

TruePnext <- function(sigma.mat, sigma.next, xi, alpha){
  sigma_s.mat <- t(t(apply(sigma.mat^(1/xi), 2,sum)))^xi
  lambda_next_s.mat <- (sigma_s.mat/sigma.next)^(1/xi)
  lambda_next_s_inv_sansV <- lambda_next_s.mat^(-1)
  lambda_next_s_inv_sansV <- t(as.numeric(lambda_next_s_inv_sansV))
  lambda_next <- V(lambda_next_s_inv_sansV,alpha)
  TruePrecord <- 1/(1+lambda_next)
  return(TruePrecord)
}

TrueP1toT <- function(sigma.mat, sigma.next, xi, alpha){
  TruePrecord <- matrix(NA, ncol=(dim(sigma.mat)[2]+1), nrow=dim(sigma.mat)[1])
  TrueLambdas <- matrix(NA, ncol=(dim(sigma.mat)[2]+1), nrow=dim(sigma.mat)[1])
  for (i in 1: dim(sigma.mat)[2]){
    sigma_maxst <- cumsum((sigma.mat[,i]^(1/xi)))^xi
    lambda_st <- (sigma_maxst/ c(sigma.mat[-1,1],sigma.next))^(1/xi)
    P_record <- 1/(1+lambda_st)
    TruePrecord[,i] <- P_record
    TrueLambdas[,i] <- lambda_st
  }
  lambdaMatMatinv <- TrueLambdas^(-1)
  lambda_Vt <- V(lambdaMatMatinv[,1:dim(sigma.mat)[2]],alpha)
  Truerecord_all <- 1/(1+lambda_Vt)
  TrueLambdas[,(dim(sigma.mat)[2]+1)] <- lambda_Vt
  TruePrecord[,(dim(sigma.mat)[2]+1)] <- Truerecord_all
  return(list(TruePrecord,TrueLambdas))
}

# Functions for estimated probability
HatGmaxs <- function(dat.mat, theta, sigma.mat, xi , h1, t){
  theta <- theta[1]
  Gmaxhat.mat <- matrix(NA, ncol=dim(sigma.mat)[2], nrow=dim(sigma.mat)[1])
  Gmax.mat <- matrix(NA, ncol=dim(sigma.mat)[2], nrow=dim(sigma.mat)[1])
  for (i in 1: dim(sigma.mat)[2]){
    TrueGmax <- Gmaxtrue_untilt(theta=theta, sigma_t=sigma.mat[,i], ksi=xi)
    Ghats <- G_estimator(t, dat.mat[,i], h1[i])
    for (j in 2: dim(sigma.mat)[1]){
      Gmax.mat[j,i] <- TrueGmax(dat.mat[j,1],(j-1))
      
      Gmaxhats_tj <- Gmax_glissant_estimator(t=j,Ghats)
      Gmaxhat_s <- Gmaxhats_tj(dat.mat[j,1]) 
      Gmaxhat.mat[j,i] <- Gmaxhat_s
    }
    
  }
  Gmaxhat.mat <-as.matrix(Gmaxhat.mat[-1,])
  Gmax.mat <- as.matrix(Gmax.mat[-1,])
  return(list(Gmaxhat.mat,Gmax.mat))
}

HatEGmaxs <- function(h2, t, Gmaxhat.mat, Gmax.mat){
  EGmaxhat.mat <- matrix(NA, ncol=dim(Gmaxhat.mat)[2],nrow=dim(Gmaxhat.mat)[1])
  EGmax.mat <- matrix(NA, ncol=dim(Gmax.mat)[2],nrow=dim(Gmax.mat)[1])
  for (i in 1:dim(Gmaxhat.mat)[2]){
    Kij <- outer(t[-1],t[-1],function(zz,z) dEpan((zz - z) / h2[i])) 
    W <- Kij / rowSums(Kij)
    EGmax.mat[,i] <- W %*% as.numeric(Gmax.mat[,i])
    EGmaxhat.mat[,i] <- W %*% as.numeric(Gmaxhat.mat[,i])
  }
  return(list(EGmaxhat.mat,EGmax.mat))
}

HatPnext <- function(EGmaxhat.mat,EGmax.mat){
  HWpred.mat <- rep(0, dim(EGmax.mat)[2])
  HWhatpred.mat <-rep(0, dim(EGmaxhat.mat)[2])
  Lambdapred.mat <- rep(0, dim(EGmax.mat)[2])
  HatLambdapred.mat <-rep(0, dim(EGmaxhat.mat)[2])
  for (i in 1:dim(EGmax.mat)[2]){
    HWGmax_s <- HoltWinters(ts(EGmax.mat[,i], frequency = 1), gamma = FALSE)
    HWGmaxpred_s <- predict(HWGmax_s , n.ahead = 1, prediction.interval = TRUE)
    HWpred.mat[i] <- HWGmaxpred_s[,1]
    
    HWGmaxhat_s <- HoltWinters(ts(EGmaxhat.mat[,i], frequency = 1), gamma = FALSE)
    HWGmaxhatpred_s <- predict(HWGmaxhat_s, n.ahead = 1, prediction.interval = TRUE)
    HWhatpred.mat[i] <- HWGmaxhatpred_s[,1]
    
    
    Lambdapred.mat[i] <- 1/(HWGmaxpred_s[,1]/(1-HWGmaxpred_s[,1]))
    HatLambdapred.mat[i] <- 1/(HWGmaxhatpred_s[,1]/(1-HWGmaxhatpred_s[,1]))
  }
  return(list(HWhatpred.mat,HWpred.mat,HatLambdapred.mat,Lambdapred.mat ))
}

LambdaSfunc <- function(Lambdapred.mat,alpha){
  Lambdapred.mat <- as.matrix(Lambdapred.mat, ncol=length(Lambdapred.mat) )
  lambdaMatMatinv_from_EGmax_hat <- Lambdapred.mat^(-1)
  lambda_next_from_EGmax_hat <- V(t(lambdaMatMatinv_from_EGmax_hat),alpha)
  Pnext_hat <- 1/(1+lambda_next_from_EGmax_hat)
  return(Pnext_hat)
}


##############################################################################################################################################################################################################################################################

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

TruePnext <- function(sigma.mat, xi, alpha, dependence.type = "log",siteid = 0){
  # Objective : compute true records probability at Tnext, from theoretical parameters
  # input : sigma.mat (matrix of all scale parameters), xi(shape parameter) and type of dependence and it's parameter
  #         sitsid is a value or a vector with the "locations of interest" (where we want to know records probability)
  #         if siteid = 0, all locations are locations of interest
  # output : list of 5 arrays.
  #          for the 3-dimension arrays we have 1d = time (rows) , 2d=sites of comparison (cols), 2d=sites of interest
  
  #############################################################################
  # code: 
  
  Tobs <- (dim(sigma.mat)[1] -1)
  nbsites <- dim(sigma.mat)[2]
  sigma.obs <- sigma.mat[1:Tobs,]
  sigma.next <- sigma.mat[Tobs+1,]
  sigma.obs_s <- t(t(apply(sigma.obs^(1/xi), 2,sum)))^xi 
  
  if (siteid[1] == 0){ # site id = 0 means all sites
    
    #arrays
    interest_sites <- c(1:nbsites)
    P.next <- array(NA,dim = c(1,1,nbsites))
    lambda.next <- array(NA,dim = c(1,1,nbsites))
    PMarginal.next <- array(NA,dim = c(1,nbsites,nbsites))
    lambdaMarginal.next <- array(NA,dim = c(1,nbsites,nbsites))
    
    # reste
    for (i in 1:nbsites){
      lambda.next_s <- (sigma.obs_s/sigma.next[i])^(1/xi)
      lambdaMarginal.next[1,,i] <- lambda.next_s
      PMarginal.next[1,,i] <- 1/(1+lambda.next_s)
      lambda.next_s_inv_sansV <- lambda.next_s^(-1)
      if (dependence.type == "log"){
        lambda.next_s_inv_sansV <- t(as.numeric(lambda.next_s_inv_sansV))
        lambda.nextV <- V(lambda.next_s_inv_sansV,alpha)
        lambda.next[,,i] <- lambda.nextV 
        Precord <- 1/(1+lambda.nextV)
        P.next[,,i] <- Precord
      }
    }
  }
  
  if (siteid[1] != 0){
    sites.so_length <- length(siteid)
    
    #arrays
    interest_sites <- siteid
    P.next <- array(NA,dim = c(1,1,sites.so_length))
    lambda.next <- array(NA,dim = c(1,1,sites.so_length))
    PMarginal.next <- array(NA,dim = c(1,nbsites,sites.so_length))
    lambdaMarginal.next <- array(NA,dim = c(1,nbsites,sites.so_length))
    
    # reste
    for (i in 1:length(siteid)){
      site <- siteid[i]
      lambda.next_s <- (sigma.obs_s/sigma.next[site])^(1/xi)
      lambdaMarginal.next[1,,i] <- lambda.next_s
      PMarginal.next[1,,i] <- 1/(1+lambda.next_s)
      lambda.next_s_inv_sansV <- lambda.next_s^(-1)
      if (dependence.type == "log"){
        lambda.next_s_inv_sansV <- t(as.numeric(lambda.next_s_inv_sansV))
        lambda.nextV <- V(lambda.next_s_inv_sansV,alpha)
        lambda.next[,,i] <- lambda.nextV
        Precord <- 1/(1+lambda.nextV)
      }
      P.next[,,i] <- Precord
    }
  }
  return(list("index_SitesInteret" = interest_sites,
              "ProbabiliteMarginal.next_SitesInteret"=PMarginal.next,
              "lambdaMarginal.next_SitesInteret"=lambdaMarginal.next,
              "ProbabiliteTotal.next_SitesInteret"=P.next,
              "lambdaTotal.next_SitesInteret"=lambda.next))
} 

TrueP1toT <- function(sigma.mat, sigma.next, xi, alpha){
  TruePrecord <- matrix(NA, ncol=(dim(sigma.mat)[2]+1), nrow=dim(sigma.mat)[1])
  TrueLambdas <- matrix(NA, ncol=(dim(sigma.mat)[2]+1), nrow=dim(sigma.mat)[1])
  for (i in 1: dim(sigma.mat)[2]){
    sigma_maxst <- cumsum((sigma.mat[,i]^(1/xi)))^xi
    lambda_st <- (sigma_maxst/ c(sigma.mat[-1,1],sigma_next))^(1/xi)
    P_record <- 1/(1+lambda_st)
    TruePrecord[,i] <- P_record
    TrueLambdas[,i] <- lambda_st
  }
  lambdaMatMatinv <- TrueLambdas^(-1)
  lambda_Vt <- V(lambdaMatMatinv[,1:dim(sigma.mat)[2]],alpha)
  Truerecord_all <- 1/(1+lambda_Vt)
  TrueLambdas[,(dim(sigma.mat)[2]+1)] <- lambda_Vt
  TruePrecord[,(dim(sigma.mat)[2]+1)] <- Truerecord_all
  return(list(TruePrecord,TrueLambdas))
}

TrueP2toTnext <- function(sigma.mat, xi, alpha, dependence.type = "log", siteid = 0){
  # Objective : compute true records probability at from t=2 to t=Ttoday=Tnext-1 from theoretical parameters
  # input : sigma.mat (matrix of all scale parameters), xi(shape parameter) and type of dependence and it's parameter
  #         sitsid is a value or a vector with the "locations of interest" (where we want to know records probability)
  #         if siteid = 0, all locations are locations of interest
  # output : list of 5 arrays.
  #          for the 3-dimension arrays we have 1d = time (rows) , 2d=sites of comparison (cols), 2d=sites of interest
  
  #############################################################################
  # code: 
  
  nbsites <- dim(sigma.mat)[2]
  Tnext <- dim(sigma.mat)[1]
  Tobs <- Tnext - 1
  sigma.obs <- sigma.mat[1:Tobs,]
  sigma.next <- as.numeric(sigma.mat[Tnext,])
  sigmaMax.obs <- (apply(sigma.obs^(1/xi),2,cumsum))^xi # ncol = nbsites , nrow = Tobs
  
  nb.evaluations <- dim(sigma.mat)[1]-1
  tt<-c(1:Tobs)
  
  if (siteid[1] == 0){
    
    # vecteurs a retourner
    interest_sites <- c(1:nbsites)
    # a rendre dans une liste finale
    lambda.2toTnext_marginal <- array(NA,dim = c(Tobs,nbsites,nbsites))
    # a rendre dans une liste finale
    P.2toTnext_marginal <- array(NA,dim = c(Tobs,nbsites,nbsites))
    # a rendre dans une liste finale
    lambda.2toTnext <- matrix(NA,nrow=Tobs, ncol=nbsites)
    # a rendre dans une liste finale
    P.2toTnext <- matrix(NA,nrow=Tobs, ncol=nbsites)
    
    for (i in 1:nbsites){
      lambda.2toTnext_marginal_indx <- (sigmaMax.obs/ sigma.mat[-1, i])^(1/xi) # matrix a sauvgarder
      P.2toTnext_marginal_indx <- 1/(1+lambda.2toTnext_marginal_indx) # matrix a sauvgarder
      
      lambda.2toTnext_marginal[,,i]<- lambda.2toTnext_marginal_indx
      P.2toTnext_marginal[,,i]<- P.2toTnext_marginal_indx 
      
      if (dependence.type == "log"){
        lambdaInv.2toTnext_marginal_indx  <- 1/(lambda.2toTnext_marginal_indx*tt)
        #lambdaInv.2toTnext_marginal_indx <- lambda.2toTnext_marginal_indx^(-1) # <- faire un array de vecterys 
        lambda.2toTnext_indx <- V(lambdaInv.2toTnext_marginal_indx,alpha) # <- mettre dans la liste
        P.2toTnext_indx <- 1/(1+lambda.2toTnext_indx) # un vecteur de la pro de t=2 a t=Tnext <- faire un array de vecterys 
        
        lambda.2toTnext[,i] <- lambda.2toTnext_indx
        P.2toTnext[,i] <- P.2toTnext_indx 
      }
    }
  }
  
  if (siteid[1] != 0){
    sites.so_length <- length(siteid)
    
    interest_sites <- siteid
    # a rendre dans une liste finale
    lambda.2toTnext_marginal <- array(NA,dim = c(p_timesteps,nbsites,sites.so_length))
    # a rendre dans une liste finale
    P.2toTnext_marginal <- array(NA,dim = c(p_timesteps,nbsites,sites.so_length))
    # a rendre dans une liste finale
    lambda.2toTnext <- array(NA,dim = c(p_timesteps,1,sites.so_length))
    # a rendre dans une liste finale
    P.2toTnext <- array(NA,dim = c(p_timesteps,1,sites.so_length))
    
    for (i in 1:sites.so_length){
      site <- siteid[i]
      lambda.2toTnext_marginal_indx <- (sigmaMax.obs/ sigma.mat[-1, site])^(1/xi) # matrix a sauvgarder
      P.2toTnext_marginal_indx <- 1/(1+lambda.2toTnext_marginal_indx) # matrix a sauvgarder
      lambda.2toTnext_marginal[,,i]<- lambda.2toTnext_marginal_indx
      P.2toTnext_marginal[,,i]<- P.2toTnext_marginal_indx 
      if (dependence.type == "log"){
        #lambdaInv.2toTnext_marginal_indx <- lambda.2toTnext_marginal_indx^(-1) # <- faire un array de vecterys 
        lambdaInv.2toTnext_marginal_indx  <- 1/(lambda.2toTnext_marginal_indx*tt)
        lambda.2toTnext_indx <- V(lambdaInv.2toTnext_marginal_indx,alpha) # <- mettre dans la liste
        P.2toTnext_indx <- 1/(1+lambda.2toTnext_indx) # un vecteur de la pro de t=2 a t=Tnext <- faire un array de vecterys 
        lambda.2toTnext[,1,i] <- lambda.2toTnext_indx
        P.2toTnext[,1,i] <- P.2toTnext_indx 
      }
    }
  }
  
  return(list("index_SitesInteret" = interest_sites,
              "lambdaMarginal.2toTnext_SitesInteret"=lambda.2toTnext_marginal,
              "ProbabiliteMarginal.2toTnext_SitesInteret"=P.2toTnext_marginal,
              "lambdaTotal.2toTnext_SitesInteret"=lambda.2toTnext,
              "ProbabiliteTotal.2toTnext_SitesInteret"=P.2toTnext))
} 

Plot.ProbLambda <- function(ArrayProbLambda, siteid=1, site.other = 0, plot.type = "total",
                            par1=1,par2=2,div.screen="TRUE"){
  # Objective : Plot true theoretical marginal and aggregated record probabilities
  # Input : Array coming from function TrueP2toTnext()
  # Output: plots of marginal Prob & lambda or Total Prob & lambda
  # !!! this 
  
  #############################################################################
  # code:
  
  # siteid must be 1 value
  if(div.screen=="TRUE"){par(mfrow=c(par1,par2), mar = c(5, 5.3, 2, 2))}
  
  if(plot.type == "total"){
    
    if(siteid[1]==0){
      siteid <- as.numeric(ArrayProbLambda$index_SitesInteret)
    }
    
    if(length(siteid) == 1){
      proba_total.vec <- as.numeric(ArrayProbLambda$ProbabiliteTotal.2toTnext_SitesInteret[,siteid])
      lambda_total.vec <- as.numeric(ArrayProbLambda$lambdaTotal.2toTnext_SitesInteret[,siteid])
      x_axis <- c(1:length(proba_total.vec))
      plot(x_axis, proba_total.vec, type = "l", col='black',xlab="t from 2 to T+1",ylab= TeX(r'($P_t^{s_o}$)'),
           main=bquote("Record prob at site "~ .(siteid)),lwd=2)
      plot(x_axis, lambda_total.vec, type = "l", col='black',xlab="t from 2 to T+1",ylab= TeX(r'($\lambda_t^{s_o}$)'),
           main=bquote("Associated" ~ lambda ~ "at site "~ .(siteid)),lwd=2)
    }
    
    if(length(siteid) > 1){
      proba_total.mat <- ArrayProbLambda$ProbabiliteTotal.2toTnext_SitesInteret[,siteid]
      lambda_total.mat <- ArrayProbLambda$lambdaTotal.2toTnext_SitesInteret[,siteid]
      
      x_axis <- c(1:dim(proba_total.mat)[1])
      
      for( i in 1:length(siteid)){
        iplot <- siteid[i]
        plot(x_axis, proba_total.mat[,i], type = "l", col='black',xlab="t from 2 to T+1",ylab= TeX(r'($P_t^{s_o}$)'),
             main=bquote("Record prob at site "~ .(iplot)),lwd=2)
        plot(x_axis, lambda_total.mat[,i], type = "l", col='black',xlab="t from 2 to T+1",ylab= TeX(r'($\lambda_t^{s_o}$)'),
             main=bquote("Associated" ~ lambda ~ "at site "~ .(iplot)),lwd=2)
      }
    }
  }
  
  
  if(plot.type == "marginal"){
    site.other <- as.numeric(site.other)
    if(site.other[1]==0){
      site.other <- as.numeric(ArrayProbLambda$index_SitesInteret)
    }
    
    if(siteid[1]==0){
      siteid <- as.numeric(ArrayProbLambda$index_SitesInteret)
    }
    
    if(length(site.other) == 1){
      
      if(length(siteid) == 1){
        proba_marginal.vec <- as.numeric(ArrayProbLambda$ProbabiliteMarginal.2toTnext_SitesInteret[,site.other,siteid])
        lambda_marginal.vec <- as.numeric(ArrayProbLambda$lambdaMarginal.2toTnext_SitesInteret[,site.other,siteid])
        x_axis <- c(1:length(proba_marginal.vec))
        plot(x_axis, proba_marginal.vec, type = "l", col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($P_t^{s,s_o}$)'),
             main=bquote("Marginal record prob at site "~ .(siteid)~" w.r.t site " ~ .(site.other)),lwd=2)
        plot(x_axis, lambda_marginal.vec, type = "l", col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($\lambda_t^{s,s_o}$)'),
             main=bquote("Associated marginal" ~ lambda ~ "at site "~ .(siteid)~" w.r.t site " ~ .(site.other)),lwd=2)
      }
      if(length(siteid) > 1){
        for( i in 1:length(siteid)){
          proba_marginal.mat <- ArrayProbLambda$ProbabiliteMarginal.2toTnext_SitesInteret[,site.other,siteid]
          lambda_marginal.mat <- ArrayProbLambda$lambdaMarginal.2toTnext_SitesInteret[,site.other,siteid]
          x_axis <- c(1:dim(proba_marginal.mat)[1])
          iplot <- siteid[i]
          plot(x_axis, proba_marginal.mat[,i], type = "l", col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($P_t^{s,s_o}$)'),
               main=bquote("Marginal record prob at site "~ .(iplot)~" w.r.t site " ~ .(site.other)),lwd=2)
          plot(x_axis, lambda_marginal.mat[,i], type = "l", col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($\lambda_t^{s,s_o}$)'),
               main=bquote("Associated marginal " ~ lambda ~ "at site "~ .(iplot)~" w.r.t site " ~ .(site.other)),lwd=2)
        }
      }
    }
    
    
    if(length(site.other) > 1){
      
      if(div.screen=="TRUE"){par(mfrow=c(par1,par2), mar = c(5, 5.3, 4, 2))}
      
      if(length(siteid) == 1){
        
        proba_marginal.mat <- ArrayProbLambda$ProbabiliteMarginal.2toTnext_SitesInteret[,site.other,siteid]
        lambda_marginal.mat <- ArrayProbLambda$lambdaMarginal.2toTnext_SitesInteret[,site.other,siteid]
        x_axis <- c(1:dim(ArrayProbLambda$ProbabiliteMarginal.2toTnext_SitesInteret)[1])
        for( j in 1:length(site.other)){
          jplot <- site.other[j]
          plot(x_axis,proba_marginal.mat[,j],type='l', col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($P_t^{s,s_o}$)'),
               main=bquote(atop("Marginal Record prob at site " ~ .(siteid), " w.r.t site " ~ .(jplot))))
          plot(x_axis,lambda_marginal.mat[,j],type='l', col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($\lambda_t^{s,s_o}$)'),
               main=bquote(atop("Associated marginal " ~ lambda ~ "at site "~ .(siteid)," w.r.t site " ~ .(jplot))))
        }
      }
      
      if(length(siteid) > 1){
        proba_marginal.array <- ArrayProbLambda$ProbabiliteMarginal.2toTnext_SitesInteret[,site.other,siteid]
        lambda_marginal.array <- ArrayProbLambda$lambdaMarginal.2toTnext_SitesInteret[,site.other,siteid]
        x_axis <- c(1:dim(ArrayProbLambda$ProbabiliteMarginal.2toTnext_SitesInteret)[1])
        for (i in 1:length(siteid)){
          iplot <- siteid[i]
          for (j in 1:length(site.other)){
            jplot <- site.other[j]
            plot(x_axis,proba_marginal.array[,j],type='l', col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($P_t^{s,s_o}$)'),
                 main=bquote(atop("Marginal Record prob at site "~ .(iplot)," w.r.t site " ~ .(jplot))))
            plot(x_axis,proba_marginal.array[,j],type='l', col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($P_t^{s,s_o}$)'),
                 main=bquote(atop("Associated marginal "~ lambda ~ "at site "~ .(iplot)," w.r.t site " ~ .(jplot))))
          }
        }
      }
    }
  }
  
  
  
}

PlotData.separate <- function(Data, siteid = 0, par1 = 4, par2 = 2 , add.max=FALSE , max.type="all" ){ # if  siteid = 0 we plot data in all sites 
  par(mfrow=c(par1,par2))
  timesteps <- dim(Data)[1]
  nbsites <- dim(Data)[2]
  if (length(siteid)==1){
    if (siteid>0){
      par1 <- 1
      par2 <- 1 
      plot(c(1:timesteps),Data[,siteid],ylim=c(min(Data),max(Data)),
           xlab="time",ylab='observations', main=bquote("Observations at site "~ .(siteid)))
      if(add.max==TRUE){
        max.traj <- apply(X=Data, MARGIN=1, FUN=max)
        points(c(1:timesteps),max.traj,col='red',pch=4)
        legend('bottomright',legend="max all sites",col="red",pch=4, bg="transparent")
      }
      #text(c(1:timesteps)[10],max(Data[1:(timesteps/3)])-1, cex=1.3,col="blue",paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
    }
    
    
    
    
    
    if (siteid==0){
      ylim.max = max(Data)
      ylim.min = min(Data)
      if(add.max==FALSE){
        for(i in 1:nbsites){
          plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="time",ylab='observations',
               main=bquote("Observations at site "~ .(i)))
          #text(c(1:timesteps)[10],max(Data[1:(timesteps/3)])-1,cex=1.3,col="blue", paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
        }
      }
      
      if(add.max==TRUE){
        max.traj <- apply(X=Data, MARGIN=1, FUN=max)
        for(i in 1:nbsites){
          plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="time",ylab='observations',
               main=bquote("Observations at site "~ .(i)))
          points(c(1:timesteps),max.traj,col='red',pch=4)
          #text(c(1:timesteps)[10],max(Data[1:(timesteps/3)])-1,cex=1.3,col="blue", paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
          legend('bottomright',legend="max all sites",col="red",pch=4, bg="transparent")
        }
      }
    }
  }
  
  if (length(siteid)>1){
    
    nb.plots <- length(siteid)
    ylim.max = max(Data)
    ylim.min = min(Data)
    
    if(add.max==FALSE){
      for(i in siteid ){
        plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="time",ylab='observations',
             main=bquote("Observations at site "~ .(i)))
        #text(c(1:timesteps)[10],max(Data[1:(timesteps/3)])-1, cex=1.3,col='blue',paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
      }
    }
    
    if(add.max==TRUE){
      if(max.type=="all"){
        max.traj <- apply(X=Data, MARGIN=1, FUN=max)
        colmax <- "red"
        legend.text <- "max all sites" }
      
      if(max.type=="not.all"){
        max.traj <- apply(X=Data[,c(siteid)], MARGIN=1, FUN=max)
        colmax <- "green"
        legend.text <- "max showed sites"}
      
      for(i in siteid ){
        plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="time",ylab='observations',
             main=bquote("Observations at site "~ .(i)))
        points(c(1:timesteps),max.traj,col=colmax,pch=4)
        #text(c(1:timesteps)[10],max(Data[1:(timesteps/3)])-1, cex=1.3,col='blue',paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
        legend('bottomright',legend=legend.text,col=colmax,pch=4, bg="transparent")
      }
    }
    
  }
}

PlotData.compare <- function(Data, site.interest = 1, site.others = 0,  max.plot = FALSE){
  
  if(length(site.interest) != 1){
    stop("We must observe only 1 site of interest")}
  
  nbtimesteps <- dim(Data)[1]
  timesteps <- c(1:nbtimesteps)
  nbsites <- dim(Data)[2]
  
  site.id <- site.interest
  others.id <- site.others
  col.points.max <- "green"
  if(length(site.others) == 1){
    if(site.others == 0){
      others.id <- c(1:nbsites)
      col.points.max <- "red"}
  }
  
  if(max.plot == FALSE){par(mfrow=c(1,1), mar = c(5, 5.3, 4, 2))}
  if(max.plot == TRUE){par(mfrow=c(1,2), mar = c(5, 5.3, 4, 2))}
  
  ylim.max = max(Data)
  ylim.min = min(Data)
  
  sitesother.names <- paste(others.id,collapse=" ")
  
  plot(timesteps,Data[,site.id],ylim=c(ylim.min,ylim.max),xlab="time",ylab='observations',
       main=bquote(atop("Observations at site of interest "~ .(site.id),"w.r.t to other sites " ~ .(sitesother.names))))
  for(i in others.id ){points(timesteps,Data[,i],col='gray')}
  points(timesteps,Data[,site.id],col="black")
  # text(c(1:nbtimesteps)[13],max(Data)-2, cex=1.3,col='blue',paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
  legend('bottomright',legend=c("site interest","other sites"),col=c("black",'gray'),pch=1, bg="transparent")
  
  if(max.plot == TRUE){
    plot(timesteps,Data[,site.id],ylim=c(ylim.min,ylim.max),xlab="time",ylab='observations',
         main=bquote(atop("Observations at site of interest "~ .(site.id),"w.r.t to other sites " ~ .(sitesother.names))))
    for(i in others.id ){points(timesteps,Data[,i],col='gray')}
    points(timesteps,Data[,site.id],col="black")
    max.traj <- apply(X=Data[,c(site.id,others.id)], MARGIN=1, FUN=max)
    points(timesteps,max.traj,pch=4,col=col.points.max)
    # text(c(1:nbtimesteps)[13],max(Data)-2, cex=1.3,col='blue',paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
    legend('bottomright',legend=c("site interest","other sites","max plotted sites"),col=c("black",'gray',col.points.max),pch=c(1,1,4), bg="transparent")
  }
}

GmaxY.estimation.function.step1 <- function(Data, site.interest = 0, site.others = 0, bandwidth= 35, plot=FALSE,
                                            plot.site.interest=1,plot.site.others=1){
  nbtimesteps <- dim(Data)[1]
  timesteps <- seq(1,nbtimesteps,by=1)
  nbsites <- dim(Data)[2]
  
  if(length(site.interest) > nbsites | length(site.others) > nbsites ){
    stop("the number of sites of study exceeds the number of columns of data ")}
  
  if(site.interest[1] == 0){site.interest <- seq(1,dim(Data)[2],by=1)}
  
  if(site.others[1] == 0){site.others <- seq(1,dim(Data)[2],by=1)}
  
  nbsites.interest <- length(site.interest)
  nbsites.others <- length(site.others)
  
  if (length(bandwidth) == 1){bandwidth <- rep(bandwidth,nbsites.others)}
  
  if(length(bandwidth)!= nbsites.others){
    stop("Bandwight must be a value or a vector of length nbsites.others")}
  
  array.Gmaxhat.mat <- array(NA,c((dim(Data)[1]-1),nbsites.others,nbsites.interest))
  
  for (i.dim3 in 1:nbsites.interest){
    i.interest <- site.interest[i.dim3]
    
    Gmaxhat.mat <- matrix(NA,ncol=nbsites.others,nrow=dim(Data)[1])
    
    for (i.dim2 in 1:nbsites.others){
      i.others <- site.others[i.dim2]
      Ghat <- G_estimator(timesteps, Data[,i.others], bandwidth[i.dim2])
      
      for (i.dim1 in 2:nbtimesteps){
        Gmaxhat_ti.func <- Gmax_glissant_estimator(t = i.dim1,Ghat)
        Gmaxhat <- Gmaxhat_ti.func(Data[i.dim1,i.interest]) 
        Gmaxhat.mat[i.dim1,i.dim2] <- Gmaxhat
      }
    }
    Gmaxhat.mat <-as.matrix(Gmaxhat.mat[-1,])
    array.Gmaxhat.mat[,,i.dim3] <- Gmaxhat.mat
  }
  
  id.1.plot <- which(as.numeric(site.interest) == plot.site.interest)
  id.2.plot <- which(as.numeric(site.others) == plot.site.others)
  
  if (plot==TRUE){
    par(mfrow=c(1,1), mar = c(5, 5.3, 4, 2))
    plot(timesteps[-1],array.Gmaxhat.mat[,id.2.plot,id.1.plot],xlab='t from 2 to T',
         ylab=TeX(r'($\hat{G}_{\max,t -1}^{s}(Y_{t}^{s_o})$)'),
         main=bquote("GmaxHat_s(Y_so) with so =  "~ .(id.1.plot) ~ "and s = " ~ .(id.2.plot)),col='red')
  }
  
  return(list("Array.GmaxHatY" = array.Gmaxhat.mat))
}

plot.separate.function.step.1 <- function(Array.step1, third.dimension = 1, second.dimension=0, par1=1, par2=1){
  par(mfrow=c(par1,par2), mar = c(5, 5.3, 4, 2))
  timesteps <- dim(Array.step1)[1]
  nbsites.interest <- dim(Array.step1)[3]
  nbsites.others <- dim(Array.step1)[2]
  
  if(second.dimension[1] == 0){second.dimension <- seq(1,nbsites.others,by=1)}
  
  if(length(second.dimension)==1){
    plot(c(1:timesteps),Array.step1[,second.dimension,third.dimension],ylim=c(min(Array.step1[,,third.dimension]),max(Array.step1[,,third.dimension])),
         xlab="t from 2 to T",ylab=TeX(r'($\hat{G}_{\max,t -1}^{s}(Y_{t}^{s_o})$)'), 
         main=bquote("GmaxHatY at 3d :"~ .(third.dimension) ~ " and 2d :"~ .(second.dimension)),col='red')
  }
  
  if(length(second.dimension)>1){
    for(i in 1:length(second.dimension)){
      plot(c(1:timesteps),Array.step1[,second.dimension[i],third.dimension],ylim=c(min(Array.step1[,,third.dimension]),max(Array.step1[,,third.dimension])),
           xlab="t from 2 to T",ylab=TeX(r'($\hat{G}_{\max,t -1}^{s}(Y_{t}^{s_o})$)'), 
           main=bquote("GmaxHatY at 3d :"~ .(third.dimension) ~ " and 2d :"~ .(second.dimension[i])),col="red")
    }
  }
}

EGmaxY.estimation.function.step2 <- function(Array.GmaxY, bandwidth= 45){
  nbtimesteps <- dim(Array.GmaxY)[1]
  nbsites.others <- dim(Array.GmaxY)[2]
  nbsites.interest <- dim(Array.GmaxY)[3]
  
  tt <- c(1:nbtimesteps)
  # timesteps <- seq(1,nbtimesteps,by=1)
  
  if (length(bandwidth) == 1){bandwidth <- rep(bandwidth,nbsites.interest)}
  
  if(length(bandwidth)!= nbsites.interest){
    stop("Bandwight must be a value or a vector of length nbsites.interest, i.e dim(Array.GmaxY)[3] ")}
  
  array.EGmaxHatY.mat <- array(NA,c(nbtimesteps,nbsites.others,nbsites.interest))
  
  for(i.interest in 1:nbsites.interest){
    #Gmaxhat_i.mat <- as.matrix(Array.GmaxY[,,i.interest])
    
    EGmaxhat_i.mat <- matrix(NA, ncol=nbsites.others,nrow=nbtimesteps)
    
    for (i.others in 1:nbsites.others){
      Kij <- outer(tt,tt,function(zz,z) dEpan((zz - z) / bandwidth[i.interest]))
      W <- Kij / rowSums(Kij)
      EGmaxhat_i.mat[,i.others] <- W %*% as.numeric(Array.GmaxY[,i.others,i.interest])
    }
    
    array.EGmaxHatY.mat[,,i.interest] <- EGmaxhat_i.mat
    
  }
  return(list("Array.EGmaxHatY" = array.EGmaxHatY.mat))
}

Kernel.marginal.estimation.step1step2 <- function(Data, site.interest = 0, site.others = 0, bandwidth.step1 = 35,bandwidth.step2 = 45){
  Array.step1 <- GmaxY.estimation.function.step1(Data = Data, site.interest = site.interest, site.others = site.others, bandwidth= bandwidth.step1)$Array.GmaxHatY
  Array.E.kernel.estimation <- EGmaxY.estimation.function.step2(Array.GmaxY= Array.step1, bandwidth= bandwidth.step2)$Array.EGmaxHatY
  Array.lambda.kernel.estimation <- (1-Array.E.kernel.estimation)/Array.E.kernel.estimation
  return(list("Marginal.probability.estimation" = Array.E.kernel.estimation,"Marginal.lambda.estimation" = Array.lambda.kernel.estimation))
} 

Marginal.aggregation.function.step3 <- function(Marginal.Array, dependence.type="log", dep.param = 0.3){
  
  marginal.lambdaHat <- (1-Marginal.Array)/Marginal.Array # Array
  
  nbtimesteps <- dim(marginal.lambdaHat)[1]
  nbsites.others <- dim(marginal.lambdaHat)[2]
  nbsites.interest <- dim(marginal.lambdaHat)[3]
  
  mat.aggregated.probability <- matrix(NA,ncol=nbsites.interest,nrow=nbtimesteps)
  mat.aggregated.lambda <- matrix(NA,ncol=nbsites.interest,nrow=nbtimesteps)
  
  if (dependence.type=="log"){
    for (i in 1:nbsites.interest){
      marginal.lambda.site.i <- as.matrix(marginal.lambdaHat[,,i])
      marginal.lambda.site.i.invSansV <- marginal.lambda.site.i^(-1)
      aggregated.lambda.site.i <- V(marginal.lambda.site.i.invSansV,dep.param)
      mat.aggregated.lambda[,i] <- aggregated.lambda.site.i
      aggregated.probability.site.i <- 1/(1+aggregated.lambda.site.i)
      mat.aggregated.probability[,i] <- aggregated.probability.site.i
    }
  }
  
  return(list("Kernel.aggregated.probability" = mat.aggregated.probability,"Kernel.aggregated.lambda" = mat.aggregated.lambda))
  
}

Estimation.RecordProbability <- function(Data, site.interest = 0, site.others = 0, 
                                         bandwidth.step1 = 35,bandwidth.step2 = 45,
                                         dependence.type="log", dep.param = 0.3){
  
  Array.step1 <- GmaxY.estimation.function.step1(Data = Data, site.interest = site.interest, site.others = site.others, bandwidth= bandwidth.step1)$Array.GmaxHatY
  Array.E.kernel.estimation <- EGmaxY.estimation.function.step2(Array.GmaxY= Array.step1, bandwidth= bandwidth.step2)$Array.EGmaxHatY
  
  Array.marginal.lambdaHat <- (1-Array.E.kernel.estimation)/Array.E.kernel.estimation 
  
  Aggregated.Kernel.estimation <- Marginal.aggregation.function.step3(Marginal.Array = Array.E.kernel.estimation, dependence.typ = dependence.type, dep.param = dep.param)
  
  Aggregated.probability <- Aggregated.Kernel.estimation$Kernel.aggregated.probability
  Aggregated.lambda <- Aggregated.Kernel.estimation$Kernel.aggregated.lambda
  
  return(list("Aggregated.probability.estimation" = Aggregated.probability,
              "Aggregated.lambda.estimation" = Aggregated.lambda,
              "Marginal.probability.estimation" = Array.E.kernel.estimation,
              "Marginal.lambda.estimation" = Array.marginal.lambdaHat))
} 

Compare.estimation.and.forecast.plot.func <- function(List.Forecast.object, List.Kernel.Estimation, Listes.True.quantities, type.plot="total", siteid = 0 , site.others = 1 , par1 = 3,par2 = 3 ){
  
  
  indexes.sites <- Listes.True.quantities$index_SitesInteret
  
  if(type.plot=="total"){
    
    par(mfrow=c(par1,par2), mar = c(5, 5.3, 2, 2))
    
    # lambda <- Listes.True.quantities$lambdaTotal.2toTnext_SitesInteret
    probability <- Listes.True.quantities$ProbabiliteTotal.2toTnext_SitesInteret
    nbtimesteps <- dim(probability)[1]
    HW.probability.evolution <- List.Forecast.object$HW.fit
    HW.probability.Tnext <- List.Forecast.object$HW.forecast
    probability.kernel <- List.Kernel.Estimation$Aggregated.probability.estimation
    # lambda.kernel <- List.Kernel.Estimation$Aggregated.lambda.estimation
    if(is.matrix(HW.probability.evolution) == FALSE){
      stop("Please use the aggregated forecast object as imput")
    }
    
    if(siteid[1]==0){
      siteid <- as.numeric(indexes.sites)
    }
    
    ylim.max <- max(probability,HW.probability.evolution)
    
    
    for (i in 1:length(siteid)){
      site.interest <- siteid[i]
      plot(c(1:nbtimesteps), probability[,site.interest], col='black', lwd=2, ylab=TeX(r'($P_{t}^{s_o}$)'), xlab="t from 2 to T+1",
           main=bquote("Record prob at site "~ .(site.interest)),ylim=c(0,ylim.max),xlim=c(1,nbtimesteps),type="l")
      points(c(4:nbtimesteps),HW.probability.evolution[,site.interest],col="green",pch=19)
      points(c(2:nbtimesteps),probability.kernel[,site.interest],pch=4,col="red")
      points(nbtimesteps,HW.probability.Tnext[site.interest],col="blue")
    }
    
  }
  
  
  if(type.plot=="marginal"){
    
    par(mfrow=c(par1,par2), mar = c(5, 5.3, 4, 2))
    
    # lambda <- Listes.True.quantities$lambdaMarginal.2toTnext_SitesInteret
    probability <- Listes.True.quantities$ProbabiliteMarginal.2toTnext_SitesInteret
    nbtimesteps <- dim(probability)[1]
    HW.probability.evolution <- List.Forecast.object$HW.fit
    HW.probability.Tnext <- List.Forecast.object$HW.forecast
    probability.kernel <- List.Kernel.Estimation$Marginal.probability.estimation
    # lambda.kernel <- List.Kernel.Estimation$Marginal.lambda.estimation
    
    if(is.array(HW.probability.evolution) == FALSE){
      stop("Please use the marginal forecast object as imput")
    }
    
    if(siteid[1]==0){
      siteid <- as.numeric(indexes.sites)
    }
    
    if(site.others[1]==0){
      site.others <- as.numeric(indexes.sites)
    }
    
    ylim.max <- max(probability,HW.probability.evolution)
    
    for(i in 1:length(siteid)){
      for (j in 1:length(site.others)){
        site.interest <- siteid[i]
        jplot <- site.others[j]
        plot(c(1:nbtimesteps), probability[,jplot,site.interest], col='blue', lwd=2, ylab=TeX(r'($P_{t}^{s,s_o}$)'), xlab="t from 2 to T+1",
             main=bquote(atop("Marginal Record prob at site " ~ .(site.interest), " w.r.t site " ~ .(jplot))),ylim=c(0,ylim.max),xlim=c(1,nbtimesteps),type="l")
        points(c(4:nbtimesteps),HW.probability.evolution[,jplot,site.interest],col="green",pch=19)
        points(c(2:nbtimesteps),probability.kernel[,jplot,site.interest],pch=4,col="red")
        points(nbtimesteps,HW.probability.Tnext[jplot,site.interest],col="blue")
      }
    }
    
  }
  
}


Compare.estimation.plot.func <- function(List.Kernel.Estimation, Listes.True.quantities, type.plot="total",
                                         plot.lambda = FALSE, siteid = 0, site.other = 0, par1 = 4, par2 = 2 ){
  
  if(siteid[1]==0){
    siteid <- as.numeric(Listes.True.quantities$index_SitesInteret)
  }
  
  if(type.plot=="total"){
    par(mfrow=c(par1,par2), mar = c(5, 5.3, 2, 2))
    
    Estimated.probability <- List.Kernel.Estimation$Aggregated.probability.estimation
    Estimated.lambda <- List.Kernel.Estimation$Aggregated.lambda.estimation
    True.probability <- Listes.True.quantities$ProbabiliteTotal.2toTnext_SitesInteret[-nrow(Listes.True.quantities$ProbabiliteTotal.2toTnext_SitesInteret), ]
    True.lambda <- Listes.True.quantities$lambdaTotal.2toTnext_SitesInteret[-nrow(Listes.True.quantities$ProbabiliteTotal.2toTnext_SitesInteret), ]
    nbtimesteps <- dim(True.lambda)[1]
    ylim.max.p <- max(True.probability,Estimated.probability)
    
    if(length(siteid)==1){
      plot(c(1:nbtimesteps),True.probability[,siteid],type="l", col="black",lwd=2,
           main=bquote("Record prob at site "~ .(siteid)),ylim=c(0,ylim.max.p),
           ylab=TeX(r'($P_t^{s_o}$)'),xlab="t from 2 to T")
      points(c(1:nbtimesteps),Estimated.probability[,siteid],col='red')
      if(plot.lambda == TRUE){
        ylim.max.l <- max(True.lambda, Estimated.lambda)
        plot(c(1:nbtimesteps),True.lambda[,siteid],type="l", col="black",lwd=2,
             main=bquote("Record prob at site "~ .(siteid)),ylim=c(0,ylim.max.l),
             ylab=TeX(r'($\lambda_t^{s_o}$)'),xlab="t from 2 to T")
        points(c(1:nbtimesteps),Estimated.lambda[,siteid],col='red')
      }
    }
    
    if(length(siteid)>1){
      for (i in 1:length(siteid)){
        siteid.i <- siteid[i]
        plot(c(1:nbtimesteps),True.probability[,siteid.i],type="l", col="black",lwd=2,
             main=bquote("Record prob at site "~ .(siteid.i)),ylim=c(0,ylim.max.p),
             ylab=TeX(r'($P_t^{s_o}$)'),xlab="t from 2 to T")
        points(c(1:nbtimesteps),Estimated.probability[,siteid.i],col='red')
        if(plot.lambda == TRUE){
          ylim.max.l <- max(True.lambda, Estimated.lambda)
          plot(c(1:nbtimesteps),True.lambda[,siteid.i],type="l", col="black",lwd=2,
               main=bquote("Associtaed lambda at site "~ .(siteid.i)), ylim=c(0,ylim.max.l),
               ylab=TeX(r'($\lambda_t^{s_o}$)'),xlab="t from 2 to T")
          points(c(1:nbtimesteps),Estimated.lambda[,siteid.i],col='red')
        }
      }
    }
  }
  
  
  if(type.plot=="marginal"){
    par(mfrow=c(par1,par2), mar = c(5, 5.3, 4, 2))
    if(site.other[1]==0){
      site.other <- as.numeric(Listes.True.quantities$index_SitesInteret)
    }
    Estimated.probability <- List.Kernel.Estimation$Marginal.probability.estimation
    Estimated.lambda <- List.Kernel.Estimation$Marginal.lambda.estimation
    True.probability <- Listes.True.quantities$ProbabiliteMarginal.2toTnext_SitesInteret[-(dim(Listes.True.quantities$ProbabiliteMarginal.2toTnext_SitesInteret)[1]),,]
    True.lambda <- Listes.True.quantities$lambdaMarginal.2toTnext_SitesInteret[-(dim(Listes.True.quantities$lambdaMarginal.2toTnext_SitesInteret)[1]),,]
    nbtimesteps <- dim(True.lambda)[1]
    ylim.max.p <- max(True.probability,Estimated.probability)
    # a[-(dim(a)[1]),,]
    if(length(siteid)==1){
      if(length(site.other)==1){
        plot(c(1:nbtimesteps),True.probability[,site.other,siteid],type="l", col="blue",lwd=2,
             main=bquote(atop("Marginal record prob at site "~ .(siteid)," w.r.t site " ~ .(site.other))),
             ylim=c(0,ylim.max.p), ylab=TeX(r'($P_t^{s,s_o}$)'),xlab="t from 2 to T")
        points(c(1:nbtimesteps),Estimated.probability[,site.other,siteid],col='red')
        if(plot.lambda == TRUE){
          ylim.max.l <- max(True.lambda, Estimated.lambda)
          plot(c(1:nbtimesteps),True.lambda[,site.other,siteid],type="l", col="blue",lwd=2,ylim=c(0,ylim.max.l),
               main=bquote(atop("Associated marginal "~ lambda ~ "at site "~ .(siteid)," w.r.t site " ~ .(site.other))),
               ylab=TeX(r'($\lambda_t^{s,s_o}$)'),xlab="t from 2 to T")
          points(c(1:nbtimesteps),Estimated.lambda[,site.other,siteid],col='red')
        }
      }
      if(length(site.other)>1){
        for (j in 1:length(site.other)){
          other.j <- site.other[j]
          plot(c(1:nbtimesteps),True.probability[,other.j,siteid],type="l", col="blue",lwd=2,
               main=bquote(atop("Marginal record prob at site "~ .(siteid)," w.r.t site " ~ .(other.j))),
               ylim=c(0,ylim.max.p), ylab=TeX(r'($P_t^{s,s_o}$)'),xlab="t from 2 to T")
          points(c(1:nbtimesteps),Estimated.probability[,other.j,siteid],col='red')
          if(plot.lambda == TRUE){
            ylim.max.l <- max(True.lambda, Estimated.lambda)
            plot(c(1:nbtimesteps),True.lambda[,other.j,siteid],type="l", col="blue",lwd=2,ylim=c(0,ylim.max.l),
                 main=bquote(atop("Associated marginal "~ lambda ~ "at site "~ .(siteid)," w.r.t site " ~ .(other.j))),
                 ylab=TeX(r'($\lambda_t^{s,s_o}$)'),xlab="t from 2 to T")
            points(c(1:nbtimesteps),Estimated.lambda[,other.j,siteid],col='red')
          }
        }
      }
    }
    
    if(length(siteid)>1){stop("please chose a site of interest")}
  }
}

Forecast.Holt.Winters.func <- function(input.traj){
  
  if(is.matrix(input.traj) == TRUE){ # on travaille avec les probas/lambdas aggregés.
    nbtimesteps <- dim(input.traj)[1]
    nbsites.interest <- dim(input.traj)[2]
    HW_2toT.fit.all <- matrix(NA,nrow=(nbtimesteps-2),ncol=nbsites.interest)
    HW_Tnext.fit.all <- rep(NA,nbsites.interest)
    for (i in 1:nbsites.interest){
      HW_2toT <- HoltWinters(ts(input.traj[,i], frequency = 1), gamma = FALSE)
      HW_2toT.fit  <- HW_2toT$fitted[,1]# green lines 
      HW_Tnext.fit <- predict(HW_2toT, n.ahead = 1, prediction.interval = TRUE)
      HW_2toT.fit.all[,i] <- HW_2toT.fit
      HW_Tnext.fit.all[i] <- HW_Tnext.fit[,1]
    }
  }
  
  if(is.matrix(input.traj) == FALSE){ # on travaille avec les probas/lambdas marginals
    
    nbtimesteps <- dim(input.traj)[1]
    nbsites.others <- dim(input.traj)[2]
    nbsites.interest <- dim(input.traj)[3]
    HW_2toT.fit.all <- array(NA,c((nbtimesteps-2),nbsites.others,nbsites.interest))
    HW_Tnext.fit.all <- matrix(NA,nrow=nbsites.others,ncol=nbsites.interest)
    
    for (i in 1:nbsites.interest){
      for (j in 1:nbsites.others){
        HW_2toT <- stats::HoltWinters(ts(input.traj[,j,i], frequency = 1), gamma = FALSE)
        HW_2toT.fit  <- HW_2toT$fitted[,1]
        HW_Tnext.fit <- predict(HW_2toT, n.ahead = 1, prediction.interval = TRUE)
        HW_2toT.fit.all[,j,i] <- HW_2toT.fit
        HW_Tnext.fit.all[j,i] <- HW_Tnext.fit[,1]
      }
    }
  }
  return(list("HW.fit" = HW_2toT.fit.all,"HW.forecast" = HW_Tnext.fit.all))
}


Compare_stationnary_fonction <- function(Estimated.probability2toTnext.mat, sigma.mat, xi, site.interest = 0, dependence.type="log", dep.param){
  nbtimesteps <- dim(Estimated.probability2toTnext.mat)[1]
  timesteps <- seq(1,nbtimesteps,by=1)
  nbsites <- dim(Estimated.probability2toTnext.mat)[2]
  
  sigma_stationary.vec <- as.numeric(sigma.mat[1,])
  
  if(length(site.interest) > nbsites){
    stop("the number of sites of study exceeds the number of columns of data ")}
  
  if(site.interest[1] == 0){site.interest <- seq(1,dim(Estimated.probability)[2],by=1)} # all the sites where we did the estimation are the interest sites
  
  p.stationary.mat <- matrix(NA, ncol=length(site.interest), nrow = nbtimesteps)
  ratio_site.mat <- matrix(NA, ncol=length(site.interest), nrow = nbtimesteps)
  
  HAT.p.stationary.mat <- matrix(NA, ncol=length(site.interest), nrow = nbtimesteps)
  HAT.ratio_site.mat <- matrix(NA, ncol=length(site.interest), nrow = nbtimesteps)
  
  t_i_minus1 <- (seq(from=1, to=(nbtimesteps), by=1))
  
  for (i in site.interest){
    lambda_stationary.mat <- matrix(NA, ncol=nbsites, nrow = nbtimesteps)
    sigma_intermediate <- (sigma_stationary.vec/sigma_stationary.vec[i])^(1/xi) 
    lambda.stationary.marginal.siteid <- outer(t_i_minus1,sigma_intermediate,FUN = "*")
    
    if (dependence.type == "log"){
      lambda.stationary.marginal.siteid.invSansV <- lambda.stationary.marginal.siteid^(-1)
      lam.aggregated.stationary.siteid <- V(lambda.stationary.marginal.siteid.invSansV,dep.param)
      p.stationary.siteid <- 1/(1+lam.aggregated.stationary.siteid)
      
      HAT.p.stationary.siteid <- (nbsites)/((t_i_minus1)^(dep.param))
      
    }
    
    p.stationary.mat[,i] <- as.numeric(p.stationary.siteid)
    ratio_site.i <- as.numeric(Estimated.probability2toTnext.mat[,i]/p.stationary.siteid)
    ratio_site.mat[,i] <- ratio_site.i # each colum is a site of interest 
    
    HAT.p.stationary.mat[,i] <- as.numeric(HAT.p.stationary.siteid)
    HAT.ratio_site.i <- as.numeric(Estimated.probability2toTnext.mat[,i]/HAT.p.stationary.siteid)
    HAT.ratio_site.mat[,i] <- HAT.ratio_site.i # each colum is a site of interest 
    
  }
  
  return(list("ratio" = ratio_site.mat,
              "p.stationary" = p.stationary.mat,
              "HAT.ratio" = ratio_site.mat,
              "HAT.p.stationary" = p.stationary.mat))
}



############################################################################
find_records_func <- function(Data,marginal=FALSE, same.site = TRUE , other.site = 1){
  nbsites <- dim(Data)[2]
  nbtimesteps <- dim(Data)[1]
  timesteps <- 1:nbtimesteps
  
  record.value <- c()
  record.time <- c()
  record.site <- c()
  
  if (marginal==FALSE){
    marginalRecords.mat <- apply(X=t(Data), MARGIN=1, FUN=cummax)
    TotalRecords.vec <- apply(X=t(marginalRecords.mat), MARGIN=2, FUN=max)
    for (s in 1:nbsites){
      ok <- TotalRecords.vec == Data[,s]
      ok[1] <- FALSE
      if (any(ok)==TRUE){
        record.s.value.vec <- Data[ok,s]
        record.s.time.vec <- timesteps[ok]
        record.value <- append(record.value,record.s.value.vec)
        record.time <- append(record.time,record.s.time.vec)
        record.site <- append(record.site,rep(s,length(record.s.time.vec)))
      }
    }
  }
  
  if (marginal==TRUE){
    marginalRecords.mat <- apply(X=t(Data), MARGIN=1, FUN=cummax)
    if (same.site == TRUE){
      for (s in 1:nbsites){
        ok <- marginalRecords.mat[,s] == Data[,s]
        ok[1] <- FALSE
        if (any(ok)==TRUE){
          record.s.value.vec <- Data[ok,s]
          record.s.time.vec <- timesteps[ok]
          record.value <- append(record.value,record.s.value.vec)
          record.time <- append(record.time,record.s.time.vec)
          record.site <- append(record.site,rep(s,length(record.s.time.vec)))
        }
      }
    }
    
    
    if (same.site == FALSE){
      for (s in 1:nbsites){
        ok <- marginalRecords.mat[,other.site] <= Data[,s]
        ok[1] <- FALSE
        if (any(ok)==TRUE){
          record.s.value.vec <- Data[ok,s]
          record.s.time.vec <- timesteps[ok]
          record.value <- append(record.value,record.s.value.vec)
          record.time <- append(record.time,record.s.time.vec)
          record.site <- append(record.site,rep(s,length(record.s.time.vec)))
        }
      }
    }
    
    
  }
  
  record.table <- cbind(record.time,record.value,record.site)
  
  return(record.table)
}

PlotObservedRecord.separate <- function(Data, siteid = 0, par1 = 4, par2 = 2 , 
                                        add.max=FALSE , max.type="all", add.record.total = FALSE, add.record.marginal = FALSE ){ # if  siteid = 0 we plot data in all sites 
  par(mfrow=c(par1,par2))
  timesteps <- dim(Data)[1]
  nbsites <- dim(Data)[2]
  
  # If we chose to plot observations in one single site
  if (length(siteid)==1){
    if (siteid>0){
      par1 <- 1
      par2 <- 1 
      plot(c(1:timesteps),Data[,siteid],ylim=c(min(Data),max(Data)),
           xlab="time",ylab='observations', main=bquote("Observations at site "~ .(siteid)))
      if(add.record.marginal == TRUE){
        MarginalRecord.matrix <- find_records_func(Data,marginal=TRUE)
        MarginalRecord.siteid <- MarginalRecord.matrix[MarginalRecord.matrix[,3] == siteid,]
        points(MarginalRecord.siteid[,1],MarginalRecord.siteid[,2],col="black", bg="wheat2",pch=21)
      }
      if(add.record.total == TRUE){
        TotalRecord.matrix <- find_records_func(Data)
        TotalRecord.siteid <- TotalRecord.matrix[TotalRecord.matrix[,3] == siteid,]
        points(TotalRecord.siteid[,1],TotalRecord.siteid[,2],col="black", bg="blue",pch=21)
      }
      if(add.max==TRUE){
        max.traj <- apply(X=Data, MARGIN=1, FUN=max)
        points(c(1:timesteps),max.traj,col='red',pch=4)
        if(add.record.total==TRUE){
          points(TotalRecord.siteid[,1],TotalRecord.siteid[,2],col="blue",pch=4)
        }
      }
      if(add.max==TRUE & add.record.total==TRUE & add.record.marginal == TRUE){
        legend('bottomright',legend=c("max all sites", "record all sites", "marginal record"),col=c("red","blue","gray"),pch=c(4,4,16), bg="transparent")
      }
      if(add.max==TRUE & add.record.total==TRUE & add.record.marginal == FALSE){
        legend('bottomright',legend=c("max all sites", "record all sites"),col=c("red","blue"),pch=c(4,4), bg="transparent")
      }
      if(add.max==TRUE & add.record.total==FALSE & add.record.marginal == TRUE){
        legend('bottomright',legend=c("max all sites","marginal record"),col=c("red","wheat2"),pch=c(4,16), bg="transparent")
      }
      if(add.max==TRUE & add.record.total==FALSE & add.record.marginal == FALSE){
        legend('bottomright',legend="max all sites",col="red",pch=4, bg="transparent")
      }
    }
    
    if (siteid==0){
      ylim.max = max(Data)
      ylim.min = min(Data)
      
      if(add.record.marginal == TRUE){
        MarginalRecord.matrix <- find_records_func(Data,marginal=TRUE)
      }
      if(add.record.total == TRUE){
        TotalRecord.matrix <- find_records_func(Data)
      }
      
      for(i in 1:nbsites){
        plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="time",ylab='observations',
             main=bquote("Observations at site "~ .(i)))
        if(add.record.marginal == TRUE){
          MarginalRecord.i <- MarginalRecord.matrix[MarginalRecord.matrix[,3] == i,]
          points(MarginalRecord.i[,1],MarginalRecord.i[,2],col="black", bg="wheat2",pch=21)
        }
        if(add.record.total == TRUE){
          if (i %in%  TotalRecord.matrix[,3] == TRUE ){
            TotalRecord.i <- TotalRecord.matrix[TotalRecord.matrix[,3] == i,]
            TotalRecord.i <- as.matrix(TotalRecord.i)
            if(dim(TotalRecord.i)[2]==1){TotalRecord.i <- t(TotalRecord.i)}
            points(as.numeric(TotalRecord.i[,1]),as.numeric(TotalRecord.i[,2]),col="black", bg="blue",pch=21)
          }
        }
        
        if(add.max==TRUE){
          max.traj <- apply(X=Data, MARGIN=1, FUN=max)
          points(c(1:timesteps),max.traj,col='red',pch=4)
          if(add.record.total==TRUE){
            if (i %in%  TotalRecord.matrix[,3] == TRUE ){
              TotalRecord.i <- TotalRecord.matrix[TotalRecord.matrix[,3] == i,]
              TotalRecord.i <- as.matrix(TotalRecord.i)
              if(dim(TotalRecord.i)[2]==1){TotalRecord.i <- t(TotalRecord.i)}
              points(as.numeric(TotalRecord.i[,1]),as.numeric(TotalRecord.i[,2]),col="blue",pch=4)
            }
          }
        }
        if(add.max==TRUE & add.record.total==TRUE & add.record.marginal == TRUE){
          legend('bottomright',legend=c("max all sites", "record all sites", "marginal record"),col=c("red","blue","wheat2"),pch=c(4,4,16), bg="transparent")
        }
        if(add.max==TRUE & add.record.total==TRUE & add.record.marginal == FALSE){
          legend('bottomright',legend=c("max all sites", "record all sites"),col=c("red","blue"),pch=c(4,4), bg="transparent")
        }
        if(add.max==TRUE & add.record.total==FALSE & add.record.marginal == TRUE){
          legend('bottomright',legend=c("max all sites","marginal record"),col=c("red","wheat2"),pch=c(4,16), bg="transparent")
        }
        if(add.max==TRUE & add.record.total==FALSE & add.record.marginal == FALSE){
          legend('bottomright',legend="max all sites",col="red",pch=4, bg="transparent")
        }
      }
    }
  }
  
  
  if (length(siteid)>1){
    
    nb.plots <- length(siteid)
    ylim.max = max(Data)
    ylim.min = min(Data)
    
    # plotted record will also depend on "max.type"
    if(add.record.marginal == TRUE){
      MarginalRecord.matrix <- find_records_func(Data,marginal=TRUE)
    }
    if(add.record.total == TRUE){
      if(max.type=="all"){
        TotalRecord.matrix <- find_records_func(Data)
      }
      if(max.type=="not.all"){
        TotalRecord.matrix <- find_records_func(Data[,siteid])
        nbcol <- 1:dim(Data[,siteid])[2]
      }
    }
    
    colrecord <- "blue"
    
    
    if(add.max==TRUE){
      if(max.type=="all"){
        max.traj <- apply(X=Data, MARGIN=1, FUN=max)
        colmax <- "red"
        colrecord <- "blue"
        legendmax.text <- "max all sites" 
        legendrecord.text <- "record all sites"}
      
      if(max.type=="not.all"){
        max.traj <- apply(X=Data[,c(siteid)], MARGIN=1, FUN=max)
        colmax <- "green"
        colrecord <- "cyan"
        legendmax.text <- "max showed sites"
        legendrecord.text <- "record showed sites"}
      
      j <- 1
      for(i in siteid ){
        plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="time",ylab='observations',
             main=bquote("Observations at site "~ .(i)))
        if(add.record.marginal == TRUE){
          MarginalRecord.i <- MarginalRecord.matrix[MarginalRecord.matrix[,3] == i,]
          points(MarginalRecord.i[,1],MarginalRecord.i[,2],col="black", bg="gray",pch=21)
        }
        if(add.record.total == TRUE){
          TotalRecord.i <- TotalRecord.matrix[TotalRecord.matrix[,3] == j,]
          points(TotalRecord.i[,1],TotalRecord.i[,2],col="black", bg=colrecord,pch=21)
        }
        if(add.max==TRUE){
          points(c(1:timesteps),max.traj,col=colmax,pch=4)
          if(add.record.total==TRUE){
            points(TotalRecord.i[,1],TotalRecord.i[,2],col=colrecord,pch=4)
          }
        }
        
        if(add.max==TRUE & add.record.total==TRUE & add.record.marginal == TRUE){
          legend('bottomright',legend=c(legendmax.text, legendrecord.text, "marginal record"),col=c("red",colrecord,"gray"),pch=c(4,4,16), bg="transparent")
        }
        if(add.max==TRUE & add.record.total==TRUE & add.record.marginal == FALSE){
          legend('bottomright',legend=c(legendmax.text, legendrecord.text),col=c(colmax,colrecord),pch=c(4,4), bg="transparent")
        }
        if(add.max==TRUE & add.record.total==FALSE & add.record.marginal == TRUE){
          legend('bottomright',legend=c(legendmax.text,"marginal record"),col=c(colmax,"gray"),pch=c(4,16), bg="transparent")
        }
        if(add.max==TRUE & add.record.total==FALSE & add.record.marginal == FALSE){
          legend('bottomright',legend=legendmax.text,col=colmax,pch=4, bg="transparent")
        } 
        
        if(add.record.total==TRUE & add.record.marginal == TRUE){
          legend('bottomright',legend=c("total record","marginal record"),col=c("blue",'gray'),pch=16, bg="transparent")
        } 
        
        j <- j+1 
      }
      
    }
    
  }
}


GiYi.func <- function(Data,h){
  nbsites <- dim(Data)[2]
  nbtimesteps <- dim(Data)[1]
  timesteps <- c(1:nbtimesteps)
  GiYi.matrix <- matrix(NA,ncol=nbsites,nrow=nbtimesteps)
  if(length(h)==1){
    h <- rep(h,nbsites)}
  if(length(h) != nbsites){
    stop("You are not using the correct Bandwidth length, a value or a vector of the legth of number of colums in matrix ")}
  for (s in 1:nbsites){
    Ghat.kernel <- G_estimator(timesteps, Data[,s], h[s])
    for (ti in 1:nbtimesteps){
      GiYi <- Ghat.kernel(t0 = ti, x = Data[ti,s])
      GiYi.matrix[ti,s] <- GiYi
    }
  }
  return(GiYi.matrix)
}

haversine_distance <- function(lon1, lat1, lon2, lat2) {
  # Earth radius in kilometers
  R <- 6371
  
  # Convert degrees to radians
  lon1 <- lon1 * (pi / 180)
  lat1 <- lat1 * (pi / 180)
  lon2 <- lon2 * (pi / 180)
  lat2 <- lat2 * (pi / 180)
  
  # Haversine formula
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  distance <- R * c
  
  return(distance)
}

################################## # FUNCTION POUR ESTIMER DEP QUI TOURNER ICI
Version2_V.applied.lambda.estim.func<- function(GiYi.mat, lambda , use.timesteps = 0 ){ 
  if(is.numeric(lambda)==TRUE){lambda <- as.matrix(lambda,nrow=1)}
  nb.evaluations <- dim(lambda)[1]
  vw <- rep(NA,nb.evaluations)
  cw <- rep(NA,nb.evaluations)
  Aw <- rep(NA,nb.evaluations)
  Vw <- rep(NA,nb.evaluations)
  lambdaV <- rep(NA,nb.evaluations)
  
  #zi.mat  <- -1/(log((lambda)^(-1))) 
  # zi.mat  <- (-log(lambda))^(-1)
  # zi.mat  <- (lambda)^(-1) # maybe it is not like that, I don't remember, LOOK FOR IT LATER
  #zi.mat  <- 1/log(lambda)
  #zi.mat  <- 1/exp(-lambda)
  #zi.mat  <- exp(-lambda)
  tt<-1:nb.evaluations
  #zi.mat  <- 1/(lambda*tt)
  zi.mat  <- 1/lambda
  
  w.mat <- t(apply(zi.mat,1, function(x) x/sum(x)))
  # w.mat <- zi.mat/(apply(zi.mat,1,sum)) # maybe it is not like that, I don't remember, LOOK FOR IT LATER
  
  GiYi.inform  <- GiYi.mat
  #GiYi.inform  <- exp(-1/GiYi.mat)
  # GiYi.inform  <- exp(GiYi.mat)
  
  # use.timesteps = 0 means we will use all the information to construct V (we may want to use only de central time steps)
  if (use.timesteps[1] != 0 ){
    GiYi.inform <- GiYi.mat[use.timesteps,] } # truncate the matrix used
  nbsites <- dim(GiYi.inform)[2]
  nbcolumns <- dim(GiYi.inform)[1]
  
  for(ti in 1:nb.evaluations){ 
    wi.eval <- as.numeric(w.mat[ti,])
    FiYiwi <- t(t(GiYi.inform)^(1/wi.eval))
    Mean.FiYiwi<- rowMeans(FiYiwi) # a vector of length nbcolumns
    Max.FiYiwi<- apply(FiYiwi,1,max) # a vector of length nbcolumns
    vw.ti <- mean(Max.FiYiwi - Mean.FiYiwi)
    vw[ti] <- vw.ti
    cw.ti <- mean(wi.eval/(1+wi.eval))
    cw[ti] <- cw.ti
    Aw.ti <- (vw.ti + cw.ti)/ (1 - vw.ti - cw.ti)
    Aw[ti] <- Aw.ti
    sum.zi.ti <- sum(zi.mat[ti,]^(-1)) # wi = xi/(sumxi)
    lambdaV[ti] <- sum.zi.ti * Aw.ti # lambdaV.ti
  }
  lambdaV <- lambdaV*tt
  return(list("lambdaV"=lambdaV,"A"=Aw))
}

#############################@
import_SQRT_TX_infromation <- function(min_years = 55, plot_map=FALSE ,get_notused=FALSE, plot_all=FALSE, just_plot =FALSE ){
  # library(dplyr) # <- needs package dlypr 
  
  # import data 
  stations_corse_names <- c("CAP PERTUSATO","AJACCIO","CALVI","BASTIA")
  SQR_TX.Information <- read.csv("/Users/pgonzale/Downloads/Liste_SQR_TX_metro.csv",sep = ";",skip=2) %>%
    mutate(year.start = floor(date_debut_SQR.YYYYMMDD./ 10000), month.start = floor((date_debut_SQR.YYYYMMDD. %% 10000) / 100), day.start = date_debut_SQR.YYYYMMDD. %% 100) %>% 
    select(-date_debut_SQR.YYYYMMDD.) %>% mutate(year.end = floor(date_fin_SQR.YYYYMMDD./ 10000), month.end = floor((date_fin_SQR.YYYYMMDD. %% 10000) / 100), day.end = date_fin_SQR.YYYYMMDD. %% 100) %>% select(-date_fin_SQR.YYYYMMDD.) %>%  
    filter(!(nom_usuel %in% stations_corse_names))
  
  # min_years is the minimum length for the time series we will use. but some site might have more data
  
  FirstYear_data <- min(SQR_TX.Information$year.start) # first first observation year of all data (normally 1951)
  LastYear_data <- max(SQR_TX.Information$year.end) # latest last observation year of all data (normally 2022)
  
  first_year_filter <- LastYear_data-min_years # latest first observation year of all data (if min_years = 55 -> it is 1967)
  
  stations_notused_names <- subset(SQR_TX.Information, year.start > first_year_filter , select = nom_usuel)[[1]]
  nbstations_notused <- length(stations_notused_names) # number of removed stations (if min_years = 55 -> it is 40)
  
  # information of the kept stations 
  StudyZone_SQR_TX.Information <- SQR_TX.Information %>%  filter(!(nom_usuel %in% stations_notused_names))
  nbsitesused <- nrow(StudyZone_SQR_TX.Information) # number of kept stations (if min_years = 55 -> it is 100)
  
  if (plot_map==TRUE & get_notused==FALSE ){
    
    if (plot_all==FALSE){
      addmaptext <-paste(nbsitesused,"stations")
      france_map <- map_data("france")
      france_stations <- ggplot() +
        geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey50', size = 1) +
        coord_quickmap() +
        geom_point(data = StudyZone_SQR_TX.Information, aes(x = longitude..., y = latitude...), fill="seagreen1", col = "grey2", size = 1,pch=21) + 
        geom_text(aes(x = -2, y = 42,label = addmaptext),stat = "unique") +
        labs(title = "Meteorological stations")
      print(france_stations)
    }
    
    if (plot_all==TRUE){
      addmaptext <-paste(nbsitesused,"stations")
      france_map <- map_data("france")
      france_stations_used <- ggplot() +
        geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey50', size = 1) +
        coord_quickmap() +
        geom_point(data = StudyZone_SQR_TX.Information, aes(x = longitude..., y = latitude...), fill="seagreen1", col = "grey2", size = 1,pch=21) + 
        geom_text(aes(x = -2, y = 42,label = addmaptext),stat = "unique")# +
      #labs(title = "Meteorological stations")
      #print(france_stations_used)
      nballstations <- nbsitesused + nbstations_notused
      addmaptext_all <-paste(nballstations,"stations")
      france_map <- map_data("france")
      france_stations_all <- ggplot() +
        geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey50', size = 1) +
        coord_quickmap() +
        geom_point(data = SQR_TX.Information, aes(x = longitude..., y = latitude...), fill="seagreen1", col = "grey2", size = 1,pch=21) + 
        geom_text(aes(x = -2, y = 42,label = addmaptext_all),stat = "unique") #+
      #labs(title = "Meteorological stations")
      #print(france_stations_all)
      library(ggpubr)
      plotmap <-ggarrange(
        france_stations_used, france_stations_all, labels = c("A", "B"),
        common.legend = TRUE, legend = "bottom"
      )
      print(plotmap)
    }
  }
  
  if (get_notused==TRUE ){
    # information of the removed stations 
    NotusedZone_SQR_TX.Information <- SQR_TX.Information %>%  filter((nom_usuel %in% stations_notused_names))
  }
  
  if (plot_map==TRUE & get_notused==TRUE){
    if (plot_all==FALSE){
      addmaptext <-paste(nbsitesused,"stations used ") # blue the used , red not used
      france_map <- map_data("france")
      france_stations_used<- ggplot() +
        geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey50', size = 1) +
        coord_quickmap() +
        geom_point(data = SQR_TX.Information, aes(x = longitude..., y = latitude...), fill="red", col = "grey2", size = 1,pch=21) + 
        geom_point(data = StudyZone_SQR_TX.Information, aes(x = longitude..., y = latitude...), fill="seagreen1", col = "grey2", size = 1,pch=21)+
        geom_text(aes(x = -1, y = 42,label = addmaptext),stat = "unique") +
        labs(title = "Statios Meteo France")
      print(france_stations_used)
    }
    if (plot_all==TRUE){
      addmaptext <-paste(nbsitesused,"stations used ") # blue the used , red not used
      france_map <- map_data("france")
      france_stations_used <- ggplot() +
        geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey50', size = 1) +
        coord_quickmap() +
        geom_point(data = StudyZone_SQR_TX.Information, aes(x = longitude..., y = latitude...), fill="seagreen1", col = "grey2", size = 1,pch=21)+
        geom_text(aes(x = -1, y = 42,label = addmaptext),stat = "unique") # +
      #labs(title = "Statios Meteo France")
      # print(plotmap)
      nballstations <- nbsitesused + nbstations_notused
      addmaptext_all <-paste(nballstations,"stations")
      france_map <- map_data("france")
      france_stations_all <- ggplot() +
        geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey50', size = 1) +
        coord_quickmap() +
        geom_point(data = SQR_TX.Information, aes(x = longitude..., y = latitude...), fill="red", col = "grey2", size = 1,pch=21) + 
        geom_point(data = StudyZone_SQR_TX.Information, aes(x = longitude..., y = latitude...), fill="seagreen1", col = "grey2", size = 1,pch=21)+
        geom_text(aes(x = -2, y = 42,label = addmaptext_all),stat = "unique") #+
      #labs(title = "Weather stations")
      #print(france_stations_all)
      plotmap <-ggarrange(
        france_stations_used, france_stations_all, labels = c("A", "B"),
        common.legend = TRUE, legend = "bottom"
      )
      print(plotmap)
    }
  }
  
  if (get_notused==FALSE  & just_plot==FALSE  ){
    return(StudyZone_SQR_TX.Information)
  }
  if (get_notused==TRUE & just_plot==FALSE ){
    return(list("StudyZone_SQR_TX.Information" = StudyZone_SQR_TX.Information,
                "NotusedZone_SQR_TX.Information " = NotusedZone_SQR_TX.Information ))
  }
}
import_files <- function(data_frame, folder_path,print_missing=FALSE){
  file_names <-data_frame$nom_fichier
  nbsitesused <- length(file_names) 
  list.dataframes.usedsites <- vector("list", length = nbsitesused)
  
  for (i in seq_along(file_names)){ 
    file_path <- paste(folder_path,file_names[i],sep="")
    datafile <- read.csv(file_path,sep = ";",skip=8)
    datafile <- datafile %>% mutate(year = floor(AAAAMMJJ/ 10000), month = floor((AAAAMMJJ %% 10000) / 100), day = AAAAMMJJ %% 100) %>% select(-AAAAMMJJ)
    site.firstyear <- min(datafile$year)
    
    dataMAX  <- datafile %>% group_by(year) %>% summarize(maxTX = max(VALEUR,na.rm = TRUE))
    nbyears_dataMAX <- dim(dataMAX)[1] 
    
    if(!any(dplyr::filter(datafile, year == site.firstyear)$month<= 8)){ # not use first year if it begins after summer
      datafile <- datafile %>%  filter(!(year == site.firstyear))
    } 
    
    summer.dataframe <- dplyr::filter(datafile, month <= 8 & month >= 6)
    summer.observations <- summer.dataframe$VALEUR
    
    if (any(is.na(summer.observations))){
      if (print_missing==TRUE){
        print("missing values on summer")
      }
    }
    
    if(any(!is.na(summer.observations)) == FALSE){ # if all values are missing 
      if (print_missing==TRUE){
        print("all summer vales are mising, site erased") # here we jaut filter column , but we should make the difference between first and other years
      }
      dataMAX$VALUER <- rep(NA,nbyears_dataMAX)
    }
    list.dataframes.usedsites[[i]] <- dataMAX # we get a list containing all the sites (col1=year,col2=maxvalue)
  }
  return(list.dataframes.usedsites)
}

Get_observation.matrix <- function(data_frame, folder_path){
  FirstYear_data <- min(data_frame$year.start)
  LastYear_data <- max(data_frame$year.end)
  nbsitesused <- nrow(data_frame)
  Usedsites.max.summary_matrix <- matrix(ncol=nbsitesused,nrow=(LastYear_data -FirstYear_data)+1)
  rownames(Usedsites.max.summary_matrix) <- FirstYear_data:LastYear_data
  colnames(Usedsites.max.summary_matrix) <- 1:nbsitesused
  list_of_dataframes <- import_files(data_frame = StudyZone_dataframe, folder_path= folder_path)
  for(i in 1:nbsitesused){
    mat1 <- as.matrix(list_of_dataframes[[i]])
    mat.obs <- mat1[,-1]
    rownames(mat1) <- mat1[,1]
    Usedsites.max.summary_matrix[rownames(mat1),i] <- mat.obs 
  }
  return(Usedsites.max.summary_matrix)
}

Get_GiYi.matrix <- function(Usedsites_max.matrix,bandwidth.vec){
  nbsitesused <- ncol(Usedsites_max.matrix)
  nbyears <- nrow(Usedsites_max.matrix)
  FirstYear_data <- as.numeric(rownames(Usedsites_max.matrix)[1])
  LastYear_data <- as.numeric(rownames(Usedsites_max.matrix)[nbyears])
  Usedsites.GiYi.matrix <- matrix(ncol=nbsitesused,nrow=(LastYear_data -FirstYear_data)+1)
  rownames(Usedsites.GiYi.matrix) <- FirstYear_data:LastYear_data
  colnames(Usedsites.GiYi.matrix) <- 1:nbsitesused
  if (length(bandwidth.vec)==1){bandwidth.vec <- rep(bandwidth.vec,nbsitesused)}
  for(i in 1:nbsitesused){
    yearsAll <- rownames(Usedsites_max.matrix)
    Yi.allyears <- Usedsites_max.matrix[,i]
    Yi.and.years <- cbind(yearsAll,Yi.allyears)
    # Yi <- Yi[!is.na(Yi)]
    Yi.and.years <- Yi.and.years[!is.na(Yi.and.years[,2]),]
    observationYears <- Yi.and.years[,1]
    Yi <- as.numeric(Yi.and.years[,2])
    GiYi <- GiYi.func(Data = as.matrix(Yi), h = bandwidth.vec[i])
    GiYi.and.years <- cbind(observationYears,as.numeric(GiYi))
    GiYi.obs <- as.numeric(GiYi.and.years[,-1])
    rownames(GiYi.and.years) <- GiYi.and.years[,1]
    Usedsites.GiYi.matrix[rownames(GiYi.and.years),i] <- GiYi.obs 
  }
  return(Usedsites.GiYi.matrix)
}

Get_biv_madogram <- function(GiYi_matrix){
  nbsitesused <- ncol(GiYi_matrix)
  sites_seq <- 1:nbsitesused
  mat.dij <- matrix(NA,nrow=nbsitesused,ncol=nbsitesused)
  mat.startyear.ij <- matrix(NA,nrow=nbsitesused,ncol=nbsitesused)
  for ( i in sites_seq){ # create matrix of dissimilarities
    for ( j in sites_seq){
      site1site2 <- GiYi_matrix[,c(i,j)]
      site1site2 <- as.matrix(na.omit(site1site2))
      mat.startyear.ij[i,j] <- rownames(site1site2)[1]
      dij <- mean(abs(site1site2[,1] - site1site2[,2]))/2
      mat.dij[i,j] <- dij # la diagonal doit être 0 et j,i = i,j
    }
  }
  return(list("dissimilarity" = mat.dij,
              "years" = mat.startyear.ij))
}

library(cluster)
Get_clusters <- function(k=8,diss_ij.matrix, StudyDataframe_info){
  nbsitesused <- nrow(diss_ij.matrix)
  pamClusters <- pam(x = diss_ij.matrix , k=k, diss=TRUE)
  cluster_id <- pamClusters$clustering
  StudyZone.Information_clusteid <- cbind(StudyDataframe_info, cluster_id )
  colorCluster <- rep(NA,nbsitesused)
  sites_seq <- 1:nbsitesused
  colorVec <- rainbow(k, s = 1, v = 1)
  if(k==8){
    colorVec <- c("#FF0000","#FFA700FF","#FFFF00","#00FF40", "#00FFFF", "#0040FF", "#8000FF", "#FF00BF")
  }
  if(k==7){
    colorVec <- c("#FF0000","#FFFF00","#00FF40", "#00FFFF", "#0040FF", "#8000FF", "#FF00BF")
  }
  if(k==6){
    colorVec <- c("#FF0000","#FFFF00","#00FF40", "#00FFFF", "#8000FF", "#FF00BF")
  }
  if(k==5){
    colorVec <- c("#FF0000","#FFFF00","#00FF40", "#00FFFF", "#8000FF")
  }
  for (i in sites_seq){
    clusterid <- StudyZone.Information_clusteid$cluster_id[i]
    colorCluster[i] <- colorVec[clusterid]
  }
  StudyZone.Information_clusteid_color <- cbind(StudyZone.Information_clusteid,colorCluster)
  medoids_id_site <- pamClusters$id.med
  StudyZone_medoids.Information_clusteid_color <- StudyZone.Information_clusteid_color %>% slice(medoids_id_site)
  StudyZone_medoids.Information_clusteid_color <- cbind(StudyZone_medoids.Information_clusteid_color,medoids_id_site)
  return(list("all_stations" =StudyZone.Information_clusteid_color,
              "only_medoids" = StudyZone_medoids.Information_clusteid_color))
}

plotmap_stations_partition <- function(StudyZone_Info.dataframe, StudyMedoids_Info.dataframe, k = 8) {
  colorVec <- c("#FF0000", "#FFA700FF", "#FFFF00", "#00FF40", "#00FFFF", "#0040FF", "#8000FF", "#FF00BF")
  #addmaptext <-paste(nbsitesused,"stations")
  france_map <- map_data("france")
  map_clustered <- ggplot() +
    geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey60', size = 1) +
    coord_quickmap() +
    geom_point(data = StudyZone_Info.dataframe, aes(x = longitude..., y = latitude..., fill = factor(colorCluster)), color = "grey2", pch = 21, size = 1) +
    scale_fill_manual(values = colorVec, name = "Clusters", labels = 1:k) +
    geom_point(data = StudyMedoids_Info.dataframe, aes(x = longitude..., y = latitude...), size = 2, shape = 8) +
    theme(legend.position = "bottom")+
    labs(title = "Weather stations")
  print(map_clustered)
}

get_color_name <- function(color_hex) {
  library(plotrix)
  out = vector(length = length(color_hex))
  
  for(i in seq_along(color_hex)) {
    out[i] <- color.id(color_hex[i])[1]
  }
  out
  
}

Get_dataframeInfo_cluster <- function(k.id = ClusterID,StudyZone.dataframe,StudyMedoids.dataframe,dissimilaity.matrix){
  StudyZone_cluster_kid.Information <- StudyZone.dataframe %>%  filter((cluster_id == k.id))
  Medoid_cluster_kid.Information<- slice(StudyMedoids.dataframe, k.id)
  medoid_kid<- Medoid_cluster_kid.Information$medoids_id_site
  indexsites_inCluster <- which(StudyZone.dataframe$cluster_id==k.id)
  StudyZone_cluster_kid.Information <- cbind(StudyZone_cluster_kid.Information,indexsites_inCluster)
  dij_to_medoid <- as.numeric(dissimilaity.matrix[medoid_kid,indexsites_inCluster])
  StudyZone_cluster_kid.Information <- cbind(StudyZone_cluster_kid.Information,dij_to_medoid)
  StudyZone_cluster_kid.Information <- arrange(StudyZone_cluster_kid.Information, dij_to_medoid)
  return(StudyZone_cluster_kid.Information)
}

plotmap_stations_cluster <- function(Cluster.Info,label=FALSE,label.dij=FALSE){
  medoid.Info <- Cluster.Info[1,]
  #medoid.Info <- Cluster.Info[,14]
  color.kid <- factor(Cluster.Info$colorCluster)[1]
  color.id <- get_color_name(color.kid)
  id.kid <- medoid.Info$cluster_id
  if (label==FALSE & label.dij==FALSE){
    france_map <- map_data("france")
    mapcluster <- ggplot() +
      geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey60', size = 1) +
      coord_quickmap() +
      geom_point(data = Cluster.Info, aes(x = longitude..., y = latitude...), fill= color.id , color = "grey2", pch = 21, size = 1) +
      geom_point(data =medoid.Info, aes(x = longitude..., y = latitude...), size = 2, shape = 8)+
      theme(legend.position="none") +
      labs(title = paste("Study Zone",id.kid ))
  }
  
  if (label==TRUE & label.dij==FALSE){
    library(ggpubr)
    library(ggrepel)
    france_map <- map_data("france")
    mapcluster <- ggplot() +
      geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey60', size = 1) +
      coord_quickmap() +
      geom_point(data = Cluster.Info, aes(x = longitude..., y = latitude...), fill= color.id , color = "grey2", pch = 21, size = 1) +
      geom_point(data = medoid.Info, aes(x = longitude..., y = latitude...), size = 2, shape = 8)+
      geom_text_repel(data=Cluster.Info, aes(x = longitude..., y = latitude..., label=1:nrow(Cluster.Info),fontface = "bold"),size = 2.9, color="blue")+
      theme(legend.position="none")+
      labs(title = paste("Study Zone",id.kid ))
  }
  
  
  if (label==TRUE & label.dij==TRUE){
    library(ggpubr)
    library(ggrepel)
    france_map <- map_data("france")
    mapcluster <- ggplot() +
      geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey60', size = 1) +
      coord_quickmap() +
      geom_point(data = Cluster.Info, aes(x = longitude..., y = latitude...), fill= color.id , color = "grey2", pch = 21, size = 1) +
      geom_point(data = medoid.Info, aes(x = longitude..., y = latitude...), size = 2, shape = 8)+
      geom_text_repel(data=Cluster.Info, aes(x = longitude..., y = latitude..., label=round(dij_to_medoid,3),fontface = "bold"),size = 2.8, color="blue")+
      theme(legend.position="none")+
      labs(title = paste("Study Zone",id.kid ))
  }
  print(mapcluster)
}

plot_dij_spread_intraCluster <- function(Cluster.Info, add.lines=FALSE, add.map=FALSE){
  color.kid <- factor(Cluster.Info$colorCluster)[1]
  color.id <- get_color_name(color.kid)
  if (add.lines == FALSE & add.map==FALSE){
    par(mfrow = c(1, 1))
    plot(1:nrow(Cluster.Info),Cluster.Info$dij_to_medoid,
         xlab= "stations", ylab="Dissimilarity",main = "Dissimilarity  w.r.t Medoid",
         col="grey2", bg=color.id,pch=21)
  }
  if (add.lines == TRUE & add.map==FALSE){
    par(mfrow = c(1, 1))
    plot(1:nrow(Cluster.Info),Cluster.Info$dij_to_medoid,
         xlab= "stations", ylab="Dissimilarity",main = "Dissimilarity  w.r.t Medoid",
         col="grey2", bg=color.id,pch=21)
    mm <- median(Cluster.Info$dij_to_medoid[-1])
    mad <- median(abs(Cluster.Info$dij_to_medoid[-1] - mm))
    xline.mm <- mad+mm
    abline(h=xline.mm, col="grey")
    yline.mm <- which(Cluster.Info$dij_to_medoid >= xline.mm)[1]
    abline(v=yline.mm, col="grey")
  }
  
  if (add.lines == FALSE & add.map==TRUE){
    #par(mfrow=c(1,2), mar = c(5, 5.3, 4, 2))
    library(gridExtra)
    medoid.Info <- StudyZone_cluster_kid.Information[1,]
    #medoid.Info <- Cluster.Info[,14]
    color.kid <- medoid.Info$colorCluster
    color.id <- get_color_name(color.kid)
    id.kid <- medoid.Info$cluster_id
    france_map <- map_data("france")
    mapcluster <- ggplot() +
      geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey60', size = 1) +
      coord_quickmap() +
      geom_point(data = Cluster.Info, aes(x = longitude..., y = latitude...),fill= color.id, color = "grey2", pch = 21, size = 1) +
      geom_point(data = medoid.Info, aes(x = longitude..., y = latitude...), size = 2, shape = 8)+
      geom_text_repel(data=Cluster.Info, aes(x = longitude..., y = latitude..., label=1:nrow(Cluster.Info),fontface = "bold"),size = 2.9, color="blue")+
      theme(legend.position="none")
    
    plot.dij <-  ggplot(Cluster.Info, aes(x = factor(1:nrow(Cluster.Info)), y = dij_to_medoid)) +
      geom_point(color = "grey2", fill= color.id, shape = 21) +
      labs(x = "stations", y = "Dissimilarity") +
      theme(legend.position="none")
    
    plist <- list(mapcluster,plot.dij)
    grid.arrange(grobs = plist, nrow = 2)
    
    
  }
  
  
  if (add.lines == TRUE & add.map==TRUE){
    #par(mfrow=c(1,2), mar = c(5, 5.3, 4, 2))
    #library(gridExtra)
    medoid.Info <- StudyZone_cluster_kid.Information[1,]
    #medoid.Info <- Cluster.Info[,14]
    
    id.kid <- medoid.Info$cluster_id
    france_map <- map_data("france")
    mapcluster <- ggplot() +
      geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey60', size = 1) +
      coord_quickmap() +
      geom_point(data = Cluster.Info, aes(x = longitude..., y = latitude...),fill= color.id, color = "grey2", pch = 21, size = 1) +
      geom_point(data = medoid.Info, aes(x = longitude..., y = latitude...), size = 2, shape = 8)+
      geom_text_repel(data=Cluster.Info, aes(x = longitude..., y = latitude..., label=1:nrow(Cluster.Info),fontface = "bold"),size = 2.9, color="blue")+
      theme(legend.position="none")
    
    # baseplot.dij <- plot(1:nrow(Cluster.Info),Cluster.Info$dij_to_medoid,
    #      xlab= "stations", ylab="Dissimilarity",
    #      col="grey2", bg=factor(Cluster.Info$colorCluster[1]),pch=21)
    mm <- median(Cluster.Info$dij_to_medoid[-1])
    mad <- median(abs(Cluster.Info$dij_to_medoid[-1] - mm))
    xline.mm <- mad+mm
    # abline(h=xline.mm, col="grey")
    yline.mm <- which(Cluster.Info$dij_to_medoid >= xline.mm)[1]
    # abline(v=yline.mm, col="grey")
    
    # Create ggplot
    labelpoints <- StudyZone_cluster_kid.Information$indexsites_inCluster
    plot.dij <-  ggplot(Cluster.Info, aes(x = factor(1:nrow(Cluster.Info)), y = dij_to_medoid)) +
      geom_point(color = "grey2", fill= color.id,shape = 21) +
      geom_text(label=labelpoints, nudge_x = 0.0025, nudge_y = -0.005,size =3, max.overlaps=Inf)+
      geom_hline(yintercept = xline.mm, color = "grey") +
      geom_vline(xintercept = yline.mm, color = "grey") +
      labs(x = "stations", y = "Dissimilarity") +
      theme(legend.position="none")
    
    
    plist <- list(mapcluster,plot.dij)
    grid.arrange(grobs = plist, nrow = 2)
  }
  
}

Get_observations.StudyZone <- function(obs.matrix,StudyZone_cluster_kid.Information, remove.years=FALSE){
  obs_studyzone.matrix <- obs.matrix[,StudyZone_cluster_kid.Information$indexsites_inCluster]
  if (remove.years==TRUE){
    last.first.year <- max(StudyZone_cluster_kid.Information$year.start)
    obs_studyzone.matrix <- obs_studyzone.matrix[!rownames(obs_studyzone.matrix) < last.first.year, ]
  }
  return(obs_studyzone.matrix)
}

####################################@

MapVisualise_location <- function(Cluster.Info,station.id, add.medoid=TRUE,add.others=TRUE){
  medoid.Info <- StudyZone_cluster_kid.Information[1,]
  color.kid <- factor(Cluster.Info$colorCluster)[1]
  color.id <- get_color_name(color.kid)
  id.kid <- medoid.Info$cluster_id
  Highligth_station.Info <- dplyr::filter(Cluster.Info, indexsites_inCluster  %in% station.id)
  if (add.medoid==TRUE & add.others==TRUE){
    france_map <- map_data("france")
    mapcluster <- ggplot() +
      geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey60', size = 1) +
      coord_quickmap() +
      geom_point(data = Cluster.Info, aes(x = longitude..., y = latitude...), fill="grey", color = "grey2", pch = 21, size = 1) +
      geom_point(data = Highligth_station.Info , aes(x = longitude..., y = latitude..., ),fill=color.id,color = "grey2", pch = 21, size = 1.5)+
      geom_point(data =medoid.Info, aes(x = longitude..., y = latitude...), size = 2, shape = 8)+
      theme(legend.position="none") +
      labs(title = paste("Study Zone",id.kid ))
  }
  if (add.medoid==TRUE & add.others==FALSE){
    france_map <- map_data("france")
    mapcluster <- ggplot() +
      geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey60', size = 1) +
      coord_quickmap() +
      geom_point(data = Highligth_station.Info , aes(x = longitude..., y = latitude...),fill=color.id,color = "grey2", pch = 21, size = 1.5)+
      geom_point(data =medoid.Info, aes(x = longitude..., y = latitude...), size = 2, shape = 8)+
      theme(legend.position="none") +
      labs(title = paste("Study Zone",id.kid ))
  }
  
  if (add.medoid==FALSE & add.others==TRUE){
    france_map <- map_data("france")
    mapcluster <- ggplot() +
      geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey60', size = 1) +
      coord_quickmap() +
      geom_point(data = Cluster.Info, aes(x = longitude..., y = latitude...), fill="grey", color = "grey2", pch = 21, size = 1) +
      geom_point(data = Highligth_station.Info , aes(x = longitude..., y = latitude...),fill=color.id,color = "grey2", pch = 21, size = 1.5)+
      theme(legend.position="none") +
      labs(title = paste("Study Zone",id.kid ))
  }
  
  if (add.medoid==FALSE & add.others==FALSE){
    france_map <- map_data("france")
    mapcluster <- ggplot() +
      geom_polygon(data = france_map, aes(x = long, y = lat, group = group), fill = NA, col = 'grey60', size = 1) +
      coord_quickmap() +
      geom_point(data = Highligth_station.Info , aes(x = longitude..., y = latitude...),fill=color.id,color = "grey2", pch = 21, size = 1.5)+
      theme(legend.position="none") +
      labs(title = paste("Study Zone",id.kid ))
  }
  print(mapcluster)
}

PlotData.separate.VersionSQRTTMAX <- function(Data,Studyzone.Info, siteid = 0, par1 = 4, par2 = 2 , add.max=FALSE , max.type="all" ){ # if  siteid = 0 we plot data in all sites 
  par(mfrow=c(par1,par2))
  timesteps <- dim(Data)[1]
  nbsites <- dim(Data)[2]
  indexsites_zone <- as.numeric(Studyzone.Info$indexsites_inCluster)
  startYears.vec <- as.numeric(Studyzone.Info$year.start)
  maxyear <- max(startYears.vec)
  plotaxis.vec <- c(1,seq(10,40,10),timesteps)
  if (length(siteid)==1){
    if (siteid>0){
      indexsitesplot <- indexsites_zone[siteid]
      par1 <- 1
      par2 <- 1 
      plot(c(1:timesteps),Data[,siteid],ylim=c(min(Data),max(Data)),xaxt = "n",
           xlab="year",ylab='T°', main=bquote("T° at site "~ .(indexsitesplot)))
      if(add.max==TRUE){
        max.traj <- apply(X=Data, MARGIN=1, FUN=max)
        points(c(1:timesteps),max.traj,col='red',pch=4)
        axis(1, at= plotaxis.vec, labels= plotaxis.vec + (maxyear-1),cex.axis=.5)
        legend('bottomright',legend="max all sites",col="red",pch=4, bg="transparent")
      }
      #text(c(1:timesteps)[10],max(Data[1:(timesteps/3)])-1, cex=1.3,col="blue",paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
    }
    
    if (siteid==0){
      ylim.max = max(Data)
      ylim.min = min(Data)
      if(add.max==FALSE){
        for(i in 1:nbsites){
          indexsitesplot <- indexsites_zone[i]
          plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="year",ylab='T°',xaxt = "n",
               main=bquote("T° at site "~ .(indexsitesplot)))
          axis(1, at= plotaxis.vec, labels= plotaxis.vec + (maxyear-1),cex.axis=.5)
          #text(c(1:timesteps)[10],max(Data[1:(timesteps/3)])-1,cex=1.3,col="blue", paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
        }
      }
      
      if(add.max==TRUE){
        max.traj <- apply(X=Data, MARGIN=1, FUN=max)
        for(i in 1:nbsites){
          indexsitesplot <- indexsites_zone[i]
          plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="year",ylab='T°',xaxt = "n",
               main=bquote("T° at site "~ .(indexsitesplot)))
          points(c(1:timesteps),max.traj,col='red',pch=4)
          axis(1, at= plotaxis.vec, labels= plotaxis.vec + (maxyear-1),cex.axis=.5)
          #text(c(1:timesteps)[10],max(Data[1:(timesteps/3)])-1,cex=1.3,col="blue", paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
          legend('bottomright',legend="max all sites",col="red",pch=4, bg="transparent")
        }
      }
    }
  }
  
  if (length(siteid)>1){
    nb.plots <- length(siteid)
    ylim.max = max(Data)
    ylim.min = min(Data)
    
    if(add.max==FALSE){
      for(i in siteid ){
        indexsitesplot <- indexsites_zone[i]
        plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="year",ylab='T°',xaxt = "n",
             main=bquote("T° at site "~ .(indexsitesplot)))
        axis(1, at= plotaxis.vec, labels= plotaxis.vec + (maxyear-1),cex.axis=.5)
        #text(c(1:timesteps)[10],max(Data[1:(timesteps/3)])-1, cex=1.3,col='blue',paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
      }
    }
    
    if(add.max==TRUE){
      if(max.type=="all"){
        max.traj <- apply(X=Data, MARGIN=1, FUN=max)
        colmax <- "red"
        legend.text <- "max all sites" }
      
      if(max.type=="not.all"){
        max.traj <- apply(X=Data[,c(siteid)], MARGIN=1, FUN=max)
        colmax <- "green"
        legend.text <- "max showed sites"}
      
      for(i in siteid ){
        indexsitesplot <- indexsites_zone[i]
        plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="year",ylab='T°',xaxt = "n",
             main=bquote("T° at site "~ .(indexsitesplot)))
        points(c(1:timesteps),max.traj,col=colmax,pch=4)
        axis(1, at= plotaxis.vec, labels= plotaxis.vec + (maxyear-1),cex.axis=.5)
        #text(c(1:timesteps)[10],max(Data[1:(timesteps/3)])-1, cex=1.3,col='blue',paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
        legend('bottomright',legend=legend.text,col=colmax,pch=4, bg="transparent")
      }
    }
    
  }
}

PlotData.compare.VersionSQRTTMAX <- function(Data, Studyzone.Info, site.interest = 1, site.others = 0,  max.plot = FALSE){
  nbtimesteps <- dim(Data)[1]
  nbsites <- dim(Data)[2]
  indexsites_zone <- as.numeric(Studyzone.Info$indexsites_inCluster)
  startYears.vec <- as.numeric(Studyzone.Info$year.start)
  maxyear <- max(startYears.vec)
  plotaxis.vec <- c(1,seq(10,40,10),nbtimesteps)
  if(length(site.interest) != 1){
    stop("We must observe only 1 site of interest")}
  
  nbtimesteps <- dim(Data)[1]
  timesteps <- c(1:nbtimesteps)
  nbsites <- dim(Data)[2]
  plotaxis.vec <- c(1,seq(10,40,10),nbtimesteps)
  
  site.id <- site.interest
  others.id <- site.others
  col.points.max <- "green"
  if(length(site.others) == 1){
    if(site.others == 0){
      others.id <- c(1:nbsites)
      col.points.max <- "red"}
  }
  
  if(max.plot == FALSE){par(mfrow=c(1,1), mar = c(5, 5.3, 4, 2))}
  if(max.plot == TRUE){par(mfrow=c(1,2), mar = c(5, 5.3, 4, 2))}
  
  ylim.max = max(Data)
  ylim.min = min(Data)
  
  sitesother.names <- paste(others.id,collapse=" ")
  
  indexsitesplot <- indexsites_zone[site.interest]
  plot(timesteps,Data[,site.id],ylim=c(ylim.min,ylim.max),xlab="year",ylab='T°',xaxt = "n",
       main=bquote(atop("T ° at station "~ .(indexsitesplot),"w.r.t to Zone ")))
  for(i in others.id ){points(timesteps,Data[,i],col='gray')}
  points(timesteps,Data[,site.id],col="black")
  axis(1, at= plotaxis.vec, labels= plotaxis.vec + (maxyear-1),cex.axis=.5)
  # text(c(1:nbtimesteps)[13],max(Data)-2, cex=1.3,col='blue',paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
  legend('bottomright',legend=c("station of study","other stations"),col=c("black",'gray'),pch=1, bg="transparent")
  
  if(max.plot == TRUE){
    indexsitesplot <- indexsites_zone[site.interest]
    plot(timesteps,Data[,site.id],ylim=c(ylim.min,ylim.max),xlab="year",ylab='T°',xaxt = "n",
         main=bquote(atop("T° at stationn "~ .(indexsitesplot),"w.r.t to Zone")))
    for(i in others.id ){points(timesteps,Data[,i],col='gray')}
    points(timesteps,Data[,site.id],col="black")
    max.traj <- apply(X=Data[,c(site.id,others.id)], MARGIN=1, FUN=max)
    points(timesteps,max.traj,pch=4,col=col.points.max)
    axis(1, at= plotaxis.vec, labels= plotaxis.vec + (maxyear-1),cex.axis=.5)
    # text(c(1:nbtimesteps)[13],max(Data)-2, cex=1.3,col='blue',paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
    legend('bottomright',legend=c("station of study","other stations","max stations"),col=c("black",'gray',col.points.max),pch=c(1,1,4), bg="transparent")
  }
}

find_records_func.VersionSQRTTMAX  <- function(Data,Studyzone.Info,marginal=FALSE, same.site = TRUE , other.site = 1){
  nbsites <- dim(Data)[2]
  nbtimesteps <- dim(Data)[1]
  timesteps <- 1:nbtimesteps
  yearsnames <- as.numeric(rownames(Data))
  indexstations <- as.numeric(Studyzone.Info$indexsites_inCluster)
  record.value <- c()
  record.time <- c()
  record.site <- c()
  
  if (marginal==FALSE){
    marginalRecords.mat <- apply(X=t(Data), MARGIN=1, FUN=cummax)
    TotalRecords.vec <- apply(X=t(marginalRecords.mat), MARGIN=2, FUN=max)
    for (s in 1:nbsites){
      ok <- TotalRecords.vec == Data[,s]
      ok[1] <- FALSE
      if (any(ok)==TRUE){
        record.s.value.vec <- Data[ok,s]
        record.s.time.vec <- timesteps[ok]
        record.value <- append(record.value,record.s.value.vec)
        record.time <- append(record.time,record.s.time.vec)
        record.site <- append(record.site,rep(s,length(record.s.time.vec)))
      }
    }
    record.time <- yearsnames[as.numeric(record.time)]
    record.site <- indexstations[as.numeric(record.site)]
  }
  
  if (marginal==TRUE){
    marginalRecords.mat <- apply(X=t(Data), MARGIN=1, FUN=cummax)
    if (same.site == TRUE){
      for (s in 1:nbsites){
        ok <- marginalRecords.mat[,s] == Data[,s]
        ok[1] <- FALSE
        if (any(ok)==TRUE){
          record.s.value.vec <- Data[ok,s]
          record.s.time.vec <- timesteps[ok]
          record.value <- append(record.value,record.s.value.vec)
          record.time <- append(record.time,record.s.time.vec)
          record.site <- append(record.site,rep(s,length(record.s.time.vec)))
        }
      }
      record.time <- yearsnames[as.numeric(record.time)]
      record.site <- indexstations[as.numeric(record.site)]
    }
    
    
    if (same.site == FALSE){
      for (s in 1:nbsites){
        ok <- marginalRecords.mat[,other.site] <= Data[,s]
        ok[1] <- FALSE
        if (any(ok)==TRUE){
          record.s.value.vec <- Data[ok,s]
          record.s.time.vec <- timesteps[ok]
          record.value <- append(record.value,record.s.value.vec)
          record.time <- append(record.time,record.s.time.vec)
          record.site <- append(record.site,rep(s,length(record.s.time.vec)))
        }
      }
      record.time <- yearsnames[as.numeric(record.time)]
      record.site <- indexstations[as.numeric(record.site)]
    }
    
    
  }
  
  record.table <- cbind(record.time,record.value,record.site)
  
  sorted_record.table <- record.table[order(record.table[, 1]), ]
  return(sorted_record.table)
}

PlotObservedRecord.separate.VersionSQRTTMAX <- function(Data, Studyzone.Info, siteid = 0, par1 = 4, par2 = 2 ,add.max=FALSE , max.type="all", add.record.total = FALSE, add.record.marginal = FALSE,add.line=FALSE ){ # if  siteid = 0 we plot data in all sites 
  par(mfrow=c(par1,par2))
  timesteps <- dim(Data)[1]
  nbsites <- dim(Data)[2]
  indexsites_zone <- as.numeric(Studyzone.Info$indexsites_inCluster)
  startYears.vec <- as.numeric(Studyzone.Info$year.start)
  maxyear <- max(startYears.vec)
  plotaxis.vec <- c(1,seq(10,40,10),timesteps)
  # If we chose to plot observations in one single site
  if (length(siteid)==1){
    if (siteid>0){
      par1 <- 1
      par2 <- 1 
      indexsitesplot <- indexsites_zone[siteid]
      plot(c(1:timesteps),Data[,siteid],ylim=c(min(Data),max(Data)),xaxt = "n",
           xlab="year",ylab='T°', main=bquote("Records at site "~ .(indexsitesplot)))
      axis(1, at= plotaxis.vec, labels= plotaxis.vec + (maxyear-1),cex.axis=.5)
      if(add.record.marginal == TRUE){
        MarginalRecord.matrix <- find_records_func(Data,marginal=TRUE)
        MarginalRecord.siteid <- MarginalRecord.matrix[MarginalRecord.matrix[,3] == siteid,]
        if(add.line==TRUE){
          abline(v=MarginalRecord.siteid[,1], col="gray",lty=2)
        }
        points(MarginalRecord.siteid[,1],MarginalRecord.siteid[,2],col="black", bg="wheat2",pch=21)
        
      }
      if(add.record.total == TRUE){
        TotalRecord.matrix <- find_records_func(Data)
        TotalRecord.siteid <- TotalRecord.matrix[TotalRecord.matrix[,3] == siteid,]
        if(add.line==TRUE){
          abline(v=TotalRecord.siteid[,1], col="gray",lty=2)
        }
        points(TotalRecord.siteid[,1],TotalRecord.siteid[,2],col="black", bg="blue",pch=21)
        
      }
      if(add.max==TRUE){
        max.traj <- apply(X=Data, MARGIN=1, FUN=max)
        points(c(1:timesteps),max.traj,col='red',pch=4)
        if(add.record.total==TRUE){
          points(TotalRecord.siteid[,1],TotalRecord.siteid[,2],col="blue",pch=4)
        }
      }
      if(add.max==TRUE & add.record.total==TRUE & add.record.marginal == TRUE){
        legend('bottomright',legend=c("max all sites", "record all sites", "marginal record"),col=c("red","blue","gray"),pch=c(4,4,16), bg="transparent")
      }
      if(add.max==TRUE & add.record.total==TRUE & add.record.marginal == FALSE){
        legend('bottomright',legend=c("max all sites", "record all sites"),col=c("red","blue"),pch=c(4,4), bg="transparent")
      }
      if(add.max==TRUE & add.record.total==FALSE & add.record.marginal == TRUE){
        legend('bottomright',legend=c("max all sites","marginal record"),col=c("red","wheat2"),pch=c(4,16), bg="transparent")
      }
      if(add.max==TRUE & add.record.total==FALSE & add.record.marginal == FALSE){
        legend('bottomright',legend="max all sites",col="red",pch=4, bg="transparent")
      }
    }
    
    if (siteid==0){
      ylim.max = max(Data)
      ylim.min = min(Data)
      
      if(add.record.marginal == TRUE){
        MarginalRecord.matrix <- find_records_func(Data,marginal=TRUE)
      }
      if(add.record.total == TRUE){
        TotalRecord.matrix <- find_records_func(Data)
      }
      
      for(i in 1:nbsites){
        indexsitesplot <- indexsites_zone[i]
        plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="year",ylab='T°',xaxt = "n",
             main=bquote("Records at site "~ .(indexsitesplot)))
        axis(1, at= plotaxis.vec, labels= plotaxis.vec + (maxyear-1),cex.axis=.5)
        if(add.record.marginal == TRUE){
          MarginalRecord.i <- MarginalRecord.matrix[MarginalRecord.matrix[,3] == i,]
          if(add.line==TRUE){
            abline(v=MarginalRecord.i[,1], col="gray",lty=2)
          }
          points(MarginalRecord.i[,1],MarginalRecord.i[,2],col="black", bg="wheat2",pch=21)
          
        }
        if(add.record.total == TRUE){
          if (i %in%  TotalRecord.matrix[,3] == TRUE ){
            TotalRecord.i <- TotalRecord.matrix[TotalRecord.matrix[,3] == i,]
            TotalRecord.i <- as.matrix(TotalRecord.i)
            if(dim(TotalRecord.i)[2]==1){TotalRecord.i <- t(TotalRecord.i)}
            if(add.line==TRUE){
              abline(v=TotalRecord.i[,1], col="gray",lty=2)
            }
            points(as.numeric(TotalRecord.i[,1]),as.numeric(TotalRecord.i[,2]),col="black", bg="blue",pch=21)
            
            
          }
          
        }
        
        if(add.max==TRUE){
          max.traj <- apply(X=Data, MARGIN=1, FUN=max)
          points(c(1:timesteps),max.traj,col='red',pch=4)
          if(add.record.total==TRUE){
            if (i %in%  TotalRecord.matrix[,3] == TRUE ){
              TotalRecord.i <- TotalRecord.matrix[TotalRecord.matrix[,3] == i,]
              TotalRecord.i <- as.matrix(TotalRecord.i)
              if(dim(TotalRecord.i)[2]==1){TotalRecord.i <- t(TotalRecord.i)}
              points(as.numeric(TotalRecord.i[,1]),as.numeric(TotalRecord.i[,2]),col="blue",pch=4)
            }
          }
        }
        if(add.max==TRUE & add.record.total==TRUE & add.record.marginal == TRUE){
          legend('bottomright',legend=c("max all sites", "record all sites", "marginal record"),col=c("red","blue","wheat2"),pch=c(4,4,16), bg="transparent")
        }
        if(add.max==TRUE & add.record.total==TRUE & add.record.marginal == FALSE){
          legend('bottomright',legend=c("max all sites", "record all sites"),col=c("red","blue"),pch=c(4,4), bg="transparent")
        }
        if(add.max==TRUE & add.record.total==FALSE & add.record.marginal == TRUE){
          legend('bottomright',legend=c("max all sites","marginal record"),col=c("red","wheat2"),pch=c(4,16), bg="transparent")
        }
        if(add.max==TRUE & add.record.total==FALSE & add.record.marginal == FALSE){
          legend('bottomright',legend="max all sites",col="red",pch=4, bg="transparent")
        }
      }
    }
  }
  
  
  if (length(siteid)>1){
    
    nb.plots <- length(siteid)
    ylim.max = max(Data)
    ylim.min = min(Data)
    
    # plotted record will also depend on "max.type"
    if(add.record.marginal == TRUE){
      MarginalRecord.matrix <- find_records_func(Data,marginal=TRUE)
    }
    if(add.record.total == TRUE){
      if(max.type=="all"){
        TotalRecord.matrix <- find_records_func(Data)
      }
      if(max.type=="not.all"){
        TotalRecord.matrix <- find_records_func(Data[,siteid])
        nbcol <- 1:dim(Data[,siteid])[2]
      }
    }
    
    colrecord <- "blue"
    
    
    if(add.max==TRUE){
      if(max.type=="all"){
        max.traj <- apply(X=Data, MARGIN=1, FUN=max)
        colmax <- "red"
        colrecord <- "blue"
        legendmax.text <- "max all sites" 
        legendrecord.text <- "record all sites"}
      
      if(max.type=="not.all"){
        max.traj <- apply(X=Data[,c(siteid)], MARGIN=1, FUN=max)
        colmax <- "green"
        colrecord <- "cyan"
        legendmax.text <- "max showed sites"
        legendrecord.text <- "record showed sites"}
      
      j <- 1
      for(i in siteid ){
        indexsitesplot <- indexsites_zone[i]
        plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="year",ylab='T°',xaxt = "n",
             main=bquote("Records at site "~ .(indexsitesplot)))
        axis(1, at= plotaxis.vec, labels= plotaxis.vec + (maxyear-1),cex.axis=.5)
        if(add.record.marginal == TRUE){
          MarginalRecord.i <- MarginalRecord.matrix[MarginalRecord.matrix[,3] == i,]
          if(add.line==TRUE){
            abline(v=MarginalRecord.i[,1], col="gray",lty=2)
          }
          points(MarginalRecord.i[,1],MarginalRecord.i[,2],col="black", bg="gray",pch=21)
          
        }
        if(add.record.total == TRUE){
          TotalRecord.i <- TotalRecord.matrix[TotalRecord.matrix[,3] == j,]
          if(add.line==TRUE){
            abline(v=TotalRecord.i[,1], col="gray",lty=2)
          }
          points(TotalRecord.i[,1],TotalRecord.i[,2],col="black", bg=colrecord,pch=21)
          
        }
        if(add.max==TRUE){
          points(c(1:timesteps),max.traj,col=colmax,pch=4)
          if(add.record.total==TRUE){
            points(TotalRecord.i[,1],TotalRecord.i[,2],col=colrecord,pch=4)
          }
        }
        
        if(add.max==TRUE & add.record.total==TRUE & add.record.marginal == TRUE){
          legend('bottomright',legend=c(legendmax.text, legendrecord.text, "marginal record"),col=c("red",colrecord,"gray"),pch=c(4,4,16), bg="transparent")
        }
        if(add.max==TRUE & add.record.total==TRUE & add.record.marginal == FALSE){
          legend('bottomright',legend=c(legendmax.text, legendrecord.text),col=c(colmax,colrecord),pch=c(4,4), bg="transparent")
        }
        if(add.max==TRUE & add.record.total==FALSE & add.record.marginal == TRUE){
          legend('bottomright',legend=c(legendmax.text,"marginal record"),col=c(colmax,"gray"),pch=c(4,16), bg="transparent")
        }
        if(add.max==TRUE & add.record.total==FALSE & add.record.marginal == FALSE){
          legend('bottomright',legend=legendmax.text,col=colmax,pch=4, bg="transparent")
        } 
        
        if(add.record.total==TRUE & add.record.marginal == TRUE){
          legend('bottomright',legend=c("total record","marginal record"),col=c("blue",'gray'),pch=16, bg="transparent")
        } 
        
        j <- j+1 
      }
      
    }
    
  }
}


##############################

cdf_exponential <- function(x) pexp(x, rate = 1) 

compute_CRPS <- function(x, cdf_exponential) {
  n <- length(x)
  crps <- 0
  for (i in 1:n) {
    crps <- crps + (cdf_exponential(x[i]) - sum(x <= x[i]) / n)^2
  }
  crps <- crps / n
  return(crps)
}

########################

Plot.ProbMargs_so.VersionSQRTTMAX <- function(E.GmaxY,Studyzone.Info,site.interest=1,site.others=0,par1=ppar1,par2=ppar2,log.scale=FALSE,add.iid=FALSE,add.iid.scale=FALSE){
  library(latex2exp)
  timesteps <- dim(E.GmaxY)[1]
  startYears.vec <- as.numeric(Studyzone.Info$year.start)
  maxyear <- max(startYears.vec)
  plotaxis.vec <- c(1,seq(10,40,10),timesteps)
  sitesnames <- as.numeric(Studyzone.Info$indexsites_inCluster)
  E.GmaxY <- as.matrix(E.GmaxY[,,site.interest],nrow=dim(E.GmaxY)[1],ncol=dim(E.GmaxY)[2])
  par(mfrow=c(par1,par2), mar = c(5, 5.3, 4, 2))
  if(site.others==0){
    site.others <- 1:dim(E.GmaxY)[2]
  }
  
  if(log.scale==FALSE){
    if(add.iid == FALSE){
      ylim.max.p <- max(E.GmaxY)
      for(j in site.others){
        plot(c(1:dim(E.GmaxY)[1]),E.GmaxY[,j], type='l',ylim=c(0,ylim.max.p),col="red",xaxt = "n",
             main=bquote(atop("Record probability "~ .(sitesnames[site.interest])," w.r.t site " ~ .(sitesnames[j]))),
             ylab=TeX(r'($\hat{p}_t^{s_o,s}$)'), xlab="year")
        axis(1, at= plotaxis.vec, labels= plotaxis.vec + (maxyear-1),cex.axis=.5)
      }
    }
    if(add.iid == TRUE){
      p.iid <- 1/(c(2:(dim(E.GmaxY)[1]+1)))
      ylim.max.p <- max(E.GmaxY)
      if(add.iid.scale == TRUE){ylim.max.p <- max(E.GmaxY,p.iid )}
      for(j in site.others){
        plot(c(1:dim(E.GmaxY)[1]),p.iid , type='l',ylim=c(0,ylim.max.p),col="black",xaxt = "n",
             main=bquote(atop("Record probability "~ .(sitesnames[site.interest])," w.r.t site " ~ .(sitesnames[j]))),
             ylab=TeX(r'($\hat{p}_t^{s_o,s}$)'), xlab="year")
        axis(1, at= plotaxis.vec, labels= plotaxis.vec + (maxyear),cex.axis=.5)
        lines(c(1:dim(E.GmaxY)[1]),E.GmaxY[,j], type='l',col="red")
      }
    }
  }
  
  if(log.scale==TRUE){
    plotaxis.vec <- c(2,seq(10,40,10),timesteps)
    plotaxisLog.vec <- log(plotaxis.vec)
    E.GmaxY <- log(E.GmaxY)
    if(add.iid == FALSE){
      ylim.max.p <- min(E.GmaxY)
      for(j in site.others){
        plot(log(c(2:(dim(E.GmaxY)[1]+1))),E.GmaxY[,j], type='l',ylim=c(ylim.max.p,0),col="red",xaxt = "n",
             main=bquote(atop("Record probability "~ .(sitesnames[site.interest])," w.r.t site " ~ .(sitesnames[j]))),
             ylab=TeX(r'($\og(\hat{p}_t^{s_o,s})$)'), xlab=" ")
        axis(2, at= plotaxisLog.vec, labels= (plotaxis.vec) + (maxyear),cex.axis=.5)
      }
    }
    if(add.iid == TRUE){
      p.iid <- log(1/(c(2:(dim(E.GmaxY)[1]+1))))
      ylim.max.p <- min(E.GmaxY)
      if(add.iid.scale == TRUE){ylim.max.p <- min(E.GmaxY,p.iid )}
      for(j in site.others){
        plot(log(c(2:(dim(E.GmaxY)[1]+1))),p.iid, type='l',ylim=c(ylim.max.p,0),col="black",xaxt = "n",
             main=bquote(atop("Record probability "~ .(sitesnames[site.interest])," w.r.t site " ~ .(sitesnames[j]))),
             ylab=TeX(r'($\log(\hat{p}_t^{s_o,s})$)'), xlab=" ")
        axis(1, at= plotaxisLog.vec, labels= (plotaxis.vec) + (maxyear),cex.axis=.5)
        lines(log(c(2:(dim(E.GmaxY)[1]+1))),E.GmaxY[,j], type='l',col="red")
      }
    }
  }
  
  
  
}

Plot.QQplotExp1.ProbMargs_so.VersionSQRTTMAX <- function(GmaxY.MatAarray,E.GmaxY.MatAarray,sitesnames,site.interest=1,site.others=0,par1=ppar1,par2=ppar2, add.print=FALSE){
  library(latex2exp)
  #sitesnames <- as.numeric(Studyzone.Info$indexsites_inCluster)
  GmaxY.MatAarray <- as.matrix(GmaxY.MatAarray[,,site.interest],nrow=dim(GmaxY.MatAarray)[1],ncol=dim(GmaxtY.MatAarray)[2])
  E.GmaxY.MatAarray <- as.matrix(E.GmaxY.MatAarray[,,site.interest],nrow=dim(E.GmaxY.MatAarray)[1],ncol=dim(E.GmaxY.MatAarray)[2])
  if(site.others[1]==0){site.others <- 1:dim(E.GmaxY.MatAarray)[2]}
  if (add.print==TRUE){
    countZeros.so <- apply(GmaxY.MatAarray, 2, function(c)sum(c==0))
    print(rbind(sitesnames,countZeros.so))}
  # when we use "matGZ_func_bis(.)" we won't need to do that -> do that later
  MinValReplace <- min(GmaxY.MatAarray[GmaxY.MatAarray>0])
  GmaxY.MatAarray[GmaxY.MatAarray==0] <- MinValReplace
  ExpoDist.MatAarray <- -log(GmaxY.MatAarray)
  lambda.MatAarray <- (1-E.GmaxY.MatAarray)/E.GmaxY.MatAarray
  ExpUnit.MatAarray <-  ExpoDist.MatAarray * (lambda.MatAarray)^(-1)
  
  par(mfrow=c(par1,par2), mar = c(5, 5, 5, 2))
  
  for(j in site.others){
    ExpUnit.ij <- as.numeric(ExpUnit.MatAarray[,j])
    theoretical_quantiles <- qexp(ppoints(length(ExpUnit.ij)))
    qqplot(theoretical_quantiles,ExpUnit.ij, main = bquote(atop("Exponential QQ-plot",~ .(sitesnames[site.interest])~" w.r.t " ~ .(sitesnames[j]))),
           xlab = "Theoretical Exp(1) Quantiles", ylab="Empirical Exp(1) Quantiles",xlim=c(0,max(ExpUnit.MatAarray )),ylim=c(0,max(ExpUnit.MatAarray)),
           cex.main=1.1,cex.lab=1.1)
    abline(0, 1, col = "red")
  }
}

Get_CRPS_soVSzone <- function(Data,Studyzone.Info,site.interest=1,site.others=0,h1.vec,h2.vec,add.print=FALSE){
  sitesnames <- as.numeric(Studyzone.Info$indexsites_inCluster)
  siteid <- site.interest
  if(site.others==0){site.others <- 1:dim(Data)[2]}
  if (length(h1.vec)==1){h1.vec <- rep(h1.vec,dim(Data)[2])}
  if (length(h2.vec)==1){h2.vec <- rep(h2.vec,dim(Data)[2])}
  CRPS.so.zone <- rep(NA,length(site.others))
  for (j in site.others){
    YiYj_mat <- as.matrix(na.omit(Data[,c(siteid,j)]))
    Gmaxj.Yso <- GmaxY.estimation.function.step1(YiYj_mat,site.interest = 1,site.others = 2,bandwidth= h1.vec[j])$Array.GmaxHatY
    Exp.transf_Gmaxj.so <- -log(as.numeric(Gmaxj.Yso))
    P.marginalj.so <- as.numeric(EGmaxY.estimation.function.step2(Gmaxj.Yso, bandwidth= h2.vec[j])$Array.EGmaxHatY)
    Lambda.marginalj.so <- (1-P.marginalj.so )/P.marginalj.so 
    ExpUnit_j.so <- Exp.transf_Gmaxj.so * ( Lambda.marginalj.so)^(-1)
    cdf_exponential <- function(x) pexp(x, rate = 1) 
    crps_j.so <- compute_CRPS(ExpUnit_j.so, cdf_exponential)
    if(add.print==TRUE){cat(paste("station n°",sitesnames[j]," has CRPS =  ", round(crps_j.so,3),"\n"))}
    CRPS.so.zone[j] <- crps_j.so
  }
  names(CRPS.so.zone) <- sitesnames
  return(CRPS.so.zone)
}  

Get_CRPS_HoTrue <- function(E.GmaxY.MatAarray,site.interest=1,site.others=0,Nn = 2000){
  siteid <- site.interest
  if(site.others==0){site.others <- 1:dim(E.GmaxY.MatAarray)[2]}
  CRPS.Nn_HoTrue <- matrix(NA,ncol=length(site.others),nrow=Nn)
  lambda.MatAarray <- (1-E.GmaxY.MatAarray)/E.GmaxY.MatAarray
  lambdaMarg_so.mat<- as.matrix(lambda.MatAarray[,site.others,siteid],ncol=length(site.others),nrow=dim(E.GmaxY.MatAarray)[1])
  for(n in 1:Nn){
    ExpDistr.ijt <- apply(lambdaMarg_so.mat, 1:2, function(x) rexp(1,rate=(x)^(-1)))
    ExpUnit.ijt  <- ExpDistr.ijt * (lambdaMarg_so.mat)^(-1) 
    crpsHo_so.j <- apply(ExpUnit.ijt, 2, function(x) compute_CRPS(x,cdf_exponential))
    CRPS.Nn_HoTrue[n,] <- crpsHo_so.j
  }
  return(CRPS.Nn_HoTrue)
}

Eval_CRPS_HoTrue <- function(crpsHat.vec,crpsHo.mat,type.range=TRUE,add.plot = FALSE,par1,par2){
  sitesnames <- names(crpsHat.vec)
  nbsites <- length(crpsHat.vec)
  is.ok.vec <- rep(NA,nbsites)
  if(is.null(sitesnames)==TRUE){sitesnames <- 1:nbsites }
  if (type.range==FALSE){
    crps_thershold.vec <- apply(crpsHo.mat, 2, function(x) max(x))
  }
  if (type.range==TRUE){
    crps_thershold.vec <- apply(crpsHo.mat, 2, function(x) quantile(x, 0.975))  
    
    # crps_thershold_low.vec <- apply(crpsHat.vec, 2, function(x) quantile(x, 0.025))
    # within_range <- sum(crp_ho_vec >= quantile(crp_ho_vec, 0.025) & crp_ho_vec <= quantile(crp_ho_vec, 0.975)) / length(crp_ho_vec)
    # is_within_95_range <- crps_hat >= quantile(crp_ho_vec , 0.025) & crps_hat <= quantile(crp_ho_vec, 0.975)
    # cat("Proportion of values within the 95% range:", within_range, "\n")
    # cat("\nIs n_o within the 95% range?", is_within_95_range, "\n")
    
  }
  
  for(j in 1:nbsites){
    is_ok <-  crps_thershold.vec[j] >= crpsHat.vec[j]
    is.ok.vec[j] <- is_ok
  }
  
  if (add.plot==TRUE){
    par(mfrow=c(par1,par2), mar = c(4, 4, 2, 2))
    for(j in 1:nbsites){
      sitej <- sitesnames[j]
      hist(crpsHo.mat[,j], main = paste("Histogram station" , sitej),
           xlim=c(min(min(crpsHo.mat[,j]),crpsHat.vec[j]),max(max(crpsHo.mat[,j]),crpsHat.vec[j])),
           xlab='CRPS when Ho TRUE')
      abline(v=crpsHat.vec[j],col='red', lwd = 2, lty = "dashed")
      abline(v=crps_thershold.vec[j],col="blue", lwd = 2, lty = "dashed")
      text(x=(.7*crpsHat.vec[j]), y=300, col="red",paste0("crps =", round(crpsHat.vec[j], 3)))
      text(x=(.5*crps_thershold.vec[j]), y=400, col="blue", paste0("thershold=", round(crps_thershold.vec[j], 3)))
    }
  }
  names(is.ok.vec) <- sitesnames 
  return(is.ok.vec)
}

Visu_stations_Ho <- function(crpsHo.truefalse,Studyzone.Info,Data,GmaxY.MatAarray,E.GmaxY.MatAarray, Ho = FALSE, add.qqplot = TRUE,add.map = FALSE, add.obs= FALSE,par1,par2){
  if(Ho==TRUE){sitesVisu <- as.numeric(which(crpsHo.truefalse == TRUE))}
  if(Ho==FALSE){sitesVisu <- as.numeric(which(crpsHo.truefalse == FALSE))}
  
  if(add.qqplot == TRUE){
    Plot.QQplotExp1.ProbMargs_so.VersionSQRTTMAX(GmaxY.MatAarray,E.GmaxY.MatAarray,Studyzone.Info,site.interest=1,site.others=sitesVisu,par1,par2, add.print=FALSE)
  }
  
  if(add.obs == TRUE){
    par(mfrow=c(par1,par2), mar = c(4, 4, 2, 2))
    for (sitej in sitesVisu){
      PlotData.compare.VersionSQRTTMAX(Data,Studyzone.Info, site.interest = 1, site.others = sitej, max.plot = FALSE)
    }
  }
  
  if(add.map == TRUE){
    stationidHo <- as.numeric(Studyzone.Info$indexsites_inCluster)[sitesVisu]
    MapVisualise_location(Cluster.Info = Studyzone.Info
                          ,station.id = stationidHo , add.medoid=TRUE,add.others=TRUE)
  }
  
}

####################################################################################

Get.P.spatiotemp <- function(Data,EGmaxY,Studyzone.Info,site.interest=1,add.iid=TRUE,add.indep=TRUE,add.stat=TRUE,h1,par1=ppar1,par2=ppar2,use.timesteps = 0){
  
  siteid <- site.interest
  Pmarg_so_j <- as.matrix(EGmaxY[,,siteid])
  lambda_so_j <- (1 - Pmarg_so_j)/Pmarg_so_j
  FjYj <- GiYi.func(Data = Data, h = h1)
  lambdaPickands <- Version2_V.applied.lambda.estim.func(GiYi.mat = FjYj,lambda = lambda_so_j, use.timesteps = use.timesteps)
  lambda_so <- lambdaPickands$lambdaV
  Pickands_so <- lambdaPickands$A
  Precord_so <- 1/(1+lambda_so)
  nbtimesteps <- length(Precord_so)
  
  returnObject_list <- list(P.record = Precord_so , A = Pickands_so)
  
  if(add.iid==TRUE){
    tt <- c(1:nbtimesteps) +1
    Piid_so <- 1/(1+((tt-1)*(dim(FjYj)[2])))
    iidObject_list <- list(P.record.iid = Piid_so)
    returnObject_list <- c(returnObject_list,iidObject_list)
  }
  if(add.indep==TRUE){ 
    tt <- c(1:nbtimesteps)
    zi_so_j  <- 1/(lambda_so_j *tt)
    sum.zi.ti <- rowSums(zi_so_j^(-1))
    Pindep_so <- 1/(1+sum.zi.ti)
    indepObject_list <- list(P.record.indep = Pindep_so)
    returnObject_list <- c(returnObject_list,indepObject_list)
  }
  
  if((add.stat==TRUE)){
    bandwithInf <- nbtimesteps*10
    GmaxHatY.stat <- GmaxY.estimation.function.step1(Data, bandwidth = bandwithInf)$Array.GmaxHatY
    E.GmaxY.stat <- EGmaxY.estimation.function.step2(GmaxHatY.stat, bandwidth = bandwithInf)$Array.EGmaxHatY
    Pstat_so_j <- as.matrix(E.GmaxY.stat[,,siteid])
    lambdastat_so_j <- (1 - Pstat_so_j)/Pstat_so_j
    FjYj.stat <- GiYi.func(Data = Data, h = bandwithInf)
    lambdaPickands.stat <- Version2_V.applied.lambda.estim.func(GiYi.mat = FjYj.stat,lambda = lambda_so_j, use.timesteps = use.timesteps)
    lambda_stat_so <- lambdaPickands.stat$lambdaV
    Precord_stat_so <- 1/(1+lambda_stat_so)
    Pickands_stat_so <- lambdaPickands.stat$A
    statObject_list <- list(P.record.stat = Precord_stat_so,Astat = Pickands_stat_so)
    returnObject_list <- c(returnObject_list,statObject_list)
  }
  return(returnObject_list)
  
}

# Version2_V.applied.lambda.estim.func<- function(GiYi.mat, lambda , use.timesteps = 0 ){ 
#   if(is.numeric(lambda)==TRUE){lambda <- as.matrix(lambda,nrow=1)}
#   nb.evaluations <- dim(lambda)[1]
#   vw <- rep(NA,nb.evaluations)
#   cw <- rep(NA,nb.evaluations)
#   Aw <- rep(NA,nb.evaluations)
#   Vw <- rep(NA,nb.evaluations)
#   lambdaV <- rep(NA,nb.evaluations)
#  
#   #zi.mat  <- -1/(log((lambda)^(-1))) 
#   # zi.mat  <- (-log(lambda))^(-1)
#   # zi.mat  <- (lambda)^(-1) # maybe it is not like that, I don't remember, LOOK FOR IT LATER
#   #zi.mat  <- 1/log(lambda)
#   #zi.mat  <- 1/exp(-lambda)
#   #zi.mat  <- exp(-lambda)
#   tt<-1:nb.evaluations
#   zi.mat  <- 1/(lambda*tt)
#   #zi.mat  <- 1/lambda
#   
#   w.mat <- t(apply(zi.mat,1, function(x) x/sum(x)))
#   # w.mat <- zi.mat/(apply(zi.mat,1,sum)) # maybe it is not like that, I don't remember, LOOK FOR IT LATER
#   
#   GiYi.inform  <- GiYi.mat
#   #GiYi.inform  <- exp(-1/GiYi.mat)
#   # GiYi.inform  <- exp(GiYi.mat)
#   
#   # use.timesteps = 0 means we will use all the information to construct V (we may want to use only de central time steps)
#   if (use.timesteps[1] != 0 ){
#     GiYi.inform <- GiYi.mat[use.timesteps,] } # truncate the matrix used
#   nbsites <- dim(GiYi.inform)[2]
#   nbcolumns <- dim(GiYi.inform)[1]
#   
#   for(ti in 1:nb.evaluations){ 
#     wi.eval <- as.numeric(w.mat[ti,])
#     FiYiwi <- t(t(GiYi.inform)^(1/wi.eval))
#     Mean.FiYiwi<- rowMeans(FiYiwi) # a vector of length nbcolumns
#     Max.FiYiwi<- apply(FiYiwi,1,max) # a vector of length nbcolumns
#     vw.ti <- mean(Max.FiYiwi - Mean.FiYiwi)
#     vw[ti] <- vw.ti
#     cw.ti <- mean(wi.eval/(1+wi.eval))
#     cw[ti] <- cw.ti
#     Aw.ti <- (vw.ti + cw.ti)/ (1 - vw.ti - cw.ti)
#     Aw[ti] <- Aw.ti
#     sum.zi.ti <- sum(zi.mat[ti,]^(-1)) # wi = xi/(sumxi)
#     lambdaV[ti] <- sum.zi.ti * Aw.ti # lambdaV.ti
#   }
#   #lambdaV <- lambdaV*tt
#   return(list("lambdaV"=lambdaV,"A"=Aw))
# }

get_intercept <- function(x1, y1, slope) {
  # Point slope formula (y - y1) = slope(x - x1)
  y_intercept = slope * (- x1) + y1
  return(y_intercept)
}

