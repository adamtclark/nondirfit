#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, version 3.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#Author: Adam Clark, University of Minnesota, 2016

#Example:
if(FALSE) {
  source("nondirfit.R") #load R function
  
  #Make fake data:
  x1<-2+rnorm(100); x2<-7+3*x1+rnorm(100)
  plot(x1~x2)
  vardf<-data.frame(x1, x2)
  
  #Run analysis
  out<-nondirfit(vardf)
  
  #Extract parameters
  parms<-out$pars
  #Reformat for x1 = b0 + b1*x2
  (reform_parms<-parms[c("(Intercept)", "x2")]/(-parms[c("x1")]))
  abline(a=reform_parms["(Intercept)"], b=reform_parms["x2"])
  
  #Compare to linear model
  mod<-lm(x1~x2)
  coef(mod)
  abline(mod, lty=2)
  #Should match closely, but not perfectly
}



nondirfit<-function(vardf, doscale=TRUE) {
  #Runs a type II regression among all variables in vardf
  #vardf is a dataframe with columns of variables to be fit nodirectionally
  
  #Function automatically scales all variables, then returns back-tranformed parameter estimates
  #note, output$pars is in the format 0 = b1*x1+b2*x2+b3*x3...
  
  var<-vardf
  if(ncol(var)<=1 | nrow(var)<=1) {
    stop("too few variables to regress")
  }
  
  if(is.null(colnames(var))) {
    if(ncol(var)>3)
      stop("too many columns - don't even think about it")
    colnames(var)<-letters(1:ncol(var))
  }
  
  #scale data
  if(doscale) {
    sdmat<-apply(var, 2, function(x) sd(x, na.rm=T))
    mnmat<-colMeans(var, na.rm=T)
  } else {
    sdmat<-rep(1, ncol(var)); names(sdmat)=colnames(var)
    mnmat<-rep(0, ncol(var)); names(mnmat)=colnames(var)
  }
  
  varsc<-data.frame(t((t(var)-mnmat)/sdmat))
  
  vn<-names(var)
  
  #fit linear model to get initial guess of parameters
  expr<-paste("lm(", vn[1], "~", paste(vn[-1], collapse="+"), ", data=varsc)", sep="")    
  modc<-coef(eval(parse(text=expr)))
  modc[vn[1]]<-(-1)
  
  #optimize nondirectionally
  optc<-optim(par = modc, fn = optfun, gr = NULL, extraparms=list(varsc=varsc, vn=vn))$par
  optc<-optc/optc[2]
  
  #Calcuate predicted (and unscaled) values
  hatv<-matrix(nrow=nrow(varsc), ncol=ncol(varsc)); colnames(hatv)<-vn
  hatvsc<-matrix(nrow=nrow(varsc), ncol=ncol(varsc)); colnames(hatvsc)<-vn
  for(i in 1:ncol(hatv)) {
    pos<-which(names(optc)==vn[i])
    hatvsc[,i]<-(optc[-c(1, pos)]%*%t(varsc[-i])+optc[1])/(-optc[pos])
    hatv[,i]<-(hatvsc[,i]*sdmat[i])+mnmat[i]
  }
  hatv<-data.frame(hatv)
  hatvsc<-data.frame(hatvsc)
  
  #project variables onto fit line or plane
  #given ax+by+cz+...+int = 0, and point {x0, y0, z0...}
  #k = (ax0+by0+cz0+...+int)/(a^2+b^2+c^2...)
  #xnew = x0-a*k
  
  k<-(optc[c(1, match(vn, names(optc)))]%*%t(cbind(1, varsc)))/sum(optc[-1]^2)
  possnapsc<-t(t(varsc)-optc[match(vn, names(optc))]*t(matrix(ncol=length(vn), nrow=length(k), k)))
  possnap<-array(dim=dim(possnapsc))
  
  for(i in 1:ncol(possnapsc)) {
    possnap[,i]<-possnapsc[,i]*sdmat[i]+mnmat[i]
  }
  
  #Transform parms into non-scaled space
  optc_notsc<-numeric(length(optc)); names(optc_notsc)<-names(optc)
  optc_notsc[1]<-optc[1]
  for(i in 2:length(optc)) {
    scpos<-which(names(sdmat)==names(optc[i]))
    optc_notsc[i]<-optc[i]/sdmat[scpos]
    optc_notsc[1]<-optc_notsc[1]-mnmat[scpos]/sdmat[scpos]*optc[i]
  }
  
  colnames(possnap)<-colnames(var)
  
  #Get estimate of model fit
  SSres<-sum((possnapsc-varsc)^2, na.rm=T)
  SStot<-sum((t(varsc)-colMeans(varsc, na.rm=T))^2, na.rm=T) #NB - ybar ~ 0
  rsq_est<-(1-SSres/SStot)
  p<-ncol(varsc); n<-nrow(varsc)
  rsq_est_adj<-rsq_est-(1-rsq_est)*p/(n-p-1)
  
  return(list(possnap=possnap, pred=hatv, vars=var, pars=optc_notsc, rsq=list(rsq_est=rsq_est, rsq_est_adj=rsq_est_adj),
              scl=list(sdmat=sdmat, mnmat=mnmat, possnapsc=possnapsc, predsc=hatvsc, varsc=varsc, parssc=optc)))
  #regular and scaled lists
  #possnap is values snapped to tradeoff
  #hatv is predicted values
  #pars is pars for nonparametric regression (in scaled space!!)
  #scl is scaling parameters for xscal = (x-mean(x))/sd(x)
}

#The following are helper functions needed by the nondirfit function
optfun<-function(modc, extraparms) {
  k<-(modc[c(1, match(extraparms$vn, names(modc)))]%*%t(cbind(1, extraparms$varsc)))/sum(modc[-1]^2)
  possnapsc<-t(t(extraparms$varsc)-modc[match(extraparms$vn, names(modc))]*t(matrix(ncol=length(extraparms$vn), nrow=length(k), k)))
  
  diffobs<-sum((rowSums((possnapsc-extraparms$varsc)^2)), na.rm=T)
  
  return(diffobs)
}

get_rsq_nondir<-function(trout, colnums) {
  possnapsc<-trout$scl$possnapsc[,colnums]
  varsc<-trout$scl$varsc[,colnums]
  
  SSres<-sum((possnapsc-varsc)^2, na.rm=T)
  SStot<-sum((t(varsc)-colMeans(as.matrix(varsc),na.r=T))^2, na.rm=T) #NB - ybar ~ 0
  rsq_est<-(1-SSres/SStot)
  p<-ncol(as.matrix(varsc)); n<-nrow(as.matrix(varsc))
  rsq_est_adj<-rsq_est-(1-rsq_est)*p/(n-p-1)
  
  return(list(rst_est=rsq_est, rsq_est_adj=rsq_est_adj))
}
