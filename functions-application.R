#
# Packages
#
suppressMessages(library(mvtnorm))
suppressMessages(library(tmvtnorm))
suppressMessages(library(glmnet))
suppressMessages(library(dplyr))
suppressMessages(library(glmmLasso))
suppressMessages(library(sandwich))
suppressMessages(library(readr))
suppressMessages(library(statmod))
suppressMessages(library(gamlss.data))

#
# Functions
#

# The likelihood function
logver.QGauss <-function(nj,y,x,w,betas,DD,Q){
  p <-ncol(x)
  q1 <-ncol(w)
  m <-length(nj)
  out <- gauss.quad(Q,"hermite")
  nu1 <- out$nodes
  ww1 <- out$weights
  loglik <- matrix(0,m,1)
  for(j in 1:m){
    y1 <- y[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]
    x1 <- matrix(x[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=p)
    w1 <- matrix(w[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=q1)
    
    aux <- c(x1%*%betas)+(nu1)%*%t(w1)
    auxp <- pnorm(aux,log = T)
    auxp2 <- pnorm(-aux,log = T)
    res <- auxp%*%y1+auxp2%*%(1-y1)
    loglik[j]<- log(sum(ww1*exp(res)))
  }
  return(sum(loglik))
}
logverLogit.QGauss <-function(nj,y,x,w,betas,DD,Q){
  p<-ncol(x)
  q1<-ncol(w)
  m<-length(nj)
  
  out <- gauss.quad(Q,"hermite")
  nu1 <- out$nodes
  ww1 <- out$weights
  
  loglik<- matrix(0,m,1)
  
  for(j in 1:m){
    y1=y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))  ]
    x1=matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
    w1=matrix(w[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
    
    aux<-c(x1%*%betas)+(nu1)%*%t(w1)
    auxp<- plogis(aux,log.p = TRUE)
    auxp2<- plogis(-aux,log.p = TRUE)
    res<- auxp%*%y1+auxp2%*%(1-y1)
    loglik[j]<-log(sum(ww1*exp(res)))
  }
  return(sum(loglik))
}
# Bernoulli mixed models by EM algorithm
EMN_ProbitLogit <- function(nj, y, x, w, eps,iter.max,Q,type="probit"){
  p<-ncol(x)
  q1<-ncol(w)
  m<-length(nj)
  N<-sum(nj)
  criterio<-1
  betas<- solve(t(x)%*%x)%*%t(x)%*%y
  DD=diag(q1)
  teta <- c(betas,DD[upper.tri(DD, diag = T)])
  cont<-0
  
  if(type=="probit"){
    
    while(criterio>eps){
      cont<-cont+1
      print(cont)
      
      suma1<-matrix(0,p,p)
      suma2<-matrix(0,p,1)
      suma3<-matrix(0,q1,q1)
      suma4<-matrix(0,p+1,p+1)
      aux3<-matrix(0,p+1,1)
      
      for(j in 1:m){
        y1=y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))  ]
        x1=matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
        w1=matrix(w[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
        
        Omegai<-w1%*%DD%*%t(w1)+diag(nj[j])
        Deltai<-DD%*%t(w1)%*%solve(Omegai)
        gammai<-x1%*%betas
        Lambdai<-DD-DD%*%t(w1)%*%solve(Omegai)%*%w1%*%DD
        
        if(length(y1)==1){y1 <- matrix(y1)}
        
        ##Media e Variancia:
        
        Ai<-diag(1-2*y1)
        aux2<-Ai%*%Omegai%*%Ai
        aux2<-(aux2+t(aux2))/2
        aux<-mtmvnorm(mean=c(Ai%*%gammai), sigma=aux2, upper=rep(0,nj[j]))
        M1<-Ai%*%aux$tmean
        M2<-Ai%*%aux$tvar%*%Ai+ M1%*%t(M1)
        
        bi<-Deltai%*%(M1-gammai)
        b2i<-Lambdai+Deltai%*%(M2+gammai%*%t(gammai)-M1%*%t(gammai)-t(M1%*%t(gammai)))%*%t(Deltai)
        suma1<-suma1+t(x1)%*%x1
        suma2<-suma2+t(x1)%*%(M1-w1%*%bi)
        suma3<-suma3+b2i
        aux3[1:p]<-t(M1-w1%*%bi-gammai)%*%x1
        aux3[p+1]<- -0.5*sum(diag(solve(DD)%*%matrix(1,1,1)))+
          0.5*sum(diag(solve(DD)%*%matrix(1,1,1)%*%solve(DD)%*%b2i))
        suma4<-suma4+(aux3)%*%t(aux3)
        
      }
      teta1<-teta
      betas <- solve(suma1)%*%suma2 
      DD<- suma3/m
      teta <- c(betas,DD[upper.tri(DD, diag = T)])
      criterio <- (teta1-teta)%*%(teta1-teta)#sum((teta1-teta)^2)
      if(cont==iter.max){criterio=eps/10}
      logver<-logver.QGauss(nj,y,x,w,betas,(DD+t(DD))/2,Q)
      epbetas<-sqrt(diag(solve(suma4)))
    }
  }
  
  
  if(type=="logit"){
    
    c=16*sqrt(3)/(15*pi)
    while(criterio>eps){
      cont<-cont+1
      print(cont)
      
      suma1<-matrix(0,p,p)
      suma2<-matrix(0,p,1)
      suma3<-matrix(0,q1,q1)
      suma4<-matrix(0,p+1,p+1)
      aux3<-matrix(0,p+1,1)
      
      for(j in 1:m){
        y1=y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))  ]
        x1=matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
        w1=matrix(w[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
        
        x1=c*x1
        w1=c*w1
        
        Omegai<-w1%*%DD%*%t(w1)+diag(nj[j])
        Deltai<-DD%*%t(w1)%*%solve(Omegai)
        gammai<-x1%*%betas
        Lambdai<-DD-DD%*%t(w1)%*%solve(Omegai)%*%w1%*%DD
        
        if(length(y1)==1){y1 <- matrix(y1)}
        
        ##Media e Variancia:
        
        Ai<-diag(1-2*y1)
        aux2<-Ai%*%Omegai%*%Ai
        aux2<-(aux2+t(aux2))/2
        aux<-mtmvnorm(mean=c(Ai%*%gammai), sigma=aux2, upper=rep(0,nj[j]))
        M1<-Ai%*%aux$tmean
        M2<-Ai%*%aux$tvar%*%Ai+ M1%*%t(M1)
        
        bi<-Deltai%*%(M1-gammai)
        b2i<-Lambdai+Deltai%*%(M2+gammai%*%t(gammai)-M1%*%t(gammai)-t(M1%*%t(gammai)))%*%t(Deltai)
        suma1<-suma1+t(x1)%*%x1
        suma2<-suma2+t(x1)%*%(M1-w1%*%bi)
        suma3<-suma3+b2i
        aux3[1:p]<-t(M1-w1%*%bi-gammai)%*%x1
        aux3[p+1]<- -0.5*sum(diag(solve(DD)%*%matrix(1,1,1)))+
          0.5*sum(diag(solve(DD)%*%matrix(1,1,1)%*%solve(DD)%*%b2i))
        suma4<-suma4+(aux3)%*%t(aux3)
        
      }
      teta1<-teta
      betas <- solve(suma1)%*%suma2 
      DD<- suma3/m
      teta <- c(betas,DD[upper.tri(DD, diag = T)])
      criterio <- (teta1-teta)%*%(teta1-teta)#sum((teta1-teta)^2)
      if(cont==iter.max){criterio=eps/10}
      logver<-logverLogit.QGauss(nj,y,x,w,betas,(DD+t(DD))/2,Q)
      epbetas<-sqrt(diag(solve(suma4)))
    }
  }
  
  
  return(list(teta=teta,ep=epbetas, loglik=logver))
}
# PML estimation in the Bernoulli mixed model
EMGMLasso_ProbitLogit <-function(nj,y,x,w,eps,iter.max,Intercep=FALSE,folds,Q,type="probit"){
  q1 <- ncol(w)
  m <- length(nj)
  N <- sum(nj)
  criterio <- 1
  
  if(Intercep==FALSE){
    p <- ncol(x)  
    mod0 <- cv.glmnet(x,y,family=binomial(link=probit),
                      intercept=FALSE,nfolds=folds,alpha=1)
    lambda0 <- mod0$lambda.min
    betaI <- glmnet(x,y,family=binomial(link=probit),
                    intercept=FALSE,alpha=1,
                    lambda=lambda0)
    betas <- matrix(betaI$beta,nrow=p,ncol=1)  
  }else{
    p <- ncol(x)+1
    mod0 <- cv.glmnet(x,y,family=binomial(link=probit),
                      intercept=TRUE,
                      nfolds=folds,
                      alpha=1)
    lambda0 <- mod0$lambda.min
    betaI <- glmnet(x,y,family=binomial(link=probit),intercept=TRUE,
                    alpha=1,
                    lambda=lambda0)
    betas <- matrix(c(as.numeric(betaI$a0),as.numeric(betaI$beta)),p+1,1)}
  
  DD <- diag(q1)
  if(type=="probit"){
    loglik0 <- logver.QGauss(nj,y,x,w,betas,(DD+t(DD))/2,Q) 
    teta <- c(betas,DD[upper.tri(DD,diag=T)])
    cont<-0
    while(criterio>eps){
      cont<-cont+1
      suma3<-matrix(0,q1,q1)
      if(Intercep==FALSE){x.aux <- x}else{x.aux <- cbind(1,x)}
      
      zm <-NULL
      for(j in 1:m){
        y1=y[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]
        x1=matrix(x.aux[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=p)
        w1=matrix(w[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=q1)
        
        Omegai<-w1%*%DD%*%t(w1)+diag(nj[j])
        Deltai<-DD%*%t(w1)%*%solve(Omegai)
        gammai<-x1%*%betas
        Lambdai<-DD-DD%*%t(w1)%*%solve(Omegai)%*%w1%*%DD
        
        ##Media e Variancia:
        if(length(y1)==1){y1 <- matrix(y1)}
        
        Ai<-diag(1-2*y1)
        #aux<-MomemNT(Ai%*%gammai,Ai%*%Omegai%*%Ai,rep(0,nj[j]))
        aux2<-Ai%*%Omegai%*%Ai
        aux2<-(aux2+t(aux2))/2
        aux<-mtmvnorm(mean=c(Ai%*%gammai), sigma=aux2, upper=rep(0,nj[j]))
        M1<-Ai%*%aux$tmean
        M2<-Ai%*%aux$tvar%*%Ai+ M1%*%t(M1)
        
        bi<-Deltai%*%(M1-gammai)
        b2i<-Lambdai+Deltai%*%(M2+gammai%*%t(gammai)-M1%*%t(gammai)-
                                 t(M1%*%t(gammai)))%*%t(Deltai)
        suma3<- suma3+b2i
        y1m<- M1-w1%*%bi
        zm <- rbind(zm,y1m)
      }
      
      teta1<-teta
      
      if(Intercep==FALSE){
        mod01 <- cv.glmnet(x,zm,family="gaussian",intercept=FALSE,
                           nfolds=folds,
                           alpha=1)
        lambda.atual <- mod01$lambda.min
        la.eq2 <- glmnet(x,zm,family="gaussian",lambda=lambda.atual,
                         intercept=FALSE,alpha=1)
        betas <- matrix(la.eq2$beta,p,1)  
      }else{
        mod01 <- cv.glmnet(x,zm,family="gaussian",intercept=TRUE,
                           nfolds=folds,
                           alpha=1)
        lambda.atual <- mod01$lambda.min
        la.eq2 <- glmnet(x,zm,family="gaussian",lambda=lambda.atual,
                         intercept=TRUE,alpha=1)
        betas <- matrix(c(as.numeric(la.eq2$a0),as.numeric(la.eq2$beta)),p,1)}
      
      plot(mod01)
      
      DD <- suma3/m
      teta <- c(betas,DD[upper.tri(DD, diag = T)])
      
      if(Intercep==FALSE){
        loglik <- logver.QGauss(nj,y,x,w,betas,(DD+t(DD))/2,Q)  
      }else{
        loglik <- logver.QGauss(nj,y,cbind(1,x),w,betas,(DD+t(DD))/2,Q)  
      }
      
      print(loglik)
      
      #loglik1 <- loglik0
      #a.k <- (loglik-loglik1)/(loglik1-loglik)
      #loglik.ass <-  loglik1+(1/(1-a.k))*(loglik-loglik1)
      #criterio <- abs(loglik-loglik.ass)
      
      criterio <- (teta1-teta)%*%(teta1-teta)#sum((teta1-teta)^2)
      print(cont)
      if(cont==iter.max){criterio=eps/10}
    }
    
    npar <- length(c(teta))-length(which(betas == 0))
    
    if(Intercep==FALSE){
      loglik <- logver.QGauss(nj,y,x,w,betas,(DD+t(DD))/2,Q)  
    }else{
      loglik <- logver.QGauss(nj,y,cbind(1,x),w,betas,(DD+t(DD))/2,Q)  
    }
    
    AICc<- -2*loglik +2*npar
    BICc <- -2*loglik +log(N)*npar
    
    return(list(betas=betas, DD=DD, loglik=loglik, AIC=AICc, BIC=BICc, 
                lambda=lambda.atual,count=cont,model=la.eq2))
  }
  if(type=="logit"){
    loglik0 <- logverLogit.QGauss(nj,y,x,w,betas,(DD+t(DD))/2,Q) 
    teta <- c(betas,DD[upper.tri(DD,diag=T)])
    cont<-0
    c=16*sqrt(3)/(15*pi)
    while(criterio>eps){
      cont<-cont+1
      suma3<-matrix(0,q1,q1)
      if(Intercep==FALSE){x.aux <- x}else{x.aux <- cbind(1,x)}
      
      zm <-NULL
      for(j in 1:m){
        y1<-y[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]
        x1<-matrix(x.aux[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=p)
        w1<-matrix(w[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=q1)
        
        x1=c*x1
        w1=c*w1
        
        Omegai<-w1%*%DD%*%t(w1)+diag(nj[j])
        Deltai<-DD%*%t(w1)%*%solve(Omegai)
        gammai<-x1%*%betas
        Lambdai<-DD-DD%*%t(w1)%*%solve(Omegai)%*%w1%*%DD
        
        ##Media e Variancia:
        if(length(y1)==1){y1 <- matrix(y1)}
        
        Ai<-diag(1-2*y1)
        #aux<-MomemNT(Ai%*%gammai,Ai%*%Omegai%*%Ai,rep(0,nj[j]))
        aux2<-Ai%*%Omegai%*%Ai
        aux2<-(aux2+t(aux2))/2
        aux<-mtmvnorm(mean=c(Ai%*%gammai), sigma=aux2, upper=rep(0,nj[j]))
        M1<-Ai%*%aux$tmean
        M2<-Ai%*%aux$tvar%*%Ai+ M1%*%t(M1)
        
        bi<-Deltai%*%(M1-gammai)
        b2i<-Lambdai+Deltai%*%(M2+gammai%*%t(gammai)-M1%*%t(gammai)-
                                 t(M1%*%t(gammai)))%*%t(Deltai)
        suma3<- suma3+b2i
        y1m<- M1-w1%*%bi
        zm <- rbind(zm,y1m)
      }
      
      teta1<-teta
      
      if(Intercep==FALSE){
        mod01 <- cv.glmnet(x,zm,family="gaussian",intercept=FALSE,
                           nfolds=folds,
                           alpha=1)
        lambda.atual <- mod01$lambda.min
        la.eq2 <- glmnet(x,zm,family="gaussian",lambda=lambda.atual,
                         intercept=FALSE,alpha=1)
        betas <- matrix(la.eq2$beta,p,1)  
      }else{
        mod01 <- cv.glmnet(x,zm,family="gaussian",intercept=TRUE,
                           nfolds=folds,
                           alpha=1)
        lambda.atual <- mod01$lambda.min
        la.eq2 <- glmnet(x,zm,family="gaussian",lambda=lambda.atual,
                         intercept=TRUE,alpha=1)
        betas <- matrix(c(as.numeric(la.eq2$a0),as.numeric(la.eq2$beta)),p,1)}
      
      plot(mod01)
      
      DD <- suma3/m
      teta <- c(betas,DD[upper.tri(DD, diag = T)])
      
      if(Intercep==FALSE){
        loglik <- logverLogit.QGauss(nj,y,x,w,betas,(DD+t(DD))/2,Q)  
      }else{
        loglik <- logverLogit.QGauss(nj,y,cbind(1,x),w,betas,(DD+t(DD))/2,Q)  
      }
      
      print(loglik)
      
      #loglik1 <- loglik0
      #a.k <- (loglik-loglik1)/(loglik1-loglik)
      #loglik.ass <-  loglik1+(1/(1-a.k))*(loglik-loglik1)
      #criterio <- abs(loglik-loglik.ass)
      
      criterio <- (teta1-teta)%*%(teta1-teta)#sum((teta1-teta)^2)
      print(cont)
      if(cont==iter.max){criterio=eps/10}
    }
    
    npar <- length(c(teta))-length(which(betas == 0))
    
    if(Intercep==FALSE){
      loglik <- logverLogit.QGauss(nj,y,x,w,betas,(DD+t(DD))/2,Q)  
    }else{
      loglik <- logverLogit.QGauss(nj,y,cbind(1,x),w,betas,(DD+t(DD))/2,Q)  
    }
    
    AICc<- -2*loglik +2*npar
    BICc <- -2*loglik +log(N)*npar
    
    return(list(betas=betas, DD=DD, loglik=loglik, AIC=AICc, BIC=BICc, 
                lambda=lambda.atual,count=cont,model=la.eq2))
  }
  
  
}

