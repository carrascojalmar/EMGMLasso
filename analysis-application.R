#
source("functions-application.R")

#------------------------------
#  Data
#------------------------------

data("respInf")
#head(respInf)

age.p <- scale(respInf$age,center=TRUE, scale=TRUE)
age2.p <- scale(respInf$age**2,center=TRUE, scale=TRUE)

age1.p <- scale(respInf$age1,center=TRUE, scale=TRUE)
age12.p <- scale(respInf$age1**2,center=TRUE, scale=TRUE)

time2.p <- scale(respInf$time2,center=TRUE, scale=TRUE)
time22.p <- scale(respInf$time2**2,center=TRUE, scale=TRUE)

mi <- as.numeric(table(respInf$id))

x <- data.frame(respInf$female,respInf$height,respInf$cosine,
           respInf$sine,respInf$xero, age.p,age2.p,age1.p,
           age12.p,time2.p,time22.p)
colnames(x) <- c("Female","Height","Cosine","Sine","Xero",
                 "age.p","age2.p","age1.p",
                 "age12.p","time2.p","time22.p")
x <- data.matrix(x)
x[,1] <- ifelse(x[,1]==1,0,1)
x[,5] <- ifelse(x[,5]==1,0,1)
Y <- respInf$time
m <- dim(respInf)[1]
nj <- mi
w <- matrix(rep(1,sum(mi)),ncol=1,nrow=sum(mi))
DD <- matrix(1,1,1)

#------------------------------
# Estimation - probit
#------------------------------

mod0 <- EMN_ProbitLogit(nj=nj,y=Y,x=x,w=w,eps=0.000001,
                        iter.max=1000,Q=100, 
                                    type="probit")
se <- mod0$ep
z <- mod0$teta/se
pvalue <- 2 * (1 - stats::pnorm(abs(z)))
TAB <- cbind(Estimate = mod0$teta, Std.Error = se,
             z.value = z, `Pr(>|z|)` = pvalue)
row.names(TAB) <- c("Female","Height","Cosine","Sine","Xero",
                    "age.p","age2.p","age1.p",
                    "age12.p","time2.p","time22.p",
                    "lambda")
mTab <- list(Lik=mod0$lik,Coefficients = round(TAB,digits = 3))
mTab

#------------------------------
# Selection
# glmnet
#------------------------------

mod01 <- cv.glmnet(x,Y,family=binomial(link=probit),
                   intercept=FALSE,alpha=1)

mod02 <- glmnet(x,Y,family=binomial(link=probit),
                   intercept=FALSE,alpha=1,
                lambda=mod01$lambda.min)

mod02$beta
#plot(mod01)

#------------------------------
# Selecetion
# EMGMLasso
#------------------------------

mod1 <- EMGMLasso(nj=nj,y=Y,x=x,w=w,
                       eps=0.00001,
                       iter.max=1000,
                       Intercep = FALSE,
                       folds=10,
                       type="probit",
                       Q=30)
mod1

#------------------------------
# Selection
# glmmLasso
#------------------------------

ind <- rep(1:275, mi)
dados <- data.frame(Y,x,ind=factor(ind))
dados$ind <- factor(dados$ind)
colnames(dados)<- c("Y","Female","Height","Cosine","Sine","Xero",
                    "age.p","age2.p","age1.p",
                    "age12.p","time2.p","time22.p",
                    "ind")
lambda <- seq(500,0,by=-5)
family <- binomial(probit)
BIC_vec<-rep(Inf,length(lambda))

for(j in 1:length(lambda)){
  glm1 <- try(glmmLasso(Y~-1+Female+Height+Cosine+
  Sine+Xero+age.p+age2.p+age1.p+age12.p+time2.p+time22.p,
                        rnd=list(ind=~1),family=family,
                        data=dados,
                        lambda=lambda[j],
                        control=list(center=FALSE)),
              silent=TRUE)  
  if(!inherits(glm1, "try-error")){BIC_vec[j]<-glm1$bic}
}

opt <- which.min(BIC_vec)
mod2 <- glmmLasso(Y~-1+Female+Height+Cosine+
                    Sine+Xero+age.p+age2.p+age1.p+
                    age12.p+time2.p+time22.p,
                  rnd=list(ind=~1),  
                  family=family,data=dados,
                  lambda=lambda[opt],
                  control=list(center=FALSE))
summary(mod2)  


#------------------------------
# High-dimension
#------------------------------

np <- 1000
x.1aux <- runif(275*np,-1,1)
x.aux <- matrix(x.1aux,nrow=275,ncol=np,byrow = T)
z.aux <- (x.aux-apply(x.aux,2,mean))/apply(x.aux,2,sd)
mX <- matrix(NA,nrow=sum(nj),ncol=np)
for(i in 1:np){mX[,i] <- rep(z.aux[,i],nj)}
x.np <- cbind(x,mX)

#------------------------------
# Selection
# glmnet
#------------------------------

mod01.np <- cv.glmnet(x.np,Y,family=binomial(link=probit),
                      intercept=FALSE,alpha=1)
lambda.np <- mod01.np$lambda.min
mod02.np <- glmnet(x.np,Y,family=binomial(link=probit),
                   lambda=lambda.np,
                   intercept=TRUE,alpha=1)
mod02.np$beta[1:11]

#------------------------------
# Selection
# EMGMLasso
#------------------------------

mod1.np <- EMGMLasso(nj=nj,y=Y,x=x.np,w=w,
                          eps=0.000001,
                          iter.max=1000,folds=10,
                          type="probit",
                          Intercep = FALSE,Q=30)
mod1.np$beta[1:11]

#------------------------------
# Selection
# glmmLasso
#------------------------------

dados2 <- data.frame(Y,x.np,ind=factor(rep(1:275, mi)))

form.aux1 <- paste0(colnames(dados2)[2:(np+12)],
                    collapse="+")
form.aux2 <- paste("-1",form.aux1,sep="+")
form.aux3 <- paste(colnames(dados2)[1],form.aux2,sep="~")
form <- as.formula(form.aux3)

lambda <- seq(500,0,by=-5)
family <- binomial(probit)
BIC_vec<-rep(Inf,length(lambda))
for(j in 1:length(lambda)){
  glm1 <- try(glmmLasso(form,
                        rnd = list(ind=~1), 
                        data=dados2,  
                        family=family,lambda=lambda[j],
                        control=list(center=FALSE)),silent=TRUE)  
  if(!inherits(glm1, "try-error")){BIC_vec[j]<- glm1$bic}
}

opt <- which.min(BIC_vec)

mod2 <- glmmLasso(form,
                  rnd = list(ind=~1), data=dados2, 
                  family = binomial(link = "probit"), 
                  control=list(center=FALSE),
                  lambda=lambda[opt])
summary(mod2) 
