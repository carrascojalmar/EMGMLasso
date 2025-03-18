#
source("functions-application.R")

#------------------------------
#  Data
#------------------------------

data("respInf")
head(respInf,n=20)

age1 <- respInf$age1
age12 <- respInf$age1**2

time2.p <- scale(respInf$time2,center=TRUE)
time22.p <- scale(respInf$time2**2,center=TRUE)

mi <- as.numeric(table(respInf$id))

x <- data.frame(respInf$female,respInf$height,respInf$cosine,
           respInf$sine,respInf$xero,respInf$stunted,age1,age12,
           time2.p,time22.p)
colnames(x) <- c("Female","Height","Cosine","Sine","Xero","stunted",
                 "Age1","Age12","time2.p","time2.p")
x <- data.matrix(x)
Y <- respInf$time
m <- dim(respInf)[1]
nj <- mi
w <- matrix(rep(1,sum(mi)),ncol=1,nrow=sum(mi))
DD <- matrix(1,1,1)

#------------------------------
# High-dimension
#------------------------------

#set.seed(1234)
np <- 1000
x.1aux <- rnorm(275*np,0,1)
x.aux <- matrix(x.1aux,nrow=275,ncol=np,byrow = T)
z.aux <- x.aux#scale(x.aux,center=TRUE,scale = TRUE)
mX <- matrix(NA,nrow=sum(nj),ncol=np)
colnames(mX) <- paste0("x", 1:np)
for(i in 1:np){mX[,i] <- rep(z.aux[,i],nj)}
x.np <- cbind(x,mX)
colnames(x.np) <- c(colnames(x),colnames(mX))

#------------------------------
# Selection and estimation
# EMGMLasso
#------------------------------

mod4 <- EMGMLasso_ProbitLogit(nj=nj,y=Y,x=x.np,w=w,eps=1E-5,
                              iter.max=1E3,folds=7,
                              type="probit",
                              Intercep = FALSE,
                              Q=100)

#saveRDS(mod4, file = "mod4.rds")

#------------------------------
# Estimation - Probit
#------------------------------

new.x <- x.np[,which(mod4$betas != 0 & mod4$betas < 1E-5)]

mod5 <- EMN_ProbitLogit(nj=nj,y=Y,x=cbind(1,new.x),w=w,eps=1E-5,
                        iter.max=1E3,Q=1E2,type="probit")
se <- mod5$ep
z <- mod5$teta/se
pvalue <- 2 * (1 - stats::pnorm(abs(z)))
TAB <- data.frame(Estimate = mod5$teta, Std.Error = se,
                  z.value = z, pvalue = pvalue)
mTab <- list(Lik=mod5$loglik,Coefficients = round(TAB,digits = 3))
row.names(TAB) <- c("Intercept",colnames(new.x),"lambda")

#saveRDS(dplyr::filter(TAB,pvalue<0.06), file = "mod5.rds")

#------------------------------
# Selection and estimation
# glmmLasso
#------------------------------

dados2 <- data.frame(Y,x.np,ind=factor(rep(1:275, mi)))

form.aux1 <- paste0(colnames(dados2)[2:(np+11)],
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

form.aux1 <- paste0(colnames(dados2)[2:(np+11)],
                    collapse="+")
form.aux3 <- paste(colnames(dados2)[1],form.aux1,sep="~")
form <- as.formula(form.aux3)

mod6 <- glmmLasso(form,
                  rnd = list(ind=~1), data=dados2, 
                  family = binomial(link = "probit"), 
                  control=list(center=FALSE),
                  final.re=TRUE,
                  lambda=lambda[opt])

#saveRDS(summary(mod6), file = "mod6.rds")

