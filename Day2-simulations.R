
### Longitudinal Data Analysis with R Statistical Software
### R simulation studies
### Benjamin French and Bryan Shepherd
### Vanderbilt-Nigeria Biostatistics Workshop
### April 27 -- May 1, 2026

## General set up

library(ggplot2)
library(geepack)
library(lme4)

n <- 50
m <- 4
beta <- c(70, 5, 7)
sigma <- sqrt(10)
D <- sqrt(15)

x <- rep(c(0,1), each=n/2)
b0 <- rnorm(n, 0, D)

data <- data.frame(id=rep(1:n, each=m),
                   x=rep(x, each=m),
                   t=rep((1:m)-1, n),
                   b0=rep(b0, each=m),
                   e=rnorm(n*m, 0, sigma))

data <- within(data, {
  y = beta[1] + b0 + beta[2]*t*(x==0) + beta[3]*t*(x==1) + e
})

## Day 2

# Goal: Compare GEE to linear regression

m_lm <- lm(y ~ x*t, data=data)
summary(m_lm)

m_gee <- geeglm(y ~ x*t, data=data, id=id, corstr="independence")
summary(m_gee)

m_gee_e <- geeglm(y ~ x*t, data=data, id=id, corstr="exchangeable")
summary(m_gee_e)

m_gee <- geeglm(y ~ factor(x):t, data=data, id=id, corstr="independence")
summary(m_gee)

m_gee_e <- geeglm(y ~ factor(x):t, data=data, id=id, corstr="exchangeable")
summary(m_gee_e)


m_gee_a <- geeglm(y ~ x*t, data=data, id=id, corstr="ar1")
summary(m_gee_a)


### Now doing a loop
Nsims<-1000
ests.lm<-ests.gee<-ests.gee.e<-ests.gee.a<-matrix(NA,nrow=Nsims,ncol=4)
ses.lm<-ses.gee<-ses.gee.e<-ses.gee.a<-matrix(NA,nrow=Nsims,ncol=4)
for (i in 1:Nsims){
  n <- 50
  m <- 4
  beta <- c(70, 5, 7)
  sigma <- sqrt(10)
  D <- sqrt(15)
  
  x <- rep(c(0,1), each=n/2)
  b0 <- rnorm(n, 0, D)
  
  data <- data.frame(id=rep(1:n, each=m),
                     x=rep(x, each=m),
                     t=rep((1:m)-1, n),
                     b0=rep(b0, each=m),
                     e=rnorm(n*m, 0, sigma))
  
  data <- within(data, {
    y = beta[1] + b0 + beta[2]*t*(x==0) + beta[3]*t*(x==1) + e
  })
  
  
  m_lm <- lm(y ~ factor(x)*t, data=data)
  ests.lm[i,]<-summary(m_lm)$coeff[,1]
  ses.lm[i,]<-summary(m_lm)$coeff[,2]
  
  m_gee <- geeglm(y ~ factor(x)*t, data=data, id=id, corstr="independence")
  ests.gee[i,]<-summary(m_gee)$coeff[,1]
  ses.gee[i,]<-summary(m_gee)$coeff[,2]
  
  m_gee_e <- geeglm(y ~ factor(x)*t, data=data, id=id, corstr="exchangeable")
  ests.gee.e[i,]<-summary(m_gee_e)$coeff[,1]
  ses.gee.e[i,]<-summary(m_gee_e)$coeff[,2]
  
  m_gee_a <- geeglm(y ~ factor(x)*t, data=data, id=id, corstr="ar1")
  ests.gee.a[i,]<-summary(m_gee_a)$coeff[,1]
  ses.gee.a[i,]<-summary(m_gee_a)$coeff[,2]
  
}

mean(ests.lm[,3] - 1.96*ses.lm[,3] < 5 & ests.lm[,3] + 1.96*ses.lm[,3] >5 )
mean(ests.gee[,3] - 1.96*ses.gee[,3] < beta[2] & ests.gee[,3] + 1.96*ses.gee[,3] >beta[2])
mean(ests.gee.e[,3] - 1.96*ses.gee.e[,3] < beta[2] & ests.gee.e[,3] + 1.96*ses.gee.e[,3] >beta[2])
mean(ests.gee.a[,3] - 1.96*ses.gee.a[,3] < beta[2] & ests.gee.a[,3] + 1.96*ses.gee.a[,3] >beta[2])

mean(ests.lm[,4] - 1.96*ses.lm[,4] < 2 & ests.lm[,4] + 1.96*ses.lm[,4] >2)
mean(ests.gee[,4] - 1.96*ses.gee[,4] < 2 & ests.gee[,4] + 1.96*ses.gee[,4] >2)
mean(ests.gee.e[,4] - 1.96*ses.gee.e[,4] < 2 & ests.gee.e[,4] + 1.96*ses.gee.e[,4] >2)
mean(ests.gee.a[,4] - 1.96*ses.gee.a[,4] < 2 & ests.gee.a[,4] + 1.96*ses.gee.a[,4] >2)



