### Longitudinal Data Analysis with R Statistical Software
### R simulation studies
### Benjamin French and Bryan Shepherd
### Vanderbilt-Nigeria Biostatistics Workshop
### April 27 -- May 1, 2026


###################################################
## Day 4
###################################################

rm(list=ls())
library(ggplot2)
library(geepack)
library(lme4)
library(mice)
library(tidyr)
library(dplyr)

set.seed(20260430)

n <- 10000
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

m_lme <- lmer(y ~ (1 | id) + x*t, data=data)
summary(m_lme)

m_gee_e <- geeglm(y ~ x*t, data=data, id=id, corstr="exchangeable")
summary(m_gee_e)

m_gee_i <- geeglm(y ~ x*t, data=data, id=id, corstr="independence")
summary(m_gee_i)

cbind(summary(m_lme)$coeff[,1],
      summary(m_gee_e)$coeff[,1],
      summary(m_gee_i)$coeff[,1])

cbind(summary(m_lme)$coeff[,2],
      summary(m_gee_e)$coeff[,2],
      summary(m_gee_i)$coeff[,2])


# Creating MCAR data
data.mcar<-data
missing<-rbinom(n,1,0.25)
data.mcar<-data[missing==0,]

m_lme1 <- lmer(y ~ (1 | id) + x*t, data=data.mcar)
summary(m_lme1)

m_gee_e1 <- geeglm(y ~ x*t, data=data.mcar, id=id, 
                   corstr="exchangeable")
summary(m_gee_e1)

m_gee_i1 <- geeglm(y ~ x*t, data=data.mcar, id=id, 
                   corstr="independence")
summary(m_gee_i1)

cbind(summary(m_lme1)$coeff[,1],
      summary(m_gee_e1)$coeff[,1],
      summary(m_gee_i1)$coeff[,1])

cbind(summary(m_lme1)$coeff[,2],
      summary(m_gee_e1)$coeff[,2],
      summary(m_gee_i1)$coeff[,2])


# Creating MAR data (#2)
data.mar <- data 

missing<-rbinom(n,1,2*data.mar$t/10)
data.mar<-data[missing==0,]

m_lme2 <- lmer(y ~ (1 | id) + x*t, data=data.mar)
summary(m_lme2)

m_gee_e2 <- geeglm(y ~ x*t, data=data.mar, id=id, 
                   corstr="exchangeable")
summary(m_gee_e2)

m_gee_i2 <- geeglm(y ~ x*t, data=data.mar, id=id, 
                   corstr="independence")
summary(m_gee_i2)

cbind(summary(m_lme2)$coeff[,1],
      summary(m_gee_e2)$coeff[,1],
      summary(m_gee_i2)$coeff[,1])

cbind(summary(m_lme2)$coeff[,2],
      summary(m_gee_e2)$coeff[,2],
      summary(m_gee_i2)$coeff[,2])


# Creating MAR data (#3)
data.mar <- data 

data.mar <- data %>%
  group_by(id) %>%
  arrange(t) %>%
  mutate(
    y_lag = lag(y),
    # higher previous outcome → higher dropout probability
    p_drop = plogis(-8 + 0.1 * y_lag),
    drop = rbinom(n(), 1, ifelse(is.na(p_drop), 0, p_drop)),
    cumdrop = cumsum(drop)
  ) %>%
  ungroup() %>%
  mutate(
    y1 = ifelse(cumdrop>0,NA,y)
  )
#  filter(cumdrop == 0)   # once dropped, remove all future observations

# data.mar$y1<-with(data.mar, ifelse(cumdrop>0,NA,y))

data.mar <- data.frame(data.mar[order(data.mar$id, data.mar$t),])


m_lme3 <- lmer(y1 ~ (1 | id) + x*t, data=data.mar)
summary(m_lme3)

m_gee_e3 <- geeglm(y1 ~ x*t, data=data.mar, id=id, corstr="exchangeable")
summary(m_gee_e3)

m_gee_i3 <- geeglm(y1 ~ x*t, data=data.mar, id=id, corstr="independence")
summary(m_gee_i3)

cbind(summary(m_lme3)$coeff[,1],
      summary(m_gee_e3)$coeff[,1],
      summary(m_gee_i3)$coeff[,1])

cbind(summary(m_lme3)$coeff[,2],
      summary(m_gee_e3)$coeff[,2],
      summary(m_gee_i3)$coeff[,2])


#### Now using multiple imputation

#### Long format imputation:
data.mar1<-data.frame(data.mar[,c("id","x","t","y1")])
data.imputed <- mice(data.mar1, m=5, method="pmm", seed=10)

fit_gee_i <- with(data.imputed, geeglm(y1 ~ x*t, id=id, corstr="independence"))

pooled_results <- pool(fit_gee_i)

pooled_results


data.mar2 <- reshape(data.mar1, idvar="id", timevar="t", v.names="y1",
                     direction="wide")

data.imputed2 <- mice(data.mar2, m=5, method="pmm", seed=10)

imp_long <- complete(data.imputed2, action="long", include=TRUE)

imp_long_reshaped <- imp_long %>%
  pivot_longer(
    cols = starts_with("y1."), 
    names_to = "t", 
    names_prefix = "y1.",
    values_to = "y"
  ) %>%
  mutate(time = as.numeric(t)) %>%
  select(-.id) %>%
  arrange(.imp, id, t) %>%
  as.data.frame()

rownames(imp_long_reshaped) <- seq_len(nrow(imp_long_reshaped))
# Convert back to a mids object so mice functions recognize it

imp_mids_final <- as.mids(imp_long_reshaped, .imp = ".imp")

fit_gee_i2 <- with(imp_mids_final, geeglm(y ~ x*time, id=id, corstr="independence"))

summary(pool(fit_gee_i2))

