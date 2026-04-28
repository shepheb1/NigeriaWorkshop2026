rm(list=ls())
library(ggplot2)
#detach("package:rms")
#detach("package:Hmisc")

d.cd4 <- read.csv("C:/Users/rukay/Downloads/WORKSHOP 2026/d.cd4.csv",
                  stringsAsFactors = FALSE)
ord<-order(d.cd4$patient_id, d.cd4$days.art)
d.cd4<-d.cd4[ord,]
head(d.cd4)


d.vl <- read.csv("C:/Users/rukay/Downloads/WORKSHOP 2026/d.vl.csv",
                 stringsAsFactors = FALSE)
ord1<-order(d.vl$patient_id, d.vl$days.art)
d.vl<-d.vl[ord1,]
head(d.vl)

library(dplyr)
library(purrr)


#### Number of observations per person
n_per_patient_cd4 <- table(d.cd4$patient_id)
n_per_patient_cd4

d.cd4 |>
  group_by(patient_id) |>
  dplyr::summarize(m=sum(!is.na(cd4_v))) |>
  dplyr::summarize(min(m), median(m), max(m))


n_per_patient_vl <- table(d.vl$patient_id)
n_per_patient_vl
d.vl |>
  group_by(patient_id) |>
  dplyr::summarize(m=sum(!is.na(detect_vl))) |>
  dplyr::summarize(min(m), median(m), max(m))


### Computing slopes and intercepts for CD4 count
slopes_df <- d.cd4 %>%
  group_by(patient_id) %>%
  dplyr::summarize(
    intercept = coef(lm(cd4_v ~ days.art))[1],
    slope = coef(lm(cd4_v ~ days.art))[2]
  )
head(slopes_df)

### Creating a flat dataset that only has one observation per person
d.cd4.flat<-d.cd4[!duplicated(d.cd4$patient_id),]
d.flat<-merge(d.cd4.flat,slopes_df,by="patient_id")
head(d.flat)

summary(d.flat$slope)
summary(d.flat$intercept)
table(d.flat$country)

### Doing analyses for Brazil
### Seeing if the CD4 slopes depend on baseline predictors
d.b<-d.flat[d.flat$country=="Brazil",]
summary(d.b$slope)
boxplot(d.b$slope)
mod.slopes<-lm(slope ~ male_y + age + class + cd4.baseline, data=d.b)
summary(mod.slopes)

### There are some extreme outliers in the estimated slopes, so I am going to switch to a rank-based analysis
library(rms)
mod.slopes1<-orm(slope ~ male_y + age + class + cd4.baseline, data=d.b)
mod.slopes1

### Would I get similar results if I excluded the extreme outliers?
mod.slopes2<-lm(slope ~ male_y + age + class + cd4.baseline, data=d.b, subset=abs(slope)<2)
summary(mod.slopes2)
### Significance depends a little on what I define as an outlier

### Do the intercepts depend on baseline predictors?
summary(d.b$intercept)
boxplot(d.b$intercept)
mod.intcpt<-lm(intercept ~ male_y + age + class + cd4.baseline, data=d.b)
summary(mod.intcpt)

### Again, there are some extreme outliers in the estimated intercepts, so I am going to switch to a rank-based analysis
mod.intcpt1<-orm(intercept ~ male_y + age + class , data=d.b)
mod.intcpt1

### Would I get similar results if I excluded the extreme outliers?
mod.intcpt2<-lm(intercept ~ male_y + age + class, data=d.b, subset=(intercept>0&intercept<2000))
summary(mod.intcpt2)


## Boxplot of CD4 slopes for all countries
boxplot(slope ~ country,
        data = d.flat,
        main = "CD4 Slopes by Country",
        xlab = "Country",
        ylab = "CD4 slope",
        las = 2)

## Boxplot of CD4 intercepts for all countries
boxplot(intercept ~ country,
        data = d.flat,
        main = "CD4 Intercepts by Country",
        xlab = "Country",
        ylab = "CD4 intercept",
        las = 2)


library(rms)

countries <- unique(d.flat$country)

for (cty in countries) {
  
  cat("\n====================\n")
  cat("Country:", cty, "\n")
  cat("====================\n")
  
  d.b <- d.flat[d.flat$country == cty, ]
  
  ## CD4 slope
  print(summary(d.b$slope))
  boxplot(d.b$slope, main = paste("Slope -", cty))
  
  print(summary(lm(slope ~ male_y + age + class + cd4.baseline, data = d.b)))
  print(orm(slope ~ male_y + age + class + cd4.baseline, data = d.b))
  print(summary(lm(slope ~ male_y + age + class + cd4.baseline,
                   data = d.b, subset = abs(slope) < 2)))
  
  ## CD4 intercept
  print(summary(d.b$intercept))
  boxplot(d.b$intercept, main = paste("Intercept -", cty))
  
  print(summary(lm(intercept ~ male_y + age + class + cd4.baseline, data = d.b)))
  print(orm(intercept ~ male_y + age + class + cd4.baseline, data = d.b))
  print(summary(lm(intercept ~ male_y + age + class + cd4.baseline,
                   data = d.b, subset = intercept > 0 & intercept < 2000)))
}