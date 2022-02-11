# R code for "A multi-pollutant approach to estimate the causal effects of air pollution mixtures on overall mortality in a large prospective cohort of Dutch individuals"

Packages <- c("car", "mvGPS", "knitr", "cobalt", "WeightIt", "plsRglm", "beepr", "haven", "dplyr", "survival", "tidyverse", "data.table", "readxl", "bkmr", "qgraph", "gWQS", "qgcomp", "corrplot", "cluster","factoextra","gridExtra","table1","stargazer","gbm","fastDummies")
lapply(Packages, library, character.only = TRUE)

lifework <- readRDS("lifework.rds")

# Exposures correlation

names(lifework)
exposure<-lifework[,10:14]
cor(exposure, method = "spearman")

# Table 1. characteristics by clusters

set.seed(123)
#fviz_nbclust(exposure, kmeans, method = "wss")
#2-4 clusters suggested
k2 <- kmeans(exposure, centers = 2, nstart = 20)
k3 <- kmeans(exposure, centers = 3, nstart = 20)
k4 <- kmeans(exposure, centers = 4, nstart = 20)
k5 <- kmeans(exposure, centers = 5, nstart = 20)
p1 <- fviz_cluster(k2, geom = "point", data = exposure) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = exposure) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = exposure) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = exposure) + ggtitle("k = 5")
grid.arrange(p1, p2, p3, p4, nrow = 2) #3 seems like a good classifier
datacluster <- data.frame(lifework, k3$cluster)
table1(~ cohort + age + sex + education + smoking  + BMI + CVDdiagnosis + COPDdiagnosis + CAdiagnosis +
         Income + NDVI + no2 + pm25 + pm25abs + pm10 + op_dtt | k3.cluster, data=datacluster)

# Table 2. regression modeling

## A. Univariate models

lm1 <- glm(vitalstatus~no2+age+sex+smoking+BMI+CVDdiagnosis,data=lifework,family="binomial")
lm2 <- glm(vitalstatus~pm25+age+sex+smoking+BMI+CVDdiagnosis,data=lifework,family="binomial")
lm3 <- glm(vitalstatus~pm25abs+age+sex+smoking+BMI+CVDdiagnosis,data=lifework,family="binomial")
lm4 <- glm(vitalstatus~pm10+age+sex+smoking+BMI+CVDdiagnosis,data=lifework,family="binomial")
lm5 <- glm(vitalstatus~op_dtt+age+sex+smoking+BMI+CVDdiagnosis,data=lifework,family="binomial")

## B. Minimally adjusted multivariate

fit_log_1 <- glm(vitalstatus~no2+pm25+pm25abs+pm10+op_dtt+age+sex+smoking+BMI+CVDdiagnosis,data=lifework,family="binomial")
exp(cbind(OR=coef(fit_log_1), confint(fit_log_1)))
as.data.frame(vif(fit_log_1))

## C. Fully adjusted multivariate

fit_log_2 <- glm(vitalstatus~no2+pm25+pm25abs+pm10+op_dtt+age+sex+smoking+BMI+education+Income+NDVI+CVDdiagnosis+COPDdiagnosis+CAdiagnosis,data=lifework,family="binomial")
exp(cbind(OR=coef(fit_log_2), confint(fit_log_2)))
as.data.frame(vif(fit_log_2))

### Cox

table1(~ time | vitalstatus, data=lifework)
su_obj <- Surv(lifework$time, lifework$vitalstatus)
str(su_obj)
fit_km <- survfit(su_obj ~ 1, data = lifework)
print(fit_km, print.rmean = TRUE)
plot(survfit(su_obj ~ 1, data = lifework), 
     xlab = "Days", 
     ylab = "Overall survival probability")
mall_adj <- coxph(su_obj ~ no2+pm25+pm25abs+pm10+op_dtt+age+sex+smoking+BMI+CVDdiagnosis, data = lifework)
summary(mall_adj) #check non linearities

# Supplementary figure 1. relative importance (weights) from WQS

# define the names of the covariates in the mixture
exposure<- names(lifework[,10:14])
results1_adj <- gwqs(vitalstatus ~ wqs+age+sex+smoking+BMI+CVDdiagnosis, mix_name = exposure, data = lifework, q = 10, validation = 0.6,  
                     b = 100, b1_pos = T, b1_constr = F, family = "binomial", 
                     seed = 123)
results1_adj$final_weights
summary(results1_adj$fit)

## Supplemetary figure 2. H statistics (interaction importance) from random forest

set.seed(123)
LogLossBinary = function(actual, predicted, eps = 1e-15) {  
  predicted = pmin(pmax(predicted, eps), 1-eps)  
  - (sum(actual * log(predicted) + (1 - actual) * log(1 - predicted))) / length(actual)
}
smp_siz = floor(0.75*nrow(lifework)) 
train_ind = sample(seq_len(nrow(lifework)),size = smp_siz)
train =lifework[train_ind,]
test =lifework[-train_ind,]
table(train$Status_C)
table(test$Status_C)
gbm_train<-train[,c(10:15)]
# train GBM model
gbm.fit.final <- gbm(
  formula = vitalstatus ~ .,
  distribution = "bernoulli",
  data = gbm_train,
  n.trees = 15,
  interaction.depth = 3,
  shrinkage = 0.05,
  n.minobsinnode = 100,
  train.fraction = 1,
  n.cores = NULL, # will use all cores by default
  verbose = FALSE
)  
## look at variables
par(mfrow=c(1,1))
summary(
  gbm.fit.final, 
  cBars = 10,
  method = relative.influence, # also can use permutation.test.gbm
  las = 2
)
a<-matrix(ncol=5,nrow=5)
for(i in 1:5) {
  for(j in 1:5) {
    a[i,j]<-interact.gbm(gbm.fit.final, gbm_train,c(i,j))
  }
}

# Table 3. Propensity score

### single exposure propensity score using IPW, trim at 0.99 / 0.97 / 0.95 replicate for each exposure variable

W <- weightit(no2 ~ age+sex+smoking+BMI+CVDdiagnosis, 
              data = lifework, method = "ps", 
              estimand = "ATE")
w <- W$weights
w <- trim(w, at = 0.99, lower = TRUE)

lm1 <- glm(vitalstatus~no2,data=lifework, weights = w, family="binomial")

### mvGPS (excluding PM25abs), trim at 0.99 / 0.97 / 0.95
## mvGPS was delevoped by Justin Williams (https://github.com/williazo) to estimate causal effects for multivariate continuous exposures. Weights are formed assuming a multivariate normal distribution for the simultaneous exposures. The mvGPS was designed for continuous confounders and it easily accomodates binary and categorical confounders as long as they are included as dummy variables in the design matrix.
## mvGPS version 1.2.1 (2021-10-16) was used in the current study.

airpollutants <- c("no2", "pm25", "pm10", "op_dtt")
D <- lifework[airpollutants]
confounders <- c("age","sex_1","bsmoking_1","bsmoking_2","BMI","CVDdiagnosis_0")
C <- lifework[confounders]
out_mvGPS <- mvGPS(D=D, C=C, common = T, trim_w = T, trim_quantile = 0.99)
summary(out_mvGPS)
w <- out_mvGPS$w
dt <- data.frame(lifework$vitalstatus, D)
mvGPS_mod <- glm(dt$lifework ~ dt$no2+dt$pm25+dt$pm10+dt$op_dtt, data=dt, weights=w, family="binomial")
unadj_mod <- glm(dt$lifework ~ dt$no2+dt$pm25+dt$pm10+dt$op_dtt, data=dt, family="binomial")
bias_tbl <- cbind(unadj_OR=exp(unadj_mod$coefficients), unadj_95CI=exp(confint(unadj_mod)), mvGPS_OR=exp(mvGPS_mod$coefficients), mvGPS_95CI=exp(confint(mvGPS_mod)))
bias_tbl
hist(w, breaks = 100)
bal_results <- bal(model_list=c("mvGPS","PS"), D=D, C=C, common = TRUE, trim_w = TRUE, trim_quantile=0.99)
bal_summary <- bal_results$bal_metrics 
#contains overall summary statistics with respect to balance
bal_summary <-data.frame(bal_summary, ESS=c(bal_results$ess, nrow(D)))
#adding in ESS with last value representing the unweighted case
bal_summary <- bal_summary[order(bal_summary$max_cor), ]
bal_summary


# To integrate binary/categorical confounders in the mvGPS function, an alternative option that might provide more stable estimates would involve restoring the intercept into the denominator of the (full conditional) generalized propensity score (lines 133 and 148 https://github.com/williazo/mvGPS/blob/babf41bf51672b722e5b987a83d624ffcc6d91c8/R/mvGPS.R) as follows:

mvGPS <- function(D, C, common=FALSE, trim_w=FALSE, trim_quantile=0.99){
    check_result <- D_C_check(D, C, common)
    assign("D", check_result$D)
    assign("C", check_result$C)

    m <- ncol(D)
    
    for(i in seq_len(m)){
        if(i==1){
            #marginal densities factorized
            d_1 <- lm(D[, i] ~ 1)
            d_1_mu <- coef(d_1)
            d_1_sigma <- summary(d_1)$sigma
            f_d_1 <- dnorm(D[, i], mean=d_1_mu, sd=d_1_sigma)

            #generalized propensity score
            gps_d_1 <- lm(D[, i] ~ C[[i]] - 1) # change into gps_d_1 <- lm(D[, i] ~ C[[i]])
            gps_1_beta <- coef(gps_d_1)
            gps_1_Xb <- model.matrix(gps_d_1) %*% gps_1_beta
            gps_1_sigma <- summary(gps_d_1)$sigma
            f_gps_1 <- dnorm(D[, i], mean=gps_1_Xb, sd=gps_1_sigma)
        } else {
            cond_dens <- lapply(seq_len(m - 1) + 1, function(x){
                #full conditional marginal densities
                d_x <- lm(D[, x] ~ D[, seq_len(x-1)])
                d_x_beta <- coef(d_x)
                d_x_Xb <- model.matrix(d_x) %*% d_x_beta
                d_x_sigma <- summary(d_x)$sigma
                f_d_x <- dnorm(D[, x], mean=d_x_Xb, sd=d_x_sigma)

                #full conditional generalized propensity scores
                gps_x <- lm(D[, x] ~ D[, seq_len(x-1)] + C[[x]] - 1) # change into gps_x <- lm(D[, x] ~ D[, seq_len(x-1)] + C[[x]])
                gps_x_beta <- coef(gps_x)
                gps_x_Xb <- model.matrix(gps_x) %*% gps_x_beta
                gps_x_sigma <- summary(gps_x)$sigma
                f_gps_x <- dnorm(D[, x], gps_x_Xb, gps_x_sigma)

                return(list(marg=f_d_x, gps=f_gps_x))
            })
        }
    }
    cond_results <- unlist(cond_dens, recursive=FALSE)
    num_args <- cond_results[which(names(cond_results)=="marg")]
    num_args[["marg_1"]] <- f_d_1
    denom_args <- cond_results[which(names(cond_results)=="gps")]
    denom_args[["gps_1"]] <- f_gps_1
    
    score <- Reduce("*", denom_args)
    w <- Reduce("*", num_args)/score
    if(trim_w==TRUE){
        #trimming the large weights
        w <- ifelse(w<quantile(w, trim_quantile), w, quantile(w, trim_quantile))
    }
    return(list(score=score, wts=w))
}

