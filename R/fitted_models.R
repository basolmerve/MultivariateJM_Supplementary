library(readr)
library(dplyr)
library(magrittr)
library(JMbayes2)
library(splines)
library(pROC)
library(dynpred) # AUCw()


### Read Data
# Note that the dataset is accessible upon request to senior author.
# Please see "Data Availability" section in the published paper.
pd_data_long <- read_delim(file = "data/longitudinal_data84.txt", delim = "\t", col_names = TRUE, 
                           locale = locale(decimal_mark = "."))

pd_data_additional <- read_delim("data/additional_data.csv", delim = ",", locale = locale(decimal_mark = "."))

pd_data <- pd_data_long %>%
  full_join(pd_data_additional) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(preHD = factor(preHD, levels = c(0, 1), labels = c("No", "Yes")),
         highlow = factor(highlow, levels = c(1, 2), labels = c("Low", "High"))) %>%
  as.data.frame


# Two patients, with IDs 4 and 13, are randomly selected as test patients. These patients
# are excluded from complete data.
# Also, removed patients with missing data.
pd_data_train <- pd_data %>%
  filter(complete.cases(time, event, gender, age, BMI, preHD, ill_count, peritonitrate, ALB, 
                        BUN, KR, CA, P)) %>% 
  filter(!(id %in% c(13, 4)))

# Counting process
counting <- function(times, survtime){
  c(times[-1], survtime)
}

pd_data_train <- pd_data_train %>%
  group_by(id) %>%
  mutate(tstart = time, tstop = counting(time, unique(surv)),
         statusC = if_else(
           event == 0, 
           rep(0, length(time)),
           rep(c(0, 1), c(length(time) - 1, 1))
         ))


###### SURVIVAL DATA (baseline measurements only) #####
pd_data_train_baseline <- pd_data_train %>%
  filter(time == "0")

#### TEST SET
# test1: longitudinal data for the patient with id = 4
# test1_baseline: baseline data for the patient with id = 4
# test2: longitudinal data for the patient with id = 13
# test2_baseline: baseline data for the patient with id = 13

test1 <- pd_data[pd_data$id == 4, ]
test1$event <- 0  # JMbayes2 package assumes all patients are alive at the baseline.
test1_baseline <- test1 %>%
  filter(time == "0")

test2 <- pd_data[pd_data$id == 13, ]
test2_baseline <- test2 %>%
  filter(time == "0")

TestSamples <- rbind(test1, test2)

#### M1: BASELINE ----
### Cox PH model for baseline albumin levels
alb_base <- coxph(Surv(surv, event) ~  age + BMI + preHD + ill_count + 
                    peritonitrate + ALB, data = pd_data_train_baseline)
summary(alb_base)

### Prediction values of model for Test Set (IDs: 4 and 13)
prob_4 <- summary(survfit(alb_base, data = test1_baseline, start.time = 49))
prob_13 <- summary(survfit(alb_base, data = test2_baseline, start.time = 47))

### Time Dependent area under the ROC Curve values 
alb_base_auc <- AUCw(Surv(surv, event) ~  age + BMI + preHD + ill_count + 
                     peritonitrate + ALB, data = pd_data_train_baseline, width = 6)


### Cox PH model for baseline BUN levels
bun_base <- coxph(Surv(surv, event) ~  age + BMI + preHD + ill_count +
                    peritonitrate + BUN, data = pd_data_train_baseline)
summary(bun_base)

#### Prediction values of model for Test Set
prob_bun_4 <- summary(survfit(bun_base, data = test1_baseline, start.time = 49))
prob_bun_13 <- summary(survfit(bun_base, data = test2_baseline, start.time = 47))

### Time Dependent area under the ROC Curve values 
bun_base_auc <- AUCw(Surv(surv, event) ~  age + BMI + preHD + ill_count + 
             peritonitrate + BUN, data = pd_data_train_baseline, width = 6)


###### Cox PH model for baseline Creatinine levels
kr_base <- coxph(Surv(surv, event) ~  age + BMI + preHD + ill_count +
                   peritonitrate + KR, data = pd_data_train_baseline)
summary(kr_base)

#### Prediction values of model for Test Set
prob_kr_114 <- summary( survfit(kr_base, data = test1_baseline, start.time = 49))
prob_kr_13 <- summary( survfit(kr_base, data = test2_baseline, start.time = 47))

### Time Dependent area under the ROC Curve values 
kr_base_auc <- AUCw(Surv(surv, event) ~  age + BMI + preHD + ill_count + 
                    peritonitrate + KR, data = pd_data_train_baseline, width = 6)

###### Cox PH model for baseline Calcium levels
ca_base <- coxph(Surv(surv, event) ~  age + BMI + preHD + ill_count + 
                   peritonitrate + CA, data = pd_data_train_baseline)
summary(ca_base)

#### Prediction values of model for Test Set
prob_ca_4 <- summary( survfit(ca_base, data = test1_baseline, start.time = 48))
prob_ca_13 <- summary( survfit(ca_base, data = test2_baseline, start.time = 47))

### Time Dependent area under the ROC Curve values 
ca_base_auc <- AUCw(Surv(surv, event) ~  age + BMI + preHD + ill_count + 
                     peritonitrate + CA, data = pd_data_train_baseline, width = 6)

###### Cox PH model for baseline Phosporus levels
p_base <- coxph(Surv(surv, event) ~  age + BMI + preHD + ill_count + 
                  peritonitrate + P, data = pd_data_train_baseline)
summary(p_base)

### Prediction values of model for Test Set
prob_p_4 <- summary( survfit(p_base, data = test1_baseline, start.time = 48))
prob_p_13 <- summary( survfit(p_base, data = test2_baseline, start.time = 47))

### Time Dependent area under the ROC Curve values 
p_base_auc <- AUCw(Surv(surv, event) ~  age + BMI + preHD + ill_count + 
                   peritonitrate + P, data = pd_data_train_baseline, width=6)


#### M2: AVERAGED ----
pd_data_M2 <- pd_data_train %>%
  filter(!duplicated(id))

# Note: Variable name for averaged biomarker values has suffix "o", e.g., albo, buno, etc.
###### Cox PH model for Average Albumin levels
alb_ave <- pd_data_M2 %>%
  coxph(Surv(surv, event) ~  age + BMI + preHD + ill_count + peritonitrate +
          albo, data = .)
summary(alb_ave)

#### Prediction values of model for Test Set
probave_alb_4 <- summary( survfit(alb_ave, data = test1_baseline, start.time = 48))
probave_alb_13 <- summary( survfit(alb_ave, data = test2_baseline, start.time = 47))

### Time Dependent area under the ROC Curve values 

alb_ave_auc <- AUCw(Surv(surv, event) ~  age + BMI + preHD + ill_count + 
                     peritonitrate + albo, data = pd_data_M2, width = 6)


######### Cox PH model for Average BUN levels
bun_ave <- pd_data_M2 %>%
  coxph(Surv(surv, event) ~  age + BMI + preHD + ill_count + peritonitrate +
          buno, data = .)
summary(bun_ave)

#### Prediction values of model for Test Set
probave_bun_4 <- summary(survfit(bun_ave, data = test1_baseline, start.time = 48))
probave_bun_13 <- summary(survfit(bun_ave, data = test2_baseline, start.time = 47))

### Time Dependent area under the ROC Curve values 
bun_ave_auc <- AUCw(Surv(surv, event) ~  age + BMI + preHD + ill_count + 
                     peritonitrate + buno, data = pd_data_M2, width = 6)

######### Cox PH model for Average Creatinine levels
kr_ave <- pd_data_M2 %>%
  coxph(Surv(surv, event) ~  age + BMI + preHD + ill_count + peritonitrate +
          kreao, data = .)
summary(kr_ave)

#### Prediction values of model for Test Set
probave_kr_4 <- summary(survfit(kr_ave, data = test1_baseline, start.time = 48))
probave_kr_13 <- summary(survfit(kr_ave, data = test2_baseline, start.time = 47))


### Time Dependent area under the ROC Curve values 
kr_ave_auc <- AUCw(Surv(surv, event) ~  age + BMI + preHD + ill_count + 
                    peritonitrate + kreao, data = pd_data_M2, width = 6)

######### Cox PH model for Average Calcium levels
ca_ave <- pd_data_M2 %>%
  coxph(Surv(surv, event) ~  age + BMI + preHD + ill_count + peritonitrate +
          cao, data = .)
summary(ca_ave)

#### Prediction values of model for Test Set
probave_ca_4 <- summary(survfit(ca_ave, data = test1_baseline, start.time = 48))
probave_ca_13 <- summary(survfit(ca_ave, data = test2_baseline, start.time = 47))

### Time Dependent area under the ROC Curve values 
ca_ave_auc <- AUCw(Surv(surv, event) ~  age + BMI + preHD + ill_count + 
                   peritonitrate + cao, data = pd_data_M2, width = 6)


######### Cox PH model for Average Phosporus levels
p_ave <- pd_data_M2 %>%
  coxph(Surv(surv, event) ~  age + BMI + preHD + ill_count + peritonitrate +
          po, data = .)
summary(p_ave)

#### Prediction values of model for Test Set
probave_p_4 <- summary(survfit(p_ave, data = test1_baseline, start.time = 48))
probave_p_13 <- summary(survfit(p_ave, data = test2_baseline, start.time = 47))

### Time Dependent area under the ROC Curve values 
p_ave_auc <- AUCw(Surv(surv, event) ~  age + BMI + preHD + ill_count + 
                   peritonitrate + po, data = pd_data_M2, width = 6)


#### M3: TIME DEPENDENT COX PH ----
###### Cox PH model for Time Dependent Albumin levels
tdcox_ALB <- coxph(Surv(tstart, tstop, statusC) ~  age + BMI + preHD + ill_count
                   + peritonitrate + ALB, data = pd_data_train)
summary(tdcox_ALB)

#### Prediction values of model for Test Set
prob_tdALB_4 <- summary(survfit(tdcox_ALB, data = test1, start.time = 48 ))
prob_tdALB_13 <- summary(survfit(tdcox_ALB, data = test2, start.time = 47 ))

### Time Dependent area under the ROC Curve values 
alb_td_auc <- AUCw(Surv(tstart, tstop, statusC) ~  age + BMI + preHD + ill_count
                  + peritonitrate + ALB, data = pd_data_train, width = 6)

########## Cox PH model for Time Dependent BUN levels
tdcox_BUN <- coxph(Surv(tstart, tstop, statusC) ~  age + BMI + preHD + ill_count
                   + peritonitrate + BUN, data = pd_data_train)
summary(tdcox_BUN)

#### Prediction values of model for Test Set
prob_tdBUN_4 <- summary(survfit(tdcox_BUN, data = test1, start.time = 48 ))
prob_tdBUN_13 <- summary(survfit(tdcox_BUN, data = test2, start.time = 47 ))

### Time Dependent area under the ROC Curve values 
bun_td_auc <- AUCw(Surv(tstart, tstop, statusC) ~  age + BMI + preHD + ill_count
                   + peritonitrate + BUN, data = pd_data_train, width = 6)


###### Cox PH model for Time Dependent Creatinine levels
tdcox_KR <- coxph(Surv(tstart, tstop, statusC) ~  age + BMI + preHD + ill_count
                  + peritonitrate + KR, data = pd_data_train)
summary(tdcox_KR)

#### Prediction values of model for Test Set
prob_tdKR_4 <- summary(survfit(tdcox_KR, data = test1, start.time = 48))
prob_tdKR_13 <- summary(survfit(tdcox_KR, data = test2, start.time = 47))

### Time Dependent area under the ROC Curve values 
kr_td_auc <- AUCw(Surv(tstart, tstop, statusC) ~  age + BMI + preHD + ill_count
                 + peritonitrate + KR, data = pd_data_train, width = 6)


###### Cox PH model for Time Dependent Calcium levels
tdcox_CA <- coxph(Surv(tstart, tstop, statusC) ~  age + BMI + preHD + ill_count
                  + peritonitrate + CA, data = pd_data_train)
summary(tdcox_CA)

#### Prediction values of model for Test Set
prob_tdCA_4 <- summary(survfit(tdcox_CA, data = test1, start.time = 48))
prob_tdCA_13 <- summary(survfit(tdcox_CA, data = test2, start.time = 47))

### Time Dependent area under the ROC Curve values 
ca_td_auc <- AUCw(Surv(tstart, tstop, statusC) ~  age + BMI + preHD + ill_count
                 + peritonitrate + CA, data = pd_data_train, width = 6)

###### Cox PH model for Time Dependent Phosporus levels
tdcox_P <- coxph(Surv(tstart, tstop, statusC) ~  age + BMI + preHD + ill_count
                 + peritonitrate + P, data = pd_data_train)
summary(tdcox_P)

#### Prediction values of model for Test Set
prob_tdP_4 <- summary(survfit(tdcox_P, data = test1, start.time = 48))
prob_tdP_13 <- summary(survfit(tdcox_P, data = test2, start.time = 47))

### Time Dependent area under the ROC Curve values 
p_td_auc <- AUCw(Surv(tstart, tstop, statusC) ~  age + BMI + preHD + ill_count
                + peritonitrate + P, data = pd_data_train, width = 6)

#### M4: UNIVARIATE JOIN MODELLING ----
## Longitudinal Part 
model1 <- lme(ALB ~ age + highlow + time, random = ~ time | id, data = pd_data_train)
summary(model1)

model2 <- lme(BUN ~ BMI + age + highlow + ill_count + time, random = ~ time | id, data = pd_data_train)
summary(model2)

model3 <- lme(KR ~ age + ill_count + highlow + preHD + time, random = ~ time | id, data = pd_data_train)
summary(model3)

model4 <- lme(CA ~ highlow + time, random = ~ time | id, data = pd_data_train)
summary(model4)

model5 <- lme(P ~ BMI + age + highlow + time, random = ~ time | id, data = pd_data_train)
summary(model5)

## Survival Part (Cox PH) ##
suv_data_JM <- pd_data_train %>%
  filter(!duplicated(id))

survFit <- coxph(Surv(surv, event) ~  age + BMI + preHD + ill_count + peritonitrate, data = suv_data_JM)

## Univariate Joint Models 
JM_U_ALB <- jm(survFit, model1, time_var = "time")
summary(JM_U_ALB)

JM_U_BUN <- jm(survFit, model2, time_var = "time")
summary(JM_U_BUN)

JM_U_CR <- jm(survFit, model3, time_var = "time")
summary(JM_U_CR)

JM_U_CA <- jm(survFit, model4, time_var = "time")
summary(JM_U_CR)

JM_U_P <- jm(survFit, model5, time_var = "time")
summary(JM_U_P)

#### Prediction values of model for Test Set
prob1_JM_4 <- predict(object = JM_U_ALB, newdata = test1, process = "event",
                          times = c(48, 49, 51, 52, 54, 55, 56, 58, 62, 63, 66,
                                    76, 77, 80, 82,83, 87, 93, 95, 115, 130))
prob1_JM_U_4 <- 1 - prob1_JM_4$pred
   

prob1_JM_23 <- predict(object = JM_U_ALB, newdata = test2, process="event",
                          times = c(49, 51, 52, 54, 55, 56, 58, 62, 63, 66, 
                                    76, 77, 80, 82, 83, 87, 93, 95, 115, 130))
prob1_JM_U_23 <- 1 - prob1_JM_23$pred

# BUN
prob2_JM_4 <- predict(JM_U_BUN, newdata = test1, process = "event",
                      times = c(49, 51, 52, 
                                54, 55, 56, 58, 62, 63, 66, 76, 77, 80, 82,
                                83, 87, 93, 95, 115, 130))
prob2_JM_U_4 <- 1 - prob2_JM_4$pred


prob2_JM_23 <- predict(object = JM_U_BUN, newdata = test2, process = "event",
                       times = c(49, 51, 52, 54, 55, 56, 58, 62, 63, 66, 
                                76, 77, 80, 82, 83, 87, 93, 95, 115, 130))
prob2_JM_U_23 <- 1 - prob2_JM_23$pred


# CR
prob3_JM_4 <- predict(JM_U_CR, newdata = test1, process = "event",
                      times = c(49, 51,52, 54, 55, 56, 58, 62, 63, 66, 76, 77, 80, 
                                82, 83, 87, 93, 95, 115, 130))

prob3_JM_U_4 <- 1 - prob3_JM_4$pred

prob3_JM_23 <- predict(object = JM_U_CR, newdata = test2, process = "event",
                        times=c(49, 51, 52, 54, 55, 56, 58, 62, 63, 66, 76, 77, 80, 
                                82, 83, 87, 93, 95, 115, 130))

prob3_JM_U_23 <- 1 - prob3_JM_23$pred

# CA
prob4_JM_4 <- predict(JM_U_CA, newdata = test1, process = "event",
                      times = c(49, 51,52, 54, 55, 56, 58, 62, 63, 66, 76, 77, 80, 
                                82, 83, 87, 93, 95, 115, 130))

prob4_JM_U_4 <- 1 - prob4_JM_4$pred

prob4_JM_23 <- predict(object = JM_U_CA, newdata = test2, process = "event",
                       times=c(38, 41, 42, 44, 45, 46, 49, 51, 52, 54, 55, 56, 58, 
                               62, 63, 66, 76, 77, 80, 82, 83, 87, 93, 95, 115, 130))
prob4_JM_U_23 <- 1 - prob4_JM_23$pred


# P
prob5_JM_4 <- predict(JM_U_P, newdata = test1, process = "event",
                      times = c(49, 51,52, 54, 55, 56, 58, 62, 63, 66, 76, 77, 80, 
                                82, 83, 87, 93, 95, 115, 130))

prob5_JM_U_4 <- 1 - prob5_JM_4$pred

prob5_JM_23 <- predict(object = JM_U_P, newdata = test2, process="event",
                       times=c(49, 51, 52, 54, 55, 56, 58, 62, 63, 66, 76, 77, 80, 
                               82, 83, 87, 93, 95, 115, 130))
prob5_JM_U_23 <- 1 - prob5_JM_23$pred


### Time Dependent area under the ROC Curve values 
auc_U_ALB <- auc_U_BUN <- auc_U_KR <- auc_U_CA <- auc_U_P <- NULL

for (i in seq(0, 66, 6)) {
  auc_U_ALB[i] <- tvAUC(object = JM_U_ALB, newdata = pd_data_train, Tstart = i, Dt = 6.1, idVar = "id", simulate = F)
  auc_U_BUN[i] <- tvAUC(object = JM_U_BUN, newdata = pd_data_train, Tstart = i, Dt = 6.1, idVar = "id", simulate = F)
  auc_U_KR[i] <- tvAUC(object = JM_U_KR, newdata = pd_data_train, Tstart = i, Dt = 6.1, idVar = "id", simulate = F)
  auc_U_CA[i] <- tvAUC(object = JM_U_CA, newdata = pd_data_train, Tstart = i, Dt = 6.1, idVar = "id", simulate = F)
  auc_U_P[i] <- tvAUC(object = JM_U_P, newdata = pd_data_train, Tstart = i, Dt = 6.1, idVar = "id", simulate = F)
}

#### M5: MULTIVARIATE JOIN MODEL ----
## Longitudinal Part (using fitted models from M4)
LME_list_M5 <- list(model1, model2, model3, model4, model5)

## Survival Part 
survFit <- coxph(Surv(surv, event) ~  age + BMI + preHD + ill_count + peritonitrate, 
                 data = suv_data_JM, model = TRUE)

## Multivariate Joint Model 
JM_Mult <- jm(survFit, LME_list_M5, time_var = "time")
summary(JM_Mult)

#### Prediction values of model for Test SeT
prob6_JM_M_4 <- predict(JM_Mult, newdata = test1, process = "event",
                        times = c(49, 51,52, 54, 55, 56, 58, 62, 63, 66, 76, 77, 80, 
                                  82, 83, 87, 93, 95, 115, 130))
prob6_JM_4 <- 1 - prob6_JM_M_4$pred

prob6_JM_M_23 <- predict(object = JM_Mult,newdata = test2, process = "event",
                         times = c(49, 51, 52, 54, 55, 56, 58, 62, 63, 66, 76, 77, 80, 
                                   82, 83, 87, 93, 95, 115, 130))
prob6_JM_M_23 <- 1 - prob6_JM_M_23$pred


### Time Dependent area under the ROC Curve values 
auc_MJM <- NULL

for (i in seq(0, 66, 6)) {
  auc_MJM[i] <- tvAUC(object = JM_Mult, newdata = pd_data_train, Tstart = i, Dt = 6.1,
                      idVar = "id", simulate = F)
}

