
############################################################################
############################################################################
####################### ICM-MM CNV ANALYSIS ################################
############################################################################
############################################################################


############################# LOAD PACKAGES ################################


library(dplyr)
library(tidyr)
library(readr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(forcats)
library(MASS)


############################# READ IN  DATA ###############################

cnv <- read_csv("~/Documents/LINC/ukbb/ukbb.cnvs.51.ndd.csv") %>%
  filter(!is.na(f.eid)) %>%
  mutate(nCNV = rowSums(.[2:52], na.rm=T)) %>%
  mutate(anyCNV = if_else(nCNV >= 1, 1, 0)) 

multi.cnv <- cnv %>%
  filter(nCNV >1)

withdrawn <- read_csv("~/Documents/LINC/ukbb/withdrawn2024.csv", col_names = F)  

cnv <- cnv %>%
  filter(!f.eid %in% withdrawn$X1)

demo <- read_csv("~/Documents/LINC/ukbb/ukbb.demo.csv") %>%
  select(f.eid, Age, Sex, Townsend, Ethnicity, YOB) 

bmi <- read_csv("~/Documents/LINC/ukbb/anthropometry_participant.csv") %>%
  dplyr::select(eid, p21001_i0)

names(bmi) <-  c("f.eid", "BMI")

key <-  read_csv("~/Downloads/match.ids.linc.all.csv") 

pcs <-  read_table("~/Documents/LINC/ukbb/UKB_15175.cov", col_names=T) %>%
  dplyr::select(app15175, c01, c02, c03, c04, c05, c06, c07, c08, c09, c10) %>%
  inner_join(key)

demo <- demo %>% 
  left_join(pcs, by=c("f.eid" = "linc")) %>%
  filter(!f.eid %in% withdrawn$X1)

comb <- read_csv("~/Documents/LINC/ukbb/edited/linc.ehr.all.july.24.csv")
full <- read_csv("~/Documents/LINC/ukbb/edited/linc.ehr.full.july.24.csv")



sums <-  colSums(cnv[-1], na.rm=T)

any.cnv <- cnv %>%
select(f.eid, anyCNV) 


################################ PHENOTYPES #################################

################################ MULTIMORBIDITY #############################

multi <- comb %>%
  mutate(IntHyp = if_else(InternalisingAny == 1 & Hypertension == 1, 1 ,0),
         IntDys = if_else(InternalisingAny == 1 & Dyslipidemia == 1, 1 ,0),
         IntOb = if_else(InternalisingAny == 1 & Obesity == 1, 1 ,0),
         IntDia = if_else(InternalisingAny == 1 & T2D == 1, 1 ,0),
         IntKid = if_else(InternalisingAny == 1 & ChronicKidney == 1, 1 ,0)) %>%
  select(f.eid, InternalisingAny, CMAny, ICMMM, IntHyp, IntDys, IntOb, IntDia, IntKid)


multi2 <- full %>%
  mutate(IntHyp = if_else(InternalisingAny == 1 & Hypertension == 1, 1 ,0),
         IntDys = if_else(InternalisingAny == 1 & Dyslipidemia == 1, 1 ,0),
         IntOb = if_else(InternalisingAny == 1 & Obesity == 1, 1 ,0),
         IntDia = if_else(InternalisingAny == 1 & T2D == 1, 1 ,0),
         IntKid = if_else(InternalisingAny == 1 & ChronicKidney == 1, 1 ,0)) %>%
  select(f.eid, InternalisingAny, CMAny, ICMMM, IntHyp, IntDys, IntOb, IntDia, IntKid)


############################################################################
############################################################################
############################### ANALYSES ###################################
############################################################################
############################################################################


############################### ANY CNV ####################################

any.cnv.t <- any.cnv %>%
  inner_join(demo) %>%
  inner_join(multi) %>%
  mutate(Sex = if_else(Sex == "Male", 0, 1))

outs1 <- names(multi[2:9])

cond1 <- c("Any internalizing",
          "Any cardiometabolic",          
          "Any ICM-MM",
          "Any internalizing and hypertension",
          "Any internalizing and dyslipidemia",
          "Any internalizing and obesity",
          "Any internalizing and T2D",
          "Any internalizing and CKD"
)



any.cnv.reg <- lapply(outs1, function(x) glm(formula(paste(x, "~ anyCNV + Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05")) , data = any.cnv.t, family = "binomial"))
names(any.cnv.reg) <- outs1
lapply(any.cnv.reg, function(x) summary(x))
pvals <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,4])
betas <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,1])
or <- exp(betas)

confint <- lapply(any.cnv.reg, function(x){exp(confint.default(x)[2,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)

confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"

confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"

plot <- as.data.frame(cbind(cond1, or, confint.down, confint.up)) 
plot$cond.group[1] <- "int"
plot$cond.group[2] <- "cm"
plot$cond.group[3:8] <- "mm"


############################### IND CONDITIONS ####################################


any.cnv.t2 <- any.cnv %>%
  inner_join(demo) %>%
  inner_join(comb) %>%
  mutate(Sex = if_else(Sex == "Male", 0, 1))

outs2 <- names(comb[2:9])

cond2 <- c("Anxiety", "Depression", "Somatic symptom disorder", "Hypertension", "Obesity", "Dyslipidemia","T2D", "CKD")   


any.cnv.reg <- lapply(outs2, function(x) glm(formula(paste(x, "~ anyCNV + Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05")) , data = any.cnv.t2, family = "binomial"))
names(any.cnv.reg) <- outs2
lapply(any.cnv.reg, function(x) summary(x))
pvals <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,4])
betas <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,1])
or <- exp(betas)

confint <- lapply(any.cnv.reg, function(x){exp(confint.default(x)[2,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)

confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"

confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"

plot2 <- as.data.frame(cbind(cond2, or, confint.down, confint.up)) 
plot2$cond.group[1] <- "int"
plot2$cond.group[4:8] <- "cm"
names(plot2) <- c("cond1", "or", "lower", "upper", "cond.group")


plot.all <- rbind(plot, plot2) %>%
  arrange(factor(cond1, levels = c("Anxiety", "Depression", "Somatic symptom disorder", "Any internalizing",
                                  "Hypertension", "Obesity", "Dyslipidemia","T2D", "CKD" , "Any cardiometabolic",
                                  "Any internalizing and hypertension",
                                  "Any internalizing and dyslipidemia",
                                  "Any internalizing and obesity",
                                  "Any internalizing and T2D",
                                  "Any internalizing and CKD", "Any ICM-MM")))



ggplot(plot.all, aes(y=fct_rev(fct_inorder(cond1)), x=or, colour=cond.group, shape= cond.group)) +
  geom_point(size=4) +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=.3) +
  geom_vline(xintercept=1, linetype='longdash') +
  labs(x='Odds ratio', y=element_blank()) + 
  theme_classic() + labs(color  = "cond.group", shape = "cond.group") +
  scale_colour_manual(values= c(   "#990033", "#0099CC", "#006666"), labels=c("Cardiometabolic", "Internalizing","ICM-MM")) + scale_shape(labels=c("Cardiometabolic", "Internalizing","ICM-MM")) 



############################### FIRTH CORRECTION ####################################

install.packages("logistf")
library(logistf)

any.cnv.reg <-  lapply(outs1, function(x)logistf(formula = x ~ anyCNV + Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05), data = any.cnv.t)
summary(any.cnv.reg)

for (i in outs1) {
  list[[i]] <- logistf(formula = i ~ anyCNV + Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05, data = any.cnv.t)
}

reg1 <- logistf(formula = IntKid ~ anyCNV + Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05, data = any.cnv.t)
summary(reg1)

any.cnv.reg <- lapply(outs1, function(x) glm(formula(paste(x, "~ anyCNV + Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05")) , data = any.cnv.t, family = "binomial"))
names(any.cnv.reg) <- outs1
lapply(any.cnv.reg, function(x) summary(x))
pvals <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,4])
betas <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,1])
or <- exp(betas)



############################### FULL EHR DATA ####################################

any.cnv.t.full <- any.cnv %>%
  inner_join(demo) %>%
  inner_join(multi2) %>%
  mutate(Sex = if_else(Sex == "Male", 0, 1))

any.cnv.reg <- lapply(outs1, function(x) glm(formula(paste(x, "~ anyCNV + Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05")) , data = any.cnv.t.full, family = "binomial"))
names(any.cnv.reg) <- outs1
lapply(any.cnv.reg, function(x) summary(x))
pvals <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,4])
betas <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,1])
or <- exp(betas)

confint <- lapply(any.cnv.reg, function(x){exp(confint.default(x)[2,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)

confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"

confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"

plot <- as.data.frame(cbind(cond1, or, confint.down, confint.up))
plot$cond.group[1] <- "int"
plot$cond.group[2] <- "cm"
plot$cond.group[3:8] <- "mm"



ggplot(plot, aes(y=fct_rev(fct_inorder(cond1)), x=or, colour=cond.group, shape= cond.group)) +
  geom_point(size=4) +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=.3) +
  geom_vline(xintercept=1, linetype='longdash') +
  labs(x='Odds ratio', y=element_blank()) + 
  theme_classic() + theme(legend.position="none") +
  scale_colour_manual(values= c(   "#990033", "#0099CC", "#006666"), labels=c("Cardiometabolic", "Internalising","ICM-MM"))


############################### EXCLUDE 16p DEL ####################################

any.cnv.sens <- cnv[ , -c(34,36)] 
any.cnv.sens <-  any.cnv.sens %>%
  filter(!is.na(f.eid)) %>%
  mutate(nCNV = rowSums(.[2:50], na.rm=T)) %>%
  mutate(anyCNV = if_else(nCNV >= 1, 1, 0)) %>%
  select(f.eid, anyCNV)  


any.cnv.t.sens <- any.cnv.sens %>%
  inner_join(demo) %>%
  inner_join(multi) %>%
  mutate(Sex = if_else(Sex == "Male", 0, 1))


any.cnv.reg.sens <- lapply(outs1, function(x) glm(formula(paste(x, "~ anyCNV + Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05")) , data = any.cnv.t.sens, family = "binomial"))
names(any.cnv.reg.sens) <- outs1
lapply(any.cnv.reg, function(x) summary(x))
pvals <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,4])
betas <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,1])
or <- exp(betas)

confint <- lapply(any.cnv.reg, function(x){exp(confint.default(x)[2,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)

confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"

confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"

plot <- as.data.frame(cbind(cond1, or, confint.down, confint.up)) 
plot$cond.group[1] <- "int"
plot$cond.group[2] <- "cm"
plot$cond.group[3:8] <- "mm"


ggplot(plot, aes(y=fct_rev(fct_inorder(cond1)), x=or, colour=cond.group, shape=cond.group)) +
  geom_point(size=4) +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=.3) +
  geom_vline(xintercept=1, linetype='longdash') +
  labs(x='Odds ratio', y=element_blank()) +
  theme_classic() + theme(legend.position="none") +
  scale_colour_manual(values= c(   "#990033", "#0099CC", "#006666"))



########################## ETHNICITY SPECIFIC #################################


multi <- multi %>%
  select(-InternalisingAny, -CMAny)

ethn <- demo %>%
  mutate(ethn = case_when(Ethnicity == "British" | Ethnicity == "Irish" | Ethnicity == "Any other white background" ~ "White",
                          Ethnicity == "White and Black Caribbean" | Ethnicity == "White and Black African" | Ethnicity == "White and Asian" | Ethnicity == "Any other mixed background" ~ "Mixed",
                          Ethnicity == "Indian" | Ethnicity == "Pakistani" | Ethnicity == "Bangladeshi" | Ethnicity == "Any other Asian background" ~ "Asian",
                          Ethnicity == "Caribbean" | Ethnicity == "African" | Ethnicity == "Any other Black background" ~ "Black",
                          Ethnicity == "Chinese" ~ "Chinese",
                          Ethnicity == "Other ethnic group" ~ "Other"))



any.cnv.t <- any.cnv %>%
  inner_join(ethn) %>%
  inner_join(multi) %>%
  mutate(Sex = if_else(Sex == "Male", 0, 1))


white <- any.cnv.t %>%
  filter(ethn == "White")

any.cnv.reg <- lapply(outs1, function(x) glm(formula(paste(x, "~ anyCNV + Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05")) , data = white, family = "binomial"))
names(any.cnv.reg) <- outs1
lapply(any.cnv.reg, function(x) summary(x))
pvals <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,4])
betas <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,1])
or <- exp(betas)

confint <- lapply(any.cnv.reg, function(x){exp(confint.default(x)[2,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)

confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"

confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"

plot <- as.data.frame(cbind(cond, or, confint.down, confint.up)) 


black <- any.cnv.t %>%
  filter(ethn == "Black")

any.cnv.reg <- lapply(outs1, function(x) glm(formula(paste(x, "~ anyCNV + Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05")) , data = black, family = "binomial"))
names(any.cnv.reg) <- outs1
lapply(any.cnv.reg, function(x) summary(x))
pvals <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,4])
betas <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,1])
or <- exp(betas)

confint <- lapply(any.cnv.reg, function(x){exp(confint.default(x)[2,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)

confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"

confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"

plot <- as.data.frame(cbind(cond, or, confint.down, confint.up)) 

asian <- any.cnv.t %>%
  filter(ethn == "Asian" | ethn == "Chinese")

any.cnv.reg <- lapply(outs1, function(x) glm(formula(paste(x, "~ anyCNV + Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05")) , data = asian, family = "binomial"))
names(any.cnv.reg) <- outs1
pvals <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,4])
betas <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,1])
or <- exp(betas)

confint <- lapply(any.cnv.reg, function(x){exp(confint.default(x)[2,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)

confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"

confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"

plot <- as.data.frame(cbind(cond, or, confint.down, confint.up)) 


mixed <- any.cnv.t %>%
  filter(ethn == "Mixed")

any.cnv.reg <- lapply(outs1, function(x) glm(formula(paste(x, "~ anyCNV + Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05")) , data = mixed, family = "binomial"))
names(any.cnv.reg) <- outs1
pvals <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,4])
betas <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,1])
or <- exp(betas)

confint <- lapply(any.cnv.reg, function(x){exp(confint.default(x)[2,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)

confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"

confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"

plot <- as.data.frame(cbind(cond, or, confint.down, confint.up)) 

other <- any.cnv.t %>%
  filter(ethn == "Other")

any.cnv.reg <- lapply(outs1, function(x) glm(formula(paste(x, "~ anyCNV + Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05")) , data = other, family = "binomial"))
names(any.cnv.reg) <- outs1
pvals <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,4])
betas <- sapply(any.cnv.reg, function(x) summary(x)$coefficients[2,1])
or <- exp(betas)

confint <- lapply(any.cnv.reg, function(x){exp(confint.default(x)[2,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)

confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"

confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"

plot <- as.data.frame(cbind(cond, or, confint.down, confint.up)) 



######################### SEX DIFFERENCES ####################################


cnv.reg.sex <- lapply(outs1, function(x) glm(formula(paste(x, "~ Age + anyCNV + Townsend + Sex + Sex*anyCNV + c01 + c02 + c03 + c04 + c05")), family="binomial", data=any.cnv.t)) 
names(cnv.reg.sex) <- outs1
lapply(cnv.reg.sex, function(x){summary(x)})

### interaction for int+hyper, inte+t2d, any cm ###

### PLOT ###

test.cnv.fem <- any.cnv.t %>%
  filter(Sex == "1")

cnv.reg <- lapply(outs1, function(x) glm(formula(paste(x, "~ Age + anyCNV + Townsend")), family="binomial", data=test.cnv.fem)) 
names(cnv.reg) <- outs1
pval <- sapply(cnv.reg, function(x) summary(x)$coefficients[3,4])
beta <- sapply(cnv.reg, function(x) summary(x)$coefficients[3,1])
or <- exp(beta)
confint <- lapply(cnv.reg, function(x){exp(confint.default(x)[3,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"

plot.female <- as.data.frame(cbind(cond1, or, confint.down, confint.up)) %>%
  mutate(Sex = "Female")

test.cnv.m <- any.cnv.t %>%
  filter(Sex == "0")

cnv.reg <- lapply(outs1, function(x) glm(formula(paste(x, "~ Age + anyCNV + Townsend")), family="binomial", data=test.cnv.m)) 
names(cnv.reg) <- outs1
pval <- sapply(cnv.reg, function(x) summary(x)$coefficients[3,4])
beta <- sapply(cnv.reg, function(x) summary(x)$coefficients[3,1])
or <- exp(beta)
confint <- lapply(cnv.reg, function(x){exp(confint.default(x)[3,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"

plot.male <- as.data.frame(cbind(cond1, or, confint.down, confint.up)) %>%
  mutate(Sex = "Male")

plot <- rbind(plot.female, plot.male)

ggplot(plot, aes(y=fct_rev(fct_inorder(cond1)), x=or, colour=Sex, shape= Sex)) +
  geom_point(size=3, position = position_dodge(0.9)) +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=.3, position = position_dodge(0.9)) +
  geom_vline(xintercept=1, linetype='longdash') +
  labs(x='Odds ratio', y=element_blank()) +
  theme_classic() +
  scale_colour_manual(values = c("#668b8b","#8b4789"))


###################### INDIVIDUAL CONDITIONS ###################### 


cnv.reg.sex <- lapply(outs2, function(x) glm(formula(paste(x, "~ Age + anyCNV + Townsend + Sex + Sex*anyCNV + c01 + c02 + c03 + c04 + c05")), family="binomial", data=any.cnv.t2)) 
names(cnv.reg.sex) <- outs2
lapply(cnv.reg.sex, function(x){summary(x)})

### sig interaction for t2d and hypertension ###

### PLOT ###

test.cnv.fem <- any.cnv.t2 %>%
  filter(Sex == "1")

cnv.reg <- lapply(outs2, function(x) glm(formula(paste(x, "~ Age + anyCNV + Townsend")), family="binomial", data=test.cnv.fem)) 
names(cnv.reg) <- outs2
pval <- sapply(cnv.reg, function(x) summary(x)$coefficients[3,4])
beta <- sapply(cnv.reg, function(x) summary(x)$coefficients[3,1])
or <- exp(beta)
confint <- lapply(cnv.reg, function(x){exp(confint.default(x)[3,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"

plot.female <- as.data.frame(cbind(cond2, or, confint.down, confint.up)) %>%
  mutate(Sex = "Female")

test.cnv.m <- any.cnv.t2 %>%
  filter(Sex == "0")

cnv.reg <- lapply(outs2, function(x) glm(formula(paste(x, "~ Age + anyCNV + Townsend")), family="binomial", data=test.cnv.m)) 
names(cnv.reg) <- outs2
pval <- sapply(cnv.reg, function(x) summary(x)$coefficients[3,4])
beta <- sapply(cnv.reg, function(x) summary(x)$coefficients[3,1])
or <- exp(beta)
confint <- lapply(cnv.reg, function(x){exp(confint.default(x)[3,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"

plot.male <- as.data.frame(cbind(cond2, or, confint.down, confint.up)) %>%
  mutate(Sex = "Male")

plot <- rbind(plot.female, plot.male)

ggplot(plot, aes(y=fct_rev(fct_inorder(cond2)), x=or, colour=Sex, shape=Sex)) +
  geom_point(size=3, position = position_dodge(0.9)) +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=.3, position = position_dodge(0.9)) +
  geom_vline(xintercept=1, linetype='longdash') +
  labs(x='Odds ratio', y=element_blank()) +
  theme_classic() +
  scale_colour_manual(values = c("#668b8b","#8b4789"))



####################### DELETIONS VS DUPLICATIONS ##########################


deletions <- cnv %>%
  select(f.eid, ends_with("del")) %>%
  mutate(anyDel = rowSums(.[2:30], na.rm=T)) %>%
  mutate(anyDel = if_else(anyDel >= 1, 1, 0)) %>%
  select(f.eid, anyDel)

duplications <- cnv %>%
  select(f.eid, ends_with("dup")) %>%
  mutate(anyDup = rowSums(.[2:19], na.rm=T)) %>%
  mutate(anyDup = if_else(anyDup >= 1, 1, 0)) %>%
  select(f.eid, anyDup)

cnv.geno <- deletions %>%
  full_join(duplications) 

both <- cnv.geno %>%
  filter(anyDel == 1 & anyDup == 1)

cnv.geno <- cnv.geno %>%
  filter(!f.eid %in% both$f.eid) %>%
  mutate(Geno = case_when(anyDel == 1 ~ "Del",
                          anyDup == 1 ~ "Dup",
                          anyDel == 0 & anyDup == 0 ~ "Control")) %>%
  inner_join(demo) %>%
  inner_join(multi) 


geno.reg <- lapply(outs1, function(x) glm(formula(paste(x, "~ Geno + Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05")) , data = cnv.geno, family = "binomial"))
names(geno.reg) <- outs1
lapply(geno.reg, function(x) summary(x))
pvals <- as.data.frame(sapply(geno.reg, function(x) summary(x)$coefficients[2,4]))
betas <- as.data.frame(sapply(geno.reg, function(x) summary(x)$coefficients[2,1]))


### using the emmeans package to do pairwise comparisons of the marginal means ###

library(emmeans)

marginal <- emmeans(geno.reg[[2]], ~ Geno)

pairs(marginal)

### only sig for int+obesity ###

### PLOT ###

test.cnv.del <- cnv.geno %>%
  filter(!anyDup == 1) 

cnv.reg.del <- lapply(outs1, function(x) glm(formula(paste(x, "~ Age + Sex + anyDel + Townsend")), family="binomial", data=test.cnv.del)) 
names(cnv.reg.del) <- outs1
pval <- sapply(cnv.reg.del, function(x) summary(x)$coefficients[4,4])
beta <- sapply(cnv.reg.del, function(x) summary(x)$coefficients[4,1])
or <- exp(beta)
confint <- lapply(cnv.reg.del, function(x){exp(confint.default(x)[4,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"

plot.del <- as.data.frame(cbind(cond1, or, confint.down, confint.up)) %>%
  mutate(CNV = "Deletions")

test.cnv.dup <- cnv.geno %>%
  filter(!anyDel == 1) 

cnv.reg.dup <- lapply(outs1, function(x) glm(formula(paste(x, "~ Age + Sex + anyDup +  Townsend")), family="binomial", data=test.cnv.dup)) 
names(cnv.reg.dup) <- outs1
pval <- sapply(cnv.reg.dup, function(x) summary(x)$coefficients[4,4])
beta <- sapply(cnv.reg.dup, function(x) summary(x)$coefficients[4,1])
or <- exp(beta)
confint <- lapply(cnv.reg.dup, function(x){exp(confint.default(x)[4,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"
plot.dup <- as.data.frame(cbind(cond1, or, confint.down, confint.up)) %>%
  mutate(CNV = "Duplications")

plot <- rbind(plot.del, plot.dup)

ggplot(plot, aes(y=fct_rev(fct_inorder(cond1)), x=or, colour=CNV, shape=CNV)) +
  geom_point(size=3, position = position_dodge(0.9)) +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=.3, position = position_dodge(0.9)) +
  geom_vline(xintercept=1, linetype='longdash') +
  labs(x='Odds ratio', y=element_blank()) +
  theme_classic() +
  scale_colour_manual(values = c("#6959CD","#336633"))



###################### INDIVIDUAL CONDITIONS ###################### 

cnv.geno2 <- cnv.geno %>%
  filter(!f.eid %in% both$f.eid) %>%
  mutate(Geno = case_when(anyDel == 1 ~ "Del",
                          anyDup == 1 ~ "Dup",
                          anyDel == 0 & anyDup == 0 ~ "Control")) %>%
  inner_join(demo) %>%
  inner_join(comb) 


geno.reg <- lapply(outs2, function(x) glm(formula(paste(x, "~ Geno + Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05")) , data = cnv.geno2, family = "binomial"))
names(geno.reg) <- outs2
lapply(geno.reg, function(x) summary(x))
pvals <- as.data.frame(sapply(geno.reg, function(x) summary(x)$coefficients[2,4]))
betas <- as.data.frame(sapply(geno.reg, function(x) summary(x)$coefficients[2,1]))


### using the emmeans package to do pairwise comparisons of the marginal means ###

library(emmeans)

marginal <- emmeans(geno.reg[[1]], ~ Geno)

pairs(marginal)

###  sig for obesity, dyslipidemia, t2d ###

### PLOT ###

test.cnv.del <- cnv.geno2 %>%
  filter(!anyDup == 1) 

cnv.reg.del <- lapply(outs2, function(x) glm(formula(paste(x, "~ Age + Sex + anyDel + Townsend")), family="binomial", data=test.cnv.del)) 
names(cnv.reg.del) <- outs2
pval <- sapply(cnv.reg.del, function(x) summary(x)$coefficients[4,4])
beta <- sapply(cnv.reg.del, function(x) summary(x)$coefficients[4,1])
or <- exp(beta)
confint <- lapply(cnv.reg.del, function(x){exp(confint.default(x)[4,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"

plot.del <- as.data.frame(cbind(cond2, or, confint.down, confint.up)) %>%
  mutate(CNV = "Deletions")

test.cnv.dup <- cnv.geno2 %>%
  filter(!anyDel == 1) 

cnv.reg.dup <- lapply(outs2, function(x) glm(formula(paste(x, "~ Age + Sex + anyDup +  Townsend")), family="binomial", data=test.cnv.dup)) 
names(cnv.reg.dup) <- outs2
pval <- sapply(cnv.reg.dup, function(x) summary(x)$coefficients[4,4])
beta <- sapply(cnv.reg.dup, function(x) summary(x)$coefficients[4,1])
or <- exp(beta)
confint <- lapply(cnv.reg.dup, function(x){exp(confint.default(x)[4,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"
plot.dup <- as.data.frame(cbind(cond2, or, confint.down, confint.up)) %>%
  mutate(CNV = "Duplications")

plot <- rbind(plot.del, plot.dup)

ggplot(plot, aes(y=fct_rev(fct_inorder(cond2)), x=or, colour=CNV, shape=CNV)) +
  geom_point(size=3, position = position_dodge(0.9)) +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=.3, position = position_dodge(0.9)) +
  geom_vline(xintercept=1, linetype='longdash') +
  labs(x='Odds ratio', y=element_blank()) +
  theme_classic() +
  scale_colour_manual(values = c("#6959CD","#336633"))





############################## GENE DOSAGE EFFECT #################################

################################ GET GENE COUNTS ##################################

cnv.genes <- read_csv("~/Documents/LINC/CNV_genes.csv") %>%
  dplyr::select(f.eid, Genes)

genes.long <- cnv.genes %>%
  separate_longer_delim(Genes, ",")

write.table(file="~/Desktop/genes.csv", genes.long, col.names=T, row.names=F, quote=F, sep=",")
genes.long <- read_csv("~/Documents/LINC/ukbb/edited/genes.csv")

count.genes <- genes.long %>%
  count(f.eid) 

haplo.triplo <- read_csv("~/Documents/LINC/haplo.triplo.genes.collins.et.al.csv") %>%
  select(Gene, pHaplo, pTriplo)

haplo <- haplo.triplo %>%
  filter(pHaplo >= 0.86) %>%
  select(Gene)

triplo <- haplo.triplo %>%
  filter(pTriplo >= 0.94) %>%
  select(Gene)

genes.haplo <- genes.long %>%
  inner_join(haplo, by=c("Genes"="Gene"))


genes.triplo <- genes.long %>%
  inner_join(triplo, by=c("Genes"="Gene"))  

count.haplo <- genes.haplo %>%
  count(f.eid) 
write.table(file="~/Desktop/haplo.csv", count.haplo, col.names=T, row.names=F, quote=F, sep=",")

count.triplo <- genes.triplo %>%
  count(f.eid) 
write.table(file="~/Desktop/triplo.csv", count.triplo, col.names=T, row.names=F, quote=F, sep=",")

count.triplo <- read_csv("~/Documents/LINC/ukbb/edited/triplo.csv")
count.haplo <- read_csv("~/Documents/LINC/ukbb/edited/haplo.csv")


names(count.triplo) <- c("f.eid", "ntriplo")
names(count.haplo) <- c("f.eid", "nhaplo")

counts <- count.genes %>%
  full_join(count.triplo) %>%
  full_join(count.haplo) %>%
  mutate_at(vars(3:4), ~replace(., is.na(.), 0))



counts <- counts %>%
  inner_join(key, by = c("f.eid"="app14421")) %>%
  dplyr::select(linc, n, nhaplo, ntriplo)


counts.cnv <- deletions %>%
  inner_join(duplications) %>%
  left_join(counts, by=c("f.eid"="linc")) %>%
  mutate_at(vars(4:6), ~replace(., is.na(.), 0)) %>%
  mutate(nhaplo = if_else(anyDel == 1, nhaplo, 0), 
         ntriplo = if_else(anyDup == 1, ntriplo, 0))


write.table(file="~/Documents/LINC/ukbb/edited/dosage.sensitive.gene.counts.cnv.june.24.csv", counts.cnv, col.names=T, row.names=F, quote=F, sep=",")

counts.cnv <- read_csv("~/Documents/LINC/ukbb/edited/dosage.sensitive.gene.counts.cnv.june.24.csv")

####################################### TEST #########################################

count.cnv.t.hap <- any.cnv.t %>%
  inner_join(counts.cnv) %>%
  filter(anyDup == 0)

cnv.reg.hap <- lapply(outs1, function(x) glm(formula(paste(x, "~ Age + Sex + nhaplo + n + Townsend + c01 + c02 + c03 + c04 + c05")), family="binomial", data=count.cnv.t.hap)) 
names(cnv.reg.hap) <- outs1
lapply(cnv.reg.hap, function(x){summary(x)})
pval <- sapply(cnv.reg.hap, function(x) summary(x)$coefficients[4,4])
beta <- sapply(cnv.reg.hap, function(x) summary(x)$coefficients[4,1])
or <- exp(beta)
confint <- lapply(cnv.reg.hap, function(x){exp(confint.default(x)[4,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  dplyr::select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  dplyr::select(value)
names(confint.down) <- "lower"

plot.hap <- as.data.frame(cbind(cond1, or, confint.down, confint.up)) %>%
  mutate(Genes = "Haploinsufficient")

count.cnv.t.trip <- any.cnv.t %>%
  inner_join(counts.cnv) %>%
  filter(anyDel == 0)

cnv.reg.trip <- lapply(outs1, function(x) glm(formula(paste(x, "~ Age + Sex + ntriplo + n + Townsend + c01 + c02 + c03 + c04 + c05")), family="binomial", data=count.cnv.t.trip)) 
names(cnv.reg.trip) <- outs1
lapply(cnv.reg.trip, function(x){summary(x)})
pval <- sapply(cnv.reg.trip, function(x) summary(x)$coefficients[4,4])
beta <- sapply(cnv.reg.trip, function(x) summary(x)$coefficients[4,1])
or <- exp(beta)
confint <- lapply(cnv.reg.trip, function(x){exp(confint.default(x)[4,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  dplyr::select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  dplyr::select(value)
names(confint.down) <- "lower"

plot.trip <- as.data.frame(cbind(cond1, or, confint.down, confint.up)) %>%
  mutate(Genes = "Triplosensitive")

plot <- rbind(plot.hap, plot.trip)

ggplot(plot, aes(y=fct_rev(fct_inorder(cond1)), x=or, colour=Genes, shape=Genes)) +
  geom_point(size=3, position = position_dodge(0.9)) +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=.3, position = position_dodge(0.9)) +
  geom_vline(xintercept=1, linetype='longdash') +
  labs(x='Odds ratio', y=element_blank()) + 
  theme_classic() +
  scale_colour_manual(values = c("#6959CD","#336633"))




########################### INDIVIDUAL CNVS #########################################

cnv <- cnv[ , colSums(cnv, na.rm=T) >=5] 

cnv.t <- demo %>%
  inner_join(cnv) %>%
  inner_join(multi)

preds <- names(cnv.t)[20:52]

cnv.reg1 <- lapply(preds, function(x) glm(formula(paste("ICMMM ~", x, "+ Age + Sex + Townsend + c01 + c02 + c03 + c04 + c05")), family="binomial", data=cnv.t))
names(cnv.reg1) <- preds
beta1 <- as.data.frame(sapply(cnv.reg1, function(x) summary(x)$coefficients[2,1]))
colnames(beta1) <- "ICMMM"
p1 <- as.data.frame(sapply(cnv.reg1, function(x) summary(x)$coefficients[2,4]))
colnames(p1) <- "ICMMM"
confint1 <- lapply(cnv.reg1, function(x){exp(confint.default(x)[2,])}) 
confint1 <- as.data.frame(confint1)
confint1$names <- rownames(confint1)
confint1.up <- confint1[2,] %>%
  melt() %>%
  dplyr::select(value)
names(confint1.up) <- "upper"
confint1.down <- confint1[1,] %>%
  melt() %>%
  dplyr::select(value)
names(confint1.down) <- "lower"
odd1 <- exp(beta1)
names(odd1) <- "OR"
names(p1) <- "P"
r1 <- cbind(odd1, confint1.up, confint1.down, p1) %>%
  mutate(Phenotype = "Any ICM-MM")
r1$CNV <- rownames(r1)


cnv.reg2 <- lapply(preds, function(x) glm(formula(paste("IntHyp ~", x, "+ Age + Sex + c01 + c02 + c03 + c04 + c05")), family="binomial", data=cnv.t))
names(cnv.reg2) <- preds
beta2 <- as.data.frame(sapply(cnv.reg2, function(x) summary(x)$coefficients[2,1]))
colnames(beta2) <- "IntHyp"
p2 <- as.data.frame(sapply(cnv.reg2, function(x) summary(x)$coefficients[2,4]))
colnames(p2) <- "IntHyp"
confint2 <- lapply(cnv.reg2, function(x){exp(confint.default(x)[2,])}) 
confint2 <- as.data.frame(confint2)
confint2$names <- rownames(confint2)
confint2.up <- confint2[2,] %>%
  melt() %>%
  dplyr::select(value)
names(confint2.up) <- "upper"
confint2.down <- confint2[1,] %>%
  melt() %>%
  dplyr::select(value)
names(confint2.down) <- "lower"
odd2 <- exp(beta2)
names(odd2) <- "OR"
names(p2) <- "P"
r2 <- cbind(odd2, confint2.up, confint2.down, p2) %>%
  mutate(Phenotype = "Any internalizing and hypertension")
r2$CNV <- rownames(r2)


cnv.reg3 <- lapply(preds, function(x) glm(formula(paste("IntDys ~", x, "+ Age + Sex + c01 + c02 + c03 + c04 + c05")), family="binomial", data=cnv.t))
names(cnv.reg3) <- preds
beta3 <- as.data.frame(sapply(cnv.reg3, function(x) summary(x)$coefficients[2,1]))
colnames(beta3) <- "IntDys"
p3 <- as.data.frame(sapply(cnv.reg3, function(x) summary(x)$coefficients[2,4]))
colnames(p3) <- "IntDys"
confint3 <- lapply(cnv.reg3, function(x){exp(confint.default(x)[2,])}) 
confint3 <- as.data.frame(confint3)
confint3$names <- rownames(confint3)
confint3.up <- confint3[2,] %>%
  melt() %>%
  dplyr::select(value)
names(confint3.up) <- "upper"
confint3.down <- confint3[1,] %>%
  melt() %>%
  dplyr::select(value)
names(confint3.down) <- "lower"
odd3 <- exp(beta3)
names(odd3) <- "OR"
names(p3) <- "P"
r3 <- cbind(odd3, confint3.up, confint3.down, p3) %>%
  mutate(Phenotype = "Any internalizing and dyslipidemia")
r3$CNV <- rownames(r3)


cnv.reg4 <- lapply(preds, function(x) glm(formula(paste("IntOb ~", x, "+ Age + Sex + c01 + c02 + c03 + c04 + c05")), family="binomial", data=cnv.t))
names(cnv.reg4) <- preds
beta4 <- as.data.frame(sapply(cnv.reg4, function(x) summary(x)$coefficients[2,1]))
colnames(beta4) <- "IntOb"
p4 <- as.data.frame(sapply(cnv.reg4, function(x) summary(x)$coefficients[2,4]))
colnames(p4) <- "IntOb"
confint4 <- lapply(cnv.reg4, function(x){exp(confint.default(x)[2,])}) 
confint4 <- as.data.frame(confint4)
confint4$names <- rownames(confint4)
confint4.up <- confint4[2,] %>%
  melt() %>%
  dplyr::select(value)
names(confint4.up) <- "upper"
confint4.down <- confint4[1,] %>%
  melt() %>%
  dplyr::select(value)
names(confint4.down) <- "lower"
odd4 <- exp(beta4)
names(odd4) <- "OR"
names(p4) <- "P"
r4 <- cbind(odd4, confint4.up, confint4.down, p4) %>%
  mutate(Phenotype = "Any internalizing and obesity")
r4$CNV <- rownames(r4)


cnv.reg5 <- lapply(preds, function(x) glm(formula(paste("IntDia ~", x, "+ Age + Sex + c01 + c02 + c03 + c04 + c05")), family="binomial", data=cnv.t))
names(cnv.reg5) <- preds
beta5 <- as.data.frame(sapply(cnv.reg5, function(x) summary(x)$coefficients[2,1]))
colnames(beta5) <- "IntDia"
p5 <- as.data.frame(sapply(cnv.reg5, function(x) summary(x)$coefficients[2,4]))
colnames(p5) <- "IntDia"
confint5 <- lapply(cnv.reg5, function(x){exp(confint.default(x)[2,])}) 
confint5 <- as.data.frame(confint5)
confint5$names <- rownames(confint5)
confint5.up <- confint5[2,] %>%
  melt() %>%
  dplyr::select(value)
names(confint5.up) <- "upper"
confint5.down <- confint5[1,] %>%
  melt() %>%
  dplyr::select(value)
names(confint5.down) <- "lower"
odd5 <- exp(beta5)
names(odd5) <- "OR"
names(p5) <- "P"
r5 <- cbind(odd5, confint5.up, confint5.down, p5) %>%
  mutate(Phenotype = "Any internalizing and T2D")
r5$CNV <- rownames(r5)


cnv.reg6 <- lapply(preds, function(x) glm(formula(paste("IntKid ~", x, "+ Age + Sex + c01 + c02 + c03 + c04 + c05")), family="binomial", data=cnv.t))
names(cnv.reg6) <- preds
beta6 <- as.data.frame(sapply(cnv.reg6, function(x) summary(x)$coefficients[2,1]))
colnames(beta6) <- "IntKid"
p6 <- as.data.frame(sapply(cnv.reg6, function(x) summary(x)$coefficients[2,4]))
colnames(p6) <- "IntKid"
confint6 <- lapply(cnv.reg6, function(x){exp(confint.default(x)[2,])}) 
confint6 <- as.data.frame(confint6)
confint6$names <- rownames(confint6)
confint6.up <- confint6[2,] %>%
  melt() %>%
  dplyr::select(value)
names(confint6.up) <- "upper"
confint6.down <- confint6[1,] %>%
  melt() %>%
  dplyr::select(value)
names(confint6.down) <- "lower"
odd6 <- exp(beta6)
names(odd6) <- "OR"
names(p6) <- "P"
r6 <- cbind(odd6, confint6.up, confint6.down, p6) %>%
  mutate(Phenotype = "Any internalizing and CKD")
r6$CNV <- rownames(r6)


plot <- rbind(r1, r2, r3, r4, r5, r6) %>%
  filter(P < 0.0002525253)
plot$CNV <- gsub("_", " ", plot$CNV)
plot$CNV <- gsub("CNV", "", plot$CNV)


ggplot(plot, aes(x = CNV, y = OR, ymin = lower, ymax = upper,
                 col = Phenotype, fill = Phenotype)) + 
  geom_linerange(linewidth = 0.5, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_point(size = 3, shape = 21, colour = "white", stroke = 0.5,
             position = position_dodge(width = 0.5)) +
  scale_x_discrete(name = "CNV") +
  scale_y_continuous(name = "Odds ratio") + coord_flip(ylim = c(1, 30),
                                                       clip = 'off') +
  theme_minimal() +
  theme(legend.title=element_blank())



############################# CNV + PRS INTERACTION ###############################

prs.int <- read_csv("~/Documents/LINC/ukbb/edited/august.23.adjusted.prs.csv", col_names=T)  %>%
  inner_join(key, by = c("IID"="app15175")) %>%
  select(linc, anx.adjusted, mdd.adjusted)

prs.cm <- read_csv("~/Documents/LINC/ukbb/edited/oct.24.adjusted.cm.prs.csv", col_names=T) %>%
  inner_join(key, by = c("IID"="app15175")) %>%
  select(linc, ckd.adjusted, bmi.adjusted, t2d.adjusted, ldl.adjusted)

prs.sbp <- read_csv("~/Documents/LINC/ukbb/edited/apr.25.adjusted.sbp.prs.csv", col_names=T) %>%
  inner_join(key, by = c("IID"="app15175")) %>%
  select(linc, sbp.adjusted) 


multi.prs.cnv <- any.cnv %>%
  left_join(multi2) %>%
  left_join(prs.int, by=c("f.eid"="linc")) %>%
  left_join(prs.cm, by=c("f.eid"="linc")) %>%
  left_join(prs.sbp, by=c("f.eid"="linc")) %>%
  left_join(demo, by="f.eid")

test <- glm(anyCNV ~ anx.adjusted, data=multi.prs.cnv, family = "binomial")
summary(test)
confint <- exp(confint.default(test))

cnv.reg.mdd <- lapply(outs1, function(x) glm(formula(paste(x, "~ Age + anyCNV + mdd.adjusted + mdd.adjusted*anyCNV+ Townsend + Sex + c01 + c02 + c03 + c04 + c05")), family="binomial", data=multi.prs.cnv)) 
names(cnv.reg.mdd) <- outs1
lapply(cnv.reg.mdd, function(x){summary(x)})
pval <- sapply(cnv.reg.mdd, function(x) summary(x)$coefficients[12,4])
beta <- sapply(cnv.reg.mdd, function(x) summary(x)$coefficients[12,1])
or <- exp(beta)
confint <- lapply(cnv.reg.mdd, function(x){exp(confint.default(x)[12,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"
plot.mdd <- as.data.frame(cbind(cond, or, confint.down, confint.up, pval)) %>%
  mutate(PRS = "MDD")



cnv.reg.anx <- lapply(outs1, function(x) glm(formula(paste(x, "~ Age + anyCNV + anx.adjusted + anx.adjusted*anyCNV+ Townsend + Sex + c01 + c02 + c03 + c04 + c05")), family="binomial", data=multi.prs.cnv)) 
names(cnv.reg.anx) <- outs1
lapply(cnv.reg.anx, function(x){summary(x)})
pval <- sapply(cnv.reg.anx, function(x) summary(x)$coefficients[12,4])
beta <- sapply(cnv.reg.anx, function(x) summary(x)$coefficients[12,1])
or <- exp(beta)
confint <- lapply(cnv.reg.anx, function(x){exp(confint.default(x)[12,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"
plot.anx <- as.data.frame(cbind(cond, or, confint.down, confint.up, pval)) %>%
  mutate(PRS = "Anxiety")

cnv.reg.ckd <- lapply(outs1, function(x) glm(formula(paste(x, "~ Age + anyCNV + ckd.adjusted + ckd.adjusted*anyCNV+ Townsend + Sex + c01 + c02 + c03 + c04 + c05")), family="binomial", data=multi.prs.cnv)) 
names(cnv.reg.ckd) <- outs1
lapply(cnv.reg.ckd, function(x){summary(x)})
pval <- sapply(cnv.reg.ckd, function(x) summary(x)$coefficients[12,4])
beta <- sapply(cnv.reg.ckd, function(x) summary(x)$coefficients[12,1])
or <- exp(beta)
confint <- lapply(cnv.reg.ckd, function(x){exp(confint.default(x)[12,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"
plot.ckd <- as.data.frame(cbind(cond, or, confint.down, confint.up, pval)) %>%
  mutate(PRS = "CKD")

cnv.reg.bmi <- lapply(outs1, function(x) glm(formula(paste(x, "~ Age + anyCNV + bmi.adjusted + bmi.adjusted*anyCNV+ Townsend + Sex + c01 + c02 + c03 + c04 + c05")), family="binomial", data=multi.prs.cnv)) 
names(cnv.reg.bmi) <- outs1
lapply(cnv.reg.bmi, function(x){summary(x)})
pval <- sapply(cnv.reg.bmi, function(x) summary(x)$coefficients[12,4])
beta <- sapply(cnv.reg.bmi, function(x) summary(x)$coefficients[12,1])
or <- exp(beta)
confint <- lapply(cnv.reg.bmi, function(x){exp(confint.default(x)[12,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"
plot.bmi <- as.data.frame(cbind(cond, or, confint.down, confint.up, pval)) %>%
  mutate(PRS = "BMI")

cnv.reg.t2d <- lapply(outs1, function(x) glm(formula(paste(x, "~ Age + anyCNV + t2d.adjusted + t2d.adjusted*anyCNV+ Townsend + Sex + c01 + c02 + c03 + c04 + c05")), family="binomial", data=multi.prs.cnv)) 
names(cnv.reg.t2d) <- outs1
lapply(cnv.reg.t2d, function(x){summary(x)})
pval <- sapply(cnv.reg.t2d, function(x) summary(x)$coefficients[12,4])
beta <- sapply(cnv.reg.t2d, function(x) summary(x)$coefficients[12,1])
or <- exp(beta)
confint <- lapply(cnv.reg.t2d, function(x){exp(confint.default(x)[12,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"
plot.t2d <- as.data.frame(cbind(cond, or, confint.down, confint.up, pval)) %>%
  mutate(PRS = "T2D")

cnv.reg.ldl <- lapply(outs1, function(x) glm(formula(paste(x, "~ Age + anyCNV + ldl.adjusted + ldl.adjusted*anyCNV+ Townsend + Sex + c01 + c02 + c03 + c04 + c05")), family="binomial", data=multi.prs.cnv)) 
names(cnv.reg.ldl) <- outs1
lapply(cnv.reg.ldl, function(x){summary(x)})
pval <- sapply(cnv.reg.ldl, function(x) summary(x)$coefficients[12,4])
beta <- sapply(cnv.reg.ldl, function(x) summary(x)$coefficients[12,1])
or <- exp(beta)
confint <- lapply(cnv.reg.ldl, function(x){exp(confint.default(x)[12,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"
plot.ldl <- as.data.frame(cbind(cond, or, confint.down, confint.up, pval)) %>%
  mutate(PRS = "LDL")

cnv.reg.hyp <- lapply(outs, function(x) glm(formula(paste(x, "~ Age + anyCNV + Hypertension + Hypertension*anyCNV+ Townsend + Sex + c01 + c02 + c03 + c04 + c05")), family="binomial", data=multi.prs.cnv)) 
names(cnv.reg.hyp) <- outs
lapply(cnv.reg.hyp, function(x){summary(x)})

cnv.reg.sbp <- lapply(outs1, function(x) glm(formula(paste(x, "~ Age + anyCNV + sbp.adjusted + sbp.adjusted*anyCNV+ Townsend + Sex + c01 + c02 + c03 + c04 + c05")), family="binomial", data=multi.prs.cnv)) 
names(cnv.reg.sbp) <- outs1
lapply(cnv.reg.sbp, function(x){summary(x)})
pval <- sapply(cnv.reg.sbp, function(x) summary(x)$coefficients[12,4])
beta <- sapply(cnv.reg.sbp, function(x) summary(x)$coefficients[12,1])
or <- exp(beta)
confint <- lapply(cnv.reg.sbp, function(x){exp(confint.default(x)[12,])}) 
confint <- as.data.frame(confint)
confint$names <- rownames(confint)
confint.up <- confint[2,] %>%
  melt() %>%
  select(value)
names(confint.up) <- "upper"
confint.down <- confint[1,] %>%
  melt() %>%
  select(value)
names(confint.down) <- "lower"
plot.sbp <- as.data.frame(cbind(cond, or, confint.down, confint.up, pval)) %>%
  mutate(PRS = "SBP")

plot <- rbind(plot.mdd, plot.anx, plot.ckd, plot.bmi, plot.t2d, plot.sbp, plot.ldl)

write.table(file = "~/Desktop/prs.results.csv", plot, col.names = T, row.names= T, quote = F, sep=",")


