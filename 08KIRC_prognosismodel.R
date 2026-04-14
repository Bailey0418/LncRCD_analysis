setwd("D:/Project/03SCI/250710_LncRCD/Computational and Structural Biotechnology Journal/Supplementary/KIRC/Prognosis")
a <- read.csv("E:/RCD/data/cancer/KIRC/Apoptosis.csv")
b <- read.csv("E:/RCD/data/cancer/KIRC/Cuproptosis.csv")
c <- read.csv("E:/RCD/data/cancer/KIRC/Necroptosis.csv")
d <- read.csv("E:/RCD/data/cancer/KIRC/Pyroptosis.csv")
e <- rbind(a,b,c,d)
lnc <- unique(e$lncRNA)
rm(a,b,c,d,e)

library(data.table)
library(dplyr)
library(stringr)
library(tibble)
clinical <- fread("D:/Project/03SCI/250710_LncRCD/Computational and Structural Biotechnology Journal/Supplementary/KIRC/Prognosis/TCGA-KIRC.clinical.tsv.gz", data.table = FALSE)
survival <- fread("D:/Project/03SCI/250710_LncRCD/Computational and Structural Biotechnology Journal/Supplementary/KIRC/Prognosis/TCGA-KIRC.survival.tsv.gz", data.table = FALSE)
fpkm <- read.csv("E:/RCD/data/cancer/KIRC/fpkm/TCGA_KIRC_lncRNA_tum_fpkm.csv", row.names = 1)
colnames(fpkm) <- gsub("\\.", "-", colnames(fpkm))

surv2 <- survival %>%
  transmute(
    sample = as.character(sample),
    patient_id = as.character(`_PATIENT`),
    time = as.numeric(`OS.time`),
    status = as.numeric(OS)
  ) %>%
  filter(!is.na(sample), !is.na(patient_id), !is.na(time), !is.na(status)) %>%
  filter(time > 0) %>%
  distinct(patient_id, .keep_all = TRUE)
write.csv(surv2,"surv.csv")

clin2 <- clinical %>%
  transmute(
    sample = as.character(sample),
    patient_id = as.character(submitter_id),
    age = as.numeric(age_at_index.demographic),
    gender = as.character(gender.demographic),
    stage = as.character(ajcc_pathologic_stage.diagnoses),
    T = as.character(ajcc_pathologic_t.diagnoses),
    M = as.character(ajcc_pathologic_m.diagnoses),
    N = as.character(ajcc_pathologic_n.diagnoses),
    sample_type = as.character(sample_type.samples)
  ) %>%
  filter(!is.na(patient_id)) %>%
  distinct(patient_id, .keep_all = TRUE)
write.csv(clin2,"clin.csv")

expr_mat <- as.matrix(fpkm)
expr_mat <- t(expr_mat)
expr_mat <- as.data.frame(expr_mat)
expr_mat$sample <- rownames(expr_mat)
rownames(expr_mat) <- NULL
expr_mat$sample_type_code <- substr(expr_mat$sample, 14, 15)
table(expr_mat$sample_type_code)
expr_tumor <- expr_mat %>%
  filter(sample_type_code == "01")
expr_tumor$patient_id <- substr(expr_tumor$sample, 1, 12)
expr_tumor <- expr_tumor %>%
  distinct(patient_id, .keep_all = TRUE)
meta_cols <- c("sample", "sample_type_code", "patient_id")
expr_cols <- setdiff(colnames(expr_tumor), meta_cols)

dat_model <- expr_tumor %>%
  inner_join(surv2, by = "patient_id") %>%
  left_join(clin2 %>% select(-sample), by = "patient_id")

if("sample.x" %in% colnames(dat_model)) {
  dat_model$sample <- dat_model$sample.x}
dat_model <- dat_model %>%
  select(-any_of(c("sample.x", "sample.y")))

dat_model <- dat_model %>%
  filter(!is.na(time), !is.na(status), time > 0)
write.table(dat_model,file = "D:/Project/03SCI/250710_LncRCD/Computational and Structural Biotechnology Journal/Supplementary/KIRC/Prognosis/KIRC_model_data_all.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

library(survival)
library(survminer)
library(glmnet)
library(timeROC)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(rms)
dat <- dat_model %>%
  select(patient_id, sample, time, status, age, gender, stage, T, M, N, all_of(lnc))
dat <- dat %>%
  filter(!is.na(time), !is.na(status), time > 0)
dat_clin <- dat %>%
  filter(!is.na(age), !is.na(gender), !is.na(stage))
#单因素cox
set.seed(1234)
train_id <- sample(1:nrow(dat), size = floor(nrow(dat) * 0.5))
train <- dat[train_id, ]
test  <- dat[-train_id, ]
candidate_lnc <- lnc
uni_cox <- lapply(candidate_lnc, function(gene) {
  formula <- as.formula(paste0("Surv(time, status) ~ `", gene, "`"))
  fit <- coxph(formula, data = train)
  s <- summary(fit)
  data.frame(
    gene = gene,
    HR = s$coefficients[,"exp(coef)"],
    lower95 = s$conf.int[,"lower .95"],
    upper95 = s$conf.int[,"upper .95"],
    pvalue = s$coefficients[,"Pr(>|z|)"]
  )
})
uni_cox_res <- do.call(rbind, uni_cox)
uni_cox_res <- uni_cox_res[order(uni_cox_res$pvalue), ]
uni_cox_res
sig_lnc <- uni_cox_res$gene[uni_cox_res$pvalue < 0.05]
sig_lnc
write.csv(uni_cox_res,"uni_cox_res.csv")
#LASSO
x_train <- as.matrix(train[, sig_lnc, drop = FALSE])
y_train <- Surv(train$time, train$status)
set.seed(1234)
cvfit <- cv.glmnet(
  x_train,
  y_train,
  family = "cox",
  alpha = 1,
  nfolds = 10)
plot(cvfit)
fit <- glmnet(x_train, y_train, family = "cox", alpha = 1)
plot(fit, xvar = "lambda")
coef_min <- coef(cvfit, s = "lambda.min")
active_index <- which(coef_min != 0)
lasso_genes <- rownames(coef_min)[active_index]
lasso_genes
lasso_coef <- as.matrix(coef_min[active_index])
lasso_coef
#多因素cox
final_formula <- as.formula(
  paste("Surv(time, status) ~", paste(paste0("`", lasso_genes, "`"), collapse = " + ")))
final_cox <- coxph(final_formula, data = train)
summary(final_cox)
final_coef <- summary(final_cox)$coefficients
final_coef
write.csv(final_coef,"final_coef.csv")
#risk score
train$riskScore <- predict(final_cox, type = "lp", newdata = train)
test$riskScore <- predict(final_cox, type = "lp", newdata = test)
dat$riskScore <- predict(final_cox, type = "lp", newdata = dat)
cut_res <- surv_cutpoint(
  data = train[, c("time", "status", "riskScore")],
  time = "time",
  event = "status",
  variables = "riskScore"
)
cutoff <- cut_res$cutpoint[1, "cutpoint"]
cutoff
plot(cut_res, "riskScore", palette = "npg")

train$riskGroup <- ifelse(train$riskScore > cutoff, "High", "Low")
test$riskGroup  <- ifelse(test$riskScore > cutoff, "High", "Low")
dat$riskGroup   <- ifelse(dat$riskScore > cutoff, "High", "Low")
train$riskGroup <- factor(train$riskGroup, levels = c("Low", "High"))
test$riskGroup  <- factor(test$riskGroup, levels = c("Low", "High"))
dat$riskGroup   <- factor(dat$riskGroup, levels = c("Low", "High"))
#KM
fit_km_train <- survfit(Surv(time, status) ~ riskGroup, data = train)
ggsurvplot(
  fit_km_train,
  data = train,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata", # 根据分层更改风险表颜色
  linetype = "strata", # 根据分层更改线型
  surv.median.line = "hv", # 同时显示垂直和水平参考线
  ggtheme = theme_bw(), # 更改ggplot2的主题
  palette = c("#E7B800", "#2E9FDF"))#定义颜色

fit_km_test <- survfit(Surv(time, status) ~ riskGroup, data = test)
ggsurvplot(
  fit_km_test,
  data = test,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata", # 根据分层更改风险表颜色
  linetype = "strata", # 根据分层更改线型
  surv.median.line = "hv", # 同时显示垂直和水平参考线
  ggtheme = theme_bw(), # 更改ggplot2的主题
  palette = c("#E7B800", "#2E9FDF"))#定义颜色
#ROC
roc_train <- timeROC(
  T = train$time,
  delta = train$status,
  marker = train$riskScore,
  cause = 1,
  times = c(365, 1095, 1825),
  iid = TRUE)
plot(roc_train, time = 365, col = "red", title = FALSE)
plot(roc_train, time = 1095, add = TRUE, col = "blue")
plot(roc_train, time = 1825, add = TRUE, col = "green")
legend("bottomright",
       legend = c(
         paste0("1-year AUC=", round(roc_train$AUC[1], 3)),
         paste0("3-year AUC=", round(roc_train$AUC[2], 3)),
         paste0("5-year AUC=", round(roc_train$AUC[3], 3))
       ),
       col = c("red", "blue", "green"),
       lwd = 2)

roc_test <- timeROC(
  T = test$time,
  delta = test$status,
  marker = test$riskScore,
  cause = 1,
  times = c(365, 1095, 1825),
  iid = TRUE)
plot(roc_test, time = 365, col = "red", title = FALSE)
plot(roc_test, time = 1095, add = TRUE, col = "blue")
plot(roc_test, time = 1825, add = TRUE, col = "green")
legend("bottomright",
       legend = c(
         paste0("1-year AUC=", round(roc_test$AUC[1], 3)),
         paste0("3-year AUC=", round(roc_test$AUC[2], 3)),
         paste0("5-year AUC=", round(roc_test$AUC[3], 3))
       ),
       col = c("red", "blue", "green"),
       lwd = 2)
#独立预后分析
dat2 <- dat[,c("sample", "time", "status", "riskScore", "age", "gender", "stage")]
dat2 <- as.data.frame(dat2)
dat2$time <- as.numeric(dat2$time)
dat2$status <- as.numeric(dat2$status)
dat2$gender <- as.factor(dat2$gender)
dat2$stage  <- as.factor(dat2$stage)
dat$age <- as.numeric(dat$age)
vars <- c("riskScore", "age", "gender", "stage")
uni_res <- lapply(vars, function(v){
  fml <- as.formula(paste0("Surv(time, status) ~ ", v))
  fit <- coxph(fml, data = dat2)
  s <- summary(fit)
  
  data.frame(
    Variable = rownames(s$coefficients),
    HR = s$coefficients[, "exp(coef)"],
    lower95 = s$conf.int[, "lower .95"],
    upper95 = s$conf.int[, "upper .95"],
    pvalue = s$coefficients[, "Pr(>|z|)"],
    row.names = NULL
  )
})
uni_res <- do.call(rbind, uni_res)
uni_res

dat2$gender <- factor(dat2$gender)
dat2$gender <- relevel(dat2$gender, ref = "female")
dat2$stage <- factor(dat2$stage)
dat2$stage <- relevel(dat2$stage, ref = "Stage I")
multi_fit <- coxph(Surv(time, status) ~ riskScore + age + gender + stage, data = dat2)
multi_sum <- summary(multi_fit)
multi_res <- data.frame(
  Variable = rownames(multi_sum$coefficients),
  HR = multi_sum$coefficients[, "exp(coef)"],
  lower95 = multi_sum$conf.int[, "lower .95"],
  upper95 = multi_sum$conf.int[, "upper .95"],
  pvalue = multi_sum$coefficients[, "Pr(>|z|)"],
  row.names = NULL)
write.csv(multi_res, "multivariate_cox_result_yuhou.csv", row.names = FALSE)

ggforest(
  multi_fit,
  data = dat2,
  main = "Multivariate Cox regression",
  fontsize = 1.0,
  refLabel = "Reference",
  noDigits = 3)
