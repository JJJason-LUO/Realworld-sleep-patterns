setwd("D:/PheWAS")
library(tidyr)
library(dplyr)
library(viridis)
library(dplyr)
library(tidyr)
library(readr)
library(PheWAS)
library(ggplot2)
library(stringr)
library(ggrepel)
library(qqman)
library(purrr)
library(lubridate)
library(survival)
library(doParallel)
library(rms)
library(boot)
library(readxl)
library(gridExtra)
library(gfoRmula)
library(CMAverse)
library(survminer)
library(readxl)
#-----------------PheWAS ----------------------
cox_time <- read.csv("95559_cox_time_new.csv")
cox_status <- read.csv("95559_cox_status.csv")
cox_time <-cox_time[,c(2,5:2054)]
cox_status <- cox_status[,c(2:2052)]
exp_cov <- read.csv("exposure_covari.cleaned.csv")
exp_cov <- exp_cov[,c(2:29)]
exp_cov <- merge(df2[, c("ID", "session")], exp_cov, by = "ID")
covariates <- exp_cov_light[,c(1:18)]

exp_cov_light <- exp_cov %>%
  mutate(light_mean = n1_mean + n2_mean)

exp_cov_tst3 <- exp_cov_light %>%
  mutate(tst_mean = case_when(
    tst_mean < 420 ~ "0-6",        
    tst_mean <= 540 ~ "7-9",      
    TRUE ~ "9+" ))           

sa <- cox_status[,c("ID", "NS_333.1")]
exp_cov_sa <- exp_cov_light %>%
  left_join(sa, by = "ID")
exp_cov_sa <- exp_cov_sa %>% rename("SA" = NS_333.1)
exp_cov_sa$SA <- ifelse(exp_cov_sa$SA == 0, "F", "T")


combined_data <- cox_status %>%
  inner_join(cox_time, by = "ID") %>%
  inner_join(exp_cov_tst3, by = "ID")


phenotype_cols <- names(cox_status)[2:2051]  # 第2到2051列为疾病状态列

event_counts <- sapply(phenotype_cols, function(phecode) {
  sum(combined_data[[phecode]] == 1, na.rm = TRUE)
})
valid_phenotypes <- phenotype_cols[event_counts >= 20]

covariates <- setdiff(names(covariates), "ID")
stopifnot("NA" = length(covariates) > 0)
missing_covars <- setdiff(covariates, names(combined_data))
if (length(missing_covars) > 0) {
  stop("NE：", paste(missing_covars, collapse = ", "))
}

#DISCOVERY-FOCUSED
registerDoParallel(cores = 16)

df5 <- foreach(
  phecode = valid_phenotypes,
  .combine = bind_rows,
  .packages = c("dplyr", "survival")
) %dopar% {
  time_col <- paste0("time_", phecode)
  
  if (!time_col %in% names(combined_data)) {
    message("S：", phecode, "（F）")
    return(NULL)
  }
  
  # 构建公式
  formula_str <- paste(
    "Surv(", time_col, ", ", phecode, ") ~ light_mean",
    if (length(covariates) > 0) paste("+", paste(covariates, collapse = " + "))
  )
  formula <- as.formula(formula_str)
  
  tryCatch({
    fit <- coxph(formula, data = combined_data)
    summary_fit <- summary(fit)
    coef_df <- as.data.frame(summary_fit$coefficients["light_mean", , drop = FALSE])
    confint_df <- as.data.frame(summary_fit$conf.int["light_mean", , drop = FALSE])
    
    data.frame(
      Phenotype = phecode,
      
      HR = confint_df$`exp(coef)`,
      CI_lower = confint_df$`lower .95`,
      CI_upper = confint_df$`upper .95`,
      P_value = coef_df$`Pr(>|z|)`,
      N = fit$n,
      Events = sum(combined_data[[phecode]] == 1, na.rm = TRUE)
    )
  }, error = function(e) {
    message("F（", phecode, "）：", e$message)
    return(NULL)
  })
}

stopImplicitCluster()

# COMPARISON-FOCUSED
combined_data_06 <- combined_data %>%
  filter(tst_mean == "0-6" | tst_mean=="7-9")
combined_data_9p <- combined_data %>%
  filter(tst_mean == "7-9" | tst_mean=="9+")

combined_data_06 <- combined_data_06 %>%
  mutate(tst_mean_3 =if_else(tst_mean == "0-6",1,0))
combined_data_9p <- combined_data_9p %>%
  mutate(tst_mean_3 =if_else(tst_mean == "9+",1,0))

registerDoParallel(cores = 16)
df5 <- foreach(
  phecode = valid_phenotypes,
  .combine = bind_rows,
  .packages = c("dplyr", "survival")
) %dopar% {
  time_col <- paste0("time_", phecode)
  
  if (!time_col %in% names(combined_data_06)) {
    message("SK：", phecode, "（NA）")
    return(NULL)
  }
  
  # 构建公式
  formula_str <- paste(
    "Surv(", time_col, ", ", phecode, ") ~ tst_mean_3",  
    if (length(covariates) > 0) paste("+", paste(covariates, collapse = " + "))
  )
  formula <- as.formula(formula_str)
  
  tryCatch({
    fit <- coxph(formula, data = combined_data_06)
    summary_fit <- summary(fit)
    

    coef_names <- c("tst_mean_3")  
    results <- lapply(coef_names, function(name) {
      if (!name %in% rownames(summary_fit$coefficients)) {
        message("V ", name, " M ", phecode, " SK")
        return(NULL)
      }
      
      data.frame(
        Phenotype = phecode,
        HR = summary_fit$conf.int[name, "exp(coef)"],
        CI_lower = summary_fit$conf.int[name, "lower .95"],
        CI_upper = summary_fit$conf.int[name, "upper .95"],
        P_value = summary_fit$coefficients[name, "Pr(>|z|)"],
        N = fit$n,
        Events = sum(combined_data_06[[phecode]] == 1, na.rm = TRUE)
      )
    })
    
    do.call(rbind, Filter(Negate(is.null), results))
    
  }, error = function(e) {
    message("F（", phecode, "）：", e$message)
    return(NULL)
  })
}

stopImplicitCluster()

# PLOT
bonferroni_threshold <- 0.05/nrow(df5)
-log10(bonferroni_threshold)

df5 <- df5 %>%
  mutate(
    bonferroni = ifelse(P_value < bonferroni_threshold, 1, 0),
    Association = case_when(
      P_value < 0.05 & HR > 1 ~ "positive",
      P_value < 0.05 & HR < 1 ~ "negative",
      TRUE ~ "nonsig"
    )
  )
df5$log10_p <- -log10(df5$P_value)
df5$P_FDR <- p.adjust(df5$P_value, method = "fdr")

phemap <- read.csv("phenomap.csv")

phecode_info <- phemap %>% distinct(Phecode, Category, phecode_string)

plot_data <- df5 %>%
  left_join(phecode_info, by = c("Phenotype" = "Phecode")) %>%
  arrange(Category,phecode_string, Phenotype)  

summary(plot_data)
table(factor(plot_data$Association))

bonferroni_check <- subset(plot_data, log10_p > -log10(bonferroni_threshold))
p_check <- subset(plot_data, log10_p > -log10(0.05))
fdr_check <- subset(plot_data, P_FDR < 0.05) 	
table(factor(plot_data$Association))
table(factor(fdr_check$Association))
table(factor(bonferroni_check$Association))

plot_data$Association <- factor(plot_data$Association, levels = c("negative","positive", "nonsig"))

p <- ggplot(plot_data, aes(Category, log10_p))+#添加数据
  
  geom_jitter(aes(fill = Category, shape = Association, size = HR), color = "black", stroke = 1.5, width = 0.3)+#添加抖动的点，颜色根据门水平，形状根据change 大小根据LogFC
  
  scale_y_continuous(expand = c(0,0),limits = c(0,10))+#更改y轴范围
  
  scale_size_continuous(range = c(4,8))+#更改点的大小极值
  scale_shape_manual(values = c(21,1,46))+
  #scale_shape_manual(values = c(19,1,13))+#更改形状  
  
  scale_fill_manual(values = c("#B2182B", "#D6604D", "#F4A582", "#92C5DE", "#4393C3", "#2166AC",
                               "#FFDD55", "#E69F00", "#56B4E9", "#009E73","#B2182B", "#D6604D",
                               "#F4A582", "#92C5DE", "#4393C3", "#2166AC","#FFDD55"))+#更改颜色
  
  geom_hline(yintercept = -log10(0.05), size = 1, linetype = "dashed", color = "#999999")+
  
  geom_hline(yintercept = -log10(1.480613e-03), size = 1, linetype = "dashed", color = "#226FB8")+
  
  #geom_hline(yintercept = -log10(bonferroni_threshold), size = 1, linetype = "dashed", color = "#B2182B")+
  
  labs(x = "", y = "-Log10(P_value)")+
  
  theme_classic()+ #更改主题
  
  theme(legend.position = "top")+#更改图例位置  也可以是 left bottom right top  none  c(0.8,0.8)
  
  geom_text(
    data = subset(plot_data, log10_p >= 2.829559),
    aes(x = Category, y = log10_p, label = phecode_string),
    position = position_jitter(width = 0.5),  # 与geom_jitter同步抖动
    hjust = -0.5,
    vjust = -0.5,
    size = 3,
    color = "#999999"
  )  + 
  guides(
    fill = "none",          # 隐藏填充颜色的图例
    color = "none",         # 隐藏边框颜色的图例
    shape = "none",
    size = guide_legend(    # 显示 HR 的图例
      title = "HRs", 
      override.aes = list(shape = 16)  # 确保图例使用可填充形状
    )
  ) 
p

#----------------------RCS----------------------
rcs_data <- merge(cox_status[,c("ID", "CV_424")], cox_time[,c("ID","time_CV_424")])
rcs_data <- rcs_data %>%
  left_join(exp_cov, by = "ID")

dd <- datadist(rcs_data)
dd$limits["Adjust to", "tst_mean"] <- 480
options(datadist = "dd")
model <- cph(Surv(time_CV_424, CV_424) ~ rcs(tst_mean, 5) + Age + Sex + TDI + Ethnicity + Education + Employment.status + Income + ubrancity
             + Smoking + Alcohol + WatchingTV + Computerusing + Phoneusing + BMI +physicalactivity + PM2.5 + nighttime_noise, data = rcs_data)

model$stats["P"]

linear_model <- cph(Surv(time_CV_424, CV_424) ~ tst_mean + Age + Sex + TDI + Ethnicity + Education + Employment.status + Income + ubrancity
                    + Smoking + Alcohol + WatchingTV + Computerusing + Phoneusing + BMI +physicalactivity + PM2.5 + nighttime_noise, data = rcs_data)
anova_result <- anova(model)

print(anova_result)

lr_test <- lrtest(model, linear_model)
print(lr_test)

sleep_grid <- seq(min(rcs_data$tst_mean), max(rcs_data$tst_mean), length=1000)
pred <- Predict(model, tst_mean=sleep_grid, fun=exp, ref.zero = T)  

ref_index <- which.min(abs(sleep_grid - 480))
cat("refHR:", pred$yhat[ref_index], "\n")

min_hr_point <- sleep_grid[which.min(pred$yhat)]
cat("Low risk:", min_hr_point, "M\n")

mean(rcs_data$tst_mean, na.rm = T)

boot_func <- function(data, indices) {
  sample_data <- data[indices, ]
  boot_model <- update(model, data=sample_data)  
  pred <- Predict(boot_model, tst_mean=sleep_grid, fun=exp)
  sleep_grid[which.min(pred$yhat)]
}
set.seed(123)
boot_results <- boot(rcs_data, boot_func, R=500) 
boot_ci <- boot.ci(boot_results, type="perc")
cat("95% CI:", boot_ci$percent[4:5], "分钟\n")

min_risk <- min_hr_point
ci_lower <- boot_ci$percent[4]
ci_upper <- boot_ci$percent[5]
nonlinear_p <- anova_result[2,3]
global_p <- anova_result[20,3]
lr_p <- lr_test$stats[3]

ggplot(Predict(model, tst_mean, fun=exp, ref.zero = TRUE), aes(x = tst_mean, y = yhat)) +
  geom_line(color = "#92C5DE", linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#92C5DE") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1.3) +
  coord_cartesian(ylim = c(0, 6), clip = "off") +
  labs(x = "Sleep duration (min)", y = "HR(95%CI) for incidence", size = 8, caption = NULL) +
  annotate("text",
           x = min(sleep_grid), y = max(pred$upper)*0.9, 
           label = paste("Nonlinear P_value: ", format.pval(nonlinear_p, digits = 5), 
                         "\nLikelihood Ratio test P_value: ", format.pval(lr_p, digits = 5)),
           hjust = 0, vjust = 1,  
           size = 8, color = "black") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.position = "none",  
    plot.caption = element_blank()  
  )

results_df <- data.frame(
  matrix = c("min_point", 
             "95%_lower", 
             "95%_upper", 
             "nonlinear_p",
             "global_p",
             "lr_p"),
  value = c(
    round(min_risk, 1),
    round(ci_lower, 1),
    round(ci_upper, 1),
    format.pval(nonlinear_p, digits = 4, eps = 0.0001),
    format.pval(global_p, digits = 4, eps = 0.0001),
    format.pval(lr_p, digits = 4, eps = 0.0001)
  )
)

print(results_df)

#--------------------- sleep architecture------------------
exposure <- exp_cov_light[,c(1,20:21,25,27:29)]
exposure <- exposure %>%
  mutate(WHOrange = as.integer(tst_mean >= 420 & tst_mean <= 540))
exposure <- exposure[,c(1,8,2:7)]
exposure <- exposure %>%
  mutate(rem_tst = rem_mean/tst_mean, 
         n3_tst = n3_mean/tst_mean, 
         light_tst = light_mean/tst_mean,
         waso_tst = waso_mean/tst_mean, 
         irregularity_tst = sleep_irregular/tst_mean)

numeric_vars <- names(exposure)[sapply(exposure, is.numeric)]
numeric_vars <- numeric_vars[!numeric_vars %in% c("ID","WHOrange")]

summary_stats <- exposure %>%
  group_by(WHOrange) %>%
  summarise(across(
    all_of(numeric_vars),
    list(mean = ~mean(., na.rm = TRUE), 
         sd = ~sd(., na.rm = TRUE)
    ))) %>%
  pivot_wider(
    names_from = WHOrange,
    values_from = ends_with("_mean") | ends_with("_sd"),
    names_glue = "WHOrange_{.value}_{WHOrange}"
  ) %>%
  rename_with(~ gsub("_mean_", "_mean_", .x)) %>%
  rename_with(~ gsub("_sd_", "_sd_", .x))

p_values <- sapply(numeric_vars, function(var) {
  t.test(as.formula(paste(var, "~ WHOrange")), data = exposure)$p.value
}) 

final_result <- data.frame(
  variable = numeric_vars,
  who7_9_1_mean = as.numeric(summary_stats[1, grep("_mean_1", names(summary_stats))]),
  who7_9_1_sd = as.numeric(summary_stats[1, grep("_sd_1", names(summary_stats))]),
  who7_9_0_mean = as.numeric(summary_stats[1, grep("_mean_0", names(summary_stats))]),
  who7_9_0_sd = as.numeric(summary_stats[1, grep("_sd_0", names(summary_stats))]),
  p = p_values
)

print(final_result)

#-------------------CMA------------------------------
data06 <- combined_data_06 %>%
  select(ID,MB_280.13,time_MB_280.13,tst_mean_3)
data06 <- data06 %>%
  left_join(exp_cov_light, by = "ID")
data06_clean <- data06[!is.na(data06$time_MB_280.13), ]
data06_clean <- data06_clean %>%
  mutate(TREM = rem_mean/tst_mean,
         TDEEP = n3_mean/tst_mean,
         TLIGHT = light_mean/tst_mean,
         TWASO = waso_mean/tst_mean,
         TSI = sleep_irregular/tst_mean)

data9p <- combined_data_9p %>%
  select(ID,MB_280.13,time_MB_280.13,tst_mean_3)
data9p <- data9p %>%
  left_join(exp_cov_light, by = "ID")
data9p_clean <- data9p[!is.na(data9p$time_MB_280.13), ]
data9p_clean <- data9p_clean %>%
  mutate(TREM = rem_mean/tst_mean,
         TDEEP = n3_mean/tst_mean,
         TLIGHT = light_mean/tst_mean,
         TWASO = waso_mean/tst_mean,
         TSI = sleep_irregular/tst_mean)

mediator_list <- c('rem_mean','n3_mean','light_mean','waso_mean','sleep_irregular')
mediator_list <- c('TREM','TDEEP','TLIGHT','TWASO','TSI')

res_rb <- cmest(data = data9p_clean, 
                model='rb',
                exposure=c('tst_mean_3'),
                mediator = mediator_list,
                outcome = 'time_MB_280.13',
                event='MB_280.13',
                yreg = 'coxph',
                basec=c('Age', 'Sex' ),
                mreg=list('linear','linear','linear','linear','linear'),
                mval=list(0,0,0,0,0),
                EMint = T,
                astar=0,a=1,
                estimation='imputation',
                inference = 'bootstrap',
                nboot = 10,
                multimp = F)
summary(res_rb)

mediator_list <- c('TREM')
mediator_list <- c('TDEEP')
mediator_list <- c('TLIGHT')
mediator_list <- c('TSI')
mediator_list <- c('TWASO')

res_rb <- cmest(data = data9p_clean, 
                model='rb',
                exposure=c('tst_mean_3'),
                mediator = mediator_list,
                outcome = 'time_MB_280.13',
                event='MB_280.13',
                yreg = 'coxph',
                basec=c('Age', 'Sex' ),
                mreg=list('linear'),
                mval=list(0),
                EMint = T,
                astar=0,a=1,
                estimation='imputation',
                inference = 'bootstrap',
                nboot = 10,
                multimp = F)
summary(res_rb)

p <- cmdag(outcome = "time_CV_424", exposure = "tst_mean_3", mediator = c('rem_mean'), 
           basec=c('Age', 'Sex' , 'TDI' , 'Ethnicity' , 'Education' , 'Employment.status' , 'Income',
                   'ubrancity', 'Smoking' , 'Alcohol' , 'WatchingTV' , 'Computerusing' , 'Phoneusing' , 'BMI' , 'physicalactivity' ,
                   'PM2.5' , 'nighttime_noise'), postc = NULL, node = T, text_col = "white")
#-----------------------PAR------------------------------
cox_disease <- read.csv("95559_cox_time_new.csv")
cox_disease <- cox_disease[,c(2,4:2054)]
colnames(cox_disease) <- sub("^time_", "", colnames(cox_disease))
target_disease <- cox_disease %>% select(ID,refer_time, 
                                         CV_424,
                                         EM_256.2,
                                         GU_582.1,
                                         GU_594.3,
                                         GU_591,
                                         ID_003,
                                         MB_308.1,
                                         MB_280.13,
                                         NS_325,
                                         NS_350.5,
                                         SS_819)

reference <- target_disease[[2]]
cols_to_convert <- 3:13
target_disease[cols_to_convert] <- lapply(target_disease[cols_to_convert], function(x) {
  ifelse(
    !is.na(x),        
    ifelse(x == reference, 0, 1), 
    NA                  
  )
})
target_disease <- target_disease[,c(1,3:13)]

exposure <- read.csv("exposure_covari.cleaned.csv")
exposure <- exposure %>%
  mutate(
    range = case_when(
      tst_mean >= 420 & tst_mean <= 540 ~ 1,
      tst_mean < 420 ~ -1,
      TRUE ~ 0
    ))
exposure <- exposure[,c(2,30)]

par_data <- exposure %>%
  left_join(target_disease, "ID")


# 0-6（-1）vs 7-9（1)(NE)
par_results_06_vs_79 <- lapply(names(par_data)[3:ncol(par_data)], function(col_name) {
  
  temp <- par_data[par_data$range %in% c(-1, 1), c("range", col_name)]
  temp <- temp[!is.na(temp[[col_name]]), ]
  
  total_cases <- sum(temp[[col_name]] == 1)
  total_n <- nrow(temp)
  p_total <- total_cases / total_n
  
  unexposed <- temp[temp$range == 1, ]
  unexposed_cases <- sum(unexposed[[col_name]] == 1)
  unexposed_n <- nrow(unexposed)
  p_unexposed <- unexposed_cases / unexposed_n
  
  par <- (p_total - p_unexposed) / p_total * 100
  
  var_p_unexposed <- (p_unexposed * (1 - p_unexposed)) / unexposed_n
  var_p_total <- (p_total * (1 - p_total)) / total_n
  se <- sqrt(var_p_unexposed / (p_total^2) + var_p_total / (p_total^2)) * 100
  
  lower <- par - 1.96 * se
  upper <- par + 1.96 * se
  
  return(list(par = par, lower = lower, upper = upper))
  
})

par_df_06_vs_79 <- data.frame(
  variable = names(par_data)[3:ncol(par_data)],
  par = sapply(par_results_06_vs_79, function(x) x$par),
  lower_ci = sapply(par_results_06_vs_79, function(x) x$lower),
  upper_ci = sapply(par_results_06_vs_79, function(x) x$upper)
)

# 9+（0） vs 7-9（1）(refer)
par_results_79_vs_9plus <- lapply(names(par_data)[3:ncol(par_data)], function(col_name) {

  temp <- par_data[par_data$range %in% c(1, 0), c("range", col_name)]
  temp <- temp[!is.na(temp[[col_name]]), ]
  
  total_cases <- sum(temp[[col_name]] == 1)
  total_n <- nrow(temp)
  p_total <- total_cases / total_n
  
  unexposed <- temp[temp$range == 1, ]
  unexposed_cases <- sum(unexposed[[col_name]] == 1)
  unexposed_n <- nrow(unexposed)
  
  p_unexposed <- unexposed_cases / unexposed_n
  
  (p_total - p_unexposed)/p_total * 100
  
  par <- (p_total - p_unexposed) / p_total * 100
  
  var_p_unexposed <- (p_unexposed * (1 - p_unexposed)) / unexposed_n
  var_p_total <- (p_total * (1 - p_total)) / total_n
  se <- sqrt(var_p_unexposed / (p_total^2) + var_p_total / (p_total^2)) * 100
  
  lower <- par - 1.96 * se
  upper <- par + 1.96 * se
  
  return(list(par = par, lower = lower, upper = upper))
})

par_df_79_vs_9plus <- data.frame(
  variable = names(par_data)[3:ncol(par_data)],
  par = sapply(par_results_79_vs_9plus, function(x) x$par),
  lower_ci = sapply(par_results_79_vs_9plus, function(x) x$lower),
  upper_ci = sapply(par_results_79_vs_9plus, function(x) x$upper)
)

par_06 <- data.frame(
  Disease = names(par_results_06_vs_79),
  PAR_percent = unname(par_results_06_vs_79)
)

par_9plus <- data.frame(
  Disease = names(par_results_79_vs_9plus),
  PAR_percent = unname(par_results_79_vs_9plus)
)

par_df <- na.omit(par_df)

print(par_df_06_vs_79)
print(par_df_79_vs_9plus)