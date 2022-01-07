# "LOOCV"

library(data.table)
library(tidyverse)
library(survminer)
library(survival)

steps <- seq(0.02, 0.98, by = 0.02)

loo_cv <- bind_rows(lapply(steps, function(step) {
  temp <- data_in %>% filter(aneuploidy_score >= !!step)
  
  bind_rows(lapply(1:nrow(temp), function(i) {
    temp_loo <- temp[-i,]
    if(length(unique(temp_loo$timing)) < 2) return()
    s <- summary(survfit(Surv(os_time, os_event) ~ timing, data = temp_loo), times = 6, extend = TRUE)
    s2 <- survdiff(Surv(os_time, os_event) ~ timing, data = temp_loo)
    list(aneuploidy_thresh = step,
         iteration = i, 
         p_value = 1 - pchisq(s2$chisq, length(s2$n) - 1),
         one_year_surv1 = s$surv[1],
         one_year_surv2 = s$surv[2],
         surv_diff =  abs(s$surv[1] -  s$surv[2]))
  }))
}))

loocv_results_p <- bind_rows(lapply(unique(loo_cv$aneuploidy_thresh), function(x) {
  temp <- loo_cv %>% filter(aneuploidy_thresh == !!x)
  m <- mean(temp$p_value)
  t <- qt(p = 0.025, df = nrow(temp) - 1, lower.tail = F)
  se <- sd(temp$p_value) / sqrt(nrow(temp))
  marg <- t * se
  list(aneuploidy_thresh = x,
       mean_p = m,
       lower_ci = m - marg,
       upper_ci = m + marg)
}))

loocv_results_diff <- bind_rows(lapply(unique(loo_cv$aneuploidy_thresh), function(x) {
  temp <- loo_cv %>% filter(aneuploidy_thresh == !!x)
  m <- mean(temp$surv_diff)
  t <- qt(p = 0.025, df = nrow(temp) - 1, lower.tail = F)
  se <- sd(temp$surv_diff) / sqrt(nrow(temp))
  marg <- t * se
  list(aneuploidy_thresh = x,
       mean_surv_diff = m,
       lower_ci = m - marg,
       upper_ci = m + marg)
}))


loocv_results_comb <- bind_rows(loocv_results_p %>% mutate(type = "p") %>% dplyr::rename(mean = mean_p),
                                loocv_results_diff %>% mutate(type = "diff") %>% dplyr::rename(mean = mean_surv_diff))

best_thresh <- loocv_results_comb %>% dplyr::select(aneuploidy_thresh, mean, type) %>%
  spread(type, mean) %>%
  mutate(metric = diff - p)
