# "COSINR Enrichment and Survival Analyses"

library(data.table)
library(tidyverse)
library(qvalue)
library(survminer)
library(survival)
library(rms)

run_survival <- function(data, censor_var, time_var, cov = NA, start_col = 55, count_thresh = 3) {
  if(start_col >= ncol(data)) return()
  
  bind_rows(lapply(names(data)[start_col:ncol(data)], function(x) {
    if(length(levels(factor(data[,x]))) == 1) return()
    if(!is.na(count_thresh)) {
      if (sum(data[,x]) < count_thresh) return()
    } 
    
    f1 <- as.formula(paste0("Surv(", time_var, ", ", censor_var, ") ~ ", x)) 
    if(!is.na(cov)) {
      f2 <- as.formula(paste0("Surv(", time_var, ", ", censor_var, ") ~ ", x, " + ", paste(cov, collapse = " +"))) 
    } else {
      f2 <- as.formula(paste0("Surv(", time_var, ", ", censor_var, ") ~ ", x)) 
    }
    s <- survdiff(f1, data = data)
    c <- summary(coxph(f2, data = data))
    multicol_flag <- any(vif(cph(f2, data = data)) > 4)
    data.frame(event = x, 
               p_value_wald =c$coefficients[1,5],
               p_value_logrank = pchisq(s$chisq, length(s$n) - 1, lower.tail = FALSE), 
               cov = paste(cov, collapse = ", "), 
               groups = paste(rownames(s$n), collapse = ", "),
               ns = paste(s$n, collapse = ", "),
               hazard_ratio = c$coefficients[1,2],
               medians = paste(round(surv_median(survfit(f1, data = data))[,2], 2), collapse = " vs. "),
               multicoll_detected = multicol_flag)
  }))
}

get_q_vals <- function(df, p_value_col = "p_value") {
  if(!p_value_col %in% names(df)) {
    message("p_value column not found, q-values not computed")
    return(df)
  }
  if(nrow(df) < 2) {
    message("not enough observations to compute q-values")
    return(df)
  }
  df.mut <- drop_na(df, any_of(p_value_col))
  df.mut[,p_value_col][df.mut[,p_value_col] >= 1] <- 1
  df.mut[,p_value_col][df.mut[,p_value_col] <= 0] <- 0
  df.mut <- cbind(df.mut, q_value = qvalue(unlist(df.mut[,p_value_col]), pi0 = 1)$qvalues)
  names(df.mut)[names(df.mut) == "q_value"] <- gsub("^p", "q", p_value_col)
  df.mut
}

by_var <- function(df, var, f, ...) {
  df_a <- split(df, df[,var]) # split the data frame
  df_out <- lapply(df_a, f, ...)  # apply the function to the sections
  bind_rows(df_out, .id = "GROUP") # rejoin the sections
}

compare_all_events <- function(df, var, start_index) {
  if(start_index > ncol(df)) return(data.frame())
  df_out <- bind_rows(lapply(names(df)[start_index:length(names(df))], 
                             function(x) {
                               temp = df %>% dplyr::select(one_of(var, x)) %>% drop_na() %>% group_by_(var, x) %>% 
                                 dplyr::summarize(count = n(), .groups = "keep") %>% ungroup() %>% spread_(var, "count") 
                               temp = temp[order(c(unlist(abs(temp[,1]))), decreasing = T),]
                               r = unlist(temp[,1])
                               temp = as.matrix(temp[,-1])
                               rownames(temp) = r
                               temp[is.na(temp)] = 0
                               out = compare_events(temp)
                               if(class(out) == "NULL") return()
                               out = out %>% mutate(event = !!x) %>% dplyr::select(event, everything())
                             }))
  if(class(df_out) == "NULL") return()
  df_out
}