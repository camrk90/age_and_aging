library(tidyverse)
library(ggplot2)
library(ggeffects)
library(lme4)

# Simulation settings
n_individuals<- 100
beta_0<- 10
beta_within<- 1
beta_between<- -0.5
sigma<- 2

df_list <- vector("list", 1000)

for (i in seq_along(df_list)) {
  
  sim_data <- do.call(rbind, lapply(seq_len(n_individuals), function(id) {
    
    n_samples <- sample(2:5, size = 1)
    baseline_age <- runif(1, min = 5, max = 30)
    
    age_offsets <- sort(runif(n_samples, min = 0, max = 5))
    age <- baseline_age + age_offsets
    
    data.frame(
      id = id,
      sample = seq_len(n_samples),
      age = age
    )
  }))
  
  sim_data <- sim_data %>%
    group_by(id) %>%
    mutate(n = n(),
           mean_age = mean(age),
           within_age = age - mean_age) %>%
    ungroup() %>%
    mutate(y = beta_0 + beta_within * within_age + beta_between * mean_age +
             rnorm(nrow(.), mean = 0, sd = sigma))
  
  df_list[[i]] <- sim_data
}

mods<- lapply(df_list, function(x){
  
  #Generate models
  eq1<- lmer(y ~ age + (1|id), data = x)
  eq1_summary<- summary(eq1)
  eq2<- lmer(y ~ within_age + mean_age + (1|id), data = x)
  eq2_summary<- summary(eq2)
  eq3<- lmer(y ~ age + mean_age + (1|id), data = x)
  eq3_summary<- summary(eq3)
  
  mod_list<- list(eq1=eq1, eq1_summary=summary(eq1),
                  eq2=eq2, eq2_summary=summary(eq2),
                  eq3=eq3, eq3_summary=summary(eq3))
  
})
names(mods)<- 1:1000

coefs<- lapply(names(mods), function(x){
  
  mod<- mods[[x]]
  mod_df<- rbind(mod[["eq2_summary"]]$coefficients[, 1])
  mod_df
  
})

coefs<- as.data.frame(do.call(rbind, coefs))

coefs %>%
  ggplot(aes(within_age)) +
  geom_histogram(aes(fill = after_stat(x)), bins = 50) +
  geom_vline(xintercept = mean(coefs$within_age), linetype = "dashed", colour = "red") +
  scale_fill_gradient2(low = "white", mid = "purple3", high = "white", midpoint = -0.5) +
  theme_classic(base_size=6) +
  theme(legend.key.width = unit(2, 'mm'), 
        legend.key.height = unit(5, 'mm'),
        legend.position = "none") +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  xlab(expression(beta["Eq.3"])) +
  ylab("Count")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/sim_hist2.svg",
       plot = last_plot(),
       height = 25, width = 25, units = "mm")

coefs %>%
  ggplot(aes(mean_age)) +
  geom_histogram(aes(fill = after_stat(x)), bins = 50, colour = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "red") +
  scale_fill_gradient2(low = "white", mid = "grey30", high = "white", midpoint = 1) +
  theme_classic(base_size=12) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5))

eq1_coefs<- lapply(names(mods), function(x){
  
  mod<- mods[[x]]
  mod_df<- rbind(mod[["eq1_summary"]]$coefficients[, 1])
  mod_df
  
})

eq1_coefs<- as.data.frame(do.call(rbind, eq1_coefs))

eq1_coefs %>%
  ggplot(aes(age)) +
  geom_histogram(aes(fill = after_stat(x)), bins = 50, colour = "black") +
  geom_vline(xintercept = mean(eq1_coefs$age), linetype = "dashed", colour = "red") +
  scale_fill_gradient2(low = "white", mid = "steelblue3", high = "white", midpoint = mean(eq1_coefs$age)) +
  theme_classic(base_size=12) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5))

fe.chron<- predict_response(mods[["74"]][["eq1"]], "age")

fe.eq2<- predict_response(mods[["74"]][["eq2"]], "within_age")

fe.btwn<- predict_response(mods[["74"]][["eq2"]], "mean_age")

fe.eq3 <- predict_response(mods[["74"]][["eq3"]], "age")

ggplot() + 
  geom_point(data = sim_data, inherit.aes = FALSE, aes(x = age, y = y)) +
  geom_path(data = sim_data, inherit.aes = FALSE, aes(x = age, y = y, group = id)) +
  geom_line(data = fe.btwn, inherit.aes = FALSE, 
            aes(x = x, y = predicted),
            color = "grey30",
            linewidth = 1) +
  geom_line(data = fe.chron,
            inherit.aes = FALSE,
            aes(x = x,
                y = predicted),
            color = "steelblue2",
            linewidth = 1) +
  geom_line(data = fe.eq3,
            inherit.aes = FALSE,
            aes(x = x,
                y = predicted),
            color = "purple",
            linewidth = 1) +
  theme_classic(base_size=6) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  xlab("Age") +
  ylab("Response")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/sim2.svg",
       plot = last_plot(),
       height = 50, width = 50, units = "mm")

#BetaW and BetaB match

# Simulation settings
n_individuals<- 100
beta_0<- 10
beta_within<- 1
beta_between<- 1
sigma<- 2

df_list <- vector("list", 1000)

for (i in seq_along(df_list)) {
  
  sim_data <- do.call(rbind, lapply(seq_len(n_individuals), function(id) {
    
    n_samples <- sample(2:5, size = 1)
    baseline_age <- runif(1, min = 5, max = 30)
    
    age_offsets <- sort(runif(n_samples, min = 0, max = 5))
    age <- baseline_age + age_offsets
    
    data.frame(
      id = id,
      sample = seq_len(n_samples),
      age = age
    )
  }))
  
  sim_data <- sim_data %>%
    group_by(id) %>%
    mutate(n = n(),
           mean_age = mean(age),
           within_age = age - mean_age) %>%
    ungroup() %>%
    mutate(y = beta_0 + beta_within * within_age + beta_between * mean_age +
             rnorm(nrow(.), mean = 0, sd = sigma))
  
  df_list[[i]] <- sim_data
}

mods<- lapply(df_list, function(x){
  
  #Generate models
  eq1<- lmer(y ~ age + (1|id), data = x)
  eq1_summary<- summary(eq1)
  eq2<- lmer(y ~ within_age + mean_age + (1|id), data = x)
  eq2_summary<- summary(eq2)
  eq3<- lmer(y ~ age + mean_age + (1|id), data = x)
  eq3_summary<- summary(eq3)
  
  mod_list<- list(eq1=eq1, eq1_summary=summary(eq1),
                  eq2=eq2, eq2_summary=summary(eq2),
                  eq3=eq3, eq3_summary=summary(eq3))
  
})
names(mods)<- 1:1000

coefs<- lapply(names(mods), function(x){
  
  mod<- mods[[x]]
  mod_df<- rbind(mod[["eq2_summary"]]$coefficients[, 1])
  mod_df
  
})

coefs<- as.data.frame(do.call(rbind, coefs))

fe.chron<- predict_response(mods[["74"]][["eq1"]], "age")

fe.eq2<- predict_response(mods[["74"]][["eq2"]], "within_age")

fe.btwn<- predict_response(mods[["74"]][["eq2"]], "mean_age")

fe.eq3 <- predict_response(mods[["74"]][["eq3"]], "age")

ggplot() + 
  geom_point(data = sim_data, inherit.aes = FALSE, aes(x = age, y = y)) +
  geom_path(data = sim_data, inherit.aes = FALSE, aes(x = age, y = y, group = id)) +
  geom_line(data = fe.btwn, inherit.aes = FALSE, 
            aes(x = x, y = predicted),
            color = "grey30",
            linewidth = 1) +
  geom_line(data = fe.chron,
            inherit.aes = FALSE,
            aes(x = x,
                y = predicted),
            color = "steelblue2",
            linewidth = 1) +
  geom_line(data = fe.eq3,
            inherit.aes = FALSE,
            aes(x = x,
                y = predicted),
            color = "purple",
            linewidth = 1) +
  theme_classic(base_size=6) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  xlab("Age") +
  ylab("Response")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/sim1.svg",
       plot = last_plot(),
       height = 50, width = 50, units = "mm")
  

#BetaW and BetaB same direction

# Simulation settings
n_individuals<- 100
beta_0<- 10
beta_within<- 1
beta_between<- 0.5
sigma<- 2

df_list <- vector("list", 1000)

for (i in seq_along(df_list)) {
  
  sim_data <- do.call(rbind, lapply(seq_len(n_individuals), function(id) {
    
    n_samples <- sample(2:5, size = 1)
    baseline_age <- runif(1, min = 5, max = 30)
    
    age_offsets <- sort(runif(n_samples, min = 0, max = 5))
    age <- baseline_age + age_offsets
    
    data.frame(
      id = id,
      sample = seq_len(n_samples),
      age = age
    )
  }))
  
  sim_data <- sim_data %>%
    group_by(id) %>%
    mutate(n = n(),
           mean_age = mean(age),
           within_age = age - mean_age) %>%
    ungroup() %>%
    mutate(y = beta_0 + beta_within * within_age + beta_between * mean_age +
             rnorm(nrow(.), mean = 0, sd = sigma))
  
  df_list[[i]] <- sim_data
}

mods<- lapply(df_list, function(x){
  
  #Generate models
  eq1<- lmer(y ~ age + (1|id), data = x)
  eq1_summary<- summary(eq1)
  eq2<- lmer(y ~ within_age + mean_age + (1|id), data = x)
  eq2_summary<- summary(eq2)
  eq3<- lmer(y ~ age + mean_age + (1|id), data = x)
  eq3_summary<- summary(eq3)
  
  mod_list<- list(eq1=eq1, eq1_summary=summary(eq1),
                  eq2=eq2, eq2_summary=summary(eq2),
                  eq3=eq3, eq3_summary=summary(eq3))
  
})
names(mods)<- 1:1000

coefs<- lapply(names(mods), function(x){
  
  mod<- mods[[x]]
  mod_df<- rbind(mod[["eq2_summary"]]$coefficients[, 1])
  mod_df
  
})

coefs<- as.data.frame(do.call(rbind, coefs))

fe.chron<- predict_response(mods[["74"]][["eq1"]], "age")

fe.eq2<- predict_response(mods[["74"]][["eq2"]], "within_age")

fe.btwn<- predict_response(mods[["74"]][["eq2"]], "mean_age")

fe.eq3 <- predict_response(mods[["74"]][["eq3"]], "age")

ggplot() + 
  geom_point(data = sim_data, inherit.aes = FALSE, aes(x = age, y = y)) +
  geom_path(data = sim_data, inherit.aes = FALSE, aes(x = age, y = y, group = id)) +
  geom_line(data = fe.btwn, inherit.aes = FALSE, 
            aes(x = x, y = predicted),
            color = "grey30",
            linewidth = 1) +
  geom_line(data = fe.chron,
            inherit.aes = FALSE,
            aes(x = x,
                y = predicted),
            color = "steelblue2",
            linewidth = 1) +
  geom_line(data = fe.eq3,
            inherit.aes = FALSE,
            aes(x = x,
                y = predicted),
            color = "purple",
            linewidth = 1) +
  theme_classic(base_size=6) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  xlab("Age") +
  ylab("Response")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/sim3.svg",
       plot = last_plot(),
       height = 50, width = 50, units = "mm")



