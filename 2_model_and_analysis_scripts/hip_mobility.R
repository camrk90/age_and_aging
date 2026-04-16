library(tidyverse)
library(ggplot2)
library(ggeffects)
library(lme4)
library(lmerTest)
library(stargazer)

#Define functions---------------------------------------------------------------
generate_models<- function(dat, var){
  
  var <- deparse(substitute(var))
  
  f1 <- as.formula(paste(var, "~ age + individual_sex + (1|individual_code)"))
  f2 <- as.formula(paste(var, "~ within_age + between_age + individual_sex + (1|individual_code)"))
  f3 <- as.formula(paste(var, "~ age + between_age + individual_sex + (1|individual_code)"))
  
  chron <- lmer(f1, data = dat)
  eq2   <- lmer(f2, data = dat)
  eq3   <- lmer(f3, data = dat)

  coefs<- as.data.frame(rbind(summary(chron)[["coefficients"]], summary(eq2)[["coefficients"]],
                                    summary(eq3)[["coefficients"]]))
  colnames(hip_ext_table)<- c("β", "se", "df", "t", "p", "2.5%", "97.5%")
  
  coefs <- round(coefs, 3)
  
  rownames(coefs)<- c("eq1_intercept", "eq1_age", "eq1_sexM",
                      "eq2_intercept", "eq2_age.w", "eq2_age.btwn", "eq2_sexM",
                      "eq3_intercept", "eq3_age.w", "eq3_age.btwn", "eq3_sexM")
  
  conf_ints<- as.data.frame(rbind(confint(chron), 
                                  confint(eq2),
                                  confint(eq3)))
  
  conf_ints<- conf_ints[!grepl("sig", rownames(conf_ints)),]
  
  coefs<- cbind(coefs, conf_ints)
  
  fe.chron<- predict_response(chron, "age")
  #re.chron<- predict_response(chron_mod, terms = c("age", "individual_code"), type = "random")
  
  fe.eq2<- predict_response(eq2, "within_age")
  
  fe.btwn<- predict_response(eq2, "between_age")
  #re.eq2<- predict_response(eq2_mod, terms = c("within_age", "individual_code"), type = "random")
  
  fe.eq3 <- predict_response(eq3, "age")
  #re.eq3<- predict_response(eq3_mod, terms = c("age", "individual_code"), type = "random")
  
  return(list(chron_mod = chron, eq2_mod = eq2, eq3_mod = eq3, 
              coefs = coefs, fe.chron=fe.chron, fe.eq2=fe.eq2, fe.eq3=fe.eq3, fe.btwn=fe.btwn))
  
}

generate_model_plots <- function(dat, x_var, y_var, fe_chron, fe_eq3, fe_btwn) {
  
  ggplot() + 
    geom_point(data = dat,
               inherit.aes = FALSE,
               aes(x = {{ x_var }},
                   y = {{ y_var }}),
               size = 0.1) +
    geom_line(data = fe_btwn,
              inherit.aes = FALSE,
              aes(x = x,
                  y = predicted),
              color = "grey30",
              linewidth = 1) +
    geom_line(data = fe_chron,
              inherit.aes = FALSE,
              aes(x = x,
                  y = predicted),
              color = "steelblue2",
              linewidth = 1) +
    geom_line(data = fe_eq3,
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
          panel.grid.minor = element_line(color = "grey98", linewidth = 0.5))
}

#Import data--------------------------------------------------------------------
hip_flexion<- read_csv("/scratch/ckelsey4/Cayo_meth/hip_flexion.csv", col_names = T)
hip_flexion<- hip_flexion %>%
  group_by(individual_code) %>%
  filter(age > 0) %>%
  mutate(n = n()) %>%
  filter(n > 1) %>%
  mutate(between_age = mean(age),
         within_age = age - between_age) %>%
  relocate(within_age, .after = age) %>%
  relocate(between_age, .after = within_age)

hip_flexion<- hip_flexion %>%
  group_by(individual_code) %>%
  mutate(min_age = min(age))
hip_flexion$age<- round(hip_flexion$age, 0)

samples_dist_hip<- hip_flexion %>%
  ggplot(aes(x=age, y=reorder(individual_code, min_age), colour=as.factor(individual_sex))) +
  geom_path(linewidth = 0.5) +
  geom_point(colour="black", size = 0.25) +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  scale_colour_manual(values = c("red3", "pink2"), name = "Sex") +
  ylab("Individual") +
  xlab("Age") +
  theme_classic(base_size=6) +
  theme(legend.position = "none", 
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = margin(1, 1, 1, 1, "pt"))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/samples_dist_hip.svg", 
         samples_dist_hip, 
         height = 85, width = 55, units = "mm")

hip_flexion %>%
  ggplot(aes(x=age, fill=as.factor(individual_sex))) +
  geom_bar(colour='black', position = 'dodge') +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  scale_fill_manual(values = c("red3", "pink2"), name = "Sex") +
  ylab("Count") +
  xlab("Age") +
  theme_classic(base_size = 24) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(colour = "black", linewidth=3))

n_samples_morph<- hip_flexion %>%
  ggplot(aes(x=n, fill=as.factor(individual_sex))) +
  geom_bar(position = 'dodge') +
  scale_fill_manual(values = c("red3", "pink2"), name = "Sex") +
  ylab("Count") +
  xlab("N Samples") +
  theme_classic(base_size=6) +
  theme(panel.background = element_rect(colour = "black", linewidth=0.5),
        axis.line = element_line(colour = "black", linewidth = 0.25),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        legend.position = "none")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/n_samples_morph.svg", 
       n_samples_morph, 
       height = 20, width = 30, units = "mm")

#Hip Extension------------------------------------------------------------------
#Outliers under 100 degrees are removed
hip_ext<- hip_flexion %>%
  filter(hip_extension_deg > 100)

#Generate models and coefs
hip_extension<- generate_models(hip_ext, hip_extension_deg)

#Generate base plot
hip_ext_plot<- generate_model_plots(dat = hip_ext,
                                    x_var = age,
                                    y_var = hip_extension_deg,
                                    fe_chron = hip_extension$fe.chron,
                                    fe_eq3 = hip_extension$fe.eq3,
                                    fe_btwn = hip_extension$fe.btwn)
hip_ext_plot

#Add text and axis scales to identify plot values for Illustrator
hip_ext_plot<- hip_ext_plot +
  xlab("Age") +
  ylab("Hip Extension (Deg.)") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30)) +
  scale_y_continuous(breaks = seq(100, 180, 20), limits = c(100, 180))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/hip_ext_plot.svg", 
       hip_ext_plot, 
       height = 55, width = 55, units = "mm")

#Generate table
hip_ext_table<- hip_extension[['coefs']]
colnames(hip_ext_table)<- c("β", "se", "df", "t", "p", "2.5%", "97.5%")
hip_ext_table<- hip_ext_table %>%
  mutate(var = rownames(hip_ext_table)) %>%
  separate_wider_delim(cols = var, delim = "_", names = c("Eq", "Var")) %>%
  relocate(Eq, Var, .before = β)
kable(hip_ext_table, format = 'html', align = "c") %>%
  save_kable("hip_ext_table.pdf")

#Hip Internal Rotation----------------------------------------------------------
#Generate models and coefs list
hip_rotation<- generate_models(hip_flexion, hip_internal_rotation_deg)

#Generate base plot
hip_int_plot<- generate_model_plots(dat = hip_flexion,
                                        x_var = age,
                                        y_var = hip_internal_rotation_deg,
                                        fe_chron = hip_rotation$fe.chron,
                                        fe_eq3 = hip_rotation$fe.eq3,
                                        fe_btwn = hip_rotation$fe.btwn)

#Add text and axis scales to identify plot values for Illustrator
hip_int_plot<- hip_int_plot +
  xlab("Age") +
  ylab("Hip Internal Rotation (Deg.)") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30)) +
  scale_y_continuous(breaks = seq(0, 80, 20), limits = c(0, 80))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/hip_int_plot.svg", 
       hip_int_plot, 
       height = 55, width = 55, units = "mm")

#Body Weight--------------------------------------------------------------------
bw<- hip_flexion %>%
  filter(body_weight_lb < 50)

bw_mods<- generate_models(bw, body_weight_lb)

#Generate base plot
bw_plot<- generate_model_plots(dat = bw,
                               x_var = age,
                               y_var = body_weight_lb,
                               fe_chron = bw_mods$fe.chron,
                               fe_eq3 = bw_mods$fe.eq3,
                               fe_btwn = bw_mods$fe.btwn)

bw_plot<- bw_plot +
  xlab("Age") +
  ylab("Body Weight (lbs)") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/bw_plot.svg", 
       bw_plot, 
       height = 55, width = 55, units = "mm")

#Femoral Adduction--------------------------------------------------------------
femoral_adduction<- generate_models(hip_flexion, femoral_adduction_deg)

#### Extra Models ####
#Femoral Abduction--------------------------------------------------------------
femoral_abduction<- generate_models(hip_flexion, femoral_abduction_deg)

#Hip Flexion-------------------------------------------------------------------

#Hip External Rotation----------------------------------------------------------
lowlim<- mean(hip_flexion$hip_external_rotation_deg, na.rm = T) + (2*sd(hip_flexion$hip_external_rotation_deg, na.rm = T))
hip_ext<- hip_flexion %>%
  filter(hip_external_rotation_deg > -1*lowlim & hip_external_rotation_deg < lowlim)

#Generate models and coefs list
hip_rotation<- generate_models(hip_ext, hip_external_rotation_deg)

#Generate base plot
hip_int_plot<- generate_model_plots(dat = hip_ext,
                                    x_var = age,
                                    y_var = hip_external_rotation_deg,
                                    fe_chron = hip_rotation$fe.chron,
                                    fe_eq3 = hip_rotation$fe.eq3,
                                    fe_btwn = hip_rotation$fe.btwn)

hip_int_plot

#Add text and axis scales to identify plot values for Illustrator
hip_int_plot<- hip_int_plot +
  xlab("Age") +
  ylab("Hip external Rotation (Deg.)") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30)) +
  scale_y_continuous(breaks = seq(0, 80, 20), limits = c(0, 80))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/hip_int_plot.svg", 
       hip_int_plot, 
       height = 55, width = 55, units = "mm")

