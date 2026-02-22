library(tidyverse)
library(ggplot2)
library(ggeffects)
library(lme4)
#library(lmerTest)

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

hip_flexion %>%
  ggplot(aes(x=age, y=reorder(individual_code, min_age), colour=as.factor(individual_sex))) +
  geom_path(linewidth = 1.5, alpha = 0.8) +
  geom_point(colour="black", size = 0.5) +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  scale_colour_manual(values = c("red3", "pink2"), name = "Sex") +
  ylab("Individual") +
  xlab("Age") +
  theme_classic(base_size = 24) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(colour = "black", linewidth=3))
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/morph_samples.svg", height = 100, width = 87, units = "mm")

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

hip_flexion %>%
  ggplot(aes(x=n, fill=as.factor(individual_sex))) +
  geom_bar(colour='black', position = 'dodge') +
  scale_fill_manual(values = c("red3", "pink2"), name = "Sex") +
  ylab("Count") +
  xlab("N Samples") +
  theme_classic(base_size = 24) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(colour = "black", linewidth=3))

#Hip Extension------------------------------------------------------------------
#Outliers under 100 degrees are removed
hip_ext<- hip_flexion %>%
  filter(hip_extension_deg > 100)

hip_extension<- generate_models(hip_ext, hip_extension_deg)

ggplot() + 
  geom_point(data = hip_ext, 
             inherit.aes=F,
             aes(x = age,
                 y = hip_extension_deg), 
             cex = 1,
             alpha = 0.8) +
  geom_line(data = hip_extension$fe.chron,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'steelblue2',
            linewidth = 1.5) +
  geom_line(data = hip_extension$fe.eq3,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
            color = 'purple',
            linewidth =1.5) +
  geom_line(data = hip_extension$fe.btwn,
            inherit.aes=F,
            aes(x = x,
                y = predicted),
            color = 'green4',
            linewidth =1.5,
            linetype = "dashed") +
  theme_classic(base_size = 18) + 
  theme(panel.background = element_rect(colour = "black", linewidth=1)) +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.text.y  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.title = element_text(vjust = -5)) +
  theme(panel.background = element_rect(colour = "black", linewidth=3)) +
  xlab("Age") +
  ylab("Hip Extension (Deg.)") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30)) +
  scale_y_continuous(breaks = seq(100, 180, 20), limits = c(100, 180))


#Hip Internal Rotation----------------------------------------------------------
hip_rotation<- generate_models(hip_flexion, hip_internal_rotation_deg)

ggplot() + 
  geom_point(data = hip_flexion, 
             inherit.aes=F,
             aes(x = age,
                 y = hip_internal_rotation_deg), 
             cex = 1,
             alpha = 0.8) +
  geom_line(data = hip_rotation$fe.chron,
            inherit.aes=F,
            aes(x = x,
                y = predicted),
            color = 'steelblue2',
            linewidth = 1.5) +
  geom_line(data = hip_rotation$fe.eq3,
            inherit.aes=F,
            aes(x = x,
                y = predicted),
            color = 'purple',
            linewidth = 1.5) +
  geom_line(data = hip_rotation$fe.btwn,
            inherit.aes=F,
            aes(x = x,
                y = predicted),
            color = 'green4',
            linewidth =1.5,
            linetype = "dashed") +
  theme_classic(base_size = 18) + 
  theme(panel.background = element_rect(colour = "black", linewidth=1)) +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.text.y  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.title = element_text(vjust = -5)) +
  theme(panel.background = element_rect(colour = "black", linewidth=3)) +
  xlab("Age") +
  ylab("Hip Internal Rotation (Deg.)") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30)) +
  scale_y_continuous(breaks = seq(0, 80, 20), limits = c(0, 80))
  

#Femoral Abduction--------------------------------------------------------------
femoral_abduction<- generate_models(hip_flexion, femoral_abduction_deg)

ggplot() + 
  geom_point(data = hip_flexion, 
             inherit.aes=F,
             aes(x = age,
                 y = femoral_abduction_deg), 
             cex = 1,
             shape = 1,
             alpha = 0.8) +
  geom_line(data = femoral_abduction$fe.chron,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'steelblue2',
            linewidth = 1.5) +
  geom_line(data = femoral_abduction$fe.eq3,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'purple',
            linewidth = 1.5) +
  geom_line(data = femoral_abduction$fe.btwn,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
            color = 'green4',
            linewidth = 1.5,
            linetype = "dashed") +
  theme_classic(base_size = 18) + 
  theme(panel.background = element_rect(colour = "black", linewidth=1)) +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.text.y  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.title = element_text(vjust = -5)) +
  theme(panel.background = element_rect(colour = "black", linewidth=3)) +
  xlab("Age") +
  ylab("Femoral Abduction") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))

#Femoral Adduction--------------------------------------------------------------
femoral_adduction<- generate_models(hip_flexion, femoral_adduction_deg)

ggplot() + 
  geom_point(data = hip_flexion, 
             inherit.aes=F,
             aes(x = age,
                 y = femoral_adduction_deg), 
             cex = 1,
             shape = 1,
             alpha = 0.8) +
  geom_line(data = femoral_abduction$fe.chron,
            inherit.aes=F,
            aes(x = x,
                y = predicted),
            color = 'steelblue2',
            linewidth = 1.5) +
  geom_line(data = femoral_adduction$fe.eq3,
            inherit.aes=F,
            aes(x = x,
                y = predicted),
            color = 'purple',
            linewidth = 1.5) +
  geom_line(data = femoral_adduction$fe.btwn,
            inherit.aes=F,
            aes(x = x,
                y = predicted),
            color = 'green4',
            linewidth = 1.5,
            linetype = "dashed") +
  theme_classic(base_size = 18) + 
  theme(panel.background = element_rect(colour = "black", linewidth=1)) +
  theme(axis.text.x  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.text.y  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.title = element_text(vjust = -5)) +
  theme(panel.background = element_rect(colour = "black", linewidth=3)) +
  xlab("Age") +
  ylab("Femoral Adduction") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))

#Body Weight--------------------------------------------------------------------
bw<- hip_flexion %>%
  filter(body_weight_lb < 50)

bw_mods<- generate_models(bw, body_weight_lb)

ggplot() + 
  geom_point(data = bw, 
             inherit.aes=F,
             aes(x = age,
                 y = body_weight_lb), 
             cex = 1,
             alpha = 0.8) +
  geom_line(data = bw_mods$fe.chron,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'steelblue2',
              linewidth = 1.5) +
  geom_line(data = bw_mods$fe.eq3,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'purple',
              linewidth = 1.5) +
  geom_line(data = bw_mods$fe.btwn,
            inherit.aes=F,
            aes(x = x,
                y = predicted),
            color = 'green4',
            linewidth = 1.5,
            linetype = "dashed") +
  theme_classic(base_size = 18) + 
  theme(panel.background = element_rect(colour = "black", linewidth=1)) +
  theme(axis.text.x  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.text.y  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.title = element_text(vjust = -5)) +
  theme(panel.background = element_rect(colour = "black", linewidth=3)) +
  xlab("Age") +
  ylab("Body Weight (lbs)") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30)) +
  scale_y_continuous(breaks = seq(10, 35, 5), limits = c(10, 35))

#### Flexion and Internal Rotation ####
#Hip Flexion-------------------------------------------------------------------
hip_flexion_chron<- lmer(hip_flexion_deg ~ age + individual_sex + (1|individual_code), 
                         data = hip_flexion)

hip_flexion_eq3<- lmer(hip_flexion_deg ~ age + between_age + individual_sex + (1|individual_code), 
                       data = hip_flexion)

hip_flexion_within<- lmer(hip_flexion_deg ~ within_age + between_age + individual_sex + (1|individual_code), 
                          data = hip_flexion)

summary(hip_flexion_chron)[["coefficients"]]
summary(hip_flexion_eq3)[["coefficients"]]
summary(hip_flexion_within)[["coefficients"]]

#Generate predicted responses
flexion_predicted_chron<- predict_response(hip_flexion_chron, "age")
re.flexion_predicted_chron<- predict_response(hip_flexion_chron, 
                                              terms = c("age", "individual_code"), type = "random")

flexion_predicted_eq3<- predict_response(hip_flexion_eq3, "age")
re.flexion_predicted_eq3<- predict_response(hip_flexion_eq3, 
                                            terms = c("age", "individual_code"), type = "random")

flexion_predicted_within<- predict_response(hip_flexion_within, "within_age")
re.flexion_predicted_within<- predict_response(hip_flexion_within, 
                                               terms = c("within_age", "individual_code"), type = "random")

ggplot() + 
  geom_point(data = hip_flexion, 
             inherit.aes=F,
             aes(x = age,
                 y = hip_flexion_deg), 
             cex = 1,
             shape = 1,
             alpha = 0.8) +
  geom_smooth(data = flexion_predicted_chron,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'steelblue2',
              method="lm",
              se = F) +
  geom_smooth(data = flexion_predicted_eq3,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'purple',
              method="lm",
              se = F) +
  geom_smooth(data = flexion_predicted_eq3,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'green4',
              method="lm",
              se = F,
              linetype = "dashed") +
  theme_classic(base_size = 18) + 
  theme(panel.background = element_rect(colour = "black", linewidth=1)) +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.text.y  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.title = element_text(vjust = -5)) +
  theme(panel.background = element_rect(colour = "black", linewidth=3)) +
  xlab("Age") +
  ylab("Hip Flexion Deg.") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))

#Generate plots
chron_flexion<- plot_slopes(flexion_predicted_chron, re.flexion_predicted_chron, 
                            hip_flexion, age, hip_flexion_deg)

ggplot() + 
  geom_point(data = hip_flexion, 
             inherit.aes=F,
             aes(x = age,
                 y = hip_flexion_deg), 
             cex = 1,
             shape = 1,
             alpha = 0.8) +
  theme_classic(base_size = 18) + 
  theme(panel.background = element_rect(colour = "black", linewidth=1)) +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.text.y  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.title = element_text(vjust = -5)) +
  theme(panel.background = element_rect(colour = "black", linewidth=3)) +
  xlab("Age") +
  ylab("Hip Flexion Deg.") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))

chron_flexion +
  xlab("Age") +
  ylab("Hip Flexion Deg.") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))

eq3_flexion<- plot_slopes(flexion_predicted_eq3, re.flexion_predicted_eq3, 
                          hip_flexion, age, hip_flexion_deg)
eq3_flexion + 
  xlab("Age") +
  ylab("Hip Flexion Deg.") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30)) +
  ylim(0, 100)

within_flexion<- plot_slopes(flexion_predicted_within, re.flexion_predicted_within, 
                             hip_flexion, within_age, hip_flexion_deg)
within_flexion +
  xlab("Within-Ind. Centered Age") +
  ylab("Hip Flexion Deg.")

#External Hip Rotation----------------------------------------------------------
hip_ext_rotation_within<- lmer(hip_external_rotation_deg ~ within_age + between_age + individual_sex + (1|individual_code), 
                               data = hip_flexion)

hip_ext_rotation_eq3<- lmer(hip_external_rotation_deg ~ age + between_age + individual_sex + (1|individual_code), 
                            data = hip_flexion)

hip_ext_rotation_chron<- lmer(hip_external_rotation_deg ~ age + individual_sex + (1|individual_code), 
                              data = hip_flexion)

summary(hip_ext_rotation_within)[["coefficients"]]
summary(hip_ext_rotation_eq3)[["coefficients"]]
summary(hip_ext_rotation_chron)[["coefficients"]]

fm.dat.chron<- predict_response(hip_ext_rotation_chron, "age")
re.dat.chron<- predict_response(hip_ext_rotation_chron, terms = c("age", "individual_code"), type = "random")

fm.dat.within<- predict_response(hip_ext_rotation_within, "within_age")
re.dat.within<- predict_response(hip_ext_rotation_within, terms = c("within_age", "individual_code"), type = "random")

fm.dat <- predict_response(hip_ext_rotation_eq3, "age")
re.dat<- predict_response(hip_ext_rotation_eq3, terms = c("age", "individual_code"), type = "random")

chron_ext_rotation<- plot_slopes(fm.dat.chron, re.dat.chron, hip_flexion, age, hip_external_rotation_deg)
chron_ext_rotation +
  xlab("Age") +
  ylab("Hip External Rotation Deg.") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))

eq3_ext_rotation<- plot_slopes(fm.dat, re.dat, hip_flexion, age, hip_external_rotation_deg)
eq3_ext_rotation + 
  xlab("Age") +
  ylab("Hip External Rotation Deg.") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))

within_ext_rotation<- plot_slopes(fm.dat.within, re.dat.within, hip_flexion, within_age, hip_external_rotation_deg)
within_ext_rotation +
  xlab("Within-Ind. Centered Age") +
  ylab("Hip External Rotation Deg.")

