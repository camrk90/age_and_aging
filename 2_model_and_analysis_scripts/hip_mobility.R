library(tidyverse)
library(ggplot2)
library(ggeffects)
library(lme4)
#library(lmerTest)

#Define functions---------------------------------------------------------------
plot_slopes<- function(fixed, re, metadata, predictor, response){
  
  ggplot() + 
    geom_ribbon(data = fixed, 
                inherit.aes = F,
                aes(x = x,
                    y = predicted,
                    ymin = conf.low,
                    ymax = conf.high),
                fill = "grey90")  +
    geom_line(data = re, 
              aes(x = x, 
                  y = predicted, 
                  colour = group),
              alpha = .2,
              linewidth = .5) +
    geom_smooth(data = fixed,
                inherit.aes=F,
                aes(x = x,
                    y = predicted),
                color = 'black',
                alpha = .35,
                method="lm",
                se = F) +
    geom_point(data = metadata, 
               inherit.aes=F,
               aes(x = {{predictor}},
                   y = {{response}}), 
               cex = 1,
               shape = 1,
               alpha = 0.4) +
    theme_classic() + 
    theme(panel.background = element_rect(colour = "black", linewidth=1)) +
    theme(legend.position = "none") +
    theme(axis.text.x  = element_text(vjust=0.5, size=7, colour="black")) + 
    theme(axis.text.y  = element_text(vjust=0.5, size=7, colour="black")) + 
    theme(axis.title = element_text(size=7, vjust = -5))
  
}
generate_predictions<- function(chron_mod, eq2_mod, eq3_mod){
  
  coefs<- as.data.frame(rbind(summary(chron_mod)[["coefficients"]], summary(eq2_mod)[["coefficients"]],
                                    summary(eq3_mod)[["coefficients"]]))
  
  rownames(coefs)<- c("eq1_intercept", "eq1_age", "eq1_sexM",
                      "eq2_intercept", "eq2_age.w", "eq2_age.btwn", "eq2_sexM",
                      "eq3_intercept", "eq3_age.w", "eq3_age.btwn", "eq3_sexM")
  
  conf_ints<- as.data.frame(rbind(confint(chron_mod), 
                                  confint(eq2_mod),
                                  confint(eq3_mod)))
  
  conf_ints<- conf_ints[!grepl("sig", rownames(conf_ints)),]
  
  coefs<- cbind(coefs, conf_ints)
  
  fe.chron<- predict_response(chron_mod, "age")
  #re.chron<- predict_response(chron_mod, terms = c("age", "individual_code"), type = "random")
  
  fe.eq2<- predict_response(eq2_mod, "within_age")
  #re.eq2<- predict_response(eq2_mod, terms = c("within_age", "individual_code"), type = "random")
  
  fe.eq3 <- predict_response(eq3_mod, "age")
  #re.eq3<- predict_response(eq3_mod, terms = c("age", "individual_code"), type = "random")
  
  return(list(coefs = coefs, fe.chron=fe.chron, fe.eq2=fe.eq2, fe.eq3=fe.eq3))
  
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

hip_flexion %>%
  filter(hip_extension_deg > 100) %>%
  ggplot(aes(age, hip_extension_deg, colour = individual_code)) +
  geom_point(alpha = 0.3) +
  geom_path(alpha = 0.5) +
  theme(legend.position = "none")

#Hip Extension------------------------------------------------------------------
#Outliers under 100 degrees are removed
hip_ext<- hip_flexion %>%
  filter(hip_extension_deg > 100) 

hip_extension_chron<- lmer(hip_extension_deg ~ age + individual_sex + (1|individual_code), 
                           data = hip_ext)

hip_extension_within<- lmer(hip_extension_deg ~ within_age + between_age + individual_sex + (1|individual_code), 
                            data = hip_ext)

hip_extension_eq3<- lmer(hip_extension_deg ~ age + between_age + individual_sex + (1|individual_code), 
                         data = hip_ext)

hip_extension<- generate_predictions(chron_mod = hip_extension_chron, 
                                     eq2_mod = hip_extension_within, 
                                     eq3_mod = hip_extension_eq3)

df2<- hip_ext %>%
  select(individual_code, individual_sex, between_age, hip_extension_deg) %>%
  group_by(individual_code) %>%
  mutate(mean_hip_ext = mean(hip_extension_deg)) %>%
  ungroup() %>%
  distinct(individual_code, .keep_all = T)
df2$age.w_slope<- hip_extension$coefs[5,1]

x_range <- diff(range(df2$between_age))
half_dx <- 0.05 * x_range

# compute segment endpoints for each individual's slope segment
df2 <- df2 %>%
  mutate(
    dx = half_dx,
    x0 = between_age - dx,
    x1 = between_age + dx,
    y0 = mean_hip_ext - age.w_slope * dx,
    y1 = mean_hip_ext + age.w_slope * dx
  )

df2 %>%
  ggplot(aes(x=between_age, y=mean_hip_ext)) + 
  geom_segment(aes(x = x0, xend = x1, y = y0, yend = y1),
               colour = 'green4', 
               size = 0.8, lineend = "round", show.legend = TRUE) +
  geom_line(data = hip_extension$fe.chron,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
            color = 'steelblue2',
            linewidth = 1.5) +
  geom_point(cex = 1,
             alpha = 0.8) +
  theme_classic(base_size = 18) + 
  theme(panel.background = element_rect(colour = "black", linewidth=1)) +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.text.y  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.title = element_text(vjust = -5)) +
  theme(panel.background = element_rect(colour = "black", linewidth=3)) +
  xlab("Age") +
  ylab("Hip Extension (Deg.)") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))

ggplot() + 
  geom_point(data = hip_ext, 
             inherit.aes=F,
             aes(x = age,
                 y = hip_extension_deg), 
             cex = 1,
             shape = 1,
             alpha = 0.8) +
  geom_line(data = hip_extension$fe.chron,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'steelblue2') +
  geom_line(data = hip_extension$fe.eq3,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
            color = 'purple',
            linewidth =1.5) +
  geom_smooth(data = hip_extension$fe.eq3,
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
  ylab("Hip Extension (Deg.)") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))
  


#Hip Internal Rotation-------------------------------------------------------------------
hip_int_rotation_within<- lmer(hip_internal_rotation_deg ~ within_age + between_age + individual_sex + (1|individual_code), 
                           data = hip_flexion)

hip_int_rotation_eq3<- lmer(hip_internal_rotation_deg ~ age + between_age + individual_sex + (1|individual_code), 
                        data = hip_flexion)

hip_int_rotation_chron<- lmer(hip_internal_rotation_deg ~ age + individual_sex + (1|individual_code), 
                          data = hip_flexion)

summary(hip_int_rotation_chron)[["coefficients"]]
confint(hip_int_rotation_chron)
summary(hip_int_rotation_within)[["coefficients"]]
confint(hip_int_rotation_within)
summary(hip_int_rotation_eq3)[["coefficients"]]
confint(hip_int_rotation_eq3)

int_coefs<- as.data.frame(rbind(summary(hip_int_rotation_chron)[["coefficients"]], summary(hip_int_rotation_within)[["coefficients"]],
                                summary(hip_int_rotation_eq3)[["coefficients"]]))

int_coefs<- int_coefs[c("age", "within_age", "between_age", "age.1", "between_age.1"),]

fm.dat.chron<- predict_response(hip_int_rotation_chron, "age")
re.dat.chron<- predict_response(hip_int_rotation_chron, terms = c("age", "individual_code"), type = "random")

fm.dat.within<- predict_response(hip_int_rotation_within, "within_age")
re.dat.within<- predict_response(hip_int_rotation_within, terms = c("within_age", "individual_code"), type = "random")

fm.dat <- predict_response(hip_int_rotation_eq3, "age")
fm.dat_btwn <- predict_response(hip_int_rotation_eq3, "between_age")
re.dat<- predict_response(hip_int_rotation_eq3, terms = c("age", "individual_code"), type = "random")

ggplot() + 
  geom_point(data = hip_ext, 
             inherit.aes=F,
             aes(x = age,
                 y = hip_internal_rotation_deg), 
             cex = 1,
             shape = 1,
             alpha = 0.8) +
  geom_smooth(data = fm.dat.chron,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'steelblue2',
              method="lm",
              se = F,
              linewidth = 3) +
  theme_classic(base_size = 18) + 
  theme(panel.background = element_rect(colour = "black", linewidth=1)) +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.text.y  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.title = element_text(vjust = -5)) +
  theme(panel.background = element_rect(colour = "black", linewidth=3)) +
  xlab("Age") +
  ylab("Hip Internal Rotation Deg.") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))

chron_int_rot<- plot_slopes(fm.dat.chron, re.dat.chron, hip_flexion, age, hip_internal_rotation_deg)
chron_int_rot +
  xlab("Age") +
  ylab("Hip Internal Rotation Deg.") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))

eq3_int_rot<- plot_slopes(fm.dat, re.dat, hip_flexion, age, hip_internal_rotation_deg)
eq3_int_rot + 
  xlab("Age") +
  ylab("Hip Internal Rotation Deg.") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30)) 

within_int_rot<- plot_slopes(fm.dat.within, re.dat.within, hip_flexion, within_age, hip_internal_rotation_deg)
within_int_rot +
  xlab("Within-Ind. Centered Age") +
  ylab("Hip Internal Rotation Deg.")















#Femoral Abduction--------------------------------------------------------------
fem_abduction_within<- lmer(femoral_abduction_deg ~ within_age + between_age + individual_sex + (1|individual_code), 
                               data = hip_flexion)

fem_abduction_eq3<- lmer(femoral_abduction_deg ~ age + between_age + individual_sex + (1|individual_code), 
                            data = hip_flexion)

fem_abduction_chron<- lmer(femoral_abduction_deg ~ age + individual_sex + (1|individual_code), 
                              data = hip_flexion)

summary(fem_abduction_chron)[["coefficients"]]
confint(fem_abduction_chron)
summary(fem_abduction_within)[["coefficients"]]
confint(fem_abduction_within)
summary(fem_abduction_eq3)[["coefficients"]]
confint(fem_abduction_eq3)

abduc_coefs<- as.data.frame(rbind(summary(fem_abduction_chron)[["coefficients"]], summary(fem_abduction_within)[["coefficients"]],
                                  summary(fem_abduction_eq3)[["coefficients"]]))

abduc_coefs<- abduc_coefs[c("age", "within_age", "between_age", "age.1", "between_age.1"),]

fm.dat.chron<- predict_response(fem_abduction_chron, "age")
re.dat.chron<- predict_response(fem_abduction_chron, terms = c("age", "individual_code"), type = "random")

fm.dat.within<- predict_response(fem_abduction_within, "within_age")
re.dat.within<- predict_response(fem_abduction_within, terms = c("within_age", "individual_code"), type = "random")

fm.dat <- predict_response(fem_abduction_eq3, "age")
fm.dat_btwn <- predict_response(fem_abduction_eq3, "between_age")
re.dat<- predict_response(fem_abduction_eq3, terms = c("age", "individual_code"), type = "random")

ggplot() + 
  geom_point(data = hip_flexion, 
             inherit.aes=F,
             aes(x = age,
                 y = femoral_abduction_deg), 
             cex = 1,
             shape = 1,
             alpha = 0.8) +
  geom_smooth(data = fm.dat.chron,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'steelblue2',
              method="lm",
              se = F) +
  geom_smooth(data = fm.dat,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'purple',
              method="lm",
              se = F) +
  geom_smooth(data = fm.dat,
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
  ylab("Femoral Abduction") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))

#Femoral Abduction--------------------------------------------------------------
femoral_adduction_within<- lmer(femoral_adduction_deg ~ within_age + between_age + individual_sex + (1|individual_code), 
                            data = hip_flexion)

femoral_adduction_eq3<- lmer(femoral_adduction_deg ~ age + between_age + individual_sex + (1|individual_code), 
                         data = hip_flexion)

femoral_adduction_chron<- lmer(femoral_adduction_deg ~ age + individual_sex + (1|individual_code), 
                           data = hip_flexion)

summary(femoral_adduction_chron)[["coefficients"]]
confint(femoral_adduction_chron)
summary(femoral_adduction_within)[["coefficients"]]
confint(femoral_adduction_within)
summary(femoral_adduction_eq3)[["coefficients"]]
confint(femoral_adduction_eq3)

abduc_coefs<- as.data.frame(rbind(summary(femoral_adduction_chron)[["coefficients"]], summary(femoral_adduction_within)[["coefficients"]],
                                  summary(femoral_adduction_eq3)[["coefficients"]]))

abduc_coefs<- abduc_coefs[c("age", "within_age", "between_age", "age.1", "between_age.1"),]

fm.dat.chron<- predict_response(femoral_adduction_chron, "age")
re.dat.chron<- predict_response(femoral_adduction_chron, terms = c("age", "individual_code"), type = "random")

fm.dat.within<- predict_response(femoral_adduction_within, "within_age")
re.dat.within<- predict_response(femoral_adduction_within, terms = c("within_age", "individual_code"), type = "random")

fm.dat <- predict_response(femoral_adduction_eq3, "age")
fm.dat_btwn <- predict_response(femoral_adduction_eq3, "between_age")
re.dat<- predict_response(femoral_adduction_eq3, terms = c("age", "individual_code"), type = "random")

ggplot() + 
  geom_point(data = hip_flexion, 
             inherit.aes=F,
             aes(x = age,
                 y = femoral_adduction_deg), 
             cex = 1,
             shape = 1,
             alpha = 0.8) +
  geom_smooth(data = fm.dat.chron,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'steelblue2',
              method="lm",
              se = F) +
  geom_smooth(data = fm.dat,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'purple',
              method="lm",
              se = F) +
  geom_smooth(data = fm.dat,
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
  ylab("Femoral Adduction (Deg.)") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))

#Body Weight--------------------------------------------------------------------
bw<- hip_flexion %>%
  filter(body_weight_lb < 50)

body_weight_within<- lmer(body_weight_lb ~ within_age + between_age + individual_sex + (1|individual_code), 
                                data = bw)

body_weight_eq3<- lmer(body_weight_lb ~ age + between_age + individual_sex + (1|individual_code), 
                             data = bw)

body_weight_chron<- lmer(body_weight_lb ~ age + individual_sex + (1|individual_code), 
                               data = bw)

summary(body_weight_chron)[["coefficients"]]
confint(body_weight_chron)
summary(body_weight_within)[["coefficients"]]
confint(body_weight_within)
summary(body_weight_eq3)[["coefficients"]]
confint(body_weight_eq3)

bw_coefs<- as.data.frame(rbind(summary(body_weight_chron)[["coefficients"]], summary(body_weight_within)[["coefficients"]],
                                  summary(body_weight_eq3)[["coefficients"]]))

bw_coefs<- abduc_coefs[c("age", "within_age", "between_age", "age.1", "between_age.1"),]

fm.dat.chron<- predict_response(body_weight_chron, c("age", "individual_sex"))
re.dat.chron<- predict_response(body_weight_chron, terms = c("age", "individual_code"), type = "random")

fm.dat.within<- predict_response(body_weight_within, c("within_age", "individual_sex"))
re.dat.within<- predict_response(body_weight_within, terms = c("within_age", "individual_code"), type = "random")

fm.dat <- predict_response(body_weight_eq3, c("age", "individual_sex"))
fm.dat2 <- predict_response(body_weight_eq3, c("age", "individual_sex"))
fm.dat_btwn <- predict_response(body_weight_eq3, "between_age")
re.dat<- predict_response(body_weight_eq3, terms = c("age", "individual_code"), type = "random")

ggplot() + 
  geom_point(data = hip_flexion, 
             inherit.aes=F,
             aes(x = age,
                 y = femoral_adduction_deg), 
             cex = 1,
             shape = 1,
             alpha = 0.8) +
  geom_smooth(data = fm.dat.chron,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'steelblue2',
              method="lm",
              se = F) +
  geom_smooth(data = fm.dat,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'purple',
              method="lm",
              se = F) +
  geom_smooth(data = fm.dat,
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
  ylab("Body Weight (lbs)") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))

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

