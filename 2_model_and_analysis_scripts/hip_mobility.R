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

#Import data--------------------------------------------------------------------
hip_flexion<- blood_metadata<- read_csv("/scratch/ckelsey4/Cayo_meth/hip_flexion.csv", col_names = T)
hip_flexion<- hip_flexion %>%
  group_by(individual_code) %>%
  filter(age > 0) %>%
  mutate(n = n()) %>%
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
  geom_path(linewidth = 1.2, alpha = 0.8) +
  geom_point(colour="black", size = 0.5) +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  scale_colour_manual(values = c("green4", "purple4"), name = "Sex") +
  ylab("Individual") +
  xlab("Age") +
  theme_classic(base_size = 24) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/morph_samples.svg", height = 100, width = 87, units = "mm")

hip_flexion %>%
  filter(hip_extension_deg > 50) %>%
  ggplot(aes(age, hip_extension_deg, group = individual_code)) +
  geom_point(alpha = 0.3) +
  geom_path(alpha = 0.5) +
  theme(legend.position = "none")

#Hip Extension------------------------------------------------------------------
ext_range<- hip_flexion %>%
  group_by(individual_code) %>%
  summarize(min = min(hip_extension_deg), 
            max = max(hip_extension_deg))

ext_range<- ext_range %>%
  mutate(diff = min - max)

#Can't remember why I made this filter so I'm taking it out until I figure it out
#hip_ext<- hip_flexion %>%
  #filter(hip_extension_deg > 50)

hip_ext<- hip_flexion

hip_extension_chron<- lmer(hip_extension_deg ~ age + individual_sex + (1|individual_code), 
                           data = hip_ext)

hip_extension_within<- lmer(hip_extension_deg ~ within_age + between_age + individual_sex + (1|individual_code), 
                            data = hip_ext)

hip_extension_eq3<- lmer(hip_extension_deg ~ age + between_age + individual_sex + (1|individual_code), 
                         data = hip_ext)

summary(hip_extension_chron)[["coefficients"]]
summary(hip_extension_eq3)[["coefficients"]]
summary(hip_extension_within)[["coefficients"]]

chron_fm<- predict_response(hip_extension_chron, "age")
chron_re<- predict_response(hip_extension_chron, terms = c("age", "individual_code"), type = "random")

eq2_fm<- predict_response(hip_extension_within, "within_age")
eq2_between_fm<- predict_response(hip_extension_within, "between_age")
eq2_re<- predict_response(hip_extension_within, terms = c("within_age", "individual_code"), type = "random")

eq3_fm <- predict_response(hip_extension_eq3, "age")
eq3_between_fm<- predict_response(hip_extension_eq3, "between_age")
eq3_re<- predict_response(hip_extension_eq3, terms = c("age", "individual_code"), type = "random")

chron_ext<- plot_slopes(chron_fm, chron_re, hip_ext, age, hip_extension_deg)

ggplot() + 
  geom_point(data = hip_ext, 
             inherit.aes=F,
             aes(x = age,
                 y = hip_extension_deg), 
             cex = 1,
             shape = 1,
             alpha = 0.8) +
  geom_smooth(data = chron_fm,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'steelblue2',
              method="lm",
              se = F) +
  geom_smooth(data = eq3_fm,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'purple',
              method="lm",
              se = F) +
  geom_smooth(data = eq3_between_fm,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'chartreuse3',
              method="lm",
              se = F) +
  theme_classic(base_size = 18) + 
  theme(panel.background = element_rect(colour = "black", linewidth=1)) +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.text.y  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.title = element_text(vjust = -5)) +
  xlab("Age") +
  ylab("Hip Extension Deg.") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))
  
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/chron_extension.svg")

eq3_ext<- plot_slopes(eq3_fm, eq3_re, hip_ext, age, hip_extension_deg)
eq3_ext + 
  xlab("Age") +
  ylab("Hip Extension Deg.") +
  scale_x_continuous(breaks = seq(5, 30, 5), limits = c(5, 30))
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq3_extension.svg")

geom_smooth(data = eq3_between_fm,
            inherit.aes=F,
            aes(x = x,
                y = predicted),
            color = 'black',
            alpha = .35,
            method="lm",
            se = F,
            linetype = "dashed") +

eq2_ext_plot<- plot_slopes(eq2_fm, eq2_re, hip_ext, within_age, hip_extension_deg)
eq2_ext_plot +
  xlab("Within-Ind. Centered Age") +
  ylab("Hip Extension Deg.")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/within_extension.svg")

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

#Generate plots
chron_flexion<- plot_slopes(flexion_predicted_chron, re.flexion_predicted_chron, 
                            hip_flexion, age, hip_flexion_deg)
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

#Hip Internal Rotation-------------------------------------------------------------------
hip_int_rotation_within<- lmer(hip_internal_rotation_deg ~ within_age + between_age + individual_sex + (1|individual_code), 
                           data = hip_flexion)

hip_int_rotation_eq3<- lmer(hip_internal_rotation_deg ~ age + between_age + individual_sex + (1|individual_code), 
                        data = hip_flexion)

hip_int_rotation_chron<- lmer(hip_internal_rotation_deg ~ age + individual_sex + (1|individual_code), 
                          data = hip_flexion)

summary(hip_int_rotation_within)[["coefficients"]]
summary(hip_int_rotation_eq3)[["coefficients"]]
summary(hip_int_rotation_chron)[["coefficients"]]

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
              se = F) +
  geom_smooth(data = fm.dat,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'purple',
              method="lm",
              se = F) +
  geom_smooth(data = fm.dat_btwn,
              inherit.aes=F,
              aes(x = x,
                  y = predicted),
              color = 'chartreuse3',
              method="lm",
              se = F) +
  theme_classic(base_size = 18) + 
  theme(panel.background = element_rect(colour = "black", linewidth=1)) +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.text.y  = element_text(vjust=0.5, colour="black")) + 
  theme(axis.title = element_text(vjust = -5)) +
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







