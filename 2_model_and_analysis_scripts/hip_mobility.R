library(tidyverse)
library(ggplot2)
library(ggeffects)
library(lme4)

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
              size = .5) +
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
    theme(axis.text.x  = element_text(vjust=0.5, size=11, colour="black")) + 
    theme(axis.text.y  = element_text(vjust=0.5, size=11, colour="black")) + 
    theme(axis.title = element_text(size=13, vjust = -5))
  
}

#Import data--------------------------------------------------------------------
hip_flexion<- blood_metadata<- read_csv("/scratch/ckelsey4/Cayo_meth/hip_flexion.csv", col_names = T)
hip_flexion<- hip_flexion %>%
  group_by(individual_code) %>%
  mutate(n = n()) %>%
  filter(n > 1) %>%
  mutate(between_age = mean(age),
         within_age = age - between_age) %>%
  relocate(within_age, .after = age) %>%
  relocate(between_age, .after = within_age)

#Hip Extension------------------------------------------------------------------
hip_extension_chron<- lmer(hip_extension_deg ~ age + individual_sex + (1|individual_code), 
                           data = hip_flexion)

hip_extension_within<- lmer(hip_extension_deg ~ within_age + between_age + individual_sex + (1|individual_code), 
                            data = hip_flexion)

hip_extension_eq3<- lmer(hip_extension_deg ~ age + between_age + individual_sex + (1|individual_code), 
                         data = hip_flexion)

summary(hip_extension_chron)[["coefficients"]]
summary(hip_extension_within)[["coefficients"]]
summary(hip_extension_eq3)[["coefficients"]]

fm.dat.chron<- predict_response(hip_extension_chron, "age")
re.dat.chron<- predict_response(hip_extension_chron, terms = c("age", "individual_code"), type = "random")

fm.dat.within<- predict_response(hip_extension_within, "within_age")
re.dat.within<- predict_response(hip_extension_within, terms = c("within_age", "individual_code"), type = "random")

fm.dat <- predict_response(hip_extension_eq3, "age")
re.dat<- predict_response(hip_extension_eq3, terms = c("age", "individual_code"), type = "random")

chron<- plot_slopes(fm.dat.chron, re.dat.chron, hip_flexion, age, hip_extension_deg)
chron +
  xlab("Age") +
  ylab("Hip Extension Deg.") +
  xlim(6, 30)

eq3<- plot_slopes(fm.dat, re.dat, hip_flexion, age, hip_extension_deg)
eq3 + 
  xlab("Age") +
  ylab("Hip Extension Deg.") +
  xlim(6, 30)

within<- plot_slopes(fm.dat.within, re.dat.within, hip_flexion, within_age, hip_extension_deg)
within +
  xlab("Within-Ind. Centered Age") +
  ylab("Hip Extension Deg.")

#Hip Rotation-------------------------------------------------------------------
hip_rotation_within<- lmer(hip_external_rotation_deg ~ within_age + between_age + individual_sex + (1 + within_age|individual_code), 
                           data = hip_flexion)

hip_rotation_eq3<- lmer(hip_external_rotation_deg ~ age + between_age + individual_sex + (1 + age|individual_code), 
                        data = hip_flexion)

hip_rotation_chron<- lmer(hip_external_rotation_deg ~ age + individual_sex + (1|individual_code), 
                          data = hip_flexion)

summary(hip_rotation_within)[["coefficients"]]
summary(hip_rotation_eq3)[["coefficients"]]
summary(hip_rotation_chron)[["coefficients"]]