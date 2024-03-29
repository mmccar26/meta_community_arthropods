##Meta-community analysis of flying arthropods
##Bertellotti, Sommer, Schmitz, and McCary
#1 March 2023

#===load libraries======
library(ggplot2)
library(plyr)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(tidyverse)
library(Rmisc)
library(broom)
library(emmeans)
library(vegan)

#=====import data======
#relative pathname
sticky <- file.path(".", "Data", "Sticky_traps_Yale_Myers_Final.csv")
print(sticky)

#import data
stick<- read_csv(sticky)%>%
  mutate(Block = as.factor(Block),
         Trap = as.factor(Trap),
         Year = as.factor(Year))%>%
  drop_na()%>%
  group_by(Treatment, Block, Year)%>%
  summarise(across(Diptera:Total, ~ mean(.x, na.rm = TRUE)))%>%
  ungroup()

###Year x treatment analysis
##Linear mixed effects models for year x treatment for flying arthropods

##histogram not transformed
hist(stick$Total) ## needs a log transformation

##Log transformation
stick$Total<-log(stick$Total)

##Check histogram
hist(stick$Total, prob = TRUE, ylim=c(0, 2)) 
##data looks more close enough to normal, check histogram with overlay

##Histogram with normal distribution
g = stick$Total
m<-mean(g)
std<-sqrt(var(g))
hist(g, prob=TRUE, 
     xlab="Total arthropods", ylim=c(0, 2), main = NULL)
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n") 

#Run LMMs
mod0<-lmer(log(Total)~ Treatment * Year + (1|Block), data = stick)
anova(mod0, F = "Kenward-Rogers") # mirrors ground arthropod analysis; not treatment x year interaction

##Visual of model fits
##Model diagnostics (qq norm and residual checks)
par(mfrow=c(1,2)) ## set the plot matrix

##qq norm plots
qqPlot(resid(mod0), xlab="Theoretical Quantiles", ylab = "Sample Quantiles", line = "quartiles", col.lines = "black", grid = FALSE)

##residual vs fitted plot
plot(fitted(mod0), resid(mod0), xlab = "Fitted Residuals", ylab = "Residuals") #residuals vs fitted
abline(h=0)

#******
##Linear mixed effects models for both years (Simpson Diversity) -- no documented effects for flying arthropods

##to calculate Simpson Diversity
diversity_s <- as.data.frame(cbind(stick[1:3],diversity(index = "simpson", stick[4:12])))

#rename diversity variable
colnames(diversity_s)[4] ="Diversity"

##LMM to diversity treatment effects on diversity
d_simpson<-lmer(Diversity ~ Treatment * Year + (1|Block), data = diversity_s)
anova(d_simpson, F = "Kenward-Rogers")

#******
##analysis for each conducted separately for each taxonomic group for flying arthropods

#===================================================================
##Year 2020 (initial conditions of the study)
#===================================================================

#subset 2020 data
year2020<-
  stick%>%
  filter(Year == 2020)

##LMM for the effects of treatment on flying arthropods in 2020
model1<-lmer(log(Total) ~ Treatment + (1|Block), data = year2020)
anova(model1, F = "Kenward-Rogers")

model2<-lmer(log(Arthropod.other+1) ~ Treatment + (1|Block), data = year2020)
anova(model2, F = "Kenward-Rogers")

model3<-lmer(log(Diptera+1)~ Treatment + (1|Block), data = year2020)
anova(model3, F = "Kenward-Rogers")

model4<-lmer(log(Auchenorrhyncha+1) ~ Treatment+ (1|Block), data = year2020)
anova(model4, F = "Kenward-Rogers")

model5<-lmer(log(Coleoptera+1)~ Treatment+ (1|Block), data = year2020)
anova(model5, F = "Kenward-Rogers")

model6<-lmer(log(Hemiptera+1)~ Treatment+ (1|Block), data = year2020)
anova(model6, F = "Kenward-Rogers")

model7<-lmer(log(Hymenoptera+1)~ Treatment+ (1|Block), data = year2020)
anova(model7, F = "Kenward-Rogers")

model8<-lmer(log(Lepidoptera+1)~ Treatment + (1|Block), data = year2020)
anova(model8)

model9<-lmer(log(Mecoptera+1)~ Treatment + (1|Block), data = year2020)
anova(model9, F = "Kenward-Rogers")

model10<-lmer(log(Heteroptera+1)~ Treatment + (1|Block), data = year2020)
anova(model10, F = "Kenward-Rogers")

#plot sticky traps
##to plot by each individual taxon
flying<-
  year2020%>%
  pivot_longer(Diptera:Total, names_to = "Taxon", values_to = "Abundance")%>%
  mutate(Treatment= revalue(Treatment, c("Full.connection" = "Full", 
                                         "Partial.connection" = "Partial",
                                         "No.connection" = "None")))

flying%>%
  mutate(name = fct_relevel(Treatment, 
                            "Control", "Full", "Partial", "None"))%>%
  ggplot((aes(x = name, y = Abundance, color = Treatment))) +
  scale_fill_manual(values=c("gray", "white", "red", "blue"))+
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  geom_point(position = "jitter", size = 3)+
  stat_summary(fun.y = "mean", color = "black", shape = 18, size = 1)+
  xlab(NULL) +
  ylab("Activity-density") +
  scale_color_manual(values=c("darkgray", "black", "red", "blue"))+
  facet_wrap(~ Taxon, scales = "free")+
  theme(axis.text.x = element_text(colour= "black", face = "bold", size = 15),
        axis.text.y = element_text(colour= "black", face = "bold", size = 15),
        axis.line = element_line(colour = "black", size = .3),
        axis.line.x = element_line(colour = "black", size =.3),
        axis.ticks.x = element_line(colour = "black", size = 1),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size =.3),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size = 20, color = 'black', face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill= "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))

##diversity
#******
##Linear mixed effects models for 2020 (Simpson Diversity) -- no documented effects for flying arthropods

##to calculate Simpson Diversity
diversity2020 <- as.data.frame(cbind(year2020[1:3],diversity(index = "simpson", year2020[4:12])))

#rename diversity variable
colnames(diversity2020)[4] ="Diversity"

##LMM to diversity treatment effects on diversity
mod.d20<-lmer(Diversity ~ Treatment + (1|Block), data = diversity2020)
anova(mod.d20, F = "Kenward-Rogers")



#===================================================================
##Year 2021--the main results of the paper
#===================================================================

##subset by year 2021
year2021<-
  stick%>%
  filter(Year == 2021)

##LMM to test for the effects of treatment on flying arthropod abundance
model11<-lmer(log(Total+1)~ Treatment + (1|Block), data = year2021)
anova(model11, F = "Kenward-Rogers")

model12<-lmer(log(Arthropod.other+1) ~ Treatment + (1|Block), data = year2021)
anova(model12, F = "Kenward-Rogers")

model13<-lmer(log(Diptera+1)~ Treatment+ (1|Block), data = year2021)
anova(model13, F = "Kenward-Rogers")

model14<-lmer(log(Auchenorrhyncha+1)~ Treatment+ (1|Block), data = year2021)
anova(model14, F = "Kenward-Rogers")

model15<-lmer(log(Coleoptera+1)~ Treatment+ (1|Block), data = year2021)
anova(model15, F = "Kenward-Rogers")

model16<-lmer(log(Hemiptera+1)~ Treatment+ (1|Block), data = year2021)
anova(model16, F = "Kenward-Rogers")

model17<-lmer(log(Hymenoptera+1)~ Treatment+ (1|Block), data = year2021)
anova(model17, F = "Kenward-Rogers")

model18<-lmer(log(Lepidoptera+1)~ Treatment+ (1|Block), data = year2021)
anova(model18, F = "Kenward-Rogers")

model19<-lmer(log(Mecoptera+1)~ Treatment+ (1|Block), data = year2021)
anova(model19, F = "Kenward-Rogers")

model20<-lmer(log(Heteroptera+1)~ Treatment+ (1|Block), data = year2021)
anova(model20, F = "Kenward-Rogers")

#plot sticky traps
##to plot by each individual taxon
flying<-
  year2021%>%
  pivot_longer(Diptera:Total, names_to = "Taxon", values_to = "Abundance")%>%
  mutate(Treatment= revalue(Treatment, c("Full.connection" = "Full", 
                                         "Partial.connection" = "Partial",
                                         "No.connection" = "None")))

flying%>%
  mutate(name = fct_relevel(Treatment, 
                            "Control", "Full", "Partial", "None"))%>%
  ggplot((aes(x = name, y = Abundance, color = Treatment))) +
  scale_fill_manual(values=c("gray", "white", "red", "blue"))+
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  geom_point(position = "jitter", size = 3)+
  stat_summary(fun.y = "mean", color = "black", shape = 18, size = 1)+
  xlab(NULL) +
  ylab("Activity-density") +
  scale_color_manual(values=c("darkgray", "black", "red", "blue"))+
  facet_wrap(~ Taxon, scales = "free")+
  theme(axis.text.x = element_text(colour= "black", face = "bold", size = 15),
        axis.text.y = element_text(colour= "black", face = "bold", size = 15),
        axis.line = element_line(colour = "black", size = .3),
        axis.line.x = element_line(colour = "black", size =.3),
        axis.ticks.x = element_line(colour = "black", size = 1),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size =.3),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size = 20, color = 'black', face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill= "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))

##diversity
#******
##Linear mixed effects models for 2021 (Simpson Diversity) -- no documented effects

##to calculate Simpson Diversity
diversity2021 <- as.data.frame(cbind(year2021[1:3],diversity(index = "simpson", year2021[4:12])))

#rename diversity variable
colnames(diversity2021)[4] ="Diversity"

##LMM to diversity treatment effects on diversity
mod.d21<-lmer(Diversity ~ Treatment + (1|Block), data = diversity2021)
anova(mod.d21, F = "Kenward-Rogers")
