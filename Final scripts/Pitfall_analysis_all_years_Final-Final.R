##Meta-community analysis of ground arthropods
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
library(corrplot)
library(tidyr)
library(tidySEM)
library(piecewiseSEM)

#=====import data======
#relative pathname
pit <- file.path(".", "Data", "Pitfall_data_Yale_Myers_Final.csv")
print(pit)

#import data
pitfall<- read_csv(pit)%>%
  mutate(Block = as.factor(Block),
         Plot = as.factor(Plot),
         Year = as.factor(Year))%>%
  group_by(Treatment, Block, Year, Connectivity)%>%
  summarise(across(Entomobryidae:Total, ~ mean(.x, na.rm = TRUE)))%>%
  ungroup()

###Year x treatment analysis
##summarize arthropods into functional guilds
fun.groups<-
  pitfall%>%
  mutate(Fungivores = Entomobryidae+Isotomidae+Sminthuridae+Tomoceridae+Acari,
         Predators = Lycosidae+Thomisidae+Linyphiidae+Salticidae+Gnaphosidae+Araneae.other+Staphylinidae+Carabidae+Chilopoda+Opiliones,
         Detritivores = Diplopoda+Isopoda,
         Collembola = Entomobryidae+Isotomidae+Sminthuridae+Tomoceridae,
         Araneae = Lycosidae+Thomisidae+Linyphiidae+Salticidae+Gnaphosidae+Araneae.other,
         Coleoptera = Circulionidae+Staphylinidae+Carabidae+Coleoptera.other)%>%
  select(Treatment:Connectivity, Acari, Araneae, Coleoptera, Collembola, Diplopoda, Formicidae, Isopoda, Total)

##Linear mixed effects models for year x treatment
##Test of total arthropods

##histogram not transformed
hist(fun.groups$Total) ## Data appears to be a Poisson, but I'll try a log trans

##Log transformation
fun.groups$Total<-log(fun.groups$Total)

##Check histogram
hist(fun.groups$Total, prob = TRUE, ylim=c(0, .8)) 
##data looks more close enough to normal, check histogram with overlay

##Histogram with normal distribution
g = fun.groups$Total
m<-mean(g)
std<-sqrt(var(g))
hist(g, prob=TRUE, 
     xlab="Total arthropods", ylim=c(0, 0.6), main = NULL)
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n") 

#Run LMMs on total arthropod abundance
mod0<-lmer(log(Total)~ Treatment * Year + (1|Block), data = fun.groups)
anova(mod0, F = "Kenward-Rogers") #a weak year x treatment interaction

##Visual of model fits
##Model diagnostics (qq norm and residual checks)
par(mfrow=c(1,2)) ## set the plot matrix

##qq norm plots
qqPlot(resid(mod0), xlab="Theoretical Quantiles", ylab = "Sample Quantiles", line = "quartiles", col.lines = "black", grid = FALSE)

##residual vs fitted plot
plot(fitted(mod0), resid(mod0), xlab = "Fitted Residuals", ylab = "Residuals") #residuals vs fitted
abline(h=0)

#******
##Linear mixed effects models for both years (Simpson Diversity) -- no documented effects overall

##to calculate Simpson Diversity
diversity_s <- as.data.frame(cbind(pitfall[1:4],diversity(index = "simpson", pitfall[5:24])))

#rename diversity variable
colnames(diversity_s)[5] ="Diversity"

##LMM to diversity treatment effects on diversity
d_simpson<-lmer(Diversity ~ Treatment * Year + (1|Block), data = diversity_s)
anova(d_simpson, F = "Kenward-Rogers")

#******
##analysis for each conducted separately for each taxonomic group

#===================================================================
##Year 2020 (initial conditions of the study)
#===================================================================

#subset 2020 data
year2020<-
  fun.groups%>%
  filter(Year == 2020)

##Linear mixed effects models for 2020 (taxonomic abundance) -- no documented effects
mod1<-lmer(log(Acari+1)~ Treatment + (1|Block), data = year2020)
anova(mod1, F = "Kenward-Rogers")

mod2<-lmer(Araneae ~ Treatment + (1|Block), data = year2020)
anova(mod2)

mod3<-lmer(log(Coleoptera)~ Treatment + (1|Block), data = year2020)
anova(mod3, F = "Kenward-Rogers")

mod4<-lmer(log(Collembola)~ Treatment + (1|Block), data = year2020)
anova(mod4, F = "Kenward-Rogers")

mod5<-lmer(log(Diplopoda+1)~ Treatment + (1|Block), data = year2020)
anova(mod5, F = "Kenward-Rogers")

mod6<-lmer(log(Formicidae)~ Treatment + (1|Block), data = year2020)
anova(mod6, F = "Kenward-Rogers")

mod7<-lmer(log(Isopoda+1)~ Treatment + (1|Block), data = year2020)
anova(mod7, F = "Kenward-Rogers")

mod8<-lmer(log(Total)~ Treatment + (1|Block), data = year2020)
anova(mod8, F = "Kenward-Rogers")



##To plot by each individual taxon
plt<-
  year2020%>%
  pivot_longer(Acari:Total, names_to = "Taxon", values_to = "Abundance")%>%
  mutate(Treatment= revalue(Treatment, c("Full connection" = "Full", 
                                         "Partial connection" = "Partial",
                                         "No connection" = "None")))

plt%>%
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

#******
##Linear mixed effects models for 2020 (Simpson Diversity) -- we found no effects on alpha diversity in 2020

##to calculate Simpson Diversity
diversity2020 <- as.data.frame(cbind(year2020[1:4],diversity(index = "simpson", year2020[5:11])))

#rename diversity variable
colnames(diversity2020)[5] ="Diversity"

##LMM to diversity treatment effects on diversity
mod.d20<-lmer(Diversity ~ Treatment + (1|Block), data = diversity2020)
anova(mod.d20, F = "Kenward-Rogers")

#===================================================================
##Year 2021 (the main analysis of the paper)
#===================================================================

#2021 analysis
year2021<-
  fun.groups%>%
  filter(Year == 2021)

#******
##Linear mixed effects models for 2021 to test for abundance
mod9<-lmer(log(Acari+1)~ Treatment + (1|Block), data = year2021)
anova(mod9, F = "Kenward-Rogers")

mod10<-lmer(Araneae ~ Treatment + (1|Block), data = year2021)
anova(mod10)
emmeans(mod10, list(pairwise ~ Treatment), adjust = "tukey")

mod11<-lmer(log(Coleoptera)~ Treatment + (1|Block), data = year2021)
anova(mod11, F = "Kenward-Rogers")
emmeans(mod11, list(pairwise ~ Treatment), adjust = "tukey")

mod12<-lmer(log(Collembola)~ Treatment + (1|Block), data = year2021)
anova(mod12, F = "Kenward-Rogers")

mod13<-lmer(log(Diplopoda+1)~ Treatment + (1|Block), data = year2021)
anova(mod13, F = "Kenward-Rogers")

mod14<-lmer(log(Formicidae)~ Treatment + (1|Block), data = year2021)
anova(mod14, F = "Kenward-Rogers")
emmeans(mod14, list(pairwise ~ Treatment), adjust = "tukey")

mod15<-lmer(log(Isopoda)~ Treatment + (1|Block), data = year2021)
anova(mod15, F = "Kenward-Rogers")
emmeans(mod15, list(pairwise ~ Treatment), adjust = "tukey")

mod16<-lmer(log(Total)~ Treatment + (1|Block), data = year2021)
anova(mod16, F = "Kenward-Rogers")
emmeans(mod16, list(pairwise ~ Treatment), adjust = "tukey")

##To plot by each individual taxon for 2021
plt<-
  year2021%>%
  pivot_longer(Acari:Total, names_to = "Taxon", values_to = "Abundance")%>%
  mutate(Treatment= revalue(Treatment, c("Full connection" = "Full", 
                                         "Partial connection" = "Partial",
                                         "No connection" = "None")))

plt%>%
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

#******
##Linear mixed effects models for 2021 (Simpson Diversity) -- no documented effects

##to calculate Simpson Diversity
diversity2021 <- as.data.frame(cbind(year2021[1:4],diversity(index = "simpson", year2021[5:11])))

#rename diversity variable
colnames(diversity2021)[5] ="Diversity"

##LMM to diversity treatment effects on diversity
mod.d21<-lmer(log(Diversity) ~ Treatment + (1|Block), data = diversity2021)
anova(mod.d21, F = "Kenward-Rogers")

#===================================================================
##Structural equation modeling

##need to detach lmerTest
detach("package:lmerTest", unload=TRUE)

#Need to import data the full dataset 

#relative pathname
pit <- file.path(".", "Data", "Pitfall_data_Yale_Myers_Final.csv")
print(pit)

#import data
p.sem.21<- read_csv(pit)%>%
  mutate(Block = as.factor(Block),
         Plot = as.factor(Plot),
         Year = as.factor(Year),
         Microdetritivores = Entomobryidae+Isotomidae+Sminthuridae+Tomoceridae+Acari,
         Predators = Lycosidae+Thomisidae+Linyphiidae+Salticidae+Gnaphosidae+Araneae.other+Staphylinidae+Carabidae+Chilopoda+Opiliones,
         Macrodetritivores = Diplopoda+Isopoda,
         Collembola = Entomobryidae+Isotomidae+Sminthuridae+Tomoceridae,
         Araneae = Lycosidae+Thomisidae+Linyphiidae+Salticidae+Gnaphosidae+Araneae.other,
         Coleoptera = Circulionidae+Staphylinidae+Carabidae+Coleoptera.other)%>%
  filter(Year == 2021)

#### PSEM chosen model##############
#Structural equation modeling

#transform data first
p.sem.21.transformed<-
  p.sem.21%>%
  mutate(Microdetritivores = log(Microdetritivores+1),
         Macrodetritivores = log(Macrodetritivores+1),
         Predators = log(Predators+1))

model1 <- psem(lmer(Microdetritivores ~ Connectivity + (1|Block/Plot), data = p.sem.21.transformed),
               lmer(Macrodetritivores ~ Connectivity + (1|Block/Plot), data = p.sem.21.transformed),
               lmer(Predators ~ Connectivity + Microdetritivores + Macrodetritivores + (1|Block/Plot), data = p.sem.21.transformed),
               Microdetritivores %~~% Macrodetritivores)

summary(model1, standardize = "scale")
fisherC(model1)
plot(model1)
rsquared(model1)


