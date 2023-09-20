#03/09/2023
#Bushra Moussaoui
#FINAL Reanalsysis with full data set (with recount of indiv 431YMHM) (Testing for correlation between vocal measures and FoxP2 MMST/VSP Protein Expression and Expression Difference with Age)

library(warbleR)

#We are extracting vocal measures during recording block 5 since this is closest temporally to time of brain extraction. 
#This also happened to be the time block that appears to have strongest changes in vocal measures.

#Methods for prepping dataframe:
##Download dataset from lab gdrive and keep only bird avg MMST/VSP rows, name as "CellCounts_03_09_2023_Complete_Final.xlsx"
##Start with in "CellCounts_03_09_2023_Complete_Final.xlsx" C:\Users\bushr\OneDrive - New Mexico State University\NMSU - MS\Wright Lab\Thesis\AgeLearn\Neural analysis\R Data Analysis
##Make a copy, named: "Vocal_Brain_Corr_03_09_2023_Final.csv" and saved in same file path as above.

##Open "agg_dat.csv" (in:C:\Users\bushr\OneDrive - New Mexico State University\NMSU - MS\Wright Lab\Thesis\AgeLearn\Call Analysis\Final_selections) and copy rows from birds for recording block 5
###Look up bird ID in "AgeLearn_Treatments_Master.csv"spreadsheet to match the long form ID to the short one used in vocal data (example: 268YGHM = 10C)

library(ggplot2)
library(tidyverse)
library(cowplot)
library(car)
library(pwr)

#Load data
setwd("C:/Users/bushr/OneDrive - New Mexico State University/NMSU - MS/Wright Lab/Thesis/AgeLearn/Neural analysis/R Data Analysis")
voc_brain_corr <- read.csv("Vocal_Brain_Corr_03_09_2023_Final.csv", header = TRUE)
head(voc_brain_corr)
view(voc_brain_corr)

#NAs summary
##12C (Old), 2C (Old), 6C (Old), 4A (Old) had fewer than 6 calls in recording block 5. None of the new 6 birds had fewer than 6 calls.


#We can take 1 - acoustic.overlap, acoustic.overlap being calculated as the intersect over union of an individual’s beginning and ending acoustic spaces, to be able to represent vocal plasticity such that higher values representing a more plastic vocal repertoire.
    #vocal.plasticity = (1-acoustic.overlap)

voc_brain_corr$vocal.plasticity <- (1 - (voc_brain_corr$acoustic.overlap))


#Statistics approach: ANCOVA
#https://www.tutorialspoint.com/r/r_analysis_of_covariance.htm
#https://dzchilds.github.io/stats-for-bio/two-way-ancova-in-r.html

#response variable = vocal learning measure (continuous)


#2 explanatory variables
## FoxP2 MMST/VSP expression (continuous)
## Age class (categorical: Young, Old)

voc_brain_corr$treatment <- factor(voc_brain_corr$treatment, levels = c("Young", "Old"))

#Visualize data
theme_set(theme_classic())
area.change <- ggplot(voc_brain_corr, aes(x = mmst.vsp, 
y = acoustic.area, colour = treatment)) + 
geom_point(size = 6) + geom_smooth(method = "lm", se = FALSE, size = 1.5) + labs(y = "Vocal diversity", x = "FoxP2 MMSt/VSP expression", color = "Adult age") +  theme(text = element_text(size=25), axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), axis.text = element_text(colour = "black")) + scale_colour_manual(values = c("#9f79db", "#b7d851")) + coord_cartesian(xlim = c(0.25, 1.25))
area.change

ggsave("FOXP2_VocalDiversity_Correlation_Final_03_09_2023.png", width = 7, height = 6, units = "in")

overlap.self <- ggplot(voc_brain_corr, aes(x = mmst.vsp, y = vocal.plasticity, colour = treatment)) + 
geom_point(size = 6) + geom_smooth(method = "lm", se = FALSE, size = 1.5) + labs(y = "Vocal plasticity", x = "FoxP2 MMSt/VSP expression", color = "Adult age") +  theme(text = element_text(size=25), axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), axis.text = element_text(colour = "black")) + scale_colour_manual(values = c("#9f79db", "#b7d851"), name = "Adult age") + coord_cartesian(xlim = c(0.25, 1.25), ylim = c(0,1))
overlap.self

ggsave("FOXP2_VocalPlasticity_Correlation_Final__03_09_2023.png", width = 7, height = 6, units = "in")

overlap.group <- ggplot(voc_brain_corr, aes(x = mmst.vsp, 
y = acoustic.overlap.to.group, 
colour = treatment)) + 
geom_point(size = 6) + geom_smooth(method = "lm", se = FALSE, size = 1.5) + labs(y = "Vocal convergence", x = "FoxP2 MMSt/VSP expression", color = "Adult age") +  theme(text = element_text(size=25), axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), axis.text = element_text(colour = "black")) + scale_colour_manual(values = c("#9f79db", "#b7d851"), name = "Adult age")+ coord_cartesian(xlim = c(0.25, 1.25), ylim = c(0,1))
overlap.group

ggsave("FOXP2_VocalConvergence_Correlation_Final__03_03_2023.png", width = 7, height = 6, units = "in")


#models without age as a predictor (foxp2 expression does not predict any vocal measure)
model.01 <- lm(acoustic.area ~ mmst.vsp, data = voc_brain_corr)
anova(model.01)

model.02 <- lm(vocal.plasticity ~ mmst.vsp, data = voc_brain_corr)
anova(model.02)

model.03 <- lm(acoustic.overlap.to.group ~ mmst.vsp, data = voc_brain_corr)
anova(model.03)


#models with foxp2 expression & age as predictors
model.1 <- lm(acoustic.area ~ mmst.vsp * treatment, data = voc_brain_corr)
plot(model.1, add.smooth = FALSE, which = 1) #linearity?
plot(model.1, which = 2) #normality? 
shapiro.test(voc_brain_corr$acoustic.area) #Passes normality
shapiro.test(voc_brain_corr$mmst.vsp) #Passes normality
plot(model.1, add.smooth = FALSE, which = 3) #constant variance?

anova(model.1)

model.2 <- lm(vocal.plasticity ~ mmst.vsp * treatment, data = voc_brain_corr)
plot(model.2, add.smooth = FALSE, which = 1) #linearity
plot(model.2, which = 2) #normality
shapiro.test(voc_brain_corr$vocal.plasticity) #Passes normality
plot(model.2, add.smooth = FALSE, which = 3) #constant variance?

anova(model.2)

model.3 <- lm(acoustic.overlap.to.group ~ mmst.vsp * treatment, data = voc_brain_corr)
summary(model.3)
plot(model.3, add.smooth = FALSE, which = 1) #linearity
plot(model.3, which = 2) #normality
shapiro.test(voc_brain_corr$acoustic.overlap.to.group)#Passes normality
plot(model.3, add.smooth = FALSE, which = 3) #constant variance?

anova(model.3)



##################################################
#Simple difference in FoxP2 expression by brain region and adult age class

voc_brain_corr_mmst_vsp <- read.csv("Vocal_Brain_Corr_03_09_2023_Final_WithAdditionalMMSTVSP.csv", header = TRUE)
view(voc_brain_corr_mmst_vsp)
head(voc_brain_corr_mmst_vsp)
summary(voc_brain_corr_mmst_vsp)

voc_brain_corr_mmst_vsp$treatment <- factor(voc_brain_corr_mmst_vsp$treatment, levels = c("Young", "Old"))


#Boxplot: % expression FoxP2 in MMST

?geom_boxplot()
ggplot(voc_brain_corr_mmst_vsp, aes(x = treatment, y = mmst, fill = treatment)) + labs(y="FoxP2 MMSt expression", x = "Adult age") + geom_boxplot(color = "black", lwd = 1.5) + theme_classic() + theme(text = element_text(size=35), legend.position = "none", axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), axis.text = element_text(colour = "black")) + scale_fill_manual(values = c("#9f79db", "#b7d851")) + scale_x_discrete(limits = c("Young", "Old")) + geom_hline(yintercept = 1, linetype="dashed", size = 1) + coord_cartesian(ylim = c(0.2,1.25))
ggsave("FOXP2_MMSt_Age_Boxplot_03_09_2023.png", width = 6, height = 8, units = "in")

ggplot(voc_brain_corr_mmst_vsp, aes(x = treatment, y = mmst, fill = treatment)) + labs(y="FoxP2 MMSt expression", x = "Adult age") + geom_boxplot(color = "black", lwd = 1.5) + theme_classic() + theme(text = element_text(size=35), legend.position = "none", axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), axis.text = element_text(colour = "black")) + scale_fill_manual(values = c("#9f79db", "#b7d851")) + scale_x_discrete(limits = c("Young", "Old")) + geom_hline(yintercept = 1, linetype="dashed", size = 1) + coord_cartesian(ylim = c(0.2,1.25))+ geom_jitter(shape=21, fill = "#CCCCCC", stroke = 1.5, size = 2, position=position_jitter(0.2))
ggsave("FOXP2_MMSt_Age_BoxplotJitter_03_09_2023.png", width = 6, height = 8, units = "in")

#Boxplot: % expression FoxP2 in VSP
ggplot(voc_brain_corr_mmst_vsp, aes(x = treatment, y = vsp, fill = treatment)) + labs(y="FoxP2 VSP expression", x = "Adult age") + geom_boxplot(color = "black", lwd = 1.5) + theme_classic() + theme(text = element_text(size=35), legend.position = "none", axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), axis.text = element_text(colour = "black")) + scale_fill_manual(values = c("#9f79db", "#b7d851")) + scale_x_discrete(limits = c("Young", "Old")) + geom_hline(yintercept = 1, linetype="dashed", size = 1) + coord_cartesian(ylim = c(0.2,1.25))
ggsave("FOXP2_VSP_Age_Boxplot_03_09_2023.png", width = 6, height = 8, units = "in")

ggplot(voc_brain_corr_mmst_vsp, aes(x = treatment, y = vsp, fill = treatment)) + labs(y="FoxP2 VSP expression", x = "Adult age") + geom_boxplot(color = "black", lwd = 1.5) + theme_classic() + theme(text = element_text(size=35), legend.position = "none", axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), axis.text = element_text(colour = "black")) + scale_fill_manual(values = c("#9f79db", "#b7d851")) + scale_x_discrete(limits = c("Young", "Old")) + geom_hline(yintercept = 1, linetype="dashed", size = 1) + coord_cartesian(ylim = c(0.2,1.25))+ geom_jitter(shape=21, fill = "#CCCCCC", stroke = 1.5, size = 2, position=position_jitter(0.2))
ggsave("FOXP2_VSP_Age_BoxplotJitter_03_09_2023.png", width = 6, height = 8, units = "in")

#Boxplot: MMST/VSP
ggplot(voc_brain_corr_mmst_vsp, aes(x = treatment, y = mmst.vsp, fill = treatment)) + labs(y="FoxP2 MMSt/VSP expression", x = "Adult age") + geom_boxplot(color = "black", lwd = 1.5) + theme_classic() + theme(text = element_text(size=35), legend.position = "none", axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), axis.text = element_text(colour = "black")) + scale_fill_manual(values = c("#9f79db", "#b7d851")) + scale_x_discrete(limits = c("Young", "Old")) + geom_hline(yintercept = 1, linetype="dashed", size = 1) + coord_cartesian(ylim = c(0.2,1.25))
ggsave("FOXP2_MMSt.VSP_Age_Boxplot_03_09_2023.png", width = 6, height = 8, units = "in")

ggplot(voc_brain_corr_mmst_vsp, aes(x = treatment, y = mmst.vsp, fill = treatment)) + labs(y="FoxP2 MMSt/VSP expression", x = "Adult age") + geom_boxplot(color = "black", lwd = 1.5) + theme_classic() + theme(text = element_text(size=35), legend.position = "none", axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), axis.text = element_text(colour = "black")) + scale_fill_manual(values = c("#9f79db", "#b7d851")) + scale_x_discrete(limits = c("Young", "Old")) + geom_hline(yintercept = 1, linetype="dashed", size = 1) + coord_cartesian(ylim = c(0.2,1.25))+ geom_jitter(shape=21, fill = "#CCCCCC", stroke = 1.5, size = 2, position=position_jitter(0.2))
ggsave("FOXP2_MMSt.VSP_Age_BoxplotJitter_03_09_2023.png", width = 6, height = 8, units = "in")


hist(voc_brain_corr_mmst_vsp$mmst)
shapiro.test(voc_brain_corr_mmst_vsp$mmst) #passes normality
qqPlot(voc_brain_corr_mmst_vsp$mmst) #slightly deviates from normality
leveneTest(mmst ~ treatment, data = voc_brain_corr_mmst_vsp) #equal variance

hist(voc_brain_corr_mmst_vsp$vsp)
shapiro.test(voc_brain_corr_mmst_vsp$vsp) #passes normality
qqPlot(voc_brain_corr_mmst_vsp$vsp) #slightly deviates from normality
leveneTest(vsp ~ treatment, data = voc_brain_corr_mmst_vsp) #equal variance

hist(voc_brain_corr_mmst_vsp$mmst.vsp)
shapiro.test(voc_brain_corr_mmst_vsp$mmst.vsp) #passes normality
qqPlot(voc_brain_corr_mmst_vsp$mmst.vsp) #slightly deviates from normality
leveneTest(mmst.vsp ~ treatment, data = voc_brain_corr_mmst_vsp) #equal variance


t.test(mmst ~ treatment, var.equal=TRUE, data = voc_brain_corr_mmst_vsp)
#       Two Sample t-test

#data:  mmst by treatment
#t = -2.1217, df = 22, p-value = 0.04536
#alternative hypothesis: true difference in means between group Young and group Old is not equal to 0
#95 percent confidence interval:
# -0.122325596 -0.001395428
#sample estimates:
#mean in group Young   mean in group Old
#          0.3090468           0.3709074


t.test(vsp ~ treatment, var.equal=TRUE, data = voc_brain_corr_mmst_vsp)
#Two Sample t-test

#data:  vsp by treatment
#t = -0.94644, df = 22, p-value = 0.3542
#alternative hypothesis: true difference in means between group Young and group Old is not equal to 0
#95 percent confidence interval:
# -0.12424707  0.04637955
#sample estimates:
#mean in group Young   mean in group Old
#          0.4411517           0.4800855


t.test(mmst.vsp ~ treatment, var.equal=TRUE, data = voc_brain_corr_mmst_vsp)
# Two Sample t-test

#data:  mmst.vsp by treatment
#t = -0.85662, df = 22, p-value = 0.4009
#alternative hypothesis: true difference in means between group Young and group Old is not equal to 0
#95 percent confidence interval:
# -0.22828768  0.09482511
#sample estimates:
#mean in group Young   mean in group Old
#          0.7438620           0.8105933


#Also tried non-parametric tests...same results as t-test
wilcox.test(mmst ~ treatment, data = voc_brain_corr_mmst_vsp, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(vsp ~ treatment, data = voc_brain_corr_mmst_vsp, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(mmst.vsp ~ treatment, data = voc_brain_corr_mmst_vsp, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)


###### summary stats
install.packages("plotrix")
library(plotrix)

voc_brain_corr_young <- voc_brain_corr_mmst_vsp %>% filter(treatment == "Young")
voc_brain_corr_old <- voc_brain_corr_mmst_vsp %>% filter(treatment == "Old")

mean(voc_brain_corr_young$mmst)
std.error(voc_brain_corr_young$mmst)
mean(voc_brain_corr_old$mmst)
std.error(voc_brain_corr_old$mmst)

mean(voc_brain_corr_young$vsp)
std.error(voc_brain_corr_young$vsp)
mean(voc_brain_corr_old$vsp)
std.error(voc_brain_corr_old$vsp)

mean(voc_brain_corr_young$mmst.vsp)
std.error(voc_brain_corr_young$mmst.vsp)
mean(voc_brain_corr_old$mmst.vsp)
std.error(voc_brain_corr_old$mmst.vsp)


