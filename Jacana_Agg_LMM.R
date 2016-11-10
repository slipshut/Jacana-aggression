#Comparing aggression between J.spinosa and J.jacana males and females

#set working directory
setwd("~/Dropbox/Jacanas/Aggression_Assay/Aggression_Data")

# Compare all species and sexes
agg <- read.csv("~/Dropbox/Jacanas/Aggression_Assay/Aggression_Data/aggression_14_nonest.csv", header=TRUE, sep=",")
head(agg)

# Exploratory Plots
boxplot(agg$Latency02 ~ agg$Species * agg$Sex.tested, ylab = "Latency to approach 0-2m (sec)")
boxplot(agg$Latency28 ~ agg$Species * agg$Sex.tested, ylab = "Latency to approach 2-8m (sec)")
boxplot(agg$Distance ~ agg$Species * agg$Sex.tested, ylab = "Distance (m)")
boxplot(agg$Swoops ~ agg$Species * agg$Sex.tested, ylab = "Swoops")
boxplot(agg$WingSpread ~ agg$Species * agg$Sex.tested, ylab = "WingSpread")
boxplot(agg$Threats ~ agg$Species * agg$Sex.tested, ylab = "Threats")
boxplot(agg$Pecks ~ agg$Species * agg$Sex.tested, ylab = "Pecks")
boxplot(agg$Flyovers ~ agg$Species * agg$Sex.tested, ylab = "Flyovers")
boxplot(agg$Vocalizations ~ agg$Species * agg$Sex.tested, ylab = "Vocalizations")

## Subset data
NJagg <- subset(agg, Species == "J. spinosa") 
WJagg <- subset(agg, Species == "J. jacana")

NJ.F.agg <- subset(NJagg, Sex.tested == "F") 
WJ.F.agg <- subset(WJagg, Sex.tested == "F") 
NJ.M.agg <- subset(NJagg, Sex.tested == "M") 
WJ.M.agg <- subset(WJagg, Sex.tested == "M")

# Get behavior summary stats for each sex and species
library(pastecs)
stat.desc(NJ.F.agg)
stat.desc(WJ.F.agg)
stat.desc(NJ.M.agg)
stat.desc(WJ.M.agg)

write.table(stat.desc(NJ.F.agg), file="/Users/saralipshutz/Dropbox/Jacanas/Aggression_Assay/Aggression_Data/NJ.F.agg.csv", sep=",")
write.table(stat.desc(WJ.F.agg), file="/Users/saralipshutz/Dropbox/Jacanas/Aggression_Assay/Aggression_Data/WJ.F.agg.csv", sep=",")
write.table(stat.desc(NJ.M.agg), file="/Users/saralipshutz/Dropbox/Jacanas/Aggression_Assay/Aggression_Data/NJ.M.agg.csv", sep=",")
write.table(stat.desc(WJ.M.agg), file="/Users/saralipshutz/Dropbox/Jacanas/Aggression_Assay/Aggression_Data/WJ.M.agg.csv", sep=",")

# Was one group (species + sex) less likely to respond (e.g. Max Distance = 14)? **SEE Greig et al. Evolution
library(MASS)
tbl = table(agg$Group, agg$Dist_Response)
tbl
chisq.test(tbl) 
# X-squared = 27.106, df = 3, p-value = 5.593e-06
# YES

# Was one species less likely to respond (e.g. Max Distance = 14)? 
tbl.species = table(agg$Species, agg$Dist_Response)
tbl.species
chisq.test(tbl.species) 
# X-squared = 24.412, df = 1, p-value = 7.778e-07
# YES - J. jacana was significantly less likely to respond

# Was one sex less likely to respond (e.g. Max Distance = 14)?
tbl.sex = table(agg$Sex.tested, agg$Dist_Response)
tbl.sex
chisq.test(tbl.sex) 
# X-squared = 0.13384, df = 1, p-value = 0.7145


# PCA
# log transform for assumptions of multi-normality
agg <- read.csv("~/Dropbox/Jacanas/Aggression_Assay/Aggression_Data/aggression_14_nonest.csv", header=TRUE, sep=",")
log.agg <- log(agg[, 15:24]+1) # transform log (x + 1)
log.agg
pcagg1 <- prcomp(log.agg, scale = TRUE, center = TRUE) 
summary(pcagg1) # PC 1,2, and 3 have eigenvalues higher than 1 and explain 61.4% of variation (cumulative proportion)
#  PC1    PC2    PC3     
# Standard deviation     1.838 1.2071 1.1399 
# Proportion of Variance 0.338 0.1457 0.1300 
# Cumulative Proportion  0.338 0.4837 0.6136
pcagg1$rotation # loadings
write.table(pcagg1$rotation, file="~/Dropbox/Jacanas/Aggression_Assay/Aggression_Data/aggpcaloadings.csv", sep=",")
pcagg1$x # scores
write.table(pcagg1$x, file="~/Dropbox/Jacanas/Aggression_Assay/Aggression_Data/aggpca.csv", sep=",")

agg <- read.csv("~/Dropbox/Jacanas/Aggression_Assay/Aggression_Data/aggression_14_nonest.csv", header=TRUE, sep=",")


# Visualize PCA

library(ggplot2)
library(plyr)
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Use species as a factor rather than numeric - don't think this is necessary though
agg <- read.csv("~/Dropbox/Jacanas/Aggression_Assay/Aggression_Data/aggression_14_nonest.csv", header=TRUE, sep=",")
aggPC1plot <- summarySE(agg, measurevar = "PC1", groupvars = c("Sex.tested", "Species"))
aggPC1plot$Species <- factor(aggPC1plot$Species)

# Error bars represent standard error of the mean
ggplot(aggPC1plot, aes(x=Species, y=PC1, fill=Sex.tested)) + 
  geom_bar(position=position_dodge(), stat="identity") + 
  geom_errorbar(aes(ymin=PC1-se, ymax=PC1+se), width=.2, position=position_dodge(.9)) +
  xlab("Species") +
  ylab("Aggression PC1\n") +
  scale_fill_brewer(palette = "Greys") + # Legend label, use darker colors
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(colour="black",size=20), axis.text.x = element_text(colour="black",size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y=element_text(vjust=1.0))

# ggplot2 boplot
library(ggplot2) 

# Exploratory: Plot PCs by group
boxplot(agg$PC1 ~ agg$Species * agg$Sex.tested, ylab = "PC1")
wilcox.test(agg$PC1 ~ agg$Species) # W = 720.5, p-value = 7.05e-06
wilcox.test(agg$PC1 ~ agg$Sex.tested) # W = 1172, p-value = 0.07983
boxplot(agg$PC2 ~ agg$Species * agg$Sex.tested, ylab = "PC2")
wilcox.test(agg$PC2 ~ agg$Species) # W = 1374.5, p-value = 0.644
wilcox.test(agg$PC2 ~ agg$Sex.tested) # W = 1196, p-value = 0.1087
boxplot(agg$PC3 ~ agg$Species * agg$Sex.tested, ylab = "PC3")
wilcox.test(agg$PC3 ~ agg$Species) # W = 1219.5, p-value = 0.1564
wilcox.test(agg$PC3 ~ agg$Sex.tested) # W = 1665, p-value = 0.2033
plot(agg$PC2 ~ agg$PC1)


#Box Cox transform - Distance
# library(MASS)
# boxcox(Distance ~ Species, data = agg)
# distspec1.boxcox=lm(Distance^1.8 ~ Species, data = agg)
# plot(distspec1.boxcox)
# boxplot(Distance^1.8 ~ Species, data = agg)
# boxplot(Distance ~ Species, data = agg)

# Explore relationship between response variable and random effects
boxplot(PC1 ~ Trial, data=agg)
boxplot(PC1 ~ Stimulus.Set, data=agg)
# boxplot(PC1 ~ Stimulus.Set, data=WJagg)
# boxplot(PC1 ~ Stimulus.Set, data=NJagg)
boxplot(PC1 ~ Mount.ID, data = agg)
# boxplot(PC1 ~ Mount.ID, data = WJagg)
# boxplot(PC1 ~ Mount.ID, data = NJagg)
boxplot(PC1 ~ Vocal.Stimulus, data = agg)
# boxplot(PC1 ~ Vocal.Stimulus, data = WJagg)
# boxplot(PC1 ~ Vocal.Stimulus, data = NJagg)
boxplot(PC1 ~ Site, data=agg)
# boxplot(PC1 ~ Site, data=NJagg)
# boxplot(PC1, data=WJagg)

#bw_LME tutorial 1
options(contrasts=c("contr.helmert","contr.poly")) ## Setting helmert contrasts prior to fitting models



# Use lmer to make mixed model fit by maximum likelihood 
#install.packages("lme4")
library (lme4)
agg <- read.csv("~/Dropbox/Jacanas/Aggression_Assay/Aggression_Data/aggression_14_nonest.csv", header=TRUE, sep=",")

entire.model.PC1 = lmer(PC1 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
summary(entire.model.PC1)
# anova(entire.model)
## Take out random effects (REML = FALSE) ### Check Zuur
## Take out site/trial
entire.model.PC1 = lmer(PC1 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
random.site.trial.null = lmer(PC1 ~ Species + Sex.tested + Species * Sex.tested + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
anova(entire.model.PC1, random.site.trial.null) ## Keep Trial as random effect

#install.packages("AICcmodavg")
library(AICcmodavg)
AICc(entire.model.PC1)
AICc(random.site.trial.null)
## Take out mount ID
entire.model.PC1 = lmer(PC1 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
random.mount.null = lmer(PC1 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
anova(entire.model.PC1, random.mount.null) ## Get rid of mount as random effect
AICc(random.mount.null)
## Take out vocal stimulus
entire.model.PC1 = lmer(PC1 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
random.vocal.null = lmer(PC1 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Mount.ID), data = agg, REML = FALSE)
anova(entire.model.PC1, random.vocal.null) ## Get rid of vocal as random effect
AICc(random.vocal.null)

final.model.PC1 = lmer(PC1 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial), data = agg, REML = FALSE)
anova(entire.model.PC1, final.model.PC1) # final model has lower BIC
AICc(final.model.PC1)

AIC(entire.model.PC1, random.site.trial.null, random.mount.null, random.vocal.null, final.model.PC1) 
# Final model has lowest AIC


## Take out fixed effects, start with interaction
final.model.PC1 = lmer(PC1 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial), data = agg, REML = FALSE)
interaction.null = lmer(PC1 ~ Species + Sex.tested + (1|Site/Trial), data = agg, REML = FALSE)
anova(final.model.PC1, interaction.null) # Keep all predictor variables
AICc(interaction.null)

## Now that we only have 1 random effect, we can use lme instead of lmer to get significance of predictors
## LME also allows us to use weights for unequal variance if need be
final.model.PC1 <- lme(PC1 ~ Species + Sex.tested + Species * Sex.tested, random = ~ 1|Site/Trial, data = agg, method = "REML")
summary(final.model.PC1) ## All fixed effects are significant predictors of response
#Fixed effects: PC1 ~ Species + Sex.tested + Species * Sex.tested 
#                                   Value Std.Error DF   t-value p-value
#(Intercept)                    -0.8639890 0.2938791 46 -2.939947  0.0051
#Species J. spinosa              1.0520393 0.4318401 16  2.436178  0.0269
#Sex.tested M                    0.3404522 0.2596615 46  1.311138  0.1963
#Species J. spinosa:Sex.testedM  1.0368683 0.3774208 46  2.747248  0.0086


## Check residuals for normality
final.model.PC1.res <- resid(final.model.PC1, type = "normalized")
hist(final.model.PC1.res)
qqnorm(residuals(final.model.PC1, type = "normalized"), ylab = "Residuals", xlab= "Normal Scores")

## Homogeneity of variance? (Fitted values vs. residuals) Yes
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
plot(final.model.PC1, which = c(1), col = 1, add.smooth = FALSE,caption = "")
final.model.PC1.res <- resid(final.model.PC1, type="normalized")
plot(agg$Sex.tested, final.model.PC1.res, xlab = "Sex", ylab = "Residuals")
## Residuals of males response has higher variance than females
plot(agg$Species, final.model.PC1.res, xlab = "Species", ylab = "Residuals")
## J. spinosa has higher variance than J. jacana
par(op) 
## Looks like residiual variation might increase with increasing fitted values

vfSex <- varIdent(form = ~1 | Sex.tested)
vfSpecies <- varIdent(form = ~1 | Species)

############### See if lmer can handle varIdent, otherwise, figure out how to use lme with multiple random effects
############### Use gls for model without any random effects
M.gls1 <- gls(PC1 ~ Species * Sex.tested, data = agg)
M.gls2 <- gls(PC1 ~ Species * Sex.tested, weights = vfSex, data = agg)
anova(M.gls1, M.gls2) ## Accounting for heteroscedasticity in Sex makes a better model

M.gls1 <- gls(PC1 ~ Species * Sex.tested, data = agg)
M.gls3 <- gls(PC1 ~ Species * Sex.tested, weights = vfSpecies, data = agg)
anova(M.gls1, M.gls3) ## Accounting for heteroscedasticity in Species also makes a better model

vf1 <- varComb(varIdent(form = ~1 | Sex.tested), varIdent(form = ~1 | Species))

homoscedastic.final.model.PC1 <- lme(PC1 ~ Species + Sex.tested + Species * Sex.tested, random = ~ 1|Site/Trial, data = agg, method = "REML", weights = vf1)
###final.model.PC1 = lmer(PC1 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial), data = agg, REML = FALSE)
final.model.PC1 <- lme(PC1 ~ Species + Sex.tested + Species * Sex.tested, random = ~ 1|Site/Trial, data = agg, method = "REML")
AICc(final.model.PC1)
anova(homoscedastic.final.model.PC1, final.model.PC1) ## Homoscedastic model is significantly better
AICc(homoscedastic.final.model.PC1)
summary(homoscedastic.final.model.PC1) ## All fixed effects are significant predictors of response
#Fixed effects: PC1 ~ Species + Sex.tested + Species * Sex.tested 
#                                     Value Std.Error DF   t-value p-value
#(Intercept)                    -0.8716885 0.2188513 46 -3.983017  0.0002
#Species J. spinosa              1.0139656 0.3229872 16  3.139336  0.0063
#Sex.tested M                    0.3451743 0.2262797 46  1.525432  0.1340
#SpeciesJ. spinosa:Sex.testedM   1.0891194 0.3813475 46  2.855976  0.0064

anova.lme(homoscedastic.final.model.PC1, type = "marginal", adjustSigma = F)
#numDF denDF   F-value p-value
#(Intercept)            1    46 15.864428  0.0002
#Species                1    16  9.855433  0.0063
#Sex.tested             1    46  2.326944  0.1340
#Species:Sex.tested     1    46  8.156600  0.0064

anova(homoscedastic.final.model.PC1) #### Which to use?
#numDF denDF   F-value p-value
#(Intercept)            1    46  5.373091  0.0250
#Species                1    16 10.182543  0.0057
#Sex.tested             1    46 16.003536  0.0002
#Species:Sex.tested     1    46  8.156600  0.0064

## Now do this for PC2!
entire.model.PC2 = lmer(PC2 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
summary(entire.model.PC2)
AICc(entire.model.PC2)
## anova(entire.model)
## Take out random effects (REML = FALSE) ### Check Zuur
## Take out site/trial
entire.model.PC2 = lmer(PC2 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
random.site.trial.null = lmer(PC2 ~ Species + Sex.tested + Species * Sex.tested + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
anova(entire.model.PC2, random.site.trial.null) ## Get rid of site/trial as random effect
AICc(random.site.trial.null)  
## Take out mount ID
entire.model.PC2 = lmer(PC2 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
random.mount.null = lmer(PC2 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
anova(entire.model.PC2, random.mount.null) ## Not significant - get rid of mount ID as random effect
AICc(random.mount.null)  
## Take out vocal stimulus
entire.model.PC2 = lmer(PC2 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
random.vocal.null = lmer(PC2 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Mount.ID), data = agg, REML = FALSE)
anova(entire.model.PC2, random.vocal.null) ## Get rid of vocal stimulus as random effect
AICc(random.vocal.null)

final.model.PC2 = lm(PC2 ~ Species + Sex.tested + Species * Sex.tested, data = agg)
anova(entire.model.PC2, final.model.PC2) # entire model has lower AIC, BIC, difference is significant
AICc(final.model.PC2)

AIC(entire.model.PC2, random.site.trial.null, random.mount.null, random.vocal.null, final.model.PC2) 
# final model with no random effects is better than entire model with no random effects

## Take out fixed effects, start with interaction
final.model.PC2 = gls(PC2 ~ Species + Sex.tested + Species * Sex.tested, data = agg, method = "ML")
interaction.null = gls(PC2 ~ Species + Sex.tested, data = agg, method = "ML")
anova(final.model.PC2, interaction.null) # Get rid of interaction
AICc(interaction.null)

final.model.PC2 = gls(PC2 ~ Species + Sex.tested, data = agg, method = "ML")
species.null = gls(PC2 ~ Sex.tested, data = agg, method = "ML")
anova(final.model.PC2, species.null)
AICc(species.null)

final.model.PC2 = gls(PC2 ~ Species + Sex.tested, data = agg, method = "ML")
sex.null = gls(PC2 ~ Species, data = agg, method = "ML")
anova(final.model.PC2, sex.null) 
AICc(sex.null)

entire.model.PC2 <- lmer(PC2 ~ Species + Sex.tested + Species * Sex.tested +  (1|Site/Trial) + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)

AIC(entire.model.PC2, random.site.trial.null, random.mount.null, random.vocal.null, final.model.PC2)
summary(final.model.PC2) ## All fixed effects are significant predictors of response
#Value Std.Error    t-value p-value
#(Intercept)       -0.2109029 0.1982651 -1.0637421  0.2899
#SpeciesJ. spinosa  0.0375990 0.2322127  0.1619161  0.8717
#Sex.testedM        0.3799558 0.2316145  1.6404661  0.1039

## Check residuals for normality
final.model.res.PC2 <- resid(final.model.PC2, type = "normalized")
hist(final.model.res.PC2)
qqnorm(residuals(final.model.PC2, type = "normalized"), ylab = "Residuals", xlab= "Normal Scores")

## Homogeneity of variance? (Fitted values vs. residuals) Yes
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
plot(final.model.PC2, which = c(1), col = 1, add.smooth = FALSE,caption = "")
final.model.PC2.res <- resid(final.model.PC2, type="normalized")
plot(agg$Sex.tested, final.model.PC2.res, xlab = "Sex", ylab = "Residuals")
## Residuals of males response has higher variance than females
plot(agg$Species, final.model.PC2.res, xlab = "Species", ylab = "Residuals")
## J. spinosa has higher variance than J. jacana
par(op) 
## Looks like residiual variation might increase with increasing fitted values

vfSex <- varIdent(form = ~1 | Sex.tested)
vfSpecies <- varIdent(form = ~1 | Species)

############### See if lmer can handle varIdent, otherwise, figure out how to use lme with multiple random effects
M.gls1 <- gls(PC2 ~ Species * Sex.tested, data = agg)
M.gls2 <- gls(PC2 ~ Species * Sex.tested, weights = vfSex, data = agg)
anova(M.gls1, M.gls2) ## Accounting for heteroscedasticity in Sex makes a better model

M.gls1 <- gls(PC2 ~ Species * Sex.tested, data = agg)
M.gls3 <- gls(PC2 ~ Species * Sex.tested, weights = vfSpecies, data = agg)
anova(M.gls1, M.gls3) ## Accounting for heteroscedasticity in Species also makes a better model

vf1 <- varComb(varIdent(form = ~1 | Sex.tested), varIdent(form = ~1 | Species))

homoscedastic.final.model.PC2 <- gls(PC2 ~ Species + Sex.tested, data = agg, method = "REML", weights =vf1)
####Error in model.frame.default(data = agg, weights = varIdent(form = ~1 |  : variable lengths differ (found for '(weights)')
final.model.PC2 <- gls(PC2 ~ Species + Sex.tested, data = agg, method = "REML")

anova(homoscedastic.final.model.PC2, final.model.PC2) ## Homoscedastic model is better
summary(homoscedastic.final.model.PC2) ## No fixed effects are significant predictors of response

AICc(homoscedastic.final.model.PC2)
anova.lme(homoscedastic.final.model.PC2, type = "marginal", adjustSigma = F) 
anova(homoscedastic.final.model.PC2)
# Sex is a trending predictor of PC2
#Denom. DF: 105 
#           numDF   F-value p-value
#(Intercept)     1 2.1544374  0.1451
#Species         1 0.0043297  0.9477
#Sex.tested      1 2.9413491  0.0893


## Now do this for PC3!

entire.model.PC3 = lmer(PC3 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
summary(entire.model.PC3)
AICc(entire.model.PC3)
# anova(entire.model)
## Take out random effects (REML = FALSE) ### Check Zuur
## Take out site/trial
entire.model.PC3 = lmer(PC3 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
random.site.trial.null = lmer(PC3 ~ Species + Sex.tested + Species * Sex.tested + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
AICc(random.site.trial.null)
anova(entire.model.PC3, random.site.trial.null) ## Get rid of site/trial as random effect
## Take out mount ID
entire.model.PC3 = lmer(PC3 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
random.mount.null = lmer(PC3 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
AICc(random.mount.null)
anova(entire.model.PC3, random.mount.null) ## Not significant - get rid of mount ID as random effect
## Take out vocal stimulus
entire.model.PC3 = lmer(PC3 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
random.vocal.null = lmer(PC3 ~ Species + Sex.tested + Species * Sex.tested + (1|Site/Trial) + (1|Mount.ID), data = agg, REML = FALSE)
anova(entire.model.PC3, random.vocal.null) ## Get rid of vocal stimulus as random effect

final.model.PC3 = lm(PC3 ~ Species + Sex.tested + Species * Sex.tested, data = agg)
AICc(final.model.PC3)
anova(entire.model.PC3, final.model.PC3) # entire model has lower AIC, BIC, difference is significant

AIC(entire.model.PC3, random.site.trial.null, random.mount.null, random.vocal.null, final.model.PC3) 
# final model with no random effects is better than entire model with all random effects
summary(final.model.PC3)

## Take out fixed effects, start with interaction
final.model.PC3 = lm(PC3 ~ Species + Sex.tested + Species * Sex.tested, data = agg)
interaction.null = lm(PC3 ~ Species + Sex.tested, data = agg)
AICc(interaction.null)
anova(final.model.PC3, interaction.null) # Get rid of interaction

final.model.PC3 = lm(PC3 ~ Species + Sex.tested, data = agg)
AICc(final.model.PC3)
interaction.null = lm(PC3 ~  Species, data = agg)
anova(final.model.PC3, interaction.null) # Get rid of interaction

entire.model.PC3 <- lmer(PC3 ~ Species + Sex.tested + Species * Sex.tested +  (1|Site/Trial) + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE)
AICc(entire.model.PC3)
summary(entire.model.PC3) ## All fixed effects are significant predictors of response


## Check residuals for normality
final.model.res.PC3 <- resid(entire.model.PC2, type = "deviance")
hist(final.model.res.PC3)
qqnorm(residuals(entire.model.PC3, type = "deviance"), ylab = "Residuals", xlab= "Normal Scores")

## Homogeneity of variance? (Fitted values vs. residuals) Yes
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
plot(entire.model.PC3, which = c(1), col = 1, add.smooth = FALSE,caption = "")
entire.model.PC3.res <- resid(entire.model.PC3, type="deviance")
plot(agg$Sex.tested, entire.model.PC3.res, xlab = "Sex", ylab = "Residuals")
## Residuals of males response has higher variance than females
plot(agg$Species, entire.model.PC3.res, xlab = "Species", ylab = "Residuals")
## J. spinosa has higher variance than J. jacana
par(op) 
## Looks like residual variation might increase with increasing fitted values

vfSex <- varIdent(form = ~1 | Sex.tested)
vfSpecies <- varIdent(form = ~1 | Species)

############### See if lmer can handle varIdent, otherwise, figure out how to use lme with multiple random effects
M.gls1 <- gls(PC3 ~ Species * Sex.tested, data = agg)
M.gls2 <- gls(PC3 ~ Species * Sex.tested, weights = vfSex, data = agg)
anova(M.gls1, M.gls2) ## Accounting for heteroscedasticity in Sex makes a better model

M.gls1 <- gls(PC3 ~ Species * Sex.tested, data = agg)
M.gls3 <- gls(PC3 ~ Species * Sex.tested, weights = vfSpecies, data = agg)
anova(M.gls1, M.gls3) ## Accounting for heteroscedasticity in Species does NOT make a better model

vf1 <- varComb(varIdent(form = ~1 | Sex.tested), varIdent(form = ~1 | Species))

homoscedastic.final.model.PC3 <- lmer(PC3 ~ Species + Sex.tested + Species * Sex.tested +  (1|Site/Trial) + (1|Mount.ID) + (1|Vocal.Stimulus), data = agg, REML = FALSE, weights =vfSex)
####Error in model.frame.default(data = agg, weights = varIdent(form = ~1 |  : variable lengths differ (found for '(weights)')

homoscedastic.final.model.PC3 <- lme(PC3 ~ Species + Sex.tested + Species * Sex.tested, random = ~ 1|Site/Trial, data = agg, method = "REML", weights =vfSex)
AICc(homoscedastic.final.model.PC3)
final.model.PC3 <- lme(PC3 ~ Species + Sex.tested + Species * Sex.tested, random = ~ 1|Site/Trial, data = agg, method = "REML")
AICc(final.model.PC3)

anova(homoscedastic.final.model.PC3, final.model.PC3) ## Homoscedastic model is better
summary(homoscedastic.final.model.PC3) ## All fixed effects are significant predictors of response
anova.lme(homoscedastic.final.model.PC3, type = "marginal", adjustSigma = F) 
#numDF denDF   F-value p-value
#(Intercept)            1    46 0.3256052  0.5710
#Species                1    16 3.1397455  0.0954
#Sex.tested             1    46 0.0021260  0.9634
#Species:Sex.tested     1    46 1.5864847  0.2142

anova(homoscedastic.final.model.PC3)
#numDF denDF   F-value p-value
#(Intercept)            1    46 0.1360777  0.7139
#Species                1    16 1.7844578  0.2003
#Sex.tested             1    46 1.2252241  0.2741
#Species:Sex.tested     1    46 1.5864847  0.2142


