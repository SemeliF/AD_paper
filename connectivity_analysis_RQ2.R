# Load libraries
library(lme4)
library(lmerTest)
library(lsmeans)


# Read in data
data = read.csv("dataSemeli.csv", sep = ";")
data$Subject = as.factor(paste("S", 1:40, sep = ""))
colnames(data)[1] = "PairA"

# Convert PLV values to numeric
plvcols = c("PairA", "PairB", "PairC", "PairD")
data[,plvcols] = lapply(data[,plvcols], function(x) { as.numeric(gsub(",", ".", x)) })


# Put data in long format
subData = function(data, pair = "A") {
  paircol = paste("Pair", pair, sep = "")
  sub = data[,c(paircol, "Group", "Band", "Subject")]
  sub$Pair = pair
  colnames(sub)[1] = "PLV"
  return(sub)
}

subList = list()
for(pair in c("A", "B", "C", "D")) {
  subList[[pair]] = subData(data, pair)
}

data = do.call(rbind, subList)

############## THETA BAND ANALYSIS #############
#  Run linear mixed effects regression model for THETA band
mod1_theta = lmer(PLV ~ (1 |Subject) + Pair*Group, data = data[data$Band == "Theta",])
mod2_theta = lmer(PLV ~ (1 |Subject) + Pair + Group, data = data[data$Band == "Theta",])
mod3_theta = lmer(PLV ~ (1 |Subject) + Pair , data = data[data$Band == "Theta",])
mod4_theta = lmer(PLV ~ (1 |Subject) + Group, data = data[data$Band == "Theta",])
mod5_theta = lmer(PLV ~ (1 |Subject), data = data[data$Band == "Theta",])

anova(mod1_theta,mod2_theta,mod3_theta,mod4_theta,mod5_theta)

# The best model is model 2, without the interaction between Pair and Group
# However, let's keep it in the model for now (see my comments in the e-mail) and use 
# model 1. Based on the model summary of model 1, we can state that we have a main effect 
# of Pair, a mean effect of Group, but not interaction between Group and Pair. To get an 
# idea about exactly which pairs differ between both groups, we can run post-hoc 
# comparisons

# Intermezzo: Check model assumptions

# Check assumption of normality of residuals
qqnorm(resid(mod1_theta), main = "")
shapiro.test(resid(mod1_theta))

# Check assumption of constant variance
plot(predict(mod1_theta), resid(mod1_theta), xlab = "predicted values", ylab = "residuals")
abline(h = 0, lty = 2)
lines(lowess(predict(mod1_theta), resid(mod1_theta)), col = "red")


# Look at contrasts of interest through post-hoc comparisons

# Generate posthoc comparisons
pc_theta = lsmeans(mod1_theta, "Group", by = "Pair", adjust = NULL)

# Generate p-values for pairwise contrasts of interest
cs_theta = contrast(pc_theta, method = "pairwise")

# Based on the post-hoc comparisons, we conclude that there is a significant difference 
# in connectivity between AD and HC in the Theta frequency band for Pair A and Pair D, 
# and a marginally significant difference for Pair B. For Pair C the effect is not 
# significant.



############## DELTA BAND ANALYSIS #############
#  Run linear mixed effects regression model for THETA band
mod1_delta = lmer(PLV ~ (1 |Subject) + Pair*Group, data = data[data$Band == "Delta",])
mod2_delta = lmer(PLV ~ (1 |Subject) + Pair + Group, data = data[data$Band == "Delta",])
mod3_delta = lmer(PLV ~ (1 |Subject) + Pair , data = data[data$Band == "Delta",])
mod4_delta = lmer(PLV ~ (1 |Subject) + Group, data = data[data$Band == "Delta",])
mod5_delta = lmer(PLV ~ (1 |Subject), data = data[data$Band == "Delta",])

anova(mod1_delta,mod2_delta,mod3_delta,mod4_delta,mod5_delta)

# The best model is model 2, without the interaction between Pair and Group
# However, let's keep it in the model for now (see my comments in the e-mail) and use 
# model 1. Based on the model summary of model 1, we can state that we have a main effect 
# of Pair, a mean effect of Group, but not interaction between Group and Pair. To get an 
# idea about exactly which pairs differ between both groups, we can run post-hoc 
# comparisons

# Intermezzo: Check model assumptions

# Check assumption of normality of residuals
qqnorm(resid(mod1_delta), main = "")
shapiro.test(resid(mod1_delta))

# Check assumption of constant variance
plot(predict(mod1_delta), resid(mod1_delta), xlab = "predicted values", ylab = "residuals")
abline(h = 0, lty = 2)
lines(lowess(predict(mod1_delta), resid(mod1_delta)), col = "red")


# Look at contrasts of interest through post-hoc comparisons

# Generate posthoc comparisons
pc_delta = lsmeans(mod1_delta, "Group", by = "Pair", adjust = NULL)

# Generate p-values for pairwise contrasts of interest
cs_delta = contrast(pc_delta, method = "pairwise")

# Based on the post-hoc comparisons, we conclude that there is a significant difference 
# in connectivity between AD and HC in the Theta frequency band for Pair A and Pair D, 
# and a marginally significant difference for Pair B. For Pair C the effect is not 
# significant.




