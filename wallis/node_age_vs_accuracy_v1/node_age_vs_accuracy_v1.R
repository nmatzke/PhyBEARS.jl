

# START Plot to PDF, if desired
pdffn = "node_age_vs_accuracy_plots_v1.pdf"
pdf(file=pdffn, width=9, height=9)




wd = "/GitHub/PhyBEARS.jl/wallis/node_age_vs_accuracy_v1/"
setwd(wd)

fn1 = "models_table2_7000.csv"
fn2 = "plot_table2_7000.csv"

df1 = read.csv(fn1)
df2 = read.csv(fn2)

head(df1)
head(df2)

hist(df1$control_csp)

hist(df1$spread_csp)

hist(df2$CSP)




#######################################################
# Just for kicks, remove the 1s and 0s
# Also, try bin plots
#######################################################

# 1. Control_CSP

TF1 = df1$control_csp < 0.99999999
TF2 = df1$control_csp > 0.00000001

control_csp_fractions = df1$control_csp[(TF1 + TF2)==2]
control_csp_fractions_node_ages = df1$node_age[(TF1 + TF2)==2]

lm_df1_control_CSP = linear_regression_plot(x=control_csp_fractions_node_ages, y=control_csp_fractions)
title("df1$control_csp, no 0.0s and 1.0s")

# Binned plots
# install.packages("binsreg")
library(binsreg)


# Bins plot *including* 0s and 1s ("plotyrange" doesn't appear to work to set y-axis, boo!)
binsreg(y=df1$control_csp, x=df1$node_age, line=NULL, noplot=FALSE, plotyrange=c(0,1))
title("df1$control_csp")
# with 95% CI
binsreg_df1_control_csp = binsreg(y=df1$control_csp, x=df1$node_age, line=TRUE, ci=TRUE, noplot=FALSE, plotyrange=c(0,1))
title("df1$control_csp")

# Bins plot *excluding* 0s and 1s ("plotyrange" doesn't appear to work to set y-axis, boo!)
binsreg(y=control_csp_fractions, x=control_csp_fractions_node_ages, line=NULL, plotyrange=c(0,1))
title("df1$control_csp, no 0.0s and 1.0s")

# with 95% CI
binsreg(y=control_csp_fractions, x=control_csp_fractions_node_ages, line=TRUE, ci=TRUE, plotyrange=c(0,1))
title("df1$control_csp, no 0.0s and 1.0s")







# 2. spread_csp

TF1 = df1$spread_csp < 0.99999999
TF2 = df1$spread_csp > 0.00000001

spread_csp_fractions = df1$spread_csp[(TF1 + TF2)==2]
spread_csp_fractions_node_ages = df1$node_age[(TF1 + TF2)==2]

lm_df1_spread_csp = linear_regression_plot(x=spread_csp_fractions_node_ages, y=spread_csp_fractions)
title("df1$spread_csp, no 0.0s and 1.0s")


# Binned plots
# install.packages("binsreg")
library(binsreg)


# Bins plot *including* 0s and 1s ("plotyrange" doesn't appear to work to set y-axis, boo!)
binsreg(y=df1$spread_csp, x=df1$node_age, line=NULL, noplot=FALSE, plotyrange=c(0,1))
title("df1$spread_csp")
# with 95% CI
binsreg_df1_spread_csp = binsreg(y=df1$spread_csp, x=df1$node_age, line=TRUE, ci=TRUE, noplot=FALSE, plotyrange=c(0,1))
title("df1$spread_csp")


# Bins plot *excluding* 0s and 1s ("plotyrange" doesn't appear to work to set y-axis, boo!)
binsreg(y=spread_csp_fractions, x=spread_csp_fractions_node_ages, line=NULL, plotyrange=c(0,1))
title("df1$spread_csp, no 0.0s and 1.0s")

# with 95% CI
binsreg(y=spread_csp_fractions, x=spread_csp_fractions_node_ages, line=TRUE, ci=TRUE, plotyrange=c(0,1))
title("df1$spread_csp, no 0.0s and 1.0s")



# 3. df2$CSP


TF1 = df2$CSP < 0.99999999
TF2 = df2$CSP > 0.00000001

CSP_fractions = df2$CSP[(TF1 + TF2)==2]
CSP_fractions_node_ages = df2$node_age[(TF1 + TF2)==2]

lm_df2_CSP = linear_regression_plot(x=CSP_fractions_node_ages, y=CSP_fractions)
title("df2$CSP, no 0.0s and 1.0s")


# Binned plots
# install.packages("binsreg")
library(binsreg)


# Bins plot *including* 0s and 1s ("plotyrange" doesn't appear to work to set y-axis, boo!)
binsreg(y=df2$CSP, x=df2$node_age, line=NULL, noplot=FALSE, plotyrange=c(0,1))
title("df2$CSP")
# with 95% CI
binsreg_df2_CSP = binsreg(y=df2$CSP, x=df2$node_age, line=TRUE, ci=TRUE, noplot=FALSE, plotyrange=c(0,1))
title("df2$CSP")


# Bins plot *excluding* 0s and 1s ("plotyrange" doesn't appear to work to set y-axis, boo!)
binsreg(y=CSP_fractions, x=CSP_fractions_node_ages, line=NULL, plotyrange=c(0,1))
title("df2$CSP, no 0.0s and 1.0s")
# with 95% CI
binsreg(y=CSP_fractions, x=CSP_fractions_node_ages, line=TRUE, ci=TRUE, plotyrange=c(0,1))
title("df2$CSP, no 0.0s and 1.0s")


# END Plot to PDF, if desired

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


summary(lm_df1_control_csp)
summary(lm_df1_spread_csp)
summary(lm_df2_CSP)


summary(binsreg_df1_control_csp)
summary(binsreg_df1_spread_csp)
summary(binsreg_df2_CSP)









#######################################################
# Logistic regression with a Generalized Linear Model (GLM), logit link function
# This means that the log-odds of a "1" are being modeled as a linear function of x=node_age
#######################################################


# 1. df1$control_csp
glm_df1_control_csp = glm(formula=control_csp~node_age, data=df1, family=binomial(link="logit"))
summary(glm_df1_control_csp)

# Coefficients:
# (Intercept)     node_age  
#      2.2318      -0.2428  

# Interpretation: GLM logistic regression means that the 1/0 y values are being predicted by logit(intercept + slope*x)

# install.packages("blm")  # for logit() and expit/inverse of logit
library(blm)

# So, when node_age = 0.0, the probability of a "1" is...
blm::expit(( 2.2318 + 0.0*-0.2428))
# 0.903069

# So, when node_age = 1.0, the probability of a "1" is...
blm::expit(( 2.2318 + 1.0*-0.2428))
# 0.8796373

# ...in other words, the log-odds of a "1" decrease by 0.2428 per unit of node age



#######################################################
# How log-odds, logit and expit work
#######################################################
prob = 0.1
blm::logit(prob)
prob = 0.2
blm::logit(prob)
prob = 0.5
blm::logit(prob)
prob = 0.6
blm::logit(prob)
prob = 0.9
blm::logit(prob)

# expit is the inverse of logit
prob = 0.9
blm::logit(prob)
blm::expit(logit(prob))


# Summaries
summary(glm_df1_control_csp)
logLik(glm_df1_control_csp)
AIC(glm_df1_control_csp)

# calculate McFadden's R-squared for GLM model with logit link function
# https://www.statology.org/glm-r-squared/
with(summary(glm_df1_control_csp), 1 - deviance/null.deviance)




# 2. df1$spread_csp
glm_df1_spread_csp = glm(formula=spread_csp~node_age, data=df1, family=binomial(link="logit"))
summary(glm_df1_spread_csp)
logLik(glm_df1_spread_csp)
AIC(glm_df1_spread_csp)

# calculate McFadden's R-squared for GLM model with logit link function
# https://www.statology.org/glm-r-squared/
with(summary(glm_df1_spread_csp), 1 - deviance/null.deviance)


# 3. df2$CSP
glm_df2_CSP_csp = glm(formula=CSP~node_age, data=df2, family=binomial(link="logit"))
summary(glm_df2_CSP_csp)
logLik(glm_df2_CSP_csp)
AIC(glm_df2_CSP_csp)

# calculate McFadden's R-squared for GLM model with logit link function
# https://www.statology.org/glm-r-squared/
with(summary(glm_df2_CSP_csp), 1 - deviance/null.deviance)



#######################################################
# NOTE: AICs and lnLs can be compared on *identical* datasets, but not different datasets
#######################################################

