#------------------------------------------------------------------------
#
# 「TARGET声明と国際医学誌の実例からやさしく学ぶTarget Trial Emulation」
#  第10回日本糖尿病・生活習慣病ヒューマンデータ学会年次学術集会
#
# 野間久史　（統計数理研究所）
#
# 2025年12月14日
#
#------------------------------------------------------------------------

install.packages(“rqlm”)

library(rqlm)
head(exdata04, 10)

data(exdata04)    # Example dataset

fit_pp <- ttemsm(
formula = Y ~ A + L1 + time + I(time^2) + trial,
data = exdata04, id = ID, weight = w_pp,
eform = TRUE, cl = 0.95 )
# Pooled logistic regression for target trial emulation

fit_pp

# Call:
# ttemsm(formula = Y ~ A + L1 + time + I(time^2) + trial, data = exdata04, 
#     id = ID, weight = w_pp, eform = TRUE, cl = 0.95)
# 
# Coefficient estimates and CIs with robust SE estimator (cluster-robust SE):
#             Estimate Robust SE exp(coef)  Lower  Upper  z value Pr(>|z|)
# (Intercept)  -2.4142    0.2358    0.0894 0.0563 0.1420 -10.2363   0.0000
# A            -0.5190    0.2105    0.5951 0.3939 0.8990  -2.4656   0.0137
# L1            0.5050    0.0956    1.6570 1.3740 1.9983   5.2849   0.0000
# time         -0.0236    0.1292    0.9767 0.7582 1.2580  -0.1827   0.8550
# I(time^2)     0.0213    0.0188    1.0216 0.9846 1.0599   1.1357   0.2561
# trial2       -0.2638    0.2728    0.7681 0.4500 1.3110  -0.9672   0.3334
# trial3       -0.3587    0.2610    0.6986 0.4188 1.1651  -1.3744   0.1693
# trial4       -0.6249    0.2955    0.5353 0.3000 0.9552  -2.1150   0.0344
# trial5       -0.6327    0.2903    0.5311 0.3007 0.9382  -2.1798   0.0293
