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
library(“rqlm”)
head(exdata04, 10)

data(exdata04)    # Example dataset

fit_pp <- ttemsm(
formula = Y ~ A + L1 + time + I(time^2) + trial,
data = exdata04, id = ID, weight = w_pp,
eform = TRUE, cl = 0.95 )
# Pooled logistic regression for target trial emulation
