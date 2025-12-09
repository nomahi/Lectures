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

install.packages(“rqlm”)        # Install from CRAN

library(rqlm)
head(exdata04, 10)

data(exdata04)    # Example dataset; Long format

# Y：アウトカム（=0, 1）
# A：治療
# L1：交絡要因　～　簡単のため、１つだけだが複数でよい（formulaに + L2 + L3 + ... とする）
# time：追跡期間　～　３次以上の項／スプラインを入れてもよい
# trial：試験（Factor）
# ID：個人のID
# w_pp：IPWの重み

ttemsm( Y ~ A + L1 + time + I(time^2) + trial,
  data    = exdata04, id = ID, weight = w_pp,
  eform   = TRUE, cl = 0.95, var.method="standard")
# Pooled logistic regression for target trial emulation

# Call:
# ttemsm(formula = Y ~ A + L1 + time + I(time^2) + trial, data = exdata04, 
#     id = ID, weight = w_pp, eform = TRUE, cl = 0.95, var.method = "standard")
# 
# Coefficient estimates and CIs with cluster-robust SE estimator (standard):
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


ttemsm( Y ~ A + L1 + time + I(time^2) + trial,
  data    = exdata04, id = ID, weight = w_pp,
  eform   = TRUE, cl = 0.95, var.method="MBN")
# Pooled logistic regression for target trial emulation
# Morel-Bokossa-Neerchaal-type corrected SE estimator is used.

# Call:
# ttemsm(formula = Y ~ A + L1 + time + I(time^2) + trial, data = exdata04, 
#     id = ID, weight = w_pp, eform = TRUE, cl = 0.95, var.method = "MBN")
# 
# Coefficient estimates and CIs with cluster-robust SE estimator (MBN):
#             Estimate Robust SE exp(coef)  Lower  Upper  z value Pr(>|z|)
# (Intercept)  -2.4142    0.2388    0.0894 0.0560 0.1428 -10.1077   0.0000
# A            -0.5190    0.2125    0.5951 0.3924 0.9025  -2.4427   0.0146
# L1            0.5050    0.0967    1.6570 1.3710 2.0027   5.2235   0.0000
# time         -0.0236    0.1307    0.9767 0.7560 1.2618  -0.1805   0.8567
# I(time^2)     0.0213    0.0190    1.0216 0.9843 1.0603   1.1237   0.2611
# trial2       -0.2638    0.2753    0.7681 0.4478 1.3174  -0.9584   0.3378
# trial3       -0.3587    0.2638    0.6986 0.4165 1.1717  -1.3595   0.1740
# trial4       -0.6249    0.2984    0.5353 0.2983 0.9607  -2.0943   0.0362
# trial5       -0.6327    0.2932    0.5311 0.2990 0.9435  -2.1582   0.0309

# Note 1: ベースライン共変量（交絡要因）の調整は、IPTW解析としてもよい。その場合は、stabilized weightの使用が推奨されている（rqlm packageでは、stabwt関数で、stabilized weightの簡単な計算が可能です）。
# Note 2: ITT, PP解析の双方で、脱落に対してのIPCW解析を行うことが推奨されています。重みの推定には、ベースライン共変量を固定して、time-fixed variablesとして推定するのが現状のスタンダードです（rqlm packageでは、stabwtlong関数で、stabilized weightの計算が可能です）。
# Note 3: PP解析では、プロトコル逸脱に対してのIPCW解析を行う必要があります。重みの推定には、時間変化する変数を、time-varying variablesとしてモデル化するのが現状のスタンダードです（rqlm packageのstabwtlong関数で、stabilized weightの計算が可能です）。
