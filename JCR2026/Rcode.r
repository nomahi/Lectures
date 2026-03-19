######################################################################
# R/RStudioによる臨床研究の統計解析：傾向スコアによる多変量解析を事例に
#
# 第70回日本リウマチ学会総会・学術集会
# 2026年4月24日
#
# 野間 久史
# 統計数理研究所
#
#
# 事例プログラム：Hydroxychloroquine と COVID-19 傾向スコア解析
#
# このスクリプトでは、以下の流れを学びます。
# 1. データの読み込み
# 2. 未調整の Kaplan-Meier 曲線と Cox 回帰
# 3. 傾向スコアの推定と IPTW の作成
# 4. 背景因子バランスの確認（数表と Love plot）
# 5. 重み付き Cox 回帰によるハザード比推定
# 6. 重み付き Kaplan-Meier 曲線の作図
# 7. tableone を用いた Table 1 の作成
#
#
# 主な変数：
#   x      : 治療群（0=No Hydroxychloroquine, 1=Hydroxychloroquine）
#   days   : 観察日数
#   cens   : イベント指標（1=挿管・死亡, 0=打ち切り）
######################################################################


############################################################
# 0. 使用するパッケージの読み込み
############################################################

install.packages(c("survival","WeightIt","cobalt","tableone","survey"), repos = "https://cloud.r-project.org")
# 必要なライブラリのインストール（初回のみ）

# survival:
#   生存時間解析の基本パッケージです。
#   Kaplan-Meier 曲線や Cox 回帰を行います。
library(survival)

# WeightIt:
#   傾向スコア重み付け（IPTW など）を行うためのパッケージです。
library(WeightIt)

# cobalt:
#   重み付け前後で背景因子のバランスが改善したかを確認するためのパッケージです。
#   bal.tab() で数表、love.plot() で図を作れます。
library(cobalt)

# tableone:
#   背景因子の Table 1 を作るためのパッケージです。
#   臨床研究でよく見られる形式の表を簡単に出力できます。
library(tableone)

# survey:
#   重み付きデータの記述統計を扱うためのパッケージです。
#   IPTW 後の Table 1 を作る際に使います。
library(survey)


############################################################
# 1. データの読み込み
############################################################

# 必要に応じて作業フォルダを指定してください。
# setwd("D:/IPTW")

# CSV ファイルを読み込みます。
excovid <- read.csv("excovid19.csv", stringsAsFactors = FALSE)

# データセットの中身を確認します。
# 各変数の型や、観測数を確認できます。
head(excovid,20)

# 治療群を factor 型に変換します。
# 図や表で 0/1 より群名のほうが見やすいためです。
excovid$trt <- factor(
  excovid$x,
  levels = c(0, 1),
  labels = c("No Hydroxychloroquine", "Hydroxychloroquine")
)


############################################################
# 2. 未調整解析：
#    交絡調整をしない場合の Kaplan-Meier 曲線と Cox 回帰
############################################################

# 治療群ごとの Kaplan-Meier 曲線を推定します。
# Surv(days, cens) は「観察日数」と「イベント発生の有無」を指定しています。
km1 <- survfit(Surv(days, cens) ~ trt, data = excovid)

# 未調整 Kaplan-Meier 曲線を描きます。
plot(
  km1,
  col = c("blue", "red"),
  lwd = 2,
  main = "Unadjusted KM",
  xlab = "Days",
  ylab = "Survival Probability",
  mark.time = FALSE
)

# 凡例を追加します。
legend(
  "bottomleft",
  legend = levels(excovid$trt),
  col = c("blue", "red"),
  lty = 1,
  lwd = 2,
  bty = "n"
)

# 未調整の Cox 回帰を行います。
# ここでは治療群 x のみを説明変数にしており、
# 背景因子による交絡調整はまだ行っていません。
ph1 <- coxph(Surv(days, cens) ~ x, data = excovid)

# 結果を表示します。
# exp(coef) がハザード比です。
summary(ph1)


############################################################
# 3. 傾向スコアの推定と IPTW の作成
############################################################

# 傾向スコアモデルに入れる背景因子（共変量）を指定します。
vars <- c(
  "age40", "bmi30", "smoke", "lungdisease", "hypertension",
  "cancer", "kidneydisease", "statin", "acearb", "steroids", "doac"
)

# weightit() を使って IPTW を作成します。
#
# x ~ ... の左辺は治療群、右辺は背景因子です。
# method = "glm" はロジスティック回帰による傾向スコア推定を意味します。
# estimand = "ATE" は平均処置効果（Average Treatment Effect）を対象にした
# 重みを作る設定です。
W.out <- weightit(
  x ~ age40 + bmi30 + smoke + lungdisease + hypertension + cancer +
    kidneydisease + statin + acearb + steroids + doac,
  data = excovid,
  method = "glm",
  estimand = "ATE"
)

# 重みの概要を確認します。
# 極端な重みがないか、ESS がどの程度かを確認できます。
summary(W.out)

# 作成した重みをデータセットに追加します。
excovid$w_ate <- W.out$weights


############################################################
# 4. tableone を用いた Table 1 の作成
############################################################

# Table 1 では、治療群間の背景因子分布を確認します。
# print(..., smd = TRUE) とすることで、
# 標準化差（SMD）もあわせて表示できます。

# まず、未調整の Table 1 を作成します。
tab_un <- CreateTableOne(
  vars = vars,
  strata = "trt",
  data = excovid
)

print(tab_un, smd = TRUE)

# 次に、IPTW 後の Table 1 を作成します。
# svydesign() で「重み付きデータ」として定義し、
# それを svyCreateTableOne() に渡します。
#
# 注意：
# この表の n は実人数ではなく、重み付き総数として表示されます。
# そのため、未調整の n と一致しなくても異常ではありません。
tab_w <- svyCreateTableOne(
  vars = vars,
  strata = "trt",
  data = survey::svydesign(
    ~1,
    weights = ~w_ate,
    data = excovid
  )
)

print(tab_w, smd = TRUE)

# 注：関数の仕様上、n の数が重み付け前のポピュレーションとは合致しないことがありますが、重み付けを行うことで、サンプルサイズが増減することはありません。
# 重み付け後の n の出力が誤りです。論文などで報告する際には、重み付け前の表の n を報告してください。


# Love plot は、背景因子ごとの SMD を図にしたものです。
# 重み付けによってバランスが改善したかを視覚的に確認できます。
love.plot(
  W.out,
  stats = "mean.diffs",
  abs = TRUE,
  thresholds = c(m = 0.1),
  var.order = "unadjusted",
  binary = "std",
  stars = "raw",
  line = TRUE
)


############################################################
# 5. IPTW 重み付き Cox 回帰によるハザード比推定
############################################################

# WeightIt の coxph_weightit() で、IPTW-Cox回帰による解析を行います。
# 重み推定を考慮した推定量を扱いやすい関数です。
ph2_w <- coxph_weightit(
  Surv(days, cens) ~ x,
  data = excovid,
  weightit = W.out
)

# transform = "exp" によって、回帰係数ではなくハザード比で表示します。
summary(ph2_w, ci = TRUE, transform = "exp")


############################################################
# 6. 重み付き Kaplan-Meier 曲線の作図
############################################################

par(mfrow=c(1,2))
# 重み付け前と後のKaplan-Meier曲線を、2枚並べて作図します。

# まず、未調整 Kaplan-Meier 曲線を描きます。
plot(
  km1,
  col = c("blue", "red"),
  lwd = 2,
  main = "Unadjusted KM",
  xlab = "Days",
  ylab = "Survival Probability",
  mark.time = FALSE
)

# 凡例を追加します。
legend(
  "bottomleft",
  legend = levels(excovid$trt),
  col = c("blue", "red"),
  lty = 1,
  lwd = 2,
  bty = "n"
)

# survfit() に weights = w_ate を与えることで、
# IPTW を反映した重み付き Kaplan-Meier 曲線を推定できます。
km_w <- survfit(
  Surv(days, cens) ~ trt,
  data = excovid,
  weights = w_ate
)

# 重み付き Kaplan-Meier 曲線を描きます。
# par(mfrow = c(1, 2)) により、右側のパネルに描かれます。
plot(
  km_w,
  col = c("blue", "red"),
  lwd = 2,
  main = "IPTW-adjusted KM",
  xlab = "Days",
  ylab = "Survival Probability",
  mark.time = FALSE
)

# 凡例を追加します。
legend(
  "bottomleft",
  legend = levels(excovid$trt),
  col = c("blue", "red"),
  lty = 1,
  lwd = 2,
  bty = "n"
)

