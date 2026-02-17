#  因果フォレスト入門
#  An Introduction to Causal Forest
#
#  野間久史（統計数理研究所）
#
#  日本臨床疫学会第8回年次学術大会
#
#  2026年2月21日
###


rm(list = ls())
set.seed(3542)  		# 乱数の再現性：参加者全員で同じ結果が出やすいように固定

install.packages("grf")			# 外部パッケージのインストール（初回のみ）
install.packages("maq")

library("grf")
library("maq") 		# For Qini curves.
library("ggplot2")

# ------------------------------------------------------------
# 0) データ読み込みと、Y（アウトカム）, W（介入）, X（共変量）の準備
# ------------------------------------------------------------

data = read.csv("https://raw.githubusercontent.com/grf-labs/grf/master/experiments/ijmpr/synthetic_data.csv")
Y = data$outcome
W = data$treatment
X = data[, -c(1, 2)]

if (is.null(colnames(X))) colnames(X) = make.names(1:ncol(X))

# 見るポイント：Y/W/X が想定どおりか（Wが0/1か、Xの列数・列名はあるか）を最初に確認
# 落とし穴：Wが {0,1} でない／欠測が多い／Xにカテゴリが混ざる（forestは基本“数値”前提）と挙動が変わる


# ============================================================
# 1) ATE（平均治療効果）を推定する
# ============================================================

# *** Estimating an average treatment effect (ATE) ***

summary(lm(Y ~ W))
# 見るポイント：lmのW係数は「交絡を無視した平均差」＝因果効果と一致するとは限らない（符号や大きさの参考程度）
# 落とし穴：この値を「因果効果」と言い切る（観察研究では特に危険：交絡で簡単にひっくり返ります）

cf.full = causal_forest(X, Y, W)
average_treatment_effect(cf.full)
# 見るポイント：ATE推定値と95%CI（標準誤差）。“交絡調整後にどれくらい変わるか” をlmと比較
# 落とし穴：ATEが有意/非有意だけで結論を作る（サンプルサイズとoverlapでSEが大きく動きます）

hist(cf.full$W.hat, xlab = "Estimated propensity scores", main = "")
# 見るポイント：propensity score が0や1に張り付いていないか（極端な値が多いと overlap が悪い）
# 落とし穴：overlapが悪いのに解析を続ける（推定が不安定・“その領域の効果”はデータが支えていない）


# ============================================================
# 2) CATE（条件付き治療効果）を推定する：個人差/異質性の推定
# ============================================================

# *** Estimating CATEs ***

train = sample(nrow(X), 0.6 * nrow(X))
test = -train
# 見るポイント：train/test split は「学習」と「評価」を分けるため（CATEは特に過学習しやすい）
# 落とし穴：同じデータで学習→同じデータで“当たり前に良い評価”をしてしまう（楽観バイアス）

cate.forest = causal_forest(X[train, ], Y[train], W[train])

X.test = X[test, ]
tau.hat.test = predict(cate.forest, X.test)$predictions

hist(tau.hat.test, xlab = "Estimated CATEs", main = "")
# 見るポイント：CATE推定値の分布（広がりがある＝異質性の“候補”だが、それだけで確定ではない）
# 落とし穴：CATEの点推定を個人レベルで断言する（点推定はノイズが大きく、ランキングで評価するのが基本）


# ------------------------------------------------------------
# 2-1) “推定CATEの分位” ごとにグループを作って、グループ別ATEを見る
# ------------------------------------------------------------

num.groups = 4
quartile = cut(tau.hat.test,
               quantile(tau.hat.test, seq(0, 1, by = 1 / num.groups)),
               labels = 1:num.groups,
               include.lowest = TRUE)
samples.by.quartile = split(seq_along(quartile), quartile)

eval.forest = causal_forest(X.test, Y[test], W[test])
# 見るポイント：評価専用forest（eval.forest）を別に作り、“学習で使った情報”の持ち込みを減らす
# 落とし穴：cate.forestのままグループATEを出す（評価が甘くなりやすい）

ate.by.quartile = lapply(samples.by.quartile, function(samples) {
  average_treatment_effect(eval.forest, subset = samples)
})

df.plot.ate = data.frame(
  matrix(unlist(ate.by.quartile), num.groups, byrow = TRUE, dimnames = list(NULL, c("estimate","std.err"))),
  group = 1:num.groups
)

ggplot(df.plot.ate, aes(x = group, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = estimate - 1.96 * std.err, ymax = estimate + 1.96 * std.err, width = 0.2)) +
  xlab("Estimated CATE quantile") +
  ylab("Average treatment effect")
# 見るポイント：右（高CATE群）ほどATEが大きい“単調な形”になっているか（ランキングが効いているサイン）
# 落とし穴：誤差バーが大きいのに形を断定する／分割1回の偶然で「異質性がある」と言い切る


# ============================================================
# 3) TOC / AUTOC による “異質性（ランキング）の評価”
# ============================================================

# *** Evaluate heterogeneity via TOC/AUTOC ***

rate.cate = rank_average_treatment_effect(
  eval.forest,
  tau.hat.test,
  q = seq(0.05, 1, length.out = 100)
)

plot(rate.cate)
# 見るポイント：TOC曲線が“上に膨らむ”ほど、上位（高CATE）を狙う価値がある（異質性の信号）
# 落とし穴：TOCの形だけで万能判断（overlap悪い・サンプル小さいと形が不安定になりやすい）

print(rate.cate)
# 見るポイント：AUTOC（estimate）とstd.err（不確実性）。estimateが0に近いなら“狙う価値が小さい”可能性
# 落とし穴：estimateの符号・大きさを文脈抜きで解釈（アウトカムの符号や定義で意味が変わります）

2 * pnorm(-abs(rate.cate$estimate / rate.cate$std.err))
# 見るポイント：RATE=0（異質性なし）に対する目安のp値（ただし “検定の勝ち負け” だけにしない）
# 落とし穴：p値で全て決める（複数比較・設計・実務上の閾値を無視しがち）


# ============================================================
# 4) 政策評価（ターゲティング）と Qini curve
# ============================================================

# *** Policy evaluation with Qini curves ****

cost = 1
qini = maq(tau.hat.test, cost, get_scores(eval.forest), R = 200)
qini.baseline = maq(tau.hat.test, cost, get_scores(eval.forest), R = 200,
                    target.with.covariates = FALSE)
# 見るポイント：qini（CATEで狙う）とbaseline（均等）を比較して、どれくらい上積みできるか
# 落とし穴：R（ブート回数）が小さすぎてCIが不安定／逆に大きすぎてハンズオン時間が破綻

max.deployment = 2000

scale_maq(qini, max.deployment) |>
  plot(ylab = "PTSD cases prevented",
       xlab = "Units held back from deployment")

scale_maq(qini.baseline, max.deployment) |>
  plot(add = TRUE, ci.args = NULL)
# 見るポイント：CATE曲線がbaselineより上にある範囲が長いほど“ターゲティングの価値”が大きい
# 落とし穴：縦軸の単位（何が何件減るのか）を曖昧にしたまま意思決定の話に飛ぶ

average_gain(scale_maq(qini, max.deployment), 500)
difference_gain(scale_maq(qini, max.deployment),
                scale_maq(qini.baseline, max.deployment), 500)
# 見るポイント：「500人分の介入」で見積もれる効果（実務上の制約に合わせた“点”で比較）
# 落とし穴：グラフ全体の印象だけで結論（実運用は特定の予算/人数点で意思決定することが多い）


# ============================================================
# 5) “どんな特徴の人が高CATE/低CATEか？” を覗く
# ============================================================

# *** Describing the fit CATE function ****

varimp.cate = variable_importance(cate.forest)
ranked.variables = order(varimp.cate, decreasing = TRUE)
top.varnames = colnames(X)[ranked.variables[1:4]]
print(top.varnames)
# 見るポイント：上位変数は「CATE予測に寄与していそう」な候補（“説明の糸口” になる）
# 落とし穴：重要度＝因果的な効果修飾因子と断定（相関・代理変数・分布差でも重要度は上がります）

low = samples.by.quartile[[1]]
high = samples.by.quartile[[num.groups]]

df.lo = data.frame(
  covariate.value = unlist(as.vector(X.test[low, top.varnames])),
  covariate.name = rep(top.varnames, each = length(low)),
  cate.estimates = "Low"
)
df.hi = data.frame(
  covariate.value = unlist(as.vector(X.test[high, top.varnames])),
  covariate.name = rep(top.varnames, each = length(high)),
  cate.estimates = "High"
)
df.plot.hist = rbind(df.lo, df.hi)

ggplot(df.plot.hist, aes(x = covariate.value, fill = cate.estimates)) +
  geom_histogram(alpha = 0.7, position = "identity") +
  facet_wrap(~ covariate.name, scales = "free", ncol = 2)
# 見るポイント：High vs Low で分布が明確にズレる変数があるか（“どう違う人が効く？” の直感）
# 落とし穴：分布差＝メカニズムの証拠と誤解（これは探索。確証には追加検証が必要）


# ============================================================
# 6) BLP（Best Linear Projection）：効果修飾を “線形” に要約
# ============================================================

# *** Best linear projections (BLP) ****

blp.vars = c("X1", "X2", "X3")
best_linear_projection(cf.full, X[, blp.vars])
# 見るポイント：係数の符号と大きさ（どの方向の人ほど効果が大きい傾向か）＋標準誤差（不確実性）
# 落とし穴：非線形なtau(x)を線形で“完全に説明できる”と思う（BLPはあくまで要約）


# ============================================================
# 7) リスク予測（risk）と CATEターゲティングの比較
# ============================================================

# *** Risk vs CATE-based targeting ***

train.hi = train[W[train] == 0]

rf.risk = regression_forest(X[train.hi, ], 1 - Y[train.hi])
risk.hat.test = predict(rf.risk, X.test)$predictions
# 見るポイント：riskは“発症確率の高さ”、CATEは“介入でどれだけ変わるか”。似て非なる指標
# 落とし穴：「高リスク＝高効果」と決めつける（予防は効くがリスクは低い、など普通に起こる）

rate.risk = rank_average_treatment_effect(
  eval.forest,
  cbind(tau.hat.test, risk.hat.test)
)
plot(rate.risk)
print(rate.risk)
# 見るポイント：AUTOC(cate) と AUTOC(risk) と、差（cate-risk）。“狙う基準”の優劣が見える
# 落とし穴：riskモデルの学習対象（ここではW=0）など前提を忘れて一般化する

rate.risk$estimate + data.frame(lower = -1.96 * rate.risk$std.err,
                                upper = 1.96 * rate.risk$std.err,
                                row.names = rate.risk$target)
# 見るポイント：差のCIが0を跨ぐか（どちらが優位かの不確実性）
# 落とし穴：点推定の差だけで勝敗を断定（CIが広いと結論が揺れます）


# ============================================================
# 8) Appendix：m(x), e(x) を別モデルで作って CATE を比較
# ============================================================

# *** Appendix: Evaluate CATE models via the AUTOC ***

Y.forest = regression_forest(X[train, ], Y[train], num.trees = 500)
Y.hat = predict(Y.forest)$predictions
# 見るポイント：m(x)（ベースライン）を別途推定する例。推定の分離で挙動がどう変わるかを観察
# 落とし穴：m/e/tau のどれを変えたのか曖昧に比較する（比較は“何を変えたか”が命）

W.forest = regression_forest(X[train, ], W[train], num.trees = 500)
W.hat = predict(W.forest)$predictions

varimp.Y = variable_importance(Y.forest)
selected.vars = which(varimp.Y >= quantile(varimp.Y, 0.75))
print(colnames(X)[selected.vars])
# 見るポイント：変数選択の一例（重要度上位のみ）。簡略化で性能が上がる/下がる可能性を体感
# 落とし穴：選択基準を固定化して“正解”と思う（データと目的で変わります）

if (length(selected.vars) <= 1) stop("You should really try and use more than just one predictor variable with forests.")

X.subset = X[, selected.vars]
cate.forest.restricted = causal_forest(X.subset[train, ], Y[train], W[train],
                                       Y.hat = Y.hat, W.hat = W.hat)

tau.hat.test.restricted = predict(cate.forest.restricted, X.test[, selected.vars])$predictions

rate.cate.compare = rank_average_treatment_effect(
  eval.forest,
  cbind(tau.hat.test, tau.hat.test.restricted)
)
print(rate.cate.compare)
# 見るポイント：AUTOCの比較で、どちらのCATEモデルが“ターゲティングに向くか”を判断
# 落とし穴：わずかな差を過大解釈（SEも見て、差のCIやp値も併記して判断する）

data.frame(
  p.value = 2 * pnorm(-abs(rate.cate.compare$estimate / rate.cate.compare$std.err)),
  target = rate.cate.compare$target
)

rate.cate.compare$estimate + data.frame(lower = -1.96 * rate.cate.compare$std.err,
                                        upper = 1.96 * rate.cate.compare$std.err,
                                        row.names = rate.cate.compare$target)
# 見るポイント：差（full - restricted）のCIが0を跨ぐか（“改善した”と言えるか）
# 落とし穴：この合成データの結果をそのまま実データへ一般化（実データはoverlap・欠測・測定誤差が効く）
