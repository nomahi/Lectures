##############################################
# 0. 準備：パッケージ読み込み ＆ データ生成
##############################################

# 必要ならインストール
# install.packages("TrialEmulation")

library(TrialEmulation)

set.seed(12345)

#---------------------------------------------
# 0-1. 教材用のシンプルな TTE データをシミュレート
#   - 個人ID: id
#   - 離散時点: period (t = 0,1,...,K)
#   - 治療: A_t (0/1)
#   - アウトカム: Y_t (イベント発生フラグ; once 1 then stay 1)
#   - センシング: C_t (ここではかなり単純)
#   - Eligible: すべての t=0 で eligible=1（簡略化）
#   - ベースライン共変量: X3, X4, age_s
#   - 時間依存共変量: X1_t, X2_t
#---------------------------------------------

simulate_simple_tte <- function(N = 1000, K = 9) {
  id     <- rep(1:N, each = K + 1)   # 各個人を (K+1) 期間分に複製
  period <- rep(0:K, times = N)      # t = 0,...,K
  
  # ベースライン共変量（個人単位）
  X3_i  <- rbinom(N, 1, 0.5)         # 例: ベースラインのバイナリ交絡
  X4_i  <- rnorm(N, 0, 1)            # 例: ベースラインの連続交絡
  age_i <- rnorm(N, 60, 8)           # 年齢
  age_s_i <- as.numeric(scale(age_i))# 標準化年齢
  
  # person-time に展開
  X3  <- rep(X3_i,  each = K + 1)
  X4  <- rep(X4_i,  each = K + 1)
  age_s <- rep(age_s_i, each = K + 1)
  
  # 時間依存共変量を初期化
  X1 <- numeric(N * (K + 1))  # バイナリ
  X2 <- numeric(N * (K + 1))  # 連続
  
  # 治療, アウトカム, センシング
  A <- numeric(N * (K + 1))
  Y <- numeric(N * (K + 1))
  C <- numeric(N * (K + 1))
  
  # 個人ごとに時系列を生成
  for (i in 1:N) {
    idx <- ((i - 1) * (K + 1) + 1):(i * (K + 1))
    
    # t = 0 の共変量・治療
    X1[idx[1]] <- rbinom(1, 1, plogis(-0.3 + 0.5 * X3_i[i]))
    X2[idx[1]] <- rnorm(1, mean = 0.3 * X4_i[i], sd = 1)
    
    # t = 0 の治療（ベースラインの new user / non-user）
    A[idx[1]] <- rbinom(1, 1, plogis(-0.2 + 0.8 * X3_i[i] + 0.5 * age_s_i[i]))
    
    # ベースライン時点でのイベント・検閲は 0 にしておく
    Y[idx[1]] <- 0
    C[idx[1]] <- 0
    
    # t >= 1 のプロセス
    for (k in 1:K) {
      j <- idx[k + 1]  # period k に対応
  
      # 前時点のデータ
      j_prev <- idx[k]
      
      # すでにイベント or 検閲なら、その後は観測しない（Y,C は 0 のまま）
      if (Y[j_prev] == 1 || C[j_prev] == 1) {
        X1[j] <- X1[j_prev]
        X2[j] <- X2[j_prev]
        A[j]  <- A[j_prev]
        Y[j]  <- Y[j_prev]
        C[j]  <- 1  # 一度検閲されたら以降は欠測扱いとしたいので 1
        next
      }
      
      # 時間依存共変量（前回 + ノイズ）
      X1[j] <- rbinom(1, 1, plogis(-0.2 + 0.7 * X1[j_prev] + 0.3 * X3_i[i]))
      X2[j] <- rnorm(1, mean = 0.7 * X2[j_prev] + 0.2 * X4_i[i], sd = 0.7)
      
      # 治療（スイッチングを許す）
      #   - 前回の治療・共変量に依存
      pA <- plogis(-1 + 1.0 * A[j_prev] + 0.7 * X1[j] + 0.4 * age_s_i[i])
      A[j] <- rbinom(1, 1, pA)
      
      # イベントハザード（pooled logistic に対応する離散ハザード）
      #   - 治療 (A)、共変量 (X1, X2, X3, X4, age_s) に依存
      lin_haz <- -3.0 + 0.8 * A[j] + 0.7 * X1[j] + 0.4 * X3_i[i] + 0.3 * age_s_i[i]
      pY <- plogis(lin_haz)
      Y[j] <- rbinom(1, 1, pY)
      
      # 検閲（治療や共変量にやや依存）
      pC <- plogis(-3.5 + 0.3 * A[j] + 0.3 * abs(X2[j]))
      C[j] <- rbinom(1, 1, pC)
    }
  }
  
  # eligible：ここでは t=0 のみ 1、それ以外 0（simple）
  eligible <- ifelse(period == 0, 1, 0)
  
  data.frame(
    id = id,
    period = period,
    A = A,
    Y = Y,
    C = C,
    eligible = eligible,
    X1 = X1,
    X2 = X2,
    X3 = X3,
    X4 = X4,
    age_s = age_s
  )
}

simdata <- simulate_simple_tte()
head(simdata)


##############################################
# 1. ITT 解析（intention-to-treat 効果）
##############################################

#---------------------------------------------
# 1-1. data_preparation()
#   - 逐次トライアルを構築し、ITT 用の拡張データを作成
#   - 依存検閲に対する IPCW（逆確率検閲重み）も計算
#   - ここでは論文と同様に estimand_type = "ITT" を指定
#     （TrialEmulation JSS 論文の Section 4.2 に対応）:contentReference[oaicite:1]{index=1}
#---------------------------------------------

prep_ITT <- data_preparation(
  data     = simdata,
  id       = "id",        # 個体ID
  period   = "period",    # 離散時間 t
  treatment= "A",         # 観察された治療
  outcome  = "Y",         # イベント指標
  eligible = "eligible",  # ターゲットトライアルへの適格性
  estimand_type = "ITT",  # ITT 効果を推定
  
  # MSM に入れるベースライン共変量（trial baseline 時点の値に固定される）
  outcome_cov = ~ X1 + X2 + X3 + X4 + age_s,
  
  # MSM の中で「割り付けられた治療」を説明変数として使う
  model_var = "assigned_treatment",
  
  # 依存検閲を IPCW で補正（単純なロジスティックモデル）
  use_censor_weights = TRUE,
  cense      = "C",
  cense_d_cov = ~ X1 + X2 + X3 + X4 + age_s,  # 分母モデル（真面目に入れる）
  cense_n_cov = ~ X3 + X4,                    # 分子モデル（簡略）
  pool_cense  = "numerator",                  # treated/untreated で分子だけプール
  glm_function = "glm",                       # 並列化は使わずシンプルに
  quiet = TRUE
)

# どんな変数ができているか確認
str(prep_ITT$data)
head(prep_ITT$data)


#---------------------------------------------
# 1-2. trial_msm()
#   - 拡張データに対して、pooled logistic 回帰で
#     MSM（marginal structural model）をフィットする
#   - 離散ハザード ~ 治療（assigned_treatment） + baseline 共変量 + 時間
#---------------------------------------------

msm_ITT <- trial_msm(
  data = prep_ITT,
  estimand_type = "ITT",
  
  # MSM の説明変数（線形予測子に入る項）
  #   - include_followup_time, include_trial_period で
  #     ベースラインハザードと trial 間の違いを表現
  outcome_cov = ~ assigned_treatment + X1 + X2 + X3 + X4 + age_s,
  
  include_followup_time = ~ followup_time + I(followup_time^2),
  include_trial_period  = ~ trial_period,
  
  glm_function = "glm"
)

summary(msm_ITT)  # 治療効果（ロジットスケール）などを確認


#---------------------------------------------
# 1-3. predict()
#   - 指定したターゲットトライアル集団について、
#     「全員 Treat」と「全員 Non-treat」の下での
#     累積発生確率（cum_inc）を予測
#   - その差が ITT 型の因果効果
#---------------------------------------------

# ターゲットトライアル集団として、
# 「最初の trial_period の baseline に適格だった人たち」を例にする
baseline_dat_ITT <- subset(
  prep_ITT$data,
  followup_time == 0 & trial_period == min(trial_period)
)

# MSM に必要な共変量だけ残す（重複行は削除）
baseline_dat_ITT <- unique(
  baseline_dat_ITT[, c("X1", "X2", "X3", "X4", "age_s", "trial_period")]
)

# 予測したい追跡時間（例: t = 0〜9）
time_points <- 0:9

pred_ITT <- predict(
  object        = msm_ITT,
  newdata       = baseline_dat_ITT,
  predict_times = time_points,
  type          = "cum_inc",  # 累積発生確率
  conf_int      = TRUE        # 95% CI も計算
)

# 結果の構造は
#   - treatment: 0/1
#   - time: 追跡時間
#   - estimate: 累積リスク
#   - lower/upper: CI
#   - diff: treated - untreated（行による）
pred_ITT

# 教材用に、例えば t=5 の ITT 効果だけ抜き出す例
ITT_t5 <- subset(pred_ITT, time == 5 & type == "difference")
ITT_t5


##############################################
# 2. Per-protocol 解析（PP 効果）
##############################################

#---------------------------------------------
# 2-1. data_preparation()（PP 用）
#   - estimand_type = "PP"
#   - 治療スイッチが起きた時点で「人工的打ち切り」を行い、
#     そのバイアスを IPTW（治療アドヒアランスに対する逆確率重み）
#     で補正する設定
#   - JSS 論文 Section 2.1.4 の実装イメージに対応 :contentReference[oaicite:2]{index=2}
#---------------------------------------------

prep_PP <- data_preparation(
  data     = simdata,
  id       = "id",
  period   = "period",
  treatment= "A",
  outcome  = "Y",
  eligible = "eligible",
  estimand_type = "PP",     # ★ここだけ変更
  
  outcome_cov = ~ X1 + X2 + X3 + X4 + age_s,
  model_var   = "assigned_treatment",
  
  # 治療スイッチに対する逆確率重み（分母・分子の共変量）
  switch_d_cov = ~ time_on_regime + I(time_on_regime^2) + X1 + X2 + X3 + X4 + age_s,
  switch_n_cov = ~ time_on_regime + I(time_on_regime^2) + X3 + X4 + age_s,
  
  # 検閲に対する逆確率重み（ITT と同様）
  use_censor_weights = TRUE,
  cense       = "C",
  cense_d_cov = ~ X1 + X2 + X3 + X4 + age_s,
  cense_n_cov = ~ X3 + X4,
  pool_cense  = "numerator",
  
  glm_function = "glm",
  quiet = TRUE
)

str(prep_PP$data)
head(prep_PP$data)


#---------------------------------------------
# 2-2. trial_msm()（PP 用）
#   - ここでは、割り付け治療の効果が「アドヒアランスを維持した場合」
#     の潜在アウトカムに対応するよう、重み付き pooled logistic MSM をフィット
#---------------------------------------------

msm_PP <- trial_msm(
  data = prep_PP,
  estimand_type = "PP",
  
  outcome_cov = ~ assigned_treatment + X1 + X2 + X3 + X4 + age_s,
  include_followup_time = ~ followup_time + I(followup_time^2),
  include_trial_period  = ~ trial_period,
  
  glm_function = "glm"
)

summary(msm_PP)


#---------------------------------------------
# 2-3. predict() で PP 効果（累積リスク差）を取得
#---------------------------------------------

baseline_dat_PP <- subset(
  prep_PP$data,
  followup_time == 0 & trial_period == min(trial_period)
)

baseline_dat_PP <- unique(
  baseline_dat_PP[, c("X1", "X2", "X3", "X4", "age_s", "trial_period")]
)

pred_PP <- predict(
  object        = msm_PP,
  newdata       = baseline_dat_PP,
  predict_times = time_points,
  type          = "cum_inc",
  conf_int      = TRUE
)

pred_PP

# t=5 における PP 効果（treat vs non-treat の累積リスク差）
PP_t5 <- subset(pred_PP, time == 5 & type == "difference")
PP_t5
