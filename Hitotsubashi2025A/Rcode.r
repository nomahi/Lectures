#------------------------------------------------------------------------
#
# 「メタアナリシス」
#  一橋大学ソーシャルデータサイエンス学部
#
# 野間久史　（統計数理研究所）
#
# 2025年11月13日
#
#------------------------------------------------------------------------



#--------------------------------------------------------------------------------
# 1. データセットとパッケージの読み込み
#--------------------------------------------------------------------------------


setwd("D:\\")

yusuf <- read.csv(file="yusuf1985.csv")			# beta-blockerの事例データの読み込み
teo <- read.csv(file="teo1991.csv")				# マグネシウム静注の事例データの読み込み

##

pkgCheck <- function(pkg){
	if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
	library(pkg, character.only = TRUE)
}

pkgCheck("meta")			# meta, metafor packages; パッケージの読み込み。もしインストールされていなかったらインストールする。
pkgCheck("metafor")





#--------------------------------------------------------------------------------
# 2. メタアナリシスの基本的な方法（Yusufらによるβブロッカーのメタアナリシス）
#--------------------------------------------------------------------------------

btb1 <- metabin(ai, n1i, ci, n2i, data=yusuf, sm="OR", studlab=trial, method="Peto")		# 2値アウトカムのメタアナリシス; PetoのOne-step法
summary(btb1)

forest(btb1, layout="RevMan5")							# フォレストプロットの作図






#--------------------------------------------------------------------------------
# 3. メタアナリシスの基本的な方法（Teoらによるマグネシウム静脈注射のメタアナリシス）
#--------------------------------------------------------------------------------

mteo1 <- metabin(x1i, n1i, x2i, n2i, data=teo, sm="OR", studlab=paste(author, year), method="Peto")			# 2値アウトカムのメタアナリシス
summary(mteo1)			# メタアナリシスのすべての結果を表示

forest(mteo1, layout="RevMan5")

mteo2 <- update(mteo1, method="Inverse")		# 固定効果モデルの解析を、逆分散法に変更
mteo3 <- update(mteo1, method="MH")				# 固定効果モデルの解析を、Mantel-Haenszel法に変更

summary(mteo2)			# メタアナリシスのすべての結果を表示
summary(mteo3)			# メタアナリシスのすべての結果を表示

mteo4 <- update(mteo3, method.tau="DL")			# 変量効果モデルの解析を、DerSimonian-Laird法に変更
mteo5 <- update(mteo3, method.tau="REML")		# 変量効果モデルの解析を、REML法に変更

summary(mteo4)			# メタアナリシスのすべての結果を表示
summary(mteo5)			# メタアナリシスのすべての結果を表示

mteo6 <- update(mteo5, sm="RR")					# 効果の指標をリスク比に変更
forest(mteo6, layout="RevMan5")

mteo7 <- update(mteo5, sm="RD")					# 効果の指標をリスク差に変更
forest(mteo7, layout="RevMan5")

mteo8 <- update(mteo5, method.predict="HTS")		# 予測区間の結果を追加
summary(mteo8)

mteo9 <- update(mteo1, method="GLMM")				# 変量効果モデルの解析を、ロジスティック混合効果モデルに変更
summary(mteo9)






#--------------------------------------------------------------------------------
# 4. ファンネルプロットと出版バイアスの検定
#--------------------------------------------------------------------------------

funnel(mteo8, studlab=TRUE)						# ファンネルプロットを作図
	   
metabias(mteo8, method.bias="Egger", k.min=7)   # Egger検定を行う（通常は、10試験以上でないと実行されない；k.min=7 とすることで、強制的に実行させている）



