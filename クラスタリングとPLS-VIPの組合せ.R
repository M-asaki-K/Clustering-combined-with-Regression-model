library(ClusterR) #クラスタリング
library(readr) # データ読み込み
library(dplyr) # データ操作一般
library(plyr)
library(assertr) # データのチェック
library(rsample)　#サンプリング
library(genalg)　#遺伝的アルゴリズム
library(pls)　#PLS
library(e1071)　#SVR
library(kernlab) #マハラノビス距離計算に使用
library(iml)　#機械学習結果の可視化パッケージ
library(devtools)　#一般的各種ツール
library(parallelDist) #並列計算ツール
library(bigmemory)　#メモリ節約ツール
library(rBayesianOptimization)　#ベイズ最適化

#並列計算の設定（すでにやってある場合は実行不要）
pkgs <- c('foreach', 'doParallel')
lapply(pkgs, require, character.only = T)
#if you want to change the number of threads for the calculation, please change the value "detectCores()"
registerDoParallel(makeCluster(detectCores()))


#計算条件変数の設定
max.num.of.clusters <- 15 #BIC最小化によるクラスタ数選択において、最大いくつまで見るか
scale <- "BIC"　#AICも選べる
validation.method <- "LOO"　#PLSの交差検証方法（"CV"または"LOO"）
Cluster.sample.num.min <- 5
VIP.lbound <- 1.5

#クラスタリング
BIC.measured.kmeans <- as.matrix(Optimal_Clusters_KMeans(multi.regression.x.train, max.num.of.clusters, criterion = scale))
cluster.num <- which(BIC.measured.kmeans == BIC.measured.kmeans[BIC.measured.kmeans == min(BIC.measured.kmeans),])
Centroid.kmeans <- KMeans_rcpp(multi.regression.x.train, cluster.num, num_init = 20, max_iters = 100, initializer = "kmeans++")
Centroid.kmeans$centroids

#学習データに対するクラスタ算出
Clusters <- as.numeric(predict_KMeans(multi.regression.x.train, CENTROIDS = Centroid.kmeans$centroids, threads = detectCores()))
multi.regression.compounds.clusters <- cbind(multi.regression.compounds.train, Clusters)
multi.regression.x.clusters <- cbind(multi.regression.x.train, Clusters)

A <- c(1:cluster.num)

#各クラスタにおけるサンプル数をベクトルAに格納
for (n in 1:cluster.num){
A[n] <-  sum(Clusters == n)
}

#最もサンプル数の多いクラスタを抽出
max.samples.num <- which(A == max(A))
max.samples.num

#サンプル数が10未満のクラスタは、最大クラスタに編入（モデルの不安定化を防ぐ）
for (m in 1:cluster.num){
    if (sum(Clusters == m) < Cluster.sample.num.min){
    Clusters[Clusters == m] <- max.samples.num
  } 
}

#テストデータに対するクラスタ算出
Clusters.test <- as.numeric(predict_KMeans(multi.regression.x.test, CENTROIDS = Centroid.kmeans$centroids, threads = detectCores()))

#学習データにおいて最大クラスタに編入されたクラスタは、テストデータにおいても同じく最大クラスタに編入
for (mt in 1:cluster.num){
  if (sum(Clusters == mt) < Cluster.sample.num.min){
    Clusters.test[Clusters.test == mt] <- max.samples.num
  } 
}

Clusters.test
Clusters

multi.regression.compounds.clusters.test <- cbind(multi.regression.compounds.test, Clusters.test)
colnames(multi.regression.compounds.clusters.test) <- c(colnames(multi.regression.compounds.test), "Clusters")
multi.regression.x.clusters.test <- cbind(multi.regression.x.test, Clusters.test)
colnames(multi.regression.x.clusters.test) <- c(colnames(multi.regression.x.test), "Clusters")

View(multi.regression.compounds.clusters.test)


#PLSによるそれぞれのクラスタに対する予測モデル構築（並列計算）
CLUSTERMODEL <- foreach(i = unique(Clusters), .combine = rbind,.packages = c("pls", "plsVarSel"))%dopar%{
#  i = 2
  #-------------------training data----------------------------------------
  multi.regression.compounds.train.clusters <- multi.regression.compounds.clusters[multi.regression.compounds.clusters[, c(ncol(multi.regression.compounds.clusters))] == i, ]
  preprocessed.y.train.clusters <- multi.regression.compounds.train.clusters[,c(1)]
  multi.regression.x.train.clusters <- multi.regression.compounds.train.clusters[,-c(1)]
  #-----------------------test data----------------------------------------
  multi.regression.compounds.test.clusters <- matrix(multi.regression.compounds.clusters.test[multi.regression.compounds.clusters.test[, c(ncol(multi.regression.compounds.clusters.test))] == i, ], ncol = ncol(multi.regression.compounds.clusters.test))
  multi.regression.compounds.test.clusters <- as.data.frame(multi.regression.compounds.test.clusters)
  colnames(multi.regression.compounds.test.clusters) <- colnames(multi.regression.compounds.clusters.test)
  preprocessed.y.test.clusters <- multi.regression.compounds.test.clusters[,c(1)]
  multi.regression.x.test.clusters <- multi.regression.compounds.test.clusters[,-c(1)]
  
  #-----------transform into data frame--------------------------
  multi.regression.compounds.train.clusters <- as.data.frame(multi.regression.compounds.train.clusters)
#  View(multi.regression.compounds.train.clusters)
  
  #------------------------------------------------------------------------------------------
  library(plsVarSel)
  
  vip.selected <- bve_pls(preprocessed.y.train.clusters, multi.regression.x.train.clusters, ncomp = ncomp.onesigma, VIP.threshold = VIP.lbound) 
#  vip.selected
  
  multi.regression.x.train.clusters.vip <- multi.regression.x.train.clusters[,c(vip.selected$bve.selection)]
  multi.regression.x.test.clusters.vip <- multi.regression.x.test.clusters[,c(vip.selected$bve.selection)]
  multi.regression.compounds.train.clusters.vip <- cbind.data.frame(preprocessed.y.train.clusters,multi.regression.x.train.clusters.vip)
  multi.regression.compounds.test.clusters.vip <- cbind.data.frame(preprocessed.y.test.clusters,multi.regression.x.test.clusters.vip)
#  View(multi.regression.compounds.train.clusters.vip)
  #------------------------PLS training------------------------------------
  library(pls)
  
  compounds.plsr.clusters <- plsr(preprocessed.y.train.clusters~., data=multi.regression.compounds.train.clusters.vip, validation= validation.method)
  summary(compounds.plsr.clusters)
    if (length(which.min(compounds.plsr.clusters$validation$PRESS)) == 0) {                 # if ( 条件式 )
    plsr.predicted.y.clusters <- matrix(data = mean(preprocessed.y.train.clusters), nrow = nrow(multi.regression.compounds.train.clusters.vip), ncol = 1)                   #  条件式が TRUE  のときに実行される部分
    plsr.predicted.y.test.clusters <- matrix(data = mean(preprocessed.y.train.clusters), nrow = nrow(multi.regression.compounds.test.clusters.vip), ncol = 1)

      } else {
    ncomp.onesigma.clusters <- which.min(compounds.plsr.clusters$validation$PRESS)                   #  条件式が FALSE のときに実行される部分
    plsr.predicted.y.clusters <- predict(compounds.plsr.clusters)[, , ncomp.onesigma.clusters]
    plsr.predicted.y.test.clusters <- predict(compounds.plsr.clusters,newdata = multi.regression.x.test.clusters.vip)[,, ncomp.onesigma.clusters]
    
  }

  train <- cbind((1:nrow(multi.regression.compounds.train.clusters.vip)), plsr.predicted.y.clusters, preprocessed.y.train.clusters)
  colnames(train) <- c("V1", "plsr.predicted.y.clusters", "preprocessed.y.train.clusters")
  if (nrow(multi.regression.compounds.test.clusters.vip) != 0){
    
    test <- cbind((1:nrow(multi.regression.compounds.test.clusters.vip)), plsr.predicted.y.test.clusters, preprocessed.y.test.clusters)
    colnames(test) <- c("V1", "plsr.predicted.y.test.clusters", "preprocessed.y.test.clusters")
  } else {
    test <- as.matrix(cbind(NA,NA,NA))
#    View(test)
    colnames(test) <- c("V1", "plsr.predicted.y.test.clusters", "preprocessed.y.test.clusters")
  }
  
  result <-merge (train, test, by = "V1", all = T)
  data.frame(result)
}

#予測-実測プロット
limit <- c(-2,4)

plot(CLUSTERMODEL[, c(3)], CLUSTERMODEL[, c(2)], xlim = limit, ylim = limit, xlab = "", ylab = "")
par(new = T)
plot(CLUSTERMODEL[, c(5)], CLUSTERMODEL[, c(4)], col = "blue", pch = 2, xlim = limit, ylim = limit, xlab = "observed value", ylab = "predicted value")
abline(a=0, b=1)

is.completes.CLUSTERMODEL <- complete.cases(CLUSTERMODEL[, c(4,5)])
R2test.CLUSTERMODEL <- 1 - sum((CLUSTERMODEL[is.completes.CLUSTERMODEL, c(4)] - CLUSTERMODEL[is.completes.CLUSTERMODEL, c(5)])^2) / sum((mean(CLUSTERMODEL[is.completes.CLUSTERMODEL, c(5)]) - CLUSTERMODEL[is.completes.CLUSTERMODEL, c(4)])^2)
R2test.CLUSTERMODEL

