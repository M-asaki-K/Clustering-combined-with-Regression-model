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
validation.method <- "CV"　#PLSの交差検証方法（"CV"または"LOO"）
train.ratio <- 0.7

#クラスタリング
BIC.measured.kmeans <- as.matrix(Optimal_Clusters_KMeans(multi.regression.x, max.num.of.clusters, criterion = scale))
cluster.num <- which(BIC.measured.kmeans == BIC.measured.kmeans[BIC.measured.kmeans == min(BIC.measured.kmeans),])
Centroid.kmeans <- KMeans_rcpp(multi.regression.x, cluster.num, num_init = 20, max_iters = 100, initializer = "kmeans++")
Centroid.kmeans$centroids

Clusters <- as.numeric(predict_KMeans(multi.regression.x, CENTROIDS = Centroid.kmeans$centroids, threads = detectCores()))
multi.regression.compounds.clusters <- cbind(multi.regression.compounds, Clusters)
multi.regression.x.clusters <- cbind(multi.regression.x, Clusters)

#PLSによるそれぞれのクラスタに対する予測モデル構築（並列計算）
CLUSTERMODEL <- foreach(i = 1:cluster.num, .combine = rbind,.packages = c("pls"))%dopar%{
#  i = 1
  compounds.clusters <- multi.regression.compounds.clusters[multi.regression.compounds.clusters[, c(ncol(multi.regression.compounds.clusters))] == i, ]
  x.clusters <- multi.regression.x.clusters[multi.regression.x.clusters[, c(ncol(multi.regression.x.clusters))] == i, ]

  #--------------------divide into test and training data----------------------
  train_size = train.ratio
  
  n = nrow(compounds.clusters)
  #------------collect the data with n*train_size from the dataset------------
  perm = sample(n, size = round(n * train_size))
  
  #-------------------training data----------------------------------------
  multi.regression.compounds.train.clusters <- compounds.clusters[perm, ]
  preprocessed.y.train.clusters <- multi.regression.compounds.train.clusters[,c(1)]
  multi.regression.x.train.clusters <- multi.regression.compounds.train.clusters[,-c(1)]
  #-----------------------test data----------------------------------------
  multi.regression.compounds.test.clusters <- compounds.clusters[-perm, ]
  preprocessed.y.test.clusters <- multi.regression.compounds.test.clusters[,c(1)]
  multi.regression.x.test.clusters <- multi.regression.compounds.test.clusters[,-c(1)]
  
  #-----------transform into data frame--------------------------
  multi.regression.compounds.train.clusters <- as.data.frame(multi.regression.compounds.train.clusters)
  
  #------------------------PLS training------------------------------------
  library(pls)
  
  compounds.plsr.clusters <- plsr(preprocessed.y~., data=multi.regression.compounds.train.clusters, validation= validation.method)
  summary(compounds.plsr.clusters)
    if (length(which.min(compounds.plsr.clusters$validation$PRESS)) == 0) {                 # if ( 条件式 )
    plsr.predicted.y.clusters <- matrix(data = mean(preprocessed.y.train.clusters), nrow = nrow(multi.regression.compounds.train.clusters), ncol = 1)                   #  条件式が TRUE  のときに実行される部分
    plsr.predicted.y.test.clusters <- matrix(data = mean(preprocessed.y.train.clusters), nrow = nrow(multi.regression.compounds.test.clusters), ncol = 1)

      } else {
    ncomp.onesigma.clusters <- which.min(compounds.plsr.clusters$validation$PRESS)                   #  条件式が FALSE のときに実行される部分
    plsr.predicted.y.clusters <- predict(compounds.plsr.clusters)[, , ncomp.onesigma.clusters]
    plsr.predicted.y.test.clusters <- predict(compounds.plsr.clusters,newdata = multi.regression.x.test.clusters)[,, ncomp.onesigma.clusters]
    
  }

  train <- cbind((1:nrow(multi.regression.compounds.train.clusters)), plsr.predicted.y.clusters, preprocessed.y.train.clusters)
  colnames(train) <- c("V1", "plsr.predicted.y.clusters", "preprocessed.y.train.clusters")

  test <- cbind((1:nrow(multi.regression.compounds.test.clusters)), plsr.predicted.y.test.clusters, preprocessed.y.test.clusters)
  colnames(test) <- c("V1", "plsr.predicted.y.test.clusters", "preprocessed.y.test.clusters")
  
  result <-merge (train, test, by = "V1", all = T)
  data.frame(result)
}

#予測-実測プロット
limit <- c(-2,8)

plot(CLUSTERMODEL[, c(3)], CLUSTERMODEL[, c(2)], xlim = limit, ylim = limit, xlab = "", ylab = "")
par(new = T)
plot(CLUSTERMODEL[, c(5)], CLUSTERMODEL[, c(4)], col = "blue", pch = 2, xlim = limit, ylim = limit, xlab = "observed value", ylab = "predicted value")
abline(a=0, b=1)

is.completes.CLUSTERMODEL <- complete.cases(CLUSTERMODEL[, c(4,5)])
R2test.CLUSTERMODEL <- 1 - sum((CLUSTERMODEL[is.completes.CLUSTERMODEL, c(4)] - CLUSTERMODEL[is.completes.CLUSTERMODEL, c(5)])^2) / sum((mean(CLUSTERMODEL[is.completes.CLUSTERMODEL, c(5)]) - CLUSTERMODEL[is.completes.CLUSTERMODEL, c(4)])^2)
R2test.CLUSTERMODEL
