library(ClusterR) #クラスタリング
library(readr) # データ読み込み
library(dplyr) # データ操作一般
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
max.num.of.clusters <- 3
scale <- "BIC"

#クラスタリング
BIC.measured.kmeans <- as.matrix(Optimal_Clusters_KMeans(multi.regression.x, max.num.of.clusters, criterion = scale))
Centroid.kmeans <- KMeans_rcpp(multi.regression.x, which(BIC.measured.kmeans == BIC.measured.kmeans[BIC.measured.kmeans == min(BIC.measured.kmeans),]), num_init = 5, max_iters = 100, initializer = "kmeans++")
Centroid.kmeans$centroids

Clusters <- as.numeric(predict_KMeans(multi.regression.x, CENTROIDS = Centroid.kmeans$centroids, threads = detectCores()))
multi.regression.compounds.clusters <- cbind(multi.regression.compounds, Clusters)
multi.regression.x.clusters <- cbind(multi.regression.x, Clusters)

TRAIN <- foreach(i = 1:max.num.of.clusters, .combine = rbind,.packages = c("pls"))%dopar%{                        

  compounds.clusters <- multi.regression.compounds.clusters[multi.regression.compounds.clusters[, c(ncol(multi.regression.compounds.clusters))] == i, ]
  x.clusters <- multi.regression.x.clusters[multi.regression.x.clusters[, c(ncol(multi.regression.x.clusters))] == i, ]
  #--------------------divide into test and training data----------------------
  train_size = 0.7
  
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
  
  compounds.plsr.clusters <- plsr(preprocessed.y~., data=multi.regression.compounds.train.clusters, validation="CV")
  summary(compounds.plsr.clusters)
  
  ncomp.onesigma.clusters <- which.min(compounds.plsr.clusters$validation$PRESS)

  plsr.predicted.y.clusters <- predict(compounds.plsr.clusters)[, , ncomp.onesigma.clusters]
  data.frame(plsr.predicted.y.clusters, preprocessed.y.train.clusters)
}

plot(TRAIN[, c(2)], TRAIN[, c(1)])

TEST <- foreach(i = 1:max.num.of.clusters, .combine = rbind,.packages = c("pls"))%dopar%{                        
  
  compounds.clusters <- multi.regression.compounds.clusters[multi.regression.compounds.clusters[, c(ncol(multi.regression.compounds.clusters))] == i, ]
  x.clusters <- multi.regression.x.clusters[multi.regression.x.clusters[, c(ncol(multi.regression.x.clusters))] == i, ]
  #--------------------divide into test and training data----------------------
  train_size = 0.7
  
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
  
  compounds.plsr.clusters <- plsr(preprocessed.y~., data=multi.regression.compounds.train.clusters, validation="CV")
  summary(compounds.plsr.clusters)
  plot(compounds.plsr.clusters, "validation")
  
  ncomp.onesigma.clusters <- which.min(compounds.plsr.clusters$validation$PRESS)
  
  plsr.predicted.y.clusters <- predict(compounds.plsr.clusters)[, , ncomp.onesigma.clusters]
  pls.predicted.y.test.clusters <- predict(compounds.plsr.clusters,newdata = multi.regression.x.test.clusters)[,, ncomp.onesigma.clusters]
  
  data.frame(pls.predicted.y.test.clusters, preprocessed.y.test.clusters)
}

plot(TEST[, c(2)], TEST[, c(1)])
  
