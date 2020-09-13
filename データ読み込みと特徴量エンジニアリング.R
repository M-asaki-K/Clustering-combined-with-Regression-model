#-----------remove some columns if needed（今回は8列目を除去）--------------
trimed.compounds <- compounds[, -c(8)]

#-----------select rows without empty cells---------
is.completes <- complete.cases(trimed.compounds)
is.completes

complete.compounds <- trimed.compounds[is.completes,]
View(complete.compounds)

#-----------select x from the dataset-----------------
x <- complete.compounds[,c(3, 4, 6:413)]
View(x)

#-----------calculate standard distribution of x------
x.sds <- apply(x, 2, sd)
x.mean <- apply(x, 2, mean)
#-----------remove columns of 0 distribution from x----
sd.is.not.0 <- x.sds != 0
x <- x[, sd.is.not.0]
View(x)

#-----------select y from the dataset（予測対象はバンドギャップとします）------------------
y <- complete.compounds[,c(5)]
y

#-----------standarization of y------------------------
preprocessed.y <- (y - mean(y)) / sd(y)
mean(preprocessed.y)
sd(preprocessed.y)

#-----------standarization of x（通常の正規化を使います）------------------------
apply(x, 2, mean)
apply(x, 2, sd)
preprocessed.x <- apply(x, 2, function(x) {(x - mean(x)) / sd(x)})
View(preprocessed.x)

#-----------x converted into data frame type for machine learning-----------
preprocessed.x <- data.frame(preprocessed.x)

