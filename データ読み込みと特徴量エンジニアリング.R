#-----------pick up the file path--------------
path <- file.choose()
path

#-----------read csv file as compounds--------------
compounds <- read.csv(path)
View(compounds)

#-----------remove some columns if needed（今回は12列目を除去）--------------
trimed.compounds <- compounds[, -c(12)]

#-----------select rows without empty cells---------
is.completes <- complete.cases(trimed.compounds)
is.completes

complete.compounds <- trimed.compounds[is.completes,]
View(complete.compounds)

#-----------character to numeric(キャラクターを代理変数に置き換え)---------------------
labels_spacegroup <- as.numeric(as.factor(complete.compounds[, c(4)]))
labels_bandstructure <- as.numeric(as.factor(complete.compounds[, c(8)]))
labels_theoretical <- as.numeric(as.factor(complete.compounds[, c(11)]))
labels_crystal <- as.numeric(as.factor(complete.compounds[, c(13)]))

complete.compounds.characters <- as.matrix(cbind(labels_spacegroup, labels_bandstructure, labels_theoretical, labels_crystal))

#-----------select x from the dataset（complete.compoundsおよび、代理変数）-----------------
x <- cbind(complete.compounds[,c(5,6,9,10,12,14:303)], complete.compounds.characters)
View(x)

#-----------calculate standard distribution of x------
x.sds <- apply(x, 2, sd)
x.mean <- apply(x, 2, mean)
#-----------remove columns of 0 distribution from x----
sd.is.not.0 <- x.sds != 0
x <- x[, sd.is.not.0]
View(x)

#-----------select y from the dataset（予測対象はバンドギャップとします）------------------
y <- complete.compounds[,c(7)]
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
