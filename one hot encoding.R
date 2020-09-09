library(mltools)
library(data.table)
library(dplyr)

#データベースの読み込み
path <- file.choose()
path

compounds <- read.csv(path)
compounds <- compounds[c(1:2000), ]
View(compounds)

#クラス変数の選択
factorclass <- c(3, 7, 10, 13)
compounds.chr <- as.data.frame(compounds[ , factorclass])
View(compounds.chr)

#factor型に変換
compounds.chr <- compounds.chr %>% mutate_if(function(col) !is.numeric(col), as.factor)

#data.tableに変換後、one hotエンコーディング
compounds.chr <- data.table(compounds.chr)
compounds.chr <- one_hot(compounds.chr, cols = "auto", sparsifyNAs = FALSE, naCols = FALSE, dropCols = TRUE, dropUnusedLevels = FALSE)
View(as.matrix(compounds.chr))

#元のデータセットと結合
compounds <- cbind(compounds[, -factorclass], compounds.chr)
View(compounds)
