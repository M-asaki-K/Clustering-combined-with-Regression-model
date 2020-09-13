#-----------compare the number of columns and rows--------------
ncol(preprocessed.x)
nrow(preprocessed.x)

#-----------pick up columns if needed---------------------------
multi.regression.x <- preprocessed.x[ , ]

#-----------definition of multi.regression.compounds (used for MLR)--------
multi.regression.compounds <- cbind(preprocessed.y, multi.regression.x)

#--------------------divide into test and training data----------------------
train_size = 0.7

n = nrow(multi.regression.compounds)
#------------collect the data with n*train_size from the dataset------------
perm = sample(n, size = round(n * train_size))

#-------------------training data----------------------------------------
multi.regression.compounds.train <- multi.regression.compounds[perm, ]
preprocessed.y.train <- multi.regression.compounds.train[,c(1)]
multi.regression.x.train <- multi.regression.compounds.train[,-c(1)]
#-----------------------test data----------------------------------------
multi.regression.compounds.test <-multi.regression.compounds[-perm, ]
preprocessed.y.test <- multi.regression.compounds.test[,c(1)]
multi.regression.x.test <- multi.regression.compounds.test[,-c(1)]

#-----------transform into data frame--------------------------
multi.regression.compounds.train <- as.data.frame(multi.regression.compounds.train)

#------------------------PLS training------------------------------------
library(pls)

compounds.plsr <- plsr(preprocessed.y~., data=multi.regression.compounds.train, validation="CV")
summary(compounds.plsr)
plot(compounds.plsr, "validation")

ncomp.onesigma <- selectNcomp(compounds.plsr, method = "randomization", plot = TRUE, ylim = c(.1, 1))
ncomp.onesigma

predict(compounds.plsr)[, , ncomp.onesigma]
plsr.predicted.y <- predict(compounds.plsr)[, , ncomp.onesigma]

#-------------------------pls test--------------------------------
pls.predicted.y.test <- predict(compounds.plsr,newdata = multi.regression.x.test)[,, ncomp.onesigma]

plot(multi.regression.compounds.train[,c(1)], plsr.predicted.y, xlim = c(-2, 6), ylim = c(-2,6), xlab = "", ylab = "")
par(new = T)
plot(preprocessed.y.test, pls.predicted.y.test, col = "blue", pch = 2, xlim = c(-2, 6), ylim = c(-2,6), xlab = "observed value", ylab = "predicted value")
abline(a=0, b=1)

R2test.PLS <- 1 - sum((preprocessed.y.test - pls.predicted.y.test)^2) / sum((mean(preprocessed.y.test) - pls.predicted.y.test)^2)
R2test.PLS

