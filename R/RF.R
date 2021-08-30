#' RF
#' @description Calculates the Genomic Estimated Breeding Value based on RF method.
#' @param X X is a matrix of marker genotype of size n×p where n are no of Individuals under study (i.e. genotype, lines) and p are number of markers.
#' @param Y Y is a vector of individuals of size n×1.
#' @param r fraction of testing data (ranges from (0 to 1)) used during model fitting (suppose if one want to use 0.9 of data for model training and remaining 0.1 for model testing so one has to define r=0.1).
#' @details This function fits model by dividing data into two part i.e. training population and testing population. Former one is used to build the models and later one for performance evaluation. The performance of model is evaluated by calculating model's prediction accuracy of testing population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value). Whole procedures is repeated 50 times and accuracy is averaged. For machine learning methods, an auto hyperparameter adjusting through grid search is implemented to achieve better prediction accuracy. In each split, the 0.9 of the training population is used to adjust hyperparameters, of which 0.9 selected to be training sets and 0.1 to be testing sets, 10 times of simple validation is implemented.
#' @return $fit RF model fitting
#'        $Accuracy prediction accuracy of the testing population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value)
#'
#' @import randomForest
#' @export RF
#' @examples library(ASGS);rm (list=ls ());data(trout);
#' X=trout[1:50,2:11];Y=as.data.frame(trout[1:50,1]);r=0.1;ASGS::RF(X,Y,r)
RF<-function(X, Y, r){
  n<-nrow(X)
  p<-ncol(X)
  m<-round(n*r)
  correl<-c()
  #library(randomForest)
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    install.packages('randomForest')}
  requireNamespace('randomForest')

  gridsearch<-function(X, Y){
    n<-nrow(X)
    p<-ncol(X)
    m<-round(n*0.1)
    correl<-c();Accuracy<-c()
    pra<-data.frame(c(100,1),c(100,3),c(100,5),c(500,1),c(500,3),c(500,5),c(1000,1),c(1000,3),c(1000,5))
    for (i in 1:9){
      for (k in 1:10){
        tst<-sample(1:n,size=m,replace=FALSE)

        XTRN<-X[-tst,] ; YTRN<-Y[-tst]
        XTST<-X[tst,] ; YTST<-Y[tst]

        Pred1<-randomForest(x=XTRN, y=YTRN, xtest=XTST, ytst=YTST,importance=T,proximity=T,ntree=pra[1,i],nodesize=pra[2,i])
        Predic<-Pred1$test$predicted
        correl[k]<-cor(as.numeric(Predic),as.numeric(YTST))
      }
      Accuracy[i]<-mean(correl, na.rm=TRUE)
    }
    pra2<-pra[,which.max(Accuracy)]
    return(c(pra2))
  }

  for (k in 1:50){
    tst<-sample(1:n,size=m,replace=FALSE)

    XTRN<-X[-tst,] ; YTRN<-Y[-tst,]
    XTST<-X[tst,] ; YTST<-Y[tst,]

    pa<-gridsearch(XTRN,YTRN)

    Pred1<-randomForest(x=XTRN, y=YTRN, xtest=XTST, ytst=YTST,importance=T,proximity=T,ntree=pa[1],nodesize=pa[2])
    Predic<-Pred1$test$predicted
    correl[k]<-cor(as.numeric(Predic),as.numeric(YTST))
  }
  Accuracy<-mean(correl, na.rm=TRUE)
  return(list("fit"=Pred1,"Accuracy"=Accuracy))
}
