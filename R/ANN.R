#' ANN
#' @description Calculates the Genomic Estimated Breeding Value based on ANN method.
#' @param X X is a matrix of marker genotype of size n×p where n are no of Individuals under study (i.e. genotype, lines) and p are number of markers.
#' @param Y Y is a vector of individuals of size n×1.
#' @param r fraction of testing data (ranges from (0 to 1)) used during model fitting (suppose if one want to use 0.9 of data for model training and remaining 0.1 for model testing so one has to define r=0.1).
#' @details This function fits model by dividing data into two part i.e. training population and testing population. Former one is used to build the models and later one for performance evaluation. The performance of model is evaluated by calculating model's prediction accuracy of both training and testing population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value). Whole procedures is repeated 10 times and accuracy is averaged. For machine learning methods, an auto hyperparameter adjusting through grid search is implemented to achieve better prediction accuracy. In each split, the 0.9 of the training population is used to adjust hyperparameters, of which 0.9 selected to be training sets and 0.1 to be testing sets, 10 times of simple validation is implemented.
#' @return $fit ANN model fitting
#'        $Accuracy prediction accuracy of the testing population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value)
#'        $Train accuracy prediction accuracy of the training population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value)
#' @import brnn
#' @export ANN
#' @examples library(ASGS);rm (list=ls ());data(trout);
#'X=trout[1:50,2:11];Y=as.data.frame(trout[1:50,1]);r=0.1;ASGS::ANN(X,Y,r)
ANN<-function(X, Y, r){
  n<-nrow(X)
  p<-ncol(X)
  m<-round(n*r)
  X<-as.matrix(X)
  if (!requireNamespace("brnn", quietly = TRUE)) {
    install.packages('brnn')}
  requireNamespace('brnn')
  for(i in 1:ncol(X)){ (X[,i]<-X[,i]-mean(X[,i]))/sd(X[,i])}
  G<-tcrossprod(X)/ncol(X)

  correl<-c();correl2<-c()

  gridsearch<-function(X, Y){
    n<-nrow(X)
    p<-ncol(X)
    m<-round(n*0.1)
    correl<-c();Accuracy<-c()
    pra<-c(1,2,3)
    for (i in 1:3){
      for (k in 1:10){
        tst<-sample(1:n,size=m,replace=FALSE)
        GTRN<-X[-tst,] ; yTRN<-Y[-tst]
        GTST<-X[tst,] ; yTST<-Y[tst]
        NN<-brnn(y=yTRN,x=GTRN,neurons=pra[i], epochs=12,verbose=F)
        Pred1<- predict(NN, newdata=GTST)
        correl[k]<-cor(Pred1,yTST)
      }
      Accuracy[i]<-mean(correl, na.rm=TRUE)
    }
    pra2<-pra[which.max(Accuracy)]
    return(pra2)
  }

  for (k in 1:10){
    tst<-sample(1:n,size=m,replace=FALSE)

    GTRN<-G[-tst,] ; yTRN<-Y[-tst,]
    GTST<-G[tst,] ; yTST<-Y[tst,]
    pa<-gridsearch(GTRN,yTRN)

    NN<-brnn(y=yTRN,x=GTRN,neurons=pa, epochs=12,verbose=F)
    Pred2<- predict(NN, newdata=GTRN)
    Pred1<- predict(NN, newdata=GTST)
    correl[k]<-cor(Pred1,yTST)
    correl2[k]<-cor(Pred2,yTRN)
  }
  Accuracy<-mean(correl, na.rm=TRUE)
  Accuracy2<-mean(correl2, na.rm=TRUE)
  return(list("fit"=NN, "Train-accuracy"=Accuracy2,"Accuracy"=Accuracy))
}
