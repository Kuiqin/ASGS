#' SVM
#' @description Calculates the Genomic Estimated Breeding Value based on SVM method.
#' @param X X is a matrix of marker genotype of size n×p where n are no of Individuals under study (i.e. genotype, lines) and p are number of markers.
#' @param Y Y is a vector of individuals of size n×1.
#' @param r fraction of testing data (ranges from (0 to 1)) used during model fitting (suppose if one want to use 0.9 of data for model training and remaining 0.1 for model testing so one has to define r=0.1).
#' @details This function fits model by dividing data into two part i.e. training population and testing population. Former one is used to build the models and later one for performance evaluation. The performance of model is evaluated by calculating model's prediction accuracy of both training and testing population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value). Whole procedures is repeated 50 times and accuracy is averaged. For machine learning methods, an auto hyperparameter adjusting through grid search is implemented to achieve better prediction accuracy. In each split, the 0.9 of the training population is used to adjust hyperparameters, of which 0.9 selected to be training sets and 0.1 to be testing sets, 10 times of simple validation is implemented.
#' @return $fit SVM model fitting
#'        $Accuracy prediction accuracy of the testing population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value)
#'        $Train accuracy prediction accuracy of the training population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value)
#' @import kernlab
#' @export SVM
#' @examples library(ASGS);rm (list=ls ());data(trout);
#' X=trout[1:50,2:11];Y=as.data.frame(trout[1:50,1]);r=0.1;ASGS::SVM(X,Y,r)
SVM<-function(X, Y, r){
  n<-nrow(X)
  p<-ncol(X)
  m<-round(n*r)

  correl<-c();correl2<-c()
  #library(kernlab)
  if (!requireNamespace("kernlab", quietly = TRUE)) {
    install.packages('kernlab')}
  requireNamespace('kernlab')

  gridsearch<-function(X, Y){
    n<-nrow(X)
    p<-ncol(X)
    m<-round(n*0.1)
    correl<-c();Accuracy<-c()
    pra<-c(0.01,0.001,0.05)
    for (i in 1:3){
      for (k in 1:10){
        tst<-sample(1:n,size=m,replace=FALSE)

        XTRN<-X[-tst,] ; YTRN<-Y[-tst]
        XTST<-X[tst,] ; YTST<-Y[tst]

        fm<-ksvm(y=YTRN,x=as.matrix(XTRN),type="eps-svr",kernel="rbfdot",scale=T,epsilon=pra[i])
        Pred1<-predict(fm,XTST,type="response")
        correl[k]<-cor(Pred1,YTST)
      }
      Accuracy[i]<-mean(correl, na.rm=TRUE)
    }
    pra2<-pra[which.max(Accuracy)]
    return(pra2)
  }

  for (k in 1:50){
    tst<-sample(1:n,size=m,replace=FALSE)

    XTRN<-X[-tst,] ; YTRN<-Y[-tst,]
    XTST<-X[tst,] ; YTST<-Y[tst,]

    pa<-gridsearch(XTRN,YTRN)

    fm<-ksvm(y=YTRN,x=as.matrix(XTRN),type="eps-svr",kernel="rbfdot",scale=T,epsilon=pa)
    Pred1<-predict(fm,XTST,type="response")
    Pred2<-predict(fm,XTRN,type="response")
    correl[k]<-cor(Pred1,YTST)
    correl2[k]<-cor(Pred2,YTRN)
  }
  Accuracy<-mean(correl, na.rm=TRUE)
  Accuracy2<-mean(correl2, na.rm=TRUE)
  return(list("fit"=fm,"Train-accuracy"=Accuracy2,"Accuracy"=Accuracy))
}
