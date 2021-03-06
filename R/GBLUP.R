#' GBLUP
#' @description Calculates the Genomic Estimated Breeding Value based on GBLUP method.
#' @param X X is a matrix of marker genotype of size n×p where n are no of Individuals under study (i.e. genotype, lines) and p are number of markers.
#' @param Y Y is a vector of individuals of size n×1.
#' @param r fraction of testing data (ranges from (0 to 1)) used during model fitting (suppose if one want to use 0.9 of data for model training and remaining 0.1 for model testing so one has to define r=0.1).
#' @details This function fits model by dividing data into two part i.e. training population and testing population. Former one is used to build the models and later one for performance evaluation. The performance of model is evaluated by calculating model's prediction accuracy of both training and testing population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value). Whole procedures is repeated 50 times and accuracy is averaged.
#' @return $fit GBLUP model fitting
#'        $Accuracy prediction accuracy of the testing population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value)
#'        $Train accuracy prediction accuracy of the training population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value)
#' @import BGLR
#' @export GBLUP
#' @examples library(ASGS);rm (list=ls ());data(trout);
#' X=trout[1:20,2:9];Y=as.data.frame(trout[1:20,1]);r=0.1;ASGS::GBLUP(X,Y,r)
GBLUP<-function(X, Y, r){
  n<-nrow(X)
  p<-ncol(X)
  m<-round(n*r)
  correl<-c();correl2<-c()
  Y<-as.vector(Y[,1])
  #library(BGLR)
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    install.packages('BGLR')}
  requireNamespace('BGLR')
  XS<-scale(X,center=TRUE,scale=TRUE)
  G<-tcrossprod(XS)/ncol(XS)
  for (i in 1:50){
    yNA<-Y
    tst<-sample(1:n,size=m,replace=FALSE)

    yNA[tst]<-NA

    ETA<-list(list(K=G,model='RKHS'))

    fm<-BGLR(y=yNA,ETA=ETA,
             nIter=12000,burnIn=2000,df0=5,S0=NULL,saveAt='GBLUP_',verbose=F)

    correl[i]<-cor(fm$yHat[tst],Y[tst])
    correl2[i]<-cor(fm$yHat[-tst],Y[-tst])
  }
  Accuracy<-mean(correl, na.rm=TRUE)
  Accuracy2<-mean(correl2, na.rm=TRUE)
  return(list("fit"=fm, "Train-accuracy"=Accuracy2,"Accuracy"=Accuracy))
}
