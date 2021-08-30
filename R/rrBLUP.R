#' rrBLUP
#' @description Calculates the Genomic Estimated Breeding Value based on rrBLUP method.
#' @param X X is a matrix of marker genotype of size n×p where n are no of Individuals under study (i.e. genotype, lines) and p are number of markers.
#' @param Y Y is a vector of individuals of size n×1.
#' @param r fraction of testing data (ranges from (0 to 1)) used during model fitting (suppose if one want to use 0.9 of data for model training and remaining 0.1 for model testing so one has to define r=0.1).
#' @details This function fits model by dividing data into two part i.e. training population and testing population. Former one is used to build the models and later one for performance evaluation. The performance of model is evaluated by calculating model's prediction accuracy of both training and testing population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value). Whole procedures is repeated 50 times and accuracy is averaged.
#' @return $Vu variance of random effect
#'        $Ve error variance
#'        $beta estimate of fixed effects
#'        $u estimate of random effects
#'        $LL maximized log likelihood
#'        $Accuracy prediction accuracy of the testing population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value)
#'        $Train accuracy prediction accuracy of the training population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value)
#' @import rrBLUP
#' @export RRBLUP
#' @examples library(ASGS);rm (list=ls ());data(trout);
#' X=trout[,2:101];Y=as.data.frame(trout[,1]);r=0.1;ASGS::RRBLUP(X,Y,r)
RRBLUP<-function(X, Y, r){
  n<-nrow(X)
  p<-ncol(X)
  m<-round(n*r)
  #library(rrBLUP)
  if (!requireNamespace("rrBLUP", quietly = TRUE)) {
    install.packages('rrBLUP')}
  requireNamespace('rrBLUP')
  correl<-c();correl2<-c()

  for (k in 1:50){
    tst<-sample(1:n,size=m,replace=FALSE)

    XTRN<-X[-tst,] ; YTRN<-Y[-tst,]
    XTST<-X[tst,] ; YTST<-Y[tst,]



    fm<-mixed.solve(y=YTRN,Z=as.matrix(XTRN))
    mu<-rep(fm$beta,length(YTST))
    Pred1<-(as.matrix(XTST)%*%as.vector(fm$u))+mu
    Pred2<-(as.matrix(XTRN)%*%as.vector(fm$u))+mu

    correl[k]<-cor(Pred1,YTST)
    correl2[k]<-cor(Pred2,YTRN)
  }
  Accuracy<-mean(correl, na.rm=TRUE)
  Accuracy2<-mean(correl2, na.rm=TRUE)
  return(list("Vu"= fm$Vu, "Ve"=fm$Ve, "beta"=fm$beta,"u"= fm$u, "LL"= fm$LL,"Train-accuracy"=Accuracy2, "Accuracy"= Accuracy))
}
