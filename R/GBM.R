#' GBM
#' @description Calculates the Genomic Estimated Breeding Value based on GBM method.
#' @param X X is a matrix of marker genotype of size n×p where n are no of Individuals under study (i.e. genotype, lines) and p are number of markers.
#' @param Y Y is a vector of individuals of size n×1.
#' @param r fraction of testing data (ranges from (0 to 1)) used during model fitting (suppose if one want to use 0.9 of data for model training and remaining 0.1 for model testing so one has to define r=0.1).
#' @details This function fits model by dividing data into two part i.e. training population and testing population. Former one is used to build the models and later one for performance evaluation. The performance of model is evaluated by calculating model's prediction accuracy of both training and testing population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value). Whole procedures is repeated 10 times and accuracy is averaged.For machine learning methods, an auto hyperparameter adjusting through grid search is implemented to achieve better prediction accuracy. In each split, the 0.9 of the training population is used to adjust hyperparameters, of which 0.9 selected to be training sets and 0.1 to be testing sets, 10 times of simple validation is implemented.
#' @return $fit GBM model fitting
#'        $Accuracy prediction accuracy of the testing population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value)
#'        $Train accuracy prediction accuracy of the training population (pearson correlation coefficient between actual phenotypic value and predicted phenotypic value)
#' @import gbm
#' @export GBM
#' @examples library(ASGS);rm (list=ls ());data(trout);
#' X=trout[,2:6];Y=as.data.frame(trout[,1]);r=0.1;ASGS::GBM(X,Y,r)
GBM<-function(X, Y, r){
  n<-nrow(X)
  p<-ncol(X)
  m<-round(n*r)
  correl<-c();correl2<-c()
  #library(gbm)
  if (!requireNamespace("gbm", quietly = TRUE)) {
    install.packages('gbm')}
  requireNamespace('gbm')

  gridsearch<-function(X, Y){
    n<-nrow(X)
    p<-ncol(X)
    m<-round(n*0.1)
    correl<-c();Accuracy<-c()
    pra<-data.frame(c(0.01,250,2),c(0.01,250,5),c(0.01,250,8),c(0.01,500,2),c(0.01,500,5),c(0.01,500,8),c(0.01,1000,2),c(0.01,1000,5),c(0.01,1000,8),c(0.05,250,2),c(0.05,250,5),c(0.05,250,8),c(0.05,500,2),c(0.05,500,5),c(0.05,500,8),c(0.05,1000,2),c(0.05,1000,5),c(0.05,1000,8),c(0.25,250,2),c(0.25,250,5),c(0.25,250,8),c(0.25,500,2),c(0.25,500,5),c(0.25,500,8),c(0.25,1000,2),c(0.25,1000,5),c(0.25,1000,8))
    for (i in 1:27){
      for (k in 1:10){
        tst<-sample(1:n,size=m,replace=FALSE)

        XS<-scale(X,center=T,scale=T)

        XTRN<-XS[-tst,] ; YTRN<-Y[-tst]
        XTST<-XS[tst,] ; YTST<-Y[tst]

        fm<-gbm.fit(x=XTRN,y=YTRN,distribution="gaussian",shrinkage=pra[1,i],n.trees=pra[2,i],interaction.depth=pra[3,i],verbose=F)

        Predic<-predict(fm,XTST)

        correl[k]<-cor(as.numeric(Predic),as.numeric(YTST))
      }
      Accuracy[i]<-mean(correl, na.rm=TRUE)
    }
    pra2<-pra[,which.max(Accuracy)]
    return(c(pra2))
  }

  for (k in 1:10){
    tst<-sample(1:n,size=m,replace=FALSE)

    XS<-scale(X,center=T,scale=T)
    XTRN<-XS[-tst,] ; YTRN<-Y[-tst,]
    XTST<-XS[tst,] ; YTST<-Y[tst,]

    pa<-gridsearch(XTRN,YTRN)

    fm<-gbm.fit(x=XTRN,y=YTRN,distribution="gaussian",shrinkage=pa[1],n.trees=pa[2],interaction.depth=pa[3],verbose=F)

    Predic<-predict(fm,XTST)
    Predic2<-predict(fm,XTRN)

    correl[k]<-cor(as.numeric(Predic),as.numeric(YTST))
    correl2[k]<-cor(as.numeric(Predic2),as.numeric(YTRN))
  }
  Accuracy<-mean(correl, na.rm=TRUE)
  Accuracy2<-mean(correl2, na.rm=TRUE)
  return(list("fit"=fm, "Train-accuracy"=Accuracy2,"Accuracy"=Accuracy))
}
