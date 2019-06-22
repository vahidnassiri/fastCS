#' estimates the parameters for all of the splits simultaneously
#' @param data a 2-column matrix with response variable as its first column 
#' and splits indicator as its second column. All clusters should 
#' be of equal size by placing NaN to make smaller clusters the 
#' same size as the largest one.
#' @param nk a vector containing the cluster sizes
#' @param ck a vector containing the number of clusters of size nk
#' @return a list with computed materials 
#' 
#' @author Vahid Nassiri
#' @export
est.CS.all <- function(data,ck,nk){
  
  Y.mat=matrix(data[,1],max(nk),dim(data)[1]/max(nk))
  split.idx.sub=matrix(data[,2],max(nk),dim(data)[1]/max(nk))[1,]
  split.matrix=matrix(0,length(ck),sum(ck))
  for (i in 1:length(ck)){
    split.matrix[i,split.idx.sub==i]=1
  }
  subj.mean=apply(Y.mat,2,sum,na.rm=T)
  split.mu.hat=(split.matrix%*%subj.mean)/(ck*nk)
  subj.mu.hat=t(split.matrix)%*%split.mu.hat
  mean.mat=matrix(rep(c(subj.mu.hat),max(nk)),max(nk),
                  length(subj.mean),byrow=T)
  Z.mat=(Y.mat-mean.mat)
  tmp1.1=matrix(rep(apply(Z.mat^2,2,sum,na.rm=T),length(ck)),
                length(ck),sum(ck),byrow=T)
  tmp1=apply(split.matrix*tmp1.1,1,sum)
  tmp2.1=matrix(rep(apply(Z.mat,2,sum,na.rm=T),length(ck)),
                length(ck),sum(ck),byrow=T)
  tmp2.2=split.matrix*tmp2.1
  Z.mat.zero=Z.mat
  Z.mat.zero[is.na(Z.mat)==TRUE]=0
  tmp2=apply(Z.mat.zero%*%t(tmp2.2),2,sum,na.rm=T)
  tmp3=1/((ck*nk)*(nk-1))
  split.sigma2.hat=tmp3* ((nk*tmp1)-tmp2)
  split.d.hat=tmp3*(tmp2-tmp1)
  split.var.mu.hat=(split.sigma2.hat + (nk*split.d.hat))/(ck*nk)
  varcomp.factor=(2*(split.sigma2.hat^2)/((ck*nk)*(nk-1)))
  split.var.d.hat=varcomp.factor*(((split.sigma2.hat^2) +
                                     ((2*(nk-1))*(split.d.hat*split.sigma2.hat)) + 
                                     ((nk*(nk-1))*(split.d.hat^2)))/(split.sigma2.hat^2))
  split.var.sigma2.hat=varcomp.factor*nk
  split.cov.d.sigma2.hat=-1*varcomp.factor
  return(list(mu.hat=split.mu.hat,sigma2.hat=split.sigma2.hat,
              d.hat=split.d.hat,var.mu.hat=split.var.mu.hat,
              var.d.hat=split.var.d.hat,
              var.sigma2.hat=split.var.sigma2.hat,
              cov.d.sigma2.hat=split.cov.d.sigma2.hat))
}