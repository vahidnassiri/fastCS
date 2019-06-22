#' combines result from subsamples using the approximate optimal 
#' weights together with their proper variances
#' @param nk a vector containing the cluster sizes
#' @param ck a vector containing the number of clusters of size nk
#' @param mu.split.est a vector containing the \eqn{\widehat{\mu}}'s from all sub-samples.
#' @param sigma2.split.est a vector containing the \eqn{\widehat{\sigma}^2}'s from all sub-samples
#' @param d.split.est a vector containing the \eqn{\widehat{d}}'s from all sub-samples
#' @param sigma2.split.var a vector containing the \eqn{\mathrm{Var}(\widehat{\sigma}^2)}'s from all sub-samples
#' @param d.split.var a vector containing the \eqn{\mathrm{Var}(\widehat{d})}'s from all sub-samples
#' @return a list with computed materials 
#' 
#' @author Vahid Nassiri
#' @export
approx.optimal.CS <- function(nk,ck,mu.split.est,sigma2.split.est,d.split.est,
            sigma2.split.var,d.split.var){
  library(magic)
  num.split=length(ck)
  #Calculating approximated optimal weights
  W=NULL
  for (i in 1:num.split){
    V1=(2*(sigma2.split.est[i]^2))/(ck[i]*(nk[i]-1))
    V2=(-1)*((2*(sigma2.split.est[i]^2))/(ck[i]*nk[i]*(nk[i]-1)))
    V3=(2/(ck[i]*nk[i]))*((sigma2.split.est[i]^2
                           /(nk[i]-1))+((2*d.split.est[i])
                                        *sigma2.split.est[i])+(nk[i]*(d.split.est[i]^2)))
    W[[i]]=solve(matrix(c(V1,V2,V2,V3),2,2))
  }
  V.total=apply(simplify2array(W),c(1,2),sum)
  W.inv=solve(V.total)
  W.opt=NULL
  sigma2.d=rbind(sigma2.split.est,d.split.est)
  sigma2.d.est=matrix(0,dim(sigma2.d)[1],dim(sigma2.d)[2])
  for (i in 1:num.split){
    W.opt=W.inv %*% W[[i]]
    sigma2.d.est[,i]=W.opt%*% sigma2.d[,i]
  }
  varcomp.est=apply(sigma2.d.est,1,sum)
  varcomp.var=W.inv
  # Calculating proper Variance
  A1=((sigma2.split.est^2) + (nk*d.split.est))^3
  A2=(sigma2.split.est*d.split.est)*((2*sigma2.split.est)
                                     *(nk*d.split.est))
  A3=d.split.est*((sigma2.split.est+(nk*d.split.est))^2)
  A=A1-A2-A3
  dev.w.sigma2=NULL
  dev.w.d=NULL
  C3.tmp=NULL
  C4.tmp=NULL
  for (i in 1:num.split){
    tmp1=(-1*ck[i]*nk[i])/A1[i]
    tmp2=matrix(c(A[i] ,1 ,1, nk[i]),2,2)
    dev.w.sigma2[[i]]=tmp1*tmp2
    dev.w.d[[i]]=tmp1 * matrix(c(1,nk[i],nk[i],(nk[i]^2)),2,2)
    
    II=as.matrix(diag(2)-(W.inv%*%W[[i]]))
    Theta=c(sigma2.split.est[i],d.split.est[i])
    Rho=II%*%Theta
    C3.tmp[[i]]=adiag(dev.w.d[[i]],dev.w.sigma2[[i]])
    C4.tmp[[i]]= kronecker (diag(2),Rho)
  }
  C1.tmp=t(kronecker(c(1,1),diag(2)))
  C2.tmp=kronecker(diag(2),W.inv)
  proper.var=NULL
  for (i in 1:num.split){
    C=((C1.tmp%*%C2.tmp)%*%C3.tmp[[i]])%*%C4.tmp[[i]]
    AA=W.inv%*%(W[[i]])
    B=AA+C
    VV=diag(c(d.split.var[i],sigma2.split.var[i]))
    proper.var[[i]]=(B%*%VV)%*%t(B)
  }
  Proper.var=apply(simplify2array(proper.var),c(1,2),sum)
  
  # Approximate weights for mu
  Am=(ck*nk)/(sigma2.split.est+(nk*d.split.est))
  
  w.mu=Am/sum(Am)
  mu.est=sum(w.mu*mu.split.est)
  mu.var=1/sum(Am)
  # Proper variance for mu
  proper.var.mu1=mu.var
  TMP2=rep(0,num.split)
  TMP.1=rep(0,num.split)
  for (k in 1:num.split){
    TMP1=2*ck[k]*(nk[k]^2)
    for (m in 1:num.split){
      TMP2[m]=Am[m]*((mu.split.est[k]-mu.split.est[m])^2)
    }
    TMP.1[k]=TMP1*sum(TMP2)
  }
  TMP.2=sum(Am)^4
  proper.var.mu=sum(TMP.1)/TMP.2
  Proper.var.mu=(1/sum(Am))+ proper.var.mu
  return(list(mu.est=mu.est,varcomp.est=varcomp.est,mu.var=mu.var,
              varcomp.var=varcomp.var,proper.var.mu=Proper.var.mu,
              proper.var.varcomp=Proper.var))
}