#' computes the parameters estimates with their variances within each split
#' @param n an integer indicating the common cluster size within each split
#' @param C an integer indicating the number of clusters within each split
#' @param Y a vector containing the response values corresponding to each split
#' @return a list with computed materials 
#' 
#' @author Vahid Nassiri
#' @export
est.CS.for <- function(n,C,Y){
  y.matrix=matrix(Y,n,C)
  mu.hat=mean(Y)
  Z=Y-mu.hat
  Z.matrix=matrix(Z,n,C)
  tmp2=rep(0,C)
  tmp1=rep(0,C)
  J=matrix(1,n,n)
  for (i in 1:C){
    tmp1[i]=t(Z.matrix[,i])%*% Z.matrix[,i]
    tmp2[i]= (t(Z.matrix[,i])%*%J)%*%Z.matrix[,i]
  }
  tmp3=1/((C*n)*(n-1))
  sigma2.hat=tmp3*((n*sum(tmp1))-sum(tmp2))
  d.hat=tmp3*(sum(tmp2)-sum(tmp1))
  var.mu.hat=(sigma2.hat+(n*d.hat))/(C*n)
  cov.varcomp=(2*(sigma2.hat^2)/((C*n)*(n-1)))*matrix(
    c(n,-1,-1,(((sigma2.hat^2) + ((2*(n-1))*(d.hat*sigma2.hat)) + 
                  ((n*(n-1))*(d.hat^2)))/(sigma2.hat^2))),2,2)
  return(list(mu.hat=mu.hat,d.hat=d.hat,sigma2.hat=sigma2.hat,
              var.mu.hat=var.mu.hat,cov.varcomp=cov.varcomp))
}
