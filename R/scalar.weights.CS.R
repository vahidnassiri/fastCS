#' combine results from different sub-samples using scalar weights
#' @param nk a vector containing the cluster sizes 
#' @param ck a vector containing the number of clusters of size nk
#' @param mu.split.est a vector containing the \eqn{\widehat{\mu}}'s from all sub-samples.
#' @param sigma2.split.est a vector containing the \eqn{\widehat{\sigma}^2}'s from all sub-samples
#' @param d.split.est a vector containing the \eqn{\widehat{d}}'s from all sub-samples
#' @param mu.split.var a vector containing the variance of \eqn{\widehat{\mu}}'s from all sub-samples
#' @param sigma2.split.var a vector containing the \eqn{\mathrm{Var}(\widehat{\sigma}^2)}'s from all sub-samples
#' @param d.split.var a vector containing the \eqn{\mathrm{Var}(\widehat{d})}'s from all sub-samples
#' @return a list with computed materials 
#'  
#' @author Vahid Nassiri
#' @export
scalar.weights.CS <- function(nk,ck,mu.split.est,sigma2.split.est,d.split.est,
            mu.split.var,sigma2.split.var,d.split.var){
  
  ak=(ck*nk)/(sigma2.split.est + (nk*d.split.est))
  w.mu=ak/sum(ak)
  bk=ck*(nk-1)
  w.sigma2=bk/sum(bk)	
  gk=(ck*nk)/
    (((sigma2.split.est^2)/(nk-1))
     +((2*sigma2.split.est)*d.split.est)+(nk*(d.split.est^2)))
  w.d=gk/sum(gk)
  
  mu.scalar=c(sum(w.mu*mu.split.est), sum(mu.split.var* (w.mu^2)))
  d.scalar=c(sum(w.d*d.split.est), sum(d.split.var* (w.d^2)))
  sigma2.scalar=c(sum(w.sigma2*sigma2.split.est),
                  sum(sigma2.split.var* (w.sigma2^2)))
  param.scalar=rbind(mu.scalar,sigma2.scalar,d.scalar)
  colnames(param.scalar)=c("Est.","Var")
  rownames(param.scalar)=c("mu","sigma2","d")
  return(param.scalar)
}