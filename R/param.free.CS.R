#' computes the three parameter free combining rules: Prop., Equal, and Appr.sc.
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
param.free.CS <- function(nk,ck,mu.split.est,sigma2.split.est,
                          d.split.est,mu.split.var,sigma2.split.var,
                          d.split.var){
  num.split=length(ck)
  # Calculating parameter free weights
  Equal=rep(1/num.split,num.split)
  Prop=ck/sum(ck)
  Appr.sc=(ck*nk)/sum(ck*nk)
  # mu
  mu.eq=c(sum(Equal*mu.split.est),sum((Equal^2)*mu.split.var))
  mu.prop=c(sum(Prop*mu.split.est),sum((Prop^2)*mu.split.var))
  mu.appr.sc=mu.prop
  mu=rbind(mu.eq,mu.prop,mu.appr.sc)
  colnames(mu)=c("Est","Var")
  rownames(mu)=c("Equal","Prop","Appr.sc.")
  # Sigma2
  sigma2.eq=c(sum(Equal*sigma2.split.est),
              sum((Equal^2)*sigma2.split.var))
  sigma2.prop=c(sum(Prop*sigma2.split.est),
                sum((Prop^2)*sigma2.split.var))
  sigma2.appr.sc=c(sum(Appr.sc*sigma2.split.est),
                   sum((Appr.sc^2)*sigma2.split.var))
  sigma2=rbind(sigma2.eq,sigma2.prop,sigma2.appr.sc)
  colnames(sigma2)=c("Est","Var")
  rownames(sigma2)=c("Equal","Prop","Appr.sc.")
  # d
  d.eq=c(sum(Equal*d.split.est),sum((Equal^2)*d.split.var))
  d.prop=c(sum(Prop*d.split.est),sum((Prop^2)*d.split.var))
  d.appr.sc=d.prop
  d=rbind(d.eq,d.prop,d.appr.sc)
  colnames(d)=c("Est","Var")
  rownames(d)=c("Equal","Prop","Appr.sc.")
  
  return(list(mu=mu,sigma2=sigma2,d=d))
}