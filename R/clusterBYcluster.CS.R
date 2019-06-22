#' computes weighted, two stage, unbiased two stage methods for cluster by cluster analyses
#' @param nk a vector containing the cluster sizes
#' @param ck a vector containing the number of clusters of size nk
#' @param Data a 3-column matrix with first column the subject, second 
#' column the response variable, and third column the split, which show which 
#' observation  belongs to which sub-sample
#' @return a list with computed materials 
#' 
#' @author Vahid Nassiri
#' @export
clusterBYcluster.CS <- function(nk,ck,Data){
	# Data should be a 3-column matrix with first column the subject,
	# second column the response and third column the split indexes
	# which show which observation belongs to which sub-sample.
	num.split=length(ck)
	Var.varcomp=NULL
	mu.split.est=rep(0,num.split)
	mu.split.var=matrix(0,num.split,3)
	d.split.est=matrix(0,num.split,3)
	sigma2.split.est=matrix(0,num.split,3)
	d.split.var=matrix(0,num.split,3)
	sigma2.split.var=matrix(0,num.split,3)
	
	for (k in 1:num.split){
		
		# Making data for each cluster
		split.data=Data[Data[,3]==k,]
		n=nk[k]
		N=ck[k]
		# Computing t^2
		data.matrix=matrix(split.data[,2],n,N)
		mu.hat=sum(apply(data.matrix,2,sum))/(prod(dim(data.matrix)))
		t2=sum((apply(data.matrix,2,mean)-mu.hat)^2)/(dim(data.matrix)[2])
		#Computing S^2
		mean.vec=apply(data.matrix,2,mean)
		SS=rep(0,dim(data.matrix)[2])
		for (i in 1:dim(data.matrix)[2]){
			SS[i]=sum((data.matrix[,i]-mean.vec[i])^2)
		}
		s2=sum(SS)/prod(dim(data.matrix))
		# Computing s^2* and t^2*
		s2.star=(dim(data.matrix)[1]/(dim(data.matrix)[1]-1))*s2
		t2.star=(dim(data.matrix)[2]/(dim(data.matrix)[2]-1))*t2
		# Computing \widehat{\sigma}^2 and \widehat{d}
		Z=data.matrix-mu.hat
		J=matrix(1,dim(Z)[1],dim(Z)[1])
		ZZ=rep(0,dim(Z)[2])
		ZJZ=rep(0,dim(Z)[2])
		for (i in 1:dim(Z)[2]){
			ZZ[i]=t(Z[,i])%*% Z[,i]
			ZJZ[i]=(t(Z[,i])%*%J)%*%Z[,i]
		}
		d.hat=(1/(N*n*(n-1))) * (sum(ZJZ)-sum(ZZ))
		sigma2.hat=(1/(N*n*(n-1)))* ((n*sum(ZZ))-sum(ZJZ))
		
		# computing the variance of \widehat{\mu} \widehat{\sigma}^2 and \widehat{d}
		
		var.mu=(sigma2.hat+ (n*d.hat))/(N*n)
		var1=(2*(sigma2.hat^2))/(N*(n-1))
		var12=(-1)*(2*(sigma2.hat^2))/(N*n*(n-1))
		var2=(2/(n*N))*(((sigma2.hat^2)/(n-1)) + (2*sigma2.hat*d.hat)
					+ (n*(d.hat^2)))
		var.varcomp=matrix(c(var1,var12,var12,var2),2,2)
		
		# computing the variance of \widehat{\mu} \widehat{\sigma}^2 and \widehat{d}
		#based on the two stage approach
		
		var.mu.2stage=(s2+ (n*t2))/(N*n)
		var.s2=(2*(n-1)*(s2^2))/(N*(n^2))
		var.t2=(2*((N-1)^2))/((N^2)*n) * (((s2^2)/n) + (2*s2+t2 ) 
					+ (n*(t2^2) ))
		
		# computing the variance of \widehat{\mu} \widehat{\sigma}^2 and \widehat{d}
		#based on the unbiased two stage approach
		
		var.mu.2stage.unbiased=(s2.star+ (n*t2.star))/(N*n)
		var.s2.star=(2*(s2.star^2))/ (N*(n-1))
		var.t2.star= (2/(N*n)) * (((s2.star^2)/n)+(2*s2.star*t2.star)
					+(n*(t2.star^2)))
		
		# Saving the results
		
		mu.split.est[k]=mu.hat
		Var.varcomp[[k]]=var.varcomp
		mu.split.var[k,]=c(var.mu,var.mu.2stage,var.mu.2stage.unbiased)
		d.split.est[k,]=c(d.hat,t2,t2.star)
		sigma2.split.est[k,]=c(sigma2.hat,s2,s2.star)
		d.split.var[k,]=c(var.varcomp[2,2],var.t2,var.t2.star)
		sigma2.split.var[k,]=c(var.varcomp[1,1],var.s2,var.s2.star)
	}
	colnames(mu.split.var)=c("Weighted","TwoStage",
			"TwoStageUnbiased")
	colnames(d.split.var)=c("Weighted","TwoStage",
			"TwoStageUnbiased")
	colnames(sigma2.split.var)=c("Weighted","TwoStage",
			"TwoStageUnbiased")
	colnames(d.split.est)=c("Weighted","TwoStage",
			"TwoStageUnbiased")
	colnames(sigma2.split.est)=c("Weighted","TwoStage","TwoStageUnbiased")
	
	return(list(mu.split.est=mu.split.est,mu.split.var=mu.split.var,
					sigma2.split.est=sigma2.split.est, sigma2.split.var=sigma2.split.var,
					d.split.est=d.split.est,d.split.var=d.split.var, 
					Var.varcomp=Var.varcomp))}