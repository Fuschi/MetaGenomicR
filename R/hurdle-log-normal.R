#'@author Alessandro Fuschi
#'
#'@importFrom truncnorm dtruncnorm
#'@importFrom MASS fitdistr
#'
#'@export
hlnorm.mle <- function(x, log=TRUE, warning.silent=TRUE){

	if(any(x<0)) stop("find negative values in x")
	if(!is.vector(x)) stop("x must be a vector")

	if(log) x <- log(x+1)

	x1 <- x[x>0]

	n <- length(x)
	n1 <- length(x1)
	n0 <- n - n1

	phi <- n0/n

	param.trnorm.mle <- NULL
	try(param.trnorm.mle <- MASS::fitdistr(x=x1, densfun=truncnorm::dtruncnorm,
										   start=list(mean=mean(x1),sd=sd(x1)),
										   a=0,b=Inf),
		silent=warning.silent
	)

	if(!is.null(param.trnorm.mle) && phi!=0){
		loglik <- n0*log(phi) + n1*log(1-phi) + param.trnorm.mle$loglik
	}else if (!is.null(param.trnorm.mle) && phi==0){
		loglik <- param.trnorm.mle$loglik
	}else{
		loglik <- NA
	}


	if(!is.null(param.trnorm.mle)){
		return(c("phi"    =phi,
				 "meanlog"=as.numeric(param.trnorm.mle$estimate[1]),
				 "sdlog"  =as.numeric(param.trnorm.mle$estimate[2]),
				 "loglik" =loglik,
				 "success"= 1)
		)
	} else {
		return(c("phi"    =phi,
				 "meanlog"=NA,
				 "sdlog"  =NA,
				 "loglik" =NA,
				 "success"=0)
		)
	}

}




#'@author Alessandro Fuschi
#'
#'@importFrom truncnorm dtruncnorm
#'
#'@export
dhlnorm <- function(x, phi, meanlog, sdlog){

	#CHECK ARGUMENTS
	if(phi<0 || phi>1) stop("phi must be in range [0,1]")
	if(sdlog<0) stop("sdlog must be greater than 0")

	y <- (1-phi)*dtruncnorm(xi,a=0,mean=mean,sd=sd)
	y[1] <- phi

	return(y)
}




#'@author Alessandro Fuschi
#'
#'@importFrom truncnorm rtruncnorm
#'@importFrom stats rbinom
#'
#'@export
rhlnorm <- function(N, phi, meanlog, sdlog){

	if(ceiling(N)!=N) stop("N must be an integer")
	if(N<1) stop("N must be a positive integer number")
	if(phi<0 || phi>1) stop("phi must be a number in [0,1]")
	if(sdlog<0) stop("sdlog must be a positive number")

	ans <- stats::rbinom(n=N, size=1, prob=1-phi)
	m <- length(ans[ans>0])
	ans[ans==1] <- truncnorm::rtruncnorm(n=m, a=0, b=Inf,
										 mean=meanlog,sd=sdlog)

	ans <- exp(ans) - 1
	ans[intersect(which(ans>0),which(ans<1))] <- 1
	ans <- round(ans)

	return(ans)
}
