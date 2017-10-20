require(leaps)

part<-function(d, outcomeVariable, splineTerm, additionalVars = NULL, K){

  #variable assignment

  x <- as.matrix(d[splineTerm])
    colnames(x) <- 'x'
  x.2 <- x^2
    colnames(x.2) <- 'x.2'
  x.3 <- x^3
    colnames(x.3) <- 'x.3'

  y <- as.matrix(d[outcomeVariable])
  L <- floor(length(y)/K)

  #Min/Max algorithm to find absolute maximum deviate in the kth partition

  get_x <- function(k){
    me <- mean(y[(L*(k-1)+1):(k*L)])
    ym <- y[(L*(k-1)+1):(k*L)] - me
    mc <- abs(max(ym))
    l <- do.call('c',
                 lapply(1:L, function(i)
                   ifelse(ym[i] == mc, x[(L*(k-1)+1):(k*L)][i], 0)))
    return(l[l!=0])}

  pops <- do.call('c',
                  lapply(1:K, function(i)
                    get_x(i)))

  #matrix of potential spline knots

  TT <- matrix(nrow = length(x), ncol =K)
  for (i in 1:length(x))
    for (k in 1:K)
      TT[i,k] <- ifelse(0 > (x[i]-pops[k])^3 , 0 ,(x[i]-pops[k])^3)
      bigT <- data.frame(TT)
      colnames(bigT) <- paste("X", pops,sep='')

  #model selection process

  if(is.null(additionalVars) == TRUE){
    X <- cbind(x,x.2,x.3,bigT)
  }
  else{
    X <- cbind(x,x.2,x.3,bigT,d[additionalVars])
  }

  r <- regsubsets(X,y, method="exhaustive", nbest=1, nmax=(K+4))
  cfs <- coef(r,id = match(max(summary(r)$adjr2), summary(r)$adjr2))

  xhat <- X[intersect(names(X), setdiff(names(cfs), "(Intercept)"))]
  xhat$'(Intercept)' <- rep(1, nrow(xhat))
  fits <- as.matrix(xhat[names(cfs)]) %*% as.matrix(cfs)

 out <- list()
 out$fits <- fits
 out$xhat <- xhat
 out$coefs <- cfs
 out$adjr2 <- max(summary(r)$adjr2)

 return(out)
}


