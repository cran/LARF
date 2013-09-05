#######################################
### LARF: Local Response Functions in R
### Date: 2013-09-02
### Version 1.1
#######################################

##############################
### The high-interface
##############################

larf <- function(formula, treatment, instrument, data, outcome = "binary", method = "NLS", discrete = TRUE) {
  
  ## set up a model.frame
  if(missing(data)) data <- environment(formula)
  mf <- model.frame(formula = formula, data = data)
  
  ## extract response, covaraites, treatment, and instrument
  Y <- model.response(mf)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  D <- treatment
  Z <- instrument
  
  if(is.null(outcome)) outcome <- "binary"
  if(is.null(method)) method <- "NLS"
  if(is.null(discrete)) discrete <- TRUE
  if(outcome == "continuous") discrete <- NULL
  
  ## call default interface
  est <- larf.fit(Y, X, D, Z, outcome, method, discrete) 

  est$call <-match.call()
  est$formula <- formula

  return(est)
}


##############################
### Function to implement LARF
##############################


larf.fit <- function(Y, X, D, Z, outcome, method, discrete) {
  y <- as.matrix(Y)
  X <- as.matrix(X)
  d <- as.matrix(D)
  z <- as.matrix(Z)
  
  n <- length(y)  
  
  if(missing(outcome)) outcome <- "binary"
  if(missing(method)) method <- "NLS"
  if(missing(discrete)) discrete <- TRUE
  
  # Get the weights kappa
  eqA <- glm(z ~ X - 1, family=binomial(probit))
  gamma <- as.matrix(summary(eqA)$coefficients[,1])
  Phi <- eqA$fitted
  
  nc_ratio <- rep(0,n)
  nc_ratio[d==0 & z==1] <- 1/Phi[d==0 & z==1]
  nc_ratio[d==1 & z==0] <- 1/(1-Phi[d==1 & z==0])
  kappa <- 1 - nc_ratio
  
  #  LARF for a binary outcome
  if(outcome == "binary") {
    
    if (method != "NLS") { method <- method } # In the future, MLE may be implemented.
    
    # Get initial s of the parameters
    eqB <- glm(y ~ d + X - 1, family=binomial(probit))
    theta <- as.matrix(summary(eqB)$coefficients[,1])
    b1 <- theta
    
    ## Nonlinear least squares (Abadie 2003: equation 8)
    if (method == "NLS") {
      step <- 1
      gtol <- 0.0001
      dg <- 1
      k <- 1
      while( dg > gtol & k <= 1000){
        b0 <- b1
        u <- y - pnorm(cbind(d,X)%*%b0) # The residual
        g <- crossprod(cbind(d,X),diag(as.vector(kappa*dnorm(cbind(d,X)%*%b0)),n))
        delta <- step*solve(tcrossprod(g))%*%g%*%u # The Gauss-Newton Method
        b1 <- b0 + delta
        dg <- crossprod(delta,g)%*%crossprod(g,delta)
        k <- k + 1
      }
    }
    
    # The parameters
    b <- b1
    
    # Get the SE using theorem 4.2
    u <- y - pnorm(cbind(d,X)%*%matrix(b))
    
    dM_theta <- (cbind(d,X)%*%b)*dnorm(cbind(d,X)%*%b)*u+(dnorm(cbind(d,X)%*%b))^2
    M_theta <- crossprod(cbind(d,X),diag(as.vector(kappa*dM_theta),n))%*%cbind(d,X)
    
    derkappa <- rep(0,n)
    derkappa[z==1&d==0] <- 1/(Phi[z==1&d==0])^2
    derkappa[z==0&d==1] <- (-1)/(1-Phi[z==0&d==1])^2
    
    M_gamma <- -crossprod(cbind(d,X),diag(as.vector(dnorm(cbind(d,X)%*%b)*u*matrix(derkappa)*dnorm(X%*%gamma)),n))%*%X
    
    lambda <- rep(0,n)
    lambda[z==0] <- (-1)*dnorm(X[z==0,]%*%gamma)/(1-pnorm(X[z==0,]%*%gamma))
    lambda[z==1] <- dnorm(X[z==1,]%*%gamma)/pnorm(X[z==1,]%*%gamma)
    
    H_gamma <- t(X)%*%diag(lambda^2,n)%*%X
    
    # Get the variance
    S_1 <- crossprod(cbind(d,X),diag(as.vector((dnorm(cbind(d,X)%*%b)*u*kappa)^2),n))%*%cbind(d,X)
    S_2 <- M_gamma%*%tcrossprod(solve(H_gamma),(M_gamma))
    S_3 <- M_gamma%*%tcrossprod(solve(H_gamma),X)%*%diag(as.vector(lambda*kappa*(-u)*dnorm(cbind(d,X)%*%matrix(b))),n)%*%cbind(d,X)
    S <- S_1 + S_2 + S_3 + t(S_3)
    
    # The variance and SE
    V <- solve(M_theta,tol=1e-21) %*% S %*% solve(M_theta,tol=1e-21)
    se <- as.matrix(sqrt(diag(V)))
    
    wbar <- apply(diag(kappa,n)%*%cbind(d,X),2,mean) / mean(kappa)
    
    db <- as.vector(dnorm(crossprod(matrix(wbar),b)))*b
    
    G <- as.vector(dnorm(crossprod(matrix(wbar),b)))*(diag(dim(cbind(d,X))[2])-as.vector(crossprod(matrix(wbar),b))*b%*%wbar)
    
    dV <- G%*%V%*%t(G)
    dse <- as.matrix(sqrt(diag(dV)))
    
    cat <- floor((colSums(cbind(d,X)==1 | cbind(d,X)==0))/n) 
    
    # The coefficient for categorical covaraites are discrete changes, holding other variables at their means
    if (discrete == TRUE) {
      for(i in 1:length(cat)){
        if (cat[i]==1) { 
          wbar1 <- wbar
          wbar1[i] <- 1
          wbar0 <- wbar
          wbar0[i] <- 0
          db[i] <- pnorm(crossprod(matrix(wbar1),b)) - pnorm(crossprod(matrix(wbar0),b))
          g <- as.vector(dnorm(crossprod(matrix(wbar1),b)))*wbar1 - as.vector(dnorm(crossprod(matrix(wbar0),b)))*wbar0
          dse[i] <- sqrt(g%*%V%*%matrix(g))
        }
      }
      
      # Report output
      b<-as.matrix(b)
      se<-as.matrix(se)
      db<-as.matrix(db)
      dse<-as.matrix(dse)
      V<-as.matrix(V)
      fitted.values <- pnorm(as.matrix(cbind(D,X))%*%b)
      residuals <- Y - fitted.values
      
      call<-match.call()
      out <- list(call = call, coefficients = b, StdErr = se, MargEff = db, MargStdErr = dse, vcov = V, 
                  fitted.values  = fitted.values, residuals = residuals)
      
      rownames(out$coefficients)[1]<-c("Treatment")
      rownames(out$StdErr)<-rownames(out$coefficients)
      rownames(out$MargEff)<-rownames(out$coefficients)
      rownames(out$MargStdErr)<-rownames(out$coefficients)
      rownames(out$vcov)<-rownames(out$coefficients)
    } else {
      b<-as.matrix(b)
      se<-as.matrix(se)
      V<-as.matrix(V)
      fitted.values <- pnorm(as.matrix(cbind(D,X))%*%b)
      residuals <- Y - fitted.values
      
      call<-match.call()
      out <- list(call = call, coefficients = b, StdErr = se, vcov = V, 
                  fitted.values  = fitted.values, residuals = residuals)
      
      rownames(out$coefficients)[1]<-c("Treatment")
      rownames(out$StdErr)<-rownames(out$coefficients)
      rownames(out$vcov)<-rownames(out$coefficients)        
    }
    
    class(out)<-"larf"
    return(out)   
  } # End of the binary LARF
  
  
  if(outcome == "continuous") {
    
    # WLS  of the coefficients
    DK <- diag(as.vector(kappa), n)
    b <- solve(crossprod(cbind(d,X), DK) %*% cbind(d,X)) %*% crossprod(cbind(d,X), DK) %*% Y
    
    # Get the variance
    u <- y-(cbind(d,X)%*%matrix(b))
    
    M_theta <- crossprod(cbind(d,X),diag(as.vector(kappa),n))%*%cbind(d,X)
    
    derkappa <- replicate(n,0)
    derkappa[z==1&d==0] <- 1/(Phi[z==1&d==0])^2
    derkappa[z==0&d==1] <- (-1)/(1-Phi[z==0&d==1])^2
    
    M_gamma <- -crossprod(cbind(d,X),diag(as.vector(u*matrix(derkappa)*dnorm(X%*%gamma)),n))%*%X
    
    lambda <- replicate(n,0)
    lambda[z==0] <- (-1)*dnorm(X[z==0,]%*%gamma)/(1-pnorm(X[z==0,]%*%gamma))
    lambda[z==1] <- dnorm(X[z==1,]%*%gamma)/pnorm(X[z==1,]%*%gamma)
    H_gamma <- crossprod(X,diag(lambda^2,n))%*%X
    
    # 1st part of "a^2+b^2+2a*b" variance
    S_1 <- crossprod(cbind(d,X),diag(as.vector((u*kappa)^2),n))%*%cbind(d,X)
    
    # 2nd part
    S_2 <- M_gamma%*%tcrossprod(solve(H_gamma),M_gamma)
    
    # 3rd part, AKA, the a*b part
    S_3 <- M_gamma%*%tcrossprod(solve(H_gamma),X)%*%diag(as.vector(lambda*kappa*(-u)),n)%*%cbind(d,X)
    
    # Add them togather
    S <- S_1 + S_2 + S_3 + t(S_3)
    
    # The variance and SE
    V <- solve(M_theta,tol=1e-21)%*%S%*%solve(M_theta,tol=1e-21)
    se <- sqrt(diag(V))
    
    # Report output
    b<-as.matrix(b)
    se<-as.matrix(se)
    V<-as.matrix(V)
    fitted.values <- as.matrix(cbind(D,X)%*%b)
    residuals <- Y - fitted.values
    
    call<-match.call()
    out <- list(call = call, coefficients = b, StdErr = se, vcov = V, 
                fitted.values  = fitted.values, residuals = residuals)
    
    rownames(out$coefficients)[1]<-c("Treatment")
    rownames(out$StdErr)<-rownames(out$coefficients)
    rownames(out$vcov)<-rownames(out$coefficients)
    
    class(out)<-"larf"
    return(out)
  } # End of the linear LARF
}

##############################
### Generic methods
##############################

print.larf <- function(x, digits = 3, ...)
{  
  cat("Call:\n")
  print(x$call)
  
  est <- cbind(x$coefficients, x$StdErr)
  colnames(est) <- c("Coefficients", "Standard Errors")
  
  cat("\nEstimates:\n")
  print.default(format(est, digits = digits), quote = FALSE)
}


summary.larf <- function(object, ...)
{  
  if (is.null(object$call$outcome)) object$call$outcome <- "binary" 
  if (is.null(object$call$discrete)) object$call$discrete <- TRUE
  if (object$call$outcome == "continuous") object$call$discrete <- NULL
  
  if (object$call$outcome == "binary" & !is.null(object$call$discrete)) {
    TAP<-cbind(Estimate = coef(object),
               StdErr = object$StdErr,
               MargEff=object$MargEff,
               MargStdErr=object$MargStdErr)
    
    colnames(TAP)<-c("Estimate","StdErr","MargEff","MargStdErr")
    res <- list(call=object$call, coefficients=TAP)   
    
  } else {
    TAP<-cbind(Estimate = coef(object),
               StdErr = object$StdErr)
    
    colnames(TAP)<-c("Estimate","StdErr")
    res <- list(call=object$call, coefficients=TAP)   
  }
  
  class(res) <- "summary.larf"
  return(res)
}


print.summary.larf <- function(x, digits = 3, ...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\n")
  print.default(format(x$coefficients, digits = digits), quote = FALSE)
}


vcov.larf <- function(object, ...)
{  
  colnames(object$vcov) <- rownames(object$vcov)
  cat("\nCovariance Matrix for the Parameters:\n")
  object$vcov
}


predict.larf <- function(object, newCov, newTreatment, ...) {
  if(is.null(object$call$outcome)) {
    object$predicted.values <- as.matrix(cbind(newTreatment,newCov)%*%coef(object))
  } 
  else {
    object$predicted.values <- pnorm(as.matrix(cbind(newTreatment,newCov))%*%coef(object))    
  }
  colnames(object$predicted.values) <- "predicted.values"
  return(object$predicted.values )
}



