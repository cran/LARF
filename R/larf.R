#######################################
### LARF: Local Response Functions in R
### Date: 2013-02-22
### Version 1.0
#######################################


##############################
### Fucntion to implement LARF
##############################


larf <- function(Y, D, X, Z, outcome = "binary", intercept = TRUE, method = "NLS", discrete = TRUE) {

  y <- as.matrix(Y)
  d <- as.matrix(D)
  X <- as.matrix(X)
  z <- as.matrix(Z)
  
  n <- length(y)
  
  # The defaul adds an intercept
  if(intercept == 1) X <- cbind(rep(1,n),X)
    
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
        g <- t(cbind(d,X))%*%diag(as.vector(kappa*dnorm(cbind(d,X)%*%b0)),n)
        delta <- step*solve(g%*%t(g))%*%g%*%u # The Gauss-Newton Method
        b1 <- b0 + delta
        dg <- t(delta)%*%g%*%t(g)%*%delta
        k <- k + 1
      }
    }
    
    # The parameters
    b <- b1
    
    # Get the SE using theorem 4.2
    u <- y - pnorm(cbind(d,X)%*%matrix(b))
    
    dM_theta <- (cbind(d,X)%*%b)*dnorm(cbind(d,X)%*%b)*u+(dnorm(cbind(d,X)%*%b))^2
    M_theta <- t(cbind(d,X))%*%diag(as.vector(kappa*dM_theta),n)%*%cbind(d,X)
    
    derkappa <- rep(0,n)
    derkappa[z==1&d==0] <- 1/(Phi[z==1&d==0])^2
    derkappa[z==0&d==1] <- (-1)/(1-Phi[z==0&d==1])^2
    
    M_gamma <- -t(cbind(d,X))%*%diag(as.vector(dnorm(cbind(d,X)%*%b)*u*matrix(derkappa)*dnorm(X%*%gamma)),n)%*%X
    
    lambda <- rep(0,n)
    lambda[z==0] <- (-1)*dnorm(X[z==0,]%*%gamma)/(1-pnorm(X[z==0,]%*%gamma))
    lambda[z==1] <- dnorm(X[z==1,]%*%gamma)/pnorm(X[z==1,]%*%gamma)
    
    H_gamma <- t(X)%*%diag(lambda^2,n)%*%X
    
    # Get the variance
    S_1 <- t(cbind(d,X))%*%diag(as.vector((dnorm(cbind(d,X)%*%b)*u*kappa)^2),n)%*%cbind(d,X)
    S_2 <- M_gamma%*%solve(H_gamma)%*%t(M_gamma)
    S_3 <- M_gamma%*%solve(H_gamma)%*%t(X)%*%diag(as.vector(lambda*kappa*(-u)*dnorm(cbind(d,X)%*%matrix(b))),n)%*%cbind(d,X)
    S <- S_1 + S_2 + S_3 + t(S_3)
    
    # The variance and SE
    V <- solve(M_theta,tol=1e-21) %*% S %*% solve(M_theta,tol=1e-21)
    se <- as.matrix(sqrt(diag(V)))
    
    wbar <- apply(diag(kappa,n)%*%cbind(d,X),2,mean) / mean(kappa)
    
    db <- as.vector(dnorm(t(matrix(wbar))%*%b))*b
    
    G <- as.vector(dnorm(t(matrix(wbar))%*%b))*(diag(dim(cbind(d,X))[2])-as.vector(t(matrix(wbar))%*%b)*b%*%wbar)
    
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
          db[i] <- pnorm(t(matrix(wbar1))%*%b) - pnorm(t(matrix(wbar0))%*%b)
          g <- as.vector(dnorm(t(matrix(wbar1))%*%b))*wbar1 - as.vector(dnorm(t(matrix(wbar0))%*%b))*wbar0
          dse[i] <- sqrt(g%*%V%*%matrix(g))
        }
      }
        
      # Report output
      out <- cbind(b, se, db, dse)
      colnames(out) <- c( "b", "se", "db", "dse")
    } else {
      out <- cbind(b, se)
      colnames(out) <- c( "b", "se")
    }
    rownames(out)[1] <- "D"
    return(out)    
  } # End of the binary LARF
  
  
  if(outcome == "continuous") {
  
    # WLS  of the coefficients
    DK <- diag(as.vector(kappa), n)
    b <- solve(t(cbind(d,X)) %*% DK %*% cbind(d,X)) %*% t(cbind(d,X)) %*% DK %*% Y
    
    # Get the variance
    u <- y-(cbind(d,X)%*%matrix(b))

    M_theta <- t(cbind(d,X))%*%diag(as.vector(kappa),n)%*%cbind(d,X)
    
    derkappa <- replicate(n,0)
    derkappa[z==1&d==0] <- 1/(Phi[z==1&d==0])^2
    derkappa[z==0&d==1] <- (-1)/(1-Phi[z==0&d==1])^2

    M_gamma <- -t(cbind(d,X))%*%diag(as.vector(u*matrix(derkappa)*dnorm(X%*%gamma)),n)%*%X
    
    lambda <- replicate(n,0)
    lambda[z==0] <- (-1)*dnorm(X[z==0,]%*%gamma)/(1-pnorm(X[z==0,]%*%gamma))
    lambda[z==1] <- dnorm(X[z==1,]%*%gamma)/pnorm(X[z==1,]%*%gamma)
    H_gamma <- t(X)%*%diag(lambda^2,n)%*%X
    
    # 1st part of "a^2+b^2+2a*b" variance
    S_1 <- t(cbind(d,X))%*%diag(as.vector((u*kappa)^2),n)%*%cbind(d,X)
    
    # 2nd part
    S_2 <- M_gamma%*%solve(H_gamma)%*%t(M_gamma)
    
    # 3rd part, AKA, the a*b part
    S_3 <- M_gamma%*%solve(H_gamma)%*%t(X)%*%diag(as.vector(lambda*kappa*(-u)),n)%*%cbind(d,X)
    
    # Add them togather
    S <- S_1 + S_2 + S_3 + t(S_3)
    
    # The variance and SE
    V <- solve(M_theta,tol=1e-21)%*%S%*%solve(M_theta,tol=1e-21)
    se <- sqrt(diag(V))
    
    # Report output
    out <- cbind(b, se)
    colnames(out) <- c( "b", "se")
    rownames(out)[1] <- "D"
    return(out)
  } # End of the linear LARF
  
}